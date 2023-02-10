//
// Created by magnus on 1/31/23.
//


#include "CameraCalibration.h"
#include "Logger.h"
#include "stb_image.h"
#include <fstream>
#include <vector>
#include <iterator>
#include "opencv2/opencv.hpp"

#include "yaml-cpp/yaml.h"


void CameraCalibration::readResponseFunction(const std::filesystem::path &CRFPath) {
    // load the inverse response function
    std::ifstream f(CRFPath.c_str());
    std::string line;
    std::getline(f, line);
    std::istringstream l1i(line);
    std::vector<float> GInvvec = std::vector<float>(std::istream_iterator<float>(l1i), std::istream_iterator<float>());
    if (GInvvec.size() != 256) {
        Log::Logger::getInstance()->info(
                "PhotometricUndistorter: invalid format! got {} entries in first line, expected 256!",
                (int) GInvvec.size());
    }
    for (int i = 0; i < 256; i++)
        responseFunc.gInv[i] = GInvvec[i];

    float min = responseFunc.gInv[0];
    float max = responseFunc.gInv[255];
    for (float &i: responseFunc.gInv) i = 255.0f * (i - min) / (max - min);            // make it to 0..255 => 0..255.

    bool isGood = true;
    for (int i = 0; i < 255; i++) {
        if (responseFunc.gInv[i + 1] <= responseFunc.gInv[i]) {
            Log::Logger::getInstance()->info(
                    "PhotometricUndistorter: G invalid! it has to be strictly increasing, but it isn't!");
            isGood = false;
        }
    }
    if (isGood)
        Log::Logger::getInstance()->info("Successfully loaded camera inverse response function from: {}",
                                         CRFPath.c_str());
}

void CameraCalibration::readVignetteMap(const std::filesystem::path &vignetteImagePath) {

    int width = 0, height = 0, channels = 0;
    unsigned short *pixels = stbi_load_16(vignetteImagePath.c_str(), &width, &height, &channels, 0);
    if (!pixels)
        Log::Logger::getInstance()->error("Failed to load image {}", vignetteImagePath.string());
    if (channels > 1) {
        Log::Logger::getInstance()->error("Loaded vignette image with too many channels: {}", channels);

    }
    vignetteMap.d.resize(width * height * channels);
    vignetteMap.inv.resize(width * height * channels);
    vignetteMap.d.insert(vignetteMap.d.begin(), &pixels[0], &pixels[width * height * channels]);
    vignetteMap.width = width;
    vignetteMap.height = height;
    float maxValue = *max_element(std::begin(vignetteMap.d), std::end(vignetteMap.d)); // C++11
    for (int i = 0; i < vignetteMap.width * vignetteMap.height; i++) {
        vignetteMap.d[i] = vignetteMap.d[i] / maxValue;
        vignetteMap.inv[i] = 1.0f / vignetteMap.d[i];
    }
    /*
    cv::Mat im_color;
    cv::Mat vignetteScaled(vignetteMap.height, vignetteMap.width, CV_32FC1, vignetteMap.d.data());
    vignetteScaled.convertTo(vignetteScaled, CV_8UC1, 255);
    cv::applyColorMap(vignetteScaled, im_color, cv::COLORMAP_HOT);
    cv::imshow("ColorMap", im_color);
    */
    stbi_image_free(pixels);
    Log::Logger::getInstance()->info("Successfully read photometric calibration: {}", vignetteImagePath.c_str());
}

uint32_t CameraCalibration::setGlobalCalib(uint32_t w, uint32_t h, const Mat33 &K) {
    int wlvl = w;
    int hlvl = h;
    int pyrLevelsUsed = 1;
    while (wlvl % 2 == 0 && hlvl % 2 == 0 && wlvl * hlvl > 5000 && pyrLevelsUsed < PYR_LEVELS) {
        wlvl /= 2;
        hlvl /= 2;
        pyrLevelsUsed++;
    }
    Log::Logger::getInstance()->info("using pyramid levels 0 to {}. coarsest resolution: {} x {}!",
                                     pyrLevelsUsed - 1, wlvl, hlvl);
    if (wlvl > 100 && hlvl > 100) {
        Log::Logger::getInstance()->info("\n\n===============WARNING!===================\n "
                                         "using not enough pyramid levels.\n"
                                         "Consider scaling to a resolution that is a multiple of a power of 2.\n");
    }
    if (pyrLevelsUsed < 3) {
        Log::Logger::getInstance()->info("\\n===============WARNING!===================\\n \"\n"
                                         "\"I need higher resolution.\\n\"\n"
                                         "\"I will probably segfault.\\n");

        printf("\n");
    }

    wG[0] = w;
    hG[0] = h;
    KG[0] = K;
    fxG[0] = K(0, 0);
    fyG[0] = K(1, 1);
    cxG[0] = K(0, 2);
    cyG[0] = K(1, 2);
    KiG[0] = KG[0].inverse();
    fxiG[0] = KiG[0](0, 0);
    fyiG[0] = KiG[0](1, 1);
    cxiG[0] = KiG[0](0, 2);
    cyiG[0] = KiG[0](1, 2);

    for (int level = 1; level < pyrLevelsUsed; ++level) {
        wG[level] = w >> level;
        hG[level] = h >> level;

        fxG[level] = fxG[level - 1] * 0.5f;
        fyG[level] = fyG[level - 1] * 0.5f;
        cxG[level] = (cxG[0] + 0.5f) / ((int) 1 << level) - 0.5f;
        cyG[level] = (cyG[0] + 0.5f) / ((int) 1 << level) - 0.5f;

        KG[level] << fxG[level], 0.0, cxG[level], 0.0, fyG[level], cyG[level], 0.0, 0.0, 1.0;    // synthetic
        KiG[level] = KG[level].inverse();

        fxiG[level] = KiG[level](0, 0);
        fyiG[level] = KiG[level](1, 1);
        cxiG[level] = KiG[level](0, 2);
        cyiG[level] = KiG[level](1, 2);
    }
    return pyrLevelsUsed;

}


void CameraCalibration::processFrame(std::shared_ptr<VO::Frame> frame, float exposure_time, float factor) const{
    int wh = frame->inWidth * frame->inHeight;
    {
        for (int i = 0; i < wh; i++) {
            frame->dataRaw[i] = responseFunc.gInv[frame->pixels[i]];
        }

        for (int i = 0; i < wh; i++) {
            frame->dataRaw[i] *= vignetteMap.inv[i];
        }
        frame->abExposure = exposure_time;
        frame->timestamp = 0;
    }
}

void CameraCalibration::rectify(const float* in_data, float* out_data) const{

    for (int idx = outWidth * outHeight - 1; idx >= 0; idx--) {
        // get interp. values
        float xx = remapX[idx];
        float yy = remapY[idx];

        if (xx < 0)
            out_data[idx] = 0;
        else {
            // get integer and rational parts
            int xxi = xx;
            int yyi = yy;
            xx -= xxi;
            yy -= yyi;
            float xxyy = xx * yy;

            // get array base pointer
            const float *src = in_data + xxi + yyi * inWidth;

            // interpolate (bilinear)
            out_data[idx] = xxyy * src[1 + inWidth]
                            + (yy - xxyy) * src[inWidth]
                            + (xx - xxyy) * src[1]
                            + (1 - xx - yy + xxyy) * src[0];
        }
    }
}
void
CameraCalibration::undistort(std::shared_ptr<VO::Frame> frame, float exposure, double timestamp, float factor) const{
    if (frame->inWidth != inWidth || frame->inHeight != inHeight) {
        Log::Logger::getInstance()->error("Undistort::undistort: wrong image size ({}x{} instead of {}x{}) \n",
                                          frame->width, frame->height, inWidth, inHeight);
        exit(1);
    }

    frame->height = outHeight;
    frame->width = outWidth;
    frame->dataf.resize(frame->width * frame->height);
    frame->dataRawRectified.resize(frame->width * frame->height);

    rectify(frame->dataRaw.data(), frame->dataRawRectified.data());
    processFrame(frame, exposure, factor);
    rectify(frame->dataRaw.data(), frame->dataf.data());

}

bool CameraCalibration::readIntrinsicCalibration(const std::filesystem::path &intrinsicsFilePath) {

    std::map<std::string, std::vector<float> > data;

    YAML::Node config = YAML::LoadFile(intrinsicsFilePath.c_str());

    inputRes = config["input_resolution"].as<std::vector<size_t>>();
    kData = config["K"].as<std::vector<double>>();
    kDataNew = config["K_new"].as<std::vector<double>>();
    outputRes = config["output_resolution"].as<std::vector<size_t>>();

    inWidth = inputRes[0];
    inHeight = inputRes[1];
    outWidth = outputRes[0];
    outHeight = outputRes[1];

    // rescale and substract 0.5 offset.
    // the 0.5 is because I'm assuming the calibration is given such that the pixel at (0,0)
    // contains the integral over intensity over [0,0]-[1,1], whereas I assume the pixel (0,0)
    // to contain a sample of the intensity ot [0,0], which is best approximated by the integral over
    // [-0.5,-0.5]-[0.5,0.5]. Thus, the shift by -0.5.
    kData[0] = kData[0] * inWidth;
    kData[1] = kData[1] * inHeight;
    kData[2] = kData[2] * inWidth - 0.5;
    kData[3] = kData[3] * inHeight - 0.5;

    Log::Logger::getInstance()->info("Loaded camera matrix: fx:{} fy:{} cx:{} cy:{} omega:{}", kData[0], kData[1], kData[2], kData[3], kData[4]);

    kDataNew[0] = kDataNew[0] * outWidth;
    kDataNew[1] = kDataNew[1] * outHeight;
    kDataNew[2] = kDataNew[2] * outWidth - 0.5;
    kDataNew[3] = kDataNew[3] * outHeight - 0.5;

    Log::Logger::getInstance()->info("Loaded Rectified camera matrix: fx:{} fy:{} cx:{} cy:{} omega:{}", kDataNew[0], kDataNew[1], kDataNew[2], kDataNew[3], kDataNew[4]);

    K = Eigen::Matrix<double, 3, 3>::Zero();
    K.setIdentity();
    K(0, 0) = kDataNew[0];
    K(1, 1) = kDataNew[1];
    K(0, 2) = kDataNew[2];
    K(1, 2) = kDataNew[3];



    remapX.resize(outputRes[0] * outputRes[1]);
    remapY.resize(outputRes[1] * outputRes[0]);

    //makeOptimalK_crop();

    for (int y = 0; y < outHeight; y++)
        for (int x = 0; x < outWidth; x++) {
            remapX[x + y * outWidth] = x;
            remapY[x + y * outWidth] = y;
        }

    distortCoordinates(remapX.data(), remapY.data(), remapX.data(), remapY.data(), outHeight * outWidth);


    for (int y = 0; y < outHeight; y++)
        for (int x = 0; x < outWidth; x++) {
            // make rounding resistant.
            float ix = remapX[x + y * outWidth];
            float iy = remapY[x + y * outWidth];

            if (ix == 0) ix = 0.001;
            if (iy == 0) iy = 0.001;
            if (ix == inWidth - 1) ix = inWidth - 1.001;
            if (iy == inHeight - 1) ix = inHeight - 1.001;

            if (ix > 0 && iy > 0 && ix < inWidth - 1 && iy < inWidth - 1) {
                remapX[x + y * outWidth] = ix;
                remapY[x + y * outWidth] = iy;
            } else {
                remapX[x + y * outWidth] = -1;
                remapY[x + y * outWidth] = -1;
            }
        }


    Log::Logger::getInstance()->info("Successfully read intrinsic calibration file: {}", intrinsicsFilePath.c_str());
    return true;
}

void CameraCalibration::setFramePyramidInfo(VO::Frame *frame) const {
    frame->K = K;
    for (int i = 0; i < pyrLevelsUsed; ++i) {
        frame->KLevel[i] = KG[i];
        frame->KiLevel[i] = KiG[i];
        frame->fxLevel[i] = fxG[i];
        frame->fyLevel[i] = fyG[i];
        frame->cxLevel[i] = cxG[i];
        frame->cyLevel[i] = cyG[i];
        frame->widthLevel[i] = wG[i];
        frame->heightLevel[i] = hG[i];
    }
}


void CameraCalibration::makeOptimalK_crop() {
    printf("finding CROP optimal new model!\n");
    K.setIdentity();

    // 1. stretch the center lines as far as possible, to get initial coarse quess.
    float *tgX = new float[100000];
    float *tgY = new float[100000];
    float minX = 0;
    float maxX = 0;
    float minY = 0;
    float maxY = 0;

    for (int x = 0; x < 100000; x++) {
        tgX[x] = (x - 50000.0f) / 10000.0f;
        tgY[x] = 0;
    }
    distortCoordinates(tgX, tgY, tgX, tgY, 100000);
    for (int x = 0; x < 100000; x++) {
        if (tgX[x] > 0 && tgX[x] < inputRes[0] - 1) {
            if (minX == 0) minX = (x - 50000.0f) / 10000.0f;
            maxX = (x - 50000.0f) / 10000.0f;
        }
    }
    for (int y = 0; y < 100000; y++) {
        tgY[y] = (y - 50000.0f) / 10000.0f;
        tgX[y] = 0;
    }
    distortCoordinates(tgX, tgY, tgX, tgY, 100000);
    for (int y = 0; y < 100000; y++) {
        if (tgY[y] > 0 && tgY[y] < inputRes[1] - 1) {
            if (minY == 0) minY = (y - 50000.0f) / 10000.0f;
            maxY = (y - 50000.0f) / 10000.0f;
        }
    }
    delete[] tgX;
    delete[] tgY;

    minX *= 1.01;
    maxX *= 1.01;
    minY *= 1.01;
    maxY *= 1.01;


    printf("initial range: x: %.4f - %.4f; y: %.4f - %.4f!\n", minX, maxX, minY, maxY);



    // 2. while there are invalid pixels at the border: shrink square at the side that has invalid pixels,
    // if several to choose from, shrink the wider dimension.
    bool oobLeft = true, oobRight = true, oobTop = true, oobBottom = true;
    int iteration = 0;
    while (oobLeft || oobRight || oobTop || oobBottom) {
        oobLeft = oobRight = oobTop = oobBottom = false;
        for (int y = 0; y < outputRes[1]; y++) {
            remapX[y * 2] = minX;
            remapX[y * 2 + 1] = maxX;
            remapY[y * 2] = remapY[y * 2 + 1] = minY + (maxY - minY) * (float) y / ((float) outputRes[1] - 1.0f);
        }
        distortCoordinates(remapX.data(), remapY.data(), remapX.data(), remapY.data(), 2 * outputRes[1]);
        for (int y = 0; y < outputRes[1]; y++) {
            if (!(remapX[2 * y] > 0 && remapX[2 * y] < inputRes[0] - 1))
                oobLeft = true;
            if (!(remapX[2 * y + 1] > 0 && remapX[2 * y + 1] < inputRes[0] - 1))
                oobRight = true;
        }


        for (int x = 0; x < outputRes[0]; x++) {
            remapY[x * 2] = minY;
            remapY[x * 2 + 1] = maxY;
            remapX[x * 2] = remapX[x * 2 + 1] = minX + (maxX - minX) * (float) x / ((float) outputRes[0] - 1.0f);
        }
        distortCoordinates(remapX.data(), remapY.data(), remapX.data(), remapY.data(), 2 * outputRes[0]);


        for (int x = 0; x < outputRes[0]; x++) {
            if (!(remapY[2 * x] > 0 && remapY[2 * x] < inputRes[1] - 1))
                oobTop = true;
            if (!(remapY[2 * x + 1] > 0 && remapY[2 * x + 1] < inputRes[1] - 1))
                oobBottom = true;
        }


        if ((oobLeft || oobRight) && (oobTop || oobBottom)) {
            if ((maxX - minX) > (maxY - minY))
                oobBottom = oobTop = false;    // only shrink left/right
            else
                oobLeft = oobRight = false; // only shrink top/bottom
        }

        if (oobLeft) minX *= 0.995;
        if (oobRight) maxX *= 0.995;
        if (oobTop) minY *= 0.995;
        if (oobBottom) maxY *= 0.995;

        iteration++;


        printf("iteration %05d: range: x: %.4f - %.4f; y: %.4f - %.4f!\n", iteration, minX, maxX, minY, maxY);
        if (iteration > 500) {
            printf("FAILED TO COMPUTE GOOD CAMERA MATRIX - SOMETHING IS SERIOUSLY WRONG. ABORTING \n");
            exit(1);
        }
    }

    K(0, 0) = ((float) outputRes[0] - 1.0f) / (maxX - minX);
    K(1, 1) = ((float) outputRes[1] - 1.0f) / (maxY - minY);
    K(0, 2) = -minX * K(0, 0);
    K(1, 2) = -minY * K(1, 1);

}


void CameraCalibration::distortCoordinates(float *in_x, float *in_y, float *out_x, float *out_y, int n) const {
    float dist = kData[4];
    float d2t = 2.0f * tan(dist / 2.0f);



    // current camera parameters
    float fx = kData[0];
    float fy = kData[1];
    float cx = kData[2];
    float cy = kData[3];


    float ofx = K(0, 0);
    float ofy = K(1, 1);
    float ocx = K(0, 2);
    float ocy = K(1, 2);

    for (int i = 0; i < n; i++) {
        float x = in_x[i];
        float y = in_y[i];
        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;

        float r = sqrtf(ix * ix + iy * iy);
        float fac = (r == 0 || dist == 0) ? 1 : atanf(r * d2t) / (dist * r);

        ix = fx * fac * ix + cx;
        iy = fy * fac * iy + cy;

        out_x[i] = ix;
        out_y[i] = iy;
    }
}
