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
#include "CalibrationYaml.h"

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
        Log::Logger::getInstance()->info("Successfully loaded camera inverse response function from: {}", CRFPath.c_str());
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

void CameraCalibration::applyPhotometricCalibration(std::shared_ptr<VO::Frame> *frame) {
    for (int i = 0; i < (frame->get()->width * frame->get()->height); i++) {
        frame->get()->dataf[i] =
                responseFunc.gInv[static_cast<unsigned char>(frame->get()->pixels[i])] * vignetteMap.inv[i];
    }
}

uint32_t CameraCalibration::setGlobalCalib(uint32_t w, uint32_t h, const Eigen::Matrix3f &K)
    {
        int wlvl=w;
        int hlvl=h;
        int pyrLevelsUsed=1;
        while(wlvl%2==0 && hlvl%2==0 && wlvl*hlvl > 5000 && pyrLevelsUsed < PYR_LEVELS)
        {
            wlvl /=2;
            hlvl /=2;
            pyrLevelsUsed++;
        }
        Log::Logger::getInstance()->info("using pyramid levels 0 to {}. coarsest resolution: {} x {}!",
               pyrLevelsUsed-1, wlvl, hlvl);
        if(wlvl>100 && hlvl > 100)
        {
            Log::Logger::getInstance()->info("\n\n===============WARNING!===================\n "
                   "using not enough pyramid levels.\n"
                   "Consider scaling to a resolution that is a multiple of a power of 2.\n");
        }
        if(pyrLevelsUsed < 3)
        {
            Log::Logger::getInstance()->info("\\n===============WARNING!===================\\n \"\n"
                                             "\"I need higher resolution.\\n\"\n"
                                             "\"I will probably segfault.\\n");

            printf("\n");
        }

        wG[0] = w;
        hG[0] = h;
        KG[0] = K;
        fxG[0] = K(0,0);
        fyG[0] = K(1,1);
        cxG[0] = K(0,2);
        cyG[0] = K(1,2);
        KiG[0] = KG[0].inverse();
        fxiG[0] = KiG[0](0,0);
        fyiG[0] = KiG[0](1,1);
        cxiG[0] = KiG[0](0,2);
        cyiG[0] = KiG[0](1,2);

        for (int level = 1; level < pyrLevelsUsed; ++ level)
        {
            wG[level] = w >> level;
            hG[level] = h >> level;

            fxG[level] = fxG[level-1] * 0.5f;
            fyG[level] = fyG[level-1] * 0.5f;
            cxG[level] = (cxG[0] + 0.5f) / ((int)1<<level) - 0.5f;
            cyG[level] = (cyG[0] + 0.5f) / ((int)1<<level) - 0.5f;

            KG[level]  << fxG[level], 0.0, cxG[level], 0.0, fyG[level], cyG[level], 0.0, 0.0, 1.0;	// synthetic
            KiG[level] = KG[level].inverse();

            fxiG[level] = KiG[level](0,0);
            fyiG[level] = KiG[level](1,1);
            cxiG[level] = KiG[level](0,2);
            cyiG[level] = KiG[level](1,2);
        }
        return pyrLevelsUsed;

}

bool CameraCalibration::readIntrinsicCalibration(const std::filesystem::path& intrinsicsFilePath) {

    std::ifstream inFile, exFile;
    std::map<std::string, std::vector<float> > data;

    inFile.open (intrinsicsFilePath.c_str ());

    if (!inFile) {
        Log::Logger::getInstance()->error("Failed to open {} for reading", intrinsicsFilePath.c_str());
        return false;
    }

    parseYaml (inFile, data);

    inFile.close ();

    if (data["M1"].size () != 3 * 3 ||
        (data["D1"].size () != 5 && data["D1"].size () != 8) ||
        data["M2"].size () != 3 * 3 ||
        (data["D2"].size () != 5 && data["D2"].size () != 8)) {
        Log::Logger::getInstance()->error("intrinsic matrices incomplete in {}", intrinsicsFilePath.c_str());
        return false;
    }

    K = Eigen::Matrix3f::Zero();
    K(0, 0) = data["M1"][0];
    K(1, 1) = data["M1"][4];
    K(0, 2) = data["M1"][2];
    K(1, 2) = data["M1"][5];

    Log::Logger::getInstance()->info("Successfully read intrinsic calibration file: {}", intrinsicsFilePath.c_str());

    return true;
}

void CameraCalibration::setFramePyramidInfo(VO::Frame *frame) const{
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


