//
// Created by magnus on 1/28/23.
//

#include <opencv2/opencv.hpp>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION

#include <stb_image.h>

#include "FrameClass.h"
#include "Logger.h"

std::shared_ptr<VO::Frame> VO::FrameClass::getNextFrame(const CameraCalibration *calibration) {
    m_PreviousFrame = m_CurrentFrame;
    if (m_FrameID >= m_datasetImagesCount) {
        return nullptr;
    }
    // read raw image data
    std::shared_ptr<VO::Frame> frame = readFrame();
    calibration->undistort(frame, frame->abExposure, 0, 1);
    // Make image pyramids
    VO::FrameClass::makeImagePyramid(frame.get(), calibration);
    m_FrameID++;
    return frame;
}

std::shared_ptr<VO::Frame> VO::FrameClass::readFrame() {
    std::shared_ptr<Frame> frame = std::make_shared<Frame>();
    // if from dataset
    std::filesystem::path filePath = m_FileNames.at(m_FrameID);
    int width = 0, height = 0, channels = 0;
    unsigned char *pixels = stbi_load(filePath.c_str(), &width, &height, &channels, 0);
    if (!pixels)
        Log::Logger::getInstance()->error("Failed to load image {}", filePath.string());


    frame->pixels.resize(width * height * channels);
    frame->pixels.insert(frame->pixels.begin(), &pixels[0], &pixels[width * height * channels]);
    // if from camera
    //todo optimize
    for (auto pix: frame->pixels) {
        frame->dataRaw.emplace_back(pix);
    }
    frame->inHeight = height;
    frame->inWidth = width;
    frame->id = m_FrameID;
    stbi_image_free(pixels);


    std::string line;
    std::getline(infile, line);  // this does the checking!
    std::istringstream iss(line);
        float value;
        int skip = 0;
        while (iss >> value)
        {
            if (skip != 2)
                skip++;
            else{
                frame->abExposure = value;
                break;
            }
        }
    Log::Logger::getInstance()->info("Loaded image {} into memory with exposure {}", filePath.string(), frame->abExposure);
    return frame;
}

void VO::FrameClass::makeImagePyramid(Frame *frame, const CameraCalibration *calibration) {
    calibration->setFramePyramidInfo(frame);
    frame->pyramid.resize(calibration->pyrLevelsUsed);
    frame->pyramidNoPhotometric.resize(calibration->pyrLevelsUsed);
    frame->absSquaredGrad.resize(calibration->pyrLevelsUsed);
    // populate first level

    for (auto pix: frame->dataf) {
        frame->pyramid[0].emplace_back(pix, pix, pix);
        frame->color.emplace_back(pix, pix, pix);
    }

    for (int lvl = 0; lvl < calibration->pyrLevelsUsed; ++lvl) {
        int width = frame->widthLevel[lvl], height = frame->heightLevel[lvl];
        frame->absSquaredGrad[lvl].resize(height * width);

        if (lvl > 0) {
            frame->pyramid[lvl].resize(height * width);
            int wlm1 = frame->widthLevel[lvl - 1];
            Eigen::Vector3f *dI_lm = frame->pyramid[lvl - 1].data();
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {

                    frame->pyramid[lvl][x + y * width][0] = 0.25f * (dI_lm[2 * x + 2 * y * wlm1][0] +
                                                                                  dI_lm[2 * x + 1 + 2 * y * wlm1][0] +
                                                                                  dI_lm[2 * x + 2 * y * wlm1 + wlm1][0] +
                                                                                  dI_lm[2 * x + 1 + 2 * y * wlm1 + wlm1][0]);
                }
            }
        }

        for (int idx = width; idx < (width * (height - 1)); idx++) {
            float dx = 0.5f * (frame->pyramid[lvl][idx + 1][0] - frame->pyramid[lvl][idx - 1][0]);
            float dy = 0.5f * (frame->pyramid[lvl][idx + width][0] - frame->pyramid[lvl][idx - width][0]);
            if (!std::isfinite(dx)) dx = 0;
            if (!std::isfinite(dy)) dy = 0;
            frame->pyramid[lvl][idx][1] = dx;
            frame->pyramid[lvl][idx][2] = dy;

            frame->absSquaredGrad[lvl][idx] = (dx * dx) +( dy * dy);
            // With gamma calibration
            float color = (frame->pyramid[lvl][idx][0]);
            int c = color + 0.5f;
            if (c < 5) c = 5;
            if (c > 250) c = 250;
            float gw = calibration->responseFunc.gInv[c + 1] - calibration->responseFunc.gInv[c];
            frame->absSquaredGrad[lvl][idx] *= gw * gw;    // convert to gradient of original color space (before removing response).

        }
    }
    /*
    for (int i = 0; i < frame->width * frame->height; ++i){
        if (i % 1000){
            Log::Logger::getInstance()->info("{} {} {}",  frame->color[i].x(), frame->color[i].y(), frame->color[i].z());
            Log::Logger::getInstance()->info("{} {} {}",  frame->pyramid[0][i].x(), frame->pyramid[0][i].y(), frame->pyramid[0][i].z());
        }
    }
*/

    Log::Logger::getInstance()->info("Finished creating image pyramids");
}