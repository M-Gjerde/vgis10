//
// Created by magnus on 1/28/23.
//

#include <iostream>
#include "FrameClass.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Logger.h"

std::shared_ptr<VO::Frame> VO::FrameClass::getNextFrame() {
    m_PreviousFrame = m_CurrentFrame;
    if (m_FrameID >= m_datasetImagesCount){
        return nullptr;
    }
    std::shared_ptr<VO::Frame> frame = readFrame();
    std::cout << "Use count: " << frame.use_count() << "\n";
    return frame;
}

std::shared_ptr<VO::Frame> VO::FrameClass::readFrame() {
    std::shared_ptr<Frame> frame = std::make_shared<Frame>();
    // if from dataset
    std::filesystem::path filePath = m_FileNames.at(m_FrameID);
    int width = 0, height = 0, channels = 0;
    unsigned char *img = stbi_load(filePath.c_str(), &width, &height, &channels, 0);
    if (!img)
        Log::Logger::getInstance()->error("Failed to load image {}", filePath.string());

    // if from camera
    frame->height = 1920;
    frame->width = 1080;
    frame->pixels = 0;
    return frame;
}
