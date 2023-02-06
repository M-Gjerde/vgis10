//
// Created by magnus on 1/28/23.
//

#ifndef VGIS10_FRAMECLASS_H
#define VGIS10_FRAMECLASS_H

#include <filesystem>
#include <utility>
#include <vector>
#include "Frame.h"
#include "CameraCalibration.h"

namespace VO {
    class FrameClass {
    public:
        explicit FrameClass(std::filesystem::path pathToFolder) : m_FolderPath(std::move(pathToFolder)) {
            // Create dictionary of files using a vector. If from dataset sorted in its timeseries
            std::vector<std::filesystem::path> filesInDirectory;
            std::copy(std::filesystem::directory_iterator(m_FolderPath), std::filesystem::directory_iterator(),
                      std::back_inserter(filesInDirectory));
            std::sort(filesInDirectory.begin(), filesInDirectory.end());
            m_FileNames = filesInDirectory;
            m_datasetImagesCount = m_FileNames.size();
        }

        void setDataSource() {

        }

        std::shared_ptr<Frame> getNextFrame(const CameraCalibration* calibration);

    private:
        std::filesystem::path m_FolderPath{};
        std::vector<std::filesystem::path> m_FileNames{};
        uint32_t m_datasetImagesCount = 0;

        uint32_t m_FrameID = 35;
        Frame *m_CurrentFrame{};
        Frame *m_PreviousFrame{};


        std::shared_ptr<VO::Frame> readFrame();

        static void makeImagePyramid(Frame *frame, const CameraCalibration* calibration);
    };
};

#endif //VGIS10_FRAMECLASS_H
