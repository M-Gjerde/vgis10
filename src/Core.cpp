//
// Created by magnus on 1/28/23.
//

#include <iostream>
#include <opencv2/opencv.hpp>
#include "Core.h"
#include "PointSelection.h"
#include "TrackerInitializer.h"

namespace VO {
    Core::Core() {
        //m_FrameClass = std::make_unique<FrameClass>("/home/magnus/CLionProjects/vgis10/test_images/images");
        m_CamCal = std::make_unique<CameraCalibration>("/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/vignette.png", "/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/pcalib.txt");
        m_FrameClass = std::make_unique<FrameClass>("/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/images");
        //m_CamCal = std::make_unique<CameraCalibration>("/home/magnus/CLionProjects/vgis10/test_images/vignetteSmoothed.png", "/home/magnus/CLionProjects/vgis10/test_images/pcalib.txt");

        m_Tracker = std::make_unique<Tracker>(m_CamCal.get());
        m_IsRunning = true;
    }

    bool Core::spin() {

        if (!m_Tracker->initialized){
            std::shared_ptr<Frame> firstFrame = m_FrameClass->getNextFrame(m_CamCal.get());

            // Set initial frame
            if (!m_Tracker->firstFrame) {
                CandidatePointInfo info;
                m_Tracker->firstFrame = firstFrame;
                initializeCandidatePoints(m_Tracker->firstFrame, info, m_CamCal.get(), &m_Tracker->points);
            }
            // Start tracking w.r.t initial frame
            while (!m_Tracker->initializerTrackFrame(m_FrameClass->getNextFrame(m_CamCal.get()), m_CamCal.get())){
            }

            m_Tracker->initializeFromInitializer();
            m_IsRunning = false;
            m_Tracker->takeTrackedFrame(m_Tracker->secondFrame, true);
        }

        Log::Logger::getInstance()->spinNumber++;

        return m_IsRunning;
    }

    void Core::cleanUp() {

    }


}