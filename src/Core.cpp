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

        m_Tracker = std::make_unique<Tracker>(PYR_LEVELS);
        m_IsRunning = true;
    }

    bool Core::spin() {
        Log::Logger::getInstance()->frameNumber++;
        std::shared_ptr<Frame> frame = m_FrameClass->getNextFrame(m_CamCal.get());
        if (!frame) {
            m_IsRunning = false;
            return m_IsRunning;
        }

        if (!m_Tracker->initialized){
            // Set initial frame
            if (!m_Tracker->firstFrame) {
                CandidatePointInfo info;
                m_Tracker->firstFrame = frame;
                initializeCandidatePoints(m_Tracker->firstFrame, info, m_CamCal.get(), &m_Tracker->points);
            }
            // Start tracking w.r.t initial frame
            while (!initializerTrackFrame(m_FrameClass->getNextFrame(m_CamCal.get()), m_CamCal.get(), &m_Tracker->points, &m_Tracker->info,
                                         m_Tracker->firstFrame)){

                // Initialize from initializer
                initializeFromInitializer(frame);
                // Delivertracker frame

                // Make KF or make nonKF


                m_Tracker->initialized = true;
                int loop;

            }

            //m_Tracker->initialized = true;
        }

        /*
        m_CamCal->applyPhotometricCalibration(&frame);
        cv::Mat img(frame->height, frame->width, CV_32FC1, frame->data.data());
        cv::Mat original(frame->height, frame->width, CV_8UC1, frame->pixels.data());
        img.convertTo(img, CV_8UC1, 1);
        cv::imshow("Corrected", img);
        cv::imshow("Original", original);
        cv::waitKey(0);
         */


        return m_IsRunning;
    }

    void Core::cleanUp() {

    }


}