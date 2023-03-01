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
        m_CamCal = std::make_unique<CameraCalibration>(
                "/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/vignette.png",
                "/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/pcalib.txt");
        m_FrameClass = std::make_unique<FrameClass>(
                "/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/images");
        //m_CamCal = std::make_unique<CameraCalibration>("/home/magnus/CLionProjects/vgis10/test_images/vignetteSmoothed.png", "/home/magnus/CLionProjects/vgis10/test_images/pcalib.txt");

        m_Tracker = std::make_unique<Tracker>(m_CamCal.get());
        m_IsRunning = true;
    }

    bool Core::spin() {

        std::shared_ptr<Frame> frame = m_FrameClass->getNextFrame(m_CamCal.get());
        frame->shell->marginalizedAt = m_Tracker->allFrameHistory.size();
        frame->shell->camToWorld = SE3(); 		// no lock required, as fh is not used anywhere yet.
        frame->shell->aff_g2l = AffLight(0,0);
        frame->shell->marginalizedAt = frame->shell->id = m_Tracker->allFrameHistory.size();
        frame->shell->incoming_id = m_FrameClass->m_FrameID;
        m_Tracker->allFrameHistory.push_back(frame->shell);
        if (!m_Tracker->initialized) {
            std::shared_ptr<Frame> firstFrame = frame;

            // Set initial frame
            if (!m_Tracker->firstFrame) {
                Log::Logger::getInstance()->info("Initializing with first frame");
                CandidatePointInfo info;
                m_Tracker->firstFrame = firstFrame;
                initializeCandidatePoints(m_Tracker->firstFrame, info, m_CamCal.get(), &m_Tracker->points);
            } else {
                if (m_Tracker->initializerTrackFrame(frame, m_CamCal.get())) {
                    m_Tracker->initializeFromInitializer();
                    m_Tracker->takeTrackedFrame(m_Tracker->secondFrame, true);
                    m_IsRunning = true;
                } else {
                    frame->shell->poseValid = false;

                }
            }


        } else {
            // Do front-end Operation
            // =========================== SWAP tracking reference?. =========================
            if (m_Tracker->coarseTracker_forNewKF->refFrameID > m_Tracker->coarseTracker->refFrameID) {
                Log::Logger::getInstance()->info("Swapping tracking reference");
                CoarseTracker *tmp = m_Tracker->coarseTracker;
                m_Tracker->coarseTracker = m_Tracker->coarseTracker_forNewKF;
                m_Tracker->coarseTracker_forNewKF = tmp;
            }


            Vec4 tres = m_Tracker->trackNewCoarse(frame);
            if (!std::isfinite((double) tres[0]) || !std::isfinite((double) tres[1]) ||
                !std::isfinite((double) tres[2]) || !std::isfinite((double) tres[3])) {
                printf("Initial Tracking failed: LOST!\n");
                return false;
            }

            bool needToMakeKF;

            Vec2 refToFh = AffLight::fromToVecExposure(m_Tracker->coarseTracker->lastRef->abExposure,
                                                       frame->abExposure,
                                                       m_Tracker->coarseTracker->lastRef_aff_g2l,
                                                       frame->aff_g2l());

            // BRIGHTNESS CHECK
            needToMakeKF = m_Tracker->allFrameHistory.size() == 1 ||
                           setting_kfGlobalWeight * setting_maxShiftWeightT * sqrtf((double) tres[1]) /
                           (m_Tracker->calibration->wG[0] + m_Tracker->calibration->hG[0]) +
                           setting_kfGlobalWeight * setting_maxShiftWeightR * sqrtf((double) tres[2]) /
                           (m_Tracker->calibration->wG[0] + m_Tracker->calibration->hG[0]) +
                           setting_kfGlobalWeight * setting_maxShiftWeightRT * sqrtf((double) tres[3]) /
                           (m_Tracker->calibration->wG[0] + m_Tracker->calibration->hG[0]) +
                           setting_kfGlobalWeight * setting_maxAffineWeight * fabs(logf((float) refToFh[0])) > 1 ||
                           2 * m_Tracker->coarseTracker->firstCoarseRMSE < tres[0];


            m_Tracker->takeTrackedFrame(frame, needToMakeKF);
        }

        cv::waitKey(0);
        Log::Logger::getInstance()->spinNumber++;

        return m_IsRunning;
    }

    void Core::cleanUp() {

    }


}