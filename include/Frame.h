//
// Created by magnus on 1/28/23.
//

#ifndef VGIS10_FRAME_H
#define VGIS10_FRAME_H

#include <iostream>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "Logger.h"
#include "Types.h"

#include "Settings.h"
#include "FramePose.h"

#include "PointHessian.h"

namespace VO {

    struct FrameToFramePrecalc
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        // precalc values
        Mat33f PRE_RTll;
        Mat33f PRE_KRKiTll;
        Mat33f PRE_RKiTll;
        Mat33f PRE_RTll_0;

        Vec2f PRE_aff_mode;
        float PRE_b0_mode{};

        Vec3f PRE_tTll;
        Vec3f PRE_KtTll;
        Vec3f PRE_tTll_0;

        float distanceLL{};


        inline ~FrameToFramePrecalc() = default;
        inline FrameToFramePrecalc() = default;
        void set(const VO::FramePose& host, const VO::FramePose& target, float fxl , float fyl, float cyl, float cxl){


            SE3 leftToLeft_0 = target.get_worldToCam_evalPT() * host.get_worldToCam_evalPT().inverse();
            PRE_RTll_0 = (leftToLeft_0.rotationMatrix()).cast<float>();
            PRE_tTll_0 = (leftToLeft_0.translation()).cast<float>();

            SE3 leftToLeft = target.PRE_worldToCam * host.PRE_camToWorld;
            PRE_RTll = (leftToLeft.rotationMatrix()).cast<float>();
            PRE_tTll = (leftToLeft.translation()).cast<float>();
            distanceLL = leftToLeft.translation().norm();

            Mat33f K = Mat33f::Zero();
            K(0,0) = fxl;
            K(1,1) = fyl;
            K(0,2) = cxl;
            K(1,2) = cyl;
            K(2,2) = 1;
            PRE_KRKiTll = K * PRE_RTll * K.inverse();
            PRE_RKiTll = PRE_RTll * K.inverse();
            PRE_KtTll = K * PRE_tTll;


            PRE_aff_mode = AffLight::fromToVecExposure(host.abExposure, target.abExposure, host.aff_g2l(), target.aff_g2l()).cast<float>();
            PRE_b0_mode = host.aff_g2l_0().b;
        }
    };

    struct Frame {
        Frame() {
            //Log::Logger::getInstance()->info("Created Frame Object");
        }

        ~Frame() {
            //Log::Logger::getInstance()->info("Destroyed Frame Object");
        }

        int width = 0, inWidth = 0;
        int height = 0, inHeight = 0;

        std::vector<unsigned char> pixels{};

        std::vector<float> dataRaw{}; // uncalibrated
        std::vector<float> dataf{};   // Calibrated
        std::vector<float> dataRawRectified{};   // Calibrated

        std::vector<Eigen::Vector3f> color{};
        std::vector<std::vector<float>> absSquaredGrad{};
        std::vector<std::vector<Eigen::Vector3f>> pyramidNoPhotometric{};
        std::vector<std::vector<Eigen::Vector3f>> pyramid{};

        std::vector<float> histogram{};
        std::vector<float> histogramSmoothed{};
        std::vector<int> gradientHistogram{};
        int histStep = 0;

        // Pose stuff
        FramePose pose;
        FramePose trackingRef;
        std::vector<FrameToFramePrecalc> targetPrecalc;

        // Points stuff
        std::vector<PointHessian> pointHessians;                // contains all ACTIVE points.
        std::vector<PointHessian> pointHessiansMarginalized;   	// contains all MARGINALIZED points (= fully marginalized, usually because point went OOB.)
        std::vector<PointHessian> pointHessiansOut;		        // contains all OUTLIER points (= discarded.).

        // Calibration stuff (global)
        int widthLevel[PYR_LEVELS]{}, heightLevel[PYR_LEVELS]{};
        float fxLevel[PYR_LEVELS]{}, fyLevel[PYR_LEVELS]{}, cxLevel[PYR_LEVELS]{}, cyLevel[PYR_LEVELS]{};
        float fxiLevel[PYR_LEVELS]{}, fyiLevel[PYR_LEVELS]{}, cxiLevel[PYR_LEVELS]{}, cyiLevel[PYR_LEVELS]{};
        Mat33 KLevel[PYR_LEVELS], KiLevel[PYR_LEVELS], K;

        float abExposure = 1; // exposure time in ms. // TODO set in image reader class
        float timestamp = 0;
        uint32_t id = 0;            // Overall frame number into the system
        uint32_t trackingID = 0;    // Number of frames tracked

        bool flaggedForMarginalization = false;
        /**@brief  Quick display functil for this frame. Displays whatever is in the dataf */
        void display(bool blocking = true, const std::string &windowName = "display") {
            cv::Mat img = cv::Mat(height, width, CV_32F, dataf.data());
            img.convertTo(img, CV_8UC1);
            cv::imshow(windowName, img);

            if (blocking)
                cv::waitKey(0);
            else
                cv::waitKey(30);
        }

        void displayRaw(bool blocking = true) {
            cv::Mat img(inHeight, inWidth, CV_32F, dataRaw.data());
            img.convertTo(img, CV_8UC1);
            cv::imshow("displayRaw", img);

            if (blocking)
                cv::waitKey(0);
            else
                cv::waitKey(30);
        }

        void displayPyramid(size_t lvl, bool blocking = true, const std::string &windowName = "display") {
            cv::Mat img = cv::Mat(heightLevel[lvl], widthLevel[lvl], CV_32FC3, pyramid[lvl].data());
            img.convertTo(img, CV_8UC1);
            cv::imshow(windowName, img);

            if (blocking)
                cv::waitKey(0);
            else
                cv::waitKey(30);
        }

        void displayGradient(size_t lvl, bool blocking = true) {
            cv::Mat img = cv::Mat(heightLevel[lvl], widthLevel[lvl], CV_32FC1, absSquaredGrad[lvl].data());
            img.convertTo(img, CV_8UC1);
            cv::imshow("display", img);

            if (blocking)
                cv::waitKey(0);
            else
                cv::waitKey(30);
        }
    };


    static Eigen::Vector3f
    getInterpolatedElement33(const Eigen::Vector3f *const mat, const float x, const float y, const int width) {
        int ix = (int) x;
        int iy = (int) y;
        float dx = x - ix;
        float dy = y - iy;
        float dxdy = dx * dy;
        const Eigen::Vector3f *bp = mat + ix + iy * width;

        return dxdy * *(const Eigen::Vector3f *) (bp + 1 + width)
               + (dy - dxdy) * *(const Eigen::Vector3f *) (bp + width)
               + (dx - dxdy) * *(const Eigen::Vector3f *) (bp + 1)
               + (1 - dx - dy + dxdy) * *(const Eigen::Vector3f *) (bp);
    }


    static float
    getInterpolatedElement31(const Eigen::Vector3f *const mat, const float x, const float y, const int width) {
        int ix = (int) x;
        int iy = (int) y;
        float dx = x - ix;
        float dy = y - iy;
        float dxdy = dx * dy;
        const Eigen::Vector3f *bp = mat + ix + iy * width;


        return dxdy * (*(const Eigen::Vector3f *) (bp + 1 + width))[0]
               + (dy - dxdy) * (*(const Eigen::Vector3f *) (bp + width))[0]
               + (dx - dxdy) * (*(const Eigen::Vector3f *) (bp + 1))[0]
               + (1 - dx - dy + dxdy) * (*(const Eigen::Vector3f *) (bp))[0];
    }


    static Eigen::Vector3f
    getInterpolatedElement33BiLin(const Eigen::Vector3f *const mat, const float x, const float y, const int width) {
        int ix = (int) x;
        int iy = (int) y;
        const Eigen::Vector3f *bp = mat + ix + iy * width;

        float tl = (*(bp))[0];
        float tr = (*(bp + 1))[0];
        float bl = (*(bp + width))[0];
        float br = (*(bp + width + 1))[0];

        float dx = x - ix;
        float dy = y - iy;
        float topInt = dx * tr + (1 - dx) * tl;
        float botInt = dx * br + (1 - dx) * bl;
        float leftInt = dy * bl + (1 - dy) * tl;
        float rightInt = dy * br + (1 - dy) * tr;

        return Eigen::Vector3f(
                dx * rightInt + (1 - dx) * leftInt,
                rightInt - leftInt,
                botInt - topInt);
    }
};


#endif //VGIS10_FRAME_H
