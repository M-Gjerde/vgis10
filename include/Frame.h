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
#include "CalibHessian.h"

class EFFrame;
class ImmaturePoint;

namespace VO {

    class FrameToFramePrecalc;
    struct Frame {
        Frame() {
            Log::Logger::getInstance()->info("Created Frame Object {}", static_cast<void*>(this));
            shell = new VO::FramePose();
        }

        ~Frame() {
            Log::Logger::getInstance()->info("Destroyed Frame Object: {}", static_cast<void*>(this));
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

        // Points stuff
        std::vector<PointHessian*> pointHessians;                // contains all ACTIVE points.
        std::vector<PointHessian*> pointHessiansMarginalized;   	// contains all MARGINALIZED points (= fully marginalized, usually because point went OOB.)
        std::vector<PointHessian*> pointHessiansOut;		        // contains all OUTLIER points (= discarded.).
        std::vector<ImmaturePoint*> immaturePoints;     // contains all OUTLIER points (= discarded.).

        // EnergyStuff
        float frameEnergyTH = 8*8*patternNum;
        EFFrame* efFrame;

        // Pose stuff
        FramePose* shell;

        // Calibration stuff (global)
        int widthLevel[PYR_LEVELS]{}, heightLevel[PYR_LEVELS]{};
        float fxLevel[PYR_LEVELS]{}, fyLevel[PYR_LEVELS]{}, cxLevel[PYR_LEVELS]{}, cyLevel[PYR_LEVELS]{};
        float fxiLevel[PYR_LEVELS]{}, fyiLevel[PYR_LEVELS]{}, cxiLevel[PYR_LEVELS]{}, cyiLevel[PYR_LEVELS]{};
        Mat33 KLevel[PYR_LEVELS], KiLevel[PYR_LEVELS], K;

        float abExposure = 1; // exposure time in ms. // TODO set in image reader class
        float timestamp = 0;
        uint32_t frameID = 0-1;           // incremental ID for keyframes only!
        uint32_t trackingID = 0;

        bool flaggedForMarginalization = false;


        // Tracking info
        Mat66 nullspaces_pose;
        Mat42 nullspaces_affine;
        Vec6 nullspaces_scale;

        // variable info.
        SE3 worldToCam_evalPT;
        Vec10 state_zero;
        Vec10 state_scaled;
        Vec10 state;	// [0-5: worldToCam-leftEps. 6-7: a,b]
        Vec10 step;
        Vec10 step_backup;
        Vec10 state_backup;
        // precalc values
        SE3 PRE_worldToCam;
        SE3 PRE_camToWorld;
        std::vector<FrameToFramePrecalc,Eigen::aligned_allocator<FrameToFramePrecalc>> targetPrecalc;

        EIGEN_STRONG_INLINE const SE3 &get_worldToCam_evalPT() const {return worldToCam_evalPT;}
        EIGEN_STRONG_INLINE const Vec10 &get_state_zero() const {return state_zero;}
        EIGEN_STRONG_INLINE const Vec10 &get_state() const {return state;}
        EIGEN_STRONG_INLINE const Vec10 &get_state_scaled() const {return state_scaled;}
        EIGEN_STRONG_INLINE const Vec10 get_state_minus_stateZero() const {return get_state() - get_state_zero();}

        inline Vec6 w2c_leftEps() const {return get_state_scaled().head<6>();}
        inline AffLight aff_g2l() const {return AffLight(get_state_scaled()[6], get_state_scaled()[7]);}
        inline AffLight aff_g2l_0() const {return AffLight(get_state_zero()[6]*SCALE_A, get_state_zero()[7]*SCALE_B);}

        void setStateZero(const Vec10 &state_zero);
        void setState(const Vec10 &state)
        {

            this->state = state;
            state_scaled.segment<3>(0) = SCALE_XI_TRANS * state.segment<3>(0);
            state_scaled.segment<3>(3) = SCALE_XI_ROT * state.segment<3>(3);
            state_scaled[6] = SCALE_A * state[6];
            state_scaled[7] = SCALE_B * state[7];
            state_scaled[8] = SCALE_A * state[8];
            state_scaled[9] = SCALE_B * state[9];

            PRE_worldToCam = SE3::exp(w2c_leftEps()) * get_worldToCam_evalPT();
            PRE_camToWorld = PRE_worldToCam.inverse();
            //setCurrentNullspace();
        };
        void setStateScaled(const Vec10 &state_scaled)
        {

            this->state_scaled = state_scaled;
            state.segment<3>(0) = SCALE_XI_TRANS_INVERSE * state_scaled.segment<3>(0);
            state.segment<3>(3) = SCALE_XI_ROT_INVERSE * state_scaled.segment<3>(3);
            state[6] = SCALE_A_INVERSE * state_scaled[6];
            state[7] = SCALE_B_INVERSE * state_scaled[7];
            state[8] = SCALE_A_INVERSE * state_scaled[8];
            state[9] = SCALE_B_INVERSE * state_scaled[9];

            PRE_worldToCam = SE3::exp(w2c_leftEps()) * get_worldToCam_evalPT();
            PRE_camToWorld = PRE_worldToCam.inverse();
            //setCurrentNullspace();
        };
        void setEvalPT(const SE3 &worldToCam_evalPT, const Vec10 &state)
        {

            this->worldToCam_evalPT = worldToCam_evalPT;
            setState(state);
            setStateZero(state);
        };



        void setEvalPT_scaled(const SE3 &worldToCam_evalPT, const AffLight &aff_g2l)
        {
            Vec10 initial_state = Vec10::Zero();
            initial_state[6] = aff_g2l.a;
            initial_state[7] = aff_g2l.b;
            this->worldToCam_evalPT = worldToCam_evalPT;
            setStateScaled(initial_state);
            setStateZero(this->get_state());
        };

        Vec10 getPrior()
        {
            Vec10 p =  Vec10::Zero();
            if(frameID==0)
            {
                p.head<3>() = Vec3::Constant(setting_initialTransPrior);
                p.segment<3>(3) = Vec3::Constant(setting_initialRotPrior);

                p[6] = setting_initialAffAPrior;
                p[7] = setting_initialAffBPrior;
            }
            else
            {
                if(setting_affineOptModeA < 0)
                    p[6] = setting_initialAffAPrior;
                else
                    p[6] = setting_affineOptModeA;

                if(setting_affineOptModeB < 0)
                    p[7] = setting_initialAffBPrior;
                else
                    p[7] = setting_affineOptModeB;
            }
            p[8] = setting_initialAffAPrior;
            p[9] = setting_initialAffBPrior;
            return p;
        }


        inline Vec10 getPriorZero()
        {
            return Vec10::Zero();
        }

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



    struct FrameToFramePrecalc {
        VO::Frame *host;    // defines row
        VO::Frame *target;    // defines column

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

        int hostTrackingID = 0;
        int targetTrackingID = 0;

        ~FrameToFramePrecalc() = default;

        FrameToFramePrecalc() = default;

        void set(VO::Frame *host, VO::Frame *target, CalibHessian *HCalib) {

            this->host = host;
            this->target = target;

            SE3 leftToLeft_0 = target->get_worldToCam_evalPT() * host->get_worldToCam_evalPT().inverse();
            PRE_RTll_0 = (leftToLeft_0.rotationMatrix()).cast<float>();
            PRE_tTll_0 = (leftToLeft_0.translation()).cast<float>();

            SE3 leftToLeft = target->PRE_worldToCam * host->PRE_camToWorld;
            PRE_RTll = (leftToLeft.rotationMatrix()).cast<float>();
            PRE_tTll = (leftToLeft.translation()).cast<float>();
            distanceLL = leftToLeft.translation().norm();

            Mat33f K = Mat33f::Zero();
            K(0, 0) = HCalib->fxl();
            K(1, 1) = HCalib->fyl();
            K(0, 2) = HCalib->cxl();
            K(1, 2) = HCalib->cyl();
            K(2, 2) = 1;
            PRE_KRKiTll = K * PRE_RTll * K.inverse();
            PRE_RKiTll = PRE_RTll * K.inverse();
            PRE_KtTll = K * PRE_tTll;

            PRE_aff_mode = AffLight::fromToVecExposure(host->abExposure, target->abExposure, host->aff_g2l(), target->aff_g2l()).cast<float>();
            PRE_b0_mode = host->aff_g2l_0().b;
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
