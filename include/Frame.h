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

static int staticPattern[8][2] = {
        {0,  -2},
        {-1, -1},
        {1,  -1},
        {-2, 0},
        {0,  0},
        {2,  0},
        {-1, 1},
        {0,  2}};

#define PYR_LEVELS 5

// parameters controlling pixel selection
#define setting_minGradHistCut 0.5f
#define setting_minGradHistAdd 7
#define setting_gradDownweightPerLevel 0.75
#define setting_selectDirectionDistribution true

#define setting_desiredPointDensity  2000 // aimed total points in the active window.

#define setting_gradient_block_32 32
#define setting_outlierTH 12*12                    // higher -> less strict
#define setting_outlierTHSumComponent 50*50        // higher -> less strong gradient-based reweighting .
#define setting_huberTH 9


// Parameters controlling adaptive energy threshold computation
#define setting_overallEnergyTHWeight 1
#define setting_frameEnergyTHConstWeight 0.5
#define setting_frameEnergyTHN 0.7f
#define setting_frameEnergyTHFacMedian 1.5
#define setting_coarseCutoffTH 20

/* initial hessian values to fix unobservable dimensions / priors on affine lighting parameters.
 */
#define setting_idepthFixPrior 50*50
#define setting_idepthFixPriorMargFac 600*600
#define setting_initialRotPrior 1e11
#define setting_initialTransPrior 1e10
#define setting_initialAffBPrior 1e14
#define setting_initialAffAPrior 1e14
#define setting_initialCalibHessian 5e9

#define setting_affineOptModeA 1e12 //-1: fix. >=0: optimize (with prior, if > 0).
#define setting_affineOptModeB 1e8  //-1: fix. >=0: optimize (with prior, if > 0).


/* require some minimum number of residuals for a point to become valid */
#define   setting_minGoodActiveResForMarg 3
#define   setting_minGoodResForMarg 4

#define SOLVER_REMOVE_POSEPRIOR (int)32
#define SOLVER_FIX_LAMBDA (int)128
#define SOLVER_ORTHOGONALIZE_X_LATER (int)2048
#define setting_solverMode SOLVER_FIX_LAMBDA | SOLVER_ORTHOGONALIZE_X_LATER

#define patternPadding 2
#define patternNum 8
#define patternP staticPattern

namespace VO {
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

        int widthLevel[PYR_LEVELS]{}, heightLevel[PYR_LEVELS]{};
        float fxLevel[PYR_LEVELS]{}, fyLevel[PYR_LEVELS]{}, cxLevel[PYR_LEVELS]{}, cyLevel[PYR_LEVELS]{};
        float fxiLevel[PYR_LEVELS]{}, fyiLevel[PYR_LEVELS]{}, cxiLevel[PYR_LEVELS]{}, cyiLevel[PYR_LEVELS]{};
        Mat33 KLevel[PYR_LEVELS], KiLevel[PYR_LEVELS], K;

        float abExposure = 1; // exposure time in ms. // TODO set in image reader class
        float timestamp = 0;
        uint32_t id{};

        // precalc values
        SE3 PRE_worldToCam;
        SE3 PRE_camToWorld;
        float frameEnergyTH = 8*8*patternNum; // set dynamically depending on tracking residual

        Mat66 nullspaces_pose;
        Mat42 nullspaces_affine;
        Vec6 nullspaces_scale;
        // variable info.
        SE3 worldToCam_evalPT;
        Vec10 state_zero;
        Vec10 state_scaled;
        Vec10 state;	// [0-5: worldToCam-leftEps. 6-7: a,b]

        EIGEN_STRONG_INLINE const SE3 &get_worldToCam_evalPT() const {return worldToCam_evalPT;}
        EIGEN_STRONG_INLINE const Vec10 &get_state_zero() const {return state_zero;}
        EIGEN_STRONG_INLINE const Vec10 &get_state() const {return state;}
        EIGEN_STRONG_INLINE const Vec10 &get_state_scaled() const {return state_scaled;}
        EIGEN_STRONG_INLINE const Vec10 get_state_minus_stateZero() const {return get_state() - get_state_zero();}

        inline Vec6 w2c_leftEps() const {return get_state_scaled().head<6>();}
        inline AffLight aff_g2l() const {return AffLight(get_state_scaled()[6], get_state_scaled()[7]);}
        inline AffLight aff_g2l_0() const {return AffLight(get_state_zero()[6]*SCALE_A, get_state_zero()[7]*SCALE_B);}

        void setStateZero(const Vec10 &state_zero){
            assert(state_zero.head<6>().squaredNorm() < 1e-20);

            this->state_zero = state_zero;


            for(int i=0;i<6;i++)
            {
                Vec6 eps; eps.setZero(); eps[i] = 1e-3;
                SE3 EepsP = Sophus::SE3<double>::exp(eps);
                SE3 EepsM = Sophus::SE3<double>::exp(-eps);
                SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT() * EepsP) * get_worldToCam_evalPT().inverse();
                SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT() * EepsM) * get_worldToCam_evalPT().inverse();
                nullspaces_pose.col(i) = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);
            }
            //nullspaces_pose.topRows<3>() *= SCALE_XI_TRANS_INVERSE;
            //nullspaces_pose.bottomRows<3>() *= SCALE_XI_ROT_INVERSE;

            // scale change
            SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT());
            w2c_leftEps_P_x0.translation() *= 1.00001;
            w2c_leftEps_P_x0 = w2c_leftEps_P_x0 * get_worldToCam_evalPT().inverse();
            SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT());
            w2c_leftEps_M_x0.translation() /= 1.00001;
            w2c_leftEps_M_x0 = w2c_leftEps_M_x0 * get_worldToCam_evalPT().inverse();
            nullspaces_scale = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);


            nullspaces_affine.setZero();
            nullspaces_affine.topLeftCorner<2,1>()  = Vec2(1,0);
            assert(abExposure > 0);
            nullspaces_affine.topRightCorner<2,1>() = Vec2(0, expf(aff_g2l_0().a)*abExposure);
        };

        inline void setState(const Vec10 &state)
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
        inline void setStateScaled(const Vec10 &state_scaled)
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
        inline void setEvalPT(const SE3 &worldToCam_evalPT, const Vec10 &state)
        {

            this->worldToCam_evalPT = worldToCam_evalPT;
            setState(state);
            setStateZero(state);
        };

        inline void setEvalPT_scaled(const SE3 &worldToCam_evalPT, const AffLight &aff_g2l)
        {
            Vec10 initial_state = Vec10::Zero();
            initial_state[6] = aff_g2l.a;
            initial_state[7] = aff_g2l.b;
            this->worldToCam_evalPT = worldToCam_evalPT;
            setStateScaled(initial_state);
            setStateZero(this->get_state());
        };

        inline Vec10 getPrior()
        {
            Vec10 p =  Vec10::Zero();
            if(id==0)
            {
                p.head<3>() = Vec3::Constant(setting_initialTransPrior);
                p.segment<3>(3) = Vec3::Constant(setting_initialRotPrior);
                if(setting_solverMode & SOLVER_REMOVE_POSEPRIOR) p.head<6>().setZero();

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
        void display(bool blocking = true, const std::string& windowName = "display") {
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
        void displayPyramid(size_t lvl, bool blocking = true, const std::string& windowName = "display"){
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
        int ix = (int)x;
        int iy = (int)y;
        float dx = x - ix;
        float dy = y - iy;
        float dxdy = dx*dy;
        const Eigen::Vector3f* bp = mat +ix+iy*width;

        return dxdy * *(const Eigen::Vector3f*)(bp+1+width)
               + (dy-dxdy) * *(const Eigen::Vector3f*)(bp+width)
               + (dx-dxdy) * *(const Eigen::Vector3f*)(bp+1)
               + (1-dx-dy+dxdy) * *(const Eigen::Vector3f*)(bp);
    }


    static float
    getInterpolatedElement31(const Eigen::Vector3f *const mat, const float x, const float y, const int width) {
        int ix = (int)x;
        int iy = (int)y;
        float dx = x - ix;
        float dy = y - iy;
        float dxdy = dx*dy;
        const Eigen::Vector3f* bp = mat +ix+iy*width;


        return dxdy * (*(const Eigen::Vector3f*)(bp+1+width))[0]
               + (dy-dxdy) * (*(const Eigen::Vector3f*)(bp+width))[0]
               + (dx-dxdy) * (*(const Eigen::Vector3f*)(bp+1))[0]
               + (1-dx-dy+dxdy) * (*(const Eigen::Vector3f*)(bp))[0];
    }


    static Eigen::Vector3f getInterpolatedElement33BiLin(const Eigen::Vector3f* const mat, const float x, const float y, const int width)
    {
        int ix = (int)x;
        int iy = (int)y;
        const Eigen::Vector3f* bp = mat +ix+iy*width;

        float tl = (*(bp))[0];
        float tr = (*(bp+1))[0];
        float bl = (*(bp+width))[0];
        float br = (*(bp+width+1))[0];

        float dx = x - ix;
        float dy = y - iy;
        float topInt = dx * tr + (1-dx) * tl;
        float botInt = dx * br + (1-dx) * bl;
        float leftInt = dy * bl + (1-dy) * tl;
        float rightInt = dy * br + (1-dy) * tr;

        return Eigen::Vector3f(
                dx * rightInt + (1-dx) * leftInt,
                rightInt-leftInt,
                botInt-topInt);
    }
};


#endif //VGIS10_FRAME_H
