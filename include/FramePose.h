#include "Types.h"


namespace VO {

    struct FramePose {
        // precalc values
        SE3 PRE_worldToCam;
        SE3 PRE_camToWorld;

        SE3 frameToWorld;
        SE3 frameToTrackingRef;
        AffLight affLightg2l{};
        uint32_t trackingID = 0;
        float abExposure = 1;
        Mat66 nullspaces_pose;
        Mat42 nullspaces_affine;
        Vec6 nullspaces_scale;
        // variable info.
        SE3 worldToCam_evalPT;
        Vec10 state_zero;
        Vec10 state_scaled;
        Vec10 state;    // [0-5: worldToCam-leftEps. 6-7: a,b]
        Vec10 step;

        Vec10 step_backup;
        Vec10 state_backup;

        int marginalizedAt;
        double movedByOpt;
        bool poseValid = true;

        EIGEN_STRONG_INLINE const SE3 &get_worldToCam_evalPT() const { return worldToCam_evalPT; }

        EIGEN_STRONG_INLINE const Vec10 &get_state_zero() const { return state_zero; }

        EIGEN_STRONG_INLINE const Vec10 &get_state() const { return state; }

        EIGEN_STRONG_INLINE const Vec10 &get_state_scaled() const { return state_scaled; }

        EIGEN_STRONG_INLINE const Vec10 get_state_minus_stateZero() const { return get_state() - get_state_zero(); }

        inline Vec6 w2c_leftEps() const { return get_state_scaled().head<6>(); }

        inline AffLight aff_g2l() const { return AffLight(get_state_scaled()[6], get_state_scaled()[7]); }

        inline AffLight aff_g2l_0() const {
            return AffLight(get_state_zero()[6] * SCALE_A, get_state_zero()[7] * SCALE_B);
        }

        void setStateZero(const Vec10 &state_zero) {
            assert(state_zero.head<6>().squaredNorm() < 1e-20);

            this->state_zero = state_zero;


            for (int i = 0; i < 6; i++) {
                Vec6 eps;
                eps.setZero();
                eps[i] = 1e-3;
                SE3 EepsP = Sophus::SE3<double>::exp(eps);
                SE3 EepsM = Sophus::SE3<double>::exp(-eps);
                SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT() * EepsP) * get_worldToCam_evalPT().inverse();
                SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT() * EepsM) * get_worldToCam_evalPT().inverse();
                nullspaces_pose.col(i) = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log()) / (2e-3);
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
            nullspaces_scale = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log()) / (2e-3);


            nullspaces_affine.setZero();
            nullspaces_affine.topLeftCorner<2, 1>() = Vec2(1, 0);
            assert(abExposure > 0);
            nullspaces_affine.topRightCorner<2, 1>() = Vec2(0, expf(aff_g2l_0().a) * abExposure);
        };

        inline void setState(const Vec10 &state) {

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

        inline void setStateScaled(const Vec10 &state_scaled) {

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

        inline void setEvalPT(const SE3 &worldToCam_evalPT, const Vec10 &state) {

            this->worldToCam_evalPT = worldToCam_evalPT;
            setState(state);
            setStateZero(state);
        };

        inline void setEvalPT_scaled(const SE3 &worldToCam_evalPT, const AffLight &aff_g2l) {
            Vec10 initial_state = Vec10::Zero();
            initial_state[6] = aff_g2l.a;
            initial_state[7] = aff_g2l.b;
            this->worldToCam_evalPT = worldToCam_evalPT;
            setStateScaled(initial_state);
            setStateZero(this->get_state());
        };

        Vec10 getPrior() const {
            Vec10 p = Vec10::Zero();
            if(trackingID==0)
            {
                p.head<3>() = Vec3::Constant(setting_initialTransPrior);
                p.segment<3>(3) = Vec3::Constant(setting_initialRotPrior);

                p[6] = setting_initialAffAPrior;
                p[7] = setting_initialAffBPrior;
            }
            else
            {
                    p[6] = setting_affineOptModeA;
                    p[7] = setting_affineOptModeB;
            }
            p[8] = setting_initialAffAPrior;
            p[9] = setting_initialAffBPrior;
            return p;
        }

        static inline Vec10 getPriorZero() {
            return Vec10::Zero();
        }

        std::string logFrameToWorldTranslation() {
            std::string str;
            str.resize(100);
            sprintf(str.data(), "x:%f, y:%f, z:%f", frameToWorld.translation().x(), frameToWorld.translation().y(),
                    frameToWorld.translation().z());
            return str;
        }
    };
}