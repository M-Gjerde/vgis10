//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_IMMATUREPOINT_H
#define VGIS10_IMMATUREPOINT_H

#include "Frame.h"
#include "ImmaturePointTemporaryResidual.h"
#include "CalibHessian.h"
#include "Util/Enums.h"
#include "Util/Linearize.h"

class CameraCalibration;

struct ImmaturePoint {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    // static values
    ImmaturePoint(float uCoord, float vCoord, float t, const std::shared_ptr<VO::Frame> &frame);

    ~ImmaturePoint() = default;

    float color[MAX_RES_PER_POINT]{};
    float weights[MAX_RES_PER_POINT]{};
    VO::Frame* host;
    Mat22f gradH;
    Vec2f gradH_ev;
    Mat22f gradH_eig;
    int trackingID = 0;
    float energyTH{};
    float u, v;
    int idxInImmaturePoints{};
    float quality{};

    float type;
    bool remove = false;
    float my_type;

    float idepth_min{};
    float idepth_max{};

    ImmaturePointStatus lastTraceStatus;
    Vec2f lastTraceUV;
    float lastTracePixelInterval{};


    float idepth_GT;


    float getdPixdd(
            CalibHessian *HCalib,
            ImmaturePointTemporaryResidual *tmpRes,
            float idepth) {
        return 0;

    }

    float calcResidual(
            CalibHessian *HCalib, const float outlierTHSlack,
            ImmaturePointTemporaryResidual *tmpRes,
            float idepth) {

        return 0;
    }


    ImmaturePointStatus
    traceOn(std::shared_ptr<VO::Frame> frame, const Mat33f &hostToFrame_KRKi, const Vec3f &hostToFrame_Kt,
            const Vec2f &hostToFrame_affine, CalibHessian *HCalib, const CameraCalibration *calibration,
            bool debugPrint);

    double
    linearizeResidual(CalibHessian *HCalib, const float outlierTHSlack, ImmaturePointTemporaryResidual *tmpRes,
                      float &Hdd,
                      float &bd, float idepth, const VO::Frame *target);
};


#endif //VGIS10_IMMATUREPOINT_H
