//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_FRAMETOFRAMEPRECALC_H
#define VGIS10_FRAMETOFRAMEPRECALC_H


#include "CalibHessian.h"

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
    void set(const VO::FramePose& host, const VO::FramePose& target, CalibHessian * hCalib){


        SE3 leftToLeft_0 = target.get_worldToCam_evalPT() * host.get_worldToCam_evalPT().inverse();
        PRE_RTll_0 = (leftToLeft_0.rotationMatrix()).cast<float>();
        PRE_tTll_0 = (leftToLeft_0.translation()).cast<float>();

        SE3 leftToLeft = target.PRE_worldToCam * host.PRE_camToWorld;
        PRE_RTll = (leftToLeft.rotationMatrix()).cast<float>();
        PRE_tTll = (leftToLeft.translation()).cast<float>();
        distanceLL = leftToLeft.translation().norm();

        Mat33f K = Mat33f::Zero();
        K(0,0) = hCalib->fxl();
        K(1,1) = hCalib->fyl();
        K(0,2) = hCalib->cxl();
        K(1,2) = hCalib->cyl();
        K(2,2) = 1;
        PRE_KRKiTll = K * PRE_RTll * K.inverse();
        PRE_RKiTll = PRE_RTll * K.inverse();
        PRE_KtTll = K * PRE_tTll;


        PRE_aff_mode = AffLight::fromToVecExposure(host.abExposure, target.abExposure, host.aff_g2l(), target.aff_g2l()).cast<float>();
        PRE_b0_mode = host.aff_g2l_0().b;
    }
};


#endif //VGIS10_FRAMETOFRAMEPRECALC_H
