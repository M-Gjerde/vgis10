//
// Created by magnus on 3/1/23.
//

#include "Frame.h"

void VO::Frame::setStateZero(const Vec10 &state_zero) {

    assert(state_zero.head<6>().squaredNorm() < 1e-20);

    this->state_zero = state_zero;


    for(int i=0;i<6;i++)
    {
        Vec6 eps; eps.setZero(); eps[i] = 1e-3;
        SE3 EepsP = Sophus::SE3d::exp(eps);
        SE3 EepsM = Sophus::SE3d::exp(-eps);
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
    assert(ab_exposure > 0);
    nullspaces_affine.topRightCorner<2,1>() = Vec2(0, expf(aff_g2l_0().a)*abExposure);
}