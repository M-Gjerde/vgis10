//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_FRAMEFRAMEPRECALC_H
#define VGIS10_FRAMEFRAMEPRECALC_H


struct FrameFramePrecalc
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    // static values
    static int instanceCounter;
    std::shared_ptr<VO::Frame> host;	// defines row
    std::shared_ptr<VO::Frame> target;	// defines column

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


    inline ~FrameFramePrecalc() = default;
    inline FrameFramePrecalc() {host=target=0;}
    void set(std::shared_ptr<VO::Frame> host, std::shared_ptr<VO::Frame> target, const CameraCalibration* HCalib);
};


#endif //VGIS10_FRAMEFRAMEPRECALC_H
