//
// Created by magnus on 2/27/23.
//

#ifndef VGIS10_COARSETRACKER_H
#define VGIS10_COARSETRACKER_H

#include "Types.h"
#include "Frame.h"
#include "CalibHessian.h"
#include "Optimized/MatrixAccumulators.h"
#include "CameraCalibration.h"
#include "Optimized/EFPoint.h"

class CoarseTracker {

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
    CoarseTracker(int w, int h,  const CameraCalibration* calib);
    ~CoarseTracker();

    const CameraCalibration* calibration;

    bool trackNewestCoarse(
            std::shared_ptr<VO::Frame> newFrameHessian,
            SE3 &lastToNew_out, AffLight &aff_g2l_out,
            int coarsestLvl, Vec5 minResForAbort);

    void setCoarseTrackingRef(
            std::vector<std::shared_ptr<VO::Frame>> frameHessians);

    void makeK(
            CalibHessian* HCalib, const CameraCalibration *calib);

    bool debugPrint, debugPlot;

    Mat33f K[PYR_LEVELS];
    Mat33f Ki[PYR_LEVELS];
    float fx[PYR_LEVELS];
    float fy[PYR_LEVELS];
    float fxi[PYR_LEVELS];
    float fyi[PYR_LEVELS];
    float cx[PYR_LEVELS];
    float cy[PYR_LEVELS];
    float cxi[PYR_LEVELS];
    float cyi[PYR_LEVELS];
    int w[PYR_LEVELS];
    int h[PYR_LEVELS];

    std::shared_ptr<VO::Frame> lastRef;
    AffLight lastRef_aff_g2l;
    std::shared_ptr<VO::Frame> newFrame;
    int refFrameID;

    // act as pure ouptut
    Vec5 lastResiduals;
    Vec3 lastFlowIndicators;
    double firstCoarseRMSE;

    void makeCoarseDepthL0(std::vector<std::shared_ptr<VO::Frame>> frameHessians,  const CameraCalibration* calibration);
    float* idepth[PYR_LEVELS];
    float* weightSums[PYR_LEVELS];
    float* weightSums_bak[PYR_LEVELS];


    Vec6 calcResAndGS(int lvl, Mat88 &H_out, Vec8 &b_out, const SE3 &refToNew, AffLight aff_g2l, float cutoffTH);
    Vec6 calcRes(int lvl, const SE3 &refToNew, AffLight aff_g2l, float cutoffTH);
    void calcGSSSE(int lvl, Mat88 &H_out, Vec8 &b_out, const SE3 &refToNew, AffLight aff_g2l);
    void calcGS(int lvl, Mat88 &H_out, Vec8 &b_out, const SE3 &refToNew, AffLight aff_g2l);

    // pc buffers
    float* pc_u[PYR_LEVELS];
    float* pc_v[PYR_LEVELS];
    float* pc_idepth[PYR_LEVELS];
    float* pc_color[PYR_LEVELS];
    int pc_n[PYR_LEVELS];

    // warped buffers
    float* buf_warped_idepth;
    float* buf_warped_u;
    float* buf_warped_v;
    float* buf_warped_dx;
    float* buf_warped_dy;
    float* buf_warped_residual;
    float* buf_warped_weight;
    float* buf_warped_refColor;
    int buf_warped_n;


    std::vector<float*> ptrToDelete;


    dso::Accumulator9 acc;


};


#endif //VGIS10_COARSETRACKER_H
