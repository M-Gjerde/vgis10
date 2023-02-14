//
// Created by magnus on 2/7/23.
//

#ifndef VGIS10_TRACKER_H
#define VGIS10_TRACKER_H


#include <cstdlib>
#include "Point.h"

#include "sophus/geometry.hpp"
#include "Types.h"
#include "Frame.h"
#include "Optimized/MatrixAccumulators.h"
#include "CameraCalibration.h"
#include "FrameToFramePrecalc.h"
#include "CalibHessian.h"
#include "PointHessian.h"
#include "Optimized/PointFrameResidual.h"
#include "Optimized/Energy.h"
#include "CoarseDistanceMap.h"
#include "ImmaturePointTemporaryResidual.h"
#include "Util/Populate.h"

struct Info {
    SE3 thisToNext;
    AffLight thisToNext_aff;
    Eigen::DiagonalMatrix<float, 8> wM;

    float alphaK = 2.5*2.5;//*freeDebugParam1*freeDebugParam1;
    float alphaW = 150*150;//*freeDebugParam2*freeDebugParam2;
    float regWeight = 0.8;//*freeDebugParam4;
    float couplingWeight = 1;//*freeDebugParam5;


    bool snapped = false;
    int snappedAtFrame = false;
    int frameID = 0;

    // temporary buffers for H and b.
    std::vector<Vec10f> JbBuffer;			// 0-7: sum(dd * dp). 8: sum(res*dd). 9: 1/(1+sum(dd*dd))=inverse hessian entry.
    std::vector<Vec10f> JbBuffer_new;

    dso::Accumulator9 acc9;
    dso::Accumulator9 acc9SC;

    //std::vector<EnergyFunctional> ef;
    //std::vector<std::vector<PointHessian>> pointsHessian;
//    std::vector<FrameFramePrecalc,Eigen::aligned_allocator<FrameFramePrecalc>> targetPrecalc;

};


class Tracker {

public:
    explicit Tracker(const CameraCalibration * calib) :  coarseDistanceMap(calib->wG[0], calib->hG[0], calib){

        info.thisToNext = SE3();
        info.thisToNext_aff = AffLight(0,0);
        calibration = calib;
        //energyFunctional = std::make_unique<EnergyFunctional>();
        immaturePoints.resize(2);
        VO::initCalibHessian(&hCalib, calib);
    }
    // Keep first two frames in memory
    std::shared_ptr<VO::Frame> firstFrame;
    std::shared_ptr<VO::Frame> secondFrame; // Second frame after initialization


    CalibHessian hCalib;
    const CameraCalibration* calibration;
    CoarseDistanceMap coarseDistanceMap;
    // All tracking math
    //std::unique_ptr<EnergyFunctional> energyFunctional;
    std::vector<std::shared_ptr<VO::Frame>> frameHessians;
    std::vector<FrameToFramePrecalc> frameToFramePreCalc;

    std::vector<std::vector<ImmaturePoint>> immaturePoints;     // contains all OUTLIER points (= discarded.).

    PntResult points;
    bool initialized = false;
    Info info;
    float currentMinActDist = 2;

    EnergyFunctional ef;

    bool initializerTrackFrame(const std::shared_ptr<VO::Frame> &frame, const CameraCalibration *calibration);

    void initializeFromInitializer();
    void takeTrackedFrame(std::shared_ptr<VO::Frame> frame, bool needKF);
private:

    void makeKeyFrame(std::shared_ptr<VO::Frame> frame);
    void makeNonKeyFrame(std::shared_ptr<VO::Frame> frame);


    void traceNewCoarse(std::shared_ptr<VO::Frame> frame);

    void flagFramesForMarginalization(std::shared_ptr<VO::Frame> sharedPtr);

    void activatePointsMT();

    void setPrecalcValues();

    void
    activatePointsMT_Reductor(std::vector<PointHessian> *optimized, std::vector<ImmaturePoint> *toOptimize, int min,
                              int max, Vec10 *stats, int tid);

    void optimizeImmaturePoint(ImmaturePoint *point, int minObs, ImmaturePointTemporaryResidual *residuals,
                               PointHessian *phOut);
};


#endif //VGIS10_TRACKER_H
