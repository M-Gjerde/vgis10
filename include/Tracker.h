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
#include "Optimized/Energy.h"
#include "FrameFramePrecalc.h"
#include "HessianBlocks.h"

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

    std::vector<EnergyFunctional> ef;
    std::vector<std::vector<PointHessian>> pointsHessian;
    std::vector<FrameFramePrecalc,Eigen::aligned_allocator<FrameFramePrecalc>> targetPrecalc;

};

class Tracker {

public:
    explicit Tracker(size_t numPyrLevels) : points(numPyrLevels){

        info.thisToNext = SE3();
        info.thisToNext_aff = AffLight(0,0);
    }
    std::shared_ptr<VO::Frame> firstFrame;

    PntResult points;
    bool initialized = false;
    Info info;

private:

};


#endif //VGIS10_TRACKER_H
