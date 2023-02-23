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
    explicit Tracker(const CameraCalibration * calib) :  coarseDistanceMap(calib->wG[0], calib->hG[0], calib), hCalib(calib->fxG[0],
    calib->fyG[0],
    calib->cxG[0],
    calib->cyG[0]){

        info.thisToNext = SE3();
        info.thisToNext_aff = AffLight(0,0);
        calibration = calib;
        //energyFunctional = std::make_unique<EnergyFunctional>();
        immaturePoints.resize(2);
        VO::initCalibHessian(&hCalib, calib);

        ef.red = &this->treadReduce;

        calibLog = new std::ofstream();
        calibLog->open("logs/calibLog.txt", std::ios::trunc | std::ios::out);
        calibLog->precision(12);

        numsLog = new std::ofstream();
        numsLog->open("logs/numsLog.txt", std::ios::trunc | std::ios::out);
        numsLog->precision(10);

        system("rm -rf logs");
        system("mkdir logs");

        eigenAllLog = new std::ofstream();
        eigenAllLog->open("logs/eigenAllLog.txt", std::ios::trunc | std::ios::out);
        eigenAllLog->precision(10);

        eigenPLog = new std::ofstream();
        eigenPLog->open("logs/eigenPLog.txt", std::ios::trunc | std::ios::out);
        eigenPLog->precision(10);

        eigenALog = new std::ofstream();
        eigenALog->open("logs/eigenALog.txt", std::ios::trunc | std::ios::out);
        eigenALog->precision(10);

        DiagonalLog = new std::ofstream();
        DiagonalLog->open("logs/diagonal.txt", std::ios::trunc | std::ios::out);
        DiagonalLog->precision(10);

        variancesLog = new std::ofstream();
        variancesLog->open("logs/variancesLog.txt", std::ios::trunc | std::ios::out);
        variancesLog->precision(10);


        nullspacesLog = new std::ofstream();
        nullspacesLog->open("logs/nullspacesLog.txt", std::ios::trunc | std::ios::out);
        nullspacesLog->precision(10);
    }
    ~Tracker(){
        calibLog->close(); delete calibLog;
        numsLog->close(); delete numsLog;
        //errorsLog->close(); delete errorsLog;
        eigenAllLog->close(); delete eigenAllLog;
        eigenPLog->close(); delete eigenPLog;
        eigenALog->close(); delete eigenALog;
        DiagonalLog->close(); delete DiagonalLog;
        variancesLog->close(); delete variancesLog;
        nullspacesLog->close(); delete nullspacesLog;
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
    std::vector<PointFrameResidual*> activeResiduals{};

    std::vector<std::shared_ptr<VO::FramePose>> allKeyFramesHistory;

    std::vector<std::vector<ImmaturePoint>> immaturePoints;     // contains all OUTLIER points (= discarded.).
    EnergyFunctional ef;
    dso::IndexThreadReduce<Vec10> treadReduce;

    PntResult points;
    bool initialized = false;
    Info info;
    float currentMinActDist = 2;
    std::vector<float> allResVec;


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

    float optimize(int mnumOptIts);

    Vec3 linearizeAll(bool fixLinearization);

    void linearizeAll_Reductor(bool fixLinearization, std::vector<PointFrameResidual*> *toRemove, int min, int max,
                               Vec10 *stats, int tid);

    void setNewFrameEnergyTH();

    void loadSateBackup();

    void backupState(bool backupLastStep);

    bool doStepFromBackup(float stepfacC, float stepfacT, float stepfacR, float stepfacA, float stepfacD);

    void solveSystem(int iteration, double lambda);

    std::vector<VecX> getNullspaces(std::vector<VecX> &nullspaces_pose, std::vector<VecX> &nullspaces_scale,
                                    std::vector<VecX> &nullspaces_affA, std::vector<VecX> &nullspaces_affB);

    void printEigenValLine();

    std::ofstream* calibLog;
    std::ofstream* numsLog;
    std::ofstream* errorsLog;
    std::ofstream* eigenAllLog;
    std::ofstream* eigenPLog;
    std::ofstream* eigenALog;
    std::ofstream* DiagonalLog;
    std::ofstream* variancesLog;
    std::ofstream* nullspacesLog;

};


#endif //VGIS10_TRACKER_H
