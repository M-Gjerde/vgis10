//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_ENERGY_H
#define VGIS10_ENERGY_H

#include "EFResidual.h"
#include "Frame.h"
#include "Util/IndexThreadReduce.h"
#include "AccumulatedTopHessian.h"
#include "AccumulatedSCHessian.h"
#include "CalibHessian.h"
#include "EFFrame.h"

extern bool EFAdjointsValid;
extern bool EFIndicesValid;
extern bool EFDeltaValid;

class EnergyFunctional {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    EnergyFunctional();
    ~EnergyFunctional();


    void insertResidual(PointFrameResidual *r, int hostFrameID, int targetFrameID);
    void insertFrame(std::shared_ptr<VO::Frame> fh, CalibHessian* Hcalib);
    void insertPoint(PointHessian *ph, uint32_t i);

    void dropResidual(PointFrameResidual* r);
    void marginalizeFrame(EFFrame* fh);
    void removePoint(EFPoint* ph);



    void marginalizePointsF();
    void dropPointsF();
    void solveSystemF(int iteration, double lambda, CalibHessian* HCalib);
    double calcMEnergyF();
    double calcLEnergyF_MT();


    void makeIDX();

    void setDeltaF(CalibHessian* HCalib);

    void setAdjointsF(CalibHessian* Hcalib);

    std::vector<EFFrame> eFrames;

    int nPoints, nFrames, nResiduals;

    MatXX HM;
    VecX bM;

    int resInA, resInL, resInM;
    MatXX lastHS;
    VecX lastbS;
    VecX lastX;
    std::vector<VecX> lastNullspaces_forLogging;
    std::vector<VecX> lastNullspaces_pose;
    std::vector<VecX> lastNullspaces_scale;
    std::vector<VecX> lastNullspaces_affA;
    std::vector<VecX> lastNullspaces_affB;

    dso::IndexThreadReduce<Vec10>* red;


    std::map<uint64_t,
    Eigen::Vector2i,
    std::less<uint64_t>,
    Eigen::aligned_allocator<std::pair<const uint64_t, Eigen::Vector2i>>
    > connectivityMap;

    Mat88* adHost;
    Mat88* adTarget;

    Mat88f* adHostF;
    Mat88f* adTargetF;
    Mat18f* adHTdeltaF;

    VecC cPrior;
    VecCf cDeltaF;
    VecCf cPriorF;
    std::vector<EFPoint*> allPoints;
private:

    VecX getStitchedDeltaF() const;

    void resubstituteF_MT(VecX x, CalibHessian* HCalib, bool MT);
    void resubstituteFPt(const VecCf &xc, Mat18f* xAd, int min, int max, Vec10* stats, int tid);

    void accumulateAF_MT(MatXX &H, VecX &b, bool MT);
    void accumulateLF_MT(MatXX &H, VecX &b, bool MT);
    void accumulateSCF_MT(MatXX &H, VecX &b, bool MT);

    void calcLEnergyPt(int min, int max, Vec10* stats, int tid);

    void orthogonalize(VecX* b, MatXX* H);

    dso::AccumulatedTopHessianSSE* accSSE_top_L;
    dso::AccumulatedTopHessianSSE* accSSE_top_A;


    dso::AccumulatedSCHessianSSE* accSSE_bot;

    std::vector<EFPoint*> allPointsToMarg;

    float currentLambda;
};

#endif //VGIS10_ENERGY_H
