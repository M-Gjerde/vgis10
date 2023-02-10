//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_ENERGY_H
#define VGIS10_ENERGY_H


#include "Types.h"
#include "Util/IndexThreadReduce.h"
#include "Residuals.h"
#include "AccumulatedTopHessian.h"
#include "AccumulatedSCHessian.h"

enum EFPointStatus {PS_GOOD=0, PS_MARGINALIZE, PS_DROP};

class EFPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFPoint(PointHessian* d) : data(d)
    {
        takeData();
        stateFlag=EFPointStatus::PS_GOOD;
    }
    void takeData();

    PointHessian* data;



    float priorF;
    float deltaF;


    // constant info (never changes in-between).
    int idxInPoints;

    float bdSumF;
    float HdiF;
    float Hdd_accLF;
    VecCf Hcd_accLF;
    float bd_accLF;
    float Hdd_accAF;
    VecCf Hcd_accAF;
    float bd_accAF;


    EFPointStatus stateFlag;
};

class EFFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFFrame(std::shared_ptr<VO::Frame> d) : data(d)
    {
        takeData();
    }
    void takeData();


    Vec8 prior;				// prior hessian (diagonal)
    Vec8 delta_prior;		// = state-state_prior (E_prior = (delta_prior)' * diag(prior) * (delta_prior)
    Vec8 delta;				// state - state_zero.



    std::vector<EFPoint*> points;
    std::shared_ptr<VO::Frame> data;
    int idx;	// idx in frames.

    int frameID;
};


class EFResidual
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    inline EFResidual(PointFrameResidual* org, EFPoint* point_, EFFrame* host_, EFFrame* target_) :
            data(org), point(point_), host(host_), target(target_)
    {
        isLinearized=false;
        isActiveAndIsGoodNEW=false;
        J = new RawResidualJacobian();
        assert(((long)this)%16==0);
        assert(((long)J)%16==0);
    }
    inline ~EFResidual()
    {
        delete J;
    }


    void takeDataF();

    // structural pointers
    PointFrameResidual* data;
    int hostIDX, targetIDX;
    EFPoint* point;
    EFFrame* host;
    EFFrame* target;
    int idxInAll;

    RawResidualJacobian* J;

    VecNRf res_toZeroF;
    Vec8f JpJdF;


    // status.
    bool isLinearized;

    // if residual is not OOB & not OUTLIER & should be used during accumulations
    bool isActiveAndIsGoodNEW;
    inline const bool &isActive() const {return isActiveAndIsGoodNEW;}
};


class EnergyFunctional {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    EnergyFunctional();
    ~EnergyFunctional();


    EFResidual* insertResidual(PointFrameResidual* r);
    EFFrame* insertFrame(std::shared_ptr<VO::Frame> fh, CalibHessian* Hcalib);
    EFPoint* insertPoint(PointHessian* ph);

    void dropResidual(EFResidual* r);
    void marginalizeFrame(EFFrame* fh);
    void removePoint(EFPoint* ph);



    void marginalizePointsF();
    void dropPointsF();
    void solveSystemF(int iteration, double lambda, CalibHessian* HCalib );
    double calcMEnergyF();
    double calcLEnergyF_MT();


    void makeIDX();

    void setDeltaF(const CameraCalibration* HCalib );

    void setAdjointsF(const CameraCalibration* HCalib );

    std::vector<EFFrame*> frames;
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
    VecCf cDeltaF;
private:

    VecX getStitchedDeltaF() const;

    void resubstituteF_MT(VecX x, const CameraCalibration* HCalib, bool MT);
    void resubstituteFPt(const VecCf &xc, Mat18f* xAd, int min, int max, Vec10* stats, int tid);

    void accumulateAF_MT(MatXX &H, VecX &b, bool MT);
    void accumulateLF_MT(MatXX &H, VecX &b, bool MT);
    void accumulateSCF_MT(MatXX &H, VecX &b, bool MT);

    void calcLEnergyPt(int min, int max, Vec10* stats, int tid);

    void orthogonalize(VecX* b, MatXX* H);
    Mat18f* adHTdeltaF;

    Mat88f* adHostF;
    Mat88f* adTargetF;


    VecC cPrior;
    VecCf cPriorF;

    dso::AccumulatedTopHessianSSE* accSSE_top_L;
    dso::AccumulatedTopHessianSSE* accSSE_top_A;


    dso::AccumulatedSCHessianSSE* accSSE_bot;

    std::vector<EFPoint*> allPoints;
    std::vector<EFPoint*> allPointsToMarg;

    float currentLambda;
};


#endif //VGIS10_ENERGY_H
