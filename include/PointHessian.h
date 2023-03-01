//
// Created by magnus on 2/10/23.
//

#ifndef VGIS10_POINTHESSIAN_H
#define VGIS10_POINTHESSIAN_H

#include "Types.h"
#include "Optimized/PointFrameResidual.h"

class EFPoint;
class ImmaturePoint;

struct PointHessian {
    PointHessian(const ImmaturePoint* const rawPoint, CalibHessian* Hcalib);

    ~PointHessian() {
        for(unsigned int i=0;i<residuals.size();i++) delete residuals[i];
        residuals.clear();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    //EFPoint* efPoint;

    // static values
    float color[MAX_RES_PER_POINT]{};			// intensity per residual in host frame
    float weights[MAX_RES_PER_POINT]{};		// host-weights for respective residuals.

    float u,v;
    int idx{};
    float energyTH;
    bool hasDepthPrior;
    int pointIndex = 0;
    EFPoint* efPoint;
    VO::Frame* host;

    float my_type;

    float idepth_scaled{};
    float idepth_zero_scaled{};
    float idepth_zero{};
    float idepth{};
    float step{};
    float step_backup{};
    float idepth_backup{};

    float nullspaces_scale{};
    float idepth_hessian;
    float maxRelBaseline;
    int numGoodResiduals;

    std::vector<PointFrameResidual *> residuals;					    // only contains good residuals (not OOB and not OUTLIER). Arbitrary order.
    std::pair<PointFrameResidual* , ResState> lastResiduals[2]; 	// contains information about residuals to the last two (!) frames. ([0] = latest, [1] = the one before).

    bool isGood() const{
        return std::isfinite(energyTH);
    }

    enum PtStatus {ACTIVE=0, INACTIVE, OUTLIER, OOB, MARGINALIZED};
    PtStatus status;

    inline void setPointStatus(PtStatus s) {status=s;}


    inline void setIdepth(float idepth) {
        this->idepth = idepth;
        this->idepth_scaled = SCALE_IDEPTH * idepth;
    }
    inline void setIdepthScaled(float idepth_scaled) {
        this->idepth = SCALE_IDEPTH_INVERSE * idepth_scaled;
        this->idepth_scaled = idepth_scaled;
    }
    inline void setIdepthZero(float idepth) {
        idepth_zero = idepth;
        idepth_zero_scaled = SCALE_IDEPTH * idepth;
        nullspaces_scale = -(idepth*1.001 - idepth/1.001)*500;
    }


    //std::vector<PointFrameResidual*> residuals;					// only contains good residuals (not OOB and not OUTLIER). Arbitrary order.
    //std::pair<PointFrameResidual*, ResState> lastResiduals[2]; 	// contains information about residuals to the last two (!) frames. ([0] = latest, [1] = the one before).


    void release(){}

    [[nodiscard]] bool isOOB(const std::vector<VO::Frame*>& toKeep, const std::vector<VO::Frame*>& toMarg) const;
    bool isInlierNew();





};
#endif //VGIS10_POINTHESSIAN_H
