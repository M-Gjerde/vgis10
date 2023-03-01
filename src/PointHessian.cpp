//
// Created by magnus on 3/1/23.
//

#include "PointHessian.h"
#include "ImmaturePoint.h"

PointHessian::PointHessian(const ImmaturePoint *const rawPoint, CalibHessian *Hcalib) {
    host = rawPoint->host;
    hasDepthPrior = false;

    idepth_hessian = 0;
    maxRelBaseline = 0;
    numGoodResiduals = 0;

// set static values & initialization.
    u = rawPoint->u;
    v = rawPoint->v;
    assert(std::isfinite(rawPoint->idepth_max));
//idepth_init = rawPoint->idepth_GT;

    my_type = rawPoint->my_type;

    setIdepthScaled((rawPoint->idepth_max + rawPoint->idepth_min) * 0.5);
    setPointStatus(PointHessian::INACTIVE);

    int n = patternNum;
    memcpy(color, rawPoint->color, sizeof(float) * n);
    memcpy(weights, rawPoint->weights, sizeof(float) * n);
    energyTH = rawPoint->energyTH;

    efPoint = 0;


}

bool PointHessian::isOOB(const std::vector<VO::Frame *> &toKeep, const std::vector<VO::Frame *> &toMarg) const {

    int visInToMarg = 0;
    for (PointFrameResidual *r: residuals) {
        if (r->state_state != ResState::IN) continue;
        for (auto &k: toMarg)
            if (r->target == k) visInToMarg++;
    }
    if ((int) residuals.size() >= setting_minGoodActiveResForMarg &&
        numGoodResiduals > setting_minGoodResForMarg + 10 &&
        (int) residuals.size() - visInToMarg < setting_minGoodActiveResForMarg)
        return true;


    if (lastResiduals[0].second == ResState::OOB) return true;
    if (residuals.size() < 2) return false;
    if (lastResiduals[0].second == ResState::OUTLIER && lastResiduals[1].second == ResState::OUTLIER) return true;
    return false;
}


bool PointHessian::isInlierNew() {
    return (int) residuals.size() >= setting_minGoodActiveResForMarg
           && numGoodResiduals >= setting_minGoodResForMarg;
}
