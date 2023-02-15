//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFRESIDUAL_H
#define VGIS10_EFRESIDUAL_H

#include "Eigen/Eigen"
#include "Optimized/PointFrameResidual.h"

class EFResidual {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    inline EFResidual() {
        isLinearized = false;
        isActiveAndIsGoodNEW = false;
        J = new RawResidualJacobian();
        assert(((long) this) % 16 == 0);
        assert(((long) J) % 16 == 0);
    }

    inline ~EFResidual() {
        delete J;
    }

    void takeDataF(RawResidualJacobian *JHost) {
        std::swap<RawResidualJacobian *>(J, JHost);

        Vec2f JI_JI_Jd = J->JIdx2 * J->Jpdd;

        for (int i = 0; i < 6; i++)
            JpJdF[i] = J->Jpdxi[0][i] * JI_JI_Jd[0] + J->Jpdxi[1][i] * JI_JI_Jd[1];

        JpJdF.segment<2>(6) = J->JabJIdx * J->Jpdd;
    }

    // structural pointers
    int hostIDX, targetIDX;
    int idxInAll;
    RawResidualJacobian *J;
    VecNRf res_toZeroF;
    Vec8f JpJdF;

    // status.
    bool isLinearized;
    // if residual is not OOB & not OUTLIER & should be used during accumulations
    bool isActiveAndIsGoodNEW;

    inline const bool &isActive() const { return isActiveAndIsGoodNEW; }

    void applyRes(bool copyJacobians) {

    }

};

#endif //VGIS10_EFRESIDUAL_H
