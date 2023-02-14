//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFRESIDUAL_H
#define VGIS10_EFRESIDUAL_H

#include "Eigen/Eigen"
#include "Optimized/PointFrameResidual.h"
#include "EFPoint.h"
#include "EFFrame.h"

class EFResidual {
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

#endif //VGIS10_EFRESIDUAL_H
