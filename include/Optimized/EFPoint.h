//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFPOINT_H
#define VGIS10_EFPOINT_H

#include "PointHessian.h"
#include "Util/Enums.h"
#include "Settings.h"
#include "EFResidual.h"

class EFPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFPoint(PointHessian* ph)
    {
        priorF = ph->hasDepthPrior ? setting_idepthFixPrior*SCALE_IDEPTH*SCALE_IDEPTH : 0;

        deltaF = ph->idepth - ph->idepth_zero;
        stateFlag=EFPointStatus::PS_GOOD;
        data = ph;
    }

    std::vector<EFResidual> residualsAll;
    PointHessian* data;

    float priorF;
    float deltaF;


    // constant info (never changes in-between).
    int idxInPoints{};

    float bdSumF{};
    float HdiF{};
    float Hdd_accLF{};
    VecCf Hcd_accLF;
    float bd_accLF{};
    float Hdd_accAF{};
    VecCf Hcd_accAF;
    float bd_accAF{};


    EFPointStatus stateFlag;
};

#endif //VGIS10_EFPOINT_H
