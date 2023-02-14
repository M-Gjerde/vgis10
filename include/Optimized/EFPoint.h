//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFPOINT_H
#define VGIS10_EFPOINT_H

#include "PointHessian.h"

enum EFPointStatus {PS_GOOD=0, PS_MARGINALIZE, PS_DROP};


class EFPoint {
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

    // contains all residuals.

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


#endif //VGIS10_EFPOINT_H
