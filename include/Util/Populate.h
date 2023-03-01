//
// Created by magnus on 2/14/23.
//

#ifndef VGIS10_POPULATE_H
#define VGIS10_POPULATE_H

#include "CalibHessian.h"
#include "CameraCalibration.h"
#include "PointHessian.h"
#include "ImmaturePoint.h"
#include "Enums.h"

namespace VO {
    static void applyResidual(bool copyJacobians, EFResidual* efResidual, PointFrameResidual* res){

    }

    static void pointHessianFromImmaturePoint(const ImmaturePoint * rawPoint, PointHessian* ph){
        ph->hasDepthPrior=false;
        ph->idepth_hessian=0;
        ph->maxRelBaseline=0;
        ph->numGoodResiduals=0;

        // set static values & initialization.
        ph->u = rawPoint->u;
        ph->v = rawPoint->v;
        assert(std::isfinite(rawPoint->idepth_max));
        //idepth_init = rawPoint->idepth_GT;

        ph->my_type = rawPoint->type;
        ph->setIdepthScaled((rawPoint->idepth_max + rawPoint->idepth_min)*0.5);
        ph->setPointStatus(PointHessian::INACTIVE);

        int n = patternNum;
        memcpy(ph->color, rawPoint->color, sizeof(float)*n);
        memcpy(ph->weights, rawPoint->weights, sizeof(float)*n);
        ph->energyTH = rawPoint->energyTH;
    }
    static inline void initCalibHessian(CalibHessian* calibHessian, const CameraCalibration* calibration) {
        VecC initial_value = VecC::Zero();
        initial_value[0] = calibration->fxG[0];
        initial_value[1] = calibration->fyG[0];
        initial_value[2] = calibration->cxG[0];
        initial_value[3] = calibration->cyG[0];


        // WidthMinus3Global

        calibHessian->setValueScaled(initial_value);
        calibHessian->value_minus_value_zero.setZero();

        for(int i=0;i<256;i++)
            calibHessian->Binv[i] = calibHessian->B[i] = i;		// set gamma function to identity
    }

}

#endif //VGIS10_POPULATE_H
