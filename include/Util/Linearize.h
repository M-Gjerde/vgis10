//
// Created by magnus on 2/14/23.
//

#ifndef VGIS10_LINEARIZE_H
#define VGIS10_LINEARIZE_H
#include "Frame.h"
#include "PointHessian.h"
#include "CalibHessian.h"
#include "Optimized/PointFrameResidual.h"

EIGEN_STRONG_INLINE bool projectPoint(
        const float &u_pt,const float &v_pt,
        const float &idepth,
        const Mat33f &KRKi, const Vec3f &Kt,
        float &Ku, float &Kv, float wM3G, float hM3G)
{
    Vec3f ptp = KRKi * Vec3f(u_pt,v_pt, 1) + Kt*idepth;
    Ku = ptp[0] / ptp[2];
    Kv = ptp[1] / ptp[2];
    return Ku>1.1f && Kv>1.1f && Ku<wM3G && Kv<hM3G;
}


EIGEN_STRONG_INLINE float derive_idepth(
        const Vec3f &t, const float &u, const float &v,
        const int &dx, const int &dy, const float &dxInterp,
        const float &dyInterp, const float &drescale)
{
    return (dxInterp*drescale * (t[0]-t[2]*u)
            + dyInterp*drescale * (t[1]-t[2]*v))*SCALE_IDEPTH;
}

EIGEN_STRONG_INLINE bool projectPoint(
        const float &u_pt,const float &v_pt,
        const float &idepth,
        const int &dx, const int &dy,
        CalibHessian* const &HCalib,
        const Mat33f &R, const Vec3f &t,
        float &drescale, float &u, float &v,
        float &Ku, float &Kv, Vec3f &KliP, float &new_idepth, float wM3G, float hM3G)
{
    KliP = Vec3f(
            (u_pt+dx-HCalib->cxl())*HCalib->fxli(),
            (v_pt+dy-HCalib->cyl())*HCalib->fyli(),
            1);

    Vec3f ptp = R * KliP + t*idepth;
    drescale = 1.0f/ptp[2];
    new_idepth = idepth*drescale;

    if(!(drescale>0)) return false;

    u = ptp[0] * drescale;
    v = ptp[1] * drescale;
    Ku = u*HCalib->fxl() + HCalib->cxl();
    Kv = v*HCalib->fyl() + HCalib->cyl();

    return Ku>1.1f && Kv>1.1f && Ku<wM3G && Kv<hM3G;
}



static double linearize(PointHessian* point, CalibHessian* HCalib, PointFrameResidual* pointRes, const std::shared_ptr<VO::Frame>& host, const std::shared_ptr<VO::Frame>& target){
    pointRes->state_NewEnergyWithOutlier=-1;

    if(pointRes->state_state == ResState::OOB)
    { pointRes->state_NewState = ResState::OOB; return pointRes->state_energy; }
    //Log::Logger::getInstance()->info("This: {}", static_cast<void *>(host.get()));
    VO::FrameToFramePrecalc* precalc = &(host->targetPrecalc[target->trackingID]);
    float energyLeft=0;
    const Eigen::Vector3f* dIl = target->pyramid[0].data();
    //const float* const Il = target->I;
    const Mat33f &PRE_KRKiTll = precalc->PRE_KRKiTll;
    const Vec3f &PRE_KtTll = precalc->PRE_KtTll;
    const Mat33f &PRE_RTll_0 = precalc->PRE_RTll_0;
    const Vec3f &PRE_tTll_0 = precalc->PRE_tTll_0;
    const float * const color = point->color;
    const float * const weights = point->weights;

    Vec2f affLL = precalc->PRE_aff_mode;
    float b0 = precalc->PRE_b0_mode;


    Vec6f d_xi_x, d_xi_y;
    Vec4f d_C_x, d_C_y;
    float d_d_x, d_d_y;
    {
        float drescale, u, v, new_idepth;
        float Ku, Kv;
        Vec3f KliP;

        if(!projectPoint(point->u, point->v, point->idepth_zero_scaled, 0, 0,HCalib,
                         PRE_RTll_0,PRE_tTll_0, drescale, u, v, Ku, Kv, KliP, new_idepth, host->width-3, host->height - 3))
        { pointRes->state_NewState = ResState::OOB; return pointRes->state_energy; }

        pointRes->centerProjectedTo = Vec3f(Ku, Kv, new_idepth);


        // diff d_idepth
        d_d_x = drescale * (PRE_tTll_0[0]-PRE_tTll_0[2]*u)*SCALE_IDEPTH*HCalib->fxl();
        d_d_y = drescale * (PRE_tTll_0[1]-PRE_tTll_0[2]*v)*SCALE_IDEPTH*HCalib->fyl();




        // diff calib
        d_C_x[2] = drescale*(PRE_RTll_0(2,0)*u-PRE_RTll_0(0,0));
        d_C_x[3] = HCalib->fxl() * drescale*(PRE_RTll_0(2,1)*u-PRE_RTll_0(0,1)) * HCalib->fyli();
        d_C_x[0] = KliP[0]*d_C_x[2];
        d_C_x[1] = KliP[1]*d_C_x[3];

        d_C_y[2] = HCalib->fyl() * drescale*(PRE_RTll_0(2,0)*v-PRE_RTll_0(1,0)) * HCalib->fxli();
        d_C_y[3] = drescale*(PRE_RTll_0(2,1)*v-PRE_RTll_0(1,1));
        d_C_y[0] = KliP[0]*d_C_y[2];
        d_C_y[1] = KliP[1]*d_C_y[3];

        d_C_x[0] = (d_C_x[0]+u)*SCALE_F;
        d_C_x[1] *= SCALE_F;
        d_C_x[2] = (d_C_x[2]+1)*SCALE_C;
        d_C_x[3] *= SCALE_C;

        d_C_y[0] *= SCALE_F;
        d_C_y[1] = (d_C_y[1]+v)*SCALE_F;
        d_C_y[2] *= SCALE_C;
        d_C_y[3] = (d_C_y[3]+1)*SCALE_C;


        d_xi_x[0] = new_idepth*HCalib->fxl();
        d_xi_x[1] = 0;
        d_xi_x[2] = -new_idepth*u*HCalib->fxl();
        d_xi_x[3] = -u*v*HCalib->fxl();
        d_xi_x[4] = (1+u*u)*HCalib->fxl();
        d_xi_x[5] = -v*HCalib->fxl();

        d_xi_y[0] = 0;
        d_xi_y[1] = new_idepth*HCalib->fyl();
        d_xi_y[2] = -new_idepth*v*HCalib->fyl();
        d_xi_y[3] = -(1+v*v)*HCalib->fyl();
        d_xi_y[4] = u*v*HCalib->fyl();
        d_xi_y[5] = u*HCalib->fyl();
    }


    {
        pointRes->J->Jpdxi[0] = d_xi_x;
        pointRes->J->Jpdxi[1] = d_xi_y;

        pointRes->J->Jpdc[0] = d_C_x;
        pointRes->J->Jpdc[1] = d_C_y;

        pointRes->J->Jpdd[0] = d_d_x;
        pointRes->J->Jpdd[1] = d_d_y;

    }






    float JIdxJIdx_00=0, JIdxJIdx_11=0, JIdxJIdx_10=0;
    float JabJIdx_00=0, JabJIdx_01=0, JabJIdx_10=0, JabJIdx_11=0;
    float JabJab_00=0, JabJab_01=0, JabJab_11=0;

    float wJI2_sum = 0;

    for(int idx=0;idx<patternNum;idx++)
    {
        float Ku, Kv;
        if(!projectPoint(point->u+patternP[idx][0], point->v+patternP[idx][1], point->idepth_scaled, PRE_KRKiTll, PRE_KtTll, Ku, Kv, host->width - 3, host->height - 3))
        { pointRes->state_NewState = ResState::OOB; return pointRes->state_energy; }

        pointRes->projectedTo[idx][0] = Ku;
        pointRes->projectedTo[idx][1] = Kv;


        Vec3f hitColor = (VO::getInterpolatedElement33(dIl, Ku, Kv, host->width));
        float residual = hitColor[0] - (float)(affLL[0] * color[idx] + affLL[1]);

        //printf("Hitcolor %f %f %f at %f, %f\n", hitColor.x(), hitColor.y(), hitColor.z(), Ku, Kv);


        float drdA = (color[idx]-b0);
        if(!std::isfinite((float)hitColor[0]))
        { pointRes->state_NewState = ResState::OOB; return pointRes->state_energy; }


        float w = sqrtf(setting_outlierTHSumComponent / (setting_outlierTHSumComponent + hitColor.tail<2>().squaredNorm()));
        w = 0.5f*(w + weights[idx]);



        float hw = fabsf(residual) < setting_huberTH ? 1 : setting_huberTH / fabsf(residual);
        energyLeft += w*w*hw *residual*residual*(2-hw);

        {
            if(hw < 1) hw = sqrtf(hw);
            hw = hw*w;

            hitColor[1]*=hw;
            hitColor[2]*=hw;

            pointRes->J->resF[idx] = residual*hw;

            pointRes->J->JIdx[0][idx] = hitColor[1];
            pointRes->J->JIdx[1][idx] = hitColor[2];
            pointRes->J->JabF[0][idx] = drdA*hw; // ab prior thing?
            pointRes->J->JabF[1][idx] = hw;

            JIdxJIdx_00+=hitColor[1]*hitColor[1];
            JIdxJIdx_11+=hitColor[2]*hitColor[2];
            JIdxJIdx_10+=hitColor[1]*hitColor[2];

            JabJIdx_00+= drdA*hw * hitColor[1];
            JabJIdx_01+= drdA*hw * hitColor[2];
            JabJIdx_10+= hw * hitColor[1];
            JabJIdx_11+= hw * hitColor[2];

            JabJab_00+= drdA*drdA*hw*hw;
            JabJab_01+= drdA*hw*hw;
            JabJab_11+= hw*hw;


            wJI2_sum += hw*hw*(hitColor[1]*hitColor[1]+hitColor[2]*hitColor[2]);

        }
    }

    pointRes->J->JIdx2(0,0) = JIdxJIdx_00;
    pointRes->J->JIdx2(0,1) = JIdxJIdx_10;
    pointRes->J->JIdx2(1,0) = JIdxJIdx_10;
    pointRes->J->JIdx2(1,1) = JIdxJIdx_11;
    pointRes->J->JabJIdx(0,0) = JabJIdx_00;
    pointRes->J->JabJIdx(0,1) = JabJIdx_01;
    pointRes->J->JabJIdx(1,0) = JabJIdx_10;
    pointRes->J->JabJIdx(1,1) = JabJIdx_11;
    pointRes->J->Jab2(0,0) = JabJab_00;
    pointRes->J->Jab2(0,1) = JabJab_01;
    pointRes->J->Jab2(1,0) = JabJab_01;
    pointRes->J->Jab2(1,1) = JabJab_11;

    pointRes->state_NewEnergyWithOutlier = energyLeft;

    if(energyLeft > std::max<float>(host->frameEnergyTH, target->frameEnergyTH) || wJI2_sum < 2)
    {
        energyLeft = std::max<float>(host->frameEnergyTH, target->frameEnergyTH);
        pointRes->state_NewState = ResState::OUTLIER;
    }
    else
    {
        pointRes->state_NewState = ResState::IN;
    }

    pointRes->state_NewEnergy = energyLeft;
    return energyLeft;
}


#endif //VGIS10_LINEARIZE_H
