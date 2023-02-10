//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_HESSIANBLOCKS_H
#define VGIS10_HESSIANBLOCKS_H


#include "Types.h"
#include "Optimized/Energy.h"
#include "ImmaturePoint.h"
// hessian component associated with one point.

namespace VO {
    class Frame;
}


struct CalibHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static int instanceCounter;

    VecC value_zero;
    VecC value_scaled;
    VecCf value_scaledf;
    VecCf value_scaledi;
    VecC value;
    VecC step;
    VecC step_backup;
    VecC value_backup;
    VecC value_minus_value_zero;

    inline ~CalibHessian() {instanceCounter--;}

    inline CalibHessian(const CameraCalibration* calibration)
    {

        VecC initial_value = VecC::Zero();
        initial_value[0] = calibration->fxG[0];
        initial_value[1] = calibration->fyG[0];
        initial_value[2] = calibration->cxG[0];
        initial_value[3] = calibration->cyG[0];


        wM3G = calibration->wG[0] - 3, hM3G = calibration->hG[0] - 3;

        setValueScaled(initial_value);
        value_zero = value;
        value_minus_value_zero.setZero();

        instanceCounter++;
        for(int i=0;i<256;i++)
            Binv[i] = B[i] = i;		// set gamma function to identity
    };
    float wM3G, hM3G;

    // normal mode: use the optimized parameters everywhere!
    inline float& fxl() {return value_scaledf[0];}
    inline float& fyl() {return value_scaledf[1];}
    inline float& cxl() {return value_scaledf[2];}
    inline float& cyl() {return value_scaledf[3];}
    inline float& fxli() {return value_scaledi[0];}
    inline float& fyli() {return value_scaledi[1];}
    inline float& cxli() {return value_scaledi[2];}
    inline float& cyli() {return value_scaledi[3];}

    inline void setValue(const VecC &value)
    {
        // [0-3: Kl, 4-7: Kr, 8-12: l2r]
        this->value = value;
        value_scaled[0] = SCALE_F * value[0];
        value_scaled[1] = SCALE_F * value[1];
        value_scaled[2] = SCALE_C * value[2];
        value_scaled[3] = SCALE_C * value[3];

        this->value_scaledf = this->value_scaled.cast<float>();
        this->value_scaledi[0] = 1.0f / this->value_scaledf[0];
        this->value_scaledi[1] = 1.0f / this->value_scaledf[1];
        this->value_scaledi[2] = - this->value_scaledf[2] / this->value_scaledf[0];
        this->value_scaledi[3] = - this->value_scaledf[3] / this->value_scaledf[1];
        this->value_minus_value_zero = this->value - this->value_zero;
    };

    inline void setValueScaled(const VecC &value_scaled)
    {
        this->value_scaled = value_scaled;
        this->value_scaledf = this->value_scaled.cast<float>();
        value[0] = SCALE_F_INVERSE * value_scaled[0];
        value[1] = SCALE_F_INVERSE * value_scaled[1];
        value[2] = SCALE_C_INVERSE * value_scaled[2];
        value[3] = SCALE_C_INVERSE * value_scaled[3];

        this->value_minus_value_zero = this->value - this->value_zero;
        this->value_scaledi[0] = 1.0f / this->value_scaledf[0];
        this->value_scaledi[1] = 1.0f / this->value_scaledf[1];
        this->value_scaledi[2] = - this->value_scaledf[2] / this->value_scaledf[0];
        this->value_scaledi[3] = - this->value_scaledf[3] / this->value_scaledf[1];
    };


    float Binv[256];
    float B[256];


    EIGEN_STRONG_INLINE float getBGradOnly(float color)
    {
        int c = color+0.5f;
        if(c<5) c=5;
        if(c>250) c=250;
        return B[c+1]-B[c];
    }

    EIGEN_STRONG_INLINE float getBInvGradOnly(float color)
    {
        int c = color+0.5f;
        if(c<5) c=5;
        if(c>250) c=250;
        return Binv[c+1]-Binv[c];
    }
};


struct PointHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static int instanceCounter;
    EFPoint* efPoint;

    // static values
    float color[MAX_RES_PER_POINT];			// colors in host frame
    float weights[MAX_RES_PER_POINT];		// host-weights for respective residuals.



    float u,v;
    int idx;
    float energyTH;
    std::shared_ptr<VO::Frame> host;

    bool hasDepthPrior;

    float my_type;

    float idepth_scaled;
    float idepth_zero_scaled;
    float idepth_zero;
    float idepth;
    float step;
    float step_backup;
    float idepth_backup;

    float nullspaces_scale;
    float idepth_hessian;
    float maxRelBaseline;
    int numGoodResiduals;

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


    std::vector<PointFrameResidual*> residuals;					// only contains good residuals (not OOB and not OUTLIER). Arbitrary order.
    std::pair<PointFrameResidual*, ResState> lastResiduals[2]; 	// contains information about residuals to the last two (!) frames. ([0] = latest, [1] = the one before).


    void release();
    PointHessian(const ImmaturePoint* rawPoint, const CameraCalibration* Hcalib);
    inline ~PointHessian() {assert(efPoint==0); release(); instanceCounter--;}


    inline bool isOOB(const std::vector<std::shared_ptr<VO::Frame>>& toKeep, const std::vector<std::shared_ptr<VO::Frame>>& toMarg) const
    {

        int visInToMarg = 0;
        for(PointFrameResidual* r : residuals)
        {
            if(r->state_state != ResState::IN) continue;
            for(std::shared_ptr<VO::Frame> k : toMarg)
                if(r->target == k) visInToMarg++;
        }
        if((int)residuals.size() >= setting_minGoodActiveResForMarg &&
           numGoodResiduals > setting_minGoodResForMarg+10 &&
           (int)residuals.size()-visInToMarg < setting_minGoodActiveResForMarg)
            return true;





        if(lastResiduals[0].second == ResState::OOB) return true;
        if(residuals.size() < 2) return false;
        if(lastResiduals[0].second == ResState::OUTLIER && lastResiduals[1].second == ResState::OUTLIER) return true;
        return false;
    }


    inline bool isInlierNew()
    {
        return (int)residuals.size() >= setting_minGoodActiveResForMarg
               && numGoodResiduals >= setting_minGoodResForMarg;
    }

};

#endif //VGIS10_HESSIANBLOCKS_H
