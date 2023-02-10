//
// Created by magnus on 2/9/23.
//

#ifndef VGIS10_IMMATUREPOINT_H
#define VGIS10_IMMATUREPOINT_H

#include "Frame.h"


    enum ImmaturePointStatus {
        IPS_GOOD = 0,                    // traced well and good
        IPS_OOB,                    // OOB: end tracking & marginalize!
        IPS_OUTLIER,                // energy too high: if happens again: outlier!
        IPS_SKIPPED,                // traced well and good (but not actually traced).
        IPS_BADCONDITION,            // not traced because of bad condition.
        IPS_UNINITIALIZED            // not even traced once.
    };


    struct ImmaturePoint {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        // static values
        ImmaturePoint(float uCoord, float vCoord, float t, std::shared_ptr<VO::Frame> frame) : u(uCoord), v(vCoord), type(t), host(frame)  {

            gradH.setZero();
            for(int idx=0;idx<patternNum;idx++)
            {
                int dx = patternP[idx][0];
                int dy = patternP[idx][1];
                Vec3f ptc = VO::getInterpolatedElement33BiLin(frame->color.data(), u+dx, v+dy,frame->width);

                color[idx] = ptc[0];
                if(!std::isfinite(color[idx])) {energyTH=NAN; return;}
                gradH += ptc.tail<2>()  * ptc.tail<2>().transpose();
                weights[idx] = sqrtf(setting_outlierTHSumComponent / (setting_outlierTHSumComponent + ptc.tail<2>().squaredNorm()));
            }

            energyTH = patternNum*setting_outlierTH;
            energyTH *= setting_overallEnergyTHWeight*setting_overallEnergyTHWeight;

            quality=10000;
        }
        float color[MAX_RES_PER_POINT]{};
        float weights[MAX_RES_PER_POINT]{};

        Mat22f gradH;
        Vec2f gradH_ev;
        Mat22f gradH_eig;

        float energyTH{};
        float u, v;
        int idxInImmaturePoints{};
        std::shared_ptr<VO::Frame> host;
        float quality{};

        float type;

        float idepth_min{};
        float idepth_max{};

        ImmaturePointStatus lastTraceStatus;
        Vec2f lastTraceUV;
        float lastTracePixelInterval{};


    };


#endif //VGIS10_IMMATUREPOINT_H
