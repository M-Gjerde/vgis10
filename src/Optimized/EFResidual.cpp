//
// Created by magnus on 3/1/23.
//

#include "Optimized/EFResidual.h"
#include "Optimized/Energy.h"

void EFResidual::fixLinearizationF(EnergyFunctional *ef) {
    Vec8f dp = ef->adHTdeltaF[hostIDX+ef->nFrames*targetIDX];

    // compute Jp*delta
    __m128 Jp_delta_x = _mm_set1_ps(J->Jpdxi[0].dot(dp.head<6>())
                                    +J->Jpdc[0].dot(ef->cDeltaF)
                                    +J->Jpdd[0]*point->deltaF);
    __m128 Jp_delta_y = _mm_set1_ps(J->Jpdxi[1].dot(dp.head<6>())
                                    +J->Jpdc[1].dot(ef->cDeltaF)
                                    +J->Jpdd[1]*point->deltaF);
    __m128 delta_a = _mm_set1_ps((float)(dp[6]));
    __m128 delta_b = _mm_set1_ps((float)(dp[7]));

    for(int i=0;i<patternNum;i+=4)
    {
        // PATTERN: rtz = resF - [JI*Jp Ja]*delta.
        __m128 rtz = _mm_load_ps(((float*)&J->resF)+i);
        rtz = _mm_sub_ps(rtz,_mm_mul_ps(_mm_load_ps(((float*)(J->JIdx))+i),Jp_delta_x));
        rtz = _mm_sub_ps(rtz,_mm_mul_ps(_mm_load_ps(((float*)(J->JIdx+1))+i),Jp_delta_y));
        rtz = _mm_sub_ps(rtz,_mm_mul_ps(_mm_load_ps(((float*)(J->JabF))+i),delta_a));
        rtz = _mm_sub_ps(rtz,_mm_mul_ps(_mm_load_ps(((float*)(J->JabF+1))+i),delta_b));
        _mm_store_ps(((float*)&res_toZeroF)+i, rtz);
    }

    isLinearized = true;
}
