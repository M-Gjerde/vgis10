//
// Created by magnus on 2/7/23.
//

#ifndef VGIS10_TRACKERINITIALIZER_H
#define VGIS10_TRACKERINITIALIZER_H

#include <memory>

#include "Frame.h"
#include "Point.h"
#include "Optimized/MatrixAccumulators.h"
#include "ImmaturePoint.h"

namespace VO {

    static inline Vec3f calcEC(int lvl, PntResult *pnts,
                 Info *info);


    static Vec3f
    calcResAndGS(int lvl, Mat88f &H_out, Vec8f &b_out, Mat88f &H_out_sc, Vec8f &b_out_sc, const SE3 &refToNew,
                 AffLight refToNew_aff, bool plot, const std::shared_ptr<Frame> &firstFrame,
                 const std::shared_ptr<Frame> &frame,
                 Info *info, PntResult *pnts) {

        int wl = frame->widthLevel[lvl], hl = frame->heightLevel[lvl];
        Eigen::Vector3f *colorRef = firstFrame->pyramid[lvl].data();
        Eigen::Vector3f *colorNew = frame->pyramid[lvl].data();

        /*
        cv::Mat img(hl, wl, CV_32FC3, colorRef);
        img.convertTo(img, CV_8UC3);
        cv::imshow("ref", img);

        cv::Mat img2(hl, wl, CV_32FC3, colorNew);
        img2.convertTo(img2, CV_8UC3);
        cv::imshow("new", img2);
        cv::waitKey(0);
*/

        Mat33f RKi = (refToNew.rotationMatrix() * frame->KiLevel[lvl]).cast<float>();
        Vec3f t = refToNew.translation().cast<float>();
        Eigen::Vector2f r2new_aff = Eigen::Vector2f(exp(refToNew_aff.a), refToNew_aff.b);

        float fxl = frame->fxLevel[lvl];
        float fyl = frame->fyLevel[lvl];
        float cxl = frame->cxLevel[lvl];
        float cyl = frame->cyLevel[lvl];

        dso::Accumulator11 E{};
        info->acc9.initialize();
        E.initialize();

        int npts = pnts->numPointsLevel[lvl];
        Pnt *ptsl = pnts->pointsLevel[lvl].data();
        for (int i = 0; i < npts; i++) {
            Pnt *point = ptsl + i;
            point->maxstep = 1e10;
            if (!point->isGood) {
                E.updateSingle((float) (point->energy[0]));
                point->energy_new = point->energy;
                point->isGood_new = false;
                continue;
            }
            VecNRf dp0;
            VecNRf dp1;
            VecNRf dp2;
            VecNRf dp3;
            VecNRf dp4;
            VecNRf dp5;
            VecNRf dp6;
            VecNRf dp7;
            VecNRf dd;
            VecNRf r;
            info->JbBuffer_new[i].setZero();

            // sum over all residuals.
            bool isGood = true;
            float energy = 0;
            for (int idx = 0; idx < patternNum; idx++) {
                int dx = patternP[idx][0];
                int dy = patternP[idx][1];
                Vec3f pt = RKi * Vec3f(point->u + dx, point->v + dy, 1) + t * point->idepth_new;
                float u = pt[0] / pt[2];
                float v = pt[1] / pt[2];
                float Ku = fxl * u + cxl;
                float Kv = fyl * v + cyl;
                float new_idepth = point->idepth_new / pt[2];
                if (!(Ku > 1 && Kv > 1 && Ku < wl - 2 && Kv < hl - 2 && new_idepth > 0)) {
                    isGood = false;
                    break;
                }
                Vec3f hitColor = getInterpolatedElement33(colorNew, Ku, Kv, wl);
                float rlR = getInterpolatedElement31(colorRef, point->u + dx, point->v + dy, wl);

                if (!std::isfinite(rlR) || !std::isfinite((float) hitColor[0])) {
                    isGood = false;
                    break;
                }

                float residual = hitColor[0] - (r2new_aff[0] * rlR) - r2new_aff[1];
                float hw = fabs(residual) < setting_huberTH ? 1.0f : setting_huberTH / fabs(residual);
                energy += hw * residual * residual * (2.0f - hw);

                float dxdd = (t[0] - t[2] * u) / pt[2];
                float dydd = (t[1] - t[2] * v) / pt[2];

                if (hw < 1) hw = sqrtf(hw);
                float dxInterp = hw * hitColor[1] * fxl;
                float dyInterp = hw * hitColor[2] * fyl;
                dp0[idx] = new_idepth * dxInterp;
                dp1[idx] = new_idepth * dyInterp;
                dp2[idx] = -new_idepth * (u * dxInterp + v * dyInterp);
                dp3[idx] = -u * v * dxInterp - (1 + v * v) * dyInterp;
                dp4[idx] = (1 + u * u) * dxInterp + u * v * dyInterp;
                dp5[idx] = -v * dxInterp + u * dyInterp;
                dp6[idx] = -hw * r2new_aff[0] * rlR;
                dp7[idx] = -hw * 1;
                dd[idx] = dxInterp * dxdd + dyInterp * dydd;
                r[idx] = hw * residual;

                float maxstep = 1.0f / Vec2f(dxdd * fxl, dydd * fyl).norm();
                if (maxstep < point->maxstep) point->maxstep = maxstep;

                // immediately compute dp*dd' and dd*dd' in JbBuffer1.
                info->JbBuffer_new[i][0] += dp0[idx] * dd[idx];
                info->JbBuffer_new[i][1] += dp1[idx] * dd[idx];
                info->JbBuffer_new[i][2] += dp2[idx] * dd[idx];
                info->JbBuffer_new[i][3] += dp3[idx] * dd[idx];
                info->JbBuffer_new[i][4] += dp4[idx] * dd[idx];
                info->JbBuffer_new[i][5] += dp5[idx] * dd[idx];
                info->JbBuffer_new[i][6] += dp6[idx] * dd[idx];
                info->JbBuffer_new[i][7] += dp7[idx] * dd[idx];
                info->JbBuffer_new[i][8] += r[idx] * dd[idx];
                info->JbBuffer_new[i][9] += dd[idx] * dd[idx];
            }
            if (!isGood || energy > point->outlierTH * 20) {
                E.updateSingle((float) (point->energy[0]));
                point->isGood_new = false;
                point->energy_new = point->energy;
                continue;
            }
            // add into energy.
            E.updateSingle(energy);
            point->isGood_new = true;
            point->energy_new[0] = energy;
            // update Hessian matrix.
            for (int i = 0; i + 3 < patternNum; i += 4) {
                info->acc9.updateSSE(
                        _mm_load_ps(((float *) (&dp0)) + i),
                        _mm_load_ps(((float *) (&dp1)) + i),
                        _mm_load_ps(((float *) (&dp2)) + i),
                        _mm_load_ps(((float *) (&dp3)) + i),
                        _mm_load_ps(((float *) (&dp4)) + i),
                        _mm_load_ps(((float *) (&dp5)) + i),
                        _mm_load_ps(((float *) (&dp6)) + i),
                        _mm_load_ps(((float *) (&dp7)) + i),
                        _mm_load_ps(((float *) (&r)) + i));
            }

            for (int i = ((patternNum >> 2) << 2); i < patternNum; i++) {
                info->acc9.updateSingle(
                        (float) dp0[i], (float) dp1[i], (float) dp2[i], (float) dp3[i],
                        (float) dp4[i], (float) dp5[i], (float) dp6[i], (float) dp7[i],
                        (float) r[i]);
            }

        }
        E.finish();
        info->acc9.finish();

        // calculate alpha energy, and decide if we cap it.
        dso::Accumulator11 EAlpha;
        EAlpha.initialize();
        for (int i = 0; i < npts; i++) {
            Pnt *point = ptsl + i;
            if (!point->isGood_new) {
                E.updateSingle((float) (point->energy[1]));
            } else {
                point->energy_new[1] = (point->idepth_new - 1) * (point->idepth_new - 1);
                E.updateSingle((float) (point->energy_new[1]));
            }
        }
        EAlpha.finish();
        float alphaEnergy = info->alphaW * (EAlpha.A + refToNew.translation().squaredNorm() * npts);

        // compute alpha opt.
        float alphaOpt;
        if (alphaEnergy > info->alphaK * npts) {
            alphaOpt = 0;
            alphaEnergy = info->alphaK * npts;
        } else {
            alphaOpt = info->alphaW;
        }

        info->acc9SC.initialize();
        for (int i = 0; i < npts; i++) {
            Pnt *point = ptsl + i;
            if (!point->isGood_new)
                continue;

            point->lastHessian_new = info->JbBuffer_new[i][9];

            info->JbBuffer_new[i][8] += alphaOpt * (point->idepth_new - 1);
            info->JbBuffer_new[i][9] += alphaOpt;

            if (alphaOpt == 0) {
                info->JbBuffer_new[i][8] += info->couplingWeight * (point->idepth_new - point->iR);
                info->JbBuffer_new[i][9] += info->couplingWeight;
            }

            info->JbBuffer_new[i][9] = 1 / (1 + info->JbBuffer_new[i][9]);
            info->acc9SC.updateSingleWeighted(
                    (float) info->JbBuffer_new[i][0], (float) info->JbBuffer_new[i][1],
                    (float) info->JbBuffer_new[i][2], (float) info->JbBuffer_new[i][3],
                    (float) info->JbBuffer_new[i][4], (float) info->JbBuffer_new[i][5],
                    (float) info->JbBuffer_new[i][6], (float) info->JbBuffer_new[i][7],
                    (float) info->JbBuffer_new[i][8], (float) info->JbBuffer_new[i][9]);
        }
        info->acc9SC.finish();

        //printf("nelements in H: {}, in E: {}, in Hsc: {} / 9!\n", (int)acc9.num, (int)E.num, (int)acc9SC.num*9);
        H_out = info->acc9.H.topLeftCorner<8, 8>();// / acc9.num;
        b_out = info->acc9.H.topRightCorner<8, 1>();// / acc9.num;
        H_out_sc = info->acc9SC.H.topLeftCorner<8, 8>();// / acc9.num;
        b_out_sc = info->acc9SC.H.topRightCorner<8, 1>();// / acc9.num;

        H_out(0, 0) += alphaOpt * (float) npts;
        H_out(1, 1) += alphaOpt * (float) npts;
        H_out(2, 2) += alphaOpt * (float) npts;

        Vec3f tlog = refToNew.log().head<3>().cast<float>();
        b_out[0] += tlog[0] * alphaOpt * (float) npts;
        b_out[1] += tlog[1] * alphaOpt * (float) npts;
        b_out[2] += tlog[2] * alphaOpt * (float) npts;


        return Vec3f(E.A, alphaEnergy, E.num);
    }

    static inline Vec3f calcEC(int lvl, PntResult *pnts,
                 Info *info) {
        if (!info->snapped) return Vec3f(0, 0, pnts->numPointsLevel[lvl]);
        dso::AccumulatorX<2> E;
        E.initialize();
        int npts = pnts->numPointsLevel[lvl];
        for (int i = 0; i < npts; i++) {
            Pnt *point = pnts->pointsLevel[lvl].data() + i;
            if (!point->isGood_new) continue;
            float rOld = (point->idepth - point->iR);
            float rNew = (point->idepth_new - point->iR);
            E.updateNoWeight(Vec2f(rOld * rOld, rNew * rNew));

            //printf("%f %f %f!\n", point->idepth, point->idepth_new, point->iR);
        }
        E.finish();

        //printf("ER: %f %f %f!\n", couplingWeight*E.A1m[0], couplingWeight*E.A1m[1], (float)E.num.numIn1m);
        return Vec3f(info->couplingWeight * E.A1m[0], info->couplingWeight * E.A1m[1], E.num);
    }
};
#endif //VGIS10_TRACKERINITIALIZER_H
