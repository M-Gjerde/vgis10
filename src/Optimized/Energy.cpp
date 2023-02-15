//
// Created by magnus on 2/13/23.
//

#include "Optimized/Energy.h"
#include "Optimized/EFFrame.h"


bool EFAdjointsValid = false;
bool EFIndicesValid = false;
bool EFDeltaValid = false;


EnergyFunctional::EnergyFunctional() {
    adHost = 0;
    adTarget = 0;


    red = 0;

    adHostF = 0;
    adTargetF = 0;
    adHTdeltaF = 0;

    nFrames = nResiduals = nPoints = 0;

    HM = MatXX::Zero(CPARS, CPARS);
    bM = VecX::Zero(CPARS);


    accSSE_top_L = new dso::AccumulatedTopHessianSSE();
    accSSE_top_A = new dso::AccumulatedTopHessianSSE();
    accSSE_bot = new dso::AccumulatedSCHessianSSE();

    resInA = resInL = resInM = 0;
    currentLambda = 0;
}

EnergyFunctional::~EnergyFunctional() {
    /*
    for (EFFrame *f: frames) {
        for (EFPoint *p: f->points) {
            for (EFResidual *r: p->residualsAll) {
                r->data->efResidual = 0;
                delete r;
            }
            p->data->efPoint = 0;
            delete p;
        }
        f->data->efFrame = 0;
        delete f;
    }
*/

    if (adHost != 0) delete[] adHost;
    if (adTarget != 0) delete[] adTarget;


    if (adHostF != 0) delete[] adHostF;
    if (adTargetF != 0) delete[] adTargetF;
    if (adHTdeltaF != 0) delete[] adHTdeltaF;


    delete accSSE_top_L;
    delete accSSE_top_A;
    delete accSSE_bot;
}

void EnergyFunctional::insertFrame(std::shared_ptr<VO::Frame> fh, CalibHessian *Hcalib) {
    EFFrame eff(fh);

    eff.idx = eFrames.size();
    eFrames.push_back(eff);

    nFrames++;

    assert(HM.cols() == 8 * nFrames + CPARS - 8);
    bM.conservativeResize(8 * nFrames + CPARS);
    HM.conservativeResize(8 * nFrames + CPARS, 8 * nFrames + CPARS);
    bM.tail<8>().setZero();
    HM.rightCols<8>().setZero();
    HM.bottomRows<8>().setZero();

    EFIndicesValid = false;
    EFAdjointsValid = false;
    EFDeltaValid = false;

    setAdjointsF(Hcalib);
    //makeIDX();


    for (EFFrame &fh2: eFrames) {
        connectivityMap[(((uint64_t) eff.frameID) << 32) + ((uint64_t) fh2.frameID)] = Eigen::Vector2i(0, 0);
        if (fh2.frameID != eff.frameID)
            connectivityMap[(((uint64_t) fh2.frameID) << 32) + ((uint64_t) eff.frameID)] = Eigen::Vector2i(0, 0);
    }
}

void EnergyFunctional::insertResidual(PointFrameResidual *r, int hostFrameID, int targetFrameID) {

    auto& points = eFrames[hostFrameID].points;
    points[r->pointIndex].residualsAll.emplace_back();
    points[r->pointIndex].residualsAll.back().idxInAll = points[r->pointIndex].residualsAll.size() - 1;
    r->efResidual = &points[r->pointIndex].residualsAll.back();

    if( r->pointIndex >= points.size()) {
        Log::Logger::getInstance()->info("Residual point index was larger than ePoints list, {} into {}", r->pointIndex,
                                         points.size());
    return;
    }

    connectivityMap[(((uint64_t) hostFrameID) << 32) + ((uint64_t) targetFrameID)][0]++;
    nResiduals++;
}

void EnergyFunctional::insertPoint(PointHessian *ph, uint32_t i) {
    EFPoint efp(ph);
    efp.idxInPoints = ph->pointIndex;
    eFrames[i].points.push_back(efp);

    nPoints++;
    EFIndicesValid = false;
}

void EnergyFunctional::setAdjointsF(CalibHessian *Hcalib) {
    if (adHost != 0) delete[] adHost;
    if (adTarget != 0) delete[] adTarget;
    adHost = new Mat88[nFrames * nFrames];
    adTarget = new Mat88[nFrames * nFrames];

    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            VO::Frame *host = eFrames[h].data.get();
            VO::Frame *target = eFrames[t].data.get();

            SE3 hostToTarget = target->pose.get_worldToCam_evalPT() * host->pose.get_worldToCam_evalPT().inverse();

            Mat88 AH = Mat88::Identity();
            Mat88 AT = Mat88::Identity();

            AH.topLeftCorner<6, 6>() = -hostToTarget.Adj().transpose();
            AT.topLeftCorner<6, 6>() = Mat66::Identity();


            Vec2f affLL = AffLight::fromToVecExposure(host->abExposure, target->abExposure, host->pose.aff_g2l_0(),
                                                      target->pose.aff_g2l_0()).cast<float>();
            AT(6, 6) = -affLL[0];
            AH(6, 6) = affLL[0];
            AT(7, 7) = -1;
            AH(7, 7) = affLL[0];

            AH.block<3, 8>(0, 0) *= SCALE_XI_TRANS;
            AH.block<3, 8>(3, 0) *= SCALE_XI_ROT;
            AH.block<1, 8>(6, 0) *= SCALE_A;
            AH.block<1, 8>(7, 0) *= SCALE_B;
            AT.block<3, 8>(0, 0) *= SCALE_XI_TRANS;
            AT.block<3, 8>(3, 0) *= SCALE_XI_ROT;
            AT.block<1, 8>(6, 0) *= SCALE_A;
            AT.block<1, 8>(7, 0) *= SCALE_B;

            adHost[h + t * nFrames] = AH;
            adTarget[h + t * nFrames] = AT;
        }
    cPrior = VecC::Constant(setting_initialCalibHessian);


    if (adHostF != 0) delete[] adHostF;
    if (adTargetF != 0) delete[] adTargetF;
    adHostF = new Mat88f[nFrames * nFrames];
    adTargetF = new Mat88f[nFrames * nFrames];

    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            adHostF[h + t * nFrames] = adHost[h + t * nFrames].cast<float>();
            adTargetF[h + t * nFrames] = adTarget[h + t * nFrames].cast<float>();
        }

    cPriorF = cPrior.cast<float>();


    EFAdjointsValid = true;

}

void EnergyFunctional::makeIDX() {
    for (unsigned int idx = 0; idx < eFrames.size(); idx++)
        eFrames[idx].idx = idx;

    allPoints.clear();

    for (EFFrame& f: eFrames)
        for (EFPoint& p: f.points) {
            allPoints.push_back(&p);
        }


    EFIndicesValid = true;
}

void EnergyFunctional::dropResidual(PointFrameResidual * r) {
    assert(r == p->residualsAll[r->idxInAll]);
    auto& points = eFrames[r->trackingID].points;
    EFPoint* p = &points[r->pointIndex];

    int removeAt = -1;
    for(int i = 0; i < p->residualsAll.size(); ++i){
        if ( p->residualsAll[i].idxInAll == r->efResidual->idxInAll)
            removeAt = i;
    }
    if (removeAt != -1){
        auto it = p->residualsAll.begin();
        std::advance(it, removeAt);
        p->residualsAll.erase(it);
    }


    connectivityMap[(((uint64_t) r->host->trackingID) << 32) + ((uint64_t) r->target->trackingID)][0]--;
    nResiduals--;
    r->efResidual = 0;
    delete r;
}

double EnergyFunctional::calcLEnergyF_MT() {
        assert(EFDeltaValid);
        assert(EFAdjointsValid);
        assert(EFIndicesValid);

        double E = 0;
        for (EFFrame& f: eFrames)
            E += f.delta_prior.cwiseProduct(f.prior).dot(f.delta_prior);
        Log::Logger::getInstance()->info("Calculating EnergyF_MT {}", E);

        E += cDeltaF.cwiseProduct(cPriorF).dot(cDeltaF);

        red->reduce(boost::bind(&EnergyFunctional::calcLEnergyPt,
                                this, _1, _2, _3, _4), 0, allPoints.size(), 50);

        Log::Logger::getInstance()->info("Calculating EnergyF_MT {}", E);
        return E + red->stats[0];

}

void EnergyFunctional::calcLEnergyPt(int min, int max, Vec10 *stats, int tid) {

    dso::Accumulator11 E;
    E.initialize();
    VecCf dc = cDeltaF;

    for (int i = min; i < max; i++) {
        EFPoint *p = allPoints[i];
        float dd = p->deltaF;

        for (EFResidual& r: p->residualsAll) {
            if (!r.isLinearized || !r.isActive()) continue;

            Mat18f dp = adHTdeltaF[r.hostIDX + nFrames * r.targetIDX];
            RawResidualJacobian *rJ = r.J;



            // compute Jp*delta
            float Jp_delta_x_1 = rJ->Jpdxi[0].dot(dp.head<6>())
                                 + rJ->Jpdc[0].dot(dc)
                                 + rJ->Jpdd[0] * dd;

            float Jp_delta_y_1 = rJ->Jpdxi[1].dot(dp.head<6>())
                                 + rJ->Jpdc[1].dot(dc)
                                 + rJ->Jpdd[1] * dd;

            __m128 Jp_delta_x = _mm_set1_ps(Jp_delta_x_1);
            __m128 Jp_delta_y = _mm_set1_ps(Jp_delta_y_1);
            __m128 delta_a = _mm_set1_ps((float) (dp[6]));
            __m128 delta_b = _mm_set1_ps((float) (dp[7]));

            for (int i = 0; i + 3 < patternNum; i += 4) {
                // PATTERN: E = (2*res_toZeroF + J*delta) * J*delta.
                __m128 Jdelta = _mm_mul_ps(_mm_load_ps(((float *) (rJ->JIdx)) + i), Jp_delta_x);
                Jdelta = _mm_add_ps(Jdelta, _mm_mul_ps(_mm_load_ps(((float *) (rJ->JIdx + 1)) + i), Jp_delta_y));
                Jdelta = _mm_add_ps(Jdelta, _mm_mul_ps(_mm_load_ps(((float *) (rJ->JabF)) + i), delta_a));
                Jdelta = _mm_add_ps(Jdelta, _mm_mul_ps(_mm_load_ps(((float *) (rJ->JabF + 1)) + i), delta_b));

                __m128 r0 = _mm_load_ps(((float *) &r.res_toZeroF) + i);
                r0 = _mm_add_ps(r0, r0);
                r0 = _mm_add_ps(r0, Jdelta);
                Jdelta = _mm_mul_ps(Jdelta, r0);
                E.updateSSENoShift(Jdelta);
            }
            for (int i = ((patternNum >> 2) << 2); i < patternNum; i++) {
                float Jdelta = rJ->JIdx[0][i] * Jp_delta_x_1 + rJ->JIdx[1][i] * Jp_delta_y_1 +
                               rJ->JabF[0][i] * dp[6] + rJ->JabF[1][i] * dp[7];
                E.updateSingleNoShift((float) (Jdelta * (Jdelta + 2 * r.res_toZeroF[i])));
            }
        }
        E.updateSingle(p->deltaF * p->deltaF * p->priorF);
    }
    E.finish();
    (*stats)[0] += E.A;
}

double EnergyFunctional::calcMEnergyF() {
    assert(EFDeltaValid);
    assert(EFAdjointsValid);
    assert(EFIndicesValid);

    VecX delta = getStitchedDeltaF();
    return delta.dot(2 * bM + HM * delta);
}

VecX EnergyFunctional::getStitchedDeltaF() const {
    VecX d = VecX(CPARS + nFrames * 8);
    d.head<CPARS>() = cDeltaF.cast<double>();
    for (int h = 0; h < nFrames; h++) d.segment<8>(CPARS + 8 * h) = eFrames[h].delta;
    return d;
}

void EnergyFunctional::setDeltaF(CalibHessian *HCalib) {
    if (adHTdeltaF != 0) delete[] adHTdeltaF;
    adHTdeltaF = new Mat18f[nFrames * nFrames];
    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            int idx = h + t * nFrames;
            adHTdeltaF[idx] =
                    eFrames[h].data->pose.get_state_minus_stateZero().head<8>().cast<float>().transpose() * adHostF[idx]
                    + eFrames[t].data->pose.get_state_minus_stateZero().head<8>().cast<float>().transpose() *
                      adTargetF[idx];
        }

    cDeltaF = HCalib->value_minus_value_zero.cast<float>();
    for (EFFrame& f: eFrames) {
        f.delta = f.data->pose.get_state_minus_stateZero().head<8>();
        f.delta_prior = (f.data->pose.get_state() - f.data->pose.getPriorZero()).head<8>();

        for (EFPoint& p: f.points)
             p.deltaF = p.data->idepth - p.data->idepth_zero;
    }

    EFDeltaValid = true;
}
