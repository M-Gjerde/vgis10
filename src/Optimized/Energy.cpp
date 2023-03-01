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
    auto*  eff = new EFFrame(fh);
    eff->idx = eFrames.size();
    eFrames.push_back(eff);
    fh->efFrame = eff;

    nFrames++;
    //fh->efFrame = eff;

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
    makeIDX();


    for (EFFrame *fh2: eFrames) {
        connectivityMap[(((uint64_t) eff->frameID) << 32) + ((uint64_t) fh2->frameID)] = Eigen::Vector2i(0, 0);
        if (fh2 != eff)
            connectivityMap[(((uint64_t) fh2->frameID) << 32) + ((uint64_t) eff->frameID)] = Eigen::Vector2i(0, 0);
    }
}

void EnergyFunctional::insertResidual(PointFrameResidual *r) {
    EFResidual *efr = new EFResidual(r, r->point->efPoint, r->host->efFrame, r->target->efFrame);

    efr->idxInAll = r->point->efPoint->residualsAll.size();
    r->point->efPoint->residualsAll.push_back(efr);

    connectivityMap[(((uint64_t) r->host->trackingID) << 32) + ((uint64_t) r->target->trackingID)][0]++;
    nResiduals++;
    r->efResidual = efr;

}

void EnergyFunctional::insertPoint(PointHessian *ph) {
    EFPoint* efp = new EFPoint(ph, ph->host->efFrame);
    efp->idxInPoints = ph->pointIndex;
    ph->host->efFrame->points.push_back(efp);
    nPoints++;
    ph->efPoint = efp;
    EFIndicesValid = false;
}

void EnergyFunctional::setAdjointsF(CalibHessian *Hcalib) {
    if (adHost != 0) delete[] adHost;
    if (adTarget != 0) delete[] adTarget;
    adHost = new Mat88[nFrames * nFrames];
    adTarget = new Mat88[nFrames * nFrames];

    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            VO::Frame *host = eFrames[h]->data.get();
            VO::Frame *target = eFrames[t]->data.get();

            SE3 hostToTarget = target->pose->get_worldToCam_evalPT() * host->pose->get_worldToCam_evalPT().inverse();

            Mat88 AH = Mat88::Identity();
            Mat88 AT = Mat88::Identity();

            AH.topLeftCorner<6, 6>() = -hostToTarget.Adj().transpose();
            AT.topLeftCorner<6, 6>() = Mat66::Identity();


            Vec2f affLL = AffLight::fromToVecExposure(host->abExposure, target->abExposure, host->pose->aff_g2l_0(),
                                                      target->pose->aff_g2l_0()).cast<float>();
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
        eFrames[idx]->idx = idx;

    allPoints.clear();

    for (EFFrame *f: eFrames) {
        for (EFPoint *p: f->points) {
            allPoints.push_back(p);
            for (EFResidual *r: p->residualsAll) {
                r->hostIDX = r->host->idx;
                r->targetIDX = r->target->idx;
            }
        }
    }
    EFIndicesValid = true;
}

void EnergyFunctional::dropResidual(EFResidual *r) {
    EFPoint *p = r->point;
    assert(r == p->residualsAll[r->idxInAll]);

    p->residualsAll[r->idxInAll] = p->residualsAll.back();
    p->residualsAll[r->idxInAll]->idxInAll = r->idxInAll;
    p->residualsAll.pop_back();


    connectivityMap[(((uint64_t) r->host->frameID) << 32) + ((uint64_t) r->target->frameID)][0]--;
    nResiduals--;
    r->data->efResidual = 0;
    delete r;
}

double EnergyFunctional::calcLEnergyF_MT() {
    assert(EFDeltaValid);
    assert(EFAdjointsValid);
    assert(EFIndicesValid);

    double E = 0;
    for (EFFrame *f: eFrames) {
        E += f->delta_prior.cwiseProduct(f->prior).dot(f->delta_prior);
    }

    E += cDeltaF.cwiseProduct(cPriorF).dot(cDeltaF);

    red->reduce(boost::bind(&EnergyFunctional::calcLEnergyPt,
                            this, _1, _2, _3, _4), 0, allPoints.size(), 50);
    return E + red->stats[0];

}

void EnergyFunctional::calcLEnergyPt(int min, int max, Vec10 *stats, int tid) {

    dso::Accumulator11 E{};
    E.initialize();
    VecCf dc = cDeltaF;

    for (int i = min; i < max; i++) {
        EFPoint *p = allPoints[i];
        float dd = p->deltaF;

        for (EFResidual *r: p->residualsAll) {
            if (!r->isLinearized || !r->isActive()) continue;

            Mat18f dp = adHTdeltaF[r->hostIDX + nFrames * r->targetIDX];
            if (r->targetIDX != 1 || r->hostIDX != 0)
                Log::Logger::getInstance()->error("targetIDX or host idx is not right {} | {}", r->hostIDX, r->targetIDX);

            RawResidualJacobian * rJ = r->J;



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

                __m128 r0 = _mm_load_ps(((float *) &r->res_toZeroF) + i);
                r0 = _mm_add_ps(r0, r0);
                r0 = _mm_add_ps(r0, Jdelta);
                Jdelta = _mm_mul_ps(Jdelta, r0);
                E.updateSSENoShift(Jdelta);
            }
            for (int i = ((patternNum >> 2) << 2); i < patternNum; i++) {
                float Jdelta = rJ->JIdx[0][i] * Jp_delta_x_1 + rJ->JIdx[1][i] * Jp_delta_y_1 +
                               rJ->JabF[0][i] * dp[6] + rJ->JabF[1][i] * dp[7];
                E.updateSingleNoShift((float) (Jdelta * (Jdelta + 2 * r->res_toZeroF[i])));
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
    for (int h = 0; h < nFrames; h++)
        d.segment<8>(CPARS + 8 * h) = eFrames[h]->delta;
    return d;
}

void EnergyFunctional::setDeltaF(CalibHessian *HCalib) {
    if (adHTdeltaF != 0) delete[] adHTdeltaF;
    adHTdeltaF = new Mat18f[nFrames * nFrames];
    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            int idx = h + t * nFrames;
            adHTdeltaF[idx] =
                    eFrames[h]->data->pose->get_state_minus_stateZero().head<8>().cast<float>().transpose() * adHostF[idx]
                    + eFrames[t]->data->pose->get_state_minus_stateZero().head<8>().cast<float>().transpose() *
                      adTargetF[idx];
        }

    cDeltaF = HCalib->value_minus_value_zero.cast<float>();
    for (EFFrame *f: eFrames) {
        f->delta = f->data->pose->get_state_minus_stateZero().head<8>();
        f->delta_prior = (f->data->pose->get_state() - f->data->pose->getPriorZero()).head<8>();

        for (EFPoint *p: f->points)
            p->deltaF = p->data->idepth - p->data->idepth_zero;
    }

    EFDeltaValid = true;
}

void EnergyFunctional::solveSystemF(int iteration, double lambda, CalibHessian *HCalib) {
    lambda = 1e-5;

    assert(EFDeltaValid);
    assert(EFAdjointsValid);
    assert(EFIndicesValid);

    MatXX HL_top, HA_top, H_sc;
    VecX bL_top, bA_top, bM_top, b_sc;


    accumulateAF_MT(HA_top, bA_top, false);
    accumulateLF_MT(HL_top, bL_top, false);
    accumulateSCF_MT(H_sc, b_sc, false);

    bM_top = (bM + HM * getStitchedDeltaF());
    MatXX HFinal_top;
    VecX bFinal_top;


    HFinal_top = HL_top + HM + HA_top;
    bFinal_top = bL_top + bM_top + bA_top - b_sc;

    lastHS = HFinal_top - H_sc;
    lastbS = bFinal_top;

    for (int i = 0; i < 8 * nFrames + CPARS; i++) HFinal_top(i, i) *= (1 + lambda);
    HFinal_top -= H_sc * (1.0f / (1 + lambda));


    VecX x;
    VecX SVecI = (HFinal_top.diagonal() + VecX::Constant(HFinal_top.cols(), 10)).cwiseSqrt().cwiseInverse();
    MatXX HFinalScaled = SVecI.asDiagonal() * HFinal_top * SVecI.asDiagonal();
    x = SVecI.asDiagonal() *
        HFinalScaled.ldlt().solve(SVecI.asDiagonal() * bFinal_top);//  SVec.asDiagonal() * svd.matrixV() * Ub;

    if (iteration >= 2 ) {
        VecX xOld = x;
        orthogonalize(&x, 0);
    }

    lastX = x;

    //resubstituteF(x, HCalib);
    currentLambda = lambda;
    resubstituteF_MT(x, HCalib, false);
    currentLambda = 0;

}

void EnergyFunctional::accumulateAF_MT(MatXX &H, VecX &b, bool MT) {
    accSSE_top_A->setZero(nFrames);
    for (EFFrame *f: eFrames)
        for (EFPoint *p: f->points)
            accSSE_top_A->addPoint<0>(p, this);
    accSSE_top_A->stitchDoubleMT(red, H, b, this, false, false);
    resInA = accSSE_top_A->nres[0];
}

void EnergyFunctional::accumulateLF_MT(MatXX &H, VecX &b, bool MT) {
    accSSE_top_L->setZero(nFrames);
    for (EFFrame *f: eFrames)
        for (EFPoint *p: f->points)
            accSSE_top_L->addPoint<1>(p, this);
    accSSE_top_L->stitchDoubleMT(red, H, b, this, true, false);
    resInL = accSSE_top_L->nres[0];
}

void EnergyFunctional::accumulateSCF_MT(MatXX &H, VecX &b, bool MT) {
    accSSE_bot->setZero(nFrames);
    for (EFFrame *f: eFrames)
        for (EFPoint *p: f->points)
            accSSE_bot->addPoint(p, true);
    accSSE_bot->stitchDoubleMT(red, H, b, this, false);
}


void EnergyFunctional::resubstituteF_MT(VecX x, CalibHessian *HCalib, bool MT) {
    assert(x.size() == CPARS + nFrames * 8);

    VecXf xF = x.cast<float>();
    HCalib->step = -x.head<CPARS>();

    Mat18f *xAd = new Mat18f[nFrames * nFrames];
    VecCf cstep = xF.head<CPARS>();
    for (EFFrame *h: eFrames) {
        h->data->pose->step.head<8>() = -x.segment<8>(CPARS + 8 * h->idx);
        h->data->pose->step.tail<2>().setZero();

        for (EFFrame *t: eFrames)
            xAd[nFrames * h->idx + t->idx] =
                    xF.segment<8>(CPARS + 8 * h->idx).transpose() * adHostF[h->idx + nFrames * t->idx]
                    + xF.segment<8>(CPARS + 8 * t->idx).transpose() * adTargetF[h->idx + nFrames * t->idx];
    }

    resubstituteFPt(cstep, xAd, 0, allPoints.size(), 0, 0);

    delete[] xAd;
}


void EnergyFunctional::resubstituteFPt(const VecCf &xc, Mat18f *xAd, int min, int max, Vec10 *stats, int tid) {
    for (int k = min; k < max; k++) {
        EFPoint *p = allPoints[k];

        int ngoodres = 0;
        for (EFResidual *r: p->residualsAll) if (r->isActive()) ngoodres++;
        if (ngoodres == 0) {
            p->data->step = 0;
            continue;
        }
        float b = p->bdSumF;
        b -= xc.dot(p->Hcd_accAF + p->Hcd_accLF);

        for (EFResidual *r: p->residualsAll) {
            if (!r->isActive()) continue;
            b -= xAd[r->hostIDX * nFrames + r->targetIDX] * r->JpJdF;
        }

        p->data->step = -b * p->HdiF;
        assert(std::isfinite(p->data->step));
    }
}

void EnergyFunctional::orthogonalize(VecX *b, MatXX *H) {
//	VecX eigenvaluesPre = H.eigenvalues().real();
//	std::sort(eigenvaluesPre.data(), eigenvaluesPre.data()+eigenvaluesPre.size());
//	std::cout << "EigPre:: " << eigenvaluesPre.transpose() << "\n";
    // decide to which nullspaces to orthogonalize.
    std::vector<VecX> ns;
    ns.insert(ns.end(), lastNullspaces_pose.begin(), lastNullspaces_pose.end());
    ns.insert(ns.end(), lastNullspaces_scale.begin(), lastNullspaces_scale.end());
//	if(setting_affineOptModeA <= 0)
//		ns.insert(ns.end(), lastNullspaces_affA.begin(), lastNullspaces_affA.end());
//	if(setting_affineOptModeB <= 0)
//		ns.insert(ns.end(), lastNullspaces_affB.begin(), lastNullspaces_affB.end());
    // make Nullspaces matrix
    MatXX N(ns[0].rows(), ns.size());
    for (unsigned int i = 0; i < ns.size(); i++)
        N.col(i) = ns[i].normalized();
    // compute Npi := N * (N' * N)^-1 = pseudo inverse of N.
    Eigen::JacobiSVD<MatXX> svdNN(N, Eigen::ComputeThinU | Eigen::ComputeThinV);

    VecX SNN = svdNN.singularValues();
    double minSv = 1e10, maxSv = 0;
    for (int i = 0; i < SNN.size(); i++) {
        if (SNN[i] < minSv) minSv = SNN[i];
        if (SNN[i] > maxSv) maxSv = SNN[i];
    }
    for (int i = 0; i < SNN.size(); i++) {
        if (SNN[i] > setting_solverModeDelta * maxSv)
            SNN[i] = 1.0 / SNN[i];
        else SNN[i] = 0;
    }
    MatXX Npi = svdNN.matrixU() * SNN.asDiagonal() * svdNN.matrixV().transpose();    // [dim] x 9.
    MatXX NNpiT = N * Npi.transpose();    // [dim] x [dim].
    MatXX NNpiTS = 0.5 * (NNpiT + NNpiT.transpose());    // = N * (N' * N)^-1 * N'.
    if (b != 0) *b -= NNpiTS * *b;
    if (H != 0) *H -= NNpiTS * *H * NNpiTS;
//	std::cout << std::setprecision(16) << "Orth SV: " << SNN.reverse().transpose() << "\n";

//	VecX eigenvaluesPost = H.eigenvalues().real();
//	std::sort(eigenvaluesPost.data(), eigenvaluesPost.data()+eigenvaluesPost.size());
//	std::cout << "EigPost:: " << eigenvaluesPost.transpose() << "\n";
}

void EnergyFunctional::dropPointsF() {
    for (EFFrame *f: eFrames) {
        for (int i = 0; i < (int) f->points.size(); i++) {
            EFPoint *p = f->points[i];
            if (p->stateFlag == EFPointStatus::PS_DROP) {
                removePoint(p);
                i--;
            }
        }
    }

    EFIndicesValid = false;
    makeIDX();
}

void EnergyFunctional::removePoint(EFPoint *p) {
    for (EFResidual *r: p->residualsAll)
        dropResidual(r);

    EFFrame *h = p->host;
    h->points[p->idxInPoints] = h->points.back();
    h->points[p->idxInPoints]->idxInPoints = p->idxInPoints;
    h->points.pop_back();

    nPoints--;
    p->data->efPoint = 0;

    EFIndicesValid = false;

    delete p;
}

void EnergyFunctional::marginalizeFrame(EFFrame *fh) {

    assert(EFDeltaValid);
    assert(EFAdjointsValid);
    assert(EFIndicesValid);

    assert((int) fh->points.size() == 0);
    int ndim = nFrames * 8 + CPARS - 8;// new dimension
    int odim = nFrames * 8 + CPARS;// old dimension


//	VecX eigenvaluesPre = HM.eigenvalues().real();
//	std::sort(eigenvaluesPre.data(), eigenvaluesPre.data()+eigenvaluesPre.size());
//



    if ((int) fh->idx != (int) eFrames.size() - 1) {
        int io = fh->idx * 8 + CPARS;    // index of frame to move to end
        int ntail = 8 * (nFrames - fh->idx - 1);
        assert((io + 8 + ntail) == nFrames * 8 + CPARS);

        Vec8 bTmp = bM.segment<8>(io);
        VecX tailTMP = bM.tail(ntail);
        bM.segment(io, ntail) = tailTMP;
        bM.tail<8>() = bTmp;

        MatXX HtmpCol = HM.block(0, io, odim, 8);
        MatXX rightColsTmp = HM.rightCols(ntail);
        HM.block(0, io, odim, ntail) = rightColsTmp;
        HM.rightCols(8) = HtmpCol;

        MatXX HtmpRow = HM.block(io, 0, 8, odim);
        MatXX botRowsTmp = HM.bottomRows(ntail);
        HM.block(io, 0, ntail, odim) = botRowsTmp;
        HM.bottomRows(8) = HtmpRow;
    }


//	// marginalize. First add prior here, instead of to active.
    HM.bottomRightCorner<8, 8>().diagonal() += fh->prior;
    bM.tail<8>() += fh->prior.cwiseProduct(fh->delta_prior);



//	std::cout << std::setprecision(16) << "HMPre:\n" << HM << "\n\n";


    VecX SVec = (HM.diagonal().cwiseAbs() + VecX::Constant(HM.cols(), 10)).cwiseSqrt();
    VecX SVecI = SVec.cwiseInverse();


//	std::cout << std::setprecision(16) << "SVec: " << SVec.transpose() << "\n\n";
//	std::cout << std::setprecision(16) << "SVecI: " << SVecI.transpose() << "\n\n";

    // scale!
    MatXX HMScaled = SVecI.asDiagonal() * HM * SVecI.asDiagonal();
    VecX bMScaled = SVecI.asDiagonal() * bM;

    // invert bottom part!
    Mat88 hpi = HMScaled.bottomRightCorner<8, 8>();
    hpi = 0.5f * (hpi + hpi);
    hpi = hpi.inverse();
    hpi = 0.5f * (hpi + hpi);

    // schur-complement!
    MatXX bli = HMScaled.bottomLeftCorner(8, ndim).transpose() * hpi;
    HMScaled.topLeftCorner(ndim, ndim).noalias() -= bli * HMScaled.bottomLeftCorner(8, ndim);
    bMScaled.head(ndim).noalias() -= bli * bMScaled.tail<8>();

    //unscale!
    HMScaled = SVec.asDiagonal() * HMScaled * SVec.asDiagonal();
    bMScaled = SVec.asDiagonal() * bMScaled;

    // set.
    HM = 0.5 * (HMScaled.topLeftCorner(ndim, ndim) + HMScaled.topLeftCorner(ndim, ndim).transpose());
    bM = bMScaled.head(ndim);

    // remove from vector, without changing the order!
    for (unsigned int i = fh->idx; i + 1 < eFrames.size(); i++) {
        eFrames[i] = eFrames[i + 1];
        eFrames[i]->idx = i;
    }
    eFrames.pop_back();
    nFrames--;
    fh->data->efFrame = 0;

    assert((int) frames.size() * 8 + CPARS == (int) HM.rows());
    assert((int) frames.size() * 8 + CPARS == (int) HM.cols());
    assert((int) frames.size() * 8 + CPARS == (int) bM.size());
    assert((int) frames.size() == (int) nFrames);




//	VecX eigenvaluesPost = HM.eigenvalues().real();
//	std::sort(eigenvaluesPost.data(), eigenvaluesPost.data()+eigenvaluesPost.size());

//	std::cout << std::setprecision(16) << "HMPost:\n" << HM << "\n\n";

//	std::cout << "EigPre:: " << eigenvaluesPre.transpose() << "\n";
//	std::cout << "EigPost: " << eigenvaluesPost.transpose() << "\n";

    EFIndicesValid = false;
    EFAdjointsValid = false;
    EFDeltaValid = false;

    makeIDX();
}


void EnergyFunctional::marginalizePointsF() {
    assert(EFDeltaValid);
    assert(EFAdjointsValid);
    assert(EFIndicesValid);


    allPointsToMarg.clear();
    for (EFFrame *f: eFrames) {
        for (int i = 0; i < (int) f->points.size(); i++) {
            EFPoint *p = f->points[i];
            if (p->stateFlag == EFPointStatus::PS_MARGINALIZE) {
                p->priorF *= setting_idepthFixPriorMargFac;
                for (EFResidual *r: p->residualsAll)
                    if (r->isActive())
                        connectivityMap[(((uint64_t) r->host->frameID) << 32) +
                                        ((uint64_t) r->target->frameID)][1]++;
                allPointsToMarg.push_back(p);
            }
        }
    }

    accSSE_bot->setZero(nFrames);
    accSSE_top_A->setZero(nFrames);
    for (EFPoint *p: allPointsToMarg) {
        accSSE_top_A->addPoint<2>(p, this);
        accSSE_bot->addPoint(p, false);
        removePoint(p);
    }
    MatXX M, Msc;
    VecX Mb, Mbsc;
    accSSE_top_A->stitchDouble(M, Mb, this, false, false);
    accSSE_bot->stitchDouble(Msc, Mbsc, this);

    resInM += accSSE_top_A->nres[0];

    MatXX H = M - Msc;
    VecX b = Mb - Mbsc;

    float setting_margWeightFac = 0.5f * 0.5f;

    HM += setting_margWeightFac * H;
    bM += setting_margWeightFac * b;

    EFIndicesValid = false;
    makeIDX();
}

