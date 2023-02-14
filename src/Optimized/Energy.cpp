//
// Created by magnus on 2/13/23.
//

#include "Optimized/Energy.h"


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

EFFrame *EnergyFunctional::insertFrame(std::shared_ptr<VO::Frame> fh, CalibHessian *Hcalib) {
    EFFrame *eff = new EFFrame(fh);
    eff->idx = frames.size();
    frames.push_back(eff);

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
    makeIDX();


    for (EFFrame *fh2: frames) {
        connectivityMap[(((uint64_t) eff->frameID) << 32) + ((uint64_t) fh2->frameID)] = Eigen::Vector2i(0, 0);
        if (fh2 != eff)
            connectivityMap[(((uint64_t) fh2->frameID) << 32) + ((uint64_t) eff->frameID)] = Eigen::Vector2i(0, 0);
    }

    return eff;
}

EFResidual *EnergyFunctional::insertResidual(PointFrameResidual *r) {
    return nullptr;
}

EFPoint *EnergyFunctional::insertPoint(PointHessian *ph) {
    return nullptr;
}

void EnergyFunctional::setAdjointsF(CalibHessian *Hcalib) {
    if (adHost != 0) delete[] adHost;
    if (adTarget != 0) delete[] adTarget;
    adHost = new Mat88[nFrames * nFrames];
    adTarget = new Mat88[nFrames * nFrames];

    for (int h = 0; h < nFrames; h++)
        for (int t = 0; t < nFrames; t++) {
            VO::Frame *host;// = frames[h]->data;
            VO::Frame *target;// = frames[t]->data;

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
    for (unsigned int idx = 0; idx < frames.size(); idx++)
        frames[idx]->idx = idx;

    allPoints.clear();

    for (EFFrame *f: frames)
        for (EFPoint *p: f->points) {
            allPoints.push_back(p);
            /*
            for (EFResidual *r: p->residualsAll) {
                r->hostIDX = r->host->idx;
                r->targetIDX = r->target->idx;
            }
             */
        }


    EFIndicesValid = true;
}
