//
// Created by magnus on 2/15/23.
//

#include "Tracker.h"
#include "Util/Linearize.h"

float Tracker::optimize(int mnumOptIts) {

    if (frameHessians.size() < 2) return 0;
    if (frameHessians.size() < 3) mnumOptIts = 20;
    if (frameHessians.size() < 4) mnumOptIts = 15;


    // get statistics and active residuals.

    activeResiduals.clear();
    int numPoints = 0;
    int numLRes = 0;
    for(auto& fh : frameHessians)
        for(PointHessian* ph : fh->pointHessians)
        {
            for(PointFrameResidual* r : ph->residuals)
            {
                if(!r->efResidual->isLinearized)
                {
                    activeResiduals.push_back(r);
                    r->resetOOB();
                }
                else
                    numLRes++;
            }
            numPoints++;
        }


    Log::Logger::getInstance()->info("OPTIMIZE With {} pts, {} active res, {} lin res!", ef.nPoints,
                                     activeResiduals.size(), numLRes);


    Vec3 lastEnergy = linearizeAll(false);
    Log::Logger::getInstance()->info("LastEnergy: {}", lastEnergy[0]);
    double lastEnergyL = ef.calcLEnergyF_MT();
    double lastEnergyM = ef.calcMEnergyF();
    Log::Logger::getInstance()->info("Last Energy L: {} and M: {}", lastEnergyL, lastEnergyM);

    applyRes_Reductor(true,0,activeResiduals.size(),0,0);

    Log::Logger::getInstance()->info("INITIAL ERROR \t A({})=(AV {:.3f}). Num: A({}) + M({}); ab {} {}! ",
                                     lastEnergy[0],
                                     sqrtf((float) (lastEnergy[0] / (patternNum * ef.resInA))),
                                     ef.resInA,
                                     ef.resInM,
                                     frameHessians.back()->aff_g2l().a,
                                     frameHessians.back()->aff_g2l().b);

    double lambda = 1e-1;
    float stepsize = 1;
    VecX previousX = VecX::Constant(CPARS + 8 * frameHessians.size(), NAN);
    for (int iteration = 0; iteration < mnumOptIts; iteration++) {
        // solve!

        backupState(iteration != 0);
        //solveSystemNew(0);
        solveSystem(iteration, lambda);
        double incDirChange = (1e-20 + previousX.dot(ef.lastX)) / (1e-20 + previousX.norm() * ef.lastX.norm());
        previousX = ef.lastX;

        bool canbreak = doStepFromBackup(stepsize, stepsize, stepsize, stepsize, stepsize);

        // eval new energy!
        Vec3 newEnergy = linearizeAll(false);
        double newEnergyL = ef.calcLEnergyF_MT();
        double newEnergyM = ef.calcMEnergyF();


        Log::Logger::getInstance()->info("{} {}. OldEnergy: {} {} {}, NewEnergy {} {} {} \t (L {:.2f}, dir {:.2f}, ss {:.1f}): \t",
                                         (newEnergy[0] + newEnergy[1] + newEnergyL + newEnergyM <
                                          lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM) ? "ACCEPT"
                                                                                                     : "REJECT",
                                         iteration,
                                         lastEnergy[0],
                                         lastEnergyL,
                                         lastEnergyM,
                                         newEnergy[0],
                                         newEnergyL,
                                         newEnergyM,
                                         log10(lambda),
                                         incDirChange,
                                         stepsize);

        Log::Logger::getInstance()->info("A({})=(AV {:.3f}). Num: A({}) + M({}); ab {} {}! ",
                                         newEnergy[0],
                                         sqrtf((float) (newEnergy[0] / (patternNum * ef.resInA))),
                                         ef.resInA,
                                         ef.resInM,
                                         frameHessians.back()->aff_g2l().a,
                                         frameHessians.back()->aff_g2l().b);

        if (setting_forceAceptStep || (newEnergy[0] + newEnergy[1] + newEnergyL + newEnergyM <
                                       lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM)) {


            applyRes_Reductor(true,0,activeResiduals.size(),0,0);


            lastEnergy = newEnergy;
            lastEnergyL = newEnergyL;
            lastEnergyM = newEnergyM;

            lambda *= 0.25;
        } else {
            loadSateBackup();
            lastEnergy = linearizeAll(false);
            lastEnergyL = ef.calcLEnergyF_MT();
            lastEnergyM = ef.calcMEnergyF();
            lambda *= 1e2;
        }


        if (canbreak && iteration >= setting_minOptIterations) {
            Log::Logger::getInstance()->info("Breaking GN optimization. Itration: {}, canBreak: {}", iteration, canbreak);
            break;
        }
    }


    Vec10 newStateZero = Vec10::Zero();
    newStateZero.segment<2>(6) = frameHessians.back()->get_state().segment<2>(6);

    frameHessians.back()->setEvalPT(frameHessians.back()->PRE_worldToCam,
                                         newStateZero);
    EFDeltaValid = false;
    EFAdjointsValid = false;
    ef.setAdjointsF(&hCalib);
    setPrecalcValues();


    lastEnergy = linearizeAll(true);
    Log::Logger::getInstance()->info("Forxed Linearization. Last energy: {}", lastEnergy[0]);


    if (!std::isfinite((double) lastEnergy[0]) || !std::isfinite((double) lastEnergy[1]) ||
        !std::isfinite((double) lastEnergy[2])) {
        Log::Logger::getInstance()->error("!!!!KF Tracking failed: LOST!!!!!!\n");
        //isLost = true;
    }


    //statistics_lastFineTrackRMSE = sqrtf((float) (lastEnergy[0] / (patternNum * ef->resInA)));

    /*
    if (calibLog != 0) {
        (*calibLog) << hCalib.value_scaled.transpose() <<
                    " " << frameHessians.back()->get_state_scaled().transpose() <<
                    " " << sqrtf((float) (lastEnergy[0] / (patternNum * ef->resInA))) <<
                    " " << ef->resInM << "\n";
        calibLog->flush();
    }
*/
    {
        //boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
        for (auto &fh: frameHessians) {
            fh->shell->camToWorld = fh->PRE_camToWorld; // TODO double check
            fh->shell->aff_g2l = fh->aff_g2l();
        }
    }


    //debugPlotTracking();

    return sqrtf((float) (lastEnergy[0] / (patternNum * ef.resInA)));
}

void Tracker::solveSystem(int iteration, double lambda) {
    ef.lastNullspaces_forLogging = getNullspaces(
            ef.lastNullspaces_pose,
            ef.lastNullspaces_scale,
            ef.lastNullspaces_affA,
            ef.lastNullspaces_affB);

    ef.solveSystemF(iteration, lambda, &hCalib);
}


std::vector<VecX> Tracker::getNullspaces(
        std::vector<VecX> &nullspaces_pose,
        std::vector<VecX> &nullspaces_scale,
        std::vector<VecX> &nullspaces_affA,
        std::vector<VecX> &nullspaces_affB) {
    nullspaces_pose.clear();
    nullspaces_scale.clear();
    nullspaces_affA.clear();
    nullspaces_affB.clear();


    int n = CPARS + frameHessians.size() * 8;
    std::vector<VecX> nullspaces_x0_pre;
    for (int i = 0; i < 6; i++) {
        VecX nullspace_x0(n);
        nullspace_x0.setZero();
        for (auto &fh: frameHessians) {
            nullspace_x0.segment<6>(CPARS + fh->trackingID * 8) = fh->nullspaces_pose.col(i);
            nullspace_x0.segment<3>(CPARS + fh->trackingID * 8) *= SCALE_XI_TRANS_INVERSE;
            nullspace_x0.segment<3>(CPARS + fh->trackingID * 8 + 3) *= SCALE_XI_ROT_INVERSE;
        }
        nullspaces_x0_pre.push_back(nullspace_x0);
        nullspaces_pose.push_back(nullspace_x0);
    }
    for (int i = 0; i < 2; i++) {
        VecX nullspace_x0(n);
        nullspace_x0.setZero();
        for (auto &fh: frameHessians) {
            nullspace_x0.segment<2>(CPARS + fh->trackingID * 8 + 6) = fh->nullspaces_affine.col(i).head<2>();
            nullspace_x0[CPARS + fh->trackingID * 8 + 6] *= SCALE_A_INVERSE;
            nullspace_x0[CPARS + fh->trackingID * 8 + 7] *= SCALE_B_INVERSE;
        }
        nullspaces_x0_pre.push_back(nullspace_x0);
        if (i == 0) nullspaces_affA.push_back(nullspace_x0);
        if (i == 1) nullspaces_affB.push_back(nullspace_x0);
    }

    VecX nullspace_x0(n);
    nullspace_x0.setZero();
    for (auto &fh: frameHessians) {
        nullspace_x0.segment<6>(CPARS + fh->trackingID * 8) = fh->nullspaces_scale;
        nullspace_x0.segment<3>(CPARS + fh->trackingID * 8) *= SCALE_XI_TRANS_INVERSE;
        nullspace_x0.segment<3>(CPARS + fh->trackingID * 8 + 3) *= SCALE_XI_ROT_INVERSE;
    }
    nullspaces_x0_pre.push_back(nullspace_x0);
    nullspaces_scale.push_back(nullspace_x0);

    return nullspaces_x0_pre;
}

// applies step to linearization point.
bool Tracker::doStepFromBackup(float stepfacC, float stepfacT, float stepfacR, float stepfacA, float stepfacD) {
//	float meanStepC=0,meanStepP=0,meanStepD=0;
//	meanStepC += Hcalib.step.norm();

    Vec10 pstepfac;
    pstepfac.segment<3>(0).setConstant(stepfacT);
    pstepfac.segment<3>(3).setConstant(stepfacR);
    pstepfac.segment<4>(6).setConstant(stepfacA);


    float sumA = 0, sumB = 0, sumT = 0, sumR = 0, sumID = 0, numID = 0;

    float sumNID = 0;

    hCalib.setValue(hCalib.value_backup + stepfacC*hCalib.step);
    for(auto& fh : frameHessians)
    {
        fh->setState(fh->state_backup + pstepfac.cwiseProduct(fh->step));
        sumA += fh->step[6]*fh->step[6];
        sumB += fh->step[7]*fh->step[7];
        sumT += fh->step.segment<3>(0).squaredNorm();
        sumR += fh->step.segment<3>(3).squaredNorm();

        for(PointHessian* ph : fh->pointHessians)
        {
            ph->setIdepth(ph->idepth_backup + stepfacD*ph->step);
            sumID += ph->step*ph->step;
            sumNID += fabsf(ph->idepth_backup);
            numID++;

            ph->setIdepthZero(ph->idepth_backup + stepfacD*ph->step);
        }
    }

    sumA /= frameHessians.size();
    sumB /= frameHessians.size();
    sumR /= frameHessians.size();
    sumT /= frameHessians.size();
    sumID /= numID;
    sumNID /= numID;


    Log::Logger::getInstance()->info("STEPS: A {:.1f}; B {:.1f}; R {:.1f}; T {:.1f}. \t",
                                     sqrtf(sumA) / (0.0005 * setting_thOptIterations),
                                     sqrtf(sumB) / (0.00005 * setting_thOptIterations),
                                     sqrtf(sumR) / (0.00005 * setting_thOptIterations),
                                     sqrtf(sumT) * sumNID / (0.00005 * setting_thOptIterations));


    EFDeltaValid = false;
    setPrecalcValues();


    return sqrtf(sumA) < 0.0005*setting_thOptIterations &&
           sqrtf(sumB) < 0.00005*setting_thOptIterations &&
           sqrtf(sumR) < 0.00005*setting_thOptIterations &&
           sqrtf(sumT)*sumNID < 0.00005*setting_thOptIterations;
//	printf("mean steps: %f %f %f!\n",
//			meanStepC, meanStepP, meanStepD);
}

// sets linearization point.
void Tracker::backupState(bool backupLastStep) {
    hCalib.value_backup = hCalib.value;
    for(auto& fh : frameHessians)
    {
        fh->state_backup = fh->get_state();
        for(PointHessian* ph : fh->pointHessians)
            ph->idepth_backup = ph->idepth;
    }
}

// sets linearization point.
void Tracker::loadSateBackup() {
    hCalib.setValue(hCalib.value_backup);
    for(auto& fh : frameHessians)
    {
        fh->setState(fh->state_backup);
        for(PointHessian* ph : fh->pointHessians)
        {
            ph->setIdepth(ph->idepth_backup);

            ph->setIdepthZero(ph->idepth_backup);
        }

    }


    EFDeltaValid=false;
    setPrecalcValues();
}

void
Tracker::linearizeAll_Reductor(bool fixLinearization, std::vector<PointFrameResidual *> *toRemove, int min, int max,
                               Vec10 *stats, int tid) {

    for(int k=min;k<max;k++)
    {
        PointFrameResidual* r = activeResiduals[k];
        (*stats)[0] += r->linearize(&hCalib);

        if(fixLinearization)
        {
            r->applyRes(true);

            if(r->efResidual->isActive())
            {
                if(r->isNew)
                {
                    PointHessian* p = r->point;
                    Vec3f ptp_inf = r->host->targetPrecalc[r->target->trackingID].PRE_KRKiTll * Vec3f(p->u,p->v, 1);	// projected point assuming infinite depth.
                    Vec3f ptp = ptp_inf + r->host->targetPrecalc[r->target->trackingID].PRE_KtTll*p->idepth_scaled;	// projected point with real depth.
                    float relBS = 0.01*((ptp_inf.head<2>() / ptp_inf[2])-(ptp.head<2>() / ptp[2])).norm();	// 0.01 = one pixel.


                    if(relBS > p->maxRelBaseline)
                        p->maxRelBaseline = relBS;

                    p->numGoodResiduals++;
                }
            }
            else
            {
                toRemove[tid].push_back(activeResiduals[k]);
            }
        }
    }
}

void Tracker::setNewFrameEnergyTH() {

    // collect all residuals and make decision on TH.
    allResVec.clear();
    allResVec.reserve(activeResiduals.size() * 2);
    auto newFrame = frameHessians.back();

    for (PointFrameResidual *r: activeResiduals)
        if (r->state_NewEnergyWithOutlier >= 0 && r->target->trackingID == newFrame->trackingID) {
            allResVec.push_back(r->state_NewEnergyWithOutlier);

        }

    if (allResVec.size() == 0) {
        newFrame->frameEnergyTH = 12 * 12 * patternNum;
        Log::Logger::getInstance()->info("allResVec should not be empty but it was anyways...");
        return;        // should never happen, but lets make sure.
    }


    int nthIdx = setting_frameEnergyTHN * allResVec.size();

    assert(nthIdx < (int) allResVec.size());
    assert(setting_frameEnergyTHN < 1);

    std::nth_element(allResVec.begin(), allResVec.begin() + nthIdx, allResVec.end());
    float nthElement = sqrtf(allResVec[nthIdx]);


    newFrame->frameEnergyTH = nthElement*setting_frameEnergyTHFacMedian;
    newFrame->frameEnergyTH = 26.0f*setting_frameEnergyTHConstWeight + newFrame->frameEnergyTH*(1-setting_frameEnergyTHConstWeight);
    newFrame->frameEnergyTH = newFrame->frameEnergyTH*newFrame->frameEnergyTH;
    newFrame->frameEnergyTH *= setting_overallEnergyTHWeight*setting_overallEnergyTHWeight;

    Log::Logger::getInstance()->info("nthIndex and nthElement {}, {}, new frame threshold: {}", nthIdx, nthElement, newFrame->frameEnergyTH);
}

Vec3 Tracker::linearizeAll(bool fixLinearization) {
    double lastEnergyP = 0;
    double lastEnergyR = 0;
    double num = 0;

    std::vector<PointFrameResidual *> toRemove[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) toRemove[i].clear();

    Vec10 stats;
    stats.setZero();
    linearizeAll_Reductor(fixLinearization, toRemove, 0, activeResiduals.size(), &stats, 0);
    lastEnergyP = stats[0];

    setNewFrameEnergyTH();

    if (fixLinearization) {


        for (PointFrameResidual *r: activeResiduals) {
            PointHessian* ph = r->point;
            if(ph->lastResiduals[0].first == r)
                ph->lastResiduals[0].second = r->state_state;
            else if(ph->lastResiduals[1].first == r)
                ph->lastResiduals[1].second = r->state_state;


        }

        int nResRemoved = 0;
        for (int i = 0; i < NUM_THREADS; i++) {
            for (PointFrameResidual *r: toRemove[i]) {
                PointHessian* ph = r->point;

                if(ph->lastResiduals[0].first == r)
                    ph->lastResiduals[0].first=0;
                else if(ph->lastResiduals[1].first == r)
                    ph->lastResiduals[1].first=0;

                for(unsigned int k=0; k<ph->residuals.size();k++)
                    if(ph->residuals[k] == r)
                    {
                        ef.dropResidual(r->efResidual);
                        deleteOut<PointFrameResidual>(ph->residuals,k);
                        nResRemoved++;
                        break;
                    }
            }
        }
        Log::Logger::getInstance()->info("FINAL LINEARIZATION: removed {} / {} residuals!", nResRemoved,
                                         (int) activeResiduals.size());
    }

    return Vec3(lastEnergyP, lastEnergyR, num);
}

void Tracker::removeOutliers() {
    int numPointsDropped=0;
    for(auto& fh : frameHessians)
    {
        for(unsigned int i=0;i<fh->pointHessians.size();i++)
        {
            PointHessian* ph = fh->pointHessians[i];
            if(ph==0) {
                std::cout << "Poiint hessian " << std::endl;
                continue;

            }

            if(ph->residuals.size() == 0)
            {
                fh->pointHessiansOut.push_back(ph);
                ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                fh->pointHessians[i] = fh->pointHessians.back();
                fh->pointHessians.pop_back();
                i--;
                numPointsDropped++;
            }
        }
    }
    ef.dropPointsF();
}

