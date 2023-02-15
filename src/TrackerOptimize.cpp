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
    int idx = 0;

    for (auto &frameHessian: frameHessians) {
        auto *fh = frameHessian.get();

        for (auto &pointHessian: fh->pointHessians) {
            auto *ph = &pointHessian;

            for (int i = 0; i < ph->residuals.size(); ++i) {
                auto *r = &ph->residuals[i];

                for (int m = 0; m < ef.eFrames[fh->trackingID].points[ph->pointIndex].residualsAll.size(); ++m) {
                    auto *res = &ef.eFrames[fh->trackingID].points[ph->pointIndex].residualsAll[m];
                    if (ph->pointIndex == r->pointIndex && !res->isLinearized) {
                        activeResiduals.push_back(r);
                    } else
                        numLRes++;
                }

            }

        }
        numPoints++;
    }


    Log::Logger::getInstance()->info("OPTIMIZE {} pts, {} active res, {} lin res!", ef.nPoints,
                                     activeResiduals.size(), numLRes);


    Vec3 lastEnergy = linearizeAll(false);
    Log::Logger::getInstance()->info("LastEnergy: {}", lastEnergy[0]);
    double lastEnergyL = ef.calcLEnergyF_MT();
    double lastEnergyM = ef.calcMEnergyF();
    Log::Logger::getInstance()->info("Last Energy L: {} and M: {}", lastEnergyL, lastEnergyM);
/*

applyRes_Reductor(true,0,activeResiduals.size(),0,0);


if(!setting_debugout_runquiet)
{
    printf("Initial Error       \t");
    printOptRes(lastEnergy, lastEnergyL, lastEnergyM, 0, 0, frameHessians.back()->aff_g2l().a, frameHessians.back()->aff_g2l().b);
}

debugPlotTracking();



double lambda = 1e-1;
float stepsize=1;
VecX previousX = VecX::Constant(CPARS+ 8*frameHessians.size(), NAN);
for(int iteration=0;iteration<mnumOptIts;iteration++)
{
    // solve!
    backupState(iteration!=0);
    //solveSystemNew(0);
    solveSystem(iteration, lambda);
    double incDirChange = (1e-20 + previousX.dot(ef->lastX)) / (1e-20 + previousX.norm() * ef->lastX.norm());
    previousX = ef->lastX;


    if(std::isfinite(incDirChange) && (setting_solverMode & SOLVER_STEPMOMENTUM))
    {
        float newStepsize = exp(incDirChange*1.4);
        if(incDirChange<0 && stepsize>1) stepsize=1;

        stepsize = sqrtf(sqrtf(newStepsize*stepsize*stepsize*stepsize));
        if(stepsize > 2) stepsize=2;
        if(stepsize <0.25) stepsize=0.25;
    }

    bool canbreak = doStepFromBackup(stepsize,stepsize,stepsize,stepsize,stepsize);







    // eval new energy!
    Vec3 newEnergy = linearizeAll(false);
    double newEnergyL = calcLEnergy();
    double newEnergyM = calcMEnergy();




    if(!setting_debugout_runquiet)
    {
        printf("%s %d (L %.2f, dir %.2f, ss %.1f): \t",
               (newEnergy[0] +  newEnergy[1] +  newEnergyL + newEnergyM <
                lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM) ? "ACCEPT" : "REJECT",
               iteration,
               log10(lambda),
               incDirChange,
               stepsize);
        printOptRes(newEnergy, newEnergyL, newEnergyM , 0, 0, frameHessians.back()->aff_g2l().a, frameHessians.back()->aff_g2l().b);
    }

    if(setting_forceAceptStep || (newEnergy[0] +  newEnergy[1] +  newEnergyL + newEnergyM <
                                  lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM))
    {

        if(multiThreading)
            treadReduce.reduce(boost::bind(&FullSystem::applyRes_Reductor, this, true, _1, _2, _3, _4), 0, activeResiduals.size(), 50);
        else
            applyRes_Reductor(true,0,activeResiduals.size(),0,0);

        lastEnergy = newEnergy;
        lastEnergyL = newEnergyL;
        lastEnergyM = newEnergyM;

        lambda *= 0.25;
    }
    else
    {
        loadSateBackup();
        lastEnergy = linearizeAll(false);
        lastEnergyL = calcLEnergy();
        lastEnergyM = calcMEnergy();
        lambda *= 1e2;
    }


    if(canbreak && iteration >= setting_minOptIterations) break;
}



Vec10 newStateZero = Vec10::Zero();
newStateZero.segment<2>(6) = frameHessians.back()->get_state().segment<2>(6);

frameHessians.back()->setEvalPT(frameHessians.back()->PRE_worldToCam,
                                newStateZero);
EFDeltaValid=false;
EFAdjointsValid=false;
ef->setAdjointsF(&Hcalib);
setPrecalcValues();




lastEnergy = linearizeAll(true);




if(!std::isfinite((double)lastEnergy[0]) || !std::isfinite((double)lastEnergy[1]) || !std::isfinite((double)lastEnergy[2]))
{
    printf("KF Tracking failed: LOST!\n");
    isLost=true;
}


statistics_lastFineTrackRMSE = sqrtf((float)(lastEnergy[0] / (patternNum*ef->resInA)));

if(calibLog != 0)
{
    (*calibLog) << Hcalib.value_scaled.transpose() <<
                " " << frameHessians.back()->get_state_scaled().transpose() <<
                " " << sqrtf((float)(lastEnergy[0] / (patternNum*ef->resInA))) <<
                " " << ef->resInM << "\n";
    calibLog->flush();
}

{
    boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
    for(FrameHessian* fh : frameHessians)
    {
        fh->shell->camToWorld = fh->PRE_camToWorld;
        fh->shell->aff_g2l = fh->aff_g2l();
    }
}




debugPlotTracking();

return sqrtf((float)(lastEnergy[0] / (patternNum*ef->resInA)));
*/
    return 12;
}


void
Tracker::linearizeAll_Reductor(bool fixLinearization, std::vector<PointFrameResidual *> *toRemove, int min, int max,
                               Vec10 *stats, int tid) {
    for (int k = min; k < max; k++) {
        PointFrameResidual *r = activeResiduals[k];
        PointHessian *ph = &frameHessians[frameHessians.size() - 2]->pointHessians[k];
        auto host = frameHessians[frameHessians.size() - 2];
        auto target = frameHessians[frameHessians.size() - 2];
        (*stats)[0] += linearize(ph, &hCalib, r, host, target);
        EFResidual *efResidual = r->efResidual;

        if (fixLinearization) {
            VO::applyResidual(true, efResidual, r);

            if (efResidual->isActive()) {
                if (r->isNew) {
                    PointHessian *p = ph;
                    Vec3f ptp_inf = host->targetPrecalc[target->trackingID].PRE_KRKiTll *
                                    Vec3f(p->u, p->v, 1);    // projected point assuming infinite depth.
                    Vec3f ptp = ptp_inf + host->targetPrecalc[target->trackingID].PRE_KtTll *
                                          p->idepth_scaled;    // projected point with real depth.
                    float relBS = 0.01 * ((ptp_inf.head<2>() / ptp_inf[2]) -
                                          (ptp.head<2>() / ptp[2])).norm();    // 0.01 = one pixel.


                    if (relBS > p->maxRelBaseline)
                        p->maxRelBaseline = relBS;

                    p->numGoodResiduals++;
                }
            } else {
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
        if (r->state_NewEnergyWithOutlier >= 0 && r->trackingID == newFrame->trackingID) {
            allResVec.push_back(r->state_NewEnergyWithOutlier);

        }

    if (allResVec.size() == 0) {
        newFrame->frameEnergyTH = 12 * 12 * patternNum;
        return;        // should never happen, but lets make sure.
    }


    int nthIdx = setting_frameEnergyTHN * allResVec.size();

    assert(nthIdx < (int) allResVec.size());
    assert(setting_frameEnergyTHN < 1);

    std::nth_element(allResVec.begin(), allResVec.begin() + nthIdx, allResVec.end());
    float nthElement = sqrtf(allResVec[nthIdx]);


    newFrame->frameEnergyTH = nthElement * setting_frameEnergyTHFacMedian;
    newFrame->frameEnergyTH =
            26.0f * setting_frameEnergyTHConstWeight + newFrame->frameEnergyTH * (1 - setting_frameEnergyTHConstWeight);
    newFrame->frameEnergyTH = newFrame->frameEnergyTH * newFrame->frameEnergyTH;
    newFrame->frameEnergyTH *= setting_overallEnergyTHWeight * setting_overallEnergyTHWeight;



//
//	int good=0,bad=0;
//	for(float f : allResVec) if(f<newFrame->frameEnergyTH) good++; else bad++;
//	printf("EnergyTH: mean %f, median %f, result %f (in %d, out %d)! \n",
//			meanElement, nthElement, sqrtf(newFrame->frameEnergyTH),
//			good, bad);

    Log::Logger::getInstance()->info("Set new frame threshold {}", newFrame->frameEnergyTH);
}

Vec3 Tracker::linearizeAll(bool fixLinearization) {
    double lastEnergyP = 0;
    double lastEnergyR = 0;
    double num = 0;


    std::vector<PointFrameResidual *> toRemove[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) toRemove[i].clear();


    Vec10 stats;
    linearizeAll_Reductor(fixLinearization, toRemove, 0, activeResiduals.size(), &stats, 0);
    lastEnergyP = stats[0];


    setNewFrameEnergyTH();


    if (fixLinearization) {


        for (PointFrameResidual *r: activeResiduals) {
            PointHessian *ph = r->pointH;
            if (ph->lastResiduals[0].first == r)
                ph->lastResiduals[0].second = r->state_state;
            else if (ph->lastResiduals[1].first == r)
                ph->lastResiduals[1].second = r->state_state;


        }

        int nResRemoved = 0;
        for (int i = 0; i < NUM_THREADS; i++) {
            for (PointFrameResidual *r: toRemove[i]) {
                PointHessian *ph = r->pointH;

                if (ph->lastResiduals[0].first == r)
                    ph->lastResiduals[0].first = 0;
                else if (ph->lastResiduals[1].first == r)
                    ph->lastResiduals[1].first = 0;

                for (unsigned int k = 0; k < ph->residuals.size(); k++)
                    if (ph->residuals[k].pointIndex == r->pointIndex) {
                        ef.dropResidual(r);

                        auto it = ph->residuals.begin();
                        std::advance(it, k);
                        ph->residuals.erase(it);

                        nResRemoved++;
                        break;
                    }
            }
        }
        Log::Logger::getInstance()->info("FINAL LINEARIZATION: removed %d / %d residuals!\n", nResRemoved,
                                         (int) activeResiduals.size());
    }

    return Vec3(lastEnergyP, lastEnergyR, num);
}

