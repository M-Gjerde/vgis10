//
// Created by magnus on 2/7/23.
//

#include "Tracker.h"

#include "TrackerInitializer.h"
#include "PointHessian.h"
#include "Optimized/PointFrameResidual.h"
#include <sophus/so3.hpp>
#include "sophus/se3.hpp"

void Tracker::initializeFromInitializer() {
    firstFrame->trackingID = 0;
    ef.insertFrame(firstFrame, &hCalib);
    firstFrame->frameID = 0;
    firstFrame->trackingID = 0;

    // RescaleFactor
    float sumID = 1e-5, numID = 1e-5;
    for (int i = 0; i < points.numPointsLevel[0]; ++i) {
        sumID += points.pointsLevel[0][i].iR;
        numID++;
    }
    float rescaleFactor = 1 / (sumID / numID);
    float keepPercentage = setting_desiredPointDensity / (float) points.numPointsLevel[0];
    Log::Logger::getInstance()->info("Initialization: keep {:.2f}% (need {}, have {})!", 100 * keepPercentage,
                                     (int) (setting_desiredPointDensity), points.numPointsLevel[0]);

    int ptIdx = 0;
    for (int i = 0; i < points.numPointsLevel[0]; i++) {
        // Only run this full loop keepPercentage times of number of points at level 0 but choose to do so randomly
        if (rand() / (float) RAND_MAX > keepPercentage) continue;
        Pnt *point = points.pointsLevel[0].data() + i;

        ImmaturePoint *pt = new ImmaturePoint(point->u + 0.5f, point->v + 0.5f, point->my_type, firstFrame);

        if (!std::isfinite(pt->energyTH)) continue;
        pt->idepth_max = pt->idepth_min = 1;

        PointHessian *ph = new PointHessian(pt, &hCalib);
        //VO::pointHessianFromImmaturePoint(&pt, ph);k
        if (!std::isfinite(ph->energyTH)) continue;
        ph->setIdepthScaled(point->iR * rescaleFactor);
        ph->setIdepthZero(ph->idepth);
        ph->hasDepthPrior = true;
        ph->setPointStatus(PointHessian::ACTIVE);
        ph->pointIndex = ptIdx;
        ptIdx++;
        firstFrame->pointHessians.push_back(ph);
        ef.insertPoint(firstFrame->pointHessians.back());
    }

    SE3 firstToNew = info.thisToNext;
    firstToNew.translation() /= rescaleFactor;
    frameHessians.emplace_back(firstFrame);

    firstFrame->shell->camToWorld = SE3();
    firstFrame->shell->aff_g2l = AffLight(0,0);
    firstFrame->setEvalPT_scaled(firstFrame->shell->camToWorld.inverse(),firstFrame->shell->aff_g2l);
    firstFrame->shell->trackingRef=0;
    firstFrame->shell->camToTrackingRef = SE3();


    secondFrame->shell->camToWorld = firstToNew.inverse();
    secondFrame->shell->aff_g2l = AffLight(0,0);
    secondFrame->setEvalPT_scaled(secondFrame->shell->camToWorld.inverse(),secondFrame->shell->aff_g2l);
    secondFrame->shell->trackingRef = firstFrame->shell;
    secondFrame->shell->camToTrackingRef = firstToNew.inverse();
    setPrecalcValues();

    initialized = true;

    Log::Logger::getInstance()->info("INITIALIZE FROM INITIALIZER ({} pts)!", (int) firstFrame->pointHessians.size());
}


bool Tracker::initializerTrackFrame(const std::shared_ptr<VO::Frame> &frame, const CameraCalibration *calibration) {
    int maxIterations[] = {5, 5, 10, 30, 50};
    secondFrame = frame;

    info.JbBuffer.resize(frame->width * frame->height);
    info.JbBuffer_new.resize(frame->width * frame->height);
    info.wM.diagonal()[0] = info.wM.diagonal()[1] = info.wM.diagonal()[2] = SCALE_XI_ROT;
    info.wM.diagonal()[3] = info.wM.diagonal()[4] = info.wM.diagonal()[5] = SCALE_XI_TRANS;
    info.wM.diagonal()[6] = SCALE_A;
    info.wM.diagonal()[7] = SCALE_B;

    if (!info.snapped) {
        info.thisToNext.translation().setZero();
        for (int lvl = 0; lvl < calibration->pyrLevelsUsed; lvl++) {
            int npts = points.numPointsLevel[lvl];
            Pnt *ptsl = points.pointsLevel[lvl].data();
            for (int i = 0; i < npts; i++) {
                ptsl[i].iR = 1;
                ptsl[i].idepth_new = 1;
                ptsl[i].lastHessian = 0;
            }
        }
    }

    SE3 refToNew_current = info.thisToNext;
    AffLight refToNew_aff_current = info.thisToNext_aff;

    if (firstFrame->abExposure > 0 && frame->abExposure > 0)
        refToNew_aff_current = AffLight(logf(frame->abExposure / firstFrame->abExposure),
                                        0); // coarse approximation.

    Vec3f latestRes = Vec3f::Zero();
    for (int lvl = calibration->pyrLevelsUsed - 1; lvl >= 0; lvl--) {

        Mat88f H, Hsc;
        Vec8f b, bsc;
        Vec3f resOld = VO::calcResAndGS(lvl, H, b, Hsc, bsc, refToNew_current,
                                        refToNew_aff_current, false, firstFrame,
                                        frame, &info, &points);

        if (H.hasNaN()) {
            std::cout << H << std::endl;
            throw std::runtime_error("H has NAN");
        }
        if (b.hasNaN()) {
            std::cout << b << std::endl;
            throw std::runtime_error("b has NAN");
        }


        float lambda = 0.1;
        float eps = 1e-4;
        int fails = 0;


        Log::Logger::getInstance()->info(
                "lvl {}, it {} (l={}) {}: {:.3f}+{:.5f} -> {:.3f}+{:.5f} ({:.3f}->{:.3f}) (|inc| = {})!",
                lvl,
                0,
                lambda,
                "INITIA",
                sqrtf((float) (resOld[0] / resOld[2])),
                sqrtf((float) (resOld[1] / resOld[2])),
                sqrtf((float) (resOld[0] / resOld[2])),
                sqrtf((float) (resOld[1] / resOld[2])),
                (resOld[0] + resOld[1]) / resOld[2],
                (resOld[0] + resOld[1]) / resOld[2],
                0.0f);
        //std::cout << refToNew_current.log().transpose() << " AFF " << refToNew_aff_current.vec().transpose()<< "\n";


        int iteration = 0;
        while (true) {
            Mat88f Hl = H;
            for (int i = 0; i < 8; i++) Hl(i, i) *= (1 + lambda);
            Hl -= Hsc * (1 / (1 + lambda));
            Vec8f bl = b - bsc * (1 / (1 + lambda));

            Hl = info.wM * Hl * info.wM * (0.01f / (frame->widthLevel[lvl] * frame->heightLevel[lvl]));
            bl = info.wM * bl * (0.01f / (frame->widthLevel[lvl] * frame->heightLevel[lvl]));

            Vec8f inc;
            inc.head<6>() = -(info.wM.toDenseMatrix().topLeftCorner<6, 6>() *
                              (Hl.topLeftCorner<6, 6>().ldlt().solve(bl.head<6>())));
            inc.tail<2>().setZero();

            SE3 refToNew_new = SE3::exp(inc.head<6>().cast<double>()) * refToNew_current;
            AffLight refToNew_aff_new = refToNew_aff_current;
            refToNew_aff_new.a += inc[6];
            refToNew_aff_new.b += inc[7];

            Mat88f H_new, Hsc_new;
            Vec8f b_new, bsc_new;
            Vec3f resNew = VO::calcResAndGS(lvl, H_new, b_new, Hsc_new, bsc_new, refToNew_new, refToNew_aff_new, false,
                                            firstFrame,
                                            frame, &info, &points);

            if (H_new.hasNaN()) {
                std::cout << H_new << std::endl;
                throw std::runtime_error("H_new has NAN");
            }
            if (b_new.hasNaN()) {
                std::cout << b_new << std::endl;
                throw std::runtime_error("b_new has NAN");
            }

            Vec3f regEnergy = VO::calcEC(lvl, &points, &info);

            float eTotalNew = (resNew[0] + resNew[1] + regEnergy[1]);
            float eTotalOld = (resOld[0] + resOld[1] + regEnergy[0]);

            bool accept = eTotalOld > eTotalNew;

            Log::Logger::getInstance()->info(
                    "lvl {}, it {} (l={}) {}: {:.5f} + {:.5f} + {:.5f} -> {:.5f} + {:.5f} + {:.5f} ({:.2f}->{:.2f}) (|inc| = {})! \t",
                    lvl, iteration, lambda,
                    (accept ? "ACCEPT" : "REJECT"),
                    sqrtf((float) (resOld[0] / resOld[2])),
                    sqrtf((float) (regEnergy[0] / regEnergy[2])),
                    sqrtf((float) (resOld[1] / resOld[2])),
                    sqrtf((float) (resNew[0] / resNew[2])),
                    sqrtf((float) (regEnergy[1] / regEnergy[2])),
                    sqrtf((float) (resNew[1] / resNew[2])),
                    eTotalOld / resNew[2],
                    eTotalNew / resNew[2],
                    inc.norm());
            //std::cout << refToNew_new.log().transpose() << " AFF " << refToNew_aff_new.vec().transpose() << "\n";

            if (accept) {

                if (resNew[1] == info.alphaK * points.numPointsLevel[lvl])
                    info.snapped = true;
                H = H_new;
                b = b_new;
                Hsc = Hsc_new;
                bsc = bsc_new;
                resOld = resNew;
                refToNew_aff_current = refToNew_aff_new;
                refToNew_current = refToNew_new;
                //applyStep(lvl);
                //optReg(lvl);
                lambda *= 0.5;
                fails = 0;
                if (lambda < 0.0001) lambda = 0.0001;
            } else {
                fails++;
                lambda *= 4;
                if (lambda > 10000) lambda = 10000;
            }
            bool quitOpt = false;

            if (!(inc.norm() > eps) || iteration >= maxIterations[lvl] || fails >= 2) {
                Mat88f H, Hsc;
                Vec8f b, bsc;

                quitOpt = true;
            }


            if (quitOpt) break;
            iteration++;
        }
        latestRes = resOld;
    };

    info.thisToNext = refToNew_current;
    info.thisToNext_aff = refToNew_aff_current;

    info.frameID++;
    if (!info.snapped) info.snappedAtFrame = 0;

    if (info.snapped && info.snappedAtFrame == 0)
        info.snappedAtFrame = info.frameID;

    Log::Logger::getInstance()->info("Snapped: {}, FrameID:  {}, snappedAtFrame {}", info.snapped, info.frameID,
                                     info.snappedAtFrame);
    return info.snapped && info.frameID > info.snappedAtFrame + 5;
}

void Tracker::takeTrackedFrame(const std::shared_ptr<VO::Frame> &frame, bool needKF) {
    if (needKF) {
        makeKeyFrame(frame);
    } else {
        makeNonKeyFrame(frame);
    }
}

void Tracker::makeKeyFrame(std::shared_ptr<VO::Frame> frame) {
    // =========================== UPDATE POSE =========================
    assert(fh->shell->trackingRef != 0);

    frame->shell->camToWorld =  frame->shell->trackingRef->camToWorld * frame->shell->camToTrackingRef;
    frame->setEvalPT_scaled(frame->shell->camToWorld.inverse(),frame->shell->aff_g2l);
    // traceNewCoarse
    traceNewCoarse(frame);
    // =========================== Flag Frames to be Marginalized. =========================
    flagFramesForMarginalization(frame);
    // =========================== add New Frame to Hessian Struct. =========================


    frame->trackingID = frameHessians.size();
    frameHessians.emplace_back(frame);
    frame->frameID = allKeyFramesHistory.size();
    allKeyFramesHistory.emplace_back(frame->shell);
    ef.insertFrame(frame, &hCalib);
    setPrecalcValues();

    // =========================== add new residuals for old points =========================
    for (auto &fh1: frameHessians) {
        if (fh1 == frame) continue;
        for (PointHessian *ph: fh1->pointHessians) {
            auto *r = new PointFrameResidual(ph, fh1.get(), frame.get());
            r->setState(ResState::IN);
            ph->residuals.push_back(r);
            ef.insertResidual(r);
            ph->lastResiduals[1] = ph->lastResiduals[0];
            ph->lastResiduals[0] = std::pair<PointFrameResidual *, ResState>(r, ResState::IN);
        }
    }


    // =========================== Activate Points (& flag for marginalization). =========================
    activatePointsMT();
    ef.makeIDX();

    // =========================== OPTIMIZE ALL =========================
    frame->frameEnergyTH = frameHessians.back()->frameEnergyTH;
    float rmse = optimize(setting_maxOptIterations);
    Log::Logger::getInstance()->info("RMSE: {}", rmse);

    // =========================== Figure Out if INITIALIZATION FAILED =========================
    if (allKeyFramesHistory.size() <= 4) {
        if (allKeyFramesHistory.size() == 2 && rmse > 20 * benchmark_initializerSlackFactor) {
            Log::Logger::getInstance()->info("I THINK INITIALIZATION FAILED! Resetting");
            initialized = false;
            return;
        }
        if (allKeyFramesHistory.size() == 3 && rmse > 13 * benchmark_initializerSlackFactor) {
            Log::Logger::getInstance()->info("I THINK INITIALIZATION FAILED! Resetting");
            initialized = false;
            return;
        }
        if (allKeyFramesHistory.size() == 4 && rmse > 9 * benchmark_initializerSlackFactor) {
            Log::Logger::getInstance()->info("I THINK INITIALIZATION FAILED! Resetting");
            initialized = false;
            return;
        }
    }

    // =========================== REMOVE OUTLIER =========================
    removeOutliers();

    coarseTracker_forNewKF->makeK(&hCalib, calibration);
    coarseTracker_forNewKF->setCoarseTrackingRef(frameHessians);
    // =========================== (Activate-)Marginalize Points =========================
    flagPointsForRemoval();
    ef.dropPointsF();
    getNullspaces(
            ef.lastNullspaces_pose,
            ef.lastNullspaces_scale,
            ef.lastNullspaces_affA,
            ef.lastNullspaces_affB);
    ef.marginalizePointsF();

    // =========================== add new Immature points & new residuals =========================
    Log::Logger::getInstance()->info("Making new traces on frame obj: {}", reinterpret_cast<void *>(frame.get()));
    makeNewTraces(frame, 0);

    // =========================== Marginalize Frames =========================
    for (unsigned int i = 0; i < frameHessians.size(); i++)
        if (frameHessians[i]->flaggedForMarginalization) {
            marginalizeFrame(frameHessians[i].get());
            i = 0;
        }
}


void Tracker::makeNonKeyFrame(std::shared_ptr<VO::Frame> frame) {
    // =========================== UPDATE POSE =========================
    // needs to be set by mapping thread. no lock required since we are in mapping thread.
    {
        assert(fh->shell->trackingRef != 0);
        frame->shell->camToWorld =  frame->shell->trackingRef->camToWorld * frame->shell->camToTrackingRef;
        frame->setEvalPT_scaled(frame->shell->camToWorld.inverse(),frame->shell->aff_g2l);
    }

    traceNewCoarse(frame);
}


void Tracker::marginalizeFrame(VO::Frame *frame) {
    // marginalize or remove all this frames points.

    assert((int) frame->pointHessians.size() == 0);


    ef.marginalizeFrame(frame->efFrame);

    // drop all observations of existing points in that frame.

    for (auto &fh: frameHessians) {
        if (fh.get() == frame) continue;

        for (PointHessian *ph: fh->pointHessians) {
            for (unsigned int i = 0; i < ph->residuals.size(); i++) {
                PointFrameResidual *r = ph->residuals[i];
                if (r->target == frame) {
                    if (ph->lastResiduals[0].first == r)
                        ph->lastResiduals[0].first = 0;
                    else if (ph->lastResiduals[1].first == r)
                        ph->lastResiduals[1].first = 0;


                    ef.dropResidual(r->efResidual);
                    deleteOut<PointFrameResidual>(ph->residuals, i);
                    break;
                }
            }
        }
    }


    std::vector<VO::Frame *> v;
    v.push_back(frame);
    {

        // for(IOWrap::Output3DWrapper* ow : outputWrapper)
        //     ow->publishKeyframes(v, true, &Hcalib);
    }


    frame->shell->marginalizedAt = frameHessians.back()->shell->id;
    frame->shell->movedByOpt = frame->w2c_leftEps().norm();

    for(int i = 0; auto& f : frameHessians){
        if (f.get() == frame){
            frameHessians.erase(frameHessians.begin() + i);
        }
        i++;
    }
    for (unsigned int i = 0; i < frameHessians.size(); i++)
        frameHessians[i]->trackingID = i;


    setPrecalcValues();
    ef.setAdjointsF(&hCalib);
}

void Tracker::makeNewTraces(std::shared_ptr<VO::Frame> newFrame, float *gtDepth) {
    //int numPointsTotal = makePixelStatus(newFrame->dI, selectionMap, wG[0], hG[0], setting_desiredDensity);
    int numPointsTotal = pixelSelector->makeMaps(newFrame.get(), selectionMap, setting_desiredImmatureDensity);

    newFrame->pointHessians.reserve(numPointsTotal * 1.2f);
    //fh->pointHessiansInactive.reserve(numPointsTotal*1.2f);
    newFrame->pointHessiansMarginalized.reserve(numPointsTotal * 1.2f);
    newFrame->pointHessiansOut.reserve(numPointsTotal * 1.2f);

    std::shared_ptr <VO::Frame> newFrameShared(newFrame);

    for (int y = patternPadding + 1; y < calibration->hG[0] - patternPadding - 2; y++)
        for (int x = patternPadding + 1; x < calibration->wG[0] - patternPadding - 2; x++) {
            int i = x + y * calibration->wG[0];
            if (selectionMap[i] == 0) continue;

            ImmaturePoint *impt = new ImmaturePoint(x, y, selectionMap[i], newFrameShared);
            if (!std::isfinite(impt->energyTH)) delete impt;
            else newFrame->immaturePoints.push_back(impt);

        }
    //printf("MADE %d IMMATURE POINTS!\n", (int)newFrame->immaturePoints.size());

}


Vec4 Tracker::trackNewCoarse(std::shared_ptr<VO::Frame> fh) {

    assert(allFrameHistory.size() > 0);
    // set pose initialization.

    //for(IOWrap::Output3DWrapper* ow : outputWrapper)
    //    ow->pushLiveFrame(fh);



    VO::Frame *lastF = coarseTracker->lastRef.get();

    AffLight aff_last_2_l = AffLight(0, 0);

    std::vector<SE3, Eigen::aligned_allocator<SE3>> lastF_2_fh_tries;
    if (allFrameHistory.size() == 2)
        for (unsigned int i = 0; i < lastF_2_fh_tries.size(); i++) lastF_2_fh_tries.push_back(SE3());
    else {
        auto &slast = allFrameHistory[allFrameHistory.size() - 2];
        auto &sprelast = allFrameHistory[allFrameHistory.size() - 3];
        SE3 slast_2_sprelast;
        SE3 lastF_2_slast;
        {    // lock on global pose consistency!
            slast_2_sprelast = sprelast->camToWorld.inverse() * slast->camToWorld;
            lastF_2_slast = slast->camToWorld.inverse() * lastF->shell->camToWorld;
            aff_last_2_l = slast->aff_g2l;
        }
        SE3 fh_2_slast = slast_2_sprelast;// assumed to be the same as fh_2_slast.


        // get last delta-movement.
        lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast);    // assume constant motion.
        lastF_2_fh_tries.push_back(
                fh_2_slast.inverse() * fh_2_slast.inverse() * lastF_2_slast);    // assume double motion (frame skipped)
        lastF_2_fh_tries.push_back(SE3::exp(fh_2_slast.log() * 0.5).inverse() * lastF_2_slast); // assume half motion.
        lastF_2_fh_tries.push_back(lastF_2_slast); // assume zero motion.
        lastF_2_fh_tries.push_back(SE3()); // assume zero motion FROM KF.


        // just try a TON of different initializations (all rotations). In the end,
        // if they don't work they will only be tried on the coarsest level, which is super fast anyway.
        // also, if tracking rails here we loose, so we really, really want to avoid that.
        for (float rotDelta = 0.02; rotDelta < 0.05; rotDelta++) {
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, 0, 0),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, rotDelta, 0),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, 0, rotDelta),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, 0, 0),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, -rotDelta, 0),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, 0, -rotDelta),
                                           Vec3(0, 0, 0)));            // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, rotDelta, 0),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, 0, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, rotDelta, 0),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, -rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, 0, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, -rotDelta, 0),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, 0, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, -rotDelta, 0),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, 0, -rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, 0, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, -rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, -rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, -rotDelta, rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, -rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, -rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, rotDelta, -rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast *
                                       SE3(Eigen::Quaternion<double>(1, rotDelta, rotDelta, rotDelta),
                                           Vec3(0, 0, 0)));    // assume constant motion.
        }

        if (!slast->poseValid || !sprelast->poseValid || !lastF->shell->poseValid) {
            lastF_2_fh_tries.clear();
            lastF_2_fh_tries.push_back(SE3());
        }
    }


    Vec3 flowVecs = Vec3(100, 100, 100);
    SE3 lastF_2_fh = SE3();
    AffLight aff_g2l = AffLight(0, 0);


    // as long as maxResForImmediateAccept is not reached, I'll continue through the options.
    // I'll keep track of the so-far best achieved residual for each level in achievedRes.
    // If on a coarse level, tracking is WORSE than achievedRes, we will not continue to save time.


    Vec5 achievedRes = Vec5::Constant(NAN);
    bool haveOneGood = false;
    bool trackingIsGood = false;
    int tryIterations = 0;
    for (unsigned int i = 0; i < lastF_2_fh_tries.size(); i++) {
        AffLight aff_g2l_this = aff_last_2_l;
        SE3 lastF_2_fh_this = lastF_2_fh_tries[i];
        trackingIsGood = coarseTracker->trackNewestCoarse(
                fh, lastF_2_fh_this, aff_g2l_this,
                calibration->pyrLevelsUsed - 1,
                achievedRes);    // in each level has to be at least as good as the last try.
        tryIterations++;

        if (i != 0) {
            printf("RE-TRACK ATTEMPT %d with initOption %d and start-lvl %d (ab %f %f): %f %f %f %f %f -> %f %f %f %f %f \n",
                   i,
                   i, calibration->pyrLevelsUsed - 1,
                   aff_g2l_this.a, aff_g2l_this.b,
                   achievedRes[0],
                   achievedRes[1],
                   achievedRes[2],
                   achievedRes[3],
                   achievedRes[4],
                   coarseTracker->lastResiduals[0],
                   coarseTracker->lastResiduals[1],
                   coarseTracker->lastResiduals[2],
                   coarseTracker->lastResiduals[3],
                   coarseTracker->lastResiduals[4]);
        }


        // do we have a new winner?
        Log::Logger::getInstance()->info("Tracking. is good? {} Last Residual: {} achievedRes {}", trackingIsGood,
                                         coarseTracker->lastResiduals[0], achievedRes[0]);
        if (trackingIsGood && std::isfinite((float) coarseTracker->lastResiduals[0]) &&
            !(coarseTracker->lastResiduals[0] >= achievedRes[0])) {
            //printf("take over. minRes %f -> %f!\n", achievedRes[0], coarseTracker->lastResiduals[0]);
            flowVecs = coarseTracker->lastFlowIndicators;
            aff_g2l = aff_g2l_this;
            lastF_2_fh = lastF_2_fh_this;
            haveOneGood = true;
        }

        // take over achieved res (always).
        if (haveOneGood) {
            for (int i = 0; i < 5; i++) {
                if (!std::isfinite((float) achievedRes[i]) || achievedRes[i] >
                                                              coarseTracker->lastResiduals[i])    // take over if achievedRes is either bigger or NAN.
                    achievedRes[i] = coarseTracker->lastResiduals[i];
            }
        }


        if (haveOneGood && achievedRes[0] < lastCoarseRMSE[0] * setting_reTrackThreshold)
            break;

    }

    if (!haveOneGood) {
        printf("BIG ERROR! tracking failed entirely. Take predictred pose and hope we may somehow recover.\n");
        flowVecs = Vec3(0, 0, 0);
        aff_g2l = aff_last_2_l;
        lastF_2_fh = lastF_2_fh_tries[0];
    }

    lastCoarseRMSE = achievedRes;

    fh->shell->camToTrackingRef = lastF_2_fh.inverse();
    fh->shell->trackingRef = lastF->shell;
    fh->shell->aff_g2l = aff_g2l;
    fh->shell->camToWorld = fh->shell->trackingRef->camToWorld * fh->shell->camToTrackingRef;

    if (coarseTracker->firstCoarseRMSE < 0)
        coarseTracker->firstCoarseRMSE = achievedRes[0];


    Log::Logger::getInstance()->info("Coarse Tracker tracked ab =  {}  {} (exp  {}). Res  {}!", aff_g2l.a, aff_g2l.b, fh->abExposure, achievedRes[0]);


    return Vec4(achievedRes[0], flowVecs[0], flowVecs[1], flowVecs[2]);
}

void Tracker::traceNewCoarse(std::shared_ptr<VO::Frame> frame) {
    int trace_total = 0, trace_good = 0, trace_oob = 0, trace_out = 0, trace_skip = 0, trace_badcondition = 0, trace_uninitialized = 0;
    Mat33f K = Mat33f::Identity();
    K(0, 0) = hCalib.fxl();
    K(1, 1) = hCalib.fyl();
    K(0, 2) = hCalib.cxl();
    K(1, 2) = hCalib.cyl();

    for (auto &host: frameHessians)        // go through all active frames
    {

        SE3 hostToNew = frame->PRE_worldToCam * host->PRE_camToWorld;
        Mat33f KRKi = K * hostToNew.rotationMatrix().cast<float>() * K.inverse();
        Vec3f Kt = K * hostToNew.translation().cast<float>();


        Vec2f aff = AffLight::fromToVecExposure(host->abExposure, frame->abExposure, host->aff_g2l(), frame->aff_g2l()).cast<float>();

       Log::Logger::getInstance()->info("Tracing points. Frame id: {}, tracking id: {}, Immature points: {} ", host->frameID, host->trackingID,  host->immaturePoints.size());

        for (ImmaturePoint *ph: host->immaturePoints) {
            ph->traceOn(frame, KRKi, Kt, aff, &hCalib, calibration, true);
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_GOOD) trace_good++;
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_BADCONDITION) trace_badcondition++;
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_OOB) trace_oob++;
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_OUTLIER) trace_out++;
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_SKIPPED) trace_skip++;
            if (ph->lastTraceStatus == ImmaturePointStatus::IPS_UNINITIALIZED) trace_uninitialized++;
            trace_total++;
        }
    }

}

void Tracker::flagFramesForMarginalization(std::shared_ptr<VO::Frame> sharedPtr) {
    if (setting_minFrameAge > setting_maxFrames) {
        for (int i = setting_maxFrames; i < (int) frameHessians.size(); i++) {
            std::shared_ptr<VO::Frame> fh = frameHessians[i - setting_maxFrames];
            fh->flaggedForMarginalization = true;
        }
        return;
    }


    int flagged = 0;

    // marginalize all frames that have not enough points.
    for (int i = 0; i < (int) frameHessians.size(); i++) {
        std::shared_ptr<VO::Frame> fh = frameHessians[i];

        int in = fh->pointHessians.size() + fh->immaturePoints.size();
        size_t out = fh->pointHessiansMarginalized.size() + fh->pointHessiansOut.size();


        Vec2 refToFh=AffLight::fromToVecExposure(frameHessians.back()->abExposure, fh->abExposure,
                                                 frameHessians.back()->aff_g2l(), fh->aff_g2l());

        if ((in < setting_minPointsRemaining * (in + out) ||
             fabs(logf((float) refToFh[0])) > setting_maxLogAffFacInWindow)
            && ((int) frameHessians.size()) - flagged > setting_minFrames) {
            /*
             Log::Logger::getInstance()->info("MARGINALIZE frame {}, as only {}{} points remaining (%'d %'d %'d %'d). VisInLast %'d / %'d. traces %d, activated %d!\n",
                    fh->frameID, in, in+out,
                    (int)fh->pointHessians.size(), (int)fh->immaturePoints.size(),
                    (int)fh->pointHessiansMarginalized.size(), (int)fh->pointHessiansOut.size(),
                    visInLast, outInLast,
                    fh->statistics_tracesCreatedForThisFrame, fh->statistics_pointsActivatedForThisFrame);
             */
            fh->flaggedForMarginalization = true;
            flagged++;
        } else {
//			printf("May Keep frame %d, as %'d/%'d points remaining (%'d %'d %'d %'d). VisInLast %'d / %'d. traces %d, activated %d!\n",
//					fh->frameID, in, in+out,
//					(int)fh->pointHessians.size(), (int)fh->immaturePoints.size(),
//					(int)fh->pointHessiansMarginalized.size(), (int)fh->pointHessiansOut.size(),
//					visInLast, outInLast,
//					fh->statistics_tracesCreatedForThisFrame, fh->statistics_pointsActivatedForThisFrame);
        }
    }

    // marginalize one.
    if ((int) frameHessians.size() - flagged >= setting_maxFrames) {
        double smallestScore = 1;
        std::shared_ptr<VO::Frame> toMarginalize = nullptr;
        std::shared_ptr<VO::Frame> latest = frameHessians.back();


        for (std::shared_ptr<VO::Frame> &fh: frameHessians) {
            if(fh->frameID > latest->frameID-setting_minFrameAge || fh->frameID == 0) continue;
            //if(fh==frameHessians.front() == 0) continue;

            double distScore = 0;
            for (VO::FrameToFramePrecalc ffh: fh->targetPrecalc) {
                if(ffh.target->frameID > latest->frameID-setting_minFrameAge+1 || ffh.target == ffh.host) continue;
                distScore += 1/(1e-5+ffh.distanceLL);

            }
            distScore *= -sqrtf(fh->targetPrecalc.back().distanceLL);


            if (distScore < smallestScore) {
                smallestScore = distScore;
                toMarginalize = fh;
            }
        }

        Log::Logger::getInstance()->info("MARGINALIZE frame {}, as it is the closest (score {:.2f})!",
                                         toMarginalize->frameID, smallestScore);
        toMarginalize->flaggedForMarginalization = true;
        flagged++;
    }

//	printf("FRAMES LEFT: ");
//	for(FrameHessian* fh : frameHessians)
//		printf("%d ", fh->frameID);
//	printf("\n");
}


void Tracker::activatePointsMT() {
    if (ef.nPoints < setting_desiredPointDensity * 0.66)
        currentMinActDist -= 0.8;
    if (ef.nPoints < setting_desiredPointDensity * 0.8)
        currentMinActDist -= 0.5;
    else if (ef.nPoints < setting_desiredPointDensity * 0.9)
        currentMinActDist -= 0.2;
    else if (ef.nPoints < setting_desiredPointDensity)
        currentMinActDist -= 0.1;

    if (ef.nPoints > setting_desiredPointDensity * 1.5)
        currentMinActDist += 0.8;
    if (ef.nPoints > setting_desiredPointDensity * 1.3)
        currentMinActDist += 0.5;
    if (ef.nPoints > setting_desiredPointDensity * 1.15)
        currentMinActDist += 0.2;
    if (ef.nPoints > setting_desiredPointDensity)
        currentMinActDist += 0.1;

    if (currentMinActDist < 0) currentMinActDist = 0;
    if (currentMinActDist > 4) currentMinActDist = 4;


    Log::Logger::getInstance()->info("SPARSITY:  MinActDist {} (need {} points, have {} points)!",
                                     currentMinActDist, (int) (setting_desiredPointDensity), ef.nPoints);


    std::shared_ptr<VO::Frame> newestHs = frameHessians.back();

    // make dist map.
    coarseDistanceMap.makeK(&hCalib);
    coarseDistanceMap.makeDistanceMap(frameHessians, newestHs);

    //coarseTracker->debugPlotDistMap("distMap");

    std::vector<ImmaturePoint *> toOptimize;
    toOptimize.reserve(20000);


    for (auto &host: frameHessians)        // go through all active frames
    {
        if (host == newestHs) continue;

        SE3 fhToNew = newestHs->PRE_worldToCam * host->PRE_camToWorld;
        Mat33f KRKi = (coarseDistanceMap.K[1] * fhToNew.rotationMatrix().cast<float>() * coarseDistanceMap.Ki[0]);
        Vec3f Kt = (coarseDistanceMap.K[1] * fhToNew.translation().cast<float>());


        for (unsigned int i = 0; i < host->immaturePoints.size(); i += 1) {
            ImmaturePoint *ph = host->immaturePoints[i];
            ph->idxInImmaturePoints = i;

            // delete points that have never been traced successfully, or that are outlier on the last trace.
            if (!std::isfinite(ph->idepth_max) || ph->lastTraceStatus == IPS_OUTLIER) {
//				immature_invalid_deleted++;
                // remove point.
                delete ph;
                host->immaturePoints[i] = 0;
                continue;
            }

            // can activate only if this is true.
            bool canActivate = (ph->lastTraceStatus == IPS_GOOD
                                || ph->lastTraceStatus == IPS_SKIPPED
                                || ph->lastTraceStatus == IPS_BADCONDITION
                                || ph->lastTraceStatus == IPS_OOB)
                               && ph->lastTracePixelInterval < 8
                               && ph->quality > setting_minTraceQuality
                               && (ph->idepth_max + ph->idepth_min) > 0;


            // if I cannot activate the point, skip it. Maybe also delete it.
            if (!canActivate) {
                // if point will be out afterwards, delete it instead.
                if (ph->host->flaggedForMarginalization || ph->lastTraceStatus == IPS_OOB) {
//					immature_notReady_deleted++;
                    delete ph;
                    host->immaturePoints[i] = 0;
                }
//				immature_notReady_skipped++;
                continue;
            }


            // see if we need to activate point due to distance map.
            Vec3f ptp = KRKi * Vec3f(ph->u, ph->v, 1) + Kt * (0.5f * (ph->idepth_max + ph->idepth_min));
            int u = ptp[0] / ptp[2] + 0.5f;
            int v = ptp[1] / ptp[2] + 0.5f;

            if ((u > 0 && v > 0 && u < calibration->wG[1] && v < calibration->hG[1])) {

                float dist = coarseDistanceMap.fwdWarpedIDDistFinal[u + calibration->wG[1] * v] +
                             (ptp[0] - floorf((float) (ptp[0])));

                if (dist >= currentMinActDist * ph->my_type) {
                    coarseDistanceMap.addIntoDistFinal(u, v);
                    toOptimize.push_back(ph);
                }
            } else {
                delete ph;
                host->immaturePoints[i] = 0;
            }
        }
    }


//	printf("ACTIVATE: %d. (del %d, notReady %d, marg %d, good %d, marg-skip %d)\n",
//			(int)toOptimize.size(), immature_deleted, immature_notReady, immature_needMarg, immature_want, immature_margskip);

    std::vector<PointHessian *> optimized;
    optimized.resize(toOptimize.size());


    activatePointsMT_Reductor(&optimized, &toOptimize, 0, toOptimize.size(), 0, 0);

    Log::Logger::getInstance()->info("Points to optimize: {}", toOptimize.size());

    for (unsigned k = 0; k < toOptimize.size(); k++) {
        PointHessian *newpoint = optimized[k];
        ImmaturePoint *ph = toOptimize[k];

        if (newpoint != 0 && newpoint != (PointHessian *) ((long) (-1))) {
            newpoint->host->immaturePoints[ph->idxInImmaturePoints] = 0;
            newpoint->host->pointHessians.push_back(newpoint);
            ef.insertPoint(newpoint);
            for (PointFrameResidual *r: newpoint->residuals)
                ef.insertResidual(r);
            assert(newpoint->efPoint != 0);
            delete ph;
        } else if (newpoint == (PointHessian *) ((long) (-1)) || ph->lastTraceStatus == IPS_OOB) {
            delete ph;
            ph->host->immaturePoints[ph->idxInImmaturePoints] = 0;
        } else {
            assert(newpoint == 0 || newpoint == (PointHessian *) ((long) (-1)));
        }
    }

    for (auto &host: frameHessians) {
        //std::cout << "activated pointHessians size: " <<  host->pointHessians.size() << std::endl;
        for (int i = 0; i < (int) host->immaturePoints.size(); i++) {
            if (host->immaturePoints[i] == 0) {
                host->immaturePoints[i] = host->immaturePoints.back();
                host->immaturePoints.pop_back();
                i--;
            }
        }
    }
}

void Tracker::setPrecalcValues() {
    for (auto &fh: frameHessians) {
        fh->targetPrecalc.resize(frameHessians.size());

        for (unsigned int i = 0; i < frameHessians.size(); i++)
            fh->targetPrecalc[i].set(fh.get(), frameHessians[i].get(), &hCalib);
    }

    ef.setDeltaF(&hCalib);
}


void Tracker::activatePointsMT_Reductor(
        std::vector<PointHessian *> *optimized,
        std::vector<ImmaturePoint *> *toOptimize,
        int min, int max, Vec10 *stats, int tid) {
    ImmaturePointTemporaryResidual *tr = new ImmaturePointTemporaryResidual[frameHessians.size()];
    for (int k = min; k < max; k++) {
        (*optimized)[k] = optimizeImmaturePoint((*toOptimize)[k], 1, tr);
    }
    delete[] tr;
}

PointHessian *Tracker::optimizeImmaturePoint(
        ImmaturePoint *point, int minObs,
        ImmaturePointTemporaryResidual *residuals) {
    int nres = 0;
    for (auto &fh: frameHessians) {
        if (fh.get() != point->host) {
            residuals[nres].state_NewEnergy = residuals[nres].state_energy = 0;
            residuals[nres].state_NewState = ResState::OUTLIER;
            residuals[nres].state_state = ResState::IN;
            residuals[nres].target = fh.get();
            nres++;
        }
    }
    assert(nres == ((int) frameHessians.size()) - 1);

    bool print = false;//rand()%50==0;

    float lastEnergy = 0;
    float lastHdd = 0;
    float lastbd = 0;
    float currentIdepth = (point->idepth_max + point->idepth_min) * 0.5f;


    for (int i = 0; i < nres; i++) {
        lastEnergy += point->linearizeResidual(&hCalib, 1000, residuals + i, lastHdd, lastbd, currentIdepth,
                                               (residuals + i)->target);
        residuals[i].state_state = residuals[i].state_NewState;
        residuals[i].state_energy = residuals[i].state_NewEnergy;
    }

    if (!std::isfinite(lastEnergy) || lastHdd < setting_minIdepthH_act) {
        if (print)
            printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
                   nres, lastHdd, lastEnergy);
        return 0;
    }

    if (print)
        printf("Activate point. %d residuals. H=%f. Initial Energy: %f. Initial Id=%f\n",
               nres, lastHdd, lastEnergy, currentIdepth);

    float lambda = 0.1;
    for (int iteration = 0; iteration < setting_GNItsOnPointActivation; iteration++) {
        float H = lastHdd;
        H *= 1 + lambda;
        float step = (1.0 / H) * lastbd;
        float newIdepth = currentIdepth - step;

        float newHdd = 0;
        float newbd = 0;
        float newEnergy = 0;
        for (int i = 0; i < nres; i++)
            newEnergy += point->linearizeResidual(&hCalib, 1, residuals + i, newHdd, newbd, newIdepth,
                                                  (residuals + i)->target);

        if (!std::isfinite(lastEnergy) || newHdd < setting_minIdepthH_act) {
            if (print)
                printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
                       nres,
                       newHdd,
                       lastEnergy);
            return 0;
        }

        if (print)
            printf("%s %d (L %.2f) %s: %f -> %f (idepth %f)!\n",
                   (true || newEnergy < lastEnergy) ? "ACCEPT" : "REJECT",
                   iteration,
                   log10(lambda),
                   "",
                   lastEnergy, newEnergy, newIdepth);

        if (newEnergy < lastEnergy) {
            currentIdepth = newIdepth;
            lastHdd = newHdd;
            lastbd = newbd;
            lastEnergy = newEnergy;
            for (int i = 0; i < nres; i++) {
                residuals[i].state_state = residuals[i].state_NewState;
                residuals[i].state_energy = residuals[i].state_NewEnergy;
            }

            lambda *= 0.5;
        } else {
            lambda *= 5;
        }

        if (fabsf(step) < 0.0001 * currentIdepth)
            break;
    }

    if (!std::isfinite(currentIdepth)) {
        printf("MAJOR ERROR! point idepth is nan after initialization (%f).\n", currentIdepth);
        return (PointHessian *) ((long) (-1));        // yeah I'm like 99% sure this is OK on 32bit systems.
    }


    int numGoodRes = 0;
    for (int i = 0; i < nres; i++)
        if (residuals[i].state_state == ResState::IN) numGoodRes++;

    if (numGoodRes < minObs) {
        if (print) printf("OptPoint: OUTLIER!\n");
        return (PointHessian *) ((long) (-1));        // yeah I'm like 99% sure this is OK on 32bit systems.
    }


    PointHessian *p = new PointHessian(point, &hCalib);
    if (!std::isfinite(p->energyTH)) {
        delete p;
        return (PointHessian *) ((long) (-1));
    }


    p->lastResiduals[0].first = 0;
    p->lastResiduals[0].second = ResState::OOB;
    p->lastResiduals[1].first = 0;
    p->lastResiduals[1].second = ResState::OOB;
    p->setIdepthZero(currentIdepth);
    p->setIdepth(currentIdepth);
    p->setPointStatus(PointHessian::ACTIVE);

    for (int i = 0; i < nres; i++)
        if (residuals[i].state_state == ResState::IN) {
            PointFrameResidual *r = new PointFrameResidual(p, p->host, residuals[i].target);
            r->state_NewEnergy = r->state_energy = 0;
            r->state_NewState = ResState::OUTLIER;
            r->setState(ResState::IN);
            p->residuals.push_back(r);

            if (r->target == frameHessians.back().get()) {
                p->lastResiduals[0].first = r;
                p->lastResiduals[0].second = ResState::IN;
            } else if (r->target == (frameHessians.size() < 2 ? 0 : frameHessians[frameHessians.size() - 2]).get()) {
                p->lastResiduals[1].first = r;
                p->lastResiduals[1].second = ResState::IN;
            }
        }
    return p;
}

void Tracker::printEigenValLine() {
    if (ef.lastHS.rows() < 12) return;


    MatXX Hp = ef.lastHS.bottomRightCorner(ef.lastHS.cols() - CPARS, ef.lastHS.cols() - CPARS);
    MatXX Ha = ef.lastHS.bottomRightCorner(ef.lastHS.cols() - CPARS, ef.lastHS.cols() - CPARS);
    int n = Hp.cols() / 8;
    assert(Hp.cols() % 8 == 0);

    // sub-select
    for (int i = 0; i < n; i++) {
        MatXX tmp6 = Hp.block(i * 8, 0, 6, n * 8);
        Hp.block(i * 6, 0, 6, n * 8) = tmp6;

        MatXX tmp2 = Ha.block(i * 8 + 6, 0, 2, n * 8);
        Ha.block(i * 2, 0, 2, n * 8) = tmp2;
    }
    for (int i = 0; i < n; i++) {
        MatXX tmp6 = Hp.block(0, i * 8, n * 8, 6);
        Hp.block(0, i * 6, n * 8, 6) = tmp6;

        MatXX tmp2 = Ha.block(0, i * 8 + 6, n * 8, 2);
        Ha.block(0, i * 2, n * 8, 2) = tmp2;
    }

    VecX eigenvaluesAll = ef.lastHS.eigenvalues().real();
    VecX eigenP = Hp.topLeftCorner(n * 6, n * 6).eigenvalues().real();
    VecX eigenA = Ha.topLeftCorner(n * 2, n * 2).eigenvalues().real();
    VecX diagonal = ef.lastHS.diagonal();

    std::sort(eigenvaluesAll.data(), eigenvaluesAll.data() + eigenvaluesAll.size());
    std::sort(eigenP.data(), eigenP.data() + eigenP.size());
    std::sort(eigenA.data(), eigenA.data() + eigenA.size());

    int nz = std::max(100, setting_maxFrames * 10);

    {
        VecX ea = VecX::Zero(nz);
        ea.head(eigenvaluesAll.size()) = eigenvaluesAll;
        (*eigenAllLog) << frameHessians.back()->frameID << " " << ea.transpose() << "\n";
        eigenAllLog->flush();
    }
    {
        VecX ea = VecX::Zero(nz);
        ea.head(eigenA.size()) = eigenA;
        (*eigenALog) << frameHessians.back()->frameID << " " << ea.transpose() << "\n";
        eigenALog->flush();
    }
    {
        VecX ea = VecX::Zero(nz);
        ea.head(eigenP.size()) = eigenP;
        (*eigenPLog) << frameHessians.back()->frameID << " " << ea.transpose() << "\n";
        eigenPLog->flush();
    }

    {
        VecX ea = VecX::Zero(nz);
        ea.head(diagonal.size()) = diagonal;
        (*DiagonalLog) << frameHessians.back()->frameID << " " << ea.transpose() << "\n";
        DiagonalLog->flush();
    }

    {
        VecX ea = VecX::Zero(nz);
        ea.head(diagonal.size()) = ef.lastHS.inverse().diagonal();
        (*variancesLog) << frameHessians.back()->frameID << " " << ea.transpose() << "\n";
        variancesLog->flush();
    }

    std::vector<VecX> &nsp = ef.lastNullspaces_forLogging;
}

void Tracker::applyRes_Reductor(bool b, int min, int max, Vec10 *stats, int tid) {
    for (int k = min; k < max; k++)
        activeResiduals[k]->applyRes(true);
}

void Tracker::flagPointsForRemoval() {
    assert(EFIndicesValid);

    std::vector<VO::Frame *> fhsToKeepPoints;
    std::vector<VO::Frame *> fhsToMargPoints;

    //if(setting_margPointVisWindow>0)
    {
        for (int i = ((int) frameHessians.size()) - 1; i >= 0 && i >= ((int) frameHessians.size()); i--)
            if (!frameHessians[i]->flaggedForMarginalization) fhsToKeepPoints.push_back(frameHessians[i].get());

        for (int i = 0; i < (int) frameHessians.size(); i++)
            if (frameHessians[i]->flaggedForMarginalization) fhsToMargPoints.push_back(frameHessians[i].get());
    }



    //ef->setAdjointsF();
    //ef->setDeltaF(&Hcalib);
    int flag_oob = 0, flag_in = 0, flag_inin = 0, flag_nores = 0;

    for (auto &host: frameHessians)        // go through all active frames
    {
        for (unsigned int i = 0; i < host->pointHessians.size(); i++) {
            PointHessian *ph = host->pointHessians[i];
            if (ph == 0) continue;

            if (ph->idepth_scaled < 0 || ph->residuals.size() == 0) {
                host->pointHessiansOut.push_back(ph);
                ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                host->pointHessians[i] = 0;
                flag_nores++;
            } else if (ph->isOOB(fhsToKeepPoints, fhsToMargPoints) || host->flaggedForMarginalization) {
                flag_oob++;
                if (ph->isInlierNew()) {
                    flag_in++;
                    int ngoodRes = 0;
                    for (PointFrameResidual *r: ph->residuals) {
                        r->resetOOB();
                        r->linearize(&hCalib);
                        r->efResidual->isLinearized = false;
                        r->applyRes(true);
                        if (r->efResidual->isActive()) {
                            r->efResidual->fixLinearizationF(&ef);
                            ngoodRes++;
                        }
                    }
                    if (ph->idepth_hessian > setting_minIdepthH_marg) {
                        flag_inin++;
                        ph->efPoint->stateFlag = EFPointStatus::PS_MARGINALIZE;
                        host->pointHessiansMarginalized.push_back(ph);
                    } else {
                        ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                        host->pointHessiansOut.push_back(ph);
                    }


                } else {
                    host->pointHessiansOut.push_back(ph);
                    ph->efPoint->stateFlag = EFPointStatus::PS_DROP;


                    //printf("drop point in frame %d (%d goodRes, %d activeRes)\n", ph->host->idx, ph->numGoodResiduals, (int)ph->residuals.size());
                }

                host->pointHessians[i] = 0;
            }
        }


        for (int i = 0; i < (int) host->pointHessians.size(); i++) {
            if (host->pointHessians[i] == 0) {
                host->pointHessians[i] = host->pointHessians.back();
                host->pointHessians.pop_back();
                i--;
            }
        }
    }

}
