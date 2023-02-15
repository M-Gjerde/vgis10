//
// Created by magnus on 2/7/23.
//

#include "Tracker.h"

#include "TrackerInitializer.h"
#include "PointHessian.h"
#include "Optimized/PointFrameResidual.h"

void Tracker::initializeFromInitializer() {
    frameToFramePreCalc.resize(2);
    // Add First Frame

    firstFrame->trackingID = 0;
    ef.insertFrame(firstFrame, &hCalib);

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

        ImmaturePoint pt(point->u + 0.5f, point->v + 0.5f, point->my_type, firstFrame); // why add 0.5f?
        if (!std::isfinite(pt.energyTH)) continue;
        pt.idepth_max = pt.idepth_min = 1;

        PointHessian ph;
        VO::pointHessianFromImmaturePoint(&pt, &ph);
        if (!std::isfinite(ph.energyTH)) continue;
        ph.setIdepthScaled(point->iR * rescaleFactor);
        ph.setIdepthZero(ph.idepth);
        ph.hasDepthPrior = true;
        ph.setPointStatus(PointHessian::ACTIVE);
        ph.pointIndex = ptIdx;
        ptIdx++;
        firstFrame->pointHessians.push_back(ph);
        ef.insertPoint(&firstFrame->pointHessians.back(), firstFrame->trackingID);
    }

    SE3 firstToNew = info.thisToNext;
    firstToNew.translation() /= rescaleFactor;

    frameHessians.emplace_back(firstFrame);
    setPrecalcValues();

    VO::FramePose pose;
    pose.frameToWorld = SE3();
    pose.affLightg2l = AffLight(0, 0);
    pose.trackingRefFrameID = -1;
    pose.frameToTrackingRef = SE3();
    pose.setEvalPT_scaled(pose.frameToWorld, pose.affLightg2l);
    frameToFramePreCalc[pose.trackingID].set(pose, pose, &hCalib);
    firstFrame->pose = pose;

    VO::FramePose pose2;
    pose2.frameToWorld = firstToNew.inverse();
    pose2.affLightg2l = AffLight(0, 0);
    pose2.trackingRefFrameID = 0;
    pose2.frameToTrackingRef = firstToNew.inverse();
    pose2.setEvalPT_scaled(pose2.frameToWorld, pose2.affLightg2l);
    frameToFramePreCalc[pose2.trackingID].set(pose, pose2, &hCalib);
    secondFrame->pose = pose2;



    initialized = true;

    Log::Logger::getInstance()->info("INITIALIZE FROM INITIALIZER ({} pts)!", (int) firstFrame->pointHessians.size());
}


bool Tracker::initializerTrackFrame(const std::shared_ptr<VO::Frame> &frame, const CameraCalibration *calibration) {
    int maxIterations[] = {5, 5, 10, 30, 50};
    secondFrame = frame;
    // frame->displayPyramid(0, false, "frame");
    //firstFrame->displayPyramid(0,true, "firstFrame");

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

void Tracker::takeTrackedFrame(std::shared_ptr<VO::Frame> frame, bool needKF) {
    if (needKF) {
        makeKeyFrame(frame);
    } else {
        makeNonKeyFrame(frame);
    }
}

void Tracker::makeKeyFrame(std::shared_ptr<VO::Frame> frame) {
    // =========================== UPDATE POSE =========================

    frame->pose.frameToWorld = frame->trackingRef.frameToWorld * frame->pose.frameToTrackingRef;
    frame->pose.setEvalPT_scaled(frame->pose.frameToWorld, frame->pose.affLightg2l);
    // traceNewCoarse
    Log::Logger::getInstance()->info("Tracing new coarse");
    traceNewCoarse(frame);
    // =========================== Flag Frames to be Marginalized. =========================
    Log::Logger::getInstance()->info("Flagging frames for marginalization");
    flagFramesForMarginalization(frame);
    // =========================== add New Frame to Hessian Struct. =========================


    frame->trackingID = frameHessians.size();
    frameHessians.emplace_back(frame);
    setPrecalcValues();
    ef.insertFrame(frame, &hCalib);

    // =========================== add new residuals for old points =========================
    for (auto &fh1: frameHessians) {
        if (fh1 == frame) continue;
        for (auto &ph: fh1->pointHessians) {
            PointFrameResidual r(ph.pointIndex);
            r.host = fh1.get();
            r.target = frame.get();
            r.pointH = &ph;
            r.setState(ResState::IN);
            ph.residuals.push_back(r);
            ef.insertResidual(&ph.residuals.back(), fh1->trackingID, frame->trackingID);
            ph.lastResiduals[1] = ph.lastResiduals[0];
            ph.lastResiduals[0] = std::pair<PointFrameResidual *, ResState>(&ph.residuals.back(), ResState::IN);
        }
    }

    // =========================== Activate Points (& flag for marginalization). =========================
    activatePointsMT();
    // =========================== OPTIMIZE ALL =========================
    float rmse = optimize(setting_maxOptIterations);
    Log::Logger::getInstance()->info("RMSE: {}", rmse);

    // =========================== Figure Out if INITIALIZATION FAILED =========================
    // =========================== REMOVE OUTLIER =========================
    // =========================== (Activate-)Marginalize Points =========================
    // =========================== add new Immature points & new residuals =========================
    // =========================== Marginalize Frames =========================
/*
    for (auto &host: frameHessians) {
        cv::Mat img(calibration->hG[0], calibration->wG[0], CV_32FC3, host->pyramid[0].data());
        img.convertTo(img, CV_8UC3);

        std::cout << "pointHessianss: " << host->pointHessians.size() << std::endl;
        for (PointHessian &ph: host->pointHessians) {
            cv::circle(img, cv::Point2f(ph.u, ph.v), 1,
                       cv::Scalar(0, 0, 255));

        }

        cv::imshow("PointHessian Activated: " + std::to_string(host->trackingID), img);
    }
    cv::waitKey(0);
*/
}


void Tracker::makeNonKeyFrame(std::shared_ptr<VO::Frame> frame) {
    // =========================== UPDATE POSE =========================

}

void Tracker::traceNewCoarse(std::shared_ptr<VO::Frame> frame) {
    int trace_total = 0, trace_good = 0, trace_oob = 0, trace_out = 0, trace_skip = 0, trace_badcondition = 0, trace_uninitialized = 0;
    Mat33f K = Mat33f::Identity();
    K(0, 0) = hCalib.fxl();
    K(1, 1) = hCalib.fyl();
    K(0, 2) = hCalib.cxl();
    K(1, 2) = hCalib.cyl();

    for (size_t i = 0; i < frameHessians.size(); ++i) {
        SE3 newPose = frameHessians[i]->pose.PRE_worldToCam * frameHessians[i]->pose.PRE_camToWorld;
        Mat33f KRKi = K * newPose.rotationMatrix().cast<float>() * K.inverse();
        Vec3f Kt = K * newPose.translation().cast<float>();
        Vec2f aff = AffLight::fromToVecExposure(frameHessians[i]->abExposure, frame->abExposure,
                                                frameHessians[i]->pose.aff_g2l(),
                                                frameHessians[frame->trackingID]->pose.aff_g2l()).cast<float>();

        Log::Logger::getInstance()->info("Tracking frame {} --> {}", frameHessians[i]->trackingID, frame->trackingID);
        if (i >= immaturePoints.size())
            continue;
        for (auto &ip: immaturePoints[i]) {
            ip.traceOn(frame, KRKi, Kt, aff, &hCalib, calibration);
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_GOOD) trace_good++;
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_BADCONDITION) trace_badcondition++;
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_OOB) trace_oob++;
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_OUTLIER) trace_out++;
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_SKIPPED) trace_skip++;
            if (ip.lastTraceStatus == ImmaturePointStatus::IPS_UNINITIALIZED) trace_uninitialized++;
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
        Log::Logger::getInstance()->info("{}", fh->trackingID);

        size_t in = fh->pointHessians.size() + immaturePoints[fh->trackingID].size();
        size_t out = fh->pointHessiansMarginalized.size() + fh->pointHessiansOut.size();

        Vec2 refToFh = AffLight::fromToVecExposure(frameHessians.back()->abExposure, fh->abExposure,
                                                   frameHessians.back()->pose.aff_g2l(),
                                                   frameHessians[i]->pose.aff_g2l());


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
            // TODO check if I'm using correct id
            if (fh->id > latest->id - setting_minFrameAge || fh->id == 0) continue;
            //if(fh==frameHessians.front() == 0) continue;

            double distScore = 0;
            for (std::shared_ptr<VO::Frame> &ffh: frameHessians) {
                Log::Logger::getInstance()->info("{}", ffh->trackingID);
                // TODO check if I'm using correct id and check frameToFramePrecalc is correct
                if (ffh->id > latest->id - setting_minFrameAge + 1 || ffh == fh) continue;
                distScore += 1 / (1e-5 + frameToFramePreCalc[ffh->trackingID].distanceLL);

            }
            distScore *= -sqrtf(frameToFramePreCalc.back().distanceLL);


            if (distScore < smallestScore) {
                smallestScore = distScore;
                toMarginalize = fh;
            }
        }

        Log::Logger::getInstance()->info("MARGINALIZE frame {}, as it is the closest (score {:.2f})!",
                                         toMarginalize->id, smallestScore);
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


    std::shared_ptr<VO::Frame> latestFrame = frameHessians.back();

    // make dist map.
    coarseDistanceMap.makeK(&hCalib);
    coarseDistanceMap.makeDistanceMap(frameHessians, latestFrame);

    //coarseTracker->debugPlotDistMap("distMap");

    std::vector<ImmaturePoint> toOptimize;
    toOptimize.reserve(20000);


    for (size_t i = 0; i < frameHessians.size(); ++i)        // go through all active frames
    {
        if (frameHessians[i] == latestFrame) continue;

        SE3 fhToNew = latestFrame->pose.PRE_worldToCam * frameHessians[i]->pose.PRE_camToWorld;
        Mat33f KRKi = (coarseDistanceMap.K[1] * fhToNew.rotationMatrix().cast<float>() * coarseDistanceMap.Ki[0]);
        Vec3f Kt = (coarseDistanceMap.K[1] * fhToNew.translation().cast<float>());


        for (unsigned int j = 0; j < immaturePoints[i].size(); j += 1) {
            ImmaturePoint &ph = immaturePoints[i][j];
            ph.idxInImmaturePoints = j;

            // delete points that have never been traced successfully, or that are outlier on the last trace.
            if (!std::isfinite(ph.idepth_max) || ph.lastTraceStatus == IPS_OUTLIER) {
//				immature_invalid_deleted++;
                // remove point.
                auto it = immaturePoints[i].begin();
                std::advance(it, j);
                immaturePoints[i].erase(it);
                continue;
            }

            // can activate only if this is true.
            bool canActivate = (ph.lastTraceStatus == IPS_GOOD
                                || ph.lastTraceStatus == IPS_SKIPPED
                                || ph.lastTraceStatus == IPS_BADCONDITION
                                || ph.lastTraceStatus == IPS_OOB)
                               && ph.lastTracePixelInterval < 8
                               && ph.quality > setting_minTraceQuality
                               && (ph.idepth_max + ph.idepth_min) > 0;


            // if I cannot activate the point, skip it. Maybe also delete it.
            if (!canActivate) {
                // if point will be out afterwards, delete it instead.
                if (ph.host->flaggedForMarginalization || ph.lastTraceStatus == IPS_OOB) {
//					immature_notReady_deleted++;
                    auto it = immaturePoints[i].begin();
                    std::advance(it, j);
                    immaturePoints[i].erase(it);
                }
//				immature_notReady_skipped++;
                continue;
            }


            // see if we need to activate point due to distance map.
            Vec3f ptp = KRKi * Vec3f(ph.u, ph.v, 1) + Kt * (0.5f * (ph.idepth_max + ph.idepth_min));
            int u = ptp[0] / ptp[2] + 0.5f;
            int v = ptp[1] / ptp[2] + 0.5f;

            if ((u > 0 && v > 0 && u < calibration->wG[1] && v < calibration->hG[1])) {

                float dist = coarseDistanceMap.fwdWarpedIDDistFinal[u + calibration->wG[1] * v] +
                             (ptp[0] - floorf((float) (ptp[0])));

                if (dist >= currentMinActDist * ph.type) {
                    coarseDistanceMap.addIntoDistFinal(u, v);
                    toOptimize.push_back(ph);
                }
            } else {
                auto it = immaturePoints[i].begin();
                std::advance(it, j);
                immaturePoints[i].erase(it);
            }
        }
    }

    Log::Logger::getInstance()->info("Points to optimize: {}", toOptimize.size());

//	printf("ACTIVATE: %d. (del %d, notReady %d, marg %d, good %d, marg-skip %d)\n",
//			(int)toOptimize.size(), immature_deleted, immature_notReady, immature_needMarg, immature_want, immature_margskip);

    std::vector<PointHessian> optimized;
    optimized.resize(toOptimize.size());

    activatePointsMT_Reductor(&optimized, &toOptimize, 0, toOptimize.size(), 0, 0);

    Log::Logger::getInstance()->info("Points to optimize: {}", toOptimize.size());

    for (unsigned k = 0; k < toOptimize.size(); k++) {
        PointHessian &newpoint = optimized[k];
        ImmaturePoint &ph = toOptimize[k];

        if (newpoint.isGood()) {

            //newpoint->host->immaturePoints[ph.idxInImmaturePoints]=0;
            //newpoint->host->pointHessians.push_back(newpoint);
            //ef.insertPoint(newpoint);
            for (PointFrameResidual &r: newpoint.residuals)
                ef.insertResidual(&r, 0, 0);
            assert(newpoint->efPoint != 0);

            auto it = toOptimize.begin();
            std::advance(it, k);
            toOptimize.erase(it);
        } else if (ph.lastTraceStatus == IPS_OOB) {
            auto it = toOptimize.begin();
            std::advance(it, k);
            toOptimize.erase(it);

            it = immaturePoints[k].begin();
            std::advance(it, ph.idxInImmaturePoints);
            immaturePoints[k].erase(it);
        } else {
            assert(newpoint == 0 || newpoint == (PointHessian *) ((long) (-1)));
        }

    }


    for (int i = 0; i < immaturePoints.size(); ++i) {
        for (int j = 0; i < immaturePoints[i].size(); ++j) {
            if (immaturePoints[i][j].remove) {
                auto it = immaturePoints[i].begin();
                std::advance(it, j);
                immaturePoints[i].erase(it);
            }
        }
    }
}

void Tracker::setPrecalcValues() {
    for (auto &fh: frameHessians) {
        fh->targetPrecalc.resize(frameHessians.size());

        for (unsigned int i = 0; i < frameHessians.size(); i++)
            fh->targetPrecalc[i].set(fh->pose, frameHessians[i]->pose, hCalib.fxl(), hCalib.fyl(), hCalib.cxl(),
                                     hCalib.cyl());
    }

    ef.setDeltaF(&hCalib);
}


void Tracker::activatePointsMT_Reductor(
        std::vector<PointHessian> *optimized,
        std::vector<ImmaturePoint> *toOptimize,
        int min, int max, Vec10 *stats, int tid) {
    ImmaturePointTemporaryResidual *tr = new ImmaturePointTemporaryResidual[frameHessians.size()];
    for (int k = min; k < max; k++) {
        optimizeImmaturePoint(&(*toOptimize)[k], 1, tr, &(*optimized)[k]);
    }
    delete[] tr;
}

void Tracker::optimizeImmaturePoint(ImmaturePoint *point, int minObs, ImmaturePointTemporaryResidual *residuals,
                                    PointHessian *phOut) {
    int nres = 0;
    std::vector<std::shared_ptr<VO::Frame>> targetFrames;

    for (auto &fh: frameHessians) {
        if (fh != point->host) {
            residuals[nres].state_NewEnergy = residuals[nres].state_energy = 0;
            residuals[nres].state_NewState = ResState::OUTLIER;
            residuals[nres].state_state = ResState::IN;
            residuals[nres].trackingID = fh->trackingID;
            targetFrames.push_back(fh);
            nres++;
        }
    }
    assert(nres == ((int) frameHessians.size()) - 1);

    bool print = true;//rand()%50==0;

    float lastEnergy = 0;
    float lastHdd = 0;
    float lastbd = 0;
    float currentIdepth = (point->idepth_max + point->idepth_min) * 0.5f;

    PointHessian p;

    for (int i = 0; i < nres; i++) {
        lastEnergy += point->linearizeResidual(&hCalib, 1000, residuals + i, lastHdd, lastbd, currentIdepth,
                                               targetFrames[i]);
        residuals[i].state_state = residuals[i].state_NewState;
        residuals[i].state_energy = residuals[i].state_NewEnergy;
    }

    if (!std::isfinite(lastEnergy) || lastHdd < setting_minIdepthH_act) {
        if (print)
            printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
                   nres, lastHdd, lastEnergy);
        return;
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
            newEnergy += point->linearizeResidual(&hCalib, 1, residuals + i, newHdd, newbd, newIdepth, targetFrames[i]);

        if (!std::isfinite(lastEnergy) || newHdd < setting_minIdepthH_act) {
            if (print)
                printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
                       nres,
                       newHdd,
                       lastEnergy);
            return;
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
        return;        // yeah I'm like 99% sure this is OK on 32bit systems.
    }


    int numGoodRes = 0;
    for (int i = 0; i < nres; i++)
        if (residuals[i].state_state == ResState::IN) numGoodRes++;

    if (numGoodRes < minObs) {
        if (print) printf("OptPoint: OUTLIER!\n");
        return;        // yeah I'm like 99% sure this is OK on 32bit systems.
    }


    VO::pointHessianFromImmaturePoint(point, &p);

    if (!std::isfinite(p.energyTH)) { return; }

    p.lastResiduals[0].second = ResState::OOB;
    p.lastResiduals[1].second = ResState::OOB;
    p.setIdepthZero(currentIdepth);
    p.setIdepth(currentIdepth);
    p.setPointStatus(PointHessian::ACTIVE);

    for (int i = 0; i < nres; i++)
        if (residuals[i].state_state == ResState::IN) {
            PointFrameResidual r(p.pointIndex);
            r.state_NewEnergy = r.state_energy = 0;
            r.state_NewState = ResState::OUTLIER;
            r.setState(ResState::IN);
            p.residuals.push_back(r);

            if (r.trackingID == frameHessians.back()->trackingID) {
                *p.lastResiduals[0].first = r;
                p.lastResiduals[0].second = ResState::IN;
            } else if (r.trackingID == frameHessians[frameHessians.size() - 2]->trackingID) {
                *p.lastResiduals[1].first = r;
                p.lastResiduals[1].second = ResState::IN;
            }
        }

    if (print) printf("point activated!\n");
    *phOut = p;
}