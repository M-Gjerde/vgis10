//
// Created by magnus on 2/6/23.
//

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "PointSelection.h"
#include "CameraCalibration.h"


namespace VO {

    void initializeCandidatePoints(std::shared_ptr<Frame> frame, const CandidatePointInfo &info,
                                   const CameraCalibration *calibration,
                                   PntResult *result) {
        PointSelection sel(frame->width, frame->height);
        auto *statusMap = new float[frame->width * frame->height];
        bool *statusMapB = new bool[frame->width * frame->height];


        for(int lvl=0; lvl<calibration->pyrLevelsUsed; lvl++)
        {
            result->pointsLevel[lvl].clear();
        }


        float densities[] = {0.03, 0.05, 0.15, 0.5, 1};
        for (int lvl = 0; lvl < calibration->pyrLevelsUsed; lvl++) {
            sel.currentPotential = 3;
            int npts;
            if (lvl == 0)
                npts = sel.makeMaps(frame.get(), statusMap, densities[lvl] * frame->width * frame->height, 1, false,
                                    2);
            else
                npts = makePixelStatus(frame->pyramidNoPhotometric[lvl].data(), statusMapB, frame->widthLevel[lvl], frame->heightLevel[lvl],
                                       densities[lvl] * frame->width * frame->height);

            result->pointsLevel[lvl].resize(npts);

            // set idepth map to initially 1 everywhere.
            int wl = frame->widthLevel[lvl], hl = frame->heightLevel[lvl];
            Pnt *pl = result->pointsLevel[lvl].data();
            int nl = 0;
            for (int y = patternPadding + 1; y < hl - patternPadding - 2; y++)
                for (int x = patternPadding + 1; x < wl - patternPadding - 2; x++) {
                    //if(x==2) printf("y=%d!\n",y);
                    if ((lvl != 0 && statusMapB[x + y * wl]) || (lvl == 0 && statusMap[x + y * wl] != 0)) {
                        //assert(patternNum==9);
                        pl[nl].u = x+0.1;
                        pl[nl].v = y+0.1;
                        pl[nl].idepth = 1;
                        pl[nl].iR = 1;
                        pl[nl].isGood=true;
                        pl[nl].energy.setZero();
                        pl[nl].lastHessian=0;
                        pl[nl].lastHessian_new=0;
                        pl[nl].my_type= (lvl!=0) ? 1 : statusMap[x+y*wl];

                        Eigen::Vector3f *cpt = frame->pyramid[lvl].data() + x + y * frame->widthLevel[lvl];
                        float sumGrad2 = 0;
                        for (int idx = 0; idx < patternNum; idx++) {
                            int dx = patternP[idx][0];
                            int dy = patternP[idx][1];
                            float absgrad = cpt[dx + dy * frame->widthLevel[lvl]].tail<2>().squaredNorm();
                            sumGrad2 += absgrad;
                        }

//				float gth = setting_outlierTH * (sqrtf(sumGrad2)+setting_outlierTHSumComponent);
//				pl[nl].outlierTH = patternNum*gth*gth;
//

                        pl[nl].outlierTH = patternNum * setting_outlierTH;


                        nl++;
                        assert(nl <= npts);
                    }
                }


            result->numPointsLevel[lvl] = nl;
        }
        delete[] statusMap;
        delete[] statusMapB;
    }


    int
    PointSelection::makeMaps(Frame *frame, float *map_out, float density, int recursionsLeft, bool plot,
                             float thFactor) {
        float numHave = 0;
        float numWant = density;
        float quotia;

        int idealPotential = currentPotential;
        {

            // the number of selected pixels behaves approximately as
            // K / (pot+1)^2, where K is a scene-dependent constant.
            // we will allow sub-selecting pixels by up to a quotia of 0.25, otherwise we will re-select.

            if (frame->histogram.empty())
                makeHists(frame);
            // select!
            Eigen::Vector3i n = select(frame, map_out, currentPotential, thFactor, randomPattern);
            // sub-select!
            numHave = n[0] + n[1] + n[2];
            quotia = numWant / numHave;
            // by default we want to over-sample by 40% just to be sure.
            float K = numHave * (currentPotential + 1) * (currentPotential + 1);
            idealPotential = sqrtf(K / numWant) - 1;    // round down.
            if (idealPotential < 1)
                idealPotential = 1;
            if (recursionsLeft > 0 && quotia > 1.25 && currentPotential > 1) {
                //re-sample to get more points!
                // potential needs to be smaller
                if (idealPotential >= currentPotential)
                    idealPotential = currentPotential - 1;

                Log::Logger::getInstance()->info(
                        "PixelSelector: have {:.2f}%%, need {:.2f}%%. RESAMPLE with pot {} -> {}.",
                        100 * numHave / (float) (frame->width * frame->height),
                        100 * numWant / (float) (frame->width * frame->height),
                        currentPotential,
                        idealPotential);
                currentPotential = idealPotential;
                return makeMaps(frame, map_out, density, recursionsLeft - 1, plot, thFactor);
            } else if (recursionsLeft > 0 && quotia < 0.25) {
                //re-sample to get less points!
                if (idealPotential <= currentPotential)
                    idealPotential = currentPotential + 1;
                Log::Logger::getInstance()->info(
                        "PixelSelector: have {:.2f}%%, need {:.2f}%%. RESAMPLE with pot {} -> {}",
                        100 * numHave / (float) (frame->width * frame->height),
                        100 * numWant / (float) (frame->width * frame->height),
                        currentPotential,
                        idealPotential);
                currentPotential = idealPotential;
                return makeMaps(frame, map_out, density, recursionsLeft - 1, plot, thFactor);
            }
        }
        int numHaveSub = numHave;
        if (quotia < 0.95) {
            int wh = frame->width * frame->height;
            int rn = 0;
            unsigned char charTH = 255 * quotia;
            for (int i = 0; i < wh; i++) {
                if (map_out[i] != 0) {
                    if (randomPattern[rn] > charTH) {
                        map_out[i] = 0;
                        numHaveSub--;
                    }
                    rn++;
                }
            }
        }
        currentPotential = idealPotential;
        int w = frame->width;
        int h = frame->height;

        cv::Mat img(h, w, CV_32FC3, frame->pyramid[0].data());

        img.convertTo(img, CV_8UC3);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int i = x + y * w;
                if (map_out[i] == 1)
                    cv::circle(img, cv::Point(x, y), 1, cv::Vec3i(0, 255, 0), 1);
                else if (map_out[i] == 2)
                    cv::circle(img, cv::Point(x, y), 1, cv::Vec3i(255, 0, 0), 1);
                else if (map_out[i] == 4)
                    cv::circle(img, cv::Point(x, y), 1, cv::Vec3i(0, 0, 255), 1);
            }
        }
        cv::imshow("Point selection", img);
        //cv::waitKey(0);
        return numHaveSub;
    }


    int computeHistQuantil(int *hist, float below) {
        int th = hist[0]*below+0.5f;
        for(int i=0;i<90;i++)
        {
            th -= hist[i+1];
            if(th<0) return i;
        }
        return 90;
    }


    void makeHists(Frame *frame) {
        float *mapmax0 = frame->absSquaredGrad[0].data();

        int w = frame->width;
        int h = frame->height;

        frame->gradientHistogram.resize(100 * (1 + w / 32) * (1 + h / 32));
        frame->histogramSmoothed.resize((w / 32) * (h / 32) + 100);
        frame->histogram.resize((w / 32) * (h / 32) + 100);

        int w32 = w / 32;
        int h32 = h / 32;
        frame->histStep = w32;

        for (int y = 0; y < h32; y++) {
            for (int x = 0; x < w32; x++) {
                const float *map0 = mapmax0 + 32 * x + 32 * y * w;
                int *hist0 = frame->gradientHistogram.data();// + 50*(x+y*w32);
                memset(hist0, 0, sizeof(int) * 50);

                for (int j = 0; j < 32; j++)
                    for (int i = 0; i < 32; i++) {
                        int it = i + 32 * x;
                        int jt = j + 32 * y;
                        if (it > w - 2 || jt > h - 2 || it < 1 || jt < 1) continue;
                        int g = sqrtf(map0[i + j * w]);
                        if (g > 48) g = 48;
                        hist0[g + 1]++;
                        hist0[0]++;
                    }

                frame->histogram[x + y * w32] = computeHistQuantil(hist0, setting_minGradHistCut) + setting_minGradHistAdd;
            }
        }
        for (int y = 0; y < h32; y++) {
            for (int x = 0; x < w32; x++) {
                float sum = 0, num = 0;
                if (x > 0) {
                    if (y > 0) {
                        num++;
                        sum += frame->histogram[x - 1 + (y - 1) * w32];
                    }
                    if (y < h32 - 1) {
                        num++;
                        sum += frame->histogram[x - 1 + (y + 1) * w32];
                    }
                    num++;
                    sum += frame->histogram[x - 1 + (y) * w32];
                }

                if (x < w32 - 1) {
                    if (y > 0) {
                        num++;
                        sum += frame->histogram[x + 1 + (y - 1) * w32];
                    }
                    if (y < h32 - 1) {
                        num++;
                        sum += frame->histogram[x + 1 + (y + 1) * w32];
                    }
                    num++;
                    sum += frame->histogram[x + 1 + (y) * w32];
                }

                if (y > 0) {
                    num++;
                    sum += frame->histogram[x + (y - 1) * w32];
                }
                if (y < h32 - 1) {
                    num++;
                    sum += frame->histogram[x + (y + 1) * w32];
                }
                num++;
                sum += frame->histogram[x + y * w32];

                frame->histogramSmoothed[x + y * w32] = (sum / num) * (sum / num);

            }
        }
    }

    Eigen::Vector3i
    select(Frame *frame, float *map_out, int pot, float thFactor, const std::vector<unsigned char> &randomPattern) {

        Eigen::Vector3f const *const map0 = frame->pyramidNoPhotometric[0].data();

        float *mapmax0 = frame->absSquaredGrad[0].data();
        float *mapmax1 = frame->absSquaredGrad[1].data();
        float *mapmax2 = frame->absSquaredGrad[2].data();


        int w = frame->widthLevel[0];
        int w1 = frame->widthLevel[1];
        int w2 = frame->widthLevel[2];
        int h = frame->heightLevel[0];


        const Vec2f directions[16] = {
                Vec2f(0, 1.0000),
                Vec2f(0.3827, 0.9239),
                Vec2f(0.1951, 0.9808),
                Vec2f(0.9239, 0.3827),
                Vec2f(0.7071, 0.7071),
                Vec2f(0.3827, -0.9239),
                Vec2f(0.8315, 0.5556),
                Vec2f(0.8315, -0.5556),
                Vec2f(0.5556, -0.8315),
                Vec2f(0.9808, 0.1951),
                Vec2f(0.9239, -0.3827),
                Vec2f(0.7071, -0.7071),
                Vec2f(0.5556, 0.8315),
                Vec2f(0.9808, -0.1951),
                Vec2f(1.0000, 0.0000),
                Vec2f(0.1951, -0.9808)};

        memset(map_out, 0, w * h * sizeof(PixelSelectorStatus));


        float dw1 = setting_gradDownweightPerLevel;
        float dw2 = dw1 * dw1;


        int n3 = 0, n2 = 0, n4 = 0;
        for (int y4 = 0; y4 < h; y4 += (4 * pot))
            for (int x4 = 0; x4 < w; x4 += (4 * pot)) {
                int my3 = std::min((4 * pot), h - y4);
                int mx3 = std::min((4 * pot), w - x4);
                int bestIdx4 = -1;
                float bestVal4 = 0;
                Vec2f dir4 = directions[randomPattern[n2] & 0xF];
                for (int y3 = 0; y3 < my3; y3 += (2 * pot))
                    for (int x3 = 0; x3 < mx3; x3 += (2 * pot)) {
                        int x34 = x3 + x4;
                        int y34 = y3 + y4;
                        int my2 = std::min((2 * pot), h - y34);
                        int mx2 = std::min((2 * pot), w - x34);
                        int bestIdx3 = -1;
                        float bestVal3 = 0;
                        Vec2f dir3 = directions[randomPattern[n2] & 0xF];
                        for (int y2 = 0; y2 < my2; y2 += pot)
                            for (int x2 = 0; x2 < mx2; x2 += pot) {
                                int x234 = x2 + x34;
                                int y234 = y2 + y34;
                                int my1 = std::min(pot, h - y234);
                                int mx1 = std::min(pot, w - x234);
                                int bestIdx2 = -1;
                                float bestVal2 = 0;
                                Vec2f dir2 = directions[randomPattern[n2] & 0xF];
                                for (int y1 = 0; y1 < my1; y1 += 1)
                                    for (int x1 = 0; x1 < mx1; x1 += 1) {
                                        assert(x1 + x234 < w);
                                        assert(y1 + y234 < h);
                                        int idx = x1 + x234 + w * (y1 + y234);
                                        int xf = x1 + x234;
                                        int yf = y1 + y234;

                                        if (xf < 4 || xf >= w - 5 || yf < 4 || yf > h - 4) continue;


                                        float pixelTH0 = frame->histogramSmoothed[(xf >> 5) + (yf >> 5) *
                                                                                              frame->histStep];
                                        float pixelTH1 = pixelTH0 * dw1;
                                        float pixelTH2 = pixelTH1 * dw2;


                                        float ag0 = mapmax0[idx];
                                        if (ag0 > pixelTH0 * thFactor) {
                                            Vec2f ag0d = map0[idx].tail<2>();
                                            float dirNorm = fabsf((float) (ag0d.dot(dir2)));
                                            if (!setting_selectDirectionDistribution) dirNorm = ag0;

                                            if (dirNorm > bestVal2) {
                                                bestVal2 = dirNorm;
                                                bestIdx2 = idx;
                                                bestIdx3 = -2;
                                                bestIdx4 = -2;
                                            }
                                        }
                                        if (bestIdx3 == -2) continue;

                                        float ag1 = mapmax1[(int) (xf * 0.5f + 0.25f) + (int) (yf * 0.5f + 0.25f) * w1];
                                        if (ag1 > pixelTH1 * thFactor) {
                                            Vec2f ag0d = map0[idx].tail<2>();
                                            float dirNorm = fabsf((float) (ag0d.dot(dir3)));
                                            if (!setting_selectDirectionDistribution) dirNorm = ag1;

                                            if (dirNorm > bestVal3) {
                                                bestVal3 = dirNorm;
                                                bestIdx3 = idx;
                                                bestIdx4 = -2;
                                            }
                                        }
                                        if (bestIdx4 == -2) continue;

                                        float ag2 = mapmax2[(int) (xf * 0.25f + 0.125) +
                                                            (int) (yf * 0.25f + 0.125) * w2];
                                        if (ag2 > pixelTH2 * thFactor) {
                                            Vec2f ag0d = map0[idx].tail<2>();
                                            float dirNorm = fabsf((float) (ag0d.dot(dir4)));
                                            if (!setting_selectDirectionDistribution) dirNorm = ag2;

                                            if (dirNorm > bestVal4) {
                                                bestVal4 = dirNorm;
                                                bestIdx4 = idx;
                                            }
                                        }
                                    }

                                if (bestIdx2 > 0) {
                                    map_out[bestIdx2] = 1;
                                    bestVal3 = 1e10;
                                    n2++;
                                }
                            }

                        if (bestIdx3 > 0) {
                            map_out[bestIdx3] = 2;
                            bestVal4 = 1e10;
                            n3++;
                        }
                    }

                if (bestIdx4 > 0) {
                    map_out[bestIdx4] = 4;
                    n4++;
                }
            }


        return Eigen::Vector3i(n2, n3, n4);
    }



    template<int pot>
    int gridMaxSelection(Eigen::Vector3f* grads, bool* map_out, int w, int h, float THFac)
    {
        const float minUseGrad_pixsel = 10;
        memset(map_out, 0, sizeof(bool)*w*h);

        int numGood = 0;
        for(int y=1;y<h-pot;y+=pot)
        {
            for(int x=1;x<w-pot;x+=pot)
            {
                int bestXXID = -1;
                int bestYYID = -1;
                int bestXYID = -1;
                int bestYXID = -1;

                float bestXX=0, bestYY=0, bestXY=0, bestYX=0;

                Eigen::Vector3f* grads0 = grads+x+y*w;
                for(int dx=0;dx<pot;dx++)
                    for(int dy=0;dy<pot;dy++)
                    {
                        int idx = dx+dy*w;
                        Eigen::Vector3f g=grads0[idx];
                        float sqgd = g.tail<2>().squaredNorm();
                        float TH = THFac*minUseGrad_pixsel * (0.75f);

                        if(sqgd > TH*TH)
                        {
                            float agx = fabs((float)g[1]);
                            if(agx > bestXX) {bestXX=agx; bestXXID=idx;}

                            float agy = fabs((float)g[2]);
                            if(agy > bestYY) {bestYY=agy; bestYYID=idx;}

                            float gxpy = fabs((float)(g[1]-g[2]));
                            if(gxpy > bestXY) {bestXY=gxpy; bestXYID=idx;}

                            float gxmy = fabs((float)(g[1]+g[2]));
                            if(gxmy > bestYX) {bestYX=gxmy; bestYXID=idx;}
                        }
                    }

                bool* map0 = map_out+x+y*w;

                if(bestXXID>=0)
                {
                    if(!map0[bestXXID])
                        numGood++;
                    map0[bestXXID] = true;

                }
                if(bestYYID>=0)
                {
                    if(!map0[bestYYID])
                        numGood++;
                    map0[bestYYID] = true;

                }
                if(bestXYID>=0)
                {
                    if(!map0[bestXYID])
                        numGood++;
                    map0[bestXYID] = true;

                }
                if(bestYXID>=0)
                {
                    if(!map0[bestYXID])
                        numGood++;
                    map0[bestYXID] = true;

                }
            }
        }

        return numGood;
    }

    int makePixelStatus(Eigen::Vector3f *grads, bool *map, int w, int h, float desiredDensity, int recsLeft, float THFac) {

        int sparsityFactor = 5;

        if (sparsityFactor < 1) sparsityFactor = 1;

        int numGoodPoints;
        numGoodPoints = gridMaxSelection<5>(grads, map, w, h, THFac);

        /*
         * #points is approximately proportional to sparsityFactor^2.
         */

        float quotia = numGoodPoints / (float) (desiredDensity);

        int newSparsity = (sparsityFactor * sqrtf(quotia)) + 0.7f;


        if (newSparsity < 1) newSparsity = 1;


        float oldTHFac = THFac;
        if (newSparsity == 1 && sparsityFactor == 1) THFac = 0.5;


        if ((abs(newSparsity - sparsityFactor) < 1 && THFac == oldTHFac) ||
            (quotia > 0.8 && 1.0f / quotia > 0.8) ||
            recsLeft == 0) {

//		printf(" \n");
            //all good
            sparsityFactor = newSparsity;
            return numGoodPoints;
        } else {
//		printf(" -> re-evaluate! \n");
            // re-evaluate.
            sparsityFactor = newSparsity;
            return makePixelStatus(grads, map, w, h, desiredDensity, recsLeft - 1, THFac);
        }
    }
};