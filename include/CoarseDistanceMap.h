//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_COARSEDISTANCEMAP_H
#define VGIS10_COARSEDISTANCEMAP_H


#include "Types.h"
#include "Frame.h"
#include "CalibHessian.h"
#include "Optimized/PointFrameResidual.h"

class CoarseDistanceMap {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CoarseDistanceMap(int ww, int hh, const CameraCalibration* cal){
        fwdWarpedIDDistFinal = new float[ww*hh/4];

        bfsList1 = new Eigen::Vector2i[ww*hh/4];
        bfsList2 = new Eigen::Vector2i[ww*hh/4];

        int fac = 1 << (cal->pyrLevelsUsed-1);


        coarseProjectionGrid = new PointFrameResidual*[2048*(ww*hh/(fac*fac))];
        coarseProjectionGridNum = new int[ww*hh/(fac*fac)];
        calibration = cal;
        w[0]=h[0]=0;
    }
    ~CoarseDistanceMap(){
        delete[] coarseProjectionGrid;
        delete[] coarseProjectionGridNum;
        delete[] fwdWarpedIDDistFinal;
        delete[] bfsList1;
        delete[] bfsList2;
    }

    void makeDistanceMap(
            std::vector<std::shared_ptr<VO::Frame>> allFrames,
            std::shared_ptr<VO::Frame> frame){
        int w1 = w[1];
        int h1 = h[1];
        int wh1 = w1*h1;
        for(int i=0;i<wh1;i++)
            fwdWarpedIDDistFinal[i] = 1000;


        // make coarse tracking templates for latstRef.
        int numItems = 0;

        for(int i = 0; i < allFrames.size(); ++i)
        {
            if(frame == allFrames[i]) continue;
            SE3 fhToNew = frame->pose.PRE_worldToCam * allFrames[i]->pose.PRE_camToWorld;
            Mat33f KRKi = (K[1] * fhToNew.rotationMatrix().cast<float>() * Ki[0]);
            Vec3f Kt = (K[1] * fhToNew.translation().cast<float>());

            for(auto& ph : allFrames[i]->pointHessians)
            {
                assert(ph.status == PointHessian::ACTIVE);
                Vec3f ptp = KRKi * Vec3f(ph.u, ph.v, 1) + Kt*ph.idepth_scaled;
                int u = ptp[0] / ptp[2] + 0.5f;
                int v = ptp[1] / ptp[2] + 0.5f;
                if(!(u > 0 && v > 0 && u < w[1] && v < h[1])) continue;
                fwdWarpedIDDistFinal[u+w1*v]=0;
                bfsList1[numItems] = Eigen::Vector2i(u,v);
                numItems++;
            }
        }

        growDistBFS(numItems);
    }

    void makeInlierVotes(
            std::vector<std::shared_ptr<VO::Frame>> frameHessians);

    void makeK( CalibHessian* HCalib){

        w[0] = calibration->wG[0];
        h[0] = calibration->hG[0];

        fx[0] = HCalib->fxl();
        fy[0] = HCalib->fyl();
        cx[0] = HCalib->cxl();
        cy[0] = HCalib->cyl();

        for (int level = 1; level < calibration->pyrLevelsUsed; ++ level)
        {
            w[level] = w[0] >> level;
            h[level] = h[0] >> level;
            fx[level] = fx[level-1] * 0.5;
            fy[level] = fy[level-1] * 0.5;
            cx[level] = (cx[0] + 0.5) / ((int)1<<level) - 0.5;
            cy[level] = (cy[0] + 0.5) / ((int)1<<level) - 0.5;
        }

        for (int level = 0; level < calibration->pyrLevelsUsed; ++ level)
        {
            K[level]  << fx[level], 0.0, cx[level], 0.0, fy[level], cy[level], 0.0, 0.0, 1.0;
            Ki[level] = K[level].inverse();
            fxi[level] = Ki[level](0,0);
            fyi[level] = Ki[level](1,1);
            cxi[level] = Ki[level](0,2);
            cyi[level] = Ki[level](1,2);
        }
    }


    float* fwdWarpedIDDistFinal;

    Mat33f K[PYR_LEVELS];
    Mat33f Ki[PYR_LEVELS];
    float fx[PYR_LEVELS];
    float fy[PYR_LEVELS];
    float fxi[PYR_LEVELS];
    float fyi[PYR_LEVELS];
    float cx[PYR_LEVELS];
    float cy[PYR_LEVELS];
    float cxi[PYR_LEVELS];
    float cyi[PYR_LEVELS];
    int w[PYR_LEVELS];
    int h[PYR_LEVELS];

    void addIntoDistFinal(int u, int v){
        if(w[0] == 0) return;
        bfsList1[0] = Eigen::Vector2i(u,v);
        fwdWarpedIDDistFinal[u+w[1]*v] = 0;
        growDistBFS(1);
    }


private:

    PointFrameResidual** coarseProjectionGrid;
    int* coarseProjectionGridNum;
    Eigen::Vector2i* bfsList1;
    Eigen::Vector2i* bfsList2;
    const CameraCalibration* calibration;

    void growDistBFS(int bfsNum){
        assert(w[0] != 0);
        int w1 = w[1], h1 = h[1];
        for(int k=1;k<40;k++)
        {
            int bfsNum2 = bfsNum;
            std::swap<Eigen::Vector2i*>(bfsList1,bfsList2);
            bfsNum=0;

            if(k%2==0)
            {
                for(int i=0;i<bfsNum2;i++)
                {
                    int x = bfsList2[i][0];
                    int y = bfsList2[i][1];
                    if(x==0 || y== 0 || x==w1-1 || y==h1-1) continue;
                    int idx = x + y * w1;

                    if(fwdWarpedIDDistFinal[idx+1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x+1,y); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x-1,y); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx+w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x,y+1); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x,y-1); bfsNum++;
                    }
                }
            }
            else
            {
                for(int i=0;i<bfsNum2;i++)
                {
                    int x = bfsList2[i][0];
                    int y = bfsList2[i][1];
                    if(x==0 || y== 0 || x==w1-1 || y==h1-1) continue;
                    int idx = x + y * w1;

                    if(fwdWarpedIDDistFinal[idx+1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x+1,y); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x-1,y); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx+w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x,y+1); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x,y-1); bfsNum++;
                    }

                    if(fwdWarpedIDDistFinal[idx+1+w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+1+w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x+1,y+1); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-1+w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-1+w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x-1,y+1); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx-1-w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx-1-w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x-1,y-1); bfsNum++;
                    }
                    if(fwdWarpedIDDistFinal[idx+1-w1] > k)
                    {
                        fwdWarpedIDDistFinal[idx+1-w1] = k;
                        bfsList1[bfsNum] = Eigen::Vector2i(x+1,y-1); bfsNum++;
                    }
                }
            }
        }
    }
};


#endif //VGIS10_COARSEDISTANCEMAP_H
