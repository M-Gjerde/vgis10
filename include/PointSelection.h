//
// Created by magnus on 2/6/23.
//

#ifndef VGIS10_POINTSELECTION_H
#define VGIS10_POINTSELECTION_H

#include "Frame.h"
#include "CameraCalibration.h"

namespace VO {
    enum PixelSelectorStatus {PIXSEL_VOID=0, PIXSEL_1, PIXSEL_2, PIXSEL_3};


    typedef Eigen::Matrix<float,3,3> Mat33f;
    typedef Eigen::Matrix<float,10,3> Mat103f;
    typedef Eigen::Matrix<float,2,2> Mat22f;
    typedef Eigen::Matrix<float,3,1> Vec3f;
    typedef Eigen::Matrix<float,2,1> Vec2f;
    typedef Eigen::Matrix<float,6,1> Vec6f;

    struct CandidatePointInfo {
        uint32_t width = 32;
        uint32_t height = 32;
        uint32_t gth = 7;
        uint32_t distributionBlocks = 12;
    };

    void findCandidatePoints(std::shared_ptr<Frame> frame, const CandidatePointInfo &info,
                             const CameraCalibration *calibration);

    int computeHistQuantil(int* hist, float below);

    void makeHists(Frame* frame);
    Eigen::Vector3i
    select(Frame *frame, float *map_out, int pot, float thFactor, const std::vector<unsigned char> &randomPattern);
    int makeMaps(Frame *frame, float *map_out, float density, int recursionsLeft, bool plot, float thFactor,
                 const std::vector<unsigned char> &randomPattern);


    class PointSelection{

    public:
        PointSelection(int w, int h){

            randomPattern.resize(w * h);
            std::srand(3141592);    // want to be deterministic.
            for (int i = 0; i < w * h; i++) randomPattern[i] = rand() & 0xFF;

            currentPotential=3;

        }

        int currentPotential = 3;

        int makeMaps(Frame *frame, float *map_out, float density, int recursionsLeft, bool plot, float thFactor);
        std::vector<unsigned char> randomPattern{};

    };
};


#endif //VGIS10_POINTSELECTION_H
