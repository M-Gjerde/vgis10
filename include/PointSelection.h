//
// Created by magnus on 2/6/23.
//

#ifndef VGIS10_POINTSELECTION_H
#define VGIS10_POINTSELECTION_H

#include "Frame.h"
#include "CameraCalibration.h"
#include "Point.h"

namespace VO {
    enum PixelSelectorStatus {PIXSEL_VOID=0, PIXSEL_1, PIXSEL_2, PIXSEL_3};

    struct CandidatePointInfo {
        uint32_t width = 32;
        uint32_t height = 32;
        uint32_t gth = 7;
        uint32_t distributionBlocks = 12;
    };





    void initializeCandidatePoints(std::shared_ptr<Frame> frame, const CandidatePointInfo &info,
                                   const CameraCalibration *calibration,
                                   PntResult *result);

    int computeHistQuantil(int* hist, float below);

    void makeHists(Frame* frame);
    Eigen::Vector3i
    select(Frame *frame, float *map_out, int pot, float thFactor, const std::vector<unsigned char> &randomPattern);

    static int makePixelStatus(Eigen::Vector3f *grads, bool *map, int w, int h, float desiredDensity, int recsLeft = 5,float THFac = 1);
    template<int pot>
    int gridMaxSelection(Eigen::Vector3f* grads, bool* map_out, int w, int h, float THFac);
    //// CLASS DEFINITION
    class PointSelection{

    public:
        PointSelection(int w, int h){

            randomPattern.resize(w * h);
            std::srand(3141592);    // want to be deterministic.
            for (int i = 0; i < w * h; i++) randomPattern[i] = rand() & 0xFF;

            currentPotential=3;

        }

        int currentPotential = 3;

        int makeMaps(Frame *frame, float *map_out, float density, int recursionsLeft = 1, bool plot = false, float thFactor = 1);
        std::vector<unsigned char> randomPattern{};

    };
};


#endif //VGIS10_POINTSELECTION_H
