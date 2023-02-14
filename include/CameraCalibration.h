//
// Created by magnus on 1/31/23.
//

#ifndef VGIS10_CAMERACALIBRATION_H
#define VGIS10_CAMERACALIBRATION_H

#include <filesystem>
#include <utility>
#include <vector>
#include <Eigen/Eigen>

#include "Frame.h"
#include "CalibHessian.h"

class CameraCalibration {
public:

    struct VignetteMap{
        uint32_t width = 0, height = 0;
        std::vector<float> d;
        std::vector<float> inv;
    } vignetteMap;

    struct Response{
        Response(){
            gInv.resize(255);
            g.resize(255);
        }
        std::vector<float> g;
        std::vector<float> gInv;
    }responseFunc;


    uint32_t  pyrLevelsUsed = 0;

    CameraCalibration(const std::filesystem::path& vignetteImagePath, const std::filesystem::path& CRFPath) {
        readResponseFunction(CRFPath);
        readVignetteMap(vignetteImagePath);
        if (!readIntrinsicCalibration("./calibration/intrinsics.yml")){
            throw std::runtime_error("Failed to read intrinsic calibration");
        }
        // Use image dimmensions from vignette. It should match the images from the calibrated camera
        pyrLevelsUsed = setGlobalCalib(outWidth, outHeight, K);
    }



    void applyPhotometricCalibration(std::shared_ptr<VO::Frame> *pPtr);
    void setFramePyramidInfo(VO::Frame* frame) const;

    void undistort(std::shared_ptr<VO::Frame> frame, float exposure, double timestamp, float factor) const;



    int wG[PYR_LEVELS]{}, hG[PYR_LEVELS]{};
    float fxG[PYR_LEVELS]{}, fyG[PYR_LEVELS]{}, cxG[PYR_LEVELS]{}, cyG[PYR_LEVELS]{};
    float fxiG[PYR_LEVELS]{}, fyiG[PYR_LEVELS]{}, cxiG[PYR_LEVELS]{}, cyiG[PYR_LEVELS]{};
    Mat33 KG[PYR_LEVELS], KiG[PYR_LEVELS], K;
private:

    uint32_t setGlobalCalib(uint32_t w, uint32_t h, const Mat33 &K);
    void readResponseFunction(const std::filesystem::path& CRFPath);
    void readVignetteMap(const std::filesystem::path& vignetteImagePath);
    bool readIntrinsicCalibration(const std::filesystem::path& calibrationFilePath);


    std::vector<size_t> inputRes;
    std::vector<double> kData ;
    std::vector<double> kDataNew;
    std::vector<size_t> outputRes;

    uint32_t inWidth, inHeight;
    uint32_t outWidth, outHeight;

    std::vector<float> remapX;
    std::vector<float> remapY;

    void makeOptimalK_crop();

    void distortCoordinates(float *in_x, float *in_y, float *out_x, float *out_y, int n) const;

    void processFrame(std::shared_ptr<VO::Frame> frame, float exposure_time, float factor) const;

    void rectify(const float *in_data, float *out_data) const;
};

#endif //VGIS10_CAMERACALIBRATION_H
