#ifndef VGIS10_POINTFRAMERESIDUAL_H
#define VGIS10_POINTFRAMERESIDUAL_H



#include "vector"

#include "Types.h"
#include <iostream>
#include <fstream>
#include "Util/Enums.h"
#include "RawResidualJacobian.h"

class CalibHessian;

class EFResidual;
class PointHessian;
namespace VO {
    class Frame;
}
class PointHessian;

class EFResidual;


struct FullJacRowT
{
    Eigen::Vector2f projectedTo[MAX_RES_PER_POINT];
};

class PointFrameResidual
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EFResidual* efResidual;

    static int instanceCounter;


    ResState state_state;
    double state_energy;
    ResState state_NewState;
    double state_NewEnergy;
    double state_NewEnergyWithOutlier;


    void setState(ResState s) {state_state = s;}


    PointHessian* point;
    VO::Frame* host;
    VO::Frame* target;
    RawResidualJacobian* J;


    bool isNew;


    Eigen::Vector2f projectedTo[MAX_RES_PER_POINT];
    Vec3f centerProjectedTo;

    ~PointFrameResidual(){assert(efResidual==0); delete J;}
    PointFrameResidual();
    PointFrameResidual(PointHessian* point_, VO::Frame* host_, VO::Frame* target_) :
            point(point_),
            host(host_),
            target(target_)
    {
        efResidual=0;
        resetOOB();
        J = new RawResidualJacobian();
        assert(((long)J)%16==0);

        isNew=true;
    }
    double linearize(CalibHessian* HCalib);


    void resetOOB()
    {
        state_NewEnergy = state_energy = 0;
        state_NewState = ResState::OUTLIER;

        setState(ResState::IN);
    };
    void applyRes(bool copyJacobians);

    void debugPlot();

    void printRows(std::vector<VecX> &v, VecX &r, int nFrames, int nPoints, int M, int res);
};



#endif //VGIS10_POINTFRAMERESIDUAL_H
