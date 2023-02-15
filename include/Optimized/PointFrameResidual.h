#ifndef VGIS10_POINTFRAMERESIDUAL_H
#define VGIS10_POINTFRAMERESIDUAL_H



#include "vector"

#include "Types.h"
#include <iostream>
#include <fstream>
#include "Util/Enums.h"
#include "RawResidualJacobian.h"

class EFResidual;
class PointHessian;
namespace VO {
    class Frame;
}

struct FullJacRowT
    {
        Eigen::Vector2f projectedTo[MAX_RES_PER_POINT];
    };

    class PointFrameResidual
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        explicit PointFrameResidual(int index)
        {
            resetOOB();
            J = new RawResidualJacobian();

            pointIndex = index;
            isNew = true;
        }

        EFResidual* efResidual;
        VO::Frame* host;
        VO::Frame* target;
        PointHessian* pointH;

        ResState state_state;
        double state_energy;
        ResState state_NewState;
        double state_NewEnergy;
        double state_NewEnergyWithOutlier;
        bool isNew;

        void setState(ResState s) {state_state = s;}
        RawResidualJacobian* J;
        Eigen::Vector2f projectedTo[MAX_RES_PER_POINT];
        Vec3f centerProjectedTo;
        int trackingID = 0;
        int pointIndex = 0; // Which point this residual belongs to


        ~PointFrameResidual() {
        };


        void resetOOB()
        {
            state_NewEnergy = state_energy = 0;
            state_NewState = ResState::OUTLIER;

            setState(ResState::IN);
        };
        /*
        void applyRes( bool copyJacobians) {
            if (copyJacobians) {
                if (state_state == ResState::OOB) {
                    assert(!efResidual->isActiveAndIsGoodNEW);
                    return;    // can never go back from OOB
                }
                if (state_NewState == ResState::IN)// && )
                {
                    efResidual->isActiveAndIsGoodNEW=true;
                    efResidual->takeDataF();
                } else {
                    efResidual->isActiveAndIsGoodNEW=false;
                }
            }

            setState(state_NewState);
            state_energy = state_NewEnergy;
        }
         */
    };



#endif //VGIS10_POINTFRAMERESIDUAL_H
