//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFFRAME_H
#define VGIS10_EFFRAME_H

#include "Frame.h"
#include "EFPoint.h"

class EFFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFFrame(std::shared_ptr<VO::Frame> d) : data(d)
    {
        takeData();
    }
    void takeData(){
        /*
        prior = data->getPrior().head<8>();
        delta = data->get_state_minus_stateZero().head<8>();
        delta_prior =  (data->get_state() - data->getPriorZero()).head<8>();



//	Vec10 state_zero =  data->get_state_zero();
//	state_zero.segment<3>(0) = SCALE_XI_TRANS * state_zero.segment<3>(0);
//	state_zero.segment<3>(3) = SCALE_XI_ROT * state_zero.segment<3>(3);
//	state_zero[6] = SCALE_A * state_zero[6];
//	state_zero[7] = SCALE_B * state_zero[7];
//	state_zero[8] = SCALE_A * state_zero[8];
//	state_zero[9] = SCALE_B * state_zero[9];
//
//	std::cout << "state_zero: " << state_zero.transpose() << "\n";

*/
        frameID = data->id;
    }


    Vec8 prior;				// prior hessian (diagonal)
    Vec8 delta_prior;		// = state-state_prior (E_prior = (delta_prior)' * diag(prior) * (delta_prior)
    Vec8 delta;				// state - state_zero.



    std::vector<EFPoint*> points;
    std::shared_ptr<VO::Frame> data;
    int idx;	// idx in frames.

    int frameID;
};

#endif //VGIS10_EFFRAME_H
