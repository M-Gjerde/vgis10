//
// Created by magnus on 2/9/23.
//

#include "Optimized/Energy.h"

void EFFrame::takeData() {

    prior = data->getPrior().head<8>();
    delta = data->get_state_minus_stateZero().head<8>();
    delta_prior = (data->get_state() - data->getPriorZero()).head<8>();



//	Vec10 state_zero =  data->get_state_zero();
//	state_zero.segment<3>(0) = SCALE_XI_TRANS * state_zero.segment<3>(0);
//	state_zero.segment<3>(3) = SCALE_XI_ROT * state_zero.segment<3>(3);
//	state_zero[6] = SCALE_A * state_zero[6];
//	state_zero[7] = SCALE_B * state_zero[7];
//	state_zero[8] = SCALE_A * state_zero[8];
//	state_zero[9] = SCALE_B * state_zero[9];
//
//	std::cout << "state_zero: " << state_zero.transpose() << "\n";


    assert(data->id != -1);

    frameID = data->id;

}
