//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_EFFRAME_H
#define VGIS10_EFFRAME_H

#include "EFPoint.h"
#include "Frame.h"

class EFFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFFrame(std::shared_ptr<VO::Frame> frameIn)
    {

        prior =frameIn->pose.getPrior().head<8>();
        delta = frameIn->pose.get_state_minus_stateZero().head<8>();
        delta_prior =  (frameIn->pose.get_state() - frameIn->pose.getPriorZero()).head<8>();
        frameID = frameIn->trackingID;
        data = frameIn;
    }

    std::shared_ptr<VO::Frame> data;

    Vec8 prior;				// prior hessian (diagonal)
    Vec8 delta_prior;		// = state-state_prior (E_prior = (delta_prior)' * diag(prior) * (delta_prior)
    Vec8 delta;				// state - state_zero.


    std::vector<EFPoint> points;
    int idx;	// idx in frames.
    int frameID;
};

#endif //VGIS10_EFFRAME_H
