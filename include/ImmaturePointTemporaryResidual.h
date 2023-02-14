//
// Created by magnus on 2/13/23.
//

#ifndef VGIS10_IMMATUREPOINTTEMPORARYRESIDUAL_H
#define VGIS10_IMMATUREPOINTTEMPORARYRESIDUAL_H

#include "Util/Enums.h"

struct ImmaturePointTemporaryResidual
{
public:
    ResState state_state;
    double state_energy;
    ResState state_NewState;
    double state_NewEnergy;
    int trackingID = -1;
};

#endif //VGIS10_IMMATUREPOINTTEMPORARYRESIDUAL_H
