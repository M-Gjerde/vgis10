//
// Created by magnus on 2/14/23.
//

#ifndef VGIS10_ENUMS_H
#define VGIS10_ENUMS_H

enum ResLocation {ACTIVE=0, LINEARIZED, MARGINALIZED, NONE};

enum ResState {IN=0, OOB, OUTLIER};

enum ImmaturePointStatus {
    IPS_GOOD = 0,                    // traced well and good
    IPS_OOB,                    // OOB: end tracking & marginalize!
    IPS_OUTLIER,                // energy too high: if happens again: outlier!
    IPS_SKIPPED,                // traced well and good (but not actually traced).
    IPS_BADCONDITION,            // not traced because of bad condition.
    IPS_UNINITIALIZED            // not even traced once.
};

#endif //VGIS10_ENUMS_H
