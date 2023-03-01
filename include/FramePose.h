#include "Types.h"


namespace VO {

    struct FramePose {
        // precalc values
        int id;            // INTERNAL ID, starting at zero.
        int incoming_id;    // ID passed into DSO
        double timestamp;        // timestamp passed into DSO.

        // set once after tracking
        SE3 camToTrackingRef;
        FramePose *trackingRef;

        // constantly adapted.
        SE3 camToWorld;                // Write: TRACKING, while frame is still fresh; MAPPING: only when locked [shellPoseMutex].
        AffLight aff_g2l;
        bool poseValid;

        // statisitcs
        int statistics_outlierResOnThis;
        int statistics_goodResOnThis;
        int marginalizedAt;
        double movedByOpt;
        FramePose()
        {
            id=0;
            poseValid=true;
            camToWorld = SE3();
            timestamp=0;
            marginalizedAt=-1;
            movedByOpt=0;
            statistics_outlierResOnThis=statistics_goodResOnThis=0;
            trackingRef=0;
            camToTrackingRef = SE3();
        }

    };
}