//
// Created by magnus on 2/11/23.
//

#ifndef VGIS10_SETTINGS_H
#define VGIS10_SETTINGS_H

static int staticPattern[8][2] = {
        {0,  -2},
        {-1, -1},
        {1,  -1},
        {-2, 0},
        {0,  0},
        {2,  0},
        {-1, 1},
        {0,  2}};

#define PYR_LEVELS 5

// parameters controlling pixel selection
#define setting_minGradHistCut 0.5f
#define setting_minGradHistAdd 7
#define setting_gradDownweightPerLevel 0.75
#define setting_selectDirectionDistribution true

#define setting_desiredPointDensity  2000 // aimed total points in the active window.

#define setting_gradient_block_32 32
#define setting_outlierTH 12*12                    // higher -> less strict
#define setting_outlierTHSumComponent 50*50        // higher -> less strong gradient-based reweighting .
#define setting_huberTH 9


// Parameters controlling adaptive energy threshold computation
#define setting_overallEnergyTHWeight 1
#define setting_frameEnergyTHConstWeight 0.5
#define setting_frameEnergyTHN 0.7f
#define setting_frameEnergyTHFacMedian 1.5
#define setting_coarseCutoffTH 20

/* initial hessian values to fix unobservable dimensions / priors on affine lighting parameters.
 */
#define setting_idepthFixPrior 50*50
#define setting_idepthFixPriorMargFac 600*600
#define setting_initialRotPrior 1e11
#define setting_initialTransPrior 1e10
#define setting_initialAffBPrior 1e14
#define setting_initialAffAPrior 1e14
#define setting_initialCalibHessian 5e9

#define setting_affineOptModeA 1e12 //-1: fix. >=0: optimize (with prior, if > 0).
#define setting_affineOptModeB 1e8  //-1: fix. >=0: optimize (with prior, if > 0).


/* require some minimum number of residuals for a point to become valid */
#define   setting_minGoodActiveResForMarg 3
#define   setting_minGoodResForMarg 4

#define SOLVER_REMOVE_POSEPRIOR (int)32
#define SOLVER_FIX_LAMBDA (int)128
#define SOLVER_ORTHOGONALIZE_X_LATER (int)2048
#define setting_solverMode SOLVER_FIX_LAMBDA | SOLVER_ORTHOGONALIZE_X_LATER

#define patternPadding 2
#define patternNum 8
#define patternP staticPattern

/* settings controling initial immature point tracking */
#define setting_maxPixSearch  0.027                 // max length of the ep. line segment searched during immature point tracking. relative to image resolution.
#define setting_minTraceQuality  3
#define setting_minTraceTestRadius  2
#define setting_GNItsOnPointActivation  3
#define setting_trace_stepsize  1.0                 // stepsize for initial discrete search.
#define setting_trace_GNIterations  3               // max # GN iterations
#define setting_trace_GNThreshold  0.1              // GN stop after this stepsize.
#define setting_trace_extraSlackOnTH  1.2           // for energy-based outlier check, be slightly more relaxed by this factor.
#define setting_trace_slackInterval  1.5            // if pixel-interval is smaller than this, leave it be.
#define setting_trace_minImprovementFactor  2       // if pixel-interval is smaller than this, leave it be.

// Settings for keyframes
#define setting_minFrames   5 // min frames in window.
#define setting_maxFrames   7 // max frames in window.
#define setting_minFrameAge   1
#define setting_maxOptIterations 6 // max GN iterations.
#define setting_minOptIterations 1 // min GN iterations.
#define setting_thOptIterations 1.2 // factor on break threshold for GN iteration (larger = break earlier)

#define setting_desiredImmatureDensity 1500 // immature points per frame
#define setting_desiredPointDensity 2000 // aimed total points in the active window.
#define setting_minPointsRemaining 0.05  // marg a frame if less than X% points remain.
#define setting_maxLogAffFacInWindow 0.7 // marg a frame if factor between intensities to current frame is larger than 1/X or X.

/* some thresholds on when to activate / marginalize points */
#define setting_minIdepthH_act 100
#define setting_minIdepthH_marg 50

#endif //VGIS10_SETTINGS_H
