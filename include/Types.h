#ifndef VGIS10_TYPES_H
#define VGIS10_TYPES_H

#include "Eigen/Core"
#include "sophus/sim3.hpp"
#include "sophus/se3.hpp"

#define CPARS 4
#define MAX_RES_PER_POINT 8

#define SCALE_IDEPTH 1.0f		// scales internal value to idepth.
#define SCALE_XI_ROT 1.0f
#define SCALE_XI_TRANS 0.5f
#define SCALE_F 50.0f
#define SCALE_C 50.0f
#define SCALE_W 1.0f
#define SCALE_A 10.0f
#define SCALE_B 1000.0f
#define todouble(x) (x).cast<double>()

#define SCALE_IDEPTH_INVERSE (1.0f / SCALE_IDEPTH)
#define SCALE_XI_ROT_INVERSE (1.0f / SCALE_XI_ROT)
#define SCALE_XI_TRANS_INVERSE (1.0f / SCALE_XI_TRANS)
#define SCALE_F_INVERSE (1.0f / SCALE_F)
#define SCALE_C_INVERSE (1.0f / SCALE_C)
#define SCALE_W_INVERSE (1.0f / SCALE_W)
#define SCALE_A_INVERSE (1.0f / SCALE_A)
#define SCALE_B_INVERSE (1.0f / SCALE_B)


typedef Sophus::SE3d SE3;
typedef Sophus::Sim3d Sim3;
typedef Sophus::SO3d SO3;


typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatXX;
typedef Eigen::Matrix<double,CPARS,CPARS> MatCC;
#define MatToDynamic(x) MatXX(x)

typedef Eigen::Matrix<double,CPARS,10> MatC10;
typedef Eigen::Matrix<double,10,10> Mat1010;
typedef Eigen::Matrix<double,13,13> Mat1313;

typedef Eigen::Matrix<double,8,10> Mat810;
typedef Eigen::Matrix<double,8,3> Mat83;
typedef Eigen::Matrix<double,6,6> Mat66;
typedef Eigen::Matrix<double,5,3> Mat53;
typedef Eigen::Matrix<double,4,3> Mat43;
typedef Eigen::Matrix<double,4,2> Mat42;
typedef Eigen::Matrix<double,3,3> Mat33;
typedef Eigen::Matrix<double,2,2> Mat22;
typedef Eigen::Matrix<double,8,CPARS> Mat8C;
typedef Eigen::Matrix<double,CPARS,8> MatC8;
typedef Eigen::Matrix<float,8,CPARS> Mat8Cf;
typedef Eigen::Matrix<float,CPARS,8> MatC8f;

typedef Eigen::Matrix<double,8,8> Mat88;
typedef Eigen::Matrix<double,7,7> Mat77;

typedef Eigen::Matrix<double,CPARS,1> VecC;
typedef Eigen::Matrix<float,CPARS,1> VecCf;
typedef Eigen::Matrix<double,13,1> Vec13;
typedef Eigen::Matrix<double,10,1> Vec10;
typedef Eigen::Matrix<double,9,1> Vec9;
typedef Eigen::Matrix<double,8,1> Vec8;
typedef Eigen::Matrix<double,7,1> Vec7;
typedef Eigen::Matrix<double,6,1> Vec6;
typedef Eigen::Matrix<double,5,1> Vec5;
typedef Eigen::Matrix<double,4,1> Vec4;
typedef Eigen::Matrix<double,3,1> Vec3;
typedef Eigen::Matrix<double,2,1> Vec2;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;

typedef Eigen::Matrix<float,3,3> Mat33f;
typedef Eigen::Matrix<float,10,3> Mat103f;
typedef Eigen::Matrix<float,2,2> Mat22f;
typedef Eigen::Matrix<float,3,1> Vec3f;
typedef Eigen::Matrix<float,2,1> Vec2f;
typedef Eigen::Matrix<float,6,1> Vec6f;



typedef Eigen::Matrix<double,4,9> Mat49;
typedef Eigen::Matrix<double,8,9> Mat89;

typedef Eigen::Matrix<double,9,4> Mat94;
typedef Eigen::Matrix<double,9,8> Mat98;

typedef Eigen::Matrix<double,8,1> Mat81;
typedef Eigen::Matrix<double,1,8> Mat18;
typedef Eigen::Matrix<double,9,1> Mat91;
typedef Eigen::Matrix<double,1,9> Mat19;


typedef Eigen::Matrix<double,8,4> Mat84;
typedef Eigen::Matrix<double,4,8> Mat48;
typedef Eigen::Matrix<double,4,4> Mat44;


typedef Eigen::Matrix<float,MAX_RES_PER_POINT,1> VecNRf;
typedef Eigen::Matrix<float,12,1> Vec12f;
typedef Eigen::Matrix<float,1,8> Mat18f;
typedef Eigen::Matrix<float,6,6> Mat66f;
typedef Eigen::Matrix<float,8,8> Mat88f;
typedef Eigen::Matrix<float,8,4> Mat84f;
typedef Eigen::Matrix<float,8,1> Vec8f;
typedef Eigen::Matrix<float,10,1> Vec10f;
typedef Eigen::Matrix<float,6,6> Mat66f;
typedef Eigen::Matrix<float,4,1> Vec4f;
typedef Eigen::Matrix<float,4,4> Mat44f;
typedef Eigen::Matrix<float,12,12> Mat1212f;
typedef Eigen::Matrix<float,12,1> Vec12f;
typedef Eigen::Matrix<float,13,13> Mat1313f;
typedef Eigen::Matrix<float,10,10> Mat1010f;
typedef Eigen::Matrix<float,13,1> Vec13f;
typedef Eigen::Matrix<float,9,9> Mat99f;
typedef Eigen::Matrix<float,9,1> Vec9f;

typedef Eigen::Matrix<float,4,2> Mat42f;
typedef Eigen::Matrix<float,6,2> Mat62f;
typedef Eigen::Matrix<float,1,2> Mat12f;

typedef Eigen::Matrix<float,Eigen::Dynamic,1> VecXf;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> MatXXf;


typedef Eigen::Matrix<double,8+CPARS+1,8+CPARS+1> MatPCPC;
typedef Eigen::Matrix<float,8+CPARS+1,8+CPARS+1> MatPCPCf;
typedef Eigen::Matrix<double,8+CPARS+1,1> VecPC;
typedef Eigen::Matrix<float,8+CPARS+1,1> VecPCf;

typedef Eigen::Matrix<float,14,14> Mat1414f;
typedef Eigen::Matrix<float,14,1> Vec14f;
typedef Eigen::Matrix<double,14,14> Mat1414;
typedef Eigen::Matrix<double,14,1> Vec14;




// transforms points from one frame to another.
struct AffLight
{
    AffLight(double a_, double b_) : a(a_), b(b_) {};
    AffLight() : a(0), b(0) {};

    // Affine Parameters:
    double a,b;	// I_frame = exp(a)*I_global + b. // I_global = exp(-a)*(I_frame - b).

    static Vec2 fromToVecExposure(float exposureF, float exposureT, AffLight g2F, AffLight g2T)
    {
        if(exposureF==0 || exposureT==0)
        {
            exposureT = exposureF = 1;
            //printf("got exposure value of 0! please choose the correct model.\n");
            //assert(setting_brightnessTransferFunc < 2);
        }

        double a = exp(g2T.a-g2F.a) * exposureT / exposureF;
        double b = g2T.b - a*g2F.b;
        return Vec2(a,b);
    }

    Vec2 vec()
    {
        return Vec2(a,b);
    }
};

#endif //VGIS10_TYPES_H
