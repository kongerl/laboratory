/************全部参数和函数定义*************/

#ifndef LOCATION_ALL_H
#define LOCATION_ALL_H
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/************全部需要常数定义*************/

extern const int num;

//卫星结算常量
const double GM = 3.986004418e14;		//地球引力常数
const double We = 7.2921151467e-5; 		//地球自转速度
const double pi = 3.1415926535898;
const double A_ref_MEO = 2.79061e7;     //MEO
const double A_ref_IG_G = 4.2162200e7;  //IGSO/GEO
const double C = 299792458;             //光速

/************全部需要常数定义*************/

/************全部需要变量定义*************/
//卫星解算变量
double A0_MEO;
double Ak;
double angV;	//1
double Tk;		//2
double Mk;		//3
double Ek[10];		//4
double Vk;		//5
double faik;	//6
//对流层变量
int t;
double uk,rk,ik;
double Xk,Yk;
double doy = 152.0;
double R   = 287.054;
double g   = 9.80665;
double k1  = 77.604;
double k2  = 16.6;
double k3  = 377600.0;
double E;           //卫星高度角

//spp
double delta_tr[8];//相对论
double delta[8];//卫星钟差
/************全部需要变量定义*************/

/**************************结构体定义**************************/
//广播星历
typedef struct{
    char a[3];
    int year, month, day, hour, mintue, second;	//时间参数 minute
    double a_f0, a_f1, a_f2;					//钟差、钟速、钟漂

    double IDOE;								//轨道半径的正弦调和改正项的振幅
    double C_rs;							//参考时刻卫星平均角速度与计算值之差
    double Delta_n0;									//参考时刻的平近点角
    double M0;                                  //?????

    double C_uc;								//纬度幅角的余弦调和改正项的振幅
    double e;									//偏心率
    double C_us;								//纬度幅角的正弦调和改正项的振幅
    double Delta_A;								//长半轴的平方根

    double Toe;									//星历参考时刻
    double C_ic;								//轨道倾角的余弦调和改正项的振幅
    double OMEGA0;								//周历元零时刻计算的升交点经度
    double C_is;								//轨道倾角的正弦调和改正项的振幅

    double i0;									//参考时刻的轨道倾角
    double C_rc;								//轨道半径的余弦调和改正项的振幅
    double omega;								//近地点幅角
    double OMEGA_DOT;							//升交点赤经变化率

    double Idot;								//轨道倾角变化率
    double data_bit;
    double BDT;
    double A_DOT;

    double SV;
    double Health;
    double IGD;
    double ISC;

    double Launch_time;
    double IDOC;
    double Delta_n0_dot;
    double type;

}aaa;
aaa EPH[20];//数组长度为num内存会不够

//卫星解算坐标
typedef struct{
    double Xk;
    double Yk;
    double Zk;
}bbb;
bbb location[8];//命名重复   m级

//卫星伪距 卫星钟差 电离层延时 对流层延时
typedef struct{
    double rou;
    double delta_tu;
    double I_n;
    double T_n;
}ccc;
ccc da[8] = {24688564.000,0,0,0,
             21780791.271,0,0,0,
             25360248.071,0,0,0,
             22833624.564,0,0,0,
             23687638.830,0,0,0,
             37857941.485,0,0,0,
             36562515.356,0,0,0,
             21773519.462,0,0,0,
};
//结果
typedef struct{
    double X;
    double Y;
    double Z;
    double delta_0;
}qqq;
qqq result = {0,0,0,0};


/**************************结构体定义**************************/



int read_data();//读文件

//广播星历计算卫星坐标
double semiMajor_axis(int i);
double angular_Velocity(int i);								//1、平均角速度
double Naturalization_Time(int i);							//2、归化时间
double meanAnomaly(int i);									//3、平近点角
double eccentricAnomanly(int i);							//4、偏近点角Ek
double tureAnomaly(int i);									//5、真近点角
double argumentOfOperigee(int i);						 	//6、升交距角faik
double perturbation(int i);									//7(8)、计算摄动改正量
double coordinates(int i);									//9、计算轨道平面坐标系坐标
double MEOsate(int i);										//10、轨道位置
double GEOsate(int i);

//电离层延迟
double bdgim3(int i);


//对流层延时
double UNB3(int i);
double Matrix();
double yearAverage(double *p);
double amplitude(double *p);
double NMFdry(double *p);
double NMFZWD(double *p);
double calculate(double *p1,double *p2,double d_h,double d_nh);                  //最后计算



//spp
double satellite_clock(int num,double delta[num]);                                    //计算卫星钟差
void matrix_mul(int l,int n,double a[l][n],double b[n][l]);                         //矩阵转置
void matrix_contrary(int l,double a[l][l],double c[l][l]);                          //矩阵求逆
void matrix_ride(int l,int n,int q,double a[l][n],double b[n][q],double c[l][q]);   //矩阵相乘
double derivative(int i,double X,double Y,double Z);                                //求偏导
double count(int num,double median[4][1]);                                          //计算
void estimate(int num,double median[4][1]);                                         //判断收敛



/****************************************************************************************************************
 *   bdgim.h : BDGIM model constants, basic structures and function prototypes.
 *   BDS Broadcast Ionospheric model(BDGIM) calculation module
 *   Author: Zishen LI, Ningbo WANG
 *   e-mail: lizishen@aircas.ac.cn, wangnb@aircas.ac.cn
 *
 *   Copyright (C) 2021 by Aerospace Information Research Institute, Chinese Academy of Sciences, Beijing, China.
 *   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 *   Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 *   option) any later version.
 *
 *	 Create at  2021,	Apr. 23
 ***************************************************************************************************************/



/*********** define some basic consts **************************/

#define PI	          (4.0*atan(1.0))
#define MIN(a,b)   (((a)>(b))?(b):(a))    // minimum between a and b
#define MAX(a,b)  (((a)<(b))?(b):(a))     // maximum between a and b
#define FREQ1_BDS       1575420000.0      // BDS-3 B1C  frequency (Hz)
#define BRDPARANUM      9                 // Number of broadcast ionospheric parameters
#define PERIODNUM      13                 // Number of forecast period for every non-broadcast parameters
#define TRISERINUM    (PERIODNUM*2-1)     // Trigonometric series number
#define NONBRDNUM      17                 // Number of non-broadcast one group
#define MAXGROUP       12                 // 12 groups non-broadcast coefficient every day

/*********** BDGIM Non-Broadcast Ionospheric Parameters Struct **************************/
typedef struct{
    int degOrd[NONBRDNUM][2];		              // spheric harmonic function degree and order index
    double omiga[PERIODNUM];				      // omiga calculated from the period
    double perdTable[NONBRDNUM][PERIODNUM*2-1];	  // the array for storing the non-broadcast parameter period table for perdTable
    double nonBrdCoef[NONBRDNUM][MAXGROUP];	      // the non-broadcast bdgim parameter of the calculate day
} NonBrdIonData;

/*********** BDGIM Broadcast Ionospheric Parameters Struct **************************/
typedef struct{
    double brdIonCoef[BRDPARANUM];		// the body array for storing the bdgim broadcast parameter
    int degOrd[BRDPARANUM][2];		    // spheric harmonic function degree and order index for broadcast parameter
} BrdIonData;

/*-------------------------------- BDGIM model function ----------------------------*/

/* obtains the slant ionospheric delay in B1C using BDGIM ionospheric model */
int IonBdsBrdModel(NonBrdIonData* nonBrdData, BrdIonData* brdData, double mjd, double* sta_xyz, double* sat_xyz, double* brdPara, double* ion_delay);

/* obtains the vertical ionospheric tec using BDGIM ionospheric mode */
int VtecBrdSH(NonBrdIonData* nonBrdData, BrdIonData* brdData, double mjd, double ipp_b, double ipp_l, double* vtec);

/* transform earth-fixed latitude/longitude into sun-fixed latitude/longitude */
void EFLSFL(double mjd, double* lat, double* lon, int geomag, int sunframe, double* lat1, double* lon1);

/* ionospheric mapping functions */
double IonMapping(int type, double ipp_elev, double sat_elev, double Hion);

/* calculate XYZ, latitude, longitude and elevation of the IPP (accurate) according to user and satellite position */
int IPPBLH1(double* sta, double* sat, double Hion, double* IPPXYZ, double* IPP_B, double* IPP_L, double* IPP_E, double* sat_ele);

/* calculate latitude, longitude and elevation of the IPP (approximate) according to user latitude ,longitude and satellite elevation, azimuth */
int IPPBLH2(double lat_u, double lon_u, double hion, double sat_ele, double sat_azimuth, double* ipp_b, double* ipp_l, double* ipp_e);

int BrdCoefGroupIndex(double mjd, NonBrdIonData* nonBrdData);	// obtains non-broadcast coefficient session group according to the mjd
int CalNonBrdCoef(double mjd, NonBrdIonData* nonBrdData);		// calculate the non-broadcast BDGIM coefficients at the first compute epoch
void SetNonBrdCoefPeriod(NonBrdIonData* nonBrdData);		    // Set the period term of the non-broadcast perdTable for BDGIM model
double ASLEFU(double XLAT,double XLON,int INN,int IMM);		    // Normalized legendre polynomial
double FAKULT(int N);										    // compute the factorial of N
double Distance(double *xyz1,double *xyz2);
typedef struct
{
    double mjd;
    double daysec;
} MjdData;

/* calculate MJD from UTC --------------------------------*/
void UTC2MJD(int year, int month, int day, int hour, int min, double second, MjdData* mjdata);




#endif //LOCATION_ALL_H
