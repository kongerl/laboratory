/************ȫ�������ͺ�������*************/

#ifndef LOCATION_ALL_H
#define LOCATION_ALL_H
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/************ȫ����Ҫ��������*************/

extern const int num;

//���ǽ��㳣��
const double GM = 3.986004418e14;		//������������
const double We = 7.2921151467e-5; 		//������ת�ٶ�
const double pi = 3.1415926535898;
const double A_ref_MEO = 2.79061e7;     //MEO
const double A_ref_IG_G = 4.2162200e7;  //IGSO/GEO
const double C = 299792458;             //����

/************ȫ����Ҫ��������*************/

/************ȫ����Ҫ��������*************/
//���ǽ������
double A0_MEO;
double Ak;
double angV;	//1
double Tk;		//2
double Mk;		//3
double Ek[10];		//4
double Vk;		//5
double faik;	//6
//���������
int t;
double uk,rk,ik;
double Xk,Yk;
double doy = 152.0;
double R   = 287.054;
double g   = 9.80665;
double k1  = 77.604;
double k2  = 16.6;
double k3  = 377600.0;
double E;           //���Ǹ߶Ƚ�

//spp
double delta_tr[8];//�����
double delta[8];//�����Ӳ�
/************ȫ����Ҫ��������*************/

/**************************�ṹ�嶨��**************************/
//�㲥����
typedef struct{
    char a[3];
    int year, month, day, hour, mintue, second;	//ʱ����� minute
    double a_f0, a_f1, a_f2;					//�Ӳ���١���Ư

    double IDOE;								//����뾶�����ҵ��͸���������
    double C_rs;							//�ο�ʱ������ƽ�����ٶ������ֵ֮��
    double Delta_n0;									//�ο�ʱ�̵�ƽ�����
    double M0;                                  //?????

    double C_uc;								//γ�ȷ��ǵ����ҵ��͸���������
    double e;									//ƫ����
    double C_us;								//γ�ȷ��ǵ����ҵ��͸���������
    double Delta_A;								//�������ƽ����

    double Toe;									//�����ο�ʱ��
    double C_ic;								//�����ǵ����ҵ��͸���������
    double OMEGA0;								//����Ԫ��ʱ�̼���������㾭��
    double C_is;								//�����ǵ����ҵ��͸���������

    double i0;									//�ο�ʱ�̵Ĺ�����
    double C_rc;								//����뾶�����ҵ��͸���������
    double omega;								//���ص����
    double OMEGA_DOT;							//������ྭ�仯��

    double Idot;								//�����Ǳ仯��
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
aaa EPH[20];//���鳤��Ϊnum�ڴ�᲻��

//���ǽ�������
typedef struct{
    double Xk;
    double Yk;
    double Zk;
}bbb;
bbb location[8];//�����ظ�   m��

//����α�� �����Ӳ� �������ʱ ��������ʱ
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
//���
typedef struct{
    double X;
    double Y;
    double Z;
    double delta_0;
}qqq;
qqq result = {0,0,0,0};


/**************************�ṹ�嶨��**************************/



int read_data();//���ļ�

//�㲥����������������
double semiMajor_axis(int i);
double angular_Velocity(int i);								//1��ƽ�����ٶ�
double Naturalization_Time(int i);							//2���黯ʱ��
double meanAnomaly(int i);									//3��ƽ�����
double eccentricAnomanly(int i);							//4��ƫ�����Ek
double tureAnomaly(int i);									//5��������
double argumentOfOperigee(int i);						 	//6���������faik
double perturbation(int i);									//7(8)�������㶯������
double coordinates(int i);									//9��������ƽ������ϵ����
double MEOsate(int i);										//10�����λ��
double GEOsate(int i);

//������ӳ�
double bdgim3(int i);


//��������ʱ
double UNB3(int i);
double Matrix();
double yearAverage(double *p);
double amplitude(double *p);
double NMFdry(double *p);
double NMFZWD(double *p);
double calculate(double *p1,double *p2,double d_h,double d_nh);                  //������



//spp
double satellite_clock(int num,double delta[num]);                                    //���������Ӳ�
void matrix_mul(int l,int n,double a[l][n],double b[n][l]);                         //����ת��
void matrix_contrary(int l,double a[l][l],double c[l][l]);                          //��������
void matrix_ride(int l,int n,int q,double a[l][n],double b[n][q],double c[l][q]);   //�������
double derivative(int i,double X,double Y,double Z);                                //��ƫ��
double count(int num,double median[4][1]);                                          //����
void estimate(int num,double median[4][1]);                                         //�ж�����



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
