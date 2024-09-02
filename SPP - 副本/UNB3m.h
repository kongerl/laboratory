#ifndef SPP_UNB3M_H
#define SPP_UNB3M_H


typedef struct{
    double Xk;
    double Yk;
    double Zk;
    double longitude;
    double latitude;
    double high;
    double East;
    double North;
    double Up;
}ads;
//ads coord = { -2267749  , 5009154   , 3221290,
//              114.3573   ,30.5317    ,25.8
//};
ads coord = {4550111.541 , 41880407.096 , -1237547.731,//GEO2
             114.3572611   ,30    ,28.2
};



double doy = 152.0;
double pi  = 3.14159265358979323846;
double R   = 287.054;
double g   = 9.80665;
double k1  = 77.604;
double k2  = 16.6;
double k3  = 377600.0;

double E;           //Œ¿–«∏ﬂ∂»Ω«


double UNB(double x,double y,double z);
void Matrix();
void yearAverage(double *p);
void amplitude(double *p);

void NMFdry(double *p);
void NMFZWD(double *p);

double calculate(double *p1,double *p2,double d_h,double d_nh);

#endif //SPP_UNB3M_H
