#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
ads coord = { //-14003426.312 , 22219102.603 , -2658753.980,//GEO1
              //12773.329867  , 7131.757565 ,-22113.003365,//IGSO6
              //19790.745902  ,-7229.863015 ,-16141.244320,//MEO11

             4550111.541 , 41880407.096 , -1237547.731,//GEO2
              // -11196311.788 , 40448493.438 , -1598137.482,//IGSO6
              //15146115.254  ,-6778785.178 ,-22373553.567,//MEO11


            //-24114776.912,25117497.232,-23852197.899,
             // -847293.626,  27828168.883 ,  1572453.676,

             //23666.110582,-966.774962,14732.695897,


              114.3572611   ,30    ,28.2
};


/*********************函数声明*********************/
double Matrix();
double yearAverage(double *p);
double amplitude(double *p);

double NMFdry(double *p);
double NMFZWD(double *p);

double calculate(double *p1,double *p2,double d_h,double d_nh);                  //最后计算
/************************************************/

double doy = 152.0;
double pi  = 3.14159265358979323846;
double R   = 287.054;
double g   = 9.80665;
double k1  = 77.604;
double k2  = 16.6;
double k3  = 377600.0;

double E;           //卫星高度角

int main()
{



    double AVG[5]   = {0.0,sizeof(0.0)};
    double AMP[5]   = {0.0,sizeof(0.0)};
    double X_doy[5] = {0.0,sizeof(0.0)};

    Matrix();
    printf("E = %lf\n",E);

//    coord.latitude = coord.latitude * (pi/180);
//    coord.longitude = coord.longitude * (pi/180);
//    coord.high = coord.high * (pi/180);

    //年均变化值
    yearAverage(AVG);
    printf("年均值AVG：\n");
    for(int i = 0; i < 5; i++)
    {
        printf("%lf   ",AVG[i]);
    }

    //振幅值
    amplitude(AMP);
    printf("振幅AMP：\n");
    for(int i = 0; i < 5; i++)
    {
        printf("%lf   ",AMP[i]);
    }

    //年积日的值
    printf("\n年积日X_doy：\n");
    for(int i = 0; i < 5; i++)
    {
        X_doy[i] = AVG[i] - AMP[i] * cos( ((doy - 28) * 2 * pi) / 365.25 );
        printf("%lf   ",X_doy[i]);
    }

    //气象参考值
    double g_m;
    g_m = 9.784*(1 - 2.66e-3 * cos(2 * coord.latitude) - 2.8e-7 * coord.high);
    printf("\n圆柱体大气重力加速度g_m = %lf",g_m);

    double lamda;
    lamda = X_doy[4] + 1.0;
    printf("\nlamda = %lf",lamda);

    double Tm;
    Tm = X_doy[1] * ( 1 -( (X_doy[3] * R) / (g_m * lamda) ) );
    printf("\n水蒸气平均温度Tm = %lf",Tm);

    double d_h;//干
    d_h = pow((1 - (X_doy[3] * coord.high / X_doy[1]) ), ( g / R / X_doy[3]) );
    d_h = ( ( 1e-6 * k1 * R ) / g_m ) * X_doy[0] * d_h ;
    printf("\n干分量结果d_h = %lf",d_h);

    double d_nh;//湿
    d_nh = pow(1 - X_doy[3] / X_doy[1] , (lamda * g) / (R * X_doy[3]) - 1);
    d_nh = 1e-6 * (Tm * k2 + k3) * R * X_doy[2] * d_nh;
    d_nh = d_nh / ((g_m * lamda - X_doy[3] * R) * X_doy[1]);
    printf("\n湿分量d_nh = %lf\n",d_nh);

    printf("\n干湿分量之和 = %lf\n\n",d_h + d_nh);

    //投影函数NMF干分量
    double P_h[3];
    NMFdry(P_h);
    for(int i = 0; i < 3; i++)
    {
        printf("投影函数NMF干分量P_h[%d] = %lf    \n",i,P_h[i]);
    }

    //投影函数NMF湿分量
    double P_ht[3];
    NMFZWD(P_ht);
    for(int i = 0; i < 3; i++)
    {
        printf("投影函数NMF湿分量P_ht[%d] = %lf   \n",i,P_ht[i]);
    }

    //最后计算
    calculate(P_h,P_ht,d_h,d_nh);


    return 0;
}
//坐标转换
double Matrix()
{
    double a1 = coord.latitude * (pi/180);
    double b1 = coord.longitude * (pi/180);

    double a[3][3] = {-sin(b1)            , cos(b1)           , 0         ,
                      -sin(a1)*cos(b1) , -sin(a1)*sin(b1), cos(a1),
                       cos(a1)*cos(b1) , cos(a1)*sin(b1) , sin(a1)
    };
    double b[3];
    double c[3]    = {0,0,0};

    //结构体赋值给数组
    ads *p = &coord;
    b[0]    =   (p->Xk) - (-2267750.266);
    b[1]    =   (p->Yk) - 5009156.142;
    b[2]    =   (p->Zk) - 3221291.898;

    //矩阵计算
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            c[i] += a[i][j] * b[j];
        }
    }

    //矩阵结果赋值给结构体
    p->East     =    c[0];
    p->North    =    c[1];
    p->Up       =    c[2];

    printf("East = %lf\n",coord.East);
    printf("North = %lf\n",coord.North);
    printf("Up = %lf\n",coord.Up);



    double t;
    t =  (pow(coord.North,2) + pow(coord.East,2) + pow(coord.Up,2));
    E = asin(coord.Up/sqrt(t));

}

//年均值
double yearAverage(double *p)
{
    int count = -1;
    double t = coord.latitude;
    double LAT = coord.latitude;
    double table[5][5] = {1013.25, 299.65, 26.31, 6.3e-3 , 2.77,
                          1017.25, 294.15, 21.79, 6.05e-3, 3.15,
                          1015.75, 283.15, 11.66, 5.58e-3, 2.57,
                          1011.75, 272.15, 6.78 , 5.39e-3, 1.81,
                          1013.00, 263.65, 4.11 , 4.53e-3, 1.55
    };
    //打印表中数据
//        for(int i = 0; i < 5; i++)
//    {
//        for(int j = 0; j < 5; j++)
//        {
//            printf("%lf   ",table[i][j]);
//        }
//        printf("\n");
//    }
    while(t > 15.0)
    {
        t = t - 15;
        count++;
    }
    //printf("count = %d\n",count);
    //判断计算
    if(LAT <= 15.0)
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[0][i];
        }
    }
    else if(LAT >= 75.0)
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[4][i];
        }
    }
    else
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[count][i] +( (table[count + 1][i] - table[count][i] ) / 15 ) * (LAT - (count + 1) * 15.0);
        }
    }

}
//振幅
double amplitude(double *p)
{
    int count = -1;
    double t = coord.latitude;
    double LAT = coord.latitude;
    double table[5][5]  = {  0.00 , 0.00 , 0.00, 0.00   , 0.00,
                             -3.75, 7.00 , 8.85, 0.25e-3, 0.33,
                             -2.25, 11.00, 7.24, 0.32e-3, 0.46,
                             -1.75, 15.00, 5.36, 0.81e-3, 0.74,
                             -0.50, 14.50, 3.39, 0.62e-3, 0.30,

    };

    while(t > 15.0)
    {
        t = t - 15;
        count++;
    }
    //printf("count = %d\n",count);
    //判断计算
    if(LAT <= 15.0)
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[0][i];
        }
    }
    else if(LAT >= 75.0)
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[4][i];
        }
    }
    else
    {
        for(int i = 0; i < 5; i++)
        {
            p[i] = table[count][i] +( (table[count + 1][i] - table[count][i] ) / 15 ) * (LAT - (count + 1) * 15.0);
        }
    }
}
//参考值

//NMF干分量
double NMFdry(double *p)
{
    int count = -1;
    double t = coord.latitude;
    double LAT = coord.latitude;
    double table[5][3]  = {  1.2769934e-3, 2.9153695e-3, 6.2610505e-2,
                             1.2683230e-3, 2.9152299e-3, 6.2837393e-2,
                             1.2465397e-3, 2.9288445e-3, 6.3721774e-2,
                             1.2196049e-3, 2.9022565e-3, 6.3824265e-2,
                             1.2045996e-3, 2.9024912e-3, 6.4258455e-2
    };
    double table1[5][3] = {  0.0,          0.0,          0.0,
                             1.2709626e-5, 2.1414979e-5, 9.0128400e-5,
                             2.6523662e-5, 3.0160779e-5, 4.3497037e-5,
                             3.4000452e-5, 7.2562722e-5, 8.4795348e-4,
                             4.1202191e-5, 1.1723375e-4, 1.7037206e-3
    };
    while(t > 15.0)
    {
        t = t - 15;
        count++;
    }
    //printf("count = %d\n",count);
    //判断计算
    if(LAT < 15.0)
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = table[0][i] + table[0][i] * cos( (2*pi*(t - 28)) / 365.25);
        }
    }
    else if(LAT > 75.0)
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = table[4][i] + table[4][i] * cos( (2*pi*(t - 28)) / 365.25);
        }
    }
    else
    {
        for(int i = 0; i < 3; i++)
        {
//            printf("avg(i) = %lf    \n",table[count][i]);
//            printf("amp(i) = %lf    \n",table1[count][i]);
            p[i] = table[count][i] + ( table[count + 1][i] - table[count][i] ) * ( (LAT - (count + 1) * 15.0) / 15.0 );
            p[i] = p[i] +  ( table1[count][i] + (table1[count + 1][i] - table1[count][i]) * ( (LAT - (count + 1) * 15.0) / 15.0 ) * cos( 2 * pi * (t-28)/365.25 ) );
        }
    }

}
//NMF湿分量
double NMFZWD(double *p)
{
    int count = -1;
    double t = coord.latitude;
    double LAT = coord.latitude;
    double table[5][3]  = {  5.8021897e-4, 1.4275268e-3, 4.3472961e-2,
                             5.6794847e-4, 1.5138625e-3, 4.6729510e-2,
                             5.8118019e-4, 1.4572752e-3, 4.3908931e-2,
                             5.9727542e-4, 1.5007428e-3, 4.4626982e-2,
                             6.1641693e-4, 1.7599082e-3, 5.4736038e-2
    };

    while(t >= 15.0)
    {
        t = t - 15;
        count++;
    }
    //printf("count = %d\n",count);
    //判断计算
    if(LAT < 15.0)
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = table[0][i];
        }
    }
    else if(LAT > 75.0)
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = table[4][i];
        }
    }
    else
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = table[count][i] +( (table[count + 1][i] - table[count][i] ) / 15 ) * (LAT - (count + 1) * 15.0);
        }
    }

}
//最后计算
double calculate(double *p1,double *p2,double d_h,double d_nh)
{
    double a_ht = 2.53e-5;
    double b_ht = 5.49e-3;
    double c_ht = 1.14e-3;

    double M_he;
    double M_he1;
    double M_w;

    double result;

    M_he = 1 + p1[0] / (1 + p1[1] / (1 + p1[2]) );
    M_he = M_he / ( sin(E) + p1[0] / (sin(E) + p1[1] / (sin(E) + p1[2]) ) );
    //printf("M_he = %lf\n",M_he);

    M_he1 = 1 + a_ht / (1 + b_ht / (1 + c_ht) );
    M_he1 = M_he1 / ( sin(E) + a_ht / (sin(E) + b_ht / (sin(E) + b_ht) ) );
    M_he1 = ((1 / sin(E) - M_he1) * coord.high) / 1000;
    //printf("M_he1 = %lf\n",M_he1);

    M_he = M_he + M_he1;
    printf("干分量投影M_he = %lf\n",M_he);

    M_w = 1 + p2[0] / (1 + p2[1] / (1 + p2[2]) );
    M_w = M_w / ( sin(E) + p2[0] / (sin(E) + p2[1] / (sin(E) + p2[2]) ) );
    printf("湿分量投影M_w = %lf\n",M_w);

    result = d_h * M_he + d_nh * M_w;
    printf("result = %lf",result);

}