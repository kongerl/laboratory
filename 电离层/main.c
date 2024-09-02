#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <stddef.h>       //结构体指针赋值头文件

double pi       =   3.1425926535;
double H_ion    =   4.0e5;
double Re       =   6.378e6;


typedef struct{
    double Y;
    double M;
    double d;
    double h;
}Time;
Time time = {2023,6,1,0};

typedef struct {
    double Xk;
    double Yk;
    double Zk;
    double Longitude;           //经度
    double Latitude;            //纬度
    double High;                //高度
    double East;
    double North;
    double Up;
}Coordinate;
Coordinate User = {-2267750.266, 5009156.142, 3221291.898,
                   108.904893,34.152933 ,410,
                   0.0, 0.0, 0.0};

typedef struct{
    //角度
    double E;
    double A;
    double Fai;

    double z;                   //天顶角
    double el;                  //卫星仰角
    //地理
    double long_g;
    double fai_g;
    //地固坐标系地磁
    double long_m;
    double fai_m;
    //日固坐标系地磁
    double long_m1;
    double fai_m1;

}para;
para data;

/***********函数定义*************/
double Matrix();
double location();
double Ai_calculate(int i);
double factorial(double a,int b);                              //阶乘函数
double legendre(double n,double m);             //勒让德函数
double forecast();              //电离层延迟预报值A0
double VTEC();                                                  //穿刺点处电离层延迟
double M_function();                                         //穿刺点处投影函数Mf

/******************************/
int main() {
    double a[9] = {};

    printf("用户输入数据：%lf %lf %lf %lf %lf %lf \n",
           User.Xk,User.Yk ,User.Zk,User.Latitude,User.Longitude,User.High );
    //坐标系转换
    Matrix();
//    printf("User.East   =   %lf\n",User.East);
//    printf("User.North  =   %lf\n",User.North);
//    printf("User.Up     =   %lf\n",User.Up);
    //电离层穿刺点位置
    location();
    //Ai的计算
    Ai_calculate(5);
    //A0的计算
    forecast();

    //VTEC计算
    double V = VTEC();
    printf("VTEC = %lf",V);



    //Mf计算
//    double Mf1;
//    Mf1 = M_function();
//    printf("Mf = %lf\n",Mf1);
    //最终计算

    return 0;
}
//坐标矩阵转换，卫星X，Y，Z转变为站心E，N，U
double Matrix()
{
    double a[3][3] = {-sin(User.Longitude)                     ,   cos(User.Longitude)                     , 0                 ,
                      -sin(User.Latitude)*cos(User.Longitude), -sin(User.Latitude)*sin(User.Longitude), cos(User.Latitude),
                      cos(User.Latitude)*cos(User.Longitude),  cos(User.Latitude)*cos(User.Longitude), sin(User.Latitude)
    };
    double b[3];
    double c[3]    = {0,0,0};

    //结构体赋值给数组
    Coordinate *p = &User;
    b[0]    =   p->Longitude;
    b[1]    =   p->Latitude;
    b[2]    =   p->High;
//    for(int i = 0; i < 3; i++)
//    {
//        //b[i] = p->Longitude;
//        b[i] = *(double*)((char*)p + offsetof(Coordinate, Longitude) + i * sizeof(double));
//        printf("b[%d] = %lf\n",i,b[i]);
//    }

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

}
//电离层穿刺点位置
double location()
{
    double long_M      =   (-72.58 / 180) * pi;
    double fai_M    =   ( 80.27 / 180) * pi;

    //地心张角
    double a;
    a = (pow(User.North,2) + pow(User.East,2) ) / (pow(User.North,2) + pow(User.East,2) + pow(User.Up,2));
    data.E = acos(sqrt(a));
    data.A = atan(User.East / User.North );

    data.Fai = pi/2 - data.E - asin(Re/(Re + H_ion) * cos(data.E));
    printf("a = %lf data.E = %lf data.A = %lf data.Fai = %lf\n",a,data.E,data.A,data.Fai);

    //地理坐标系
    data.fai_g = asin( sin(User.Latitude)*cos(data.Fai) + cos(User.Latitude)*sin(data.Fai)*cos(data.A) );
    data.long_g = User.Longitude + atan(sin(data.Fai)*sin(data.A)*cos(User.Latitude) /
                                        (cos(data.Fai) - sin(User.Latitude)*sin(data.fai_g)));
    printf("地理纬度fai_g = %lf\n",data.fai_g);
    printf("地理精度long_g = %lf\n",data.long_g);
    //地固坐标系
    data.fai_m = asin(sin(fai_M)*sin(data.fai_g) + cos(fai_M)*cos(data.fai_g)*cos(data.long_g-long_M));
    data.long_m = atan(cos(data.fai_g)*sin(data.long_g - fai_M)*cos(fai_M) /
                       (sin(fai_M)*sin(data.fai_m) - sin(data.fai_g)) );
    printf("地固地磁纬度fai_m = %lf\n",data.fai_m);
    printf("地固地磁精度long_m = %lf\n",data.long_m);
    //日固坐标系
    double t;//约化儒略日
    double S_ion;

    t = 1721013.5 + 365*time.Y - floor(7/4*(time.Y + floor( (time.M + 9)/12 )) ) + time.d + time.h /24 + floor((time.M*275)/9);
    t = t - 2400000.5;

    S_ion = pi*(1-2*(t-floor(t)));
    printf("t = %lf\n",t);
    printf("S_ion = %lf\n",S_ion);
    data.fai_m1 = data.fai_m;
    data.long_m1 = data.long_m - atan(sin(S_ion - long_M) /
                                      sin(fai_M)*cos(S_ion - fai_M));
    printf("日固地磁纬度fai_m1 = %lf\n",data.fai_m1);
    printf("日固地磁精度long_m1 = %lf\n",data.long_m1);


}
//Ai的计算
double Ai_calculate(int i)
{
    int n,m;
    double delta;
    double N_nm,P_nm,A_i;

    switch(i)
    {
        case 1:
            n = 0; m = 0;
            break;
        case 2:
            n = 1; m = 0;
            break;
        case 3:
            n = 1; m = 1;
            break;
        case 4:
            n = 1; m = -1;
            break;
        case 5:
            n = 2; m = 0;
            break;
        case 6:
            n = 2; m = 1;
            break;
        case 7:
            n = 2; m = -1;
            break;
        case 8:
            n = 2; m = 2;
            break;
        case 9:
            n = 2; m = -2;
            break;
    }
    printf("n = %d,m = %d\n",n,m);

    int n1 = abs(n);
    int m1 = abs(m);

    if(m1 = 0)
        delta = 1;
    else
        delta = 0;

    N_nm = sqrt( factorial(n1-m1,1)*(2*n1+1)*(2-delta) / factorial(n1+m1,1) );
    printf("N_nm = %lf\n",N_nm);

    P_nm = legendre(n1,m1);
    printf("P_nm = %lf\n",P_nm);

    if(m >= 0)
    {
        A_i = N_nm * P_nm * sin(data.fai_m1) * cos(m*data.long_m1);
    }
    else
    {
        A_i = N_nm * P_nm * sin(data.fai_m1) * sin(-m*data.long_m1);
    }
    return A_i;
}
//阶乘函数(尾数决定几阶乘)
double factorial(double a,int b)
{
    int n = 1;
    if(b == 1)
    {
        if(a == 0) return 1;
        for(int i = 1;i <= a;i++)
        {
            n = n * i;
        }
    }
    else if(b == 2)
    {
        if(a == 0) return 1;
        for(int i = 1;i <= a;i++)
        {
            n *= 2 * i - 1;
        }
    }
    return n;

}
//勒让德函数                     三四有问题感觉！！！！！
double legendre(double n1,double m1)
{
    double n;
    if(n1 == m1 && n1 == 0 )
    {
        n = 1.0;
    }
    else if(n1 == m1)
    {
        n = factorial(2*n1 - 1,2) * pow(   (1-pow(sin(data.fai_m1),2))       ,n1/2.0);
    }
    else if(n1 == m1 + 1.0)
    {
        n = sin(data.fai_m1) * (2*m1 + 1) * legendre(m1,m1);
    }
    else
    {   //n必须要比m大!!!
        n = (2*n1-1) * sin(data.fai_m1) * legendre(n1-1,m1) - (n1+m1-1) * legendre(n1-2,m1);
        n = n / (n1-m1);
    }
    return n;
}
//电离层延迟预报值A0
double forecast()
{
    double B_j;

 return 1;
}

//穿刺点处电离层延迟
double VTEC()
{
    double a[9] = {1,1,};
    double V;
    double A0 = forecast();

    for(int i = 1; i <= 9; i++)
    {
        V = a[i] * Ai_calculate(i);
    }
    V = V + A0;
    return V;
}
//穿刺点处投影函数Mf
double M_function()
{
    double Mf;
    Mf = 1 - pow(Re*cos(data.E) / (Re+ H_ion),2);
    Mf = 1 / sqrt(Mf);
    return Mf;
}
//