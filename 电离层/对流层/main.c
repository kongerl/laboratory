#include <stdio.h>
#include <math.h>
#include <stddef.h>

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
Coordinate Xupt = {-14003.426312, 22219.102603, -2658.753980,
                   108.904893,34.152933 ,410,
                   0.0, 0.0, 0.0};

typedef struct{
    //角度
    double z;                   //天顶角
    double el;                  //卫星仰角

    double e;


}para;
para data;

/***********函数定义*************/
double Matrix();


/******************************/
int main() {


    Xupt.Longitude  =   34.152933;
    Xupt.Latitude   =   108.904893;
    Xupt.High       =   410;

    printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
           Xupt.Xk,Xupt.Yk ,Xupt.Zk,Xupt.Latitude,Xupt.Longitude,Xupt.High,Xupt.East ,Xupt.North,Xupt.Up );
    //坐标系转换
    Matrix();
    printf("Xupt.East   =   %lf\n",Xupt.East);
    printf("Xupt.North  =   %lf\n",Xupt.North);
    printf("Xupt.Up     =   %lf\n",Xupt.Up);

    //角度计算





    return 0;
}
//坐标矩阵转换
double Matrix()
{
    double a[3][3] = {-sin(Xupt.Longitude)                     ,   cos(Xupt.Longitude)                     , 0                 ,
                      -sin(Xupt.Latitude)*cos(Xupt.Longitude), -sin(Xupt.Latitude)*sin(Xupt.Longitude), cos(Xupt.Latitude),
                      cos(Xupt.Latitude)*cos(Xupt.Longitude),  cos(Xupt.Latitude)*cos(Xupt.Longitude), sin(Xupt.Latitude)
    };
    double b[3];
    double c[3]    = {0,0,0};

    //结构体赋值给数组
    Coordinate *p = &Xupt;
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
//







