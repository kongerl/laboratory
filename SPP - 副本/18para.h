//
// Created by h'sh on 2023/11/13.
//

#ifndef SPP_18PARA_H
#define SPP_18PARA_H


typedef struct{
    double user_Xk;
    double user_Yk;
    double user_Zk;
    double CGCS_Xk;
    double CGCS_Yk;
    double CGCS_Zk;
}MEO;
MEO meo;

typedef struct{
    double user_Xk;
    double user_Yk;
    double user_Zk;
    double CGCS_Xk;
    double CGCS_Yk;
    double CGCS_Zk;
}GEO;
GEO geo;
//参数定义(不用修改)



//全局变量需求



double semiMajor_axis(int a);
double angular_Velocity(int a);								//1、平均角速度
double Naturalization_Time(int a);							//2、归化时间
double meanAnomaly(int a);									//3、平近点角
double eccentricAnomanly(int a);								//4、偏近点角
double tureAnomaly(int a);									//5、真近点角
double argumentOfOperigee(int a);						 	//6、升交距角faik
double perturbation(int a);									//7(8)、计算摄动改正量
double coordinates(int a);									//9、计算轨道平面坐标系坐标
double MEOsate(int a);										//10、轨道位置
double GEOsate(int a);

#endif //SPP_18PARA_H
