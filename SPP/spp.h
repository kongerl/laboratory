//
// Created by h'sh on 2023/11/6
//.

#ifndef SPP_SPP_H
#define SPP_SPP_H

typedef struct{//广播星历
    char a[4];
    int year, month, day, hour, mintue, second;	//时间参数 minute
    double a_f0, a_f1, a_f2;					//钟差、钟速、钟漂

    double IDOE;								//轨道半径的正弦调和改正项的振幅
    double C_rs;							//参考时刻卫星平均角速度与计算值之差
    double Delta_n0;									//参考时刻的平近点角
    double M0;

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

    double Ek;
}ggg;
ggg EPH[8];

typedef struct{
    double rou;
    double delta_tu;
    double I_n;
    double T_n;
}qqq;
qqq da[8] = {24688564.000,0,0,0,
             21780791.271,0,0,0,
             25360248.071,0,0,0,
             22833624.564,0,0,0,
             23687638.830,0,0,0,
             37857941.485,0,0,0,
             36562515.356,0,0,0,
             21773519.462,0,0,0,


};

typedef struct{//卫星坐标
    double Xk;
    double Yk;
    double Zk;
}aaa;
aaa location[8] = {  11198739.076 ,  12158338.761 ,  22468347.179,
                    -12242344.789 ,  17647109.141 ,  17805220.662,
                    -15867517.469 , -3718708.129  ,  22641656.201,
                    -847293.626   ,  27828168.883 ,  1572453.676,
                    -23558951.204 ,  14818393.200 ,  2096312.828,
                    -13557460.955 ,  38501905.749 , -10070933.818,
                    -6643035.003  ,  24340177.672 ,  33822615.073,
                    -9177345.895  ,  19504580.681 ,  17736830.826
};//命名重复   km级？？？？？？？

typedef struct{//结果
    double X;
    double Y;
    double Z;
    double delta_0;
}ccc;
ccc result[100] = {};
double Spp();
int read_data();//读文件
void satellite_clock(int num,double delta[num]);//计算卫星钟差
void matrix_mul(int l,int n,double a[l][n],double b[n][l]);//矩阵转置
void matrix_contrary(int l,double a[l][l],double c[l][l]);//矩阵求逆
void matrix_ride(int l,int n,int q,double a[l][n],double b[n][q],double c[l][q]);//矩阵相乘
double derivative(int i,double X,double Y,double Z);//求偏导
double count(int num,double median[4][1]);//计算

void updata(double median[4][1]);//更新非线性方程组的根
void estimate(int num,double median[4][1]);//判断收敛

#endif //SPP_SPP_H
