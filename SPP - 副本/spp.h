//
// Created by h'sh on 2023/11/6
//.

#ifndef SPP_SPP_H
#define SPP_SPP_H


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

typedef struct{//��������
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
};//�����ظ�   km����������������

typedef struct{//���
    double X;
    double Y;
    double Z;
    double delta_0;
}ccc;
ccc result = {0,0,0,0};

int read_data();//���ļ�
void satellite_clock(int num,double delta[num]);//���������Ӳ�
void matrix_mul(int l,int n,double a[l][n],double b[n][l]);//����ת��
void matrix_contrary(int l,double a[l][l],double c[l][l]);//��������
void matrix_ride(int l,int n,int q,double a[l][n],double b[n][q],double c[l][q]);//�������
double derivative(int i,double X,double Y,double Z);//��ƫ��
double count(int num,double median[4][1]);//����

void updata(double median[4][1]);//���·����Է�����ĸ�
void estimate(int num,double median[4][1]);//�ж�����

#endif //SPP_SPP_H
