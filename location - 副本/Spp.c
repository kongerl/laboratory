#include <stdio.h>
#include "math.h"
#include "all.h"

double Spp()
{
    //���ļ�
    int num = 8;

    //�����Ӳ�

    satellite_clock(num,delta);
    printf("�����Ӳ�������\n");

    //�ж�����
    double median[4][1];
    estimate(num,median);

    //������
    printf("X      =  %lf\n",result.X);
    printf("Y      =  %lf\n",result.Y);
    printf("Z      =  %lf\n",result.Z);
    printf("delta  =  %lf\n",result.delta_0);

}


//���������Ӳ�
double satellite_clock(int num,double delta[num])
{


    double Tk[num];

    double F = -4.442807633e-10;

    for(int i = 0; i < num; i++)
    {
        delta_tr[i] = F * EPH[i].e * EPH[i].Delta_A * sin(Ek[i]);
        //    printf("����۸���ֵdelta_tr[%d]��%.20lf\n",i,delta_tr[i]);
    }

    //���ǲ������λʱ��ƫ��
    for(int i = 0; i < num; i++)
    {
        delta[i] = EPH[i].a_f0 + EPH[i].a_f1 * (-(da[i].rou)/C - 14.0) + EPH[i].a_f2 * pow(-(da[i].rou)/C - 14.0,2) + delta_tr[i] /*- da[i].rou/C*/;

       // printf("�����Ӳ�Ϊdelta[%d]��%.20lf\n",i,delta[i]);
    }

    return delta[num];
}

//����ת��
void matrix_mul(int l,int n,double a[l][n],double b[n][l])
{
    //  double a[l][n];
    //  double b[n][l];
//    printf("ԭ������Ϊ��\n");
//    for(int i = 0;i < l;i++)
//    {
//        for(int j = 0;j < n;j++)
//        {
//            printf("%lf   ",a[i][j]);
//        }
//        printf("\n");
//    }


    for(int i = 0;i < l;i++)
    {
        for(int j = 0;j < n;j++)
        {
            b[j][i] = a[i][j];
        }
    }
//    printf("ת�ú����Ϊ��\n");
//    for(int i = 0;i < n;i++)
//    {
//        for(int j = 0;j < l;j++)
//        {
//            printf("%lf   ",b[i][j]);
//        }
//        printf("\n");
//    }

}

//�������
void matrix_ride(int l,int n,int q,double a[l][n],double b[n][q],double c[l][q])
{
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < q; j++)
        {
            c[i][j] = 0;
            for(int k = 0; k < n; k++)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

//��������!!!
void matrix_contrary(int l,double a[l][l],double c[l][l])
{
    double temp[l];
    double b[l][l + l];
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            for (int k = l; k < 2 * l; k++)
            {
                b[i][j] = a[i][j];
                //�����������
                if ( (k - i) == l )
                {
                    b[i][k] = 1;
                }
                else
                {
                    b[i][k] = 0;
                }
            }
        }
    }

//    printf("��ֵ��\n");
//    for (int i = 0; i < l; i++)
//    {
//        for (int j = 0; j < 2 * l; j++)
//        {
//            printf("%lf  ", b[i][j]);
//        }
//        printf("\n");
//    }

    //��ߵ�λ��(��Ϊб�����ٵ�λ��)
    for (int i = 0; i < l; i++)
    {
        //ÿ�н����ͰѶԽ�����1
        int t = i;
        temp[i] = b[i][i];//��̬��1
        for(int a = 0; a < 2*l; a++) {
            b[t][a] = b[t][a] / temp[i];
        }
        //
        for (int j = 0; j < l; j++)
        {
            if (i == j) continue;//���Խ��߾�����

            double t1 = b[j][i];//�̶��˵���
            for (int k = 0; k < 2 * l; k++)
            {
                b[j][k] = b[j][k] - t1 * b[i][k];
            }
        }
    }

//    printf("��ߵ�λ����\n");
//    for(int i = 0;i < l;i++)
//    {
//        for(int j = 0;j < 2*l;j++)
//        {
//            printf("%lf  ",b[i][j]);
//        }
//        printf("\n");
//    }

    //��
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            c[i][j] = b[i][j + l];
        }
    }
//    for(int i = 0;i < 3;i++)
//    {
//        for(int j = 0;j < 3;j++)
//        {
//            printf("%lf     ",c[i][j]);
//        }
//        printf("\n");
//    }
}

//ƫ��
double derivative(int i,double X,double Y,double Z)
{
    double r;
    //����α��
    r= sqrt( pow(location[i].Xk - result.X, 2) + pow(location[i].Yk - result.Y, 2) + pow(location[i].Zk - result.Z, 2) );

    if(X != 0 && Y == 0 && Z == 0)
    {
        r = -(location[i].Xk - result.X) / r;
    }
    else if(X == 0 && Y != 0 && Z == 0)
    {
        r = -(location[i].Yk - result.Y) / r;
    }
    else if(X == 0 && Y == 0 && Z != 0)
    {
        r = -(location[i].Zk - result.Z) / r;
    }
    else
    {
        printf("������󣡣���\n");
    }
    return r;
}

//3���������Է�����
double count(int num,double median[4][1])
{

    double b[num][1];
    double G[num][4];
    double G1[4][num];
    double c[4][4];
    double c1[4][4];//��
    double c2[4][num];
    //��ֵG����
    for(int i = 0;i < num;i++)
    {
        G[i][0] = derivative(i,location[i].Xk,0,0);         //���ܴ���
        G[i][1] = derivative(i,0,location[i].Yk,0);
        G[i][2] = derivative(i,0,0,location[i].Zk);
        G[i][3] = 1.0;
    }
    //��ֵb����
    double r[num];
    for(int i = 0; i < num;i++)     //�޸ģ�����������������������
    {
        r[i] = sqrt( pow(location[i].Xk - result.X, 2) + pow(location[i].Yk - result.Y, 2) + pow(location[i].Zk - result.Z, 2) );
        
        b[i][0] = da[i].rou + delta[i] * C - UNB3(i) - bdgim3(i)- r[i] - result.delta_0;

    }
    //��ֵ
/*     printf("G����ֵΪ��\n");
   for(int i = 0;i < num;i++)
    {
        for(int j = 0;j < 4;j++)
        {
            printf("%lf     ",G[i][j]);
        }
        printf("\n");
    }

    printf("b����ֵΪ��\n");
    for(int i = 0;i < num;i++){
        printf("%lf     ",b[i][0]);
        printf("\n");
    }
    printf("\n");*/

    //ת��
    matrix_mul(num,4,G,G1);
/*     printf("G�����ת�þ���ֵΪ��\n");
   for(int i = 0;i < 4;i++)
    {
        for(int j = 0;j < num;j++)
        {
            printf("%lf   ",G1[i][j]);
        }
        printf("\n");
    }*/

    //���
    matrix_ride(4,num,4,G1,G,c);
/*    printf("��ˣ�\n");
    for(int i = 0;i < 4;i++)
    {
        for(int j = 0;j < 4;j++)
        {
            printf("%lf  ",c[i][j]);
        }
        printf("\n");
    }*/

    //����
    matrix_contrary(4,c,c1);
/*   printf("��˺����棺\n");
    for(int i = 0;i < 4;i++)
    {
        for(int j = 0;j < 4;j++)
        {
            printf("%lf     ",c1[i][j]);
        }
        printf("\n");
    }*/


    //�ٳ�
    matrix_ride(4,4,num,c1,G1,c2);
/* printf("�ٳˣ�\n");
    for(int i = 0;i < 4;i++)
    {
        for(int j = 0;j < num;j++)
        {
            printf("%lf  ",c2[i][j]);
        }
        printf("\n");
    }*/

    //���
    matrix_ride(4,num,1,c2,b,median);
    printf("����ֵ��\n");
    for(int i = 0;i < 4;i++)
    {
        printf("%lf  ",median[i][0]);
        printf("\n");
    }
}


//�ж��Ƿ�������������һ��
void estimate(int num, double median[4][1])
{
    double te;
    do
    {
        //����
        count(8,median);

        //�ж�����
        te = sqrt( pow(median[0][0], 2) + pow(median[1][0], 2) + pow(median[2][0], 2) );
        printf("te = %.20lf\n",te);

        //����
        result.X = result.X + median[0][0];
        result.Y = result.Y + median[1][0];
        result.Z = result.Z + median[2][0];
        result.delta_0 = result.delta_0 + median[3][0];

//        median[0][0] = result.X;
//        median[1][0] = result.Y;
//        median[2][0] = result.Z;
//        median[3][0] = result.delta_0;

    }while(te > 1e-8);

}
