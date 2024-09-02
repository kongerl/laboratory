#include "stdio.h"
#include <string.h>
#include "math.h"

#include "spp.h"


int k = 0;

double Spp()
{
    //���ļ�
    read_data();
    int num = 8;

    //�����Ӳ�
    double delta[num];
    satellite_clock(num,delta);

    //�ж�����
    double median[4][1];
    estimate(num,median);

    //������
    printf("X      =  %lf\n",result[k].X);
    printf("Y      =  %lf\n",result[k].Y);
    printf("Z      =  %lf\n",result[k].Z);
    printf("delta  =  %lf\n",result[k].delta_0);

}

//���ļ�
int read_data()
{
    FILE *file = fopen("E:\\ClionProject\\BRD400DLR_S_20231520000_01D_MN.rnx", "r");  // ���ļ��������ļ���Ϊ data.txt��

    if (file == NULL)
    {
        printf("�޷����ļ�\n");
        return 1;
    }

    char line[100];  // ����ÿ��������100���ַ�
    while (fgets(line, sizeof(line), file)) {
        if (    strstr(line, "EPH C26 CNV1") || strstr(line, "EPH C29 CNV1") || strstr(line, "EPH C30 CNV1")
             || strstr(line, "EPH C35 CNV1") || strstr(line, "EPH C36 CNV1") || strstr(line, "EPH C39 CNV1")
             || strstr(line, "EPH C40 CNV1") || strstr(line, "EPH C45 CNV1") )
        {
            fgets(line, sizeof(line), file);        //��ȡ��һ��
            if(    strstr(line, "C26 2023 06 01 00 00 00") || strstr(line, "C29 2023 06 01 00 00 00") || strstr(line, "C30 2023 06 01 00 00 00")
                || strstr(line, "C35 2023 06 01 00 00 00") || strstr(line, "C36 2023 06 01 00 00 00") || strstr(line, "C39 2023 06 01 00 00 00")
                || strstr(line, "C40 2023 06 01 00 00 00") || strstr(line, "C45 2023 06 01 00 00 00") )
            {
                do {
                    int count = 0;
                    //1
                    sscanf(line, "%s %d %d %d %d %d %d %lf %lf %lf ", EPH[count].a, &EPH[count].year, &EPH[count].month, &EPH[count].day, &EPH[count].hour,&EPH[count].mintue, &EPH[count].second,&EPH[count].a_f0, &EPH[count].a_f1, &EPH[count].a_f2);
                    printf("%s %d %d %d %d %d %d\t%.20lf \t %.20lf \t %.20lf\n", EPH[count].a, EPH[count].year, EPH[count].month, EPH[count].day, EPH[count].hour, EPH[count].mintue,EPH[count].second,EPH[count].a_f0, EPH[count].a_f1, EPH[count].a_f2);
                    fgets(line, sizeof(line), file);
                    //2
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].IDOE, &EPH[count].C_rs, &EPH[count].Delta_n0, &EPH[count].M0);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].IDOE, EPH[count].C_rs, EPH[count].Delta_n0, EPH[count].M0);
                    fgets(line, sizeof(line), file);
                    //3
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].C_uc, &EPH[count].e, &EPH[count].C_us, &EPH[count].Delta_A);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].C_uc, EPH[count].e, EPH[count].C_us, EPH[count].Delta_A);
                    fgets(line, sizeof(line), file);
                    //4
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].Toe, &EPH[count].C_ic, &EPH[count].OMEGA0, &EPH[count].C_is);
                    printf("\t\t %.20lf \t %.20lf  %.20lf \t %.20lf\n", EPH[count].Toe, EPH[count].C_ic, EPH[count].OMEGA0, EPH[count].C_is);
                    fgets(line, sizeof(line), file);
                    //5
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].i0, &EPH[count].C_rc, &EPH[count].omega, &EPH[count].OMEGA_DOT);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].i0, EPH[count].C_rc, EPH[count].omega,EPH[count].OMEGA_DOT);
                    fgets(line, sizeof(line), file);
                    //6
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].Idot, &EPH[count].data_bit, &EPH[count].BDT, &EPH[count].A_DOT);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].Idot, EPH[count].data_bit, EPH[count].BDT,EPH[count].A_DOT);
                    fgets(line, sizeof(line), file);
                    //7
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].SV, &EPH[count].Health, &EPH[count].IGD, &EPH[count].ISC);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].SV, EPH[count].Health, EPH[count].IGD, EPH[count].ISC);
                    fgets(line, sizeof(line), file);
                    //8
                    sscanf(line, " %lf %lf %lf %lf", &EPH[count].Launch_time, &EPH[count].IDOC, &EPH[count].Delta_n0_dot, &EPH[count].type);
                    printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].Launch_time, EPH[count].IDOC,EPH[count].Delta_n0_dot, EPH[count].type);
                    fgets(line, sizeof(line), file);
                    count++;
                } while (    strstr(line, "EPH C26 CNV1") || strstr(line, "EPH C29 CNV1") || strstr(line, "EPH C30 CNV1")
                          || strstr(line, "EPH C35 CNV1") || strstr(line, "EPH C36 CNV1") || strstr(line, "EPH C39 CNV1")
                          || strstr(line, "EPH C40 CNV1") || strstr(line, "EPH C45 CNV1") );
            }
        }
    }

    fclose(file);  // �ر��ļ�

    return 0;
}

//���������Ӳ�
void satellite_clock(int num,double delta[num])
{
    int t = -14;
    double delta_tr[8];

    double GM = 3.986004418e14;
    double A_ref_MEO = 2.79061e7;
    double Tk[num];
    //����۸���ֵ
    double F = -4.442807633e-10;

    for(int i = 0; i < num; i++)
    {
         Tk[i] = t - EPH[i].Toe;
        Tk[i] -= 14;
        if( Tk[i] > 3.024e5 )
        {
            Tk[i] -= 6.048e5;
        }
        else if( Tk[i] < -3.024e5)
        {
            Tk[i] += 6.048e5;
        }

        double A0_MEO = A_ref_MEO + EPH[i].Delta_A;
        double Ak = A0_MEO + EPH[i].A_DOT * Tk[i];
        double angV = sqrt(GM / pow(A0_MEO,3) ) + EPH[i].Delta_n0 + 0.5 * EPH[i].Delta_n0_dot * Tk[i];
        double Mk = EPH[i].M0 + angV * Tk[i];
        double temp,Ek;
        temp = Mk;
        for(int i = 0; i<100; i++)
        {
            double Ek = Mk + EPH[i].e * sin(temp); //��Ϊʲô��Mk����������M0�𣿣���
            temp = Ek;
        }
        EPH[i].Ek = Ek;
        delta_tr[i] = F * EPH[i].e * pow(EPH[i].Delta_A,2) * sin(EPH[i].Ek);//����۸���
    }

//    double t_sv,delta_t_sv;//�źŷ���ʱ�����ǲ������λʱ��    ���ǲ������λʱ��ƫ��
//    //�����źŷ���ʱ��   ��ô��
//    EPH.t;

    //���ǲ������λʱ��ƫ��
    for(int i = 0; i < num; i++)
    {

        delta[i] = EPH[i].a_f0 + EPH[i].a_f1 * (t - 0) + EPH[i].a_f2 * pow(t - 0,2) + delta_tr[i]  ;
    }


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
    r= sqrt( pow(location[i].Xk - result[k].X, 2) + pow(location[i].Yk - result[k].Y, 2) + pow(location[i].Zk - result[k].Z, 2) );

    if(X != 0 && Y == 0 && Z == 0)
    {
        r = -(location[i].Xk - result[k].X) / r;
    }
    else if(X == 0 && Y != 0 && Z == 0)
    {
        r = -(location[i].Yk - result[k].Y) / r;
    }
    else if(X == 0 && Y == 0 && Z != 0)
    {
        r = -(location[i].Zk - result[k].Z) / r;
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
    for(int i = 0; i < num;i++)
    {
        r[i] = sqrt( pow(location[i].Xk - result[k].X, 2) + pow(location[i].Yk - result[k].Y, 2) + pow(location[i].Zk - result[k].Z, 2) );
        b[i][0] = da[i].rou + r[i] - result[k].delta_0;
    }
    //��ֵ
    printf("G����ֵΪ��\n");
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
    printf("\n");

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
/*    printf("��˺����棺\n");
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
    printf("���\n");
    for(int i = 0;i < 4;i++)
    {
        printf("%lf  ",median[i][0]);
        printf("\n");
    }
}

//���·����Է�����ĸ�
void updata(double median[4][1])
{
    for(int i = 0; i < 4; i++)
    {
        result[k + 1].X = result[k].X + median[0][1];
        result[k + 1].Y = result[k].Y + median[1][1];
        result[k + 1].Z = result[k].Z + median[2][1];
        result[k + 1].delta_0 = result[k].delta_0 + median[3][1];
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
        //����
        updata(median);
        //�ж�����
        te = sqrt( pow(median[0][0] - result[k].X, 2) + pow(median[1][0] - result[k].Y, 2) + pow(median[2][0] - result[k].Z, 2) );
        k++;
    }while(te > 1e-8);

}
