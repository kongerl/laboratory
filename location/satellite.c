#include <stdio.h>
#include <string.h>
#include <math.h>
#include "all.h"


int locate(int i)
{
    //����黯ʱ��
    t = 0;
    t = t * 60 + 345600;
    Naturalization_Time(i);

    //0��������A
    semiMajor_axis(i);

     //1��ƽ�����ٶ� n
    angular_Velocity(i);

    //2���黯ʱ�� tk

      //3��ƽ����� Mk
     meanAnomaly(i);

      //4��ƫ����� Ek
     eccentricAnomanly(i);

     //5�������� Vk
     tureAnomaly(i);

     //6��������� faik
     argumentOfOperigee(i);

     //7���㶯������
     perturbation(i);

      //9�����ƽ������ϵ����
     coordinates(i);
     MEOsate(i);

    return 0;
}

int read_data()
{
    FILE *file = fopen("hour1520.txt", "r");  // ���ļ�

    if (file == NULL)
    {
        printf("�޷����ļ�\n");
        return 1;
    }
    int count = 0;
    char line[100];  // ����ÿ��������100���ַ�
    while (fgets(line, sizeof(line), file)) {
            if(    strstr(line, "C26 2023 06 01 00 00 00") || strstr(line, "C29 2023 06 01 00 00 00") || strstr(line, "C30 2023 06 01 00 00 00")
                   || strstr(line, "C35 2023 06 01 00 00 00") || strstr(line, "C36 2023 06 01 00 00 00") || strstr(line, "C39 2023 06 01 00 00 00")
                   || strstr(line, "C40 2023 06 01 00 00 00") || strstr(line, "C45 2023 06 01 00 00 00") )
            {
                do {
                    //1
                    sscanf(line, "%s %d %d %d %d %d %d %lf %lf %lf ", EPH[count].a, &EPH[count].year, &EPH[count].month, &EPH[count].day, &EPH[count].hour,&EPH[count].mintue, &EPH[count].second,&EPH[count].a_f0, &EPH[count].a_f1, &EPH[count].a_f2);
                    //printf("%s %d %d %d %d %d %d\t%.20lf \t %.20lf \t %.20lf\n", EPH[count].a, EPH[count].year, EPH[count].month, EPH[count].day, EPH[count].hour, EPH[count].mintue,EPH[count].second,EPH[count].a_f0, EPH[count].a_f1, EPH[count].a_f2);
                    fgets(line, sizeof(line), file);
                    //2
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].IDOE, &EPH[count].C_rs, &EPH[count].Delta_n0, &EPH[count].M0);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].IDOE, EPH[count].C_rs, EPH[count].Delta_n0, EPH[count].M0);
                    fgets(line, sizeof(line), file);
                    //3
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].C_uc, &EPH[count].e, &EPH[count].C_us, &EPH[count].Delta_A);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].C_uc, EPH[count].e, EPH[count].C_us, EPH[count].Delta_A);
                    fgets(line, sizeof(line), file);
                    //4
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].Toe, &EPH[count].C_ic, &EPH[count].OMEGA0, &EPH[count].C_is);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].Toe, EPH[count].C_ic, EPH[count].OMEGA0, EPH[count].C_is);
                    fgets(line, sizeof(line), file);
                    //5
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].i0, &EPH[count].C_rc, &EPH[count].omega, &EPH[count].OMEGA_DOT);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].i0, EPH[count].C_rc, EPH[count].omega,EPH[count].OMEGA_DOT);
                    fgets(line, sizeof(line), file);
                    //6
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].Idot, &EPH[count].data_bit, &EPH[count].BDT, &EPH[count].A_DOT);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].Idot, EPH[count].data_bit, EPH[count].BDT,EPH[count].A_DOT);
                    fgets(line, sizeof(line), file);
                    //7
                    sscanf(line, " %lf %lf %lf %lf ", &EPH[count].SV, &EPH[count].Health, &EPH[count].IGD, &EPH[count].ISC);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].SV, EPH[count].Health, EPH[count].IGD, EPH[count].ISC);
                    fgets(line, sizeof(line), file);
                    //8
                    sscanf(line, " %lf %lf %lf %lf", &EPH[count].Launch_time, &EPH[count].IDOC, &EPH[count].Delta_n0_dot, &EPH[count].type);
                    //printf("\t\t %.20lf \t %.20lf \t %.20lf \t %.20lf\n", EPH[count].Launch_time, EPH[count].IDOC,EPH[count].Delta_n0_dot, EPH[count].type);
                    fgets(line, sizeof(line), file);
                    count++;
                    //printf("COUNT = %d\n",count);
                } while (    strstr(line, "EPH C26 CNV1") || strstr(line, "EPH C29 CNV1") || strstr(line, "EPH C30 CNV1")
                             || strstr(line, "EPH C35 CNV1") || strstr(line, "EPH C36 CNV1") || strstr(line, "EPH C39 CNV1")
                             || strstr(line, "EPH C40 CNV1") || strstr(line, "EPH C45 CNV1") );
                if(count == 8) break;//hour�ļ���������
            }
        }
    fclose(file);  // �ر��ļ�
    return 0;
}


//0��������
double  semiMajor_axis(int i)
{
    if(!strcmp(EPH[i].a, "C39") || !strcmp(EPH[i].a, "C40"))
    {
        A0_MEO = A_ref_IG_G + EPH[i].Delta_A;
    }
    else
    {
        A0_MEO = A_ref_MEO + EPH[i].Delta_A;
    }

    Ak = A0_MEO + EPH[i].A_DOT * Tk;

    return 0;
}


//1�������������е�ƽ�����ٶ�n
double angular_Velocity(int i)//�㶯�������� ��n�� �����᳤����a
{
    angV = sqrt(GM / pow(A0_MEO,3) ) + EPH[i].Delta_n0 + 0.5 * EPH[i].Delta_n0_dot * Tk;

    return 0;
}

//2������黯ʱ��TK       �������
double Naturalization_Time(int i)
{

    Tk = t - EPH[i].Toe;

    Tk -= 14;
    Tk -= da[i].rou/C;

    if( Tk > 3.024e5 )
    {
        Tk -= 6.048e5;
    }
    else if( Tk < -3.024e5)
    {
        Tk += 6.048e5;
    }

    return 0;
}

//3������۲�ʱ�����ǵ�ƽ�����Mk
double meanAnomaly(int i)
{
    Mk = EPH[i].M0 + angV * Tk;

    return 0;
}

//4������ƫ�����Ek
double eccentricAnomanly(int i)
{
    double temp;
    temp = Mk;
    //�����˺�û����һ��
    for(int j = 0; j<100; j++)
    {
        Ek[i] = Mk + EPH[i].e * sin(temp); //��Ϊʲô��Mk����������M0�𣿣���
        temp = Ek[i];
    }

    return 0;
}

//5������������ Vk
double tureAnomaly(int i)
{
    double a = sqrt(1 - EPH[i].e * EPH[i].e) * sin(Ek[i]);
    double b = (cos(Ek[i]) - EPH[i].e);
    Vk = atan2( a , b );

    return 0;
}

//6�������������faik
double argumentOfOperigee(int i)
{
    faik = Vk + EPH[i].omega;

    return 0;
}

//7(8)�������㶯������
double perturbation(int i)
{
    double Detal_u,Detal_r,Detal_i;

    //�㶯��
    Detal_u = EPH[i].C_uc * cos( 2 *  faik) + EPH[i].C_us * sin( 2 * faik );
    Detal_r = EPH[i].C_rc * cos( 2 *  faik) + EPH[i].C_rs * sin( 2 * faik );
    Detal_i = EPH[i].C_ic * cos( 2 *  faik) + EPH[i].C_is * sin( 2 * faik );

    //�㶯����
    uk = faik + Detal_u;
    rk = Ak * ( 1 - EPH[i].e * cos(Ek[i]) ) + Detal_r;
    ik = EPH[i].i0 + EPH[i].Idot * Tk + Detal_i;

    return 0;
}

//9��������ƽ������ϵ����
double coordinates(int i)
{
    Xk = rk * cos( uk );
    Yk = rk * sin( uk );

    return 0;
}

//10��MEO����
double MEOsate(int i)
{
    double OMEGA_k;

    //�����㾫��
    OMEGA_k = EPH[i].OMEGA0 + (EPH[i].OMEGA_DOT - We) * Tk - We * EPH[i].Toe; //������
    // printf("�����㾫�Ȧ�k = %.20lf\n",OMEGA_k);
    //CGCS2000����λ��M��
    location[i].Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    location[i].Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    location[i].Zk = Yk * sin(ik);

}
double GEOsate(int i)
{
    double OMEGA_k;

    OMEGA_k = EPH[i].OMEGA0 + EPH[i].OMEGA_DOT * Tk - We * EPH[i].Toe;
    //printf("�����㾫�Ȧ�k = %.20lf\n",OMEGA_k);

    //��Ȼ����
    location[i].Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    location[i].Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    location[i].Zk = Yk * sin(ik);
    //CGCS2000
    location[i].Xk =  location[i].Xk * cos(We * Tk) + (location[i].Yk * cos(-5) + location[i].Zk * sin(-5) )* sin(We * Tk);
    location[i].Yk = (location[i].Yk * cos(-5) + location[i].Zk * sin(-5) )* cos(We * Tk) - location[i].Xk  * sin(We * Tk);
    location[i].Zk =  location[i].Zk * cos(-5) + location[i].Yk * sin(5);
}
