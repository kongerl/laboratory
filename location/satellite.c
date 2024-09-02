#include <stdio.h>
#include <string.h>
#include <math.h>
#include "all.h"


int locate(int i)
{
    //计算归化时间
    t = 0;
    t = t * 60 + 345600;
    Naturalization_Time(i);

    //0、长半轴A
    semiMajor_axis(i);

     //1、平均角速度 n
    angular_Velocity(i);

    //2、归化时间 tk

      //3、平近点角 Mk
     meanAnomaly(i);

      //4、偏近点角 Ek
     eccentricAnomanly(i);

     //5、真近点角 Vk
     tureAnomaly(i);

     //6、升交距角 faik
     argumentOfOperigee(i);

     //7、摄动改正量
     perturbation(i);

      //9、轨道平面坐标系坐标
     coordinates(i);
     MEOsate(i);

    return 0;
}

int read_data()
{
    FILE *file = fopen("hour1520.txt", "r");  // 打开文件

    if (file == NULL)
    {
        printf("无法打开文件\n");
        return 1;
    }
    int count = 0;
    char line[100];  // 假设每行最多包含100个字符
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
                if(count == 8) break;//hour文件中有两个
            }
        }
    fclose(file);  // 关闭文件
    return 0;
}


//0、长半轴
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


//1、计算卫星运行的平均角速度n
double angular_Velocity(int i)//摄动改正参数 △n， 长半轴长根号a
{
    angV = sqrt(GM / pow(A0_MEO,3) ) + EPH[i].Delta_n0 + 0.5 * EPH[i].Delta_n0_dot * Tk;

    return 0;
}

//2、计算归化时间TK       最后问题
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

//3、计算观测时刻卫星的平近点角Mk
double meanAnomaly(int i)
{
    Mk = EPH[i].M0 + angV * Tk;

    return 0;
}

//4、计算偏近点角Ek
double eccentricAnomanly(int i)
{
    double temp;
    temp = Mk;
    //迭代了和没迭代一样
    for(int j = 0; j<100; j++)
    {
        Ek[i] = Mk + EPH[i].e * sin(temp); //这为什么是Mk？？？不是M0吗？？？
        temp = Ek[i];
    }

    return 0;
}

//5、计算真近点角 Vk
double tureAnomaly(int i)
{
    double a = sqrt(1 - EPH[i].e * EPH[i].e) * sin(Ek[i]);
    double b = (cos(Ek[i]) - EPH[i].e);
    Vk = atan2( a , b );

    return 0;
}

//6、计算升交距角faik
double argumentOfOperigee(int i)
{
    faik = Vk + EPH[i].omega;

    return 0;
}

//7(8)、计算摄动改正量
double perturbation(int i)
{
    double Detal_u,Detal_r,Detal_i;

    //摄动量
    Detal_u = EPH[i].C_uc * cos( 2 *  faik) + EPH[i].C_us * sin( 2 * faik );
    Detal_r = EPH[i].C_rc * cos( 2 *  faik) + EPH[i].C_rs * sin( 2 * faik );
    Detal_i = EPH[i].C_ic * cos( 2 *  faik) + EPH[i].C_is * sin( 2 * faik );

    //摄动改正
    uk = faik + Detal_u;
    rk = Ak * ( 1 - EPH[i].e * cos(Ek[i]) ) + Detal_r;
    ik = EPH[i].i0 + EPH[i].Idot * Tk + Detal_i;

    return 0;
}

//9、计算轨道平面坐标系坐标
double coordinates(int i)
{
    Xk = rk * cos( uk );
    Yk = rk * sin( uk );

    return 0;
}

//10、MEO卫星
double MEOsate(int i)
{
    double OMEGA_k;

    //升交点精度
    OMEGA_k = EPH[i].OMEGA0 + (EPH[i].OMEGA_DOT - We) * Tk - We * EPH[i].Toe; //有问题
    // printf("升交点精度Ωk = %.20lf\n",OMEGA_k);
    //CGCS2000（单位：M）
    location[i].Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    location[i].Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    location[i].Zk = Yk * sin(ik);

}
double GEOsate(int i)
{
    double OMEGA_k;

    OMEGA_k = EPH[i].OMEGA0 + EPH[i].OMEGA_DOT * Tk - We * EPH[i].Toe;
    //printf("升交点精度Ωk = %.20lf\n",OMEGA_k);

    //自然惯性
    location[i].Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    location[i].Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    location[i].Zk = Yk * sin(ik);
    //CGCS2000
    location[i].Xk =  location[i].Xk * cos(We * Tk) + (location[i].Yk * cos(-5) + location[i].Zk * sin(-5) )* sin(We * Tk);
    location[i].Yk = (location[i].Yk * cos(-5) + location[i].Zk * sin(-5) )* cos(We * Tk) - location[i].Xk  * sin(We * Tk);
    location[i].Zk =  location[i].Zk * cos(-5) + location[i].Yk * sin(5);
}
