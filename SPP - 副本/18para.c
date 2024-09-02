#include "spp.h"
#include "UNB3m.c"
#include "18para.h"
#include<stdio.h>
//#include<string.h>
#include<math.h>



int locate(int a)
{int i = 0;
    i = a;


    //参考时间t
    printf("请输入参考时刻（范围为0-60，单位每分钟）：");
    scanf("%d",&t);
    t = t * 60 + 345600;
    Naturalization_Time(i);

    //0、长半轴A
    semiMajor_axis(i);
    printf("0.\n");
    printf("长半轴A = %.20lf\n",Ak );
    //1、平均角速度 n
    angular_Velocity(i);
    printf("1.\n");
    printf("平均角速度n = %.20lf\n",angV );

    //2、归化时间 tk
    printf("2.\n");
    printf("归化时间tk = %.20lf\n",Tk );

    //3、平近点角 Mk
    meanAnomaly(i);
    printf("3.\n");
    printf("平近点角Mk = %.20lf\n",Mk );

    //4、偏近点角 Ek
    eccentricAnomanly(i);
    printf("4.\n");
    printf("偏近点角Ek = %.20lf\n",Ek );

    //5、真近点角 Vk
    tureAnomaly(i);
    printf("5.\n");
    printf("真近点角Vk = %.20lf\n",Vk );

    //6、升交距角 faik
    argumentOfOperigee(i);
    printf("6.\n");
    printf("升交距角faik = %.20lf\n", faik);

    //7、摄动改正量
    perturbation(i);
    printf("7(8).\n");
    printf("升交距角uk = %.20lf\n", uk);
    printf("升交距角rk = %.20lf\n", rk);
    printf("升交距角ik = %.20lf\n", ik);

    //9、轨道平面坐标系坐标
    coordinates(i);
    printf("9.\n");
    printf("轨道平面坐标Xk = %.20lf\n", Xk);
    printf("轨道平面坐标Yk = %.20lf\n", Yk);

    //10、MEO/IGSO卫星坐标
    MEOsate(i);
    printf("10.\n");
    printf("自定义惯性系Xk = %.20lf\n", meo.user_Xk);
    printf("自定义惯性系Yk = %.20lf\n", meo.user_Yk);
    printf("自定义惯性系Zk = %.20lf\n", meo.user_Zk);

    printf("CGCS2000Xk = %.20lf\n", meo.CGCS_Xk);
    printf("CGCS2000Yk = %.20lf\n", meo.CGCS_Yk);
    printf("CGCS2000Zk = %.20lf\n", meo.CGCS_Zk);
    //GEO卫星
//	GEOsate();
//	printf("10.\n");
//	printf("自定义惯性系Xk = %.20lf\n", geo.user_Xk);
//	printf("自定义惯性系Yk = %.20lf\n", geo.user_Yk);
//	printf("自定义惯性系Zk = %.20lf\n", geo.user_Zk);
//
//	printf("CGCS2000Xk = %.20lf\n", geo.CGCS_Xk/1000);
//	printf("CGCS2000Yk = %.20lf\n", geo.CGCS_Yk/1000);
//	printf("CGCS2000Zk = %.20lf\n", geo.CGCS_Zk/1000);



    return 0;
}


//0、长半轴
double  semiMajor_axis(int a)
{
    A0_MEO = A_ref_MEO + EPH[i].Delta_A;

    Ak = A0_MEO + EPH[i].A_DOT * Tk;

    return 0;
}


//1、计算卫星运行的平均角速度n
double angular_Velocity(int a)//摄动改正参数 △n， 长半轴长根号a
{
//	double angV0;//平局角速度V

//	angV0 = sqrt(GM) / pow(EPH[i].sqrtA,3);
    angV = sqrt(GM / pow(A0_MEO,3) ) + EPH[i].Delta_n0 + 0.5 * EPH[i].Delta_n0_dot * Tk;


    return 0;
}

//2、计算归化时间TK       最后问题
double Naturalization_Time(int a)
{

    Tk = t - EPH[i].Toe;

    Tk -= 14;

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
double meanAnomaly(int a)
{
    Mk = EPH[i].M0 + angV * Tk;

    return 0;
}

//4、计算偏近点角Ek
double eccentricAnomanly(int a)
{
    double temp;
    temp = Mk;
    //迭代了和没迭代一样
    for(int i = 0; i<100; i++)
    {
        Ek = Mk + EPH[i].e * sin(temp); //这为什么是Mk？？？不是M0吗？？？
        temp = Ek;
    }

    return 0;
}

//5、计算真近点角 Vk
double tureAnomaly(int a)
{
    Vk = atan( sqrt(1 - EPH[i].e * EPH[i].e) * sin(Ek) / (cos(Ek) - EPH[i].e) );

    return 0;
}

//6、计算升交距角faik
double argumentOfOperigee(int a)
{

    faik = Vk + EPH[i].omega;

    return 0;
}

//7(8)、计算摄动改正量
double perturbation(int a)
{
    double Detal_u,Detal_r,Detal_i;


    //摄动量
    Detal_u = EPH[i].C_uc * cos( 2 *  faik) + EPH[i].C_us * sin( 2 * faik );
    Detal_r = EPH[i].C_rc * cos( 2 *  faik) + EPH[i].C_rs * sin( 2 * faik );
    Detal_i = EPH[i].C_ic * cos( 2 *  faik) + EPH[i].C_is * sin( 2 * faik );
    printf("升交距角uk摄动改正量 = %.20lf\n", Detal_u);
    printf("升交距角rk摄动改正量 = %.20lf\n", Detal_r);
    printf("升交距角ik摄动改正量 = %.20lf\n", Detal_i);

    //摄动改正
    uk = faik + Detal_u;
    rk = Ak * ( 1 - EPH[i].e * cos(Ek) ) + Detal_r;
    ik = EPH[i].i0 + Detal_i + EPH[i].Idot * Tk;

    return 0;
}

//9、计算轨道平面坐标系坐标
double coordinates(int a)
{
    Xk = rk * cos( uk );
    Yk = rk * sin( uk );

    return 0;
}

//10、MEO卫星
double MEOsate(int a)
{
    double OMEGA_k;

    //升交点精度
    OMEGA_k = EPH[i].OMEGA0 + (EPH[i].OMEGA_DOT - We) * Tk - We * EPH[i].Toe; //有问题
    printf("升交点精度Ωk = %.20lf\n",OMEGA_k);
    //CGCS2000（单位：M）
    meo.user_Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    meo.user_Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    meo.user_Zk = Yk * sin(ik);
    //CGCS2000（单位：KM）
    meo.CGCS_Xk = meo.user_Xk / 1000;
    meo.CGCS_Yk = meo.user_Yk / 1000;
    meo.CGCS_Zk = meo.user_Zk / 1000;
}
//
double GEOsate(int a)
{
    double OMEGA_k;

    OMEGA_k = EPH[i].OMEGA0 + EPH[i].OMEGA_DOT * Tk - We * EPH[i].Toe;
    printf("升交点精度Ωk = %.20lf\n",OMEGA_k);

    //自然惯性
    geo.user_Xk = Xk * cos(OMEGA_k) - Yk * cos(ik) * sin(OMEGA_k);
    geo.user_Yk = Xk * sin(OMEGA_k) + Yk * cos(ik) * cos(OMEGA_k);
    geo.user_Zk = Yk * sin(ik);
    //CGCS2000
    geo.CGCS_Xk = geo.user_Xk * cos(We * Tk) + (geo.user_Yk * cos(-5) + geo.user_Zk * sin(-5) )* sin(We * Tk);
    geo.CGCS_Yk = (geo.user_Yk * cos(-5) + geo.user_Zk * sin(-5) )* cos(We * Tk) - geo.user_Xk * sin(We * Tk);
    geo.CGCS_Zk = geo.user_Zk * cos(-5) + geo.user_Yk * sin(5);
}

