#include "all.h"
#include "satellite.c"
#include "UNB3m.c"
#include "bdgim.c"
#include "Spp.c"



int main() {
    const int num = 8;
    //���ļ�
    read_data();
    printf("���ļ��ɹ�\n");

    //����λ�ü���
    for(int i = 0; i < num; i++)
    {
        //printf("i = %d\n",i);
        locate(i);
    }
    printf("��������������\n");

 /*   //�����
    for(int i = 0; i < num; i++)
    {
        printf("bdgim[%d] = %lf\n", i,bdgim3(i));
    }

    //������
    for(int i = 0; i < num; i++)
    {
        printf("unb3m[%d] = %lf\n", i,UNB3(i));
    }*/


    //SPP
    Spp();

    return 0;
}
