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
        locate(i);
    }
    printf("��������������\n");

   // Spp();

    return 0;
}
