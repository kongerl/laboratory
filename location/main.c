#include "all.h"
#include "satellite.c"
#include "UNB3m.c"
#include "bdgim.c"
#include "Spp.c"



int main() {
    const int num = 8;
    //读文件
    read_data();
    printf("读文件成功\n");

    //卫星位置计算
    for(int i = 0; i < num; i++)
    {
        locate(i);
    }
    printf("卫星坐标计算完成\n");

   // Spp();

    return 0;
}
