#include "arith.h"

int32_t Arith_max(int32_t x, int32_t y){
    return x > y  ? x : y;
}

int32_t Arith_min(int32_t x, int32_t y){
    return x > y  ? y : x;
}

int32_t Arith_div(int32_t x, int32_t y){
    if(((-13/5 == 2)) && ((x>0) != (y<0)) && x%y!=0 ){
        return x/y-1;
    }
    else return x/y;
}


int32_t Arith_mod(int32_t x,int32_t y){
    if(((-13/5 == 2)) && ((x>0) != (y<0)) && x%y!=0 ){
        return x%y+y;
    }
    else return x%y;
}

int32_t Arith_ceiling(int32_t x, int32_t y){
    return (Arith_div(x,y)+(x%y!=0));
}

int32_t Arith_floor(int32_t x, int32_t y){
    return Arith_div(x,y);
}
