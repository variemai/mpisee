#ifndef _ARITH_H
#define _ARITH_H
#include <stdint.h>
extern int32_t Arith_max(int32_t x, int32_t y);
extern int32_t Arith_min(int32_t x, int32_t y);
extern int32_t Arith_div(int32_t x, int32_t y);
extern int32_t Arith_mod(int32_t x, int32_t y);
extern int32_t Arith_ceiling(int32_t x, int32_t y);
extern int32_t Arith_floor(int32_t x, int32_t y);
#endif
