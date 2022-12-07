#ifndef __FIXEDPOINT_H_
#define __FIXEDPOINT_H_

#include <stdint.h>
#include <stdio.h>
#include "math.h"
#include "immintrin.h"
#include <xmmintrin.h>

#define SCALAR
//#define SSE



#define MAX_VAL_8       127
#define MIN_VAL_8       -128

#define MAX_VAL_16      32767
#define MIN_VAL_16      -32768

#define MAX_VAL_32      2147483647
#define MIN_VAL_32      -2147483648

#define MAX_VAL_64      9223372036854775807
#define MIN_VAL_64      -9223372036854775808

#define FRACTION_BASE   31
#define ESTIMATE_Q20    20
#define LOG_BASE        26
#define LOG_IDX_OFFSET  22
#define POW_IDX_OFFSET  22
#define POW_BASE        26

#define Q31_TO_Q26      5

#define Q20_SHIFT       11

#define EMPTY_MASK_32   0x00000000
#define EMPTY_MASK_64   0x0000000000000000
#define ESTIMATE_VAL    0x80000
#define ONE             0xFFFFF
#define ROUNDING_BIT_32 0x80000000
#define ROUNDING_BIT_64 0x80000000

#define Q31 31
#define Q30 30
#define Q29 29
#define Q28 28
#define Q27 27
#define Q26 26
#define Q25 25
#define Q24 24
#define Q23 23
#define Q22 22
#define Q21 21
#define Q20 20
#define Q19 19
#define Q18 18
#define Q17 17
#define Q16 16
#define Q15 15
#define Q14 14
#define Q13 13

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

typedef unsigned char           my_uint8;
typedef unsigned short          my_uint16;
typedef unsigned int            my_uint32;
typedef unsigned long int       my_uint64;

typedef signed char             my_sint8;
typedef signed short            my_sint16;
typedef signed int              my_sint32;
typedef signed long long int    my_sint64;

#ifdef SCALAR
    typedef float my_float;
#endif //SCALAR

#ifdef SSE
    typedef __m128 my_float;
#endif // 

my_sint32 add32(const my_sint32 a, const my_sint32 b);
my_sint32 sub32(const my_sint32 a, const my_sint32 b);
my_sint32 mul32(const my_sint32 a, const my_sint32 b);
my_sint32 mac32(const my_sint32 a, const my_sint32 b, my_sint32 c);
my_sint64 mac64(const my_sint32 a, const my_sint32 b, const my_sint64 c);
my_sint32 msub32(const my_sint32 a, const my_sint32 b, my_sint32 c);
my_sint64 msub64(const my_sint32 a, const my_sint32 b, const my_sint64 c);
my_sint32 abs32(const my_sint32 a);
my_sint32 neg32(const my_sint32 a);
my_sint32 lsh32(const my_sint32 a, const my_sint32 b);
my_sint64 lsh64(const my_sint64 a, const my_sint64 b);
my_sint32 rsh32(const my_sint32 a, const my_sint32 b);
my_sint32 div32(const my_sint32 numenator, const my_sint32 denuminator);
my_sint32 div32_1_x(const my_sint32 denuminator, const my_sint32 Q);

my_sint64 add64(const my_sint64 a, const my_sint64 b);
my_sint64 sub64(const my_sint64 a, const my_sint64 b);
my_sint64 mul64(const my_sint64 a, const my_sint64 b, const my_sint64 shift);
my_sint32 mul32_q(const my_sint32 a, const my_sint32 b, const my_sint32 shift);

my_sint32 saturation32(my_sint32* sum, my_sint32* term);
my_sint64 saturation64(my_sint64* sum, my_sint64* term);

my_sint32 float_To_Fixed(my_float floatNum, my_uint8 shift);
my_float fixed_To_Float(my_sint32 fixedNum, my_uint8 shift);

extern my_float add_f(const my_float a, const my_float b);
extern my_float sub_f(const my_float a, const my_float b);
extern my_float mul_f(const my_float a, const my_float b);
extern my_float mac_f(const my_float a, const my_float b, const my_float c);
extern my_float msub_f(const my_float a, const my_float b, const my_float c);
extern my_float abs_f(const my_float a);
extern my_float neg_f(const my_float a);
extern my_float div_f(const my_float a, const my_float b);
extern float pow_f(const float a, const float b);

my_sint32 log2x(my_sint32 a);
my_sint32 pow2x(my_sint32 a);
my_sint32 pow10x(my_sint32 a);
my_sint32 my_pow(my_sint32 a, my_sint32 x);
my_sint32 log10x(my_sint32 a);

#endif // !__FIXEDPOINT_H_
