#ifndef __CROSSOVER_DATATYPES_H_
#define __CROSSOVER_DATATYPES_H_

#include "fixedpoint.h"

#define M_PI        3.1415926535f
#define Q           0.70710678118f
#define BAND_NUM    4u

typedef enum {
    LOW_CUTOFF = 0,
    MID_CUTOFF = 1,
    HIG_CUTOFF = 2,
    SAMPLERATE = 3

}crossover_params_enum;

typedef enum {
    PAR_1 = 0,
    PAR_2 = 1,
    SERIAL = 2,
    NORM_LOW = 0,
    NORM_HIG = 1

}crossover_states_enum;

typedef struct {
    my_float L;
    my_float R;
}tStereo;

typedef struct {
    tStereo LP;
    tStereo HP;
}band_split;

typedef struct {
    tStereo band_1;
    tStereo band_2;
    tStereo band_3;
    tStereo band_4;
}split_band;

/*
*   PARAMETERS
*/
typedef struct {
    my_float cutoff_freq[BAND_NUM - 1];
    my_float samplerate;
}crossover_params;


/*
*   COEFFICIENTS
*/
typedef struct {
    my_float k0;
    my_float k1;
    my_float k2;

    my_float negk0;
    my_float negk1;
    my_float negk2;
}band_coeffs;

typedef struct {
    band_coeffs band[BAND_NUM - 1];
}crossover_coeffs;


/*
*   STATES
*/

//states for one apf uint contains apf of 1-st and 2-nd order


typedef struct {
    tStereo in;

    tStereo cd1_ord1[2];
    tStereo cd1_delay_1[3];
    tStereo cd1_delay_2[3];

    tStereo norm_ord1[2];
    tStereo norm_delay_1[2];
    tStereo norm_delay_2[2];

    tStereo cd2_ord1[2];
    tStereo cd2_delay_1[3];
    tStereo cd2_delay_2[3];

    tStereo cd3_ord1[2];
    tStereo cd3_delay_1[3];
    tStereo cd3_delay_2[3];

}crossover_states; // size = 18*3 + 2*6 = 54 + 12 = 66 floats



#endif // !__CROSSOVER_DATATYPES_H_
