#ifndef __EQUALIZER_CONTROL_H__
#define __EQUALIZER_CONTROL_H__

#include <stdint.h>
#include <stddef.h>

#include "fixedpoint.h"

#define M_PI_eq  3.14159265358979323846


#define BANDS 10

typedef enum {
    LOWPASS = 0,
    HIGHPASS = 1,
    PEAK = 2,
    LOWSHELL = 3,
    HIGHSHELL = 4,
    DISABLE = 5,

    TOTAL_REGIMES = 6
}eFilterType;

typedef enum {
    FRIQUENCY = 0,
    GAIN = 1,
    QF = 2,
    TYPE = 3
}eParam;

typedef struct {
    float L;
    float R;
}tStereo_eq;

typedef struct {
    float  freq[BANDS];
    float  gain[BANDS];
    float   Q_f[BANDS];
    float  type[BANDS];
    float sample_rate;

}equalizer_params;      // parameters

typedef struct {
    float  b0[BANDS];
    float  b1[BANDS];
    float  b2[BANDS];

    float  a0[BANDS];
    float  a1[BANDS];
    float  a2[BANDS];

}equalizer_coeffs;      // coefficients

typedef struct {
    my_float x0[BANDS];
    my_float x1[BANDS];
    my_float x2[BANDS];
    my_float y0[BANDS];
    my_float y1[BANDS];
    my_float y2[BANDS];
}equalizer_states;


/*******************************************************************************
 * Provides with the required data sizes for parameters and coefficients.
 *   It is caller responsibility to allocate enough memory (bytes) for them.
 * 
 * @param[out] params_bytes   required data size for storing parameters
 * @param[out] coeffs_bytes   required data size for storing coefficients
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t equalizer_control_get_sizes(
    size_t*     params_bytes,
    size_t*     coeffs_bytes);


/*******************************************************************************
 * Initializes params, coeffs and states to default values for the requested SR.
 * 
 * @param[in]  params       pointer to the pre-allocated params
 * @param[in]  coeffs       pointer to the pre-allocated coeffs
 * @param[in]  sample_rate  sampling rate
 * 
 * @return 0 if effect is initialized, non-zero error code otherwise
 ******************************************************************************/
int32_t equalizer_control_initialize(
    void*       params,
    void*       coeffs,
    uint32_t    sample_rate);


/*******************************************************************************
 * Set single parameter and calculate corresponding coefficients.
 * 
 * @param[in] params    initialized params
 * @param[in] id        parameter ID
 * @param[in] value     parameter value
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t equalizer_set_parameter(
    void*       params,
    int32_t     id,
    float       value);


/*******************************************************************************
 * Calculate coefficients, corresponding to the parameters.
 * 
 * @param[in] params    initialized params
 * @param[in] id        parameter ID
 * @param[in] value     parameter value
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t equalizer_update_coeffs(
    void const* params,
    void*       coeffs);


#endif