#ifndef __EFFECT_CONTROL_H__
#define __EFFECT_CONTROL_H__

#include <stdint.h>
#include <stddef.h>

#include "compressor_4ch_control.h"
#include "crossover_control.h"
#include "equalizer_control.h"
#include "fixedpoint.h"

typedef struct {
    my_float L;
    my_float R;
}stereo;

typedef struct {
    void* band_1;
    void* band_2;
    void* band_3;
    void* band_4;
}band_proc;

typedef struct {
    void* audio;
    band_proc cross_b;  // * deleted
}audio_proc;

typedef struct {    
    comprssor_4ch_params compr_4ch_par;
    crossover_params cross_par;
    equalizer_params equal_par;
}chain_params;

typedef struct {
    compressor_4ch_coeffs compr_4ch_coef;
    crossover_coeffs cross_coef;
    equalizer_coeffs equal_coef;
}chain_coeffs;

typedef struct {
    compressor_4ch_states compr_4ch_st;
    crossover_states cross_st;
    equalizer_states equal_st;
}chain_states;

/*******************************************************************************
 * Provides with the required data sizes for parameters and coefficients.
 *   It is caller responsibility to allocate enough memory (bytes) for them.
 * 
 * @param[out] params_bytes   required data size for storing parameters
 * @param[out] coeffs_bytes   required data size for storing coefficients
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t effect_control_get_sizes(
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
int32_t effect_control_initialize(
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
int32_t effect_set_parameter(
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
int32_t effect_update_coeffs(
    void const* params,
    void*       coeffs);


#endif