#ifndef __CROSSOVER_4CH_CONTROL_H__
#define __CROSSOVER_4CH_CONTROL_H__

#include <stdint.h>
#include <stddef.h>

#include "fixedpoint.h"

#define M_e    ((float)2.71828182846)

typedef struct {
    float L;
    float R;
}tStereo_compr;

typedef struct {
    void* band_1;
    void* band_2;
    void* band_3;
    void* band_4;
}band_buf_comp;

typedef struct {
    void* audio;
    band_buf_comp cross_b;
}audio_buf_comp;

typedef struct {
    float threshold;
    float ratio;
    float tauAttack;
    float tauRelease;
    float makeUpGain;
    float samplerate;
    float tauEnvAtt;
    float tauEnvRel;

    uint32_t bpass;
}comprssor_params;

typedef struct {
    comprssor_params comp_ch_1_par;
    comprssor_params comp_ch_2_par;
    comprssor_params comp_ch_3_par;
    comprssor_params comp_ch_4_par;
}comprssor_4ch_params;


typedef struct {
    my_float threshold;

    my_float ratio;
    my_float _rep_ratio;

    my_float alphaAttack;
    my_float alphaRelease;
    my_float _one_g_att;
    my_float _one_g_rel;

    my_float attackEnv;
    my_float releaseEnv;
    my_float _one_env_att;
    my_float _one_env_rel;

    my_float makeUpGain;
    my_float samplerate;

    uint32_t bpass[4];
}compressor_4ch_coeffs;


typedef struct {
    my_float xL;
    my_float xR;
    my_float g_c;
    my_float g_s;
    my_float g_sPrev;
    my_float g_m;
    my_float envelope;
    my_float envelope_prev;
}compressor_4ch_states;

/*******************************************************************************
 * Provides with the required data sizes for parameters and coefficients.
 *   It is caller responsibility to allocate enough memory (bytes) for them.
 * 
 * @param[out] params_bytes   required data size for storing parameters
 * @param[out] coeffs_bytes   required data size for storing coefficients
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t compressor_4ch_control_get_sizes(
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
int32_t compressor_4ch_control_initialize(
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
int32_t compressor_4ch_set_parameter(
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
int32_t compressor_4ch_update_coeffs(
    void const* params,
    void*       coeffs);


#endif