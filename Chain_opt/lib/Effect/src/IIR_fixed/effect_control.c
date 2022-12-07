#ifndef __EFFECT_CONTROL_H__
#define __EFFECT_CONTROL_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"
#include <math.h>

#define M_PI  (3.141592654)
#define IIR_OFSET 30
typedef struct {
    my_sint32 a0;
    my_sint32 a1;
    my_sint32 a2;

    my_sint32 b0;
    my_sint32 b1;
    my_sint32 b2;
}tIIR_coeffs;

typedef struct {
    float freq;
    float Q_fac;
    float sampleRate;
}tIIR_params;

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
    size_t*     coeffs_bytes)
{
    *params_bytes = sizeof(tIIR_params);
    *coeffs_bytes = sizeof(tIIR_coeffs);
}


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
    uint32_t    sample_rate)
{
    tIIR_params* init_params;
    tIIR_coeffs* init_coeffs;

    init_params = (tIIR_params*)params;
    init_coeffs = (tIIR_coeffs*)coeffs;

    init_params->freq = 0.0;
    init_params->Q_fac = 0.0;
    init_params->sampleRate = 0.0;

    init_coeffs->a0 = 0;
    init_coeffs->a1 = 0;
    init_coeffs->a2 = 0;

    init_coeffs->b0 = 0;
    init_coeffs->b1 = 0;
    init_coeffs->b2 = 0;
}


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
    float       value)
{
    tIIR_params* set_params;

    set_params = (tIIR_params*)params;

    switch (id)
    {
        case 0:
        {
            set_params->freq = value;
            break;
        }

        case 1:
        {
            set_params->Q_fac = value;
            break;
        }

        case 2:
        {
            set_params->sampleRate = value;
            break;
        }

        default:
            break;
    }
}


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
    void*       coeffs)
{
    tIIR_params* update_params;
    tIIR_coeffs* update_coeffs;

    float alpha;
    float omega;

    update_params = (tIIR_params*)params;
    update_coeffs = (tIIR_coeffs*)coeffs;

    omega = (2 * M_PI * update_params->freq) / update_params->sampleRate;
    alpha = sinf(omega) / (2 * update_params->Q_fac);

    update_coeffs->a0 = float_To_Fixed(1.0 + alpha, IIR_OFSET);      
    update_coeffs->a1 = float_To_Fixed((-2.0) * cosf(omega), IIR_OFSET);
    update_coeffs->a2 = float_To_Fixed(1.0 - alpha, IIR_OFSET);

    update_coeffs->b0 = float_To_Fixed(1.0, IIR_OFSET);
    update_coeffs->b1 = float_To_Fixed((-2.0) * cosf(omega), IIR_OFSET);
    update_coeffs->b2 = float_To_Fixed(1.0, IIR_OFSET);
}


#endif
