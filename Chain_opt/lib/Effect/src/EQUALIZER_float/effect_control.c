#ifndef __EFFECT_CONTROL_H__
#define __EFFECT_CONTROL_H__

#include <stdint.h>
#include <stddef.h>
#include "math.h"

#define M_PI  3.14159265358979323846
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
    Q = 2,
    TYPE = 3
}eParam;

typedef struct {
    uint32_t Id;
    double value;

}equalizer_param_set;

typedef struct {
    double sample_rate;
    equalizer_param_set  freq[BANDS];
    equalizer_param_set  gain[BANDS];
    equalizer_param_set     Q[BANDS];
    equalizer_param_set  type[BANDS];

}equalizer_params;

typedef struct {
    float  b0[BANDS];
    float  b1[BANDS];
    float  b2[BANDS];

    float  a0[BANDS];
    float  a1[BANDS];
    float  a2[BANDS];

}equalizer_coeffs;

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
    *params_bytes = sizeof(equalizer_params);
    *coeffs_bytes = sizeof(equalizer_coeffs);
    return 0;
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
    equalizer_params *init_params;
    equalizer_coeffs *init_coeffs;
    size_t i;

    init_params = (equalizer_params*)params;
    init_coeffs = (equalizer_coeffs*)coeffs;
    init_params->sample_rate = sample_rate;
    
    for (i = 0; i < (size_t)BANDS; i++)
    {
        init_params->freq[i].Id = 0.0 + (float)i * 4.0;
        init_params->gain[i].Id = 1.0 + (float)i * 4.0;
        init_params->Q[i].Id =    2.0 + (float)i * 4.0;
        init_params->type[i].Id = 3.0 + (float)i * 4.0;

        init_params->freq[i].value = 0.0;
        init_params->gain[i].value = 0.0;
        init_params->Q[i].value = 0.0;
        init_params->type[i].value = 0.0;

        init_coeffs->b0[i] = 0.0;
        init_coeffs->b1[i] = 0.0;
        init_coeffs->b2[i] = 0.0;

        init_coeffs->a0[i] = 0.0;
        init_coeffs->a1[i] = 0.0;
        init_coeffs->a2[i] = 0.0;
    }
    
    return 0;
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
    equalizer_params *set_params;
    size_t band;
    size_t param_num;

    set_params = (equalizer_params*)params;

    band = id >> 2;
    param_num = id - (4 * band);

    if (id != 40)
    {
        switch (param_num)
        {
        case FRIQUENCY:
        {
            set_params->freq[band].value = value;
            break;
        }

        case GAIN:
        {
            set_params->gain[band].value = value;
            break;
        }

        case Q:
        {
            set_params->Q[band].value = value;
            break;
        }

        case TYPE:
        {
            set_params->type[band].value = value;
            break;
        }

        default:
            break;
        }
    }
    else 
    {
        set_params->sample_rate = value;
    }
    
    return 0;
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
    double  A[10], sn[10], cs[10], alpha[10], beta[10], omega[10];
    double b0[10], b1[10], b2[10], a0[10], a1[10], a2[10];

    equalizer_params *update_params;
    equalizer_coeffs *update_coeffs;

    update_params = (equalizer_params*)params;
    update_coeffs = (equalizer_coeffs*)coeffs;

    for (size_t i = 0; i < 10; i++)
    {
        A[i]        = pow(10, update_params->gain[i].value / 40);
        omega[i]    = 2 * M_PI * update_params->freq[i].value / update_params->sample_rate;
        sn[i]       = sin(omega[i]);
        cs[i]       = cos(omega[i]);
        alpha[i]    = sn[i] / (2 * update_params->Q[i].value);
        beta[i]     = sqrt(A[i] + A[i]);

        switch ((int)update_params->type[i].value)
        {
            case LOWPASS:
            {
                b0[i] = (1.0 - cs[i]) / 2.0;
                b1[i] = 1.0 - cs[i];
                b2[i] = (1.0 - cs[i]) / 2.0;
                a0[i] = 1.0 + alpha[i];
                a1[i] = -2.0 * cs[i];
                a2[i] = 1.0 - alpha[i];
                break;
            }
            
            case HIGHPASS:
            {
                b0[i] = (1 + cs[i]) / 2.0;
                b1[i] = -(1 + cs[i]);
                b2[i] = (1 + cs[i]) / 2.0;
                a0[i] = 1 + alpha[i];
                a1[i] = -2 * cs[i];
                a2[i] = 1 - alpha[i];
                break;
            }
            
            case PEAK:
            {
                b0[i] = 1 + (alpha[i] * A[i]);
                b1[i] = -2 * cs[i];
                b2[i] = 1 - (alpha[i] * A[i]);
                a0[i] = 1 + (alpha[i] / A[i]);
                a1[i] = -2 * cs[i];
                a2[i] = 1 - (alpha[i] / A[i]);
                break;
            }
            
            case LOWSHELL:
            {
                b0[i] = A[i] * ((A[i] + 1) - (A[i] - 1) * cs[i] + beta[i] * sn[i]);
                b1[i] = A[i] * ((A[i] - 1) - (A[i] + 1) * cs[i]) * 2;
                b2[i] = A[i] * ((A[i] + 1) - (A[i] - 1) * cs[i] - beta[i] * sn[i]);
                a0[i] = (A[i] + 1) + (A[i] - 1) * cs[i] + beta[i] * sn[i];
                a1[i] = ((A[i] - 1) + (A[i] + 1) * cs[i]) * (-2);
                a2[i] = (A[i] + 1) + (A[i] - 1) * cs[i] - beta[i] * sn[i];
                break;
            }
            
            case HIGHSHELL:
            {
                b0[i] = A[i] * ((A[i] + 1) + (A[i] - 1) * cs[i] + beta[i] * sn[i]);
                b1[i] = A[i] * ((A[i] - 1) + (A[i] + 1) * cs[i]) * (-2);
                b2[i] = A[i] * ((A[i] + 1) + (A[i] - 1) * cs[i] - beta[i] * sn[i]);
                a0[i] = (A[i] + 1) - (A[i] - 1) * cs[i] + beta[i] * sn[i];
                a1[i] = ((A[i] - 1) - (A[i] + 1) * cs[i]) * 2;
                a2[i] = (A[i] + 1) - (A[i] - 1) * cs[i] - beta[i] * sn[i];
                break;
            }
            
            case DISABLE:
            {
                b0[i] = 0;
                b1[i] = 0;
                b2[i] = 0;
                a0[i] = 0;
                a1[i] = 0;
                a2[i] = 0;
                break;
            }
           
        }

        // a0[i] /= a0[i];
        b0[i] /= a0[i];
        b1[i] /= a0[i];
        b2[i] /= a0[i];
        a1[i] /= a0[i];
        a2[i] /= a0[i];

        update_coeffs->a0[i] = (float)(a0[i] / 8);
        update_coeffs->a1[i] = (float)(a1[i] / 8);
        update_coeffs->a2[i] = (float)(a2[i] / 8);
        update_coeffs->b0[i] = (float)(b0[i] / 8);
        update_coeffs->b1[i] = (float)(b1[i] / 8);
        update_coeffs->b2[i] = (float)(b2[i] / 8);

    }

    return 0;
}


#endif
