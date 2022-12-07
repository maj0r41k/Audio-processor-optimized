#include <stdint.h>
#include <stddef.h>
#include "math.h"

#include "equalizer_control.h"

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
int32_t equalizer_control_initialize(
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

        init_params->freq[i] = 0.0;
        init_params->gain[i] = 0.0;
        init_params->Q_f[i] = 0.0;
        init_params->type[i] = 0.0;

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
int32_t equalizer_set_parameter(
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
            set_params->freq[band] = value;
            break;
        }

        case GAIN:
        {
            set_params->gain[band] = value;
            break;
        }

        case QF:
        {
            set_params->Q_f[band] = value;
            break;
        }

        case TYPE:
        {
            set_params->type[band] = value;
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
int32_t equalizer_update_coeffs(
    void const* params,
    void*       coeffs)
{
    double  A[BANDS],  sn[BANDS], cs[BANDS], alpha[BANDS], beta[BANDS], omega[BANDS];
    double  b0[BANDS], b1[BANDS], b2[BANDS], a0[BANDS],    a1[BANDS],   a2[BANDS];

    equalizer_params *update_params;
    equalizer_coeffs *update_coeffs;

    update_params = (equalizer_params*)params;
    update_coeffs = (equalizer_coeffs*)coeffs;

    for (size_t i = 0; i < 10; i++)
    {
        A[i]        = pow(10, update_params->gain[i] / 40);
        omega[i]    = 2 * M_PI_eq * update_params->freq[i] / update_params->sample_rate;
        sn[i]       = sin(omega[i]);
        cs[i]       = cos(omega[i]);
        alpha[i]    = sn[i] / (2 * update_params->Q_f[i]);
        beta[i]     = sqrt(A[i] + A[i]);

        switch ((int)update_params->type[i])
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

        update_coeffs->a0[i] = (float)(a0[i]);
        update_coeffs->a1[i] = (float)(a1[i]);
        update_coeffs->a2[i] = (float)(a2[i]);
        update_coeffs->b0[i] = (float)(b0[i]);
        update_coeffs->b1[i] = (float)(b1[i]);
        update_coeffs->b2[i] = (float)(b2[i]);

    }

    return 0;
}
