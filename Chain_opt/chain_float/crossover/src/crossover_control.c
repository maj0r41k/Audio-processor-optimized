#include <stdint.h>
#include <stddef.h>

#include "crossover_control.h"

#include "math.h"

/*******************************************************************************
 * Provides with the required data sizes for parameters and coefficients.
 *   It is caller responsibility to allocate enough memory (bytes) for them.
 * 
 * @param[out] params_bytes   required data size for storing parameters
 * @param[out] coeffs_bytes   required data size for storing coefficients
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t crossover_control_get_sizes(
    size_t*     params_bytes,
    size_t*     coeffs_bytes)
{
    *params_bytes = sizeof(crossover_params);
    *coeffs_bytes = sizeof(crossover_coeffs);
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
int32_t crossover_control_initialize(
    void*       params,
    void*       coeffs,
    uint32_t    sample_rate)
{
    crossover_params* init_params;
    crossover_coeffs* init_coeffs;
    uint32_t i;

    init_params = (crossover_params*)params;
    init_coeffs = (crossover_coeffs*)coeffs;

    init_params->cutoff_freq[LOW_CUTOFF] = 0.0;
    init_params->cutoff_freq[MID_CUTOFF] = 0.0;
    init_params->cutoff_freq[HIG_CUTOFF] = 0.0;
    init_params->samplerate = sample_rate;

    for (i = 0; i < (BAND_NUM - 1); i++)
    {
        init_coeffs->band[i].ord_1.k0 = 0.0;
        init_coeffs->band[i].ord_2.k1 = 0.0;
        init_coeffs->band[i].ord_2.k2 = 0.0;

        init_coeffs->band[i].ord_1.negk0 = 0.0;
        init_coeffs->band[i].ord_2.negk1 = 0.0;
        init_coeffs->band[i].ord_2.negk2 = 0.0;
    }
    
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
int32_t crossover_set_parameter(
    void*       params,
    int32_t     id,
    float       value)
{
    crossover_params* set_params;

    set_params = (crossover_params*)params;

    switch (id)
    {
        case 200:
        {
            set_params->cutoff_freq[LOW_CUTOFF] = value;
            break;
        }

        case 201:
        {
            set_params->cutoff_freq[MID_CUTOFF] = value;
            break;
        }

        case 202:
        {
            set_params->cutoff_freq[HIG_CUTOFF] = value;
            break;
        }

        case 210:
        {
            set_params->samplerate = value;
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
int32_t crossover_update_coeffs(
    void const* params,
    void*       coeffs)
{
    crossover_params* init_params;
    crossover_coeffs* init_coeffs;
    uint32_t i;

    float fb;
    float var1;
    float var2;

    init_params = (crossover_params*)params;
    init_coeffs = (crossover_coeffs*)coeffs;

    
    for (i = 0; i < (BAND_NUM - 1); i++)
    {
        fb = init_params->cutoff_freq[i] / (float)Q;
        var1 = (tanf(M_PI * (fb / init_params->samplerate)) - 1.0) / (tanf(M_PI * (fb / init_params->samplerate)) + 1.0);
        var2 = -cosf(2.0 * M_PI * init_params->cutoff_freq[i] / init_params->samplerate);

        init_coeffs->band[i].ord_1.k0 = (tanf(M_PI * (init_params->cutoff_freq[i] / init_params->samplerate)) - 1.0) / (tanf(M_PI * (init_params->cutoff_freq[i] / init_params->samplerate)) + 1.0);
        init_coeffs->band[i].ord_2.k1 = var2 * (1.0 - var1);
        init_coeffs->band[i].ord_2.k2 = -var1;

        init_coeffs->band[i].ord_1.negk0 = -init_coeffs->band[i].ord_1.k0;
        init_coeffs->band[i].ord_2.negk1 = -init_coeffs->band[i].ord_2.k1;
        init_coeffs->band[i].ord_2.negk2 = -init_coeffs->band[i].ord_2.k2;
    }
    
}
