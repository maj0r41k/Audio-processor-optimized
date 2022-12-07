#ifndef __EFFECT_CONTROL_H__
#define __EFFECT_CONTROL_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

#define TAP_65 65

//extern float h_fxd[TAP_65];
int h_fxd[65] = { -2424691, 3325511, -371827, -3931498, 4568168, 450345,
                            -6182322, 5407067, 1715256, -6521107, 3494636, 1307894,
                            -643, -2317085, -5721822, 17967180, -9642269, -25402224,
                            47935224, -12355763, -61051912, 84172456, -3337400, -109622720,
                            115801768, 21078928, -160792432, 131305712, 57666468, -200080928,
                            124153440, 96339144, 1933618560, 96339144, 124153440, -200080928,
                            57666468, 131305712, -160792432, 21078928, 115801768, -109622720,
                            -3337400, 84172456, -61051912, -12355763, 47935224, -25402224,
                            -9642269, 17967180, -5721822, -2317085, -643, 1307894,
                            3494636, -6521107, 1715256, 5407067, -6182322, 450345,
                            4568168, -3931498, -371827, 3325511, -2424691 };

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
    *params_bytes = 0;
    *coeffs_bytes = sizeof(h_fxd);
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
    int i;
    size_t coef_size;
    int* copy;
    int var;
    copy = h_fxd;

    coef_size = sizeof(h_fxd)/sizeof(my_sint32);

    for (i = 0; i < coef_size; i++)
    {
        var = h_fxd[i];
        *((my_sint32*)coeffs) = h_fxd[i];
        var = *((my_sint32*)coeffs);
        ((my_sint32*)coeffs)++;
        //*(((my_sint32*)coeffs)++) = h_fxd[i];           //     *(float*)coeffs) = h[i];   ((float*)coeffs)++;  
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
int32_t effect_set_parameter(
    void*       params,
    int32_t     id,
    float       value)
{
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
    return 0;
}


#endif
