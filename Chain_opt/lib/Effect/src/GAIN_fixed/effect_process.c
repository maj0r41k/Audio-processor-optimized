#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

typedef struct {
    my_sint32 L;
    my_sint32 R;
}tStereo;

/*******************************************************************************
 * Provides with the required data sizes for states. It is caller responsibility
 *   to allocate enough memory (bytes) for them.
 * 
 * @param[out] states_bytes required data size for storing states
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t effect_process_get_sizes(
    size_t*     states_bytes)
{
    *states_bytes = 0;
}


/*******************************************************************************
 * Reset internal states. Configuration remains the same.
 * 
 * @param[in] coeffs        initialized coeffs
 * @param[in] states        initialized states
 * 
 * @return 0 on success, otherwise fail
 ******************************************************************************/
int32_t effect_reset(
    void const* coeffs,
    void*       states)
{
    return 0;
}


/*******************************************************************************
 * Process all available data in the input audio buffer (in-place processing).
 *   Samples are presented in the interleaved manner [L, R, L, R, ...].
 * 
 * @param[in] coeffs        initialized coeffs
 * @param[in] states        initialized states
 * @param[in,out] audio     array of input samples in the interleaved manner
 * @param[in] samples_count count of multichannel samples (L, R pairs) in audio
 * 
 * @return 0 on success, otherwise fail
 ******************************************************************************/
int32_t effect_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    my_sint32* gain_coeffs = (float*)coeffs;
    tStereo* gain_states = (tStereo*)states;
    tStereo* gain_audio = (tStereo*)audio;

    uint32_t i;
    my_sint32 L_tmp;
    my_sint32 R_tmp;

    for (i = 0; i < samples_count; i++)
    {
        L_tmp = mul32(((tStereo*)gain_audio)[i].L, *gain_coeffs);
        R_tmp = mul32(((tStereo*)gain_audio)[i].R, *gain_coeffs);
        
        if (L_tmp & 0x80000000) //if negative
        {
            if (((L_tmp & 0x70000000) ^ 0x70000000) != 0)
            {
                L_tmp = 0x80000000;
            }
            else
            {
                L_tmp <<= 3;
            }
        }
        else  //if positive
        {
            if ((L_tmp & 0x70000000) != 0)
            {
                L_tmp = 0x7FFFFFFF;
            }
            else
            {
                L_tmp <<= 3;
            }
        }

        if (R_tmp & 0x80000000) //if negative
        {
            if (((R_tmp & 0x70000000) ^ 0x70000000) != 0)
            {
                R_tmp = 0x80000000;
            }
            else
            {
                R_tmp <<= 3;
            }
        }
        else  //if positive
        {
            if ((R_tmp & 0x70000000) != 0)
            {
                R_tmp = 0x7FFFFFFF;
            }
            else
            {
                R_tmp <<= 3;
            }
        }

        ((tStereo*)gain_audio)[i].L = L_tmp;
        ((tStereo*)gain_audio)[i].R = R_tmp;
    }
}


#endif
