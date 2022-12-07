#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>

#include "equalizer_process.h"
#include "equalizer_control.h"
#include "fixedpoint.h"


/*******************************************************************************
 * Provides with the required data sizes for states. It is caller responsibility
 *   to allocate enough memory (bytes) for them.
 * 
 * @param[out] states_bytes required data size for storing states
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t equalizer_process_get_sizes(
    size_t*     states_bytes)
{
    *states_bytes = sizeof(equalizer_states);
    return 0;
}


/*******************************************************************************
 * Reset internal states. Configuration remains the same.
 * 
 * @param[in] coeffs        initialized coeffs
 * @param[in] states        initialized states
 * 
 * @return 0 on success, otherwise fail
 ******************************************************************************/
int32_t equalizer_reset(
    void const* coeffs,
    void*       states)
{
    equalizer_states *reset_states;
    size_t i;

    reset_states = (equalizer_states*)states;
    
    for ( i = 0; i < (size_t)BANDS; i++)
    {
        reset_states->x0[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);
        reset_states->x1[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);
        reset_states->x2[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);
        reset_states->y0[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);
        reset_states->y1[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);
        reset_states->y2[i] = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);       
    }

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
int32_t equalizer_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    equalizer_coeffs *c = (equalizer_coeffs*)coeffs;
    equalizer_states *s = (equalizer_states*)states;
    tStereo_eq *a = (tStereo_eq*)audio;

    for (size_t i = 0; i < samples_count; i++)
    {
        for (size_t j = 0; j < 10; j++)
        {
            if (c->a0[j] != 0)
            {
                s->x0[j] = _mm_setr_ps(a[i].L, a[i].R, 0.0, 0.0);   // vector {L; R; 0; 0}

                s->y0[j] = mac_f(s->x0[j], _mm_setr_ps( c->b0[j], c->b0[j], 0.0, 0.0 ), s->x1[j]);

                s->x1[j] = mac_f(s->x0[j], _mm_setr_ps(c->b1[j], c->b1[j], 0.0, 0.0 ), s->x2[j]);
                s->x1[j] = msub_f(s->y0[j], _mm_setr_ps(c->a1[j], c->a1[j], 0.0, 0.0 ), s->x1[j]);

                s->x2[j] = mul_f(s->x0[j], _mm_setr_ps(c->b2[j], c->b2[j], 0.0, 0.0 ));
                s->x2[j] = msub_f(s->y0[j], _mm_setr_ps(c->a2[j], c->a2[j], 0.0, 0.0 ), s->x2[j]);

                a[i].L = s->y0[j].m128_f32[0];
                a[i].R = s->y0[j].m128_f32[1];

            }
        }
    }
    return 0;
}


#endif
