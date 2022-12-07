#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

#define BANDS 10

typedef struct {
    float L;
    float R;
}tStereo;

typedef struct {
    tStereo x0[BANDS];
    tStereo x1[BANDS];
    tStereo x2[BANDS];
    tStereo y0[BANDS];
    tStereo y1[BANDS];
    tStereo y2[BANDS];
} equalizer_states;

typedef struct coeffs_s {
    float  b0[BANDS];
    float  b1[BANDS];
    float  b2[BANDS];
    float  a0[BANDS];
    float  a1[BANDS];
    float  a2[BANDS];
} equalizer_coeffs;

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
int32_t effect_reset(
    void const* coeffs,
    void*       states)
{
    equalizer_states *reset_states;
    size_t i;

    reset_states = (equalizer_states*)states;
    
    for ( i = 0; i < (size_t)BANDS; i++)
    {
        reset_states->x0[i].L = 0;
        reset_states->x1[i].L = 0;
        reset_states->x2[i].L = 0;
        reset_states->y0[i].L = 0;
        reset_states->y1[i].L = 0;
        reset_states->y2[i].L = 0;

        reset_states->x0[i].R = 0;
        reset_states->x1[i].R = 0;
        reset_states->x2[i].R = 0;
        reset_states->y0[i].R = 0;
        reset_states->y1[i].R = 0;
        reset_states->y2[i].R = 0;
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
int32_t effect_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    equalizer_coeffs *c = (equalizer_coeffs*)coeffs;
    equalizer_states *s = (equalizer_states*)states;
    tStereo *a = (tStereo*)audio;

    for (size_t i = 0; i < samples_count; i++)
    {
        for (size_t j = 0; j < 10; j++)
        {
            if (c->a0[j] != 0)
            {

                s->x0[j].L = a[i].L;
                s->x0[j].R = a[i].R;


                s->y0[j].L = mac_f(c->b0[j], s->x0[j].L, s->x1[j].L) * 8;

                s->x1[j].L = mac_f(c->b1[j], s->x0[j].L, s->x2[j].L);
                s->x1[j].L = msub_f(c->a1[j], s->y0[j].L, s->x1[j].L);

                s->x2[j].L = mul_f(c->b2[j], s->x0[j].L);
                s->x2[j].L = msub_f(c->a2[j], s->y0[j].L, s->x2[j].L);

                a[i].L = s->y0[j].L;


                s->y0[j].R = mac_f(c->b0[j], s->x0[j].R, s->x1[j].R) * 8;

                s->x1[j].R = mac_f(c->b1[j], s->x0[j].R, s->x2[j].R);
                s->x1[j].R = msub_f(c->a1[j], s->y0[j].R, s->x1[j].R);

                s->x2[j].R = mul_f(c->b2[j], s->x0[j].R);
                s->x2[j].R = msub_f(c->a2[j], s->y0[j].R, s->x2[j].R);

                a[i].R = s->y0[j].R;
            }
        }
    }
    return 0;
}


#endif
