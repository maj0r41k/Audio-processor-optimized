#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

typedef struct {
    float threshold;
    float ratio;

    float alphaAttack;
    float alphaRelease;

    float attackEnv;
    float releaseEnv;

    float makeUpGain;
    float samplerate;
}compressor_coeffs;

typedef struct {
    float L;
    float R;
}tStereo;

typedef struct {
    tStereo x;             //
    tStereo y;

    tStereo x_dB;        //
    tStereo x_sc;       //

    tStereo g_c;

    tStereo g_s;
    tStereo g_sPrev;

    tStereo g_m;

    tStereo g_lin;

    tStereo envelope;
    tStereo envelope_prev;
}compressor_states;



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
    *states_bytes = sizeof(compressor_states);
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
    compressor_states* reset_states = (compressor_states*)states;

    reset_states->x.L = 0.0;
    reset_states->y.L = 0.0;
    reset_states->x_dB.L = 0.0;
    reset_states->x_sc.L = 0.0;
    reset_states->g_c.L = 0.0;
    reset_states->g_s.L = 0.0;
    reset_states->g_sPrev.L = 1.0;
    reset_states->g_m.L = 0.0;
    reset_states->g_lin.L = 0.0;
    reset_states->envelope.L = 0.0;
    reset_states->envelope_prev.L = 0.0;

    reset_states->x.R = 0.0;
    reset_states->y.R = 0.0;
    reset_states->x_dB.R = 0.0;
    reset_states->x_sc.R = 0.0;
    reset_states->g_c.R = 0.0;
    reset_states->g_s.R = 0.0;
    reset_states->g_sPrev.R = 1.0;
    reset_states->g_m.R = 0.0;
    reset_states->g_lin.R = 0.0;
    reset_states->envelope.R = 0.0;
    reset_states->envelope_prev.R = 0.0;


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
    compressor_coeffs* coeffs_c;
    compressor_states* states_c;
    tStereo* audio_c;
    uint32_t i;

    coeffs_c = (compressor_coeffs*)coeffs;
    states_c = (compressor_states*)states;
    audio_c = (tStereo*)audio;
    float xL_abs;
    float xR_abs;
    float axil1;

    for (i = 0; i < samples_count; i++)
    {
        states_c->x = ((tStereo*)audio_c)[i];
        // envelope of the signal

        xL_abs = fabsf(states_c->x.L);

        if (xL_abs > states_c->envelope_prev.L)              // comparison of current gain and previos gain
        {

            states_c->envelope.L = coeffs_c->attackEnv * states_c->envelope_prev.L + (1.0 - coeffs_c->attackEnv) * xL_abs;     // if current gain higher than previous -> attac
        }
        else
        {
            states_c->envelope.L = coeffs_c->releaseEnv * states_c->envelope_prev.L + (1.0 - coeffs_c->releaseEnv) * xL_abs;    // attenuate
        }

        states_c->envelope_prev.L = states_c->envelope.L;

        if (states_c->envelope.L < coeffs_c->threshold)
        {
            states_c->g_c.L = powf((states_c->envelope.L / coeffs_c->threshold), (coeffs_c->ratio - 1));
        }
        else
        {
            states_c->g_c.L = 1;
        }

        if (states_c->g_c.L <= states_c->g_sPrev.L)
        {
            states_c->g_s.L = coeffs_c->alphaAttack* states_c->g_sPrev.L + (1.0 - coeffs_c->alphaAttack)*states_c->g_c.L;
        }
        else
        {
            states_c->g_s.L = coeffs_c->alphaRelease* states_c->g_sPrev.L + (1.0 - coeffs_c->alphaRelease)*states_c->g_c.L;
        }

        states_c->g_sPrev.L = states_c->g_s.L;
        states_c->g_m.L = states_c->g_s.L * coeffs_c->makeUpGain;



        xR_abs = fabsf(states_c->x.R);

        if (xR_abs > states_c->envelope_prev.R)              // comparison of current gain and previos gain
        {

            states_c->envelope.R = coeffs_c->attackEnv * states_c->envelope_prev.R + (1.0 - coeffs_c->attackEnv) * xR_abs;     // if current gain higher than previous -> attac
        }
        else
        {
            states_c->envelope.R = coeffs_c->releaseEnv * states_c->envelope_prev.R + (1.0 - coeffs_c->releaseEnv) * xR_abs;    // attenuate
        }

        states_c->envelope_prev.R = states_c->envelope.R;

        if (states_c->envelope.R < coeffs_c->threshold)
        {
            states_c->g_c.R = powf((states_c->envelope.R / coeffs_c->threshold), (coeffs_c->ratio - 1));
        }
        else
        {
            states_c->g_c.R = 1;
        }

        if (states_c->g_c.R <= states_c->g_sPrev.R)
        {
            states_c->g_s.R = coeffs_c->alphaAttack* states_c->g_sPrev.R + (1.0 - coeffs_c->alphaAttack)*states_c->g_c.R;
        }
        else
        {
            states_c->g_s.R = coeffs_c->alphaRelease* states_c->g_sPrev.R + (1.0 - coeffs_c->alphaRelease)*states_c->g_c.R;
        }

        states_c->g_sPrev.R = states_c->g_s.R;
        states_c->g_m.R = states_c->g_s.R * coeffs_c->makeUpGain;

        ((tStereo*)audio_c)[i].L = states_c->x.L;
        ((tStereo*)audio_c)[i].R = states_c->x.L * states_c->g_m.L;

        //((tStereo*)audio_c)[i].L = states_c->x.L * states_c->g_m.L;
        //((tStereo*)audio_c)[i].R = states_c->x.R * states_c->g_m.R;
    }
}


#endif
