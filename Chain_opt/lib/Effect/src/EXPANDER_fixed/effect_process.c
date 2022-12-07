#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

/*#define X_MIN_VAL_THRESHOLD    0x800
#define MIN_DB                 0xC3CB43C0 
#define TWENTY_Q23             0xA000000
#define ONE_Q23                0x800000*/

#define ONE_Q31         0x7FFFFFFF


typedef struct {
    my_sint32 threshold;
    my_sint32 ratio;

    my_sint32 alphaAttack;
    my_sint32 alphaRelease;

    my_sint32 attackEnv;
    my_sint32 releaseEnv;

    my_sint32 makeUpGain;
    my_sint32 samplerate;
}compressor_coeffs;

typedef struct {
    my_sint32 L;
    my_sint32 R;
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

    reset_states->x.L = 0;
    reset_states->y.L = 0;
    reset_states->x_dB.L = 0;
    reset_states->x_sc.L = 0;
    reset_states->g_c.L = 0;
    reset_states->g_s.L = 0;
    reset_states->g_sPrev.L = 0x7FFFFFFF;
    reset_states->g_m.L = 0;
    reset_states->g_lin.L = 0;
    reset_states->envelope.L = 0;
    reset_states->envelope_prev.L = 0;

    reset_states->x.R = 0;
    reset_states->y.R = 0;
    reset_states->x_dB.R = 0;
    reset_states->x_sc.R = 0;
    reset_states->g_c.R = 0;
    reset_states->g_s.R = 0;
    reset_states->g_sPrev.R = 0x7FFFFFFF;
    reset_states->g_m.R = 0;
    reset_states->g_lin.R = 0;
    reset_states->envelope.R = 0;
    reset_states->envelope_prev.R = 0;
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

    my_sint32 xL_abs;
    my_sint32 xR_abs;

    my_sint32 axilL1;
    my_sint32 axilL2;

    my_sint32 axilR1;
    my_sint32 axilR2;

    for (i = 0; i < samples_count; i++)
    {
        states_c->x = ((tStereo*)audio_c)[i];

                            /*      Left channel   */   
        xL_abs = abs32(states_c->x.L);

        if (xL_abs > states_c->envelope_prev.L)
        {
            axilL1 = mul32(coeffs_c->attackEnv, states_c->envelope_prev.L);
            axilL2 = sub32(ONE_Q31, coeffs_c->attackEnv);
            axilL2 = mul32(axilL2, xL_abs);

            states_c->envelope.L = add32(axilL1, axilL2);
        }
        else
        {
            axilL1 = mul32(coeffs_c->releaseEnv, states_c->envelope_prev.L);
            axilL2 = sub32(ONE_Q31, coeffs_c->releaseEnv);
            axilL2 = mul32(axilL2, xL_abs);

            states_c->envelope.L = add32(axilL1, axilL2);
        }

        states_c->envelope_prev.L = states_c->envelope.L;

        if (states_c->envelope.L < coeffs_c->threshold)
        {
            axilL1 = div32(states_c->envelope_prev.L, coeffs_c->threshold );
            axilL2 = 0x7FFFFFFF;

            if (coeffs_c->ratio  == 0)
            {
                states_c->g_c.L = ONE_Q31;
            }
            else
            {
                for (int j = 0; j < coeffs_c->ratio; j++)
                {
                    axilL2 = mul32(axilL1, axilL2);
                }

                states_c->g_c.L = axilL2;
            }
            
        }
        else
        {
            states_c->g_c.L = ONE_Q31;
        }


        if (states_c->g_c.L <= states_c->g_sPrev.L)
        {
            axilL1 = mul32(coeffs_c->alphaAttack, states_c->g_sPrev.L);
            axilL2 = sub32(ONE_Q31, coeffs_c->alphaAttack);
            axilL2 = mul32(axilL2, states_c->g_c.L);

            states_c->g_s.L = add32(axilL1, axilL2);
        }
        else
        {
            axilL1 = mul32(coeffs_c->alphaRelease, states_c->g_sPrev.L);
            axilL2 = sub32(ONE_Q31, coeffs_c->alphaRelease);
            axilL2 = mul32(axilL2, states_c->g_c.L);

            states_c->g_s.L = add32(axilL1, axilL2);
        }

        states_c->g_sPrev.L = states_c->g_s.L;
        states_c->g_m.L = mul32_q(states_c->g_sPrev.L, coeffs_c->makeUpGain, Q27);
        ((tStereo*)audio_c)[i].L = mul32(states_c->x.L, states_c->g_m.L);


                                        /*      Right channel       */    
        xR_abs = abs32(states_c->x.R);

        if (xR_abs > states_c->envelope_prev.R)
        {
            axilR1 = mul32(coeffs_c->attackEnv, states_c->envelope_prev.R);
            axilR2 = sub32(ONE_Q31, coeffs_c->attackEnv);
            axilR2 = mul32(axilR2, xR_abs);

            states_c->envelope.R = add32(axilR1, axilR2);
        }
        else
        {
            axilR1 = mul32(coeffs_c->releaseEnv, states_c->envelope_prev.R);
            axilR2 = sub32(ONE_Q31, coeffs_c->releaseEnv);
            axilR2 = mul32(axilR2, xR_abs);

            states_c->envelope.R = add32(axilR1, axilR2);
        }

        states_c->envelope_prev.R = states_c->envelope.R;

        if (states_c->envelope.R < coeffs_c->threshold)
        {
            axilR1 = div32(states_c->envelope_prev.R, coeffs_c->threshold);
            axilR2 = 0x7FFFFFFF;

            for (int j = 0; j < coeffs_c->ratio; j++)
            {
                axilR2 = mul32(axilR1, axilR2);
            }

            states_c->g_c.R = axilR2;
        }
        else
        {
            states_c->g_c.R = ONE_Q31;
        }


        if (states_c->g_c.R <= states_c->g_sPrev.R)
        {
            axilR1 = mul32(coeffs_c->alphaAttack, states_c->g_sPrev.R);
            axilR2 = sub32(ONE_Q31, coeffs_c->alphaAttack);
            axilR2 = mul32(axilR2, states_c->g_c.R);

            states_c->g_s.R = add32(axilR1, axilR2);
        }
        else
        {
            axilR1 = mul32(coeffs_c->alphaRelease, states_c->g_sPrev.R);
            axilR2 = sub32(ONE_Q31, coeffs_c->alphaRelease);
            axilR2 = mul32(axilR2, states_c->g_c.R);

            states_c->g_s.R = add32(axilR1, axilR2);
        }

        states_c->g_sPrev.R = states_c->g_s.R;
        states_c->g_m.R = mul32_q(states_c->g_sPrev.R, coeffs_c->makeUpGain, Q27);
        ((tStereo*)audio_c)[i].R = mul32(states_c->x.R, states_c->g_m.R);


        
    }
}


#endif
