#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

typedef struct {
    my_sint32 a0;
    my_sint32 a1;
    my_sint32 a2;

    my_sint32 b0;
    my_sint32 b1;
    my_sint32 b2;
}tIIR_coeffs;

typedef struct {
    my_sint32 R;
    my_sint32 L;
}tStereo;

typedef struct {
    tStereo input;
    tStereo prev_input_1;
    tStereo prev_input_2;

    tStereo output;
    tStereo prev_output_1;
    tStereo prev_output_2;

    tStereo Q_noise;
}tIIR_states;

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
    *states_bytes = sizeof(tIIR_states);
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
    tIIR_states* reset_states;

    reset_states = (tIIR_states*)states;

    reset_states->input.L = 0;
    reset_states->input.R = 0;
    reset_states->prev_input_1.L = 0;
    reset_states->prev_input_1.R = 0;
    reset_states->prev_input_2.L = 0;
    reset_states->prev_input_2.R = 0;

    reset_states->output.L = 0;
    reset_states->output.R = 0;
    reset_states->prev_output_1.L = 0;
    reset_states->prev_output_1.R = 0;
    reset_states->prev_output_2.L = 0;
    reset_states->prev_output_2.R = 0;

    reset_states->Q_noise.L = 0;
    reset_states->Q_noise.R = 0;
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
    tIIR_coeffs* proc_coeffs;
    tIIR_states* proc_states;
    tStereo* proc_audio;
    uint32_t i;

    my_sint64 iir_accum_L;
    my_sint64 iir_accum_R;

    proc_coeffs = (tIIR_coeffs*)coeffs;
    proc_states = (tIIR_states*)states;
    proc_audio = (tStereo*)audio;

    for (i = 0; i < samples_count; i++)
    {
        proc_states->input = ((tStereo*)proc_audio)[i];

        iir_accum_L = 0;
        iir_accum_R = 0;

        iir_accum_L = mac64(proc_coeffs->b0, proc_states->input.L, iir_accum_L);
        iir_accum_L = mac64(proc_coeffs->b1, proc_states->prev_input_1.L, iir_accum_L);
        iir_accum_L = mac64(proc_coeffs->b2, proc_states->prev_input_2.L, iir_accum_L);
        iir_accum_L = msub64(proc_coeffs->a1, proc_states->prev_output_1.L, iir_accum_L);
        iir_accum_L = msub64(proc_coeffs->a2, proc_states->prev_output_2.L, iir_accum_L);
        iir_accum_L = add64(iir_accum_L, proc_states->Q_noise.L);

        iir_accum_R = mac64(proc_coeffs->b0, proc_states->input.R, iir_accum_R);
        iir_accum_R = mac64(proc_coeffs->b1, proc_states->prev_input_1.R, iir_accum_R);
        iir_accum_R = mac64(proc_coeffs->b2, proc_states->prev_input_2.R, iir_accum_R);
        iir_accum_R = msub64(proc_coeffs->a1, proc_states->prev_output_1.R, iir_accum_R);
        iir_accum_R = msub64(proc_coeffs->a2, proc_states->prev_output_2.R, iir_accum_R);
        iir_accum_R = add64(iir_accum_R, proc_states->Q_noise.R);

        proc_states->Q_noise.L = (my_sint32)iir_accum_L;
        proc_states->Q_noise.R = (my_sint32)iir_accum_R;

        iir_accum_L = lsh64(iir_accum_L, 1);
        iir_accum_L = iir_accum_L >> 32;
        proc_states->output.L = (my_sint32)iir_accum_L;

        iir_accum_R = lsh64(iir_accum_R, 1);
        iir_accum_R = iir_accum_R >> 32;
        proc_states->output.R = (my_sint32)iir_accum_R;

        proc_states->prev_input_2.L = proc_states->prev_input_1.L;
        proc_states->prev_input_1.L = proc_states->input.L;

        proc_states->prev_input_2.R = proc_states->prev_input_1.R;
        proc_states->prev_input_1.R = proc_states->input.R;

        proc_states->prev_output_2.L = proc_states->prev_output_1.L;
        proc_states->prev_output_1.L = proc_states->output.L;

        proc_states->prev_output_2.R = proc_states->prev_output_1.R;
        proc_states->prev_output_1.R = proc_states->output.R;

        ((tStereo*)proc_audio)[i] = proc_states->output;
    } 
}


#endif
