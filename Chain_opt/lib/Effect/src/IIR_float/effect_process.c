#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>

typedef struct {
    float a0;
    float a1;
    float a2;

    float b0;
    float b1;
    float b2;
}tIIR_coeffs;

typedef struct {
    float R;
    float L;
}tStereo;

typedef struct {
    tStereo input;
    tStereo prev_input_1;
    tStereo prev_input_2;

    tStereo output;
    tStereo prev_output_1;
    tStereo prev_output_2;
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

    reset_states->input.L = 0.0;
    reset_states->input.R = 0.0;
    reset_states->prev_input_1.L = 0.0;
    reset_states->prev_input_1.R = 0.0;
    reset_states->prev_input_2.L = 0.0;
    reset_states->prev_input_2.R = 0.0;

    reset_states->output.L = 0.0;
    reset_states->output.R = 0.0;
    reset_states->prev_output_1.L = 0.0;
    reset_states->prev_output_1.R = 0.0;
    reset_states->prev_output_2.L = 0.0;
    reset_states->prev_output_2.R = 0.0;
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

    proc_coeffs = (tIIR_coeffs*)coeffs;
    proc_states = (tIIR_states*)states;
    proc_audio = (tStereo*)audio;

    for (i = 0; i < samples_count; i++)
    {
        proc_states->input = ((tStereo*)proc_audio)[i];

        proc_states->output.L = (proc_states->input.L * proc_coeffs->b0) +
            (proc_states->prev_input_1.L * proc_coeffs->b1) +
            (proc_states->prev_input_2.L * proc_coeffs->b2) -
            (proc_states->prev_output_1.L * proc_coeffs->a1) -
            (proc_states->prev_output_2.L * proc_coeffs->a2);

        proc_states->output.R = (proc_states->input.R * proc_coeffs->b0) +
            (proc_states->prev_input_1.R * proc_coeffs->b1) +
            (proc_states->prev_input_2.R * proc_coeffs->b2) -
            (proc_states->prev_output_1.R * proc_coeffs->a1) -
            (proc_states->prev_output_2.R * proc_coeffs->a2);

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
