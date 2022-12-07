#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

#define TAP_65 65

typedef struct {
    my_sint32 L;
    my_sint32 R;
}tStereo;

typedef struct {
    tStereo buffer[TAP_65];
}tStates;

typedef struct {
    my_sint32 coef[TAP_65];
}tCoefficient;

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
    *states_bytes = sizeof(tStates);
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
    int i;
    size_t states_size;

    states_size = sizeof(tStates) / sizeof(my_sint32);

    for (i = 0; i < states_size; i++)
    {
        *(((my_sint32*)states)++) = 0;            // *(float*)states) = h[i];   ((float*)states)++;  
    }
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
    int i;
    int j;
    int currentIdx;
    int convIdx;
    
    my_sint32 R_conv;       // convolution value for certain channel in sample
    my_sint32 L_conv;

    my_sint32 result;

    tStates* ptrStates;
    tStereo* ptrAudio;

    ptrStates = (tStates*)states;           // pointer to circulart buffer
    ptrAudio = (tStereo*)audio;             // pointer to certain data sample which is used for updating circular buffer

    for (i = 0; i < samples_count; i++)    // for all data processing
    {
        R_conv = 0;
        L_conv = 0;

        currentIdx = TAP_65 - 1 -(i % TAP_65);
        ptrStates->buffer[currentIdx] = ptrAudio[i];

        for (j = 0; j < TAP_65; j++)        // for one multichannel procession
        {
            convIdx = (currentIdx + j) % TAP_65;

            L_conv = mac32(ptrStates->buffer[convIdx].L, ((my_sint32*)coeffs)[j], L_conv);
            R_conv = mac32(ptrStates->buffer[convIdx].R, ((my_sint32*)coeffs)[j], R_conv);
        }

        result = (my_sint32)L_conv;
        ptrAudio[i].L = result;
        result = ptrAudio[i].L;


        result = (my_sint32)R_conv;
        ptrAudio[i].R = result;
        result = ptrAudio[i].R;
    }
}


#endif
