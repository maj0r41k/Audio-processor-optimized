#ifndef __CROSSOVER_PROCESS_H__
#define __CROSSOVER_PROCESS_H__

#include <stdint.h>
#include <stddef.h>


/*******************************************************************************
 * Provides with the required data sizes for states. It is caller responsibility
 *   to allocate enough memory (bytes) for them.
 * 
 * @param[out] states_bytes required data size for storing states
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t crossover_process_get_sizes(
    size_t*     states_bytes);


/*******************************************************************************
 * Reset internal states. Configuration remains the same.
 * 
 * @param[in] coeffs        initialized coeffs
 * @param[in] states        initialized states
 * 
 * @return 0 on success, otherwise fail
 ******************************************************************************/
int32_t crossover_reset(
    void const* coeffs,
    void*       states);


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
int32_t crossover_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count);


extern my_float order_1_single(my_float *in, coef_1_ord *coeffs, my_float *delay);
extern my_float order_2_single(my_float *in, coef_2_ord *coeffs, my_float *delay1, my_float *delay2);

extern my_float order_1_dual(my_float *in, coef_1_ord *coeffs_1, my_float *delay);
extern my_float order_2_dual(my_float *in, coef_2_ord *coeffs_1, coef_2_ord *coeffs_2, my_float *delay1, my_float *delay2);
/*tStereo_cross order_1(tStereo_cross *in, coef_1_ord *coeffs, tStereo_cross *delay);
tStereo_cross order_2(tStereo_cross *in, coef_2_ord *coeffs, tStereo_cross *delay1, tStereo_cross *delay2);*/

#endif