#include <stdint.h>
#include <stddef.h>

#include "crossover_control.h"
#include "crossover_process.h"
#include "time.h"
#include "stdio.h"

/*******************************************************************************
 * Provides with the required data sizes for states. It is caller responsibility
 *   to allocate enough memory (bytes) for them.
 * 
 * @param[out] states_bytes required data size for storing states
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t crossover_process_get_sizes(
    size_t*     states_bytes)
{
    *states_bytes = sizeof(crossover_states);
}


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
    void*       states)
{ 
    crossover_states* reset_states = (crossover_states*)states;
    
    reset_states->in = _mm_set_ps1(0.0); 

    reset_states->cd1_ord1[0] = _mm_set_ps1(0.0);
    reset_states->cd1_ord1[1] = _mm_set_ps1(0.0);
    reset_states->cd1_ord1[2] = _mm_set_ps1(0.0);

    reset_states->cd1_delay_1[0] = _mm_set_ps1(0.0);
    reset_states->cd1_delay_1[1] = _mm_set_ps1(0.0);

    reset_states->cd1_delay_2[0] = _mm_set_ps1(0.0);
    reset_states->cd1_delay_2[1] = _mm_set_ps1(0.0);



    reset_states->norm_ord1 = _mm_set_ps1(0.0);
    reset_states->norm_delay_1 = _mm_set_ps1(0.0);
    reset_states->norm_delay_2 = _mm_set_ps1(0.0);
    


    reset_states->cd2_ord1[0] = _mm_set_ps1(0.0);
    reset_states->cd2_ord1[1] = _mm_set_ps1(0.0);
    reset_states->cd2_ord1[2] = _mm_set_ps1(0.0);

    reset_states->cd2_delay_1[0] = _mm_set_ps1(0.0);
    reset_states->cd2_delay_1[1] = _mm_set_ps1(0.0);

    reset_states->cd2_delay_2[0] = _mm_set_ps1(0.0);
    reset_states->cd2_delay_2[1] = _mm_set_ps1(0.0);

     

    /*reset_states->cd3_ord1[0] = _mm_set_ps1(0.0);
    reset_states->cd3_ord1[1] = _mm_set_ps1(0.0);
    reset_states->cd3_ord1[2] = _mm_set_ps1(0.0);

    reset_states->cd3_delay_1[0] = _mm_set_ps1(0.0);
    reset_states->cd3_delay_1[1] = _mm_set_ps1(0.0);

    reset_states->cd3_delay_2[0] = _mm_set_ps1(0.0);
    reset_states->cd3_delay_2[1] = _mm_set_ps1(0.0);*/


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
int32_t crossover_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    crossover_coeffs* pr_coeffs;
    crossover_states* pr_states;
    audio_buf_cross* pr_audio;
  
    my_float y_ord_1;
    my_float y_ord_2;
    my_float sum;
    my_float acc_ord_2;

    my_float __low_pass;
    my_float __high_pass;
    my_float __phase_corr;

    pr_coeffs = (crossover_coeffs*)coeffs;
    pr_states = (crossover_states*)states;
    pr_audio = (audio_buf_cross*)audio;
 
    for (my_uint32 i = 0; i < samples_count; i++)
    {
         
        pr_states->in = _mm_setr_ps(((tStereo_cross*)pr_audio->audio)[i].L, ((tStereo_cross*)pr_audio->audio)[i].R, 0.0, 0.0);
        
        y_ord_1 = order_1_single(&pr_states->in, &pr_coeffs->band[MID_CUTOFF].ord_1, &pr_states->cd1_ord1[PAR_1]);
        y_ord_2 = order_2_single(&pr_states->in, &pr_coeffs->band[MID_CUTOFF].ord_2, &pr_states->cd1_delay_1[PAR_1], &pr_states->cd1_delay_2[PAR_1]);
        _mm_store_ps(&acc_ord_2, y_ord_2); // store y_ord_2 by addres &acc_ord_2;

        sum = add_f(y_ord_1, y_ord_2);
        sum = div_f(sum, _mm_setr_ps( 2.0, 2.0, 1.0, 1.0 ));
   
        y_ord_1 = order_1_single(&sum, &pr_coeffs->band[MID_CUTOFF].ord_1, &pr_states->cd1_ord1[PAR_2]);
        y_ord_2 = order_2_single(&sum, &pr_coeffs->band[MID_CUTOFF].ord_2, &pr_states->cd1_delay_1[PAR_2], &pr_states->cd1_delay_2[PAR_2]);

        sum = add_f(y_ord_1, y_ord_2);
        sum = div_f(sum, _mm_setr_ps( 2.0, 2.0, 1.0, 1.0 ));

        _mm_store_ps(&__low_pass, sum);// store y_ord_2 by addres &acc_ord_2;

        y_ord_1 = order_1_single(&acc_ord_2, &pr_coeffs->band[MID_CUTOFF].ord_1, &pr_states->cd1_ord1[SERIAL]);
        __high_pass = sub_f(y_ord_1, __low_pass);  

        __phase_corr = _mm_setr_ps( __low_pass.m128_f32[0], __low_pass.m128_f32[1], __high_pass.m128_f32[0], __high_pass.m128_f32[1] ); // __phase_corr = {LOW.L, LOW.R, HIGH.L, HIGH.R}

        __phase_corr = order_1_dual(&__phase_corr, &pr_coeffs->band[HIG_CUTOFF].ord_1, &pr_coeffs->band[LOW_CUTOFF].ord_1, &pr_states->norm_ord1);     
        __phase_corr = order_2_dual(&__phase_corr, &pr_coeffs->band[HIG_CUTOFF].ord_2, &pr_coeffs->band[LOW_CUTOFF].ord_2, &pr_states->norm_delay_1, &pr_states->norm_delay_2);


        y_ord_1 = order_1_dual(&__phase_corr, &pr_coeffs->band[LOW_CUTOFF].ord_1, &pr_coeffs->band[HIG_CUTOFF].ord_1, &pr_states->cd2_ord1[PAR_1]);
        y_ord_2 = order_2_dual(&__phase_corr, &pr_coeffs->band[LOW_CUTOFF].ord_2, &pr_coeffs->band[HIG_CUTOFF].ord_2, &pr_states->cd2_delay_1[PAR_1], &pr_states->cd2_delay_2[PAR_1]);
        _mm_store_ps(&acc_ord_2, y_ord_2);

        sum = add_f(y_ord_1, y_ord_2);
        sum = div_f(sum, _mm_setr_ps(2.0, 2.0, 2.0, 2.0 ));
        
        y_ord_1 = order_1_dual(&sum, &pr_coeffs->band[LOW_CUTOFF].ord_1, &pr_coeffs->band[HIG_CUTOFF].ord_1, &pr_states->cd2_ord1[PAR_2]);
        y_ord_2 = order_2_dual(&sum, &pr_coeffs->band[LOW_CUTOFF].ord_2, &pr_coeffs->band[HIG_CUTOFF].ord_2, &pr_states->cd2_delay_1[PAR_2], &pr_states->cd2_delay_2[PAR_2]);

        sum = add_f(y_ord_1, y_ord_2);
        sum = div_f(sum, _mm_setr_ps(2.0, 2.0, 2.0, 2.0 ));

        //Band_1 init
        ((tStereo_cross*)pr_audio->cross_b.band_1)[i].L = sum.m128_f32[0];
        ((tStereo_cross*)pr_audio->cross_b.band_1)[i].R = sum.m128_f32[1];

        //Band_3 init
        ((tStereo_cross*)pr_audio->cross_b.band_3)[i].L = sum.m128_f32[2];
        ((tStereo_cross*)pr_audio->cross_b.band_3)[i].R = sum.m128_f32[3];

        y_ord_1 = order_1_dual(&acc_ord_2, &pr_coeffs->band[LOW_CUTOFF].ord_1, &pr_coeffs->band[HIG_CUTOFF].ord_1, &pr_states->cd2_ord1[SERIAL]);
        sum = sub_f(y_ord_1, sum);

        //Band_2 init
        ((tStereo_cross*)pr_audio->cross_b.band_2)[i].L = sum.m128_f32[0];
        ((tStereo_cross*)pr_audio->cross_b.band_2)[i].R = sum.m128_f32[1];

        //Band_4 init
        ((tStereo_cross*)pr_audio->cross_b.band_4)[i].L = sum.m128_f32[2];
        ((tStereo_cross*)pr_audio->cross_b.band_4)[i].R = sum.m128_f32[3];

        //((tStereo_cross*)pr_audio->audio)[i].L = ((tStereo_cross*)pr_audio->cross_b.band_1)[i].L;// +((tStereo_cross*)pr_audio->cross_b.band_2)[i].L + ((tStereo_cross*)pr_audio->cross_b.band_3)[i].L + ((tStereo_cross*)pr_audio->cross_b.band_4)[i].L;
        //((tStereo_cross*)pr_audio->audio)[i].R = ((tStereo_cross*)pr_audio->cross_b.band_1)[i].R;// +((tStereo_cross*)pr_audio->cross_b.band_2)[i].R + ((tStereo_cross*)pr_audio->cross_b.band_3)[i].R + ((tStereo_cross*)pr_audio->cross_b.band_4)[i].R;
    }
    
}

inline my_float order_1_single(my_float *in, coef_1_ord *coeffs, my_float *delay)
{
    my_float xh;
    my_float out;

    xh = mac_f(*delay, _mm_setr_ps(coeffs->negk0, coeffs->negk0, 0.0, 0.0 ), *in);
    out = mac_f(xh, _mm_setr_ps(coeffs->k0, coeffs->k0, 0.0, 0.0 ), *delay);
    *delay = xh;

    return out;    
}

inline my_float order_1_dual(my_float *in, coef_1_ord *coeffs_1, coef_1_ord *coeffs_2, my_float *delay)
{
    my_float xh;
    my_float out;

    xh = mac_f(*delay, _mm_setr_ps(coeffs_1->negk0, coeffs_1->negk0, coeffs_2->negk0, coeffs_2->negk0 ), *in);
    out = mac_f(xh, _mm_setr_ps(coeffs_1->k0, coeffs_1->k0, coeffs_2->k0, coeffs_2->k0 ), *delay);
    *delay = xh;

    return out;
}

inline my_float order_2_single(my_float *in, coef_2_ord *coeffs, my_float *delay1, my_float *delay2)
{
    my_float out;

    out = mac_f(*in, _mm_setr_ps(coeffs->k2, coeffs->k2, 0.0, 0.0), *delay1);
    *delay1 = mac_f(*in, _mm_setr_ps(coeffs->k1, coeffs->k1, 0.0, 0.0 ), *delay2);
    *delay1 = mac_f(out, _mm_setr_ps(coeffs->negk1, coeffs->negk1, 0.0, 0.0 ), *delay1);
    *delay2 = mac_f(out, _mm_setr_ps(coeffs->negk2, coeffs->negk2, 0.0, 0.0 ), *in);

    return out;   
}


inline my_float order_2_dual(my_float *in, coef_2_ord *coeffs_1, coef_2_ord *coeffs_2, my_float *delay1, my_float *delay2)
{
    my_float out;

    out = mac_f(*in, _mm_setr_ps(coeffs_1->k2, coeffs_1->k2, coeffs_2->k2, coeffs_2->k2 ), *delay1);
    *delay1 = mac_f(*in, _mm_setr_ps(coeffs_1->k1, coeffs_1->k1, coeffs_2->k1, coeffs_2->k1 ), *delay2);
    *delay1 = mac_f(out, _mm_setr_ps(coeffs_1->negk1, coeffs_1->negk1, coeffs_2->negk1, coeffs_2->negk1 ), *delay1);
    *delay2 = mac_f(out, _mm_setr_ps(coeffs_1->negk2, coeffs_1->negk2, coeffs_2->negk2, coeffs_2->negk2 ), *in);

    return out;
}
