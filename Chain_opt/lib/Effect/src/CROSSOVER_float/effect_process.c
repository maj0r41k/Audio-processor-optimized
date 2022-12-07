#ifndef __EFFECT_PROCESS_H__
#define __EFFECT_PROCESS_H__

#include <stdint.h>
#include <stddef.h>
#include "process_fcns.h"
#include "crossover_datatypes.h"


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
int32_t effect_reset(
    void const* coeffs,
    void*       states)
{ 
    crossover_states* reset_states = (crossover_states*)states;
    
    reset_states->cd1_ord1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd1_ord1[1] = (tStereo) { 0.0, 0.0 };

    reset_states->cd1_delay_1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd1_delay_1[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd1_delay_1[2] = (tStereo) { 0.0, 0.0 };

    reset_states->cd1_delay_2[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd1_delay_2[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd1_delay_2[2] = (tStereo) { 0.0, 0.0 };


    reset_states->norm_ord1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->norm_ord1[1] = (tStereo) { 0.0, 0.0 };

    reset_states->norm_delay_1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->norm_delay_1[1] = (tStereo) { 0.0, 0.0 };
    reset_states->norm_delay_1[2] = (tStereo) { 0.0, 0.0 };

    reset_states->norm_delay_2[0] = (tStereo) { 0.0, 0.0 };
    reset_states->norm_delay_2[1] = (tStereo) { 0.0, 0.0 };
    reset_states->norm_delay_2[2] = (tStereo) { 0.0, 0.0 };


    reset_states->cd2_ord1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd2_ord1[1] = (tStereo) { 0.0, 0.0 };

    reset_states->cd2_delay_1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd2_delay_1[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd2_delay_1[2] = (tStereo) { 0.0, 0.0 };

    reset_states->cd2_delay_2[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd2_delay_2[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd2_delay_2[2] = (tStereo) { 0.0, 0.0 };


    reset_states->cd3_ord1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd3_ord1[1] = (tStereo) { 0.0, 0.0 };

    reset_states->cd3_delay_1[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd3_delay_1[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd3_delay_1[2] = (tStereo) { 0.0, 0.0 };

    reset_states->cd3_delay_2[0] = (tStereo) { 0.0, 0.0 };
    reset_states->cd3_delay_2[1] = (tStereo) { 0.0, 0.0 };
    reset_states->cd3_delay_2[2] = (tStereo) { 0.0, 0.0 };

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
    crossover_coeffs* pr_coeffs;
    crossover_states* pr_states;
    tStereo* pr_audio;
    split_band band;

    my_float xh;
    my_float y_ord_1;
    my_float y_ord_2;
    my_float accum_1;       // store data result of 1-st order all-pass
    my_float accum_par;     // store intermidiate result of parallel filter sum

    my_float low_pass;
    my_float high_pass; 


    float compensation_1;
    float compensation_2;
    pr_coeffs = (crossover_coeffs*)coeffs;
    pr_states = (crossover_states*)states;
    pr_audio = (tStereo*)audio;

    for (my_uint32 i = 0; i < samples_count; i++)
    {
        pr_states->in = ((tStereo*)pr_audio)[i];
        

                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //          Parallel calculation 1
        //
        xh = mac_f(pr_coeffs->band[MID_CUTOFF].negk0, pr_states->cd1_ord1[PAR_1].L, pr_states->in.L);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[MID_CUTOFF].k0, pr_states->cd1_ord1[PAR_1].L);
        pr_states->cd1_ord1[PAR_1].L = xh;

        accum_1 = y_ord_1;
        
        // [in] on 1-st 2-nd order in parallel
        y_ord_2 = mac_f(pr_states->in.L, pr_coeffs->band[MID_CUTOFF].k2, pr_states->cd1_delay_1[PAR_1].L);
        pr_states->cd1_delay_1[PAR_1].L = mac_f(pr_states->in.L, pr_coeffs->band[MID_CUTOFF].k1, pr_states->cd1_delay_2[PAR_1].L);
        pr_states->cd1_delay_1[PAR_1].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk1, y_ord_2, pr_states->cd1_delay_1[PAR_1].L);
        pr_states->cd1_delay_2[PAR_1].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk2, y_ord_2, pr_states->in.L);
        
        
        // parallel connected 1 cascade all-passes accumulation
        accum_par = y_ord_1 + y_ord_2;
        accum_par /= 2.0;

                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //          Parallel calculation 1
        //
        xh = mac_f(pr_coeffs->band[MID_CUTOFF].negk0, pr_states->cd1_ord1[PAR_2].L, accum_par);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[MID_CUTOFF].k0, pr_states->cd1_ord1[PAR_2].L);
        pr_states->cd1_ord1[PAR_2].L = xh;

        // [in] on 2-st 2-nd order in parallel
        y_ord_2 = mac_f(accum_par, pr_coeffs->band[MID_CUTOFF].k2, pr_states->cd1_delay_1[PAR_2].L);
        pr_states->cd1_delay_1[PAR_2].L = mac_f(accum_par, pr_coeffs->band[MID_CUTOFF].k1, pr_states->cd1_delay_2[PAR_2].L);
        pr_states->cd1_delay_1[PAR_2].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk1, y_ord_2, pr_states->cd1_delay_1[PAR_2].L);
        pr_states->cd1_delay_2[PAR_2].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk2, y_ord_2, accum_par);

        low_pass = y_ord_1 + y_ord_2;
        low_pass /= 2.0;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //      Serial calculation
        //
        y_ord_2 = mac_f(accum_1, pr_coeffs->band[MID_CUTOFF].k2, pr_states->cd1_delay_1[SERIAL].L);
        pr_states->cd1_delay_1[SERIAL].L = mac_f(accum_1, pr_coeffs->band[MID_CUTOFF].k1, pr_states->cd1_delay_2[SERIAL].L);
        pr_states->cd1_delay_1[SERIAL].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk1, y_ord_2, pr_states->cd1_delay_1[SERIAL].L);
        pr_states->cd1_delay_2[SERIAL].L = mac_f(pr_coeffs->band[MID_CUTOFF].negk2, y_ord_2, accum_1);

        high_pass = y_ord_2 - low_pass;         // ne trogat
         
        
        //
        // low pass splitting on BAND_3 friquency
        //

        xh = mac_f(pr_coeffs->band[HIG_CUTOFF].negk0, pr_states->norm_ord1[NORM_LOW].L, low_pass);
        y_ord_1 = mac_f(xh, pr_coeffs->band[HIG_CUTOFF].k0, pr_states->norm_ord1[NORM_LOW].L);
        pr_states->norm_ord1[NORM_LOW].L = xh;           // verno
        
        y_ord_2 = mac_f(y_ord_1, pr_coeffs->band[HIG_CUTOFF].k2, pr_states->norm_delay_1[NORM_LOW].L);
        pr_states->norm_delay_1[NORM_LOW].L = mac_f(y_ord_1, pr_coeffs->band[HIG_CUTOFF].k1, pr_states->norm_delay_2[NORM_LOW].L);
        pr_states->norm_delay_1[NORM_LOW].L = mac_f(y_ord_2, pr_coeffs->band[HIG_CUTOFF].negk1, pr_states->norm_delay_1[NORM_LOW].L);
        pr_states->norm_delay_2[NORM_LOW].L = mac_f(y_ord_2, pr_coeffs->band[HIG_CUTOFF].negk2, y_ord_1);       //verno

        low_pass = y_ord_2;

                   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //          Cascade 2   low_pass [in]
        //

        xh = mac_f(pr_coeffs->band[LOW_CUTOFF].negk0, pr_states->cd2_ord1[PAR_1].L, low_pass);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[LOW_CUTOFF].k0, pr_states->cd2_ord1[PAR_1].L);
        pr_states->cd2_ord1[PAR_1].L = xh;
        accum_1 = y_ord_1;

        // [in] on 1-st 2-nd order in parallel
        y_ord_2 = mac_f(low_pass, pr_coeffs->band[LOW_CUTOFF].k2, pr_states->cd2_delay_1[PAR_1].L);
        pr_states->cd2_delay_1[PAR_1].L = mac_f(low_pass, pr_coeffs->band[LOW_CUTOFF].k1, pr_states->cd2_delay_2[PAR_1].L);
        pr_states->cd2_delay_1[PAR_1].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk1, y_ord_2, pr_states->cd2_delay_1[PAR_1].L);
        pr_states->cd2_delay_2[PAR_1].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk2, y_ord_2, low_pass);

        accum_par = y_ord_1 + y_ord_2;
        accum_par /= 2.0;

        ///////////////////////////////////

        xh = mac_f(pr_coeffs->band[LOW_CUTOFF].negk0, pr_states->cd2_ord1[PAR_2].L, accum_par);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[LOW_CUTOFF].k0, pr_states->cd2_ord1[PAR_2].L);
        pr_states->cd2_ord1[PAR_2].L = xh;

        y_ord_2 = mac_f(accum_par, pr_coeffs->band[LOW_CUTOFF].k2, pr_states->cd2_delay_1[PAR_2].L);
        pr_states->cd2_delay_1[PAR_2].L = mac_f(accum_par, pr_coeffs->band[LOW_CUTOFF].k1, pr_states->cd2_delay_2[PAR_2].L);
        pr_states->cd2_delay_1[PAR_2].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk1, y_ord_2, pr_states->cd2_delay_1[PAR_2].L);
        pr_states->cd2_delay_2[PAR_2].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk2, y_ord_2, accum_par);

        accum_par = y_ord_1 + y_ord_2;
        accum_par /= 2.0;                               // band 1

        band.band_1.L = accum_par;
        /////////////////////////////////////

        y_ord_2 = mac_f(accum_1, pr_coeffs->band[LOW_CUTOFF].k2, pr_states->cd2_delay_1[SERIAL].L);
        pr_states->cd2_delay_1[PAR_1].L = mac_f(accum_1, pr_coeffs->band[LOW_CUTOFF].k1, pr_states->cd2_delay_2[SERIAL].L);
        pr_states->cd2_delay_1[PAR_1].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk1, y_ord_2, pr_states->cd2_delay_1[PAR_1].L);
        pr_states->cd2_delay_2[PAR_1].L = mac_f(pr_coeffs->band[LOW_CUTOFF].negk2, y_ord_2,accum_1);

        band.band_2.L = y_ord_2 - band.band_1.L;



        //
        // high-pass splitting on BAND_4 friquency [in] = high_pass 
        //

        xh = mac_f(pr_coeffs->band[LOW_CUTOFF].negk0, pr_states->norm_ord1[NORM_HIG].L, high_pass);
        y_ord_1 = mac_f(xh, pr_coeffs->band[LOW_CUTOFF].k0, pr_states->norm_ord1[NORM_HIG].L);
        pr_states->norm_ord1[NORM_HIG].L = xh;           // verno

        y_ord_2 = mac_f(y_ord_1, pr_coeffs->band[LOW_CUTOFF].k2, pr_states->norm_delay_1[NORM_HIG].L);
        pr_states->norm_delay_1[NORM_HIG].L = mac_f(y_ord_1, pr_coeffs->band[LOW_CUTOFF].k1, pr_states->norm_delay_2[NORM_HIG].L);
        pr_states->norm_delay_1[NORM_HIG].L = mac_f(y_ord_2, pr_coeffs->band[LOW_CUTOFF].negk1, pr_states->norm_delay_1[NORM_HIG].L);
        pr_states->norm_delay_2[NORM_HIG].L = mac_f(y_ord_2, pr_coeffs->band[LOW_CUTOFF].negk2, y_ord_1);       //verno

        high_pass = y_ord_2;

        //      **************************
        //      Cascade 3 [in] = high_pass
        //      **************************

        xh = mac_f(pr_coeffs->band[HIG_CUTOFF].negk0, pr_states->cd3_ord1[PAR_1].L, high_pass);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[HIG_CUTOFF].k0, pr_states->cd3_ord1[PAR_1].L);
        pr_states->cd3_ord1[PAR_1].L = xh;
        accum_1 = y_ord_1;

        // [in] on 1-st 2-nd order in parallel
        y_ord_2 = mac_f(high_pass, pr_coeffs->band[HIG_CUTOFF].k2, pr_states->cd3_delay_1[PAR_1].L);
        pr_states->cd3_delay_1[PAR_1].L = mac_f(high_pass, pr_coeffs->band[HIG_CUTOFF].k1, pr_states->cd3_delay_2[PAR_1].L);
        pr_states->cd3_delay_1[PAR_1].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk1, y_ord_2, pr_states->cd3_delay_1[PAR_1].L);
        pr_states->cd3_delay_2[PAR_1].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk2, y_ord_2, high_pass);

        accum_par = y_ord_1 + y_ord_2;
        accum_par /= 2.0;

        ///////////////////////////////////
        xh = mac_f(pr_coeffs->band[HIG_CUTOFF].negk0, pr_states->cd3_ord1[PAR_2].L, accum_par);        // replace by nmac_f
        y_ord_1 = mac_f(xh, pr_coeffs->band[HIG_CUTOFF].k0, pr_states->cd3_ord1[PAR_2].L);
        pr_states->cd3_ord1[PAR_2].L = xh;

        y_ord_2 = mac_f(accum_par, pr_coeffs->band[HIG_CUTOFF].k2, pr_states->cd3_delay_1[PAR_2].L);
        pr_states->cd3_delay_1[PAR_2].L = mac_f(accum_par, pr_coeffs->band[HIG_CUTOFF].k1, pr_states->cd3_delay_2[PAR_2].L);
        pr_states->cd3_delay_1[PAR_2].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk1, y_ord_2, pr_states->cd3_delay_1[PAR_2].L);
        pr_states->cd3_delay_2[PAR_2].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk2, y_ord_2, accum_par);

        accum_par = y_ord_1 + y_ord_2;
        accum_par /= 2.0;

        band.band_3.L = accum_par;

        ///////
        y_ord_2 = mac_f(accum_1, pr_coeffs->band[HIG_CUTOFF].k2, pr_states->cd3_delay_1[SERIAL].L);
        pr_states->cd3_delay_1[SERIAL].L = mac_f(accum_1, pr_coeffs->band[HIG_CUTOFF].k1, pr_states->cd3_delay_2[SERIAL].L);
        pr_states->cd3_delay_1[SERIAL].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk1, y_ord_2, pr_states->cd3_delay_1[SERIAL].L);
        pr_states->cd3_delay_2[SERIAL].L = mac_f(pr_coeffs->band[HIG_CUTOFF].negk2, y_ord_2, accum_1);

        band.band_4.L = y_ord_2 - band.band_3.L;
       
       

        
        ((tStereo*)pr_audio)[i].L = band.band_1.L + band.band_2.L;
        ((tStereo*)pr_audio)[i].R = band.band_1.L + band.band_2.L;
        //((tStereo*)pr_audio)[i].R = band.band_3.L + band.band_4.L;


        

    }

}


#endif
