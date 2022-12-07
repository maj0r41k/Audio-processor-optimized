 #include <stdint.h>
#include <stddef.h>
#include "fixedpoint.h"

#include "compressor_4ch_control.h"
#include "compressor_4ch_process.h"

/*******************************************************************************
 * Provides with the required data sizes for states. It is caller responsibility
 *   to allocate enough memory (bytes) for them.
 * 
 * @param[out] states_bytes required data size for storing states
 * 
 * @return 0 if success, non-zero error code otherwise
 ******************************************************************************/
int32_t compressor_4ch_process_get_sizes(
    size_t*     states_bytes)
{
    *states_bytes = sizeof(compressor_4ch_states);
}


/*******************************************************************************
 * Reset internal states. Configuration remains the same.
 * 
 * @param[in] coeffs        initialized coeffs
 * @param[in] states        initialized states
 * 
 * @return 0 on success, otherwise fail
 ******************************************************************************/
int32_t compressor_4ch_reset(
    void const* coeffs,
    void*       states)
{
    compressor_4ch_states* reset_states = (compressor_4ch_states*)states;

    reset_states->xL =              _mm_set_ps1(0.0);
    reset_states->xR =              _mm_set_ps1(0.0);
    reset_states->g_c =             _mm_set_ps1(0.0);
    reset_states->g_s =             _mm_set_ps1(0.0);
    reset_states->g_sPrev =         _mm_set_ps1(1.0);
    reset_states->g_m =             _mm_set_ps1(0.0);
    reset_states->envelope =        _mm_set_ps1(0.0);
    reset_states->envelope_prev =   _mm_set_ps1(0.0);

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
int32_t compressor_4ch_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    compressor_4ch_coeffs* coeffs_c;
    compressor_4ch_states* states_c;
    audio_buf_comp* audio_c;

    uint32_t i;
    my_float max_abs;
    my_float abs_xL;
    my_float abs_xR;

    my_float _sign_drop_mask;
    my_float _cond_mask;
    my_float _att_mask;
    my_float _rel_mask;

    coeffs_c = (compressor_4ch_coeffs*)coeffs;
    states_c = (compressor_4ch_coeffs*)states;
    audio_c = (audio_buf_comp*)audio;

    _sign_drop_mask.m128_i32[0] = 0x7FFFFFFF;
    _sign_drop_mask.m128_i32[1] = 0x7FFFFFFF;
    _sign_drop_mask.m128_i32[2] = 0x7FFFFFFF;
    _sign_drop_mask.m128_i32[3] = 0x7FFFFFFF;

    _cond_mask.m128_i32[0] = 0xFFFFFFFF;
    _cond_mask.m128_i32[1] = 0xFFFFFFFF;
    _cond_mask.m128_i32[2] = 0xFFFFFFFF;
    _cond_mask.m128_i32[3] = 0xFFFFFFFF;
    

    for (i = 0; i < samples_count; i++)
    {
        abs_xL = _mm_setr_ps(((tStereo_compr*)audio_c->cross_b.band_1)[i].L, ((tStereo_compr*)audio_c->cross_b.band_2)[i].L, ((tStereo_compr*)audio_c->cross_b.band_3)[i].L, ((tStereo_compr*)audio_c->cross_b.band_4)[i].L); // store abs_vals???
        abs_xR = _mm_setr_ps(((tStereo_compr*)audio_c->cross_b.band_1)[i].R, ((tStereo_compr*)audio_c->cross_b.band_2)[i].R, ((tStereo_compr*)audio_c->cross_b.band_3)[i].R, ((tStereo_compr*)audio_c->cross_b.band_4)[i].R);

        abs_xL = _mm_and_ps(abs_xL, _sign_drop_mask); 
        abs_xR = _mm_and_ps(abs_xR, _sign_drop_mask);

        max_abs = _mm_max_ps(abs_xL, abs_xR);

        //  envelope
        _att_mask = _mm_cmpgt_ps(max_abs, states_c->envelope_prev);
        _rel_mask = _mm_andnot_ps(_att_mask, _cond_mask);

        _att_mask = _mm_and_ps(max_abs, _att_mask);
        _rel_mask = _mm_and_ps(max_abs, _rel_mask);

        // Envelope attack
        _att_mask = mul_f(_att_mask, coeffs_c->_one_env_att);
        _att_mask = mac_f(coeffs_c->attackEnv, states_c->envelope_prev, _att_mask);

        // Envelope release
        _rel_mask = mul_f(_rel_mask, coeffs_c->_one_env_rel);
        _rel_mask = mac_f(coeffs_c->releaseEnv, states_c->envelope_prev, _rel_mask);

        states_c->envelope_prev = add_f(_att_mask, _rel_mask);
        states_c->envelope = states_c->envelope_prev;

        ///////////////////////////////////////////////////////////

        // gain computer
        _att_mask = _mm_cmplt_ps(states_c->envelope, coeffs_c->threshold);
        _rel_mask = _mm_andnot_ps(_att_mask, _cond_mask);

        _att_mask = _mm_and_ps(_mm_set_ps1(1.0), _att_mask);
        _rel_mask = _mm_and_ps(states_c->envelope, _rel_mask);

        states_c->g_c = div_f(coeffs_c->threshold, states_c->envelope);

        states_c->g_c.m128_f32[0] = pow_f(states_c->g_c.m128_f32[0], coeffs_c->_rep_ratio.m128_f32[0]);
        states_c->g_c.m128_f32[1] = pow_f(states_c->g_c.m128_f32[1], coeffs_c->_rep_ratio.m128_f32[1]);
        states_c->g_c.m128_f32[2] = pow_f(states_c->g_c.m128_f32[2], coeffs_c->_rep_ratio.m128_f32[2]);
        states_c->g_c.m128_f32[3] = pow_f(states_c->g_c.m128_f32[3], coeffs_c->_rep_ratio.m128_f32[3]);

        states_c->g_c = add_f(states_c->g_c, _att_mask);

        // gain smoothing
        _att_mask = _mm_cmple_ps(states_c->g_c, states_c->g_sPrev);
        _rel_mask = _mm_andnot_ps(_att_mask, _cond_mask);

        _att_mask = _mm_and_ps(states_c->g_c, _att_mask);
        _rel_mask = _mm_and_ps(states_c->g_c, _rel_mask);

        // gain attack
        _att_mask = mul_f(_att_mask, coeffs_c->_one_g_att);
        _att_mask = mac_f(coeffs_c->alphaAttack, states_c->g_sPrev, _att_mask);

        // Envelope release
        _rel_mask = mul_f(_rel_mask, coeffs_c->_one_g_rel);
        _rel_mask = mac_f(coeffs_c->alphaRelease, states_c->g_sPrev, _rel_mask);

        states_c->g_s = add_f(_att_mask, _rel_mask);
        states_c->g_sPrev = states_c->g_s;



        ((tStereo_compr*)audio_c->audio)[i].L = states_c->g_s.m128_f32[0];
        ((tStereo_compr*)audio_c->audio)[i].R = states_c->g_s.m128_f32[1];
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
/*int32_t compressor_4ср_process(
    void const* coeffs,
    void*       states,
    void*       audio,
    size_t      samples_count)
{
    compressor_coeffs* coeffs_c;
    compressor_states* states_c;
    tStereo_compr* audio_c;
    uint32_t i;

    coeffs_c = (compressor_coeffs*)coeffs;
    states_c = (compressor_states*)states;
    audio_c = (tStereo_compr*)audio;
    float xL_abs;
    float xR_abs;


    for (i = 0; i < samples_count; i++)
    {
        states_c->x = ((tStereo_compr*)audio_c)[i];

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
            states_c->g_c.L = 1;
        }
        else
        {
            states_c->g_c.L = (coeffs_c->threshold * powf((states_c->envelope.L / coeffs_c->threshold), (1.0 / coeffs_c->ratio))) / states_c->envelope.L;
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




        ((tStereo_compr*)audio_c)[i].L = states_c->x.L * states_c->g_m.L;
        ((tStereo_compr*)audio_c)[i].R = states_c->x.R * states_c->g_m.R;
    }
}



/*
compressor_4ch_coeffs* coeffs_4ch_c;
compressor_4ch_states* states_4ch_c;
audio_buf_comp* audio_c;

uint32_t i;

coeffs_4ch_c = (compressor_4ch_coeffs*)coeffs;
states_4ch_c = (compressor_4ch_states*)states;
audio_c = (audio_buf_comp*)audio;

if (coeffs_4ch_c->comp_ch_1_coef.bpass == 0)
{
    compressor_process((compressor_coeffs*)&coeffs_4ch_c->comp_ch_1_coef, (compressor_states*)&states_4ch_c->comp_ch_1_st, (audio_buf_comp*)audio_c->band_4, samples_count);
}

if (coeffs_4ch_c->comp_ch_2_coef.bpass == 0)
{
    compressor_process((compressor_coeffs*)&coeffs_4ch_c->comp_ch_2_coef, (compressor_states*)&states_4ch_c->comp_ch_2_st, (audio_buf_comp*)audio_c->band_2, samples_count);
}

if (coeffs_4ch_c->comp_ch_3_coef.bpass == 0)
{
    compressor_process((compressor_coeffs*)&coeffs_4ch_c->comp_ch_3_coef, (compressor_states*)&states_4ch_c->comp_ch_3_st, (audio_buf_comp*)audio_c->band_3, samples_count);
}

if (coeffs_4ch_c->comp_ch_4_coef.bpass == 0)
{
    compressor_process((compressor_coeffs*)&coeffs_4ch_c->comp_ch_4_coef, (compressor_states*)&states_4ch_c->comp_ch_4_st, (audio_buf_comp*)audio_c->band_4, samples_count);
}
*/