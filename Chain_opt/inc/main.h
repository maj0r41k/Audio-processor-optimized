#ifndef __MAIN_H_
#define __MAIN_H_


//#define FIR_float
//#define IIR_float
//#define GAIN
//#define COMPRESSOR
//#define EQ
//#define CROSS
#define CHAIN

#define BUFFER_SIZE 4096u
#define BUFFER_SHIFT 12u

#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"

#include "fixedpoint.h"

#include "effect_control.h"
#include "effect_process.h"

typedef struct {
    uint8_t    riff[4];                    // RIFF string
    uint32_t   overall_size;  	            // overall size of file in bytes
    uint8_t    wave[4];                    // WAVE string
} riffHeader;

typedef struct {
    uint8_t    fmt_chunk_marker[4];        // fmt string with trailing null char
    uint32_t   length_of_fmt;              // length of the format data
    uint16_t   format_type;                // format type
    uint16_t   channels;                   // # of channels
    uint32_t   sample_rate;                // sampling rate 
    uint32_t   byterate;                   // SampleRate * NumChannels * BitsPerSample/8
    uint16_t   block_align;                // NumChannels * BitsPerSample/8
    uint16_t   bits_per_sample;            // bits per sample, 8- 8bits, 16- 16 bits
}fmtHeader;

typedef struct {
    uint8_t    data_chunk_header[4];           // DATA string or FLLR string
    uint32_t   data_size;                  // NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
}dataHeader;

typedef struct {
    riffHeader riff;
    fmtHeader format;
    dataHeader data;
}wav_t;

typedef struct {
    fmtHeader* fmt;
    dataHeader* data;
    void* audio;
    float Q;
    float freq;
    float gain;     //in dB
}effect_parameters;


typedef struct {
    //fmtHeader* fmt;
    //dataHeader* data;
    //void* audio;
    float threshold;
    float ratio;
    float tauAttack;
    float tauRelease;
    float makeUpGain;
    
    float envAtt;
    float envRel;

    uint32_t bpass;
}effect_params_compressor;

typedef struct {
    float freq;
    float gain;
    float Q;
    int type;
}equalizer;

typedef struct {
    equalizer eq[10];
    float sampleRate;

    //fmtHeader* fmt;
    //dataHeader* data;
    //void* audio;
}effect_params_equalizer;


typedef struct {
    my_float cutoff_1;
    my_float cutoff_2;
    my_float cutoff_3;

    //fmtHeader* fmt;
    //dataHeader* data;
    //void* audio;
}effect_params_crossover;

typedef struct {
    void* band_1;
    void* band_2;
    void* band_3;
    void* band_4;
}band_buf;

typedef struct {
    void* audio;
    band_buf cross_b; //* absence 
}audio_buf;

typedef struct {
    FILE* source;
    FILE* destin;
    wav_t* wav;

    audio_buf* audio;

    effect_params_equalizer equalizer;
    effect_params_crossover crossover;
    effect_params_compressor comp_ch1;
    effect_params_compressor comp_ch2;
    effect_params_compressor comp_ch3;
    effect_params_compressor comp_ch4;

}effect_parameter_chain;

void readHeader(riffHeader* ptrRIFF, fmtHeader* ptrFMT, dataHeader* ptrDATA, FILE* ptrWavFile, size_t* numRead, int* diviation);

void effect_crossover(effect_params_crossover cross);
void effect_equalizer(effect_params_equalizer eqzr);
void effect_compressor(effect_params_compressor effect_par_comp);
void effect_fir(effect_parameters effect_params);
void effect_iir(effect_parameters effect_params);
void effect_gain(effect_parameters effect_params);

void effect_chain(effect_parameter_chain* chain);

void set_params(void * params);

#endif // !__MAIN_H_



