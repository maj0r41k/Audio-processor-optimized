#include "main.h"
#include "cJSON.h"
#include "time.h"
#include "stdio.h"

int main()
{
    float begin_ms;
    float end_ms;
    float time;

    static FILE* ptrWavFile;
    FILE newWavFIR;
    FILE* ptrNewWavFIR = &newWavFIR;
    int closeFile;
    size_t numRead;
    size_t numWrite;

    char* fileAdd;
    char* newFileFIR;
    int diviation;

    wav_t wav;
    size_t audio_bytes;
    audio_buf bands;
    size_t samples_count;

    effect_parameter_chain chain;

    numRead = 0;
    numWrite = 0;
    diviation = 0;
    audio_bytes = (uint32_t)BUFFER_SIZE;

    bands.audio = malloc(audio_bytes);
    bands.cross_b.band_1 = malloc(audio_bytes);
    bands.cross_b.band_2 = malloc(audio_bytes);
    bands.cross_b.band_3 = malloc(audio_bytes);
    bands.cross_b.band_4 = malloc(audio_bytes);

    fileAdd = "../../Chain_opt/test_signal/a.wav";
    newFileFIR = "../../Chain_opt/test_signal/out/fir_opt.wav";

    ptrWavFile = fopen(fileAdd, "rb");      // Open existance .wav file  
    ptrNewWavFIR = fopen(newFileFIR, "w+b");   

    readHeader((riffHeader*)(&wav.riff), (fmtHeader*)(&wav.format), (dataHeader*)(&wav.data), ptrWavFile, &numRead, &diviation);

    //hProcess = OpenProcess(PROCESS_QUERY_INFORMATION, TRUE, 20);
    //GetProcessTimes(GetCurrentProcess(), &lpCreationTime, &lpExistTime, &lpKernelTime, &lpUserTime);
    //FileTimeToSystemTime(&lpUserTime, &begin);


    numWrite = fwrite((riffHeader*)(&wav.riff), sizeof(uint8_t), sizeof(riffHeader), ptrNewWavFIR);
    numWrite = fwrite((fmtHeader*)(&wav.format), sizeof(uint8_t), sizeof(fmtHeader), ptrNewWavFIR);
    numWrite = fwrite((dataHeader*)(&wav.data), sizeof(uint8_t), sizeof(dataHeader), ptrNewWavFIR);

  
    chain = (effect_parameter_chain) {
        (FILE*)ptrWavFile,
        (FILE*)ptrNewWavFIR,
        (wav_t*)(&wav),

        (audio_buf*)(&bands),

        (effect_params_equalizer) {
            {
                (equalizer) {20, -3, 1, 3},
                (equalizer) {250, 0, 0.5, 5},
                (equalizer) {500, -3, 0.5, 2},
                (equalizer) {1000, -3, 4, 2},
                (equalizer) {2000, -3, 20, 2},
                (equalizer) {4000, -3, 4, 2},
                (equalizer) { 8000, -3, 10, 2},
                (equalizer) {12000, -3, 8, 2},
                (equalizer) { 14000, -3, 0.1, 4},
                (equalizer) {18000, -12, 1, 0}
            },
                wav.format.sample_rate
        },

        (effect_params_crossover){500, 5000, 1000},

        (effect_params_compressor) {-3, 2, 30, 40, 0, 0, 30, 0},
        (effect_params_compressor) {-3, 2, 30, 40, 0, 0, 30, 0},
        (effect_params_compressor) {-3, 2, 30, 40, 0, 0, 30, 0},
        (effect_params_compressor) {-3, 2, 30, 40, 0, 0, 30, 0}
    };

    //samples_count = wav.data.data_size / wav.format.block_align;
    //numRead = fread(bands.audio, 1, wav.data.data_size, ptrWavFile);

#ifdef CHAIN
    effect_chain(&chain);
#endif // Equalizer

    //numWrite = fwrite(bands.audio, 1, wav.data.data_size, ptrNewWavFIR);

    closeFile = fclose(ptrWavFile);
    closeFile = fclose(ptrNewWavFIR);

    return 0;
}

void effect_chain(effect_parameter_chain* chain)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    size_t numRead;
    size_t numWrite;

    size_t i;           // counter
    size_t loops;       // number of full cycles of 4096 bytes procession
    size_t last;        // number os 

    float begin_ms;
    float end_ms;
    float time = 0.0;

    loops = chain->wav->data.data_size >> BUFFER_SHIFT;
    last = chain->wav->data.data_size - (loops << BUFFER_SHIFT);

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);
    coeffs = malloc(coeffs_bytes);
    states = malloc(states_bytes);

    effect_control_initialize(params, coeffs, chain->wav->format.sample_rate);
    effect_reset(coeffs, states);
    //
    set_params(params);
    //
    effect_update_coeffs(params, coeffs);

    clock_t beginp = clock();

    //GetProcessTimes(GetCurrentProcess(), &lpCreationTime, &lpExistTime, &lpKernelTime, &lpUserTime);
    //FileTimeToSystemTime(&lpUserTime, &begin);

    for (i = 0; i < loops; i++)
    {
        numRead = fread((chain->audio)->audio, 1, BUFFER_SIZE, chain->source);
        effect_process(coeffs, states, chain->audio, 512);
        numWrite = fwrite((chain->audio)->audio, 1, BUFFER_SIZE, chain->destin);
    }
    if (last != 0)
    {
        numRead = fread((chain->audio)->audio, 1, last, chain->source);
        effect_process(coeffs, states, chain->audio, (last/chain->wav->format.block_align));
        numWrite = fwrite((chain->audio)->audio, 1, last, chain->destin);
    }

    //GetProcessTimes(GetCurrentProcess(), &lpCreationTime, &lpExistTime, &lpKernelTime, &lpUserTime);
    //FileTimeToSystemTime(&lpUserTime, &end);

    //begin_ms = (begin.wMinute * 60000) + (begin.wSecond*1000) + begin.wMilliseconds;
    //end_ms = (end.wMinute * 60000) + (end.wSecond * 1000) + end.wMilliseconds;
    //time = (end_ms - begin_ms)/1000;

    printf("%lf\n", time);

    clock_t endp = clock();
    printf("%lf\n", (double)(endp - beginp) / CLOCKS_PER_SEC);
}

void set_params(void * params)
{
    FILE * js = fopen("../../Chain_opt/chain_float/parser.json", "r");

    fseek(js, 0, SEEK_END);
    size_t size = ftell(js);
    fseek(js, 0, SEEK_SET);

    char * buffer = malloc(size);
    memset(buffer, 0, size);C:
    fread(buffer, size, 1, js);

    const cJSON *band = NULL;
    const cJSON *bands = NULL;

    cJSON *json = cJSON_Parse(buffer);

    bands = cJSON_GetObjectItemCaseSensitive(json, "eq_params");
    cJSON_ArrayForEach(band, bands)
    {
        cJSON *id = cJSON_GetObjectItemCaseSensitive(band, "id");
        cJSON *value = cJSON_GetObjectItemCaseSensitive(band, "val");

        effect_set_parameter(params, id->valueint, (float)value->valuedouble);
    }
}



/*void effect_crossover(effect_params_crossover cross)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer

    effect_control_initialize(params, coeffs, cross.fmt->sample_rate);
    effect_reset(coeffs, states);

    effect_set_parameter(params, 0, cross.cutoff_1);
    effect_set_parameter(params, 1, cross.cutoff_2);
    effect_set_parameter(params, 2, cross.cutoff_3);
    effect_set_parameter(params, 3, cross.fmt->sample_rate);

    effect_update_coeffs(params, coeffs);

    effect_process(coeffs, states, cross.audio, (cross.data->data_size / cross.fmt->block_align));
}*/

/*void effect_equalizer(effect_params_equalizer eqzr)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    size_t i;
    size_t band;
    size_t param_num;

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer

    effect_control_initialize(params, coeffs, eqzr.sampleRate);
    effect_reset(coeffs, states);

    for (i = 0; i < 39; i++)
    {
        band = i >> 2;
        param_num = i - (4 * band);

        switch (param_num)
        {
            case 0:
            {
                effect_set_parameter(params, i, eqzr.eq[band].freq);
                break;
            }

            case 1:
            {
                effect_set_parameter(params, i, eqzr.eq[band].gain);
                break;
            }

            case 2:
            {
                effect_set_parameter(params, i, eqzr.eq[band].Q);
                break;
            }

            case 3:
            {
                effect_set_parameter(params, i, eqzr.eq[band].type);
                break;
            }

            default:
                break;
        }
    }

    effect_set_parameter(params, 40, eqzr.sampleRate);

    effect_update_coeffs(params, coeffs);

    effect_process(coeffs, states, eqzr.audio, (eqzr.data->data_size / eqzr.fmt->block_align));
}*/

/*void effect_compressor(effect_params_compressor effect_par_comp)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer


    effect_control_initialize(params, coeffs, effect_par_comp.fmt->sample_rate);
    effect_reset(coeffs, states);

    effect_set_parameter(params, 0, effect_par_comp.threshold);
    effect_set_parameter(params, 1, effect_par_comp.ratio);
    effect_set_parameter(params, 2, effect_par_comp.tauAttack);
    effect_set_parameter(params, 3, effect_par_comp.tauRelease);
    effect_set_parameter(params, 4, effect_par_comp.makeUpGain);
    effect_set_parameter(params, 5, effect_par_comp.fmt->sample_rate);
    effect_set_parameter(params, 6, effect_par_comp.envAtt);
    effect_set_parameter(params, 7, effect_par_comp.envRel);

    effect_update_coeffs(params, coeffs);

    effect_process(coeffs, states, effect_par_comp.audio, (effect_par_comp.data->data_size / effect_par_comp.fmt->block_align));
}*/

/*void effect_gain(effect_parameters effect_params)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer


    effect_control_initialize(params, coeffs, effect_params.fmt->sample_rate);
    
    effect_set_parameter(params, 0, effect_params.gain);
    effect_update_coeffs(params, coeffs);

    effect_process(coeffs, states, effect_params.audio, (effect_params.data->data_size / effect_params.fmt->block_align));
}*/

/*void effect_fir(effect_parameters effect_params)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    void* params;        
    void* coeffs;     
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer


    effect_control_initialize(params, coeffs, effect_params.fmt->sample_rate);
    effect_reset(coeffs, states);
    printf("coeffs[0] = %d, coeffs[1] = %d, coeffs[2] = %d\n", ((my_sint32*)coeffs)[0], ((my_sint32*)coeffs)[1], ((my_sint32*)coeffs)[2]);
    effect_process(coeffs, states, effect_params.audio, (effect_params.data->data_size/effect_params.fmt->block_align));
}*/

/*void effect_iir(effect_parameters effect_params)
{
    size_t params_bytes;
    size_t coeffs_bytes;
    size_t states_bytes;

    void* params;
    void* coeffs;
    void* states;

    effect_control_get_sizes(&params_bytes, &coeffs_bytes);
    effect_process_get_sizes(&states_bytes);

    params = malloc(params_bytes);            // empty
    coeffs = malloc(coeffs_bytes);            // for h
    states = malloc(states_bytes);            // for circular buffer


    effect_control_initialize(params, coeffs, effect_params.fmt->sample_rate);
    effect_reset(coeffs, states);

    printf("coeffs[0] = %d, coeffs[1] = %d, coeffs[2] = %d\n", ((my_sint32*)coeffs)[0], ((my_sint32*)coeffs)[1], ((my_sint32*)coeffs)[2]);

    effect_set_parameter(params, 0, effect_params.freq);
    effect_set_parameter(params, 1, effect_params.Q);
    effect_set_parameter(params, 2, effect_params.fmt->sample_rate);

    effect_update_coeffs(params, coeffs);

    effect_process(coeffs, states, effect_params.audio, (effect_params.data->data_size / effect_params.fmt->block_align));
}*/

void readHeader(riffHeader* ptrRIFF, fmtHeader* ptrFMT, dataHeader* ptrDATA, FILE* ptrWavFile, size_t* numRead, int* diviation)
{
    size_t offset;
    int fseek_ret;

    numRead = fread(ptrRIFF, sizeof(uint8_t), sizeof(riffHeader), ptrWavFile);
    numRead = fread(ptrFMT, sizeof(uint8_t), sizeof(fmtHeader), ptrWavFile);

    offset = ptrFMT->length_of_fmt + 20;
    fseek_ret = fseek(ptrWavFile, offset, SEEK_SET);

    numRead = fread(ptrDATA, sizeof(uint8_t), sizeof(dataHeader), ptrWavFile);

    ptrRIFF->overall_size = ptrDATA->data_size + ptrFMT->length_of_fmt + 28;

    if (ptrFMT->length_of_fmt != 16)
    {
        *diviation = ptrFMT->length_of_fmt - 16;
        ptrFMT->length_of_fmt = 16;
    }

}

