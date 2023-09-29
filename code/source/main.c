#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../headers/dictionary.h"
#include "../headers/dynamicarray.h"
#include "../headers/stft.h"
#include "math.h"
#include "evaluation.h"


#include "wav.h"
#include "activations.h"


void test() {
    WavFile wav = ReadWav(
            "C:/Users/ruanb/OneDrive/Desktop/Piano Transcripton/Piano transcription/MAPS/AkPnBcht/ISOL/NO/MAPS_ISOL_NO_F_S0_M23_AkPnBcht.wav");

    DynamicArray mono = StereoToMono(&wav, "average");

    SaveArrayToCSV("mono.csv", &mono);

    Spectrogram spec = STFT(&mono, 4096, 882, 8192, 1, 44100);

    SaveSpectrogramToCSV("testSpec.csv", &spec);
}

int main() {
    // Define the parameters

//    dictionaryTest();
//    stftTest();
//    TestActivations();
//    test();
    GraphTest();

    const char *dictionary_directory = "C:\\Users\\ruanb\\OneDrive\\Desktop\\Piano Transcripton\\Piano Transcription (C)\\data_persisted\\dictionaries\\";
    const char *piano_W[3] = {"AkPnBsdf", "AkPnStgb", "AkPnBcht"};

    Dictionary dictionaries[3];

    int beta = 1; // Replace with the actual value
    const char *init = "L1"; // Replace with the actual value
    const char *spec_type = "stft"; // Replace with the actual value
    int num_points = 4096; // Replace with the actual value
    int itmax_W = 500; // Replace with the actual value
    const char *note_intensity = "M"; // Replace with the actual value



    // Iterate through piano_W
    for (int i = 0; i < 3; i++) {
        char W_persisted_name[256]; // Adjust the size as needed

        // Create the formatted string
        snprintf(W_persisted_name, sizeof(W_persisted_name),
                 "%sconv_dict_piano_%s_beta_%d_T_10_init_%s_%s_%d_itmax_%d_intensity_%s.h5",
                 dictionary_directory,
                 piano_W[i],
                 beta,
                 init,
                 spec_type,
                 num_points,
                 itmax_W,
                 note_intensity);

        LoadDictionary(&dictionaries[i], W_persisted_name);
    }

    const char *filename = "MAPS_MUS-alb_se3_AkPnBcht.wav";

    WavFile wav = ReadWav(filename);

    DynamicArray mono = StereoToMono(&wav, "average");

    Spectrogram spec = STFT(&mono, 4096, 882, 8192, 5, 44100);
    Spectrogram filtered = HardFilterSpectrogram(&spec, 1500);

    Dictionary aaaa = HardFilterSpectrograms(&dictionaries[2], 1500);

    NormaliseDictionary(&aaaa);

    DestroySpectrogram(&spec);

    Matrix activations = ComputeActivations(&filtered, 20, 1, 0.1, &aaaa);

    Matrix notes = TranscribeNotesFromActivations(&activations, 0.08, 0.02);

    SaveMatrixToCSV("notes.csv", &notes);

    return 0;
}



