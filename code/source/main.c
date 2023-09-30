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
#include "transcription.h"


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
//    GraphTest();
//    matrixTest();

    const char *piano_W[3] = {"AkPnBsdf", "AkPnStgb", "AkPnBcht"};

    Dictionary dictionaries[3];

    // Iterate through piano_W
    for (int i = 0; i < 3; i++) {
        char W_persisted_name[256]; // Adjust the size as needed

        // Create the formatted string
        snprintf(W_persisted_name, sizeof(W_persisted_name),
                 "%sconv_dict_piano_%s_beta_%d_T_10_init_%s_%s_%d_itmax_%d_intensity_%s.h5",
                 "C:\\Users\\ruanb\\OneDrive\\Desktop\\Piano Transcripton\\Piano Transcription (C)\\data_persisted\\dictionaries\\",
                 piano_W[i],
                 1,
                 "L1",
                 "stft",
                 4096,
                 500,
                 "M");

        LoadDictionary(&dictionaries[i], W_persisted_name);
    }

    const char *filename = "MAPS_MUS-alb_se3_AkPnBcht.wav";

    WavFile wav = ReadWav(filename);

    DynamicArray mono = StereoToMono(&wav, "average");

    Spectrogram spec = STFT(&mono, 4096, 882, 8192, 5, 44100);
    Spectrogram filtered = HardFilterSpectrogram(&spec, 1500);

    double delay = GetDelay(&filtered, 0.05);

    Matrix activations = LoadMatrixFromCSV("activations.csv");

    if(activations.rows == 0 && activations.cols == 0){
        activations = GetActivationsFromFile(filename, &dictionaries[2]);
    }

    double threshold = GetThreshold(&activations);

    Matrix notes = TranscribeNotesFromActivations(&activations, threshold, 0.02);

    int finalIndex = 0;

    while (true) {
        if (notes.array[finalIndex][1] == 0) {
            break;
        }
        finalIndex++;
    }

    Matrix finalNotes = CreateMatrix(finalIndex, 2);

    double est_min = 99;

    for (int i = 0; i < finalIndex; i++) {
        if (notes.array[i][0] < est_min) {
            est_min = notes.array[i][0];
        }

        finalNotes.array[i][0] = notes.array[i][0];
        finalNotes.array[i][1] = notes.array[i][1];
    }

    for (int i = 0; i < finalIndex; i++) {
        finalNotes.array[i][0] -= est_min;
    }

    Matrix refs = LoadRefsFromFile(
            "C:\\Users\\ruanb\\OneDrive\\Desktop\\Piano Transcripton\\Piano transcription\\MAPS\\AkPnBcht\\MUS\\MAPS_MUS-alb_se3_AkPnBcht.txt",
            5 - delay);

    double ref_min = 99;

    for (int i = 0; i < refs.rows; i++) {
        if (refs.array[i][0] < ref_min) {
            ref_min = refs.array[i][0];
        }
    }

    for (int i = 0; i < refs.rows; i++) {
        refs.array[i][0] -= ref_min;
    }

    EvaluateTranscription(&refs, &finalNotes);


    return 0;
}



