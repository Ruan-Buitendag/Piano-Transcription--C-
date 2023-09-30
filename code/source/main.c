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

#define poo 0

int main() {
    // Define the parameters

//    dictionaryTest();
//    stftTest();
//    TestActivations();
//    test();
//    GraphTest();
//    matrixTest();

    if (poo) {
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

        SaveMatrixToCSV("activations.csv", &activations);
    } else {

        const char *filename = "MAPS_MUS-alb_se3_AkPnBcht.wav";

        WavFile wav = ReadWav(filename);

        DynamicArray mono = StereoToMono(&wav, "average");

        Spectrogram spec = STFT(&mono, 4096, 882, 8192, 5, 44100);
        Spectrogram filtered = HardFilterSpectrogram(&spec, 1500);

        double delay = GetDelay(&filtered, 0.05);

        Matrix activations = LoadMatrixFromCSV("activations.csv");

        double threshold = GetThreshold(&activations);

        Matrix notes = TranscribeNotesFromActivations(&activations, threshold, 0.02);

//    SaveMatrixToCSV("notes.csv", &notes);

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
    }

    return 0;
}



