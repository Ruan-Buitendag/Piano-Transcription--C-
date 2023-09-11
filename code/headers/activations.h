//
// Created by ruanb on 9/10/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
#define C_PIANO_TRANSCRIPTION_ACTIVATIONS_H

#include "spectrogram.h"
#include "dictionary.h"

double BetaDivergence(Spectrogram const *x, Spectrogram const *y, double beta);
Spectrogram MatrixMultiply(Spectrogram const *a, Spectrogram const *b);
Spectrogram SumSpectrograms(Spectrogram *a, unsigned int numSpectrograms, unsigned int axis);

Spectrogram ComputeActivations(Spectrogram const *input, unsigned int iterations, double beta, double error, Dictionary *dictionary, unsigned int t);

#endif //C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
