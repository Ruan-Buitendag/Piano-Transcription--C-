//
// Created by ruanb on 9/10/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
#define C_PIANO_TRANSCRIPTION_ACTIVATIONS_H

#include "spectrogram.h"
#include "dictionary.h"
#include "matrix.h"

double BetaDivergence(Matrix const *x, Matrix const *y, double beta);

Spectrogram
ComputeActivations(Spectrogram const *input, unsigned int iterations, double beta, double error, Dictionary *dictionary);

#endif //C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
