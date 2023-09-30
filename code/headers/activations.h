//
// Created by ruanb on 9/10/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
#define C_PIANO_TRANSCRIPTION_ACTIVATIONS_H

#include "spectrogram.h"
#include "dictionary.h"
#include "matrix.h"

double BetaDivergence(Matrix const *x, Matrix const *y, double beta);

Matrix ComputeConvolution(Dictionary const *dictionary, Matrix const *matrix2, unsigned int t);

Matrix
ComputeActivations(Spectrogram const *input, unsigned int iterations, double beta, double error, Dictionary *dictionary);

void TestActivations();

#endif //C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
