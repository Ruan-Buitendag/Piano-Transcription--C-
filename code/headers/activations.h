//
// Created by ruanb on 9/10/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
#define C_PIANO_TRANSCRIPTION_ACTIVATIONS_H

#include "spectrogram.h"

double BetaDivergence(Spectrogram const *x, Spectrogram const *y, double beta);

#endif //C_PIANO_TRANSCRIPTION_ACTIVATIONS_H
