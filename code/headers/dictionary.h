//
// Created by ruanb on 9/4/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_DICTIONARY_H
#define C_PIANO_TRANSCRIPTION_DICTIONARY_H

#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "spectrogram.h"

typedef struct {
    hsize_t shape[3];
    double*** data;
} Dictionary;

void LoadDictionary(Dictionary *dictionary, char *filename);
void AllocateDictionaryMemory(Dictionary* dictionary);
void PrintDictionary(Dictionary* dictionary);

// TODO: this function, remember dict dims are (10, 4097, 88)
Dictionary HardFilterSpectrograms(Dictionary* dictionary, unsigned int numNewRows);

Spectrogram GetSpectrogramFromDictionary(Dictionary* dictionary, unsigned int axis, unsigned int index);

void dictionaryTest();


#endif //C_PIANO_TRANSCRIPTION_DICTIONARY_H
