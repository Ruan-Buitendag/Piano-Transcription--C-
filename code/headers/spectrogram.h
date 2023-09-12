//
// Created by ruanb on 9/6/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_SPECTROGRAM_H
#define C_PIANO_TRANSCRIPTION_SPECTROGRAM_H

typedef struct {
    double **array;
    unsigned int rows;
    unsigned int cols;
} Spectrogram;

Spectrogram CreateSpectrogram(unsigned int nRows, unsigned int nCols);
void DestroySpectrogram(Spectrogram *spectrogram);

Spectrogram HardFilterSpectrogram(Spectrogram const *spectrogram, unsigned int numNewRows);
void NormaliseSpectrogram(Spectrogram *spectrogram);
Spectrogram ShiftSpectrogram(Spectrogram const *spectrogram, unsigned int numShifts);
Spectrogram Transpose(Spectrogram const *spectrogram);

void SaveSpectrogramToCSV(const char *filename, Spectrogram *spectrogram);

void spectrogramTest();

#endif //C_PIANO_TRANSCRIPTION_SPECTROGRAM_H
