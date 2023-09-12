//
// Created by ruanb on 9/6/2023.
//

#include "spectrogram.h"
#include "stdlib.h"
#include "stdio.h"

Spectrogram CreateSpectrogram(unsigned int nRows, unsigned int nCols) {
    Spectrogram spectrogram;
    spectrogram.array = (double **) malloc(nRows * sizeof(double *));

    for (int i = 0; i < nRows; i++) {
//        spectrogram.array[i] = (double *) malloc(nCols * sizeof(double));
        spectrogram.array[i] = (double *) calloc(nCols, sizeof(double));
    }

    spectrogram.rows = nRows;
    spectrogram.cols = nCols;

    return spectrogram;
}

void spectrogramTest() {
}

void NormaliseSpectrogram(Spectrogram *spectrogram) {
    double max = 0;
    for(int r = 0; r < spectrogram->rows; r++){
        for(int c = 0; c < spectrogram->cols; c++){
            if(spectrogram->array[r][c] > max)
                max = spectrogram->array[r][c];
        }
    }

    for(int r = 0; r < spectrogram->rows; r++){
        for(int c = 0; c < spectrogram->cols; c++){
            spectrogram->array[r][c] /= max;
        }
    }
}

Spectrogram HardFilterSpectrogram(Spectrogram const * spectrogram, unsigned int numNewRows) {
    Spectrogram filtered = CreateSpectrogram(numNewRows, spectrogram->cols);

    for(int r = 0; r < numNewRows; r++){
        for(int c = 0; c < spectrogram->cols; c++){
            filtered.array[r][c] = spectrogram->array[r][c];
        }
    }

    return filtered;
}

void SaveSpectrogramToCSV(const char *filename, Spectrogram* spectrogram) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "SaveSpectrogramToCSV: Error opening file");
        return;
    }

    // Write the array data to the CSV file
    for (int i = 0; i < spectrogram->rows; i++) {
        for (int j = 0; j < spectrogram->cols; j++) {
            fprintf(file, "%.6f", spectrogram->array[i][j]); // Adjust the format specifier as needed
            if (j < spectrogram->cols - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

Spectrogram ShiftSpectrogram(Spectrogram const *spectrogram, unsigned int numShifts) {
    Spectrogram shifted = CreateSpectrogram(spectrogram->rows, spectrogram->cols);

    for (int i = 0; i < spectrogram->rows; i++) {
        for (int j = 0; j < spectrogram->cols; j++) {
            if (j + numShifts < spectrogram->cols) {
                shifted.array[i][j + numShifts] = spectrogram->array[i][j];
            }
        }
    }

    for (int i = 0; i < spectrogram->rows; i++) {
        for (int j = 0; j < numShifts; j++) {
            shifted.array[i][j ] = 0;

        }
    }

    return shifted;

}

Spectrogram Transpose(const Spectrogram *spectrogram) {
    Spectrogram transposed = CreateSpectrogram(spectrogram->cols, spectrogram->rows);

    for (int i = 0; i < spectrogram->rows; i++) {
        for (int j = 0; j < spectrogram->cols; j++) {
            transposed.array[j][i] = spectrogram->array[i][j];
        }
    }

    return transposed;
}

void DestroySpectrogram(Spectrogram *spectrogram) {
    for (int i = 0; i < spectrogram->rows; i++) {
        free(spectrogram->array[i]);
    }

    free(spectrogram->array);
}
