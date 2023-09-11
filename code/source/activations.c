//
// Created by ruanb on 9/10/2023.
//

#include "activations.h"
#include "math.h"
#include "stdio.h"

double BetaDivergence(const Spectrogram *x, const Spectrogram *y, double beta) {
    Spectrogram d = CreateSpectrogram(x->rows, x->cols);

    double sum = 0;

    for(int i = 0; i < x->rows; i++){
        for(int j = 0; j < x->cols; j++){
            if(y->array[i][j] < 1e-8){
                y->array[i][j] = 1e-8;
            }
            if(beta == 1){
                d.array[i][j] = x->array[i][j]*log(x->array[i][j]/y->array[i][j]) - x->array[i][j] + y->array[i][j];

                sum += d.array[i][j];

            } else {
                fprintf(stderr, "BetaDivergence: beta != 1 not implemented yet\n");
            }
        }
    }

    return sum;

}

Spectrogram MatrixMultiply(const Spectrogram *a, const Spectrogram *b) {
    Spectrogram result = CreateSpectrogram(a->rows, b->cols);

    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            result.array[i][j] = 0;
            for (int k = 0; k < a->cols; k++) {
                result.array[i][j] += a->array[i][k] * b->array[k][j];
            }
        }
    }

    return result;
}

// TODO: optimisation is probably possible here
Spectrogram SumSpectrograms(Spectrogram *a, unsigned int numSpectrograms, unsigned int axis) {
    Spectrogram result = CreateSpectrogram(a[0].rows, a[0].cols);

    for(int r = 0; r < result.rows; r++){
        for(int c = 0; c < result.cols; c++){
            result.array[r][c] = 0;

            for(int i = 0; i < numSpectrograms; i++){
                if(axis == 0){
                    result.array[r][c] += a[i].array[r][c];
                } else if(axis == 1){
                    result.array[r][c] += a[i].array[c][r];
                }
            }
        }
    }

    return result;
}

Spectrogram ComputeActivations(const Spectrogram *input, unsigned int iterations, double beta, double maximum_error,
                               Dictionary *dictionary, unsigned int t) {
    Spectrogram activations = CreateSpectrogram(dictionary->shape[2], dictionary->shape[1]);

    double gamma = 1;

    if(beta != 1){
        fprintf(stderr, "ComputeActivations: beta != 1 not implemented yet\n");
    }

    Spectrogram convolutions[t];

    for(int i = 0; i < t; i++){
        convolutions[i] = MatrixMultiply()

    }

    double error = BetaDivergence(input, &MatrixMultiply(dictionary->data[0], &activations), beta);


    return activations;
}


