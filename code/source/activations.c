//
// Created by ruanb on 9/10/2023.
//

#include "activations.h"
#include "math.h"
#include "stdio.h"

double BetaDivergence(const Spectrogram *x, const Spectrogram *y, double beta) {
    Spectrogram d = CreateSpectrogram(x->rows, x->cols);

    double sum = 0;

    for (int i = 0; i < x->rows; i++) {
        for (int j = 0; j < x->cols; j++) {
            if (y->array[i][j] < 1e-8) {
                y->array[i][j] = 1e-8;
            }
            if (beta == 1) {
                d.array[i][j] = x->array[i][j] * log(x->array[i][j] / y->array[i][j]) - x->array[i][j] + y->array[i][j];

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

    for (int r = 0; r < result.rows; r++) {
        for (int c = 0; c < result.cols; c++) {
            for (int i = 0; i < numSpectrograms; i++) {
                if (axis == 0) {
                    result.array[r][c] += a[i].array[r][c];
                } else if (axis == 1) {
                    result.array[r][c] += a[i].array[c][r];
                }
            }
        }
    }

    return result;
}

Spectrogram ComputeActivations(const Spectrogram *X, unsigned int iterations, double beta, double maximum_error,
                               Dictionary *dictionary) {

    unsigned int r = dictionary->shape[2];
    unsigned int ncol = dictionary->shape[1];
    unsigned int T = dictionary->shape[0];

    Spectrogram activations = CreateSpectrogram(dictionary->shape[2], dictionary->shape[1]);

    double gamma = 1;

    if (beta != 1) {
        fprintf(stderr, "ComputeActivations: beta != 1 not implemented yet\n");
    }

    Spectrogram conved = ComputeConvolution(dictionary, &activations, beta, T);
    double error_int = BetaDivergence(X, &conved, beta);

//    TODO: denoms_cropped_for_end

    unsigned int iteration = 0;
    double obj, obj_prev = 0;

    while (iteration < iterations) {
        Spectrogram A = ComputeConvolution(dictionary, &activations, beta, T);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                A.array[i][i] = pow(A.array[i][j], (beta - 2));

                if (A.array[i][j] < 1e-8) {
                    A.array[i][j] = 1e-8;
                }

                A.array[i][j] *= X->array[i][j];
            }
        }

        Spectrogram X_hadamard_A_padded = CreateSpectrogram(A.rows, A.cols + T);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                X_hadamard_A_padded.array[i][j] = A.array[i][j];
            }
        }

        Spectrogram num = CreateSpectrogram(r, A.cols);

        for(int t = 0 ; t < T; t++){
            Spectrogram W_at_t = GetSpectrogramFromDictionary(dictionary, 0, t);
            Spectrogram transposed = Transpose(&W_at_t);

            Spectrogram X_hadamard_A_windowed = CreateSpectrogram(X_hadamard_A_padded.rows, ncol);

            Spectrogram multiplied = MatrixMultiply(&transposed, &X_hadamard_A_windowed);

            for(int i = 0; i < multiplied.rows; i++){
                for(int j = 0; j < multiplied.cols; j++){
                    num.array[i][j] += multiplied.array[i][j];
                }
            }

//            TODO: another few denom_cropped lines

            DestroySpectrogram(&conved);
            conved = ComputeConvolution(dictionary, &activations, beta, T);

            obj = BetaDivergence(X, &conved, beta);

            if((fabs(obj - obj_prev) / error_int ) < maximum_error){
                printf("ComputeActivations: Converged sufficiently\n");
                break;
            }

            obj_prev = obj;
            iteration++;
        }
    }

    return activations;
}

Spectrogram
ComputeConvolution(Dictionary const *dictionary, Spectrogram const *activations, double beta, unsigned int t) {
    Spectrogram convolutions[t];

    for (int i = 0; i < t; i++) {
        Spectrogram tspec = GetSpectrogramFromDictionary(dictionary, 0, i);
        Spectrogram shifted = ShiftSpectrogram(activations, i);
        convolutions[i] = MatrixMultiply(&tspec, &shifted);
    }

    Spectrogram conv_sum = SumSpectrograms(convolutions, t, 0);

    return conv_sum;
}


