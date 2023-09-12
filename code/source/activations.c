//
// Created by ruanb on 9/10/2023.
//

#include "activations.h"
#include "math.h"
#include "stdio.h"

double BetaDivergence(const Matrix *x, const Matrix *y, double beta) {
    Matrix d = CreateMatrix(x->rows, x->cols);

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


Matrix ComputeActivations(const Spectrogram *X, unsigned int iterations, double beta, double maximum_error,
                          Dictionary *dictionary) {

    unsigned int r = dictionary->shape[2];
    unsigned int ncol = dictionary->shape[1];
    unsigned int T = dictionary->shape[0];

    Matrix activations = CreateMatrix(dictionary->shape[2], dictionary->shape[1]);

    double gamma = 1;

    if (beta != 1) {
        fprintf(stderr, "ComputeActivations: beta != 1 not implemented yet\n");
    }


    Matrix conved = ComputeConvolution(dictionary, &activations, T);
    double error_int = BetaDivergence(&X->matrix, &conved, beta);

//    TODO: denoms_cropped_for_end
//    denom_all_col = np.sum(np.dot(W[t].T, np.ones([W.shape[1], ncol])) for t in
//    range(T))
//
//    denoms_cropped_for_end = [None]
//    for j in range(1, T + 1):
//    tab = np.sum(np.dot(W[i].T, np.ones(W[i].shape[0])) for i in range(j))
//    denoms_cropped_for_end.append(tab)
    Matrix convolutions[T];

    Matrix ones = CreateMatrix(dictionary->shape[1], ncol);
    FillMatrix(&ones, 1);

    for (int t = 0; t < T; t++) {
        Spectrogram tspec = GetSpectrogramFromDictionary(dictionary, 0, t);
        Matrix tspec_transposed = Transpose(&tspec.matrix);

        convolutions[t] = MatrixMultiply(&tspec_transposed, &ones);

        DestroySpectrogram(&tspec);
        DestroyMatrix(&tspec_transposed);
    }

    DestroyMatrix(&ones);

    Matrix denom_all_col = SumMatrices(convolutions, T);

//    Spectrogram denoms_cropped_for_end[T];
//
//    for (int i = 0; i < T; i++) {
//        Spectrogram tspec = GetSpectrogramFromDictionary(dictionary, 0, i);
//        Spectrogram tspec_transposed = Transpose(&tspec);
//
//        Spectrogram ones = CreateSpectrogram(dictionary->shape[0], ncol);
//        FillSpectrogram(&ones, 1);
//
//        convolutions[i] = MatrixMultiply(&tspec_transposed, &ones);
//
//        DestroySpectrogram(&tspec);
//        DestroySpectrogram(&tspec_transposed);
//        DestroySpectrogram(&ones);
//    }


    for (int i = 0; i < T; i++) {
        DestroyMatrix(&convolutions[i]);
    }

    unsigned int iteration = 0;
    double obj, obj_prev = 0;

    while (iteration < iterations) {
        Matrix A = ComputeConvolution(dictionary, &activations, T);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                A.array[i][i] = pow(A.array[i][j], (beta - 2));

                if (A.array[i][j] < 1e-8) {
                    A.array[i][j] = 1e-8;
                }

                A.array[i][j] *= X->matrix.array[i][j];
            }
        }

        Matrix X_hadamard_A_padded = CreateMatrix(A.rows, A.cols + T);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                X_hadamard_A_padded.array[i][j] = A.array[i][j];
            }
        }

        Matrix num = CreateMatrix(r, A.cols);

        DestroyMatrix(&A);

        for (int t = 0; t < T; t++) {
            Matrix W_at_t = GetMatrixFromDictionary(dictionary, 0, t);
            Matrix transposed = Transpose(&W_at_t);

            DestroyMatrix(&W_at_t);

            Matrix X_hadamard_A_windowed = CreateMatrix(X_hadamard_A_padded.rows, ncol);

            Matrix multiplied = MatrixMultiply(&transposed, &X_hadamard_A_windowed);

            DestroyMatrix(&X_hadamard_A_windowed);
            DestroyMatrix(&transposed);

            for (int i = 0; i < multiplied.rows; i++) {
                for (int j = 0; j < multiplied.cols; j++) {
                    num.array[i][j] += multiplied.array[i][j];
                }
            }

//            TODO: another few denom_cropped lines

            DestroyMatrix(&conved);
            conved = ComputeConvolution(dictionary, &activations, T);

            obj = BetaDivergence(&X->matrix, &conved, beta);

            if ((fabs(obj - obj_prev) / error_int) < maximum_error) {
                printf("ComputeActivations: Converged sufficiently\n");
                break;
            }

            obj_prev = obj;
            iteration++;
        }
    }

    return activations;
}

Matrix ComputeConvolution(const Dictionary *dictionary, const Matrix *matrix2, unsigned int t) {
    Matrix convolutions[t];

    for (int i = 0; i < t; i++) {
        Matrix tspec = GetMatrixFromDictionary(dictionary, 0, i);

        Matrix shifted = ShiftMatrix(matrix2, i);
        convolutions[i] = MatrixMultiply(&tspec, &shifted);

        DestroyMatrix(&tspec);
        DestroyMatrix(&shifted);
    }

    Matrix conv_sum = SumMatricesAlongAxis(convolutions, t, 0);

    for (int i = 0; i < t; i++) {
        DestroyMatrix(&convolutions[i]);
    }

    return conv_sum;
}