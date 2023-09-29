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
            if (x->array[i][j] < 1e-8) {
                x->array[i][j] = 1e-8;
            }
            if (beta == 1) {
                d.array[i][j] = x->array[i][j] * log(x->array[i][j] / y->array[i][j]) - x->array[i][j] + y->array[i][j];


                if (isnan(d.array[i][j])) {
                    printf("x: %f, y: %f\n", x->array[i][j], y->array[i][j]);
                }

                sum += d.array[i][j];

            } else {
                fprintf(stderr, "BetaDivergence: beta != 1 not implemented yet\n");
            }
        }
    }

    DestroyMatrix(&d);

    return sum;

}


Matrix ComputeActivations(const Spectrogram *X, unsigned int iterations, double beta, double maximum_error,
                          Dictionary *dictionary) {

    unsigned int r = dictionary->shape[2];
    unsigned int ncol = X->matrix.cols;
    unsigned int T = dictionary->shape[0];

    Matrix activations = CreateMatrix(r, ncol);

    for (int ii = 0; ii < activations.rows; ii++) {
        for (int jj = 0; jj < activations.cols; jj++) {
            activations.array[ii][jj] = (double) rand() / (double) RAND_MAX;
        }
    }

//    SaveMatrixToCSV("X.csv", &X->matrix);

//    FillMatrix(&activations, 1);


    double gamma = 1;

    if (beta != 1) {
        fprintf(stderr, "ComputeActivations: beta != 1 not implemented yet\n");
    }

    Matrix conved = ComputeConvolution(dictionary, &activations, T);


    double error_int = BetaDivergence(&X->matrix, &conved, beta);

    printf("ComputeActivations: Initial error = %f\n", error_int);

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


//    SaveMatrixToCSV("denom_all_col.csv", &denom_all_col);


    for (int i = 0; i < T; i++) {
        DestroyMatrix(&convolutions[i]);
    }

    Matrix denoms_cropped_for_end = CreateMatrix(T, 88);


    for (int j = 1; j < T + 1; j++) {

        Matrix temp = CreateMatrix(T, 88);

        for (int a = 0; a < j; a++) {
            Spectrogram a_spec = GetSpectrogramFromDictionary(dictionary, 0, a);
            for (int cc = 0; cc < a_spec.matrix.cols; cc++) {
                for (int rr = 0; rr < a_spec.matrix.rows; rr++) {
                    temp.array[a][cc] += a_spec.matrix.array[rr][cc];
                }
            }

            DestroySpectrogram(&a_spec);

        }

        for (int cc = 0; cc < 88; cc++) {
            for (int rr = 0; rr < j; rr++) {
                denoms_cropped_for_end.array[j - 1][cc] += temp.array[rr][cc];
            }
        }

        DestroyMatrix(&temp);
    }

//    SaveMatrixToCSV("denoms_cropped_for_end.csv", &denoms_cropped_for_end);


    unsigned int iteration = 0;
    double obj, obj_prev = 0;

    while (iteration < iterations) {
        Matrix A = ComputeConvolution(dictionary, &activations, T);

//        SaveMatrixToCSV("A.csv", &A);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {

                if (A.array[i][j] < 1e-8) {
                    A.array[i][j] = 1e-8;
                }

//                A.array[i][j] = pow(A.array[i][j], (beta - 2));
//                A.array[i][j] *= X->matrix.array[i][j];
//
                A.array[i][j] = X->matrix.array[i][j] * 1 / A.array[i][j];

            }
        }

//        SaveMatrixToCSV("X_hadamard_A.csv", &A);

        Matrix X_hadamard_A_padded = CreateMatrix(A.rows, A.cols + T);

        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                X_hadamard_A_padded.array[i][j] = A.array[i][j];
            }
        }

//        SaveMatrixToCSV("X_hadamard_A_padded.csv", &X_hadamard_A_padded);


        Matrix num = CreateMatrix(r, A.cols);

        DestroyMatrix(&A);

        for (int t = 0; t < T; t++) {
            Matrix W_at_t = GetMatrixFromDictionary(dictionary, 0, t);
            Matrix transposed = Transpose(&W_at_t);

            DestroyMatrix(&W_at_t);

            Matrix X_hadamard_A_windowed = CreateMatrix(X_hadamard_A_padded.rows, ncol);

            for (int i = 0; i < X_hadamard_A_padded.rows; i++) {
                for (int j = 0; j < ncol; j++) {
                    X_hadamard_A_windowed.array[i][j] = X_hadamard_A_padded.array[i][j + t];
                }
            }

            Matrix multiplied = MatrixMultiply(&transposed, &X_hadamard_A_windowed);

            DestroyMatrix(&X_hadamard_A_windowed);
            DestroyMatrix(&transposed);

            for (int i = 0; i < multiplied.rows; i++) {
                for (int j = 0; j < multiplied.cols; j++) {
                    num.array[i][j] += multiplied.array[i][j];
                }
            }
        }

        DestroyMatrix(&X_hadamard_A_padded);


//        SaveMatrixToCSV("num.csv", &num);

        for (int row = 0; row < activations.rows; row++) {
            for (int c = 0; c < ncol - T; c++) {
                if (denom_all_col.array[row][c] < 1e-8) {
                    denom_all_col.array[row][c] = 1e-8;
                }
                activations.array[row][c] *= (num.array[row][c] / denom_all_col.array[row][c]);
            }
        }

//        SaveMatrixToCSV("H1.csv", &activations);

        for (int c = (int) (ncol - T); c < ncol; c++) {
            for (int row = 0; row < activations.rows; row++) {
                activations.array[row][c] *=
                        (num.array[row][c] / denoms_cropped_for_end.array[ncol - c - 1][row]);
            }
        }

//        SaveMatrixToCSV("H2.csv", &activations);

        DestroyMatrix(&num);
        DestroyMatrix(&conved);

        conved = ComputeConvolution(dictionary, &activations, T);

        obj = BetaDivergence(&X->matrix, &conved, beta);

        printf("ComputeActivations: Iteration %d, obj = %f\n", iteration, obj);

        if ((fabs(obj - obj_prev) / error_int) < maximum_error) {
            printf("ComputeActivations: Converged sufficiently\n");
//            break;
        }

        obj_prev = obj;
        iteration++;
    }

    DestroyMatrix(&denom_all_col);
    DestroyMatrix(&denoms_cropped_for_end);

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


void TestActivations() {
    Matrix X = CreateMatrix(88, 100);

    for (int i = 0; i < X.rows; i++) {
        for (int j = 0; j < X.cols; j++) {
            X.array[i][j] = i * j % 100;
        }
    }

    Matrix notes = TranscribeNotesFromActivations(&X, 0, 0.08);
    SaveMatrixToCSV("notes.csv", &notes);
}

Matrix TranscribeNotesFromActivations(const Matrix *activations, double threshold, double time_step) {
//    Matrix notes = CreateMatrix(500, 3);
    Matrix notes = CreateMatrix(500, 2);

    int note_count = 0;

    int midi_codebook[88];

    for (int i = 0; i < 88; i++) {
        midi_codebook[i] = i + 21;
    }

    bool note_presence = false;

    SaveMatrixToCSV("unsmoothed_activation.csv", activations);

    double current_pitch = 0;
    double current_onset = 0;
    int sliding_window_width = 10;

//    TODO: midi file creation

    Matrix smoothed_activations = CreateMatrix(activations->rows, activations->cols);

    for (int i = 0; i < activations->cols; i++) {
        if (i - sliding_window_width < 0 || i + sliding_window_width + 1 > activations->cols) {
            int d = floor(fmin(i, activations->cols - i));

            for (int j = 0; j < activations->rows; j++) {
                double sum = 0;
                int count = 0;
                for (int l = 0; l < activations->rows; l++) {
                    for (int k = i - d; k < i + d + 1; k++) {
                        if ((k >= 0 && k < activations->cols) && (l >= 0 && l < activations->rows)) {
                            sum += activations->array[l][k];
                            count++;
                        }
                    }
                }
                smoothed_activations.array[j][i] = sum / count;
            }
        } else {
            for (int j = 0; j < activations->rows; j++) {
                double sum = 0;
                int count = 0;
                for (int k = i - sliding_window_width; k < i + sliding_window_width + 1; k++) {
                    if (k >= 0 && k < activations->cols) {
                        sum += activations->array[j][k];
                        count++;
                    }
                }

                smoothed_activations.array[j][i] = sum / count;
            }
        }
    }

    for (int note_index = 0; note_index < smoothed_activations.rows; note_index++) {
        if (note_presence)
            note_presence = false;

        for (int time_index = 0; time_index < smoothed_activations.cols; time_index++) {
            bool minimal_sustain_condition = (
                    activations->array[note_index][time_index] - smoothed_activations.array[note_index][time_index] >
                    threshold);

            if(minimal_sustain_condition){
                if(!note_presence){
                    current_pitch = midi_codebook[note_index];
                    current_onset = time_index * time_step;
                    note_presence = true;
                }
            }
            else{
                if(note_presence){
                    double current_offset = time_index * time_step;
                    notes.array[note_count][0] = current_onset;
//                    notes.array[note_count][1] = current_offset;
//                    notes.array[note_count][2] = current_pitch;
                    notes.array[note_count][1] = current_pitch;
                    note_count++;
//                    TODO: writing notes to midi file

                    note_presence = false;
                }
            }
        }

    }

    return notes;
}
