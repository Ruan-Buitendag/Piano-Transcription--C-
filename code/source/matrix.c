//
// Created by ruanb on 9/12/2023.
//

#include "matrix.h"
#include "stdlib.h"

Matrix CreateMatrix(unsigned int nRows, unsigned int nCols) {
    Matrix matrix;
    matrix.array = (double **) malloc(nRows * sizeof(double *));

    for (int i = 0; i < nRows; i++) {
        matrix.array[i] = (double *) calloc(nCols, sizeof(double));
    }

    matrix.rows = nRows;
    matrix.cols = nCols;

    return matrix;
}

void DestroyMatrix(Matrix *matrix) {
    if(matrix->array == NULL){
        return;
    }

    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->array[i]);
    }

    free(matrix->array);
    matrix->array = NULL;
}

Matrix MatrixMultiply(const Matrix *a, const Matrix *b) {
    Matrix result = CreateMatrix(a->rows, b->cols);

    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            for (int k = 0; k < a->cols; k++) {
                result.array[i][j] += a->array[i][k] * b->array[k][j];
            }
        }
    }

    return result;
}

Matrix SumMatrices(Matrix *a, unsigned int numMatrices) {
    Matrix result = CreateMatrix(a[0].rows, a[0].cols);

    for (int i = 0; i < numMatrices; i++) {
        for (int j = 0; j < a[i].rows; j++) {
            for (int k = 0; k < a[i].cols; k++) {
                result.array[j][k] += a[i].array[j][k];
            }
        }
    }

    return result;

}

Matrix SumMatricesAlongAxis(Matrix *a, unsigned int numMatrices, unsigned int axis) {
    // TODO: optimisation is probably possible here
    Matrix result = CreateMatrix(a[0].rows, a[0].cols);

    for (int r = 0; r < result.rows; r++) {
        for (int c = 0; c < result.cols; c++) {
            for (int i = 0; i < numMatrices; i++) {
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

Matrix Transpose(const Matrix *matrix) {
    Matrix transposed = CreateMatrix(matrix->cols, matrix->rows);

    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            transposed.array[j][i] = matrix->array[i][j];
        }
    }

    return transposed;
}

void FillMatrix(Matrix *matrix, double value) {
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            matrix->array[i][j] = value;
        }
    }
}

void NormaliseMatrix(Matrix *matrix) {
    double max = 0;
    for (int r = 0; r < matrix->rows; r++) {
        for (int c = 0; c < matrix->cols; c++) {
            if (matrix->array[r][c] > max)
                max = matrix->array[r][c];
        }
    }

    for (int r = 0; r < matrix->rows; r++) {
        for (int c = 0; c < matrix->cols; c++) {
            matrix->array[r][c] /= max;
        }
    }
}

Matrix ShiftMatrix(const Matrix *matrix, unsigned int numShifts) {
    Matrix shifted = CreateMatrix(matrix->rows, matrix->cols);

    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            if (j + numShifts < matrix->cols) {
                shifted.array[i][j + numShifts] = matrix->array[i][j];
            }
        }
    }

    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < numShifts; j++) {
            shifted.array[i][j] = 0;

        }
    }

    return shifted;
}


