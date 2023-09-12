//
// Created by ruanb on 9/12/2023.
//

#ifndef C_PIANO_TRANSCRIPTION_MATRIX_H
#define C_PIANO_TRANSCRIPTION_MATRIX_H

typedef struct MatrixStruct{
    double **array;
    unsigned int rows;
    unsigned int cols;
} Matrix;

Matrix CreateMatrix(unsigned int nRows, unsigned int nCols);
void DestroyMatrix(Matrix *matrix);

Matrix MatrixMultiply(Matrix const *a, Matrix const *b);
Matrix SumMatricesAlongAxis(Matrix *a, unsigned int numMatrices, unsigned int axis);
Matrix SumMatrices(Matrix *a, unsigned int numMatrices);
Matrix Transpose(Matrix const *matrix);

void FillMatrix(Matrix* matrix, double value);
void NormaliseMatrix(Matrix *matrix);
Matrix ShiftMatrix(Matrix const *matrix, unsigned int numShifts);


#endif //C_PIANO_TRANSCRIPTION_MATRIX_H
