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
