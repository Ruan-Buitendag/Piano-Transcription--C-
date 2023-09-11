//
// Created by ruanb on 9/11/2023.
//
#include <stdio.h>
#include <stdlib.h>

// Function to right shift matrix H by t columns
// H: activation matrix H
// t: shift number
// Returns: matrix H after shift
int** shift(int** H, int rows, int cols, int t) {
    int** H_shift = (int**)malloc(rows * sizeof(int*));

    for (int i = 0; i < rows; i++) {
        H_shift[i] = (int*)malloc(cols * sizeof(int));
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (j + t < cols) {
                H_shift[i][j + t] = H[i][j];
            }
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < t; j++) {
            H_shift[i][j ] = 0;

        }
    }

    return H_shift;
}

int main() {
    int rows = 5;
    int cols = 5;

    // Initialize matrix H with example values
    int** H = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        H[i] = (int*)malloc(cols * sizeof(int));
    }

    // Populate H with example values (you can replace these)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            H[i][j] = i * cols + j;
        }
    }

    int t = 2; // Shift by 2 columns

    // Call the shift function
    int** H_shifted = shift(H, rows, cols, t);

    // Print the shifted matrix
    printf("Shifted matrix H:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%d ", H_shifted[i][j]);
        }
        printf("\n");
    }

    // Free memory
    for (int i = 0; i < rows; i++) {
        free(H[i]);
        free(H_shifted[i]);
    }
    free(H);
    free(H_shifted);

    return 0;
}
