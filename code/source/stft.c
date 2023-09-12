//
// Created by ruanb on 9/5/2023.
//

#include "stft.h"
#include "math.h"
#include "spectrogram.h"

#include <stdio.h>
#include <complex.h>
#include "stdlib.h"

Spectrogram STFT(DynamicArray *x, int windowSize, int hopSize, int fftSize) {
//    TODO: redo padding sodat dit minder stupid is noudat jy calloc gebruik

    DynamicArray xCopy = CreateDynamicArray(x->size + 2 * windowSize);

    CopyDynamicArray(&xCopy, x);

    DynamicArray padding = CreateDynamicArray(windowSize / 2);
    FillDynamicArray(&padding, 0);

    DynamicArray beginningPaddedSignal = AppendDynamicArray(&padding, &xCopy);

    DynamicArray paddedSignal = AppendDynamicArray(&beginningPaddedSignal, &padding);

    DestroyDynamicArray(&padding);
    DestroyDynamicArray(&beginningPaddedSignal);

    int nPad = 0;

    while ((paddedSignal.size + nPad - windowSize) % hopSize != 0) {
        nPad++;
    }

    DynamicArray paddingForWindowing = CreateDynamicArray(nPad);
    FillDynamicArray(&paddingForWindowing, 0);

    DynamicArray finalPaddedSignal = AppendDynamicArray(&paddedSignal, &paddingForWindowing);

    DestroyDynamicArray(&paddedSignal);
    DestroyDynamicArray(&paddingForWindowing);

    DynamicArray window = HanningWindow(windowSize);

    double sum = Sum(&window);

    unsigned int N = windowSize;
    unsigned int L = finalPaddedSignal.size;
    unsigned int M = ((L - N) / hopSize) + 1;

    DestroyDynamicArray(&xCopy);

    Spectrogram X = CreateSpectrogram(fftSize / 2 + 1, M);

    for (int i = 0; i < M; i++) {
        DynamicArray xWindowed = CreateDynamicArray(fftSize);

        for (int j = 0; j < windowSize; j++) {
            xWindowed.array[j] = finalPaddedSignal.array[i * hopSize + j] * window.array[j];
        }

        for(int k = windowSize; k < fftSize; k++) {
            xWindowed.array[k] = 0;
        }


        DynamicArray XWindowed = CreateDynamicArray(fftSize / 2);

        rfft(xWindowed.array, fftSize, XWindowed.array);

        for (int j = 0; j < XWindowed.size ; j++) {
            X.matrix.array[j][i] = fabs(XWindowed.array[j]) / sum;
        }

        //    print the contents of xWindowed
        DestroyDynamicArray(&xWindowed);
        DestroyDynamicArray(&XWindowed);
    }

    SaveSpectrogramToCSV("spectrogram.csv", &X);

    return X;
}

DynamicArray HanningWindow(int windowSize) {
    DynamicArray window = CreateDynamicArray(windowSize);

    for (int i = 0; i < windowSize; i++) {
        window.array[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (windowSize - 1)));
    }

    return window;
}


void fft(const double *x, int N, double complex *result) {
    if (N <= 1) {
        result[0] = x[0];
        return;
    }

    // Create arrays to store the even and odd parts of the signal
    double *even = (double *) malloc((N / 2) * sizeof(double));
    double *odd = (double *) malloc((N / 2) * sizeof(double));

    if (even == NULL || odd == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < N / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    // Recursive FFT on even and odd parts
    double complex *even_result = (double complex *) malloc((N / 2) * sizeof(double complex));
    double complex *odd_result = (double complex *) malloc((N / 2) * sizeof(double complex));

    if (even_result == NULL || odd_result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    fft(even, N / 2, even_result);
    fft(odd, N / 2, odd_result);

    // Combine the results
    for (int k = 0; k < N / 2; k++) {
        double complex t = cexp(-I * 2.0 * M_PI * k / N) * odd_result[k];
        result[k] = even_result[k] + t;
        result[k + N / 2] = even_result[k] - t;
    }

    // Clean up allocated memory
    free(even);
    free(odd);
    free(even_result);
    free(odd_result);
}

void rfft(const double *x, int N, double *result) {
    complex double complex_result[N];

    fft(x, N, complex_result);

    for (int i = 0; i < N / 2; i++) {
        result[i] = cabs(complex_result[i]);
    }
}

void stftTest() {

    int sampleRate = 5000;  // Sample rate in Hz
    double frequency = 440.0;  // Frequency of the sawtooth wave in Hz
    double amplitude = 1;  // Amplitude of the wave

    int duration = 1;  // Duration of the waveform in seconds
    int numSamples = sampleRate * duration;

    double signal[numSamples];

    for (int i = 0; i < numSamples; i++) {
        double t = (double) i / sampleRate;
        double value = 2.0 * (t * frequency - floor(0.5 + t * frequency)); // Sawtooth wave formula

        // Scale the value by the amplitude
        value *= amplitude;

        // Output the sample (assuming 16-bit signed PCM)
        double sample = value;

        signal[i] = sample;
    }

    double fftr[1024];

    rfft(signal, 2048, fftr);

    for (int i = 0; i < 1024; i++) {
        printf("%f\n", fftr[i]);
    }
}







