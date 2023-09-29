//
// Created by ruanb on 9/29/2023.
//

#include <math.h>
#include "evaluation.h"

Graph CreateGraph(int maxRef, int maxEstCount) {
//    Matrix result = CreateMatrix(maxRef, maxEstCount);
    Graph result;

    result.array = (IntArray *) calloc(maxRef, sizeof(IntArray));

    for (int i = 0; i < maxRef; i++) {
        result.array[i] = CreateIntArray(maxEstCount);
    }

    result.left = maxRef;
    result.right = maxEstCount;

    return result;
}

void DestroyGraph(Graph *graph) {
    for (int i = 0; i < graph->left; i++) {
        DestroyIntArray(&graph->array[i]);
    }

    free(graph->array);
    graph->array = NULL;
}


void AddEdge(Graph *graph, int ref, int est) {
    int currIndex = graph->array[ref].size;

    graph->array[ref].array[currIndex] = est;

    graph->array[ref].size += 1;
}

IntArray CreateIntArray(int capacity) {
    IntArray result;

    result.array = (int *) calloc(capacity, sizeof(int));
    result.size = 0;
    result.capacity = capacity;

    return result;
}


void DestroyIntArray(IntArray *array) {
    free(array->array);
    array->array = NULL;
}

bool dfs(Graph *graph, int u, IntArray *visited, IntArray *matchR, IntArray *matchL) {
    for (int v_index = 0; graph->array[u].size; v_index++) {
        int v = graph->array[u].array[v_index];

        if (visited->array[v]) {
            continue;
        }

        visited->array[v] = 1;

        if (matchR->array[v] == -1 || dfs(graph, matchR->array[v], visited, matchR, matchL)) {
            matchL->array[u] = v;
            matchR->array[v] = u;

            return true;
        }
    }

    return false;
}

int MaxBipartiteMatching(Graph *graph) {
    IntArray leftMatch = CreateIntArray(graph->left);
    IntArray rightMatch = CreateIntArray(graph->right);

    for (int i = 0; i < graph->left; i++) {
        leftMatch.array[i] = -1;
    }
    for (int i = 0; i < graph->right; i++) {
        rightMatch.array[i] = -1;
    }

    int maxMatching = 0;

    for (int u = 0; u < graph->left; u++) {
        IntArray visited = CreateIntArray(graph->left);

        if (dfs(graph, u, &visited, &rightMatch, &leftMatch)) {
            maxMatching++;
        }

        DestroyIntArray(&visited);
    }

    DestroyIntArray(&leftMatch);
    DestroyIntArray(&rightMatch);

    return maxMatching;
}

void GraphTest() {
    Graph graph = CreateGraph(1000, 1000);

    int ref_i_test[] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16,
                        17, 18, 20, 21, 22, 25, 26, 27, 28, 30, 32, 33, 34,
                        35, 36, 37, 41, 42, 43, 44, 51, 52, 53, 54, 55, 57,
                        58, 59, 60, 61, 62, 63, 64, 66, 68, 69, 70, 71, 73,
                        74, 75, 76, 78, 79, 81, 82, 83, 86, 88, 90, 91, 92,
                        93, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 111,
                        112, 114, 115, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
                        129, 130, 132, 133, 134, 135, 136, 138, 139, 140, 141, 142, 144,
                        145, 147, 148, 149, 150, 151, 152, 153, 154, 156, 157, 158, 160,
                        161, 162, 163, 164, 167, 169, 171, 172, 174, 176, 177, 178, 179,
                        181, 183, 184, 185, 186, 188, 190, 192, 193, 195, 196, 198, 199,
                        200, 201, 203, 204, 205, 206, 207, 210, 213, 214, 215, 216, 218,
                        219, 221, 222, 223, 224, 225, 227, 228, 229, 230, 232, 232, 233,
                        234, 235, 236, 237, 238, 239, 240, 241, 242, 244, 246, 247, 248,
                        249, 250, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                        264, 265, 266, 268, 269, 270, 271, 272, 273, 274, 277, 278, 279,
                        280, 281, 282, 282, 283, 284, 285, 286, 287, 288, 289, 290, 292,
                        293, 294, 296, 296, 298, 299, 301, 302, 303, 304, 305, 306, 306,
                        307, 309, 309, 311, 314, 316, 317, 319, 321, 322, 323, 324, 325,
                        326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 337, 338, 339,
                        340, 341, 342, 343, 344, 345, 346, 348, 350, 351, 353, 354, 355,
                        356, 357, 358, 360, 361, 364, 365, 366, 367, 371, 372, 373, 374,
                        377, 378, 383, 384, 386, 387, 388, 389, 390, 391, 392};

    int est_i_test[] = {164, 198, 51, 99, 100, 165, 199, 152, 147, 200, 135, 127, 101,
                        153, 220, 69, 154, 221, 102, 103, 166, 201, 52, 167, 202, 155,
                        148, 203, 136, 128, 104, 156, 222, 105, 168, 204, 243, 293, 53,
                        3, 28, 4, 5, 157, 223, 289, 70, 244, 278, 169, 205, 80,
                        279, 6, 290, 280, 245, 7, 90, 224, 184, 71, 185, 206, 106,
                        170, 225, 246, 8, 29, 247, 9, 81, 10, 107, 91, 49, 108,
                        82, 54, 83, 11, 137, 30, 109, 12, 138, 171, 13, 139, 92,
                        72, 140, 110, 172, 207, 248, 294, 55, 14, 31, 15, 16, 226,
                        291, 73, 17, 250, 281, 173, 208, 0, 84, 18, 292, 32, 251,
                        19, 20, 93, 227, 186, 74, 187, 209, 174, 56, 228, 252, 21,
                        253, 85, 22, 111, 94, 50, 95, 112, 86, 57, 87, 23, 141,
                        33, 113, 142, 175, 24, 143, 96, 75, 114, 176, 254, 295, 79,
                        25, 26, 27, 149, 188, 282, 129, 34, 229, 255, 193, 194, 1,
                        115, 256, 35, 268, 58, 283, 36, 37, 269, 189, 130, 38, 230,
                        257, 195, 116, 270, 258, 39, 271, 59, 284, 40, 296, 300, 41,
                        297, 190, 285, 131, 272, 259, 42, 231, 260, 196, 117, 273, 261,
                        43, 274, 60, 61, 286, 44, 298, 301, 45, 299, 191, 287, 132,
                        275, 262, 232, 233, 197, 263, 118, 46, 62, 47, 48, 234, 264,
                        119, 234, 264, 76, 235, 177, 210, 63, 68, 64, 236, 65, 265,
                        178, 237, 276, 2, 66, 288, 77, 277, 179, 266, 88, 211, 120,
                        238, 144, 212, 180, 267, 121, 213, 239, 158, 214, 97, 192, 98,
                        215, 122, 181, 89, 123, 159, 150, 216, 145, 133, 124, 160, 240,
                        161, 241, 182, 217, 67, 125, 183, 218, 162, 151, 219};


    for (int i = 0; i < 297; i++) {
        AddEdge(&graph, ref_i_test[i], est_i_test[i]);
    }

    int maxMatching = MaxBipartiteMatching(&graph);

    printf("Max matching: %d\n", maxMatching);

    DestroyGraph(&graph);

}

// TODO: optimize for loops
void EvaluateTranscription(Matrix *ref, Matrix *est) {
    Matrix onset_distances = CreateMatrix(ref->rows, est->rows);

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            onset_distances.array[i][j] = round(fabs(ref->array[i][0] - est->array[j][0]) * 10000) / 10000;
        }
    }

    Matrix onset_hit_matrix = CreateMatrix(ref->rows, est->rows);

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            onset_hit_matrix.array[i][j] = onset_distances.array[i][j] <= 0.05;
        }
    }

    DestroyMatrix(&onset_distances);

    Matrix pitch_distances = CreateMatrix(ref->rows, est->rows);

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            pitch_distances.array[i][j] = fabs(1200 * (log2(ref->array[i][1]) - log2(est->array[j][1])));
        }
    }

    Matrix pitch_hit_matrix = CreateMatrix(ref->rows, est->rows);

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            pitch_hit_matrix.array[i][j] = pitch_distances.array[i][j] <= 50;
        }
    }

    DestroyMatrix(&pitch_distances);

    Matrix hits = CreateMatrix(ref->rows, est->rows);

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            hits.array[i][j] = onset_hit_matrix.array[i][j] * pitch_hit_matrix.array[i][j];
        }
    }

    DestroyMatrix(&onset_hit_matrix);
    DestroyMatrix(&pitch_hit_matrix);

    int ref_i[1000];
    int est_i[1000];
    int counter = 0;

    for (int i = 0; i < ref->rows; i++) {
        for (int j = 0; j < est->rows; j++) {
            if (hits.array[i][j]) {
                ref_i[counter] = i;
                est_i[counter] = j;
                counter++;
            }
        }
    }

    DestroyMatrix(&hits);

    Graph graph = CreateGraph(ref->rows, est->rows);

    for (int i = 0; i < counter; i++) {
        AddEdge(&graph, ref_i[i], est_i[i]);
    }

    int maxMatching = MaxBipartiteMatching(&graph);

//    printf("Max matching: %d\n", maxMatching);

    double precision = maxMatching / est->rows;
    double recall = maxMatching / ref->rows;
    double f1 = 2 * precision * recall / (precision + recall);

    printf("Precision: %f\n", precision);
    printf("Recall: %f\n", recall);
    printf("F1: %f\n", f1);

    DestroyGraph(&graph);

}


