#include "htk.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <string>
#include <limits>

using namespace std;

#define MATCH -3
#define INDEL 5
#define SUB 1

struct location {
    int i;
    int j;
};

struct align_data {
    float score;
    location back_pointer;
};

__global__ void align_kernel(char* seq1, char* seq2, align_data* align_matrix, int rows, int cols, int d, float INF) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = d - i;
    
    if (i < rows && j < cols && j >= 0) {
        float left = (i == 0) ? INF : align_matrix[(i - 1) * cols + j].score + INDEL;
        float top = (j == 0) ? INF : align_matrix[i * cols + (j - 1)].score + INDEL;
        float diagonal = (i == 0 || j == 0) ? INF : ((seq1[i] == seq2[j]) ? MATCH : SUB) + align_matrix[(i - 1) * cols + (j - 1)].score;

        align_data best;
        best.score = INF;

        if (left < best.score) best = {left, {i - 1, j}};
        if (top < best.score) best = {top, {i, j - 1}};
        if (diagonal < best.score) best = {diagonal, {i - 1, j - 1}};
        if (best.score == INF) best.score = 0;

        align_matrix[i * cols + j] = best;
    }
}

float align(string s1, string s2, bool banded, int align_length) {
    const float INF = numeric_limits<float>::infinity();

    string seq1 = "-" + s1.substr(0, align_length);
    string seq2 = "-" + s2.substr(0, align_length);

    int rows = seq2.length();
    int cols = seq1.length();

    align_data* align_matrix;
    cudaMallocManaged(&align_matrix, rows * cols * sizeof(align_data));

    char* d_seq1;
    char* d_seq2;

    //managed memory added
    cudaMallocManaged(&d_seq1, seq1.length() * sizeof(char));
    cudaMallocManaged(&d_seq2, seq2.length() * sizeof(char));

    cudaStream_t stream;
    cudaStreamCreate(&stream);
    // asyn added
    cudaMemcpyAsync(d_seq1, seq1.c_str(), seq1.length() * sizeof(char), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_seq2, seq2.c_str(), seq2.length() * sizeof(char), cudaMemcpyHostToDevice, stream);

    dim3 threadsPerBlock(16);
    for (int d = 0; d < rows + cols - 1; ++d) {
        int numElements = min(d + 1, min(rows, cols));
        dim3 numBlocks((numElements + threadsPerBlock.x - 1) / threadsPerBlock.x);
        align_kernel<<<numBlocks, threadsPerBlock, 0, stream>>>(d_seq1, d_seq2, align_matrix, rows, cols, d, INF);
    }

    cudaStreamSynchronize(stream);
    cudaStreamDestroy(stream);

    align_data result = align_matrix[(rows - 1) * cols + (cols - 1)];

    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(align_matrix);

    return result.score;
}

int main(int argc, char* argv[]) {
    string line;
    string seq[8];
    string filename[8] = {"BCoV-ENT", "BCoV-LUN", "BCoV-Mebus", "BCoV-Quebec",
                          "Mouse_Hepatitis", "Murine_Hepatitis1", "Murine_Hepatitis2", "Murine_Hepatitis3"};

    for (int i = 0; i < 8; i++) {
        ifstream file("./sequences/" + filename[i] + ".txt");
        while (getline(file, line)) seq[i] += line;
    }

    htkTime_start(Compute, "Doing the computation");

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if (j >= i) {
                auto result = align(seq[i], seq[j], false, 30000);
                printf("Result: %f\n", result);
            }
        }
    }

    htkTime_stop(Compute, "Doing the computation");

    return 0;
}
