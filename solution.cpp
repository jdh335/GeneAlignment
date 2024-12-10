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

float align(string s1, string s2, bool banded, int align_length){

  const float INF = numeric_limits<float>::infinity();

  string seq1 = "-" + s1.substr(0, align_length);
  string seq2 = "-" + s2.substr(0, align_length);

  int rows = seq2.length();
  int cols = seq1.length();

  align_data** align_matrix = new align_data*[rows];
  for (int i = 0; i < rows; ++i) {
    align_matrix[i] = new align_data[cols];
  }

  for (int d = 0; d < rows + cols - 1; ++d) {
    #pragma omp parallel for
    for (int i = max(0, d - cols + 1); i <= min(d, rows - 1); ++i) {
      int j = d - i;

      float left = (i==0) ? INF : align_matrix[i-1][j].score + INDEL;
      float top = (j==0) ? INF : align_matrix[i][j-1].score + INDEL;
      float diagonal = (i==0||j==0) ? INF : ((seq1[i]==seq2[j]) ?  MATCH : SUB) + align_matrix[i-1][j-1].score;

      align_data best;
      best.score = INF;

      if (left < best.score) best = {left, {i-1, j}};
      if (top < best.score) best = {top, {i, j-1}};
      if (diagonal < best.score) best = {diagonal, {i-1, j-1}};
      if (best.score == INF) best.score = 0;

      align_matrix[i][j] = best;
    }
  }

  align_data result = align_matrix[cols-1][rows-1];

  for (int i = 0; i < rows; ++i) {
    delete[] align_matrix[i];
  }
  delete[] align_matrix;


  return result.score;
}


int main(int argc, char *argv[]) {

  string line;
  string seq[8];
  string filename[8] = {"BCoV-ENT", "BCoV-LUN", "BCoV-Mebus", "BCoV-Quebec",
                        "Mouse_Hepatitis", "Murine_Hepatitis1", "Murine_Hepatitis2", "Murine_Hepatitis3"};

  // Read sequence into string from file
  for (int i = 0; i < 8; i++) {
    ifstream file("./sequences/" + filename[i] + ".txt");
    while (getline(file, line)) seq[i] += line;
  }

  htkTime_start(Compute, "Doing the computation");

  // Compute upper triangular of comparison string
  // #pragma omp parallel for collapse(2)
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      if (j >= i) {
        auto result = align(seq[i], seq[j], false, 15000);
        printf("Result: %f\n", result);
      }
    }
  }

  htkTime_stop(Compute, "Doing the computation");

  return 0;
}
