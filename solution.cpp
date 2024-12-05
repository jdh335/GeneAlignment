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

const double INF = numeric_limits<double>::infinity();

class GeneSequence {
public:

  string index_str(int i, int j) {
    return ("(" + to_string(i) + "," + to_string(j) + ")").c_str();
  }

  string seq1, seq2;
  map<string, tuple<double, int, int>> score_matrix;

  tuple<double, int, int> align(string s1, string s2, bool banded, int align_length){

    map<string, tuple<double, int, int>> new_matrix;
    score_matrix = new_matrix;
    
    seq1 = "-" + s1.substr(0, align_length);
    seq2 = "-" + s2.substr(0, align_length);

    int rows = seq2.length();
    int cols = seq1.length();

    if (!banded) {
      for (int d = 0; d < rows + cols - 1; ++d) {
        #pragma omp parallel for schedule(dynamic,4)
        for (int i = max(0, d - cols + 1); i <= min(d, rows - 1); ++i) {
          int j = d - i;
          tuple<double, int, int> score = compute_score(i, j);
          #pragma omp critical
          {
            score_matrix.insert({index_str(i,j), score});
          }
        }
      }
    } else {
      for (int j = 0; j < seq2.length(); j++) {
        for (int i = -3; i < 4; i++) {
          if ((i + j >= 0) && (i + j < seq1.length())) {
            tuple<double, int, int> score = compute_score(i+j, j);
            score_matrix.insert({index_str(i+j,j), score});
          }
        }
      }
    }

    // for (int i = 0; i < seq1.length(); i++) {
    //   for (int j = 0; j < seq2.length(); j++) {
    //     score_matrix.insert({index_str(i,j), compute_score(i, j)});
    //   }
    // }

    auto result = score_matrix[index_str(seq1.length()-1,seq2.length()-1)];

    return {get<0>(result), get<1>(result), get<2>(result)};
  }

  tuple<double, int, int> compute_score(int i, int j) {
    double score = INF;
    tuple<int, int> back_pointer;

    double left = getLeft(i,j);
    double top = getTop(i,j);
    double diagonal = getDiagonal(i,j);

    if (left < score){
      score = left;
      back_pointer = {i-1, j};
    }

    if (top < score){
      score = top;
      back_pointer = {i, j-1};
    }

    if (diagonal < score){
      score = diagonal;
      back_pointer = {i-1, j-1};
    }

    if (score == INF) score = 0;

    return {score, get<0>(back_pointer), get<1>(back_pointer)};
  }

  double getLeft(int i, int j){
    if (score_matrix.find(index_str(i-1, j)) != score_matrix.end())
      return get<0>(score_matrix[index_str(i-1, j)]) + INDEL;
    else
      return INF;
  }

  double getTop(int i, int j){
    if (score_matrix.find(index_str(i, j-1)) != score_matrix.end())
      return get<0>(score_matrix[index_str(i, j-1)]) + INDEL;
    else
      return INF;
  }

  double getDiagonal(int i, int j){
    if (score_matrix.find(index_str(i-1, j-1)) != score_matrix.end()) {
      double diagonal = get<0>(score_matrix[index_str(i-1, j-1)]);
      if (seq1[i] == seq2[j]) return diagonal + MATCH;
      else return diagonal + SUB;
    }
    else
      return INF;
  }

};

int main(int argc, char *argv[]) {

  string line;
  string seq[8];
  string filename[8] = {"BCoV-ENT", "BCoV-LUN", "BCoV-Mebus", "BCoV-Quebec","Mouse_Hepatitis",
                         "Murine_Hepatitis1", "Murine_Hepatitis2", "Murine_Hepatitis3"};

  // Read sequence into string from file
  for (int i = 0; i < 8; i++) {
    ifstream file("./sequences/" + filename[i] + ".txt");
    while (getline(file, line)) seq[i] += line;
  }

  htkTime_start(Compute, "Doing the computation");

  // Compute upper triangular of comparison string
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      if (j >= i) {
        GeneSequence gs;
        auto result = gs.align(seq[i], seq[j], false, 1000);
        printf("Result: %f\n", get<0>(result));
      }
    }
  }

  htkTime_stop(Compute, "Doing the computation");

  return 0;
}
