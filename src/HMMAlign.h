#ifndef HMMALIGN_H
#define HMMALIGN_H

#include <string>

class HMMAlign {
public:
  HMMAlign(std::string firstGene, std::string secondGene, double transmissionMatrix[5][5], double emissionMatrix[5][5]);
  ~HMMAlign();
  void print();
  void run();
private:
  void viterbi_log();
  void backtrace();
  void dataPreprocessing();
  int getIndexOfBase(char c);
  double** allocateDouble2D(int n, int m);
  void emptyDouble2D(double **matrix, int n);
  char** allocateChar2D(int n, int m);
  void emptyChar2D(char **matrix, int n);
  void proba_to_log(double p[][5]);
  void to_log(double **p, int n, int m);
  void setInFile1(std::string inFile);
  void setInFile2(std::string inFile);
  void printSolution();
  double transmissionMatrix[5][5];
  double emissionMatrix[5][5];

  std::string firstGene;
  std::string secondGene;
  std::string alignedX;
  std::string alignedY;
  int n;
  int m; 
  double **withMM;
  double **withEmitX;
  double **withEmitY;
  char **traceMM;
  char **traceX;
  char **traceY;
  const static int FROM_MM =  1;
  const static int FROM_EX =  2;
  const static int FROM_EY =  3;
};

#endif