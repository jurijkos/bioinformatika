#ifndef HMMALIGN_H
#define HMMALIGN_H

#include <string>

class HMMAlign {
public:
  HMMAlign(std::string firstGene, std::string secondGene);
  void print();
  void run();
private:
  void viterbi_log();
  void backtrace();
  void dataPreprocessing();
  int getIndexOfBase(char c);
  double** allocateDouble2D(int n, int m);
  char** allocateChar2D(int n, int m);
  void proba_to_log(double p[][5]);
  void to_log(double **p, int n, int m);
  void setInFile1(std::string inFile);
  void setInFile2(std::string inFile);
  void printSolution();
  double transmissionMatrix[5][5] = {
    {0.        ,0.16324153,0.70367022,0.13308825,0.        },
    {0.        ,0.24061506,0.25119354,0.38537014,0.12282125},
    {0.        ,0.2485594 ,0.37534149,0.35442674,0.02167237},
    {0.        ,0.26130111,0.42213701,0.30572308,0.0108388 },
    {0.        ,0.        ,0.        ,0.        ,0.        }
  };


  double emissionMatrix[5][5] = {
    {0.        ,0.25993764,0.15561077,0.5255119 ,0.05893969},
    {0.23525022,0.14544823,0.02395958,0.02000951,0.01412214},
    {0.27924401,0.01898803,0.17671419,0.03698476,0.00154984},
    {0.0808479 ,0.00125295,0.00399033,0.29688447,0.00110199},
    {0.40465787,0.02299597,0.00523322,0.023001  ,0.20776378}
  };

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
};

#endif