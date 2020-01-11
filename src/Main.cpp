#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "NeedlemanWunsch.cpp"

#define GET_INDEX(i, j, m) ((i) * (m) + (j)) 

double emissionMatrix[5][5] = {
  { 0.0, 0.625, 0.375, 0.0, 0.0 },
  { 0.0, 0.3125, 0.0625, 0.5, 0.125 },
  { 0.0, 0.08333333, 0.25, 0.58333333, 0.08333333 },
  { 0.0, 0.25862069, 0.22413793, 0.48275862, 0.03448276 },
  { 0.0, 0.0, 0.0, 0.0, 0.0 }        
};
double transmissionMatrix[5][5] = {
  {0.,         0.21875,    0.46875,    0.3125    , 0.0},
  {0.29166667, 0.32758621, 0.,         0.        , 0.0},
  {0.45833333, 0.,         0.5,        0.        , 0.0},
  {0.25,       0.,         0.,         0.17241379, 0.0},
  {0.25,       0.,         0.,         0.17241379, 0.0}
};

double** allocateDouble2D(int n, int m) {
  double **matrix = new double*[n];
  for (int i = 0; i < n; i++) {
    matrix[i] = new double[m];
    for (int j = 0; j < m; j++) {
      matrix[i][j] = 0.0;
    }
  }

   
  return matrix;
}

char** allocateChar2D(int n, int m) {
  char **matrix = new char*[n];
  for (int i = 0; i < n; i++) {
    matrix[i] = new char[m];
  }
  return matrix;
}


void viterbi() {
  

  std::ifstream infile1("g1.txt");
  std::ifstream infile2("g2.txt");


  std::string firstGene;
  std::cout << "jurij" << std::endl;
  std::string line;
  while (std::getline(infile1, line)) {
    //std::istringstream iss(line);
    ///std::cout << a << " " << b << std::endl;
    // process pair (a,b)
    firstGene += line;
  }

  std::string secondGene;
  while (std::getline(infile2, line)) {
    //std::istringstream iss(line);
    ///std::cout << a << " " << b << std::endl;
    // process pair (a,b)
    secondGene += line;
  }

  std::cout << "1: " << firstGene << std::endl;
  std::cout << "2: " << secondGene << std::endl;

  int n = firstGene.size();
  int m = secondGene.size();

  
  double **withMM = allocateDouble2D(n + 1, m + 1);
  double **emitX = allocateDouble2D(n + 1, m + 1);
  double **emitY = allocateDouble2D(n + 1, m + 1);

  char **traceMM = allocateChar2D(n + 1, m + 1);
  char **traceX = allocateChar2D(n + 1, m + 1);
  char **traceY = allocateChar2D(n + 1, m + 1);

  // alociraj pocetno stanje 


  withMM[0][0] = 1; 


  
}

int main(void)
{
  viterbi();
  //new NeedlemanWunsch();

  return 0;
}