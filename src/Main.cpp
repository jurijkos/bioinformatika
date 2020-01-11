#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
//#include "NeedlemanWunsch.h"

#define GET_INDEX(i, j, m) ((i) * (m) + (j)) 
#define FROM_MM 1
#define FROM_EX 2
#define FROM_EY 3

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

int getIndexOfBase(char c) {
  switch (c) {
    case '-':
      return 0;
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 3;
    case 'T':
      return 4;  
  }
  return -1;
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
  double **withEmitX = allocateDouble2D(n + 1, m + 1);
  double **withEmitY = allocateDouble2D(n + 1, m + 1);

  char **traceMM = allocateChar2D(n + 1, m + 1);
  char **traceX = allocateChar2D(n + 1, m + 1);
  char **traceY = allocateChar2D(n + 1, m + 1);

  // alociraj pocetno stanje 


  withMM[0][0] = 1; 

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
        char x = firstGene[i - 1];
        char y = secondGene[j - 1];
        double pxy = emissionMatrix[getIndexOfBase(x)][getIndexOfBase(y)];
        double qx = emissionMatrix[getIndexOfBase(x)][0];
        double qy = emissionMatrix[0][getIndexOfBase(y)];
        std::cout << pxy + qx + qy << std::endl;
        
        double wasMM = transmissionMatrix[3][3] * withMM[i - 1][j - 1];
        double wasEmitX = transmissionMatrix[2][3] * withEmitX[i - 1][j - 1];
        double wasEmitY = transmissionMatrix[1][3] * withEmitY[i - 1][j - 1];

        if (wasMM > wasEmitX && wasMM > wasEmitY) {
          withMM[i][j] = pxy *wasMM;
          traceMM[i][j] =  FROM_MM; 
        } else if (wasEmitX > wasMM && wasEmitX > wasEmitY) {
          withMM[i][j] = pxy * wasEmitX;
          traceMM[i][j] =  FROM_EX;
        } else {
          withMM[i][j] = pxy * wasEmitY;
          traceMM[i][j] = FROM_EY;
        }

        //emit in x
        wasMM = transmissionMatrix[3][2] * withMM[i - 1][j];
        wasEmitX = transmissionMatrix[2][2] * withEmitX[i - 1][j];
        
        if (wasMM > wasEmitX) {
          withEmitX[i][j] = qx *wasMM;
          traceX[i][j] =  FROM_MM; 
        } else {
          withEmitX[i][j] = qx * wasEmitX;
          traceX[i][j] = FROM_EX;
        }

        
        //emit in y
        wasMM = transmissionMatrix[3][1] * withMM[i][j - 1];
        wasEmitY = transmissionMatrix[1][1] * withEmitY[i][j - 1];
        
        if (wasMM > wasEmitY) {
          withEmitY[i][j] = qy * wasMM;
          traceY[i][j] =  FROM_MM; 
        } else {
          withEmitY[i][j] = qy * wasEmitY;
          traceY[i][j] = FROM_EY;
        }
    }
  }
  
}

int main(void)
{
  //viterbi();
  //new NeedlemanWunsch();
  std::cout << std::max({1, 3, 4}) << std::endl;

  return 0;
}