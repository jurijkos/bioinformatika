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
/*
double transmissionMatrix[5][5] = {
  { 0.0, 0.625, 0.375,   0.0,   0.0 },
  { 0.0, 0.312, 0.0625,  0.5,   0.125 },
  { 0.0, 0.083, 0.25,    0.583, 0.083 },
  { 0.0, 0.25,  0.224,   0.48,  0.034 },
  { 0.0, 0.0,   0.0,     0.0,   0.0 }        
};
double emissionMatrix[5][5] = {
  {0.,         0.21875,    0.46875,    0.3125    , 0.0},
  {0.29166667, 0.32758621, 0.,         0.        , 0.0},
  {0.45833333, 0.,         0.5,        0.        , 0.0},
  {0.25,       0.,         0.,         0.17241379, 0.0},
  {0.25,       0.,         0.,         0.17241379, 0.0}
};
*/

double transmissionMatrix[5][5] = {
  { 0.0, 1 /3.0,   2/3.0,  0.0,   0.0 },
  { 0.0, 2/8.0,   1/8.0,  4/8.0,   1/8.0 },
  { 0.0, 1/7.0,   3/7.0,  2/7.0,   1/7.0},
  { 0.0, 3/10.0,  2/10.0, 4/10.0,  1/10.0 },
  { 0.0, 0.0,   0.0,  0.0,   0.0 }        
};
double emissionMatrix[5][5] = {
  {0.0,  2/7.0,  2/7.0,  1/7.0,   2/7.0},
  {1/9.0, 5/32.0, 1/32.0, 1/32.0,  1/32.0},
  {3/9.0, 1/32.0, 5/32.0, 1/32.0,  1/32.0},
  {2/9.0, 1/32.0, 1/32.0, 5/32.0,  1/32.0},
  {3/9.0, 1/32.0, 1/32.0, 1/32.0,  5/32.0}
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

void pm(double **matrix, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

void printRevers(std::string &string) {
  for (int i = string.size() - 1; i >=0; i--) {
    std::cout << string[i];
  }
  std::cout << "\n";
}

void viterbi() {
  

  std::ifstream infile1("g1.txt");
  std::ifstream infile2("g2.txt");


  std::string firstGene;
  std::cout << "jurij" << std::endl;
  std::string line;
  while (std::getline(infile1, line)) {
    firstGene += line;
  }

  std::string secondGene;
  while (std::getline(infile2, line)) {
    secondGene += line;
  }

  std::cout << "1: " << firstGene << std::endl;
  std::cout << "2: " << secondGene << std::endl;
  
  int n = firstGene.size();
  int m = secondGene.size();
  std::cout << "n: " << n << " m: " << m << std::endl;
  
  double **withMM = allocateDouble2D(n + 1, m + 1);
  double **withEmitX = allocateDouble2D(n + 1, m + 1);
  double **withEmitY = allocateDouble2D(n + 1, m + 1);

  char **traceMM = allocateChar2D(n + 1, m + 1);
  char **traceX = allocateChar2D(n + 1, m + 1);
  char **traceY = allocateChar2D(n + 1, m + 1);
  
  withMM[0][0] = 1;
  // alociraj pocetno stanje 
  //pm(withMM, n + 1, m + 1);
  //pm(withEmitX, n + 1, m + 1);
  //pm(withEmitY, n + 1, m + 1);

  
  for (int i = 1; i <= 1; i++) {
    for (int j = 1; j <= m; j++) {
        char x = firstGene[i - 1];
        char y = secondGene[j - 1];
        double pxy = emissionMatrix[getIndexOfBase(x)][getIndexOfBase(y)];
        double qx = emissionMatrix[getIndexOfBase(x)][0];
        double qy = emissionMatrix[0][getIndexOfBase(y)];
        
        double wasMM = transmissionMatrix[3][3] * withMM[i - 1][j - 1];
        double wasEmitX = transmissionMatrix[2][3] * withEmitX[i - 1][j - 1];
        double wasEmitY = transmissionMatrix[1][3] * withEmitY[i - 1][j - 1];
        std::cout << "was: " << wasMM << " " << wasEmitX << " " << wasEmitY << std::endl; 
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
    pm(withMM, n + 1, m + 1);
    pm(withEmitX, n + 1, m + 1);
    pm(withEmitY, n + 1, m + 1);
  }
  
/*
  // backtrace
  std::string alignedX;
  std::string alignedY;
  
  char currentState;
  double wasMM = transmissionMatrix[3][4] * withMM[n][m];
  double wasEmitX = transmissionMatrix[2][4] * withEmitX[n][m];
  double wasEmitY = transmissionMatrix[1][4] * withEmitY[n][m];

  if (wasMM > wasEmitX && wasMM > wasEmitY) {
    currentState = FROM_MM;
  } else if (wasEmitX > wasMM && wasEmitX > wasEmitY) {
    currentState = FROM_EX;
  } else {
    currentState = FROM_EY;
  }

  int i = n, j = m;
  while (true) {
    
    if (currentState == FROM_MM) {
      alignedX += firstGene[i];
      alignedY += secondGene[j];
      currentState = traceMM[i][j];
      i--; j--;
    } else if (currentState == FROM_EX) {
      alignedX += firstGene[i];
      alignedY += '-';
      currentState = traceX[i][j];
      i--;
    } else if (currentState == FROM_EY) {
      alignedX += '-';
      alignedY += secondGene[j];
      currentState = traceY[i][j];
      j--;
    } else {
      std::cout << "ERROR: BACK TRACE" << std::endl;
    }

    if (i <= 0 && j <= 0) {
      break;
    }
  }
  */
  
 
}

int main(void)
{
  viterbi();
  //new NeedlemanWunsch();
  //std::cout << std::max({1, 3, 4}) << std::endl;

  return 0;
}


