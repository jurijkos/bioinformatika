#include <iostream>
#include "HMMAlign.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <bits/stdc++.h> 
#include <cmath>

#define FROM_MM 1
#define FROM_EX 2
#define FROM_EY 3

void HMMAlign::print() {
    std::cout << "hello" << std::endl;
}



double** HMMAlign::allocateDouble2D(int n, int m) {
  double **matrix = new double*[n];
  for (int i = 0; i < n; i++) {
    matrix[i] = new double[m];
    for (int j = 0; j < m; j++) {
      matrix[i][j] = 0.0;
    }
  }

   
  return matrix;
}

char** HMMAlign::allocateChar2D(int n, int m) {
  char **matrix = new char*[n];
  for (int i = 0; i < n; i++) {
    matrix[i] = new char[m];
  }
  return matrix;
}

int HMMAlign::getIndexOfBase(char c) {
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

void HMMAlign::pm(double **matrix, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%9.6f", matrix[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

void HMMAlign::ptm(char **matrix, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%3d", matrix[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

//ispisuje matrice dinamike i trace zajedno
void HMMAlign::pbm(double **mat, char **t, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%9.6f(%d)", mat[i][j], t[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

void HMMAlign::proba_to_log(double p[][5]) {
	for (int i = 0; i < 5; ++i) {
		for(int j = 0; j < 5; ++j) {
			if (fabs(p[i][j]) < 10e-5) p[i][j] = -INFINITY;
			else p[i][j] = std::log(p[i][j]);
		}
	}
}

void HMMAlign::to_log(double **p, int n, int m) {
	for (int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			if (fabs(p[i][j]) < 10e-5) p[i][j] = -INFINITY;
			else p[i][j] = std::log(p[i][j]);
		}
	}
}

void HMMAlign::viterbi_log(char g1FileName[], char g2FileName[]) {
  std::ifstream infile1(g1FileName);
  std::ifstream infile2(g2FileName);

  std::string firstGene;
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
  
  
  // alociraj pocetno stanje 
  
  proba_to_log(transmissionMatrix);
  proba_to_log(emissionMatrix);
  
  //to_log() se moze hardkodirati u to -INFINITY jer su svi 0.0
  to_log(withMM, n+1, m+1);
  to_log(withEmitX, n+1, m+1);
  to_log(withEmitY, n+1, m+1);
  withMM[0][0] = 1;
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
        char x = firstGene[i - 1];
        char y = secondGene[j - 1];
        
        double pxy = emissionMatrix[getIndexOfBase(x)][getIndexOfBase(y)];
        double qx = emissionMatrix[getIndexOfBase(x)][0];
        double qy = emissionMatrix[0][getIndexOfBase(y)];
        
        double wasMM = transmissionMatrix[3][3] + withMM[i - 1][j - 1];
        double wasEmitX = transmissionMatrix[2][3] + withEmitX[i - 1][j - 1];
        double wasEmitY = transmissionMatrix[1][3] + withEmitY[i - 1][j - 1];
        
        if (wasMM > wasEmitX && wasMM > wasEmitY) {
          withMM[i][j] = pxy + wasMM;
          traceMM[i][j] =  FROM_MM; 
        } else if (wasEmitX > wasMM && wasEmitX > wasEmitY) {
          withMM[i][j] = pxy + wasEmitX;
          traceMM[i][j] =  FROM_EX;
        } else {
          withMM[i][j] = pxy + wasEmitY;
          traceMM[i][j] = FROM_EY;
        }

        //emit in x
        wasMM = transmissionMatrix[3][2] + withMM[i - 1][j];
        wasEmitX = transmissionMatrix[2][2] + withEmitX[i - 1][j];
        
        if (wasMM > wasEmitX) {
          withEmitX[i][j] = qx + wasMM;
          traceX[i][j] =  FROM_MM; 
        } else {
          withEmitX[i][j] = qx + wasEmitX;
          traceX[i][j] = FROM_EX;
        }

        
        //emit in y
        wasMM = transmissionMatrix[3][1] + withMM[i][j - 1];
        wasEmitY = transmissionMatrix[1][1] + withEmitY[i][j - 1];
        
        if (wasMM > wasEmitY) {
          withEmitY[i][j] = qy + wasMM;
          traceY[i][j] =  FROM_MM; 
        } else {
          withEmitY[i][j] = qy + wasEmitY;
          traceY[i][j] = FROM_EY;
        }
    }
  }  

  // backtrace
  std::string alignedX;
  std::string alignedY;
  
  char currentState;
  double wasMM = transmissionMatrix[3][4] + withMM[n][m];
  double wasEmitX = transmissionMatrix[2][4] + withEmitX[n][m];
  double wasEmitY = transmissionMatrix[1][4] + withEmitY[n][m];

  if (wasMM > wasEmitX && wasMM > wasEmitY) {
    currentState = FROM_MM;
  } else if (wasEmitX > wasMM && wasEmitX > wasEmitY) {
    currentState = FROM_EX;
  } else {
    currentState = FROM_EY;
  }
  //printf("%d\n", currentState);

  int i = n, j = m;
  while (true) {
    printf("%d %d (%d)\n", i, j, currentState);
    if (currentState == FROM_MM) {
      alignedX += firstGene[i-1];
      alignedY += secondGene[j-1];
      currentState = traceMM[i][j];
      i--; j--;
    } else if (currentState == FROM_EX) {
      alignedX += firstGene[i-1];
      alignedY += '-';
      currentState = traceX[i][j];
      i--;
    } else if (currentState == FROM_EY) {
      alignedX += '-';
      alignedY += secondGene[j-1];
      currentState = traceY[i][j];
      j--;
    } else {
     // std::cout << "ERROR: BACK TRACE" << std::endl;
    }
    if (i <= 0 && j <= 0) {
      break;
    }
  }
  
  reverse(alignedX.begin(), alignedX.end());
  reverse(alignedY.begin(), alignedY.end());

  for (int i = 0; i < alignedX.size(); i = i + 100) {
    std::cout << std::string(alignedX, i, 100) << std::endl;
    std::cout << std::string(alignedY, i, 100) << std::endl;
    std::cout << "\n\n";
  }

  //std::cout << alignedX << std::endl;
  //std::cout << alignedY << std::endl;
}
