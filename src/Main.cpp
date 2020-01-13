#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <bits/stdc++.h> 
#include <cmath>
//#include "NeedlemanWunsch.h"

#define GET_INDEX(i, j, m) ((i) * (m) + (j)) 
#define FROM_MM 1
#define FROM_EX 2
#define FROM_EY 3

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

/*
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
*/

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
      printf("%9.6f", matrix[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

void ptm(char **matrix, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%3d", matrix[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

//ispisuje matrice dinamike i trace zajedno
void pbm(double **mat, char **t, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%9.6f(%d)", mat[i][j], t[i][j]);
    }
    std::cout << "\n";
  }
  std::cout << "\n\n\n";
}

void viterbi(char g1FileName[], char g2FileName[]) {
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
  
  withMM[0][0] = 1;
  // alociraj pocetno stanje 
  //pm(withMM, n + 1, m + 1);
  //pm(withEmitX, n + 1, m + 1);
  //pm(withEmitY, n + 1, m + 1);

  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
        char x = firstGene[i - 1];
        char y = secondGene[j - 1];
        
        double pxy = emissionMatrix[getIndexOfBase(x)][getIndexOfBase(y)];
        double qx = emissionMatrix[getIndexOfBase(x)][0];
        double qy = emissionMatrix[0][getIndexOfBase(y)];
        
        
        //emit match/missmatch pair
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
        if (i==1 && j == 2) std::cout << wasMM << std::endl;
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
      std::cout << "ERROR: BACK TRACE" << std::endl;
      break;
    }
    
    if (i <= 0 && j <= 0) {
      break;
    }
  }
  
  reverse(alignedX.begin(), alignedX.end());
  reverse(alignedY.begin(), alignedY.end());
  std::cout << alignedX << std::endl;
  std::cout << alignedY << std::endl;
}

void proba_to_log(double p[][5]) {
	for (int i = 0; i < 5; ++i) {
		for(int j = 0; j < 5; ++j) {
			if (fabs(p[i][j]) < 10e-5) p[i][j] = -INFINITY;
			else p[i][j] = std::log(p[i][j]);
		}
	}
}

void to_log(double **p, int n, int m) {
	for (int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			if (fabs(p[i][j]) < 10e-5) p[i][j] = -INFINITY;
			else p[i][j] = std::log(p[i][j]);
		}
	}
}

void viterbi_log(char g1FileName[], char g2FileName[]) {
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
  std::cout << alignedX << std::endl;
  std::cout << alignedY << std::endl;
}

int main(int argc, char *argv[])
{
  viterbi_log(argv[1],argv[2]);
  //new NeedlemanWunsch();
  //std::cout << std::max({1, 3, 4}) << std::endl;

  return 0;
}


