#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>
#include "HMMAlign.h"

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

std::string getGene(std::string fileName) {
  std::ifstream infile1(fileName);
  std::string gene;
  std::string line;
  while (std::getline(infile1, line)) {
    gene += line;
  }
  return gene;
}

void fillEmissionMatrix(std::string fileName) {
  std::ifstream infile(fileName);
  std::string line;
  int i = 0;
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    for (int j = 0; j < 5; j++) {
      ss >> emissionMatrix[i][j];
    }
    i++;
  }
}

void fillTransmissionMatrix(std::string fileName) {
  std::ifstream infile(fileName);
  std::string line;
  int i = 0;
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    for (int j = 0; j < 5; j++) {
      ss >> transmissionMatrix[i][j];
    }
    i++;
  }
}


int main(int argc, char *argv[]) {
  char ch;
  std::string emissionFilename;
  std::string transFilename;
  int emiss = 0, trans = 0;
  while ((ch = getopt(argc, argv, "e:t:")) != -1) {
    switch (ch) {
      case 'e':
        emiss = 1;
        emissionFilename = optarg;
        break;
      case 't':
        trans = 1;
        transFilename = optarg;
        break;
      default:
        printf("Kriva opcija\n");
        return -1;
    }        
  }
  if (emiss) {
    fillEmissionMatrix(emissionFilename);
  }
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << emissionMatrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  if (trans) {
    fillTransmissionMatrix(transFilename);
  }

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << transmissionMatrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  HMMAlign align(getGene(argv[optind]), getGene(argv[optind + 1]), transmissionMatrix, emissionMatrix);
  align.run();
  return 0;
}

