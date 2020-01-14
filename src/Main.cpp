#include <iostream>
#include <fstream>
#include <string>
#include "HMMAlign.h"

std::string getGene(std::string fileName) {
  std::ifstream infile1(fileName);
  std::string gene;
  std::string line;
  while (std::getline(infile1, line)) {
    gene += line;
  }
  return gene;
}

int main(int argc, char *argv[]) {
  HMMAlign align(getGene(argv[1]), getGene(argv[2]));
  align.run();
  return 0;
}

