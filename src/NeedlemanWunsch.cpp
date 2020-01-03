#include <iostream>
#include <vector>

class NeedlemanWunsch
{
private:
  std::vector<std::vector<int>> scoreMatrix;
  std::vector<std::vector<int>> tracebackMatrix;

public:
  NeedlemanWunsch();
};

NeedlemanWunsch::NeedlemanWunsch()
{
  std::cout << "Const. called" << std::endl;
}
