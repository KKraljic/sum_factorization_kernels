
#include "vectorization.h"
#include "gauss_formula.h"

#include <iostream>

int main()
{
  std::vector<double> gauss = get_gauss_points(6);
  std::vector<double> gauss_lobatto = get_gauss_lobatto_points(7);

  std::cout << "Gauss points: ";
  for (unsigned int i=0; i<gauss.size(); ++i)
    std::cout << gauss[i] << " ";
  std::cout << std::endl;

  std::cout << "Gauss-Lobatto points: ";
  for (unsigned int i=0; i<gauss_lobatto.size(); ++i)
    std::cout << gauss_lobatto[i] << " ";
  std::cout << std::endl;

  return 0;
}
