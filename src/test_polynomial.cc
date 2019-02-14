#include "gauss_formula.h"
#include "lagrange_polynomials.h"

#include <iostream>
#include <iomanip>

int main()
{
  LagrangePolynomialBasis basis(get_gauss_lobatto_points(7));

  std::cout << "Evaluation value list: " << std::endl;
  const unsigned int n_points = 20;
  for (unsigned int i=0; i<=n_points; ++i)
    {
      const double x = (double)i/n_points;
      std::cout << "x=" << std::setw(6) << x << " ";
      for (unsigned int p=0; p<=basis.degree(); ++p)
        std::cout << std::setw(11) << basis.value(p, x) << " ";
      std::cout << std::endl;
    }

  std::cout << "Evaluation derivative list: " << std::endl;
  for (unsigned int i=0; i<=n_points; ++i)
    {
      const double x = (double)i/n_points;
      std::cout << "x=" << std::setw(6) << x << " ";
      for (unsigned int p=0; p<=basis.degree(); ++p)
        std::cout << std::setw(11) << basis.derivative(p, x) << " ";
      std::cout << std::endl;
    }

  for (unsigned int i=0; i<6; ++i)
    {
      HermiteLikePolynomialBasis basis(i);

      std::cout << "Evaluation value list Hermite degree " << i << ": " << std::endl;
      const unsigned int n_points = 20;
      for (unsigned int i=0; i<=n_points; ++i)
        {
          const double x = (double)i/n_points;
          std::cout << "x=" << std::setw(6) << x << " ";
          for (unsigned int p=0; p<=basis.degree(); ++p)
            std::cout << std::setw(11) << basis.value(p, x) << " ";
          std::cout << std::endl;
        }

      std::cout << "Evaluation derivative list Hermite degree " << i << ": " << std::endl;
      for (unsigned int i=0; i<=n_points; ++i)
        {
          const double x = (double)i/n_points;
          std::cout << "x=" << std::setw(6) << x << " ";
          for (unsigned int p=0; p<=basis.degree(); ++p)
            std::cout << std::setw(11) << basis.derivative(p, x) << " ";
          std::cout << std::endl;
        }
    }

  return 0;
}
