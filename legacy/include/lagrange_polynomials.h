// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// This file is inspired by the file deal.II/source/base/polynomial.cc, see
// www.dealii.org for information about licenses. Here is the original deal.II
// license statement:
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef lagrange_polynomials_h
#define lagrange_polynomials_h

#include "gauss_formula.h"
#include <limits>
#include <vector>

class LagrangePolynomialBasis
{
public:
  LagrangePolynomialBasis(const std::vector<double> &points)
    :
    points(points)
  {
    lagrange_denominators.resize(points.size());
    for (unsigned int i=0; i<points.size(); ++i)
      {
        double denominator = 1.;
        for (unsigned int j=0; j<points.size(); ++j)
          if (j!=i)
            denominator *= points[i] - points[j];
        lagrange_denominators[i] = 1./denominator;
      }
  }

  unsigned int degree() const
  {
    return points.size()-1;
  }

  double value(const unsigned int polynomial_index,
               const double x) const
  {
    double value = lagrange_denominators[polynomial_index];
    for (unsigned int i=0; i<points.size(); ++i)
      if (polynomial_index != i)
        value *= x-points[i];
    return value;
  }

  double derivative(const unsigned int polynomial_index,
                    const double       x) const
  {
    double value = lagrange_denominators[polynomial_index];
    double derivative = 0;
    for (unsigned int i=0; i<points.size(); ++i)
      if (i != polynomial_index)
        {
          const double v = x-points[i];
          derivative = derivative * v + value;
          value *= v;
        }
    return derivative;
  }

private:
  std::vector<double> points;
  std::vector<double> lagrange_denominators;
};



class HermiteLikePolynomialBasis
{
private:
  static double find_support_point_x_star (const std::vector<double> &jacobi_roots)
  {
    // Initial guess for the support point position values: The zero turns
    // out to be between zero and the first root of the Jacobi polynomial,
    // but the algorithm is agnostic about that, so simply choose two points
    // that are sufficiently far apart.
    double guess_left = 0;
    double guess_right = 0.5;
    const unsigned int degree = jacobi_roots.size() + 3;

    // Compute two integrals of the product of l_0(x) * l_1(x)
    // l_0(x) = (x-y)*(x-jacobi_roots(0))*...*(x-jacobi_roos(degree-4))*(x-1)*(x-1)
    // l_1(x) = (x-0)*(x-jacobi_roots(0))*...*(x-jacobi_roots(degree-4))*(x-1)*(x-1)
    // where y is either guess_left or guess_right for the two integrals.
    // Note that the polynomials are not yet normalized here, which is not
    // necessary because we are only looking for the x_star where the matrix
    // entry is zero, for which the constants do not matter.
    std::vector<double> points = get_gauss_points(degree+1);
    std::vector<double> weights = get_gauss_weights(degree+1);
    double integral_left = 0, integral_right = 0;
    for (unsigned int q=0; q<weights.size(); ++q)
      {
        const double x = points[q];
        double poly_val_common = x;
        for (unsigned int j=0; j<degree-3; ++j)
          poly_val_common *= (x-jacobi_roots[j])*(x-jacobi_roots[j]);
        poly_val_common *= (x - 1.)*(x - 1.)*(x - 1.)*(x - 1.);
        integral_left += weights[q]*(poly_val_common*(x - guess_left));
        integral_right += weights[q]*(poly_val_common*(x - guess_right));
      }

    // compute guess by secant method. Due to linearity in the root x_star,
    // this is the correct position after this single step
    return guess_right - (guess_right-guess_left)/(integral_right-integral_left)*integral_right;
  }

public:
  HermiteLikePolynomialBasis(const unsigned int degree)
  {
    points.resize(degree+1);
    if (degree >= 3)
      {
        std::vector<double> jacobi_roots = zeros_of_jacobi_polynomial(degree+1, 2, 2);

        // Note that the Jacobi roots are given for the interval [-1,1], so we
        // must scale the eigenvalues to the interval [0,1] before using them
        for (unsigned int i=0; i<degree-3; ++i)
          jacobi_roots[i] = 0.5*jacobi_roots[i]+0.5;

        points[0] = 0.;
        points[1] = 0.;
        for (unsigned int m=0; m<degree-3; m++)
          points[m+2] = jacobi_roots[m];
        points[degree-1] = 1.;
        points[degree] = 1.;

        zero_one = find_support_point_x_star(jacobi_roots);

        lagrange_denominators.resize(degree+1, 1.);
        for (unsigned int m=2; m<degree-1; ++m)
          lagrange_denominators[m] = 1./this->value(m, points[m]);
        lagrange_denominators[0] = 1./this->value(0, 0.);
        lagrange_denominators[degree] = 1./this->value(degree, 1.);
        lagrange_denominators[1] = -this->derivative(0, 0.) / this->derivative(1, 0.);
        lagrange_denominators[degree-1] = -(this->derivative(degree, 1.) /
                                            this->derivative(degree-1, 1.));
      }
  }

  unsigned int degree() const
  {
    return points.size()-1;
  }

  double value(const unsigned int polynomial_index,
               const double x) const
  {
    if (points.size() == 1)
      return 1.;
    else if (points.size() == 2)
      {
        if (polynomial_index == 0)
          return 1.-x;
        else
          return x;
      }
    else if (points.size() == 3)
      {
        if (polynomial_index == 0)
          return (1.-x)*(1.-x);
        else if (polynomial_index == 1)
          return 2.*x*(1.-x);
        else
          return x*x;
      }
    else
      {
        double value = lagrange_denominators[polynomial_index];
        if (polynomial_index == 0)
          {
            value *= (x-zero_one);
            for (unsigned int i=2; i<points.size(); ++i)
              value *= x-points[i];
          }
        else if (polynomial_index == points.size()-1)
          {
            value *= (x-1.+zero_one);
            for (unsigned int i=0; i<points.size()-2; ++i)
              value *= x-points[i];
          }
        else
          {
            for (unsigned int i=0; i<points.size(); ++i)
              if (i != polynomial_index)
                value *= x-points[i];
          }
        return value;
      }
  }

  double derivative(const unsigned int polynomial_index,
                    const double       x) const
  {
    if (points.size() == 1)
      return 0;
    else if (points.size() == 2)
      {
        if (polynomial_index == 0)
          return -1.;
        else
          return 1.;
      }
    else if (points.size() == 3)
      {
        if (polynomial_index == 0)
          return -2.*(1.-x);
        else if (polynomial_index == 1)
          return 2. - 4.*x;
        else
          return 2.*x;
      }
    else
      {
        double value = lagrange_denominators[polynomial_index];
        double derivative = 0;
        if (polynomial_index == 0)
          {
            derivative = value;
            value *= (x-zero_one);
            for (unsigned int i=2; i<points.size(); ++i)
              {
                const double v = x-points[i];
                derivative = derivative * v + value;
                value *= v;
              }
          }
        else if (polynomial_index == points.size()-1)
          {
            derivative = value;
            value *= (x-1.+zero_one);
            for (unsigned int i=0; i<points.size()-2; ++i)
              {
                const double v = x-points[i];
                derivative = derivative * v + value;
                value *= v;
              }
          }
        else
          {
            for (unsigned int i=0; i<points.size(); ++i)
              if (i != polynomial_index)
                {
                  const double v = x-points[i];
                  derivative = derivative * v + value;
                  value *= v;
                }
          }
        return derivative;
      }
  }

private:
  std::vector<double> points;
  std::vector<double> lagrange_denominators;
  double zero_one;
};

#endif
