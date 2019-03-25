// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// implementation of cell and face terms without real vector access for
// advection (see evaluation_dg_advect.h for the vector access on a prototype
// Cartesian grid)
//
// Author: Martin Kronbichler, November 2018

#ifndef evaluation_dgele_advect_h
#define evaluation_dgele_advect_h

#include <mpi.h>

#include "gauss_formula.h"
#include "lagrange_polynomials.h"
#include "vectorization.h"
#include "aligned_vector.h"
#include "utilities.h"
#include "matrix_vector_kernel.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif



template <int dim, int degree, typename Number>
class EvaluationCellLaplacian
{
public:
  static constexpr unsigned int dimension = dim;
  static constexpr unsigned int n_q_points = Utilities::pow(degree+1,dim);
  static constexpr unsigned int dofs_per_cell = Utilities::pow(degree+1,dimension);

  void initialize(const unsigned int n_element_batches,
                  const bool         is_cartesian)
  {
    vector_offsets.resize(n_element_batches);
    for (unsigned int i=0; i<n_element_batches; ++i)
      vector_offsets[i] = i*dofs_per_cell;

    input_array.resize(n_element_batches * dofs_per_cell);
    output_array.resize(n_element_batches * dofs_per_cell);

    std::vector<double> jacobian = get_diagonal_jacobian();
    Number jacobian_determinant = 1.;
    for (unsigned int d=0; d<dim; ++d)
      jacobian_determinant *= jacobian[d];
    jacobian_determinant = 1./jacobian_determinant;

    jxw_data.resize(1);
    jxw_data[0] = jacobian_determinant;
    jacobian_data.resize(dim);
    for (unsigned int d=0; d<dim; ++d)
      jacobian_data[d] = jacobian[d];

    convection.resize(dim);
    for (unsigned int d=0; d<dim; ++d)
      {
        convection[d] = -1. + 1.1 * d;

        for (unsigned int e=0; e<dim; ++e)
          {
            normal_jac1[d][e] = VectorizedArray<Number>();
            normal_jac2[d][e] = VectorizedArray<Number>();
            normal_vector[d][e] = VectorizedArray<Number>();
          }
        normal_jac1[d][dim-1] = jacobian[d];
        normal_jac2[d][dim-1] = jacobian[d];
        normal_vector[d][d] = 1;
        Number determinant = 1;
        for (unsigned int e=0; e<dim; ++e)
          if (d!=e)
            determinant *= jacobian[e];
        face_jxw[d] = 1./determinant;
      }

    fill_shape_values();
  }

  std::size_t n_elements() const
  {
    return VectorizedArray<Number>::n_array_elements * vector_offsets.size();
  }

  void do_verification()
  {
    // no full verification implemented
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      input_array[i] = Number(i);
    matrix_vector_product();
    Number sum = 0.;
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      sum += output_array[i][0];
    if (!std::isfinite(sum))
      std::cout << "Wrong value: " << sum << std::endl;
  }

  void matrix_vector_product()
  {
    constexpr unsigned int nn = degree+1;
    constexpr unsigned int mid = nn/2;
    constexpr unsigned int n_lanes = VectorizedArray<Number>::n_array_elements;
    constexpr unsigned int dofs_per_face = Utilities::pow(degree+1,dim-1);
    constexpr unsigned int dofs_per_plane = Utilities::pow(degree+1,2);
    const VectorizedArray<Number> *__restrict shape_values_eo = this->shape_values_eo.begin();
    const VectorizedArray<Number> *__restrict shape_gradients_eo = this->shape_gradients_eo.begin();
    AlignedVector<VectorizedArray<Number> > scratch_data_array;
    VectorizedArray<Number> my_array[degree < 27 ? 2*dofs_per_cell : 1];
    VectorizedArray<Number> *__restrict data_ptr;
    VectorizedArray<Number> array_f[6][dofs_per_face];
    VectorizedArray<Number> merged_array[dim];
    if (degree < 27)
      data_ptr = my_array;
    else
      {
        std::abort();
      }

    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const VectorizedArray<Number> *__restrict input_ptr =
          input_array.begin();
        for (unsigned int i2=0; i2<(dim>2 ? degree+1 : 1); ++i2)
          {
            // x-direction
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, input_ptr+i2*nn*nn+i1*nn, in+i1*nn);
              }
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, in+i1, in+i1);
              }
          }

        if (dim == 3)
          {
            for (unsigned int f=4; f<6; ++f)
              {
                const unsigned int offset1 = (f==4 ? dofs_per_face*degree : 0);
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    for (unsigned int i=0; i<nn; ++i)
                      array_f[f][i1*nn+i] = input_ptr[offset1+i1*nn+i];
                    apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                      (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
                  }
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
              }

            for (unsigned int i2=0; i2<dofs_per_face; ++i2)
              {
                // face integrals in z direction
                {
                  VectorizedArray<Number> normal_speed = -normal_vector[2][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed -= normal_vector[2][d]*convection[d];

                  const VectorizedArray<Number> u_minus = data_ptr[i2];
                  const VectorizedArray<Number> u_plus  = array_f[4][i2];
                  array_f[4][i2] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_jxw[2] * face_quadrature_weight[i2];
                }
                {
                  VectorizedArray<Number> normal_speed = normal_vector[2][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed += normal_vector[2][d]*convection[d];

                  const VectorizedArray<Number> u_minus = data_ptr[degree*dofs_per_face+i2];
                  const VectorizedArray<Number> u_plus  = array_f[5][i2];
                  array_f[5][i2] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_jxw[2] * face_quadrature_weight[i2];
                }

                apply_1d_matvec_kernel<degree+1,dofs_per_face,0,true,false,VectorizedArray<Number>>
                  (shape_values_eo, data_ptr+i2, data_ptr+i2);
              }
          }

        // interpolate external x values for faces
        {
          for (unsigned int f=0; f<2; ++f)
            {
              const unsigned int offset1 = (f==0 ? degree : 0);
              for (unsigned int i1=0; i1<(dim==3 ? nn : 1); ++i1)
                {
                  for (unsigned int i=0; i<nn; ++i)
                    array_f[f][i1*nn+i] = input_ptr[offset1+(i1*nn+i)*nn];

                  // interpolate values onto quadrature points
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
                }
              for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
            }
        }
        // interpolate external y values for faces
        {
          for (unsigned int f=2; f<4; ++f)
            {
              for (unsigned int i1=0; i1<(dim>2 ? (degree+1) : 1); ++i1)
                {
                  const unsigned int base_offset1 = i1*nn*nn+(f==2 ? degree : 0)*nn;
                  for (unsigned int i2=0; i2<degree+1; ++i2)
                    array_f[f][i1*nn+i2] = input_ptr[base_offset1+i2];
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
                }
              for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
            }
        }

        const VectorizedArray<Number> *__restrict jxw_ptr =
          jxw_data.begin();
        const VectorizedArray<Number> *__restrict jacobian_ptr =
          jacobian_data.begin();
        for (unsigned int d=0; d<dim; ++d)
          merged_array[d] = jxw_ptr[0] * jacobian_ptr[d] * convection[d];

        for (unsigned int i2=0; i2<(dim==3 ? degree+1 : 1); ++i2)
          {
            const unsigned int offset = i2*dofs_per_plane;
            VectorizedArray<Number> *array_ptr = data_ptr + offset;
            VectorizedArray<Number> *array_2_ptr = data_ptr + dofs_per_cell + offset;
            const Number *quadrature_ptr = quadrature_weights.begin() + offset;

            VectorizedArray<Number> outy[dofs_per_plane];

            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                // evaluation from cell onto face
                VectorizedArray<Number> r0 = shape_values_on_face_eo[0] * (array_ptr[i]+
                                                                           array_ptr[i+nn*(nn-1)]);
                VectorizedArray<Number> r1 = shape_values_on_face_eo[nn-1] * (array_ptr[i]-
                                                                              array_ptr[i+nn*(nn-1)]);
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_values_on_face_eo[ind] * (array_ptr[i+nn*ind]+
                                                          array_ptr[i+nn*(nn-1-ind)]);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (array_ptr[i+nn*ind]-
                                                               array_ptr[i+nn*(nn-1-ind)]);
                  }
                if (nn%2 == 1)
                  r0 += shape_values_on_face_eo[mid] * array_ptr[i+nn*mid];

                // face integrals in y direction
                {
                  VectorizedArray<Number> normal_speed = -normal_vector[1][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed -= normal_vector[1][d]*convection[d];

                  const VectorizedArray<Number> u_minus = r0 + r1;
                  const VectorizedArray<Number> u_plus  = array_f[2][i1];
                  array_f[2][i1] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_quadrature_weight[i1] * face_jxw[1];
                }
                {
                  VectorizedArray<Number> normal_speed = normal_vector[1][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed += normal_vector[1][d]*convection[d];

                  const VectorizedArray<Number> u_minus = r0 - r1;
                  const VectorizedArray<Number> u_plus  = array_f[3][i1];
                  array_f[3][i1] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_quadrature_weight[i1] * face_jxw[1];
                }

                // evaluation from cell onto face
                r0 = shape_values_on_face_eo[0] * (array_ptr[i*nn]+
                                                   array_ptr[i*nn+(nn-1)]);
                r1 = shape_values_on_face_eo[nn-1] * (array_ptr[i*nn]-
                                                      array_ptr[i*nn+(nn-1)]);
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_values_on_face_eo[ind] * (array_ptr[i*nn+ind]+
                                                          array_ptr[i*nn+(nn-1-ind)]);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (array_ptr[i*nn+ind]-
                                                               array_ptr[i*nn+(nn-1-ind)]);
                  }
                if (nn%2 == 1)
                  r0 += shape_values_on_face_eo[mid] * array_ptr[i*nn+mid];

                // face integrals in x direction
                {
                  VectorizedArray<Number> normal_speed = -normal_vector[0][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed -= normal_vector[0][d]*convection[d];

                  const VectorizedArray<Number> u_minus = r0 + r1;
                  const VectorizedArray<Number> u_plus  = array_f[0][i1];
                  array_f[0][i1] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_quadrature_weight[i1] * face_jxw[0];
                }
                {
                  VectorizedArray<Number> normal_speed = normal_vector[0][0]*convection[0];
                  for (unsigned int d=1; d<dim; ++d)
                    normal_speed += normal_vector[0][d]*convection[d];

                  const VectorizedArray<Number> u_minus = r0 - r1;
                  const VectorizedArray<Number> u_plus  = array_f[1][i1];
                  array_f[1][i1] = -0.5 * ((u_minus+u_plus) * normal_speed +
                                           std::abs(normal_speed) * (u_minus-u_plus)) * face_quadrature_weight[i1] * face_jxw[0];
                }
              }

            // cell integral on quadrature points
            for (unsigned int i1=0; i1<degree+1; ++i1)
              {
                VectorizedArray<Number> outx[nn];
                for (unsigned int i=0; i<degree+1; ++i)
                  {
                    const VectorizedArray<Number> res = array_ptr[i1*nn+i] * quadrature_weights[i2*nn*nn+i1*nn+i];
                    outx[i] = res * merged_array[0];
                    outy[i1*nn+i] = res * merged_array[1];
                    if (dim == 3)
                      array_2_ptr[i1*nn+i] = res * merged_array[2];
                  }
                VectorizedArray<Number> array_face[2];
                array_face[0] = array_f[0][i2*(degree+1)+i1]+array_f[1][i2*(degree+1)+i1];
                array_face[1] = array_f[0][i2*(degree+1)+i1]-array_f[1][i2*(degree+1)+i1];
                apply_1d_matvec_kernel<nn,1,1,false,false,VectorizedArray<Number>,VectorizedArray<Number>,false,1>
                  (shape_gradients_eo, outx, array_ptr+i1*(degree+1),
                   nullptr, shape_values_on_face_eo.begin(), array_face);
              }
            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                VectorizedArray<Number> array_face[2];
                array_face[0] = array_f[2][i2*(degree+1)+i]+array_f[3][i2*(degree+1)+i];
                array_face[1] = array_f[2][i2*(degree+1)+i]-array_f[3][i2*(degree+1)+i];
                apply_1d_matvec_kernel<nn,nn,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,1>
                  (shape_gradients_eo, outy+i, array_ptr+i,
                   array_ptr+i, shape_values_on_face_eo.begin(), array_face);
              }
          }
        if (dim == 3)
          {
            for (unsigned int i2=0; i2<nn*nn; ++i2)
              {
                apply_1d_matvec_kernel<nn,nn*nn,1,false,true,VectorizedArray<Number>>
                  (shape_gradients_eo, data_ptr+dofs_per_cell+i2,
                   data_ptr+i2, data_ptr+i2);
                data_ptr[i2] += array_f[4][i2];
                data_ptr[nn*nn*degree+i2] += array_f[5][i2];

                apply_1d_matvec_kernel<nn,nn*nn,0,false,false,VectorizedArray<Number>>
                  (shape_values_eo, data_ptr+i2, data_ptr+i2);
              }
          }

        for (unsigned int i2=0; i2< (dim>2 ? degree+1 : 1); ++i2)
          {
            const unsigned int offset = i2*dofs_per_plane;
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              apply_1d_matvec_kernel<nn, nn, 0, false, false, VectorizedArray<Number>>
                (shape_values_eo, data_ptr+offset+i1, data_ptr+offset+i1);
            for (unsigned int i1=0; i1<nn; ++i1)
              apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                (shape_values_eo, data_ptr+offset+i1*nn, output_array.begin()+offset+i1*nn);
          }
      }
  }

private:

  std::vector<double> get_diagonal_jacobian() const
  {
    std::vector<double> jacobian(dim);
    for (unsigned int d=0; d<dim; ++d)
      jacobian[d] = 1.+d;
    return jacobian;
  }

  void fill_shape_values()
  {
    constexpr unsigned int n_q_points_1d = degree+1;
    constexpr unsigned int stride = (n_q_points_1d+1)/2;
    shape_values_eo.resize((degree+1)*stride);
    shape_gradients_eo.resize((degree+1)*stride);

    std::vector<double> points = get_gauss_lobatto_points(n_q_points_1d);
    std::vector<double> gauss_points = get_gauss_points(n_q_points_1d);
    std::vector<double> gauss_weights = get_gauss_weights(n_q_points_1d);
    LagrangePolynomialBasis basis(points);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis.value(i, gauss_points[q]);
          const double p2 = basis.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_values_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_values_eo[degree/2*stride+q] =
          basis.value(degree/2, gauss_points[q]);

    LagrangePolynomialBasis basis_gauss(gauss_points);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gauss.derivative(i, gauss_points[q]);
          const double p2 = basis_gauss.derivative(i, gauss_points[n_q_points_1d-1-q]);
          shape_gradients_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_gradients_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_gradients_eo[degree/2*stride+q] =
          basis_gauss.derivative(degree/2, gauss_points[q]);

    shape_values_on_face_eo.resize(2*(degree+1));
    for (unsigned int i=0; i<degree/2+1; ++i)
      {
        const double v0 = basis_gauss.value(i, 0);
        const double v1 = basis_gauss.value(i, 1);
        shape_values_on_face_eo[degree-i] = 0.5 * (v0 - v1);
        shape_values_on_face_eo[i] = 0.5 * (v0 + v1);

        const double d0 = basis_gauss.derivative(i, 0);
        const double d1 = basis_gauss.derivative(i, 1);
        shape_values_on_face_eo[degree+1+i] = 0.5 * (d0 + d1);
        shape_values_on_face_eo[degree+1+degree-i] = 0.5 * (d0 - d1);
      }

    std::vector<double> gauss_weight_1d = get_gauss_weights(n_q_points_1d);
    quadrature_weights.resize(Utilities::pow(n_q_points_1d,dim));
    if (dim == 3)
      for (unsigned int q=0, z=0; z<n_q_points_1d; ++z)
        for (unsigned int y=0; y<n_q_points_1d; ++y)
          for (unsigned int x=0; x<n_q_points_1d; ++x, ++q)
            quadrature_weights[q] = (gauss_weight_1d[z] * gauss_weight_1d[y]) *
              gauss_weight_1d[x];
    else if (dim == 2)
      for (unsigned int q=0, y=0; y<n_q_points_1d; ++y)
        for (unsigned int x=0; x<n_q_points_1d; ++x, ++q)
          quadrature_weights[q] = gauss_weight_1d[y] * gauss_weight_1d[x];
    else if (dim == 1)
      for (unsigned int q=0; q<n_q_points_1d; ++q)
        quadrature_weights[q] = gauss_weight_1d[q];
    else
      throw;

    face_quadrature_weight.resize(Utilities::pow(n_q_points_1d,dim-1));
    if (dim == 3)
      for (unsigned int q=0, y=0; y<n_q_points_1d; ++y)
        for (unsigned int x=0; x<n_q_points_1d; ++x, ++q)
          face_quadrature_weight[q] = gauss_weight_1d[y] * gauss_weight_1d[x];
    else if (dim == 2)
      for (unsigned int q=0; q<n_q_points_1d; ++q)
        face_quadrature_weight[q] = gauss_weight_1d[q];
    else
      face_quadrature_weight[0] = 1.;
  }

  AlignedVector<VectorizedArray<Number> > shape_values_eo;
  AlignedVector<VectorizedArray<Number> > shape_gradients_eo;
  AlignedVector<VectorizedArray<Number> > shape_values_on_face_eo;

  AlignedVector<Number> quadrature_weights;
  AlignedVector<Number> face_quadrature_weight;

  AlignedVector<VectorizedArray<Number> > input_array;
  AlignedVector<VectorizedArray<Number> > output_array;

  AlignedVector<VectorizedArray<Number> > jxw_data;
  AlignedVector<VectorizedArray<Number> > jacobian_data;
  std::array<std::array<VectorizedArray<Number>,dim>,dim> normal_jac1, normal_jac2, normal_vector;
  std::array<VectorizedArray<Number>,dim> face_jxw;

  AlignedVector<VectorizedArray<Number> > convection;

  std::vector<unsigned int> vector_offsets;
};


#endif
