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
    VectorizedArray<Number> my_array[degree < 27 ? (degree < 5 ? (degree+7) : (2*degree+2))*dofs_per_face : 1];
    VectorizedArray<Number> *__restrict data_ptr;
    VectorizedArray<Number> array_f[6][dofs_per_face], array_fd[6][dofs_per_face];
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
            const VectorizedArray<Number> sigma = Number((degree+1) * (degree+1)) * jacobian_data[2];
            for (unsigned int f=4; f<6; ++f)
              {
                const VectorizedArray<Number> w0 = (f==4 ? 1. : -1.)*hermite_derivative_on_face;
                const unsigned int offset1 = (f==4 ? dofs_per_face*degree : 0);
                const unsigned int offset2 = dofs_per_face * (f==4 ? degree-1 : 1);
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    for (unsigned int i=0; i<nn; ++i)
                      {
                        array_f[f][i1*nn+i] = input_ptr[offset1+i1*nn+i];
                        array_fd[f][i1*nn+i] = w0 * (input_ptr[offset2+i1*nn+i]-array_f[f][i1*nn+i]);
                      }
                    apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                      (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
                    apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                      (shape_values_eo, array_fd[f]+i1*nn, array_fd[f]+i1*nn);
                  }
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_fd[f]+i1, array_fd[f]+i1);
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients_eo, array_f[f]+i1*nn, data_ptr+nn*nn*nn+4*nn*nn+i1*nn);
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients_eo, array_f[f]+i1, data_ptr+nn*nn*nn+5*nn*nn+i1);

                const unsigned int offset3 = (f==4 ? 0 : dofs_per_face*degree);
                const unsigned int offset4 = dofs_per_face * (f==4 ? 1 : degree-1);
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    for (unsigned int i=0; i<nn; ++i)
                      {
                        data_ptr[nn*nn*nn+i1*nn+i] = data_ptr[offset3 + i1*nn+i];
                        data_ptr[nn*nn*nn+3*nn*nn+i1*nn+i] = w0 *
                          (data_ptr[offset4 + i1*nn+i] - data_ptr[offset3 + i1*nn+i]);
                      }
                    apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                      (shape_gradients_eo, data_ptr+nn*nn*nn + i1*nn,
                       data_ptr+nn*nn*nn+nn*nn+i1*nn);
                  }
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients_eo, data_ptr+nn*nn*nn + i1,
                     data_ptr+nn*nn*nn+2*nn*nn+i1);

                // face integrals in z direction
                for (unsigned int q=0; q<dofs_per_face; ++q)
                  {
                    VectorizedArray<Number> average_valgrad = data_ptr[nn*nn*nn+q] * normal_jac1[f][0];
                    average_valgrad += data_ptr[nn*nn*nn+q+2*nn*nn] * normal_jac1[2][1];
                    average_valgrad += data_ptr[nn*nn*nn+q+3*nn*nn] * normal_jac1[2][2];
                    average_valgrad += data_ptr[nn*nn*nn+q+4*nn*nn] * normal_jac2[2][0];
                    average_valgrad += data_ptr[nn*nn*nn+q+5*nn*nn] * normal_jac2[2][1];
                    average_valgrad += array_fd[f][q] * normal_jac2[2][2];
                    VectorizedArray<Number> average_value = 0.5 * (data_ptr[nn*nn*nn+q+nn*nn] + array_f[f][q]);
                    const VectorizedArray<Number> weight = -face_quadrature_weight[q] * face_jxw[f/2];
                    array_fd[f][q] = average_value * weight * normal_jac1[2][2];
                    data_ptr[nn*nn*nn+q] = average_value * weight * normal_jac1[2][0];
                    data_ptr[nn*nn*nn+q+nn*nn] = average_value * weight * normal_jac1[2][1];
                    array_f[f][q] = (average_valgrad + average_value * sigma) * weight;
                  }
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients_eo, data_ptr+nn*nn*nn+nn*nn+i1,
                     array_f[f] + i1, array_f[f] + i1);
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients_eo, data_ptr+nn*nn*nn+i1*nn,
                     array_f[f] + i1*nn, array_f[f] + i1*nn);
              }

            for (unsigned int i2=0; i2<nn; ++i2)
              {
                // interpolate in z direction
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<degree+1,nn*nn,0,true,false,VectorizedArray<Number>>
                    (shape_values_eo, data_ptr+i2*nn+i1, data_ptr+i2*nn+i1);

                // interpolate onto x faces
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    VectorizedArray<Number> r0, r1, r2, r3;
                    {
                      const VectorizedArray<Number> t0 = data_ptr[i1*nn*nn+i2*nn];
                      const VectorizedArray<Number> t1 = data_ptr[i1*nn*nn+i2*nn+nn-1];
                      r0 = shape_values_on_face_eo[0] * (t0+t1);
                      r1 = shape_values_on_face_eo[nn-1] * (t0-t1);
                      r2 = shape_values_on_face_eo[nn] * (t0-t1);
                      r3 = shape_values_on_face_eo[2*nn-1] * (t0+t1);
                    }
                    for (unsigned int ind=1; ind<mid; ++ind)
                      {
                        const VectorizedArray<Number> t0 = data_ptr[i1*nn*nn+i2*nn+ind];
                        const VectorizedArray<Number> t1 = data_ptr[i1*nn*nn+i2*nn+nn-1-ind];
                        r0 += shape_values_on_face_eo[ind] * (t0+t1);
                        r1 += shape_values_on_face_eo[nn-1-ind] * (t0-t1);
                        r2 += shape_values_on_face_eo[nn+ind] * (t0-t1);
                        r3 += shape_values_on_face_eo[2*nn-1-ind] * (t0+t1);
                      }
                    if (nn%2 == 1)
                      {
                        r0 += shape_values_on_face_eo[mid] * data_ptr[i1*nn*nn+i2*nn+mid];
                        r3 += shape_values_on_face_eo[nn+mid] * data_ptr[i1*nn*nn+i2*nn+mid];
                      }
                    array_f[0][i1*nn+i2] = r0 + r1;
                    array_f[1][i1*nn+i2] = r0 - r1;
                    array_fd[0][i1*nn+i2] = r2 + r3;
                    array_fd[1][i1*nn+i2] = r2 - r3;
                  }
              }
          }
        else
          {
            for (unsigned int i2=0; i2<nn; ++i2)
              {
                VectorizedArray<Number> r0, r1, r2, r3;
                {
                  const VectorizedArray<Number> t0 = data_ptr[i2*nn];
                  const VectorizedArray<Number> t1 = data_ptr[i2*nn+nn-1];
                  r0 = shape_values_on_face_eo[0] * (t0+t1);
                  r1 = shape_values_on_face_eo[nn-1] * (t0-t1);
                  r2 = shape_values_on_face_eo[nn] * (t0-t1);
                  r3 = shape_values_on_face_eo[2*nn-1] * (t0+t1);
                }
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    const VectorizedArray<Number> t0 = data_ptr[i2*nn+ind];
                    const VectorizedArray<Number> t1 = data_ptr[i2*nn+nn-1-ind];
                    r0 += shape_values_on_face_eo[ind] * (t0+t1);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (t0-t1);
                    r2 += shape_values_on_face_eo[nn+ind] * (t0-t1);
                    r3 += shape_values_on_face_eo[2*nn-1-ind] * (t0+t1);
                  }
                if (nn%2 == 1)
                  {
                    r0 += shape_values_on_face_eo[mid] * data_ptr[i2*nn+mid];
                    r3 += shape_values_on_face_eo[nn+mid] * data_ptr[i2*nn+mid];
                  }
                array_f[0][i2] = r0 + r1;
                array_f[1][i2] = r0 - r1;
                array_fd[0][i2] = r2 + r3;
                array_fd[1][i2] = r2 - r3;
              }
          }

        // interpolate internal y values onto faces
        for (unsigned int i1=0; i1<(dim==3?nn:1); ++i1)
          {
            for (unsigned int i2=0; i2<nn; ++i2)
              {
                VectorizedArray<Number> r0, r1, r2, r3;
                {
                  const VectorizedArray<Number> t0 = data_ptr[i1*nn*nn+i2];
                  const VectorizedArray<Number> t1 = data_ptr[i1*nn*nn+i2+(nn-1)*nn];
                  r0 = shape_values_on_face_eo[0] * (t0+t1);
                  r1 = shape_values_on_face_eo[nn-1] * (t0-t1);
                  r2 = shape_values_on_face_eo[nn] * (t0-t1);
                  r3 = shape_values_on_face_eo[2*nn-1] * (t0+t1);
                }
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    const VectorizedArray<Number> t0 = data_ptr[i1*nn*nn+i2+ind*nn];
                    const VectorizedArray<Number> t1 = data_ptr[i1*nn*nn+i2+(nn-1-ind)*nn];
                    r0 += shape_values_on_face_eo[ind] * (t0+t1);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (t0-t1);
                    r2 += shape_values_on_face_eo[nn+ind] * (t0-t1);
                    r3 += shape_values_on_face_eo[2*nn-1-ind] * (t0+t1);
                  }
                if (nn%2 == 1)
                  {
                    r0 += shape_values_on_face_eo[mid] * data_ptr[i1*nn*nn+i2+mid*nn];
                    r3 += shape_values_on_face_eo[nn+mid] * data_ptr[i1*nn*nn+i2+mid*nn];
                  }
                if (dim == 3)
                  {
                    array_f[2][i2*nn+i1] = r0 + r1;
                    array_f[3][i2*nn+i1] = r0 - r1;
                    array_fd[2][i2*nn+i1] = r2 + r3;
                    array_fd[3][i2*nn+i1] = r2 - r3;
                  }
                else
                  {
                    array_f[2][i2] = r0 + r1;
                    array_f[3][i2] = r0 - r1;
                    array_fd[2][i2] = r2 + r3;
                    array_fd[3][i2] = r2 - r3;
                  }
              }
          }

        for (unsigned int f=0; f<4; ++f)
          {
            // interpolate external values for faces
            const unsigned int stride1 = Utilities::pow(degree+1,(f/2+1)%dim);
            const unsigned int stride2 = Utilities::pow(degree+1,(f/2+2)%dim);
            const unsigned int offset1 = ((f%2==0) ? degree : 0) * Utilities::pow(degree+1,(f/2)%dim);
            const unsigned int offset2 = ((f%2==0) ? degree-1 : 1) * Utilities::pow(degree+1,(f/2)%dim);
            const VectorizedArray<Number> w0 = (f%2==0 ? 1. : -1.)*hermite_derivative_on_face;
            VectorizedArray<Number> *tmp_ptr = data_ptr + dofs_per_cell;
            for (unsigned int i2=0; i2<(dim==3 ? nn : 1); ++i2)
              {
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    tmp_ptr[i2*nn+i1] = data_ptr[offset1 + i2*stride2 + i1*stride1];
                    tmp_ptr[i2*nn+i1+dofs_per_face] = w0 *
                      (data_ptr[offset2+i2*stride2+i1*stride1]-tmp_ptr[i2*nn+i1]);
                  }
                // interpolate values onto quadrature points
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, tmp_ptr+i2*nn, tmp_ptr+i2*nn);
                // interpolate values onto quadrature points
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, tmp_ptr+dofs_per_face+i2*nn, tmp_ptr+dofs_per_face+i2*nn);
              }
            for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
              apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                (shape_values_eo, tmp_ptr+dofs_per_face+i1, tmp_ptr+dofs_per_face+i1);
            for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
              apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                (shape_values_eo, tmp_ptr+i1, tmp_ptr+i1);
            for (unsigned int i1=0; i1<(dim==3 ? nn : 1); ++i1)
              apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                (shape_gradients_eo, tmp_ptr+i1*nn, tmp_ptr+2*dofs_per_face+i1*nn);
            for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
              apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                (shape_gradients_eo, tmp_ptr+i1, tmp_ptr+3*dofs_per_face+i1);
            for (unsigned int i1=0; i1<(dim==3 ? nn : 1); ++i1)
              apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                (shape_gradients_eo, array_f[f]+i1*nn, tmp_ptr+4*dofs_per_face+i1*nn);
            for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
              apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                (shape_gradients_eo, array_f[f]+i1, tmp_ptr+5*dofs_per_face+i1);

            const VectorizedArray<Number> sigma = Number((degree+1) * (degree+1)) * jacobian_data[f/2];
            const VectorizedArray<Number> * jac1 = f%2==0 ? normal_jac1[f/2].data() : normal_jac2[f/2].data();
            const VectorizedArray<Number> * jac2 = f%2==0 ? normal_jac2[f/2].data() : normal_jac1[f/2].data();
            for (unsigned int q=0; q<dofs_per_face; ++q)
              {
                VectorizedArray<Number> average_valgrad = array_fd[f][q]  * jac1[dim-1];
                average_valgrad += tmp_ptr[dofs_per_face+q] * jac2[dim-1];
                average_valgrad += tmp_ptr[4*dofs_per_face+q] * jac1[0];
                average_valgrad += tmp_ptr[2*dofs_per_face+q] * jac2[0];
                if (dim==3)
                  {
                    average_valgrad += tmp_ptr[5*dofs_per_face+q] * jac1[1];
                    average_valgrad += tmp_ptr[3*dofs_per_face+q] * jac2[1];
                  }
                VectorizedArray<Number> average_value = 0.5 * (tmp_ptr[q] + array_f[f][q]);
                const VectorizedArray<Number> weight = -face_quadrature_weight[q] * face_jxw[f/2];
                array_fd[f][q] = average_value * weight * jac1[dim-1];
                tmp_ptr[q] = average_value * weight * jac1[0];
                if (dim==3)
                  tmp_ptr[dofs_per_face+q] = average_value * weight * jac1[1];
                array_f[f][q] = (average_valgrad + average_value * sigma) * weight;
              }
            for (unsigned int i1=0; i1<(dim==3?nn:0); ++i1)
              apply_1d_matvec_kernel<nn, nn, 1, false, true, VectorizedArray<Number>>
                (shape_gradients_eo, tmp_ptr+dofs_per_face+i1,
                 array_f[f] + i1, array_f[f] + i1);
            for (unsigned int i1=0; i1<(dim==3?nn:1); ++i1)
              apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                (shape_gradients_eo, tmp_ptr+i1*nn,
                 array_f[f] + i1*nn, array_f[f] + i1*nn);
          }

        const VectorizedArray<Number> *__restrict jxw_ptr =
          jxw_data.begin();
        const VectorizedArray<Number> *__restrict jacobian_ptr =
          jacobian_data.begin();
        for (unsigned int d=0; d<dim; ++d)
          merged_array[d] = jxw_ptr[0] * jacobian_ptr[d] * convection[d];

        if (dim==3)
          for (unsigned int i2=0; i2<nn*nn; ++i2)
            apply_1d_matvec_kernel<degree+1,nn*nn,1,true,false,VectorizedArray<Number>>
              (shape_gradients_eo, data_ptr+i2, data_ptr+nn*nn*nn+i2);

        for (unsigned int i2=0; i2<(dim==3 ? degree+1 : 1); ++i2)
          {
            const unsigned int offset = i2*dofs_per_plane;
            VectorizedArray<Number> *array_ptr = data_ptr + offset;
            VectorizedArray<Number> *array_2_ptr = data_ptr + dofs_per_cell + offset;
            const Number *quadrature_ptr = quadrature_weights.begin() + offset;

            VectorizedArray<Number> outy[dofs_per_plane];
            // y-derivative
            for (unsigned int i1=0; i1<nn; ++i1) // loop over x layers
              {
                apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients_eo, array_ptr+i1, outy+i1);
              }

            // x-derivative
            for (unsigned int i1=0; i1<degree+1; ++i1) // loop over y layers
              {
                VectorizedArray<Number> outx[nn];
                apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients_eo, array_ptr+i1*nn, outx);

                for (unsigned int i=0; i<degree+1; ++i)
                  {
                    const VectorizedArray<Number> weight =
                      make_vectorized_array(quadrature_weights[i2*nn*nn+i1*nn+i]);
                    if (dim==2)
                      {
                        VectorizedArray<Number> t0 = outy[i1*nn+i]*merged_array[2] + outx[i]*merged_array[0];
                        VectorizedArray<Number> t1 = outy[i1*nn+i]*merged_array[1] + outx[i]*merged_array[2];
                        outx[i] = t0 * weight;
                        outy[i1*nn+i] = t1 * weight;
                        }
                    else if (dim==3)
                      {
                        VectorizedArray<Number> t0 = outy[i1*nn+i]*merged_array[3]+array_2_ptr[i1*nn+i]*merged_array[4] + outx[i]*merged_array[0];
                        VectorizedArray<Number> t1 = outy[i1*nn+i]*merged_array[1]+array_2_ptr[i1*nn+i]*merged_array[5] + outx[i]*merged_array[3];
                        VectorizedArray<Number> t2 = outy[i1*nn+i]*merged_array[5]+array_2_ptr[i1*nn+i]*merged_array[2] + outx[i]*merged_array[4];
                        outx[i] = t0 * weight;
                        outy[i1*nn+i] = t1 * weight;
                        array_2_ptr[i1*nn+i] = t2 * weight;
                      }
                  }
                VectorizedArray<Number> array_face[4];
                array_face[0] = array_f[0][i2*nn+i1]+array_f[1][i2*nn+i1];
                array_face[1] = array_f[0][i2*nn+i1]-array_f[1][i2*nn+i1];
                array_face[2] = array_fd[0][i2*nn+i1]+array_fd[1][i2*nn+i1];
                array_face[3] = array_fd[0][i2*nn+i1]-array_fd[1][i2*nn+i1];
                apply_1d_matvec_kernel<nn,1,1,false,false,VectorizedArray<Number>,VectorizedArray<Number>,false,2>
                  (shape_gradients_eo, outx, array_ptr+i1*nn,
                   nullptr, shape_values_on_face_eo.begin(), array_face);
              }
            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                VectorizedArray<Number> array_face[4];
                array_face[0] = array_f[2][i1]+array_f[3][i1];
                array_face[1] = array_f[2][i1]-array_f[3][i1];
                array_face[2] = array_fd[2][i1]+array_fd[3][i1];
                array_face[3] = array_fd[2][i1]-array_fd[3][i1];
                apply_1d_matvec_kernel<nn,nn,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,2>
                  (shape_gradients_eo, outy+i, array_ptr+i,
                   array_ptr+i, shape_values_on_face_eo.begin(), array_face);
              }
          }
        if (dim == 3)
          {
            for (unsigned int i2=0; i2<nn*nn; ++i2)
              {
                VectorizedArray<Number> array_face[4];
                array_face[0] = array_f[4][i2]+array_f[5][i2];
                array_face[1] = array_f[4][i2]-array_f[5][i2];
                array_face[2] = array_fd[4][i2]+array_fd[5][i2];
                array_face[3] = array_fd[4][i2]-array_fd[5][i2];
                apply_1d_matvec_kernel<nn,nn*nn,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,2>
                  (shape_gradients_eo, data_ptr+dofs_per_cell+i2,
                   data_ptr+i2, data_ptr+i2, shape_values_on_face_eo.begin(),
                   array_face);

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

    std::vector<double> gauss_points = get_gauss_points(n_q_points_1d);
    std::vector<double> gauss_weights = get_gauss_weights(n_q_points_1d);
    HermiteLikePolynomialBasis basis(degree);
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

    hermite_derivative_on_face = basis.derivative(0, 0);

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

  VectorizedArray<Number> hermite_derivative_on_face;

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
