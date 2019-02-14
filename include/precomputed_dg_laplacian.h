// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// cell and face terms for DG Laplacian (interior penalty method) using a
// precomputed tensor-product matrix on Cartesian cell geometries (no
// geometric factors implemented)
//
// Author: Martin Kronbichler, April 2018

#ifndef precomputed_dg_laplacian_h
#define precomputed_dg_laplacian_h


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
class PrecomputedDGLaplacian
{
public:
  static constexpr unsigned int dimension = dim;
  static constexpr unsigned int n_q_points = Utilities::pow(degree+1,dim);
  static constexpr unsigned int dofs_per_cell = Utilities::pow(degree+1,dim);
  unsigned int blx;
  unsigned int bly;
  unsigned int blz;

  void initialize(const unsigned int *n_cells_in)
  {
    n_cells[0] = n_cells_in[0]/VectorizedArray<Number>::n_array_elements;
    for (unsigned int d=1; d<dim; ++d)
      n_cells[d] = n_cells_in[d];
    for (unsigned int d=dim; d<3; ++d)
      n_cells[d] = 1;

    n_blocks[2] = (n_cells[2] + blz - 1)/blz;
    n_blocks[1] = (n_cells[1] + bly - 1)/bly;
    n_blocks[0] = (n_cells[0] + blx - 1)/blx;

    input_array.resize(0);
    output_array.resize(0);
    input_array.resize_fast(n_elements() * dofs_per_cell);
    output_array.resize_fast(n_elements() * dofs_per_cell);

#pragma omp parallel
    {
#pragma omp for schedule (static) collapse(2)
      for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
        for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
          for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
            for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
              for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
                {
                  const unsigned int ii=(i*n_cells[1]+j)*n_cells[0];
                  for (std::size_t ix=dofs_per_cell*VectorizedArray<Number>::n_array_elements*(kb*blx+ii);
                       ix<(std::min(n_cells[0], (kb+1)*blx)+ii)*dofs_per_cell*VectorizedArray<Number>::n_array_elements; ++ix)
                    {
                      input_array[ix] = 1;
                      output_array[ix] = 0.;
                    }
                }
    }

    fill_shape_values();
  }

  std::size_t n_elements() const
  {
    std::size_t n_element = VectorizedArray<Number>::n_array_elements;
    for (unsigned int d=0; d<dim; ++d)
      n_element *= n_cells[d];
    return n_element;
  }

  void do_verification()
  {
    std::cout << "Verification currently not implemented!" << std::endl;
  }

  void do_inner_loop (const unsigned int start_x,
                      const unsigned int end_x,
                      const unsigned int iy,
                      const unsigned int iz)
  {
    constexpr unsigned int N = degree+1;
    constexpr unsigned int n_lanes = VectorizedArray<Number>::n_array_elements;
    constexpr unsigned int dofs_per_face = Utilities::pow(degree+1,dim-1);
    VectorizedArray<Number> array_0[dofs_per_cell], array_1[dofs_per_cell], array_2[dofs_per_cell];
    const VectorizedArray<Number> *__restrict laplace_1d_eo = this->laplace_1d_eo.begin();
    const VectorizedArray<Number> *__restrict mass_1d_eo = this->mass_1d_eo.begin();

    for (unsigned int ix=start_x; ix<end_x; ++ix)
      {
        const unsigned int ii=((iz*n_cells[1]+iy)*n_cells[0]+ix)*n_lanes;
        const VectorizedArray<Number>* src_array =
          reinterpret_cast<const VectorizedArray<Number>*>(input_array.begin()+ii*dofs_per_cell);
        VectorizedArray<Number>* dst_array =
          reinterpret_cast<VectorizedArray<Number>*>(output_array.begin()+ii*dofs_per_cell);

        // -------------------------------------------------------------
        // Laplacian in x-direction
        //VectorizedArray<Number> product_diags = inv_jac[0][0];
        //for (unsigned int d=1; d<dim; ++d)
        //product_diags *= inv_jac[d][d];
        //const VectorizedArray<Number> factor_x = inv_jac[0][0] * inv_jac[0][0] / product_diags;
        //for (unsigned int i=0; i<eo_dim; ++i)
        //laplace_1d[i] = laplace_1d_eo[i] * factor_x;
        const VectorizedArray<Number> interface_val_0 = make_vectorized_array(value_outer_1);
        const VectorizedArray<Number> interface_val_1 = make_vectorized_array(value_outer_2);
        unsigned int indices_left[n_lanes];
        for (unsigned int v=1; v<n_lanes; ++v)
          indices_left[v] = ii*dofs_per_cell+v-1;
        indices_left[0] = (ii-n_lanes)*dofs_per_cell+n_lanes-1;
        if (ix==0)
          {
            // assume periodic boundary conditions
            indices_left[0] = (ii+(n_cells[0]-1)*n_lanes)*dofs_per_cell+n_lanes-1;
          }
        unsigned int indices_right[n_lanes];
        for (unsigned int v=0; v<n_lanes-1; ++v)
          indices_right[v] = ii*dofs_per_cell+v+1;
        indices_right[n_lanes-1] = (ii+n_lanes)*dofs_per_cell;
        if (ix==n_cells[0]-1)
          {
            // assume periodic boundary conditions
            indices_right[n_lanes-1] = (ii-(n_cells[0]-1)*n_lanes)+dofs_per_cell;
          }
        for (unsigned int i1=0; i1<dofs_per_face; ++i1)
          {
            const VectorizedArray<Number> *__restrict in = src_array + i1*N;
            VectorizedArray<Number> *__restrict out = array_0 + i1*N;

            apply_1d_matvec_kernel<N, 1, 2, true, false, VectorizedArray<Number>>(laplace_1d_eo, in, out);

            VectorizedArray<Number> val1, val2;
            val1.gather(input_array.begin()+(N*(i1+1)-1)*n_lanes, indices_left);
            val2.gather(input_array.begin()+(N*(i1+1)-2)*n_lanes, indices_left);
            out[0] += val1 * interface_val_0;
            out[0] += val2 * interface_val_1;
            out[1] += val1 * interface_val_1;

            val1.gather(input_array.begin()+(N*i1)*n_lanes, indices_right);
            val2.gather(input_array.begin()+(N*i1+1)*n_lanes, indices_right);
            out[N-1] += val1 * interface_val_0;
            out[N-1] += val2 * interface_val_1;
            out[N-2] += val1 * interface_val_1;
          }

        // -------------------------------------------------------------
        // Laplacian in y-direction
        const unsigned int index_left = (iy > 0 ?
                                         (ii-n_cells[0]*n_lanes) :
                                         (ii+(n_cells[1]-1)*n_cells[0]*n_lanes)
                                         ) * dofs_per_cell;
        const unsigned int index_right = (iy < n_cells[1]-1 ?
                                          (ii+n_cells[0]*n_lanes) :
                                          (ii-(n_cells[1]-1)*n_cells[0]*n_lanes)
                                          ) * dofs_per_cell;
        for (unsigned int i1=0; i1<(dim>2 ? N : 1); ++i1)
          {
            const VectorizedArray<Number> *__restrict in = src_array + i1*N*N;
            VectorizedArray<Number> *__restrict out = array_1 + i1*N*N;
            for (unsigned int i2=0; i2<N; ++i2)
              {
                apply_1d_matvec_kernel<N, N, 2, true, false, VectorizedArray<Number>>(laplace_1d_eo, in+i2,
                                                                     out+i2);

                VectorizedArray<Number> val1, val2;
                val1.load(input_array.begin()+index_left+(i1*N*N+(N-1)*N+i2)*n_lanes);
                val2.load(input_array.begin()+index_left+(i1*N*N+(N-2)*N+i2)*n_lanes);
                out[i2] += val1 * interface_val_0;
                out[i2] += val2 * interface_val_1;
                out[N+i2] += val1 * interface_val_1;
                val1.load(input_array.begin()+index_right+(i1*N*N+i2)*n_lanes);
                val2.load(input_array.begin()+index_right+(i1*N*N+N+i2)*n_lanes);
                out[N*(N-1)+i2] += val1 * interface_val_0;
                out[N*(N-1)+i2] += val2 * interface_val_1;
                out[N*(N-2)+i2] += val1 * interface_val_1;
              }
          }

        // -------------------------------------------------------------
        // Laplacian in z-direction
        if (dim > 2)
          {
            const unsigned int index_left = (iz > 0 ?
                                             (ii-n_cells[1]*n_cells[0]*n_lanes) :
                                             (ii+(n_cells[2]-1)*n_cells[1]*n_cells[0]*n_lanes)
                                             )*dofs_per_cell;
            const unsigned int index_right = (iz < n_cells[2]-1 ?
                                              (ii+n_cells[1]*n_cells[0]*n_lanes) :
                                              (ii-(n_cells[2]-1)*n_cells[1]*n_cells[0]*n_lanes)
                                              )*dofs_per_cell;
            for (unsigned int i1=0; i1<dofs_per_face; ++i1)
              {
                const VectorizedArray<Number> *__restrict in = src_array + i1;
                VectorizedArray<Number> *__restrict out = array_2 + i1;

                apply_1d_matvec_kernel<N, N*N, 2, true, false, VectorizedArray<Number>>(laplace_1d_eo, in, out);

                VectorizedArray<Number> val1, val2;
                val1.load(input_array.begin()+index_left+(N*N*(N-1)+i1)*n_lanes);
                val2.load(input_array.begin()+index_left+(N*N*(N-2)+i1)*n_lanes);
                out[0] += val1 * interface_val_0;
                out[0] += val2 * interface_val_1;
                out[N*N] += val1 * interface_val_1;
                val1.load(input_array.begin()+index_right+(i1)*n_lanes);
                val2.load(input_array.begin()+index_right+(N*N+i1)*n_lanes);
                out[N*N*(N-1)] += val1 * interface_val_0;
                out[N*N*(N-1)] += val2 * interface_val_1;
                out[N*N*(N-2)] += val1 * interface_val_1;
              }
          }

        // mass matrices
        if (dim == 3)
          {
            // z mass matrices
            for (int i1=0; i1<N*N; ++i1)
              apply_1d_matvec_kernel<N, N*N, 2, true, false, VectorizedArray<Number>>(mass_1d_eo,
                                                                     array_0+i1, array_0+i1);
            for (int i1=0; i1<N*N; ++i1)
              apply_1d_matvec_kernel<N, N*N, 2, true, false, VectorizedArray<Number>>(mass_1d_eo,
                                                                     array_1+i1, array_1+i1);

            for (int i1=0; i1<N; ++i1)
              {
                // y mass matrices
                for (int i2=0; i2<N; ++i2)
                  apply_1d_matvec_kernel<N, N, 2, true, false, VectorizedArray<Number>>(mass_1d_eo,
                                                                       array_0+i1*N*N+i2,
                                                                       array_0+i1*N*N+i2);
                for (int i2=0; i2<N; ++i2)
                  apply_1d_matvec_kernel<N, N, 2, true, true, VectorizedArray<Number>>(mass_1d_eo,
                                                                      array_2+i1*N*N+i2,
                                                                      array_1+i1*N*N+i2,
                                                                      array_1+i1*N*N+i2);

                // x mass matrix
                for (int i2=0; i2<N; ++i2)
                  apply_1d_matvec_kernel<N, 1, 2, true, true, VectorizedArray<Number>, VectorizedArray<Number>,/*NT store = */ true>
                    (mass_1d_eo, array_1+i1*N*N+i2*N,
                     dst_array+i1*N*N+i2*N, array_0+i1*N*N+i2*N);
              }
          }
        else if (dim == 2)
          {
            // y mass matrix
            for (int i1=0; i1<N; ++i1)
              apply_1d_matvec_kernel<N, N, 2, true, false, VectorizedArray<Number>>(mass_1d_eo, array_1+i1,
                                                                   array_1+i1);

            // x mass matrix
            for (int i1=0; i1<N; ++i1)
              apply_1d_matvec_kernel<N, 1, 2, true, true, VectorizedArray<Number>, VectorizedArray<Number>, /*NT store = */ true>
                (mass_1d_eo, array_0 + i1*N, dst_array+i1*N, array_1 + i1*N);
          }
        else
          throw;
      }
  }

  void matrix_vector_product()
  {
    if (degree < 1)
      return;

#pragma omp parallel
    {
#ifdef LIKWID_PERFMON
      LIKWID_MARKER_START(("dg_laplacian_" + std::to_string(dim) +
                           "d_deg_" + std::to_string(degree)).c_str());
#endif

#pragma omp for schedule (static) collapse(2)
      for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
        for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
          for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
            for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
              for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
                do_inner_loop(kb*blx, std::min(n_cells[0], (kb+1)*blx), j, i);

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_STOP(("dg_laplacian_" + std::to_string(dim) +
                          "d_deg_" + std::to_string(degree)).c_str());
#endif
    }
  }

private:
  void fill_shape_values()
  {
    const unsigned int N = degree+1;
    AlignedVector<double> mass_1d(N*N);
    AlignedVector<double> lapl_1d(N*N);

    std::vector<double> gauss_points = get_gauss_points(N);
    std::vector<double> gauss_weight_1d = get_gauss_weights(N);
    const double penalty = (N-1)*(N)*1.;

    HermiteLikePolynomialBasis basis(degree);
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        {
          double sum_m = 0, sum_l = 0;
          for (unsigned int q=0; q<gauss_points.size(); ++q)
            {
              sum_m += (basis.value(i, gauss_points[q]) *
                        basis.value(j, gauss_points[q])
                        )* gauss_weight_1d[q];
              sum_l += (basis.derivative(i, gauss_points[q]) *
                        basis.derivative(j, gauss_points[q])
                        )* gauss_weight_1d[q];
            }
          mass_1d[i*N+j] = sum_m;
          sum_l += (1.*basis.value(i, 0.)*basis.value(j, 0.)*penalty
                    +
                    0.5*basis.derivative(i, 0.)*basis.value(j, 0.)
                    +
                    0.5*basis.derivative(j, 0.)*basis.value(i, 0.));
          sum_l += (1.*basis.value(i, 1.)*basis.value(j, 1.)*penalty
                    -
                    0.5*basis.derivative(i, 1.)*basis.value(j, 1.)
                    -
                    0.5*basis.derivative(j, 1.)*basis.value(i, 1.));
          lapl_1d[i*N+j] = sum_l;
        }
    value_outer_1 = (-1.*basis.value(0, 1.)*basis.value(N-1, 0.)*penalty
                     +
                     0.5*basis.derivative(0, 1.)*basis.value(N-1, 0.)
                     -
                     0.5*basis.derivative(N-1, 0.)*basis.value(0, 1.));
    value_outer_2 = (-1.*basis.value(0, 1.)*basis.value(N-2, 0.)*penalty
                     +
                     0.5*basis.derivative(0, 1.)*basis.value(N-2, 0.)
                     -
                     0.5*basis.derivative(N-2, 0.)*basis.value(0, 1.));

    const unsigned int stride = (N+1)/2;

    mass_1d_eo.resize(N*stride);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          mass_1d_eo[i*stride+q] = 0.5 * (mass_1d[i*N+q] + mass_1d[i*N+N-1-q]);
          mass_1d_eo[(degree-i)*stride+q] = 0.5 * (mass_1d[i*N+q] - mass_1d[i*N+N-1-q]);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        mass_1d_eo[degree/2*stride+q] = mass_1d[(N/2)*N+q];
    laplace_1d_eo.resize(N*stride);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          laplace_1d_eo[i*stride+q] = 0.5 * (lapl_1d[i*N+q] + lapl_1d[i*N+N-1-q]);
          laplace_1d_eo[(degree-i)*stride+q] = 0.5 * (lapl_1d[i*N+q] - lapl_1d[i*N+N-1-q]);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        laplace_1d_eo[degree/2*stride+q] = lapl_1d[(N/2)*N+q];
  }

  unsigned int n_cells[3];
  unsigned int n_blocks[3];
  AlignedVector<VectorizedArray<Number> > mass_1d_eo;
  AlignedVector<VectorizedArray<Number> > laplace_1d_eo;

  Number value_outer_1, value_outer_2;

  AlignedVector<Number> input_array;
  AlignedVector<Number> output_array;
};


#endif
