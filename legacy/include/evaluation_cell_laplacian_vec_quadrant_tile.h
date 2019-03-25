// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// manual implementation of cell terms for Laplacian for hardcoded
// vectorization within a cell
//
// Author: Martin Kronbichler, July 2017

#ifndef evaluation_cell_laplacian_h
#define evaluation_cell_laplacian_h

#include "gauss_formula.h"
#include "lagrange_polynomials.h"
#include "vectorization.h"
#include "aligned_vector.h"
#include "utilities.h"
#include "matrix_vector_kernel.h"

//#define COPY_ONLY_BENCHMARK
//#define DO_MASS_MATRIX
//#define DO_CONVECTION
#define READ_SINGLE_VECTOR


template <int dim, int degree, typename Number>
class EvaluationCellLaplacianVecEle
{
public:
  static_assert(degree%2==1, "Only odd degrees implemented");
  static const unsigned int dimension = dim;
  static const unsigned int n_q_points = Utilities::pow(degree+1,dim);
  static const unsigned int dofs_per_cell = Utilities::pow(degree+1,dim);

  void initialize(const unsigned int n_element_batches,
                  const bool         is_cartesian)
  {
    vector_offsets.resize(n_element_batches);
    for (unsigned int i=0; i<n_element_batches; ++i)
      vector_offsets[i] = i*dofs_per_cell/VectorizedArray<Number>::n_array_elements;

    input_array.resize(n_element_batches * dofs_per_cell/VectorizedArray<Number>::n_array_elements);
    output_array.resize(n_element_batches * dofs_per_cell/VectorizedArray<Number>::n_array_elements);

    fill_shape_values();

    std::vector<double> jacobian = get_diagonal_jacobian();
    Number jacobian_determinant = 1.;
    for (unsigned int d=0; d<dim; ++d)
      jacobian_determinant *= jacobian[d];
    jacobian_determinant = 1./jacobian_determinant;

    data_offsets.resize(n_element_batches);
    if (is_cartesian)
      {
        jacobian_data.resize(dim*(dim+1)/2);
        for (unsigned int d=0; d<dim; ++d)
          jacobian_data[d] = jacobian_determinant * jacobian[d] * jacobian[d];
        for (unsigned int i=0; i<n_element_batches; ++i)
          data_offsets[i] = 0;
      }
    else
      {
        jacobian_data.resize(n_element_batches*n_q_points*dim*(dim+1)/2);
        for (unsigned int i=0; i<n_element_batches; ++i)
          {
            data_offsets[i] = i*n_q_points;
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int d=0; d<dim; ++d)
                jacobian_data[(data_offsets[i]+q)*dim*(dim+1)/2+d]
                  = jacobian_determinant * jacobian[d] * jacobian[d] * quadrature_weights[q];
          }
      }
  }

  std::size_t n_elements() const
  {
    return vector_offsets.size();
  }

  void do_verification()
  {
    // check that the Laplacian applied to a linear function in all of the
    // directions equals to the value of the linear function at the boundary.
    std::vector<double> points = get_gauss_lobatto_points(degree+1);
    std::vector<double> gauss_points = get_gauss_points(degree+1);
    std::vector<double> gauss_weights = get_gauss_weights(degree+1);
    LagrangePolynomialBasis gll(points);
    std::vector<Number> values_1d(points.size());
    std::vector<double> jacobian = get_diagonal_jacobian();

    // compute boundary integral of basis functions
    AlignedVector<Number> boundary_integral(Utilities::pow(degree+1,dim-1));
    for (unsigned int i1=0, i=0; i1<(dim>2?degree+1:1); ++i1)
      for (unsigned int i0=0; i0<(dim>1?degree+1:1); ++i0, ++i)
        {
          Number sum = 0;
          if (dim == 3)
            for (unsigned int q1=0; q1<degree+1; ++q1)
              for (unsigned int q0=0; q0<degree+1; ++q0)
                sum += gll.value(i1, gauss_points[q1]) * gll.value(i0, gauss_points[q0]) * gauss_weights[q0] * gauss_weights[q1];
          else if (dim == 2)
            for (unsigned int q0=0; q0<degree+1; ++q0)
              sum += gll.value(i0, gauss_points[q0]) * gauss_weights[q0];
          boundary_integral[i] = sum;
        }
    std::vector<unsigned int> renumbering(dofs_per_cell);
    for (unsigned int z=0, p=0; z<(dim>2?degree+1:1); ++z)
      for (unsigned int y=0; y<(dim>1?degree+1:1)/2; ++y)
        for (unsigned int x=0; x<(degree+1)/2; ++x, ++p)
          {
            renumbering[4*p] = z*(degree+1)*(degree+1)+y*(degree+1)+x;
            renumbering[4*p+1] = z*(degree+1)*(degree+1)+y*(degree+1)+(degree-x);
            renumbering[4*p+2] = z*(degree+1)*(degree+1)+(degree-y)*(degree+1)+x;
            renumbering[4*p+3] = z*(degree+1)*(degree+1)+(degree-y)*(degree+1)+(degree-x);
          }
    {
      std::vector<unsigned int> renumbering_tmp(renumbering);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        renumbering[renumbering_tmp[i]] = i;
    }
    for (unsigned int test=0; test<dim; ++test)
      {
        for (unsigned int q=0; q<=degree; ++q)
          values_1d[q] = points[q]/jacobian[test];

        // set a linear function on each cell whose derivative will evaluate
        // to zero except at the boundary of the element
        unsigned int indices[3];
        for (unsigned int i=0; i<vector_offsets.size(); ++i)
          {
            Number *data_ptr = &input_array[vector_offsets[i]][0];
            indices[2] = 0;
            for (unsigned int p=0; indices[2]<(dim>2?degree+1:1); ++indices[2])
              for (indices[1]=0; indices[1]<(dim>1?degree+1:1)/2; ++indices[1])
                for (indices[0]=0; indices[0]<(degree+1)/2; ++indices[0], ++p)
                  {
#if 1
                    data_ptr[4*p] = values_1d[indices[test]];
                    data_ptr[4*p+1] = values_1d[test==0 ? degree-indices[test] : indices[test]];
                    data_ptr[4*p+2] = values_1d[test==1 ? degree-indices[test] : indices[test]];
                    data_ptr[4*p+3] = values_1d[test<2 ? degree-indices[test] : indices[test]];
#else
                    data_ptr[4*p] = values_1d[indices[0]]*values_1d[indices[1]];
                    data_ptr[4*p+1] = values_1d[degree-indices[0]]*values_1d[indices[1]];
                    data_ptr[4*p+2] = values_1d[indices[0]]*values_1d[degree-indices[1]];
                    data_ptr[4*p+3] = values_1d[degree-indices[0]]*values_1d[degree-indices[1]];
#endif

                  }
          }

        matrix_vector_product();

        // remove the boundary integral from the cell integrals and check the
        // error
        double boundary_factor = 1.;
        for (unsigned int d=0; d<dim; ++d)
          if (d!=test)
            boundary_factor /= jacobian[d];
        double max_error = 0;
#ifndef READ_SINGLE_VECTOR
        for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
#else
        unsigned int cell = 0;
#endif
          {
            Number *data_ptr = &output_array[vector_offsets[cell]][0];
            const unsigned int stride = test < dim-1 ? (degree+1) : 1;
            int shift = 1;
            for (unsigned int d=0; d<test; ++d)
              shift *= degree+1;
            if (test != 1)
              {
                // normal vector at left is negative, must add boundary
                // contribution
                for (unsigned int i=0; i<Utilities::pow(degree+1,dim-1); ++i)
                  data_ptr[renumbering[i*stride]] += boundary_factor * boundary_integral[i];
                // normal vector at left is positive, must subtract boundary
                // contribution
                for (unsigned int i=0; i<Utilities::pow(degree+1,dim-1); ++i)
                  data_ptr[renumbering[degree*shift + i*stride]] -= boundary_factor * boundary_integral[i];
              }
            else
              {
                for (unsigned int j=0; j<=(dim>2?degree:0); ++j)
                  for (unsigned int i=0; i<=degree; ++i)
                    {
                      const unsigned int ind = j*Utilities::pow(degree+1,dim-1) + i;
                      const unsigned int l = dim>2 ? i*(degree+1)+j : i;
                      data_ptr[renumbering[ind]] += boundary_factor * boundary_integral[l];
                      data_ptr[renumbering[degree*shift+ind]] -= boundary_factor * boundary_integral[l];
                    }
              }
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              max_error = std::max(max_error, (double)data_ptr[i]);
            //if (max_error > 1e-10)
            //  {
            //    std::cout << "error: ";
            //  for (unsigned int i=0; i<dofs_per_cell; ++i)
            //    std::cout << data_ptr[i] << " ";
            //  std::cout << std::endl;
            //  }
          }

          std::cout << "Error of integral in direction " << test << ": "
                    << max_error << std::endl;
      }
  }

  void matrix_vector_product()
  {
    if (degree < 1)
      return;

    constexpr unsigned int dofs_per_cell_z = dofs_per_cell / VectorizedArray<Number>::n_array_elements;
    AlignedVector<VectorizedArray<Number> > scratch_data_array;
    VectorizedArray<Number> my_array[degree < 27 ? 4*dofs_per_cell_z : 1];
    VectorizedArray<Number> *__restrict data_ptr;
    if (degree < 27)
      data_ptr = my_array;
    else
      {
        if (scratch_data_array.size() != 4*dofs_per_cell_z)
          scratch_data_array.resize_fast(4*dofs_per_cell_z);
        data_ptr = scratch_data_array.begin();
      }

    const bool is_cartesian = jacobian_data.size() == (dim*(dim+1))/2;
    const unsigned int nn = degree+1;
    const unsigned int mid = nn/2;
    const unsigned int offset = (nn+1)/2;
    VectorizedArray<Number> one_x, one_y;
    static_assert(VectorizedArray<Number>::n_array_elements==4,
                  "Only AVX case implemented right now");
    for (unsigned int i=0; i<VectorizedArray<Number>::n_array_elements; i+=2)
      {
        one_x[i] = 1;
        one_x[i+1] = -1;
      }
    for (unsigned int i=0; i<VectorizedArray<Number>::n_array_elements; i+=2)
      {
        one_y[i/2] = 1;
        one_y[i/2+2] = -1;
      }
    const VectorizedArray<Number> *__restrict shape_val_z  = shape_values.begin();
    const VectorizedArray<Number> *__restrict shape_grad_z = shape_gradients.begin();
    const VectorizedArray<Number> *__restrict shape_val_x  = shape_values_x.begin();
    const VectorizedArray<Number> *__restrict shape_grad_x = shape_gradients_x.begin();
    const VectorizedArray<Number> *__restrict shape_val_y  = shape_values_y.begin();
    const VectorizedArray<Number> *__restrict shape_grad_y = shape_gradients_y.begin();

    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const VectorizedArray<Number> *__restrict input_ptr =
          input_array.begin();
        VectorizedArray<Number> *__restrict output_ptr =
          output_array.begin();
        // --------------------------------------------------------------------
        // apply tensor product kernels
        for (unsigned int i2=0; i2<(dim>2 ? nn : 1); ++i2)
          {
            // x-direction
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = input_ptr[i2*mid*mid+i1*mid+i];
                    xp[i].data = _mm256_permute_pd(x.data, 0x5) + (one_x * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_val_x[col]               * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_val_x[ind*offset+col] * xp[ind];
                    data_ptr[i2*mid*mid+i1*mid+col].data = _mm256_permute_pd(r0.data, 0x5) + (one_x * r0).data;
                  }
              }

            // y-direction
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[i2*mid*mid+i*mid+i1];
                    xp[i].data = _mm256_permute2f128_pd(x.data, x.data, 0x5) + (one_y * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_val_y[col]               * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_val_y[ind*offset+col] * xp[ind];
                    data_ptr[i2*mid*mid+col*mid+i1].data = _mm256_permute2f128_pd(r0.data, r0.data, 0x5) + (one_y * r0).data;
                  }
              }
          }
        if (dim == 3)
          {
            for (unsigned int i1=0; i1<mid*mid; ++i1)
              {
                apply_1d_matvec_kernel<nn, mid*mid, 0, true, false, Number>
                  (shape_values, data_ptr+i1, data_ptr+i1);
                apply_1d_matvec_kernel<nn, mid*mid, 1, true, false, Number>
                  (shape_gradients, data_ptr+i1, data_ptr+i1+3*dofs_per_cell_z);
              }
          }
        for (unsigned int i2=0; i2<(dim>2 ? nn : 1); ++i2)
          {
            // y-direction
            const VectorizedArray<Number> one_y1 = -one_y;
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[i2*mid*mid+i*mid+i1];
                    xp[i].data = x.data + one_y1.data * _mm256_permute2f128_pd(x.data, x.data, 0x5);
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_grad_y[col]               * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grad_y[ind*offset+col] * xp[ind];
                    data_ptr[2*dofs_per_cell_z+i2*mid*mid+col*mid+i1].data = _mm256_permute2f128_pd(r0.data, r0.data, 0x5) + (one_y * r0).data;
                  }
              }

            // x-direction
            const VectorizedArray<Number> one_x1 = -one_x;
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[i2*mid*mid+i1*mid+i];
                    xp[i].data = x.data + one_x1.data * _mm256_permute_pd(x.data, 0x5);
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_grad_x[col]               * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grad_x[ind*offset+col] * xp[ind];
                    data_ptr[dofs_per_cell_z+i2*mid*mid+i1*mid+col].data = _mm256_permute_pd(r0.data, 0x5) + (one_x * r0).data;
                  }
              }

            constexpr unsigned int n_q_points = dofs_per_cell_z;
            VectorizedArray<Number>* phi_grads = data_ptr + dofs_per_cell_z + mid*mid*i2;
            for (unsigned int q=0; q<mid*mid; ++q)
              {
                if (dim==2)
                  {
                    VectorizedArray<Number> t0 = phi_grads[q]*jacobian_data[0] + phi_grads[q+n_q_points]*jacobian_data[2];
                    VectorizedArray<Number> t1 = phi_grads[q]*jacobian_data[2] + phi_grads[q+n_q_points]*jacobian_data[1];
                    phi_grads[q] = t0 * quadrature_weights[i2*mid*mid+q];
                    phi_grads[q+n_q_points] = t1 * quadrature_weights[i2*mid*mid+q];
                  }
                else if (dim==3)
                  {
                    VectorizedArray<Number> t0 = phi_grads[q]*jacobian_data[0] + phi_grads[q+n_q_points]*jacobian_data[3]+phi_grads[q+2*n_q_points]*jacobian_data[4];
                    VectorizedArray<Number> t1 = phi_grads[q]*jacobian_data[3] + phi_grads[q+n_q_points]*jacobian_data[1]+phi_grads[q+2*n_q_points]*jacobian_data[5];
                    VectorizedArray<Number> t2 = phi_grads[q]*jacobian_data[4] + phi_grads[q+n_q_points]*jacobian_data[5]+phi_grads[q+2*n_q_points]*jacobian_data[2];
                    phi_grads[q] = t0 * quadrature_weights[i2*mid*mid+q];
                    phi_grads[q+n_q_points] = t1 * quadrature_weights[i2*mid*mid+q];
                    phi_grads[q+2*n_q_points] = t2 * quadrature_weights[i2*mid*mid+q];
                  }
              }

            // x-direction
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[dofs_per_cell_z+i2*mid*mid+i1*mid+i];
                    xp[i].data = _mm256_permute_pd(x.data, 0x5) + (one_x * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_grad_x[col*offset]        * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grad_x[col*offset+ind] * xp[ind];
                    data_ptr[i2*mid*mid+i1*mid+col].data = (one_x.data*_mm256_permute_pd(r0.data, 0x5)) + r0.data;
                  }
              }

            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[2*dofs_per_cell_z+i2*mid*mid+i*mid+i1];
                    xp[i].data = _mm256_permute2f128_pd(x.data, x.data, 0x5) + (one_y * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_grad_y[col*offset]        * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grad_y[col*offset+ind] * xp[ind];
                    data_ptr[i2*mid*mid+col*mid+i1].data += one_y.data*_mm256_permute2f128_pd(r0.data, r0.data, 0x5) + r0.data;
                  }
              }
          }

        if (dim == 3)
          {
            for (unsigned int i1=0; i1<mid*mid; ++i1)
              {
                apply_1d_matvec_kernel<nn, mid*mid, 1, false, true, Number>
                  (shape_gradients, data_ptr+i1+3*dofs_per_cell_z, data_ptr+i1, data_ptr+i1);
                apply_1d_matvec_kernel<nn, mid*mid, 0, false, false, Number>
                  (shape_values, data_ptr+i1, data_ptr+i1);
              }
          }

        for (unsigned int i2=0; i2<(dim>2 ? nn : 1); ++i2)
          {
            // y-direction
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[i2*mid*mid+i*mid+i1];
                    xp[i].data = _mm256_permute2f128_pd(x.data, x.data, 0x5) + (one_y * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_val_y[col*offset]        * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_val_y[col*offset+ind] * xp[ind];
                    data_ptr[i2*mid*mid+col*mid+i1].data = _mm256_permute2f128_pd(r0.data, r0.data, 0x5) + (one_y * r0).data;
                  }
              }

            // x-direction
            for (unsigned int i1=0; i1<mid; ++i1)
              {
                VectorizedArray<Number> xp[mid];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    const VectorizedArray<Number> x = data_ptr[i2*mid*mid+i1*mid+i];
                    xp[i].data = _mm256_permute_pd(x.data, 0x5) + (one_x * x).data;
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0;
                    r0 = shape_val_x[col*offset]        * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_val_x[col*offset+ind] * xp[ind];
                    output_ptr[i2*mid*mid+i1*mid+col].data = _mm256_permute_pd(r0.data, 0x5) + (one_x * r0).data;
                  }
              }
          }
      }
  }

private:
  std::vector<double> get_diagonal_jacobian() const
  {
    std::vector<double> jacobian(dim);
    jacobian[0] = 4.;
    for (unsigned int d=1; d<dim; ++d)
      {
        double entry = (double)(d+1)/4.;
        jacobian[d] = 1./entry;
      }
    return jacobian;
  }

  void fill_shape_values()
  {
    const unsigned int n_q_points_1d = degree+1;
    const unsigned int stride = (n_q_points_1d+1)/2;
    shape_values.resize((degree+1)*stride);
    shape_gradients.resize((degree+1)*stride);
    shape_values_x.resize(stride*stride);
    shape_gradients_x.resize(stride*stride);
    shape_values_y.resize(stride*stride);
    shape_gradients_y.resize(stride*stride);

    LagrangePolynomialBasis basis_gll(get_gauss_lobatto_points(degree+1));
    std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gll.value(i, gauss_points[q]);
          const double p2 = basis_gll.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values[i*stride+q] = 0.5 * (p1 + p2);
          shape_values[(degree-i)*stride+q] = 0.5 * (p1 - p2);
          for (unsigned int v=0; v<2; ++v)
            {
              shape_values_x[i*stride+q][2*v] = shape_values[i*stride+q][0];
              shape_values_x[i*stride+q][2*v+1] = shape_values[(degree-i)*stride+q][0];
              shape_values_y[i*stride+q][v] = shape_values[i*stride+q][0];
              shape_values_y[i*stride+q][2+v] = shape_values[(degree-i)*stride+q][0];
            }
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_values[degree/2*stride+q] =
          basis_gll.value(degree/2, gauss_points[q]);

    LagrangePolynomialBasis basis_gauss(get_gauss_points(degree+1));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gauss.derivative(i, gauss_points[q]);
          const double p2 = basis_gauss.derivative(i, gauss_points[n_q_points_1d-1-q]);
          shape_gradients[i*stride+q] = 0.5 * (p1 + p2);
          shape_gradients[(degree-i)*stride+q] = 0.5 * (p1 - p2);
          for (unsigned int v=0; v<2; ++v)
            {
              shape_gradients_x[i*stride+q][2*v] = shape_gradients[i*stride+q][0];
              shape_gradients_x[i*stride+q][2*v+1] = shape_gradients[(degree-i)*stride+q][0];
              shape_gradients_y[i*stride+q][v] = shape_gradients[i*stride+q][0];
              shape_gradients_y[i*stride+q][2+v] = shape_gradients[(degree-i)*stride+q][0];
            }
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_gradients[degree/2*stride+q] =
          basis_gauss.derivative(degree/2, gauss_points[q]);

    // get quadrature weights
    std::vector<double> gauss_weight_1d = get_gauss_weights(n_q_points_1d);
    quadrature_weights.resize((dim>2?n_q_points_1d:1)*Utilities::pow(n_q_points_1d/2,2));
    if (dim == 3)
      for (unsigned int q=0, z=0; z<n_q_points_1d; ++z)
        for (unsigned int y=0; y<n_q_points_1d/2; ++y)
          for (unsigned int x=0; x<n_q_points_1d/2; ++x, ++q)
            quadrature_weights[q] = (gauss_weight_1d[z] * gauss_weight_1d[y]) *
              gauss_weight_1d[x];
    else if (dim == 2)
      for (unsigned int q=0, y=0; y<n_q_points_1d/2; ++y)
        for (unsigned int x=0; x<n_q_points_1d/2; ++x, ++q)
          quadrature_weights[q] = gauss_weight_1d[y] * gauss_weight_1d[x];
    else
      throw;
  }

  AlignedVector<VectorizedArray<Number> > shape_values;
  AlignedVector<VectorizedArray<Number> > shape_gradients;

  AlignedVector<VectorizedArray<Number> > shape_values_x;
  AlignedVector<VectorizedArray<Number> > shape_gradients_x;

  AlignedVector<VectorizedArray<Number> > shape_values_y;
  AlignedVector<VectorizedArray<Number> > shape_gradients_y;

  AlignedVector<Number> quadrature_weights;

  AlignedVector<unsigned int> vector_offsets;
  AlignedVector<VectorizedArray<Number> > input_array;
  AlignedVector<VectorizedArray<Number> > output_array;

  AlignedVector<unsigned int> data_offsets;
  AlignedVector<VectorizedArray<Number> > jacobian_data;
};


#endif
