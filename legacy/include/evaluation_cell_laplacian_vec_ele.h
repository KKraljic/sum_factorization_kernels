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

//#define COPY_ONLY_BENCHMARK
//#define DO_MASS_MATRIX
//#define DO_CONVECTION
#define READ_SINGLE_VECTOR


template <int dim, int degree, typename Number>
class EvaluationCellLaplacianVecEle
{
public:
  static const unsigned int dimension = dim;
  static const unsigned int n_q_points = Utilities::pow(degree+1,dim);
  static const unsigned int dofs_per_cell = Utilities::pow(degree+1,dim);

  void initialize(const unsigned int n_element_batches,
                  const bool         is_cartesian)
  {
    vector_offsets.resize(n_element_batches);
    for (unsigned int i=0; i<n_element_batches; ++i)
      vector_offsets[i] = i*dofs_per_cell;

    input_array.resize(n_element_batches * dofs_per_cell);
    output_array.resize(n_element_batches * dofs_per_cell);

    fill_shape_values();

    std::vector<double> jacobian = get_diagonal_jacobian();
    Number jacobian_determinant = 1.;
    for (unsigned int d=0; d<dim; ++d)
      jacobian_determinant *= jacobian[d];
    jacobian_determinant = 1./jacobian_determinant;

    data_offsets.resize(n_element_batches);
    if (is_cartesian)
      {
        jxw_data.resize(1);
        jxw_data[0] = jacobian_determinant;
        jacobian_data.resize(dim);
        for (unsigned int d=0; d<dim; ++d)
          jacobian_data[d] = jacobian[d];
        for (unsigned int i=0; i<n_element_batches; ++i)
          data_offsets[i] = 0;
      }
    else
      {
        const unsigned int my_n_q_points = quadrature_weights_vect.size();
        jxw_data.resize(n_element_batches*my_n_q_points);
        jacobian_data.resize(n_element_batches*my_n_q_points*dim*dim);
        for (unsigned int i=0; i<n_element_batches; ++i)
          {
            data_offsets[i] = i*my_n_q_points;
            for (unsigned int q=0; q<my_n_q_points; ++q)
              jxw_data[data_offsets[i]+q] =
                jacobian_determinant * quadrature_weights_vect[q];
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int d=0; d<dim; ++d)
                jacobian_data[data_offsets[i]*dim*dim+
                              q*dim*dim+d*dim+d] = jacobian[d];
          }
      }

    convection.resize(dim*dofs_per_cell*(is_cartesian ? 1 : n_element_batches));
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
    for (unsigned int test=0; test<dim; ++test)
      {
        for (unsigned int q=0; q<=degree; ++q)
          values_1d[q] = points[q]/jacobian[test];

        // set a linear function on each cell whose derivative will evaluate
        // to zero except at the boundary of the element
        unsigned int indices[3];
        for (unsigned int i=0; i<vector_offsets.size(); ++i)
          {
            Number *data_ptr = &input_array[vector_offsets[i]];
            indices[2] = 0;
            for (unsigned int p=0; indices[2]<(dim>2?degree+1:1); ++indices[2])
              for (indices[1]=0; indices[1]<(dim>1?degree+1:1); ++indices[1])
                for (indices[0]=0; indices[0]<degree+1; ++indices[0], ++p)
                  data_ptr[p] = values_1d[indices[test]];
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
            Number *data_ptr = &output_array[vector_offsets[cell]];
            const unsigned int stride = test < dim-1 ? (degree+1) : 1;
            int shift = 1;
            for (unsigned int d=0; d<test; ++d)
              shift *= degree+1;
            if (test != 1)
              {
                // normal vector at left is negative, must add boundary
                // contribution
                for (unsigned int i=0; i<Utilities::pow(degree+1,dim-1); ++i)
                  data_ptr[i*stride] += boundary_factor * boundary_integral[i];
                // normal vector at left is positive, must subtract boundary
                // contribution
                for (unsigned int i=0; i<Utilities::pow(degree+1,dim-1); ++i)
                  data_ptr[degree*shift + i*stride] -= boundary_factor * boundary_integral[i];
              }
            else
              {
                for (unsigned int j=0; j<=(dim>2?degree:0); ++j)
                  for (unsigned int i=0; i<=degree; ++i)
                    {
                      const unsigned int ind = j*Utilities::pow(degree+1,dim-1) + i;
                      const unsigned int l = dim>2 ? i*(degree+1)+j : i;
                      data_ptr[ind] += boundary_factor * boundary_integral[l];
                      data_ptr[degree*shift+ind] -= boundary_factor * boundary_integral[l];
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
    static_assert(dim == 2 || dim==3, "Only 3D implemented.");

    if (degree < 1)
      return;

    // do not exchange data or zero out, assume DG operator does not need to
    // exchange and that the loop below takes care of the zeroing

    // global data structures
#ifdef COPY_ONLY_BENCHMARK
    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const Number *__restrict input_ptr =
          input_array.begin()+vector_offsets[cell];
        Number *__restrict output_ptr =
          output_array.begin()+vector_offsets[cell];
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          output_ptr[i] += input_ptr[i] * jxw_data[0];
      }
#else
    AlignedVector<VectorizedArray<Number> > scratch_data_array;
    VectorizedArray<Number> my_array[degree < 27 ? 2*dofs_per_cell : 1];
    VectorizedArray<Number> *__restrict data_ptr;
    if (degree < 27)
      data_ptr = my_array;
    else
      {
        if (scratch_data_array.size() != 2*dofs_per_cell)
          scratch_data_array.resize_fast(2*dofs_per_cell);
        data_ptr = scratch_data_array.begin();
      }
    VectorizedArray<Number> merged_array[dim*(dim+1)/2];
    for (unsigned int d=0; d<dim*(dim+1)/2; ++d)
      merged_array[d] = VectorizedArray<Number>();

    const bool is_cartesian = jxw_data.size() == 1;
    const unsigned int nn = degree+1;
    const unsigned int mid = nn/2;
    const unsigned int offset = (nn+1)/2;
    const VectorizedArray<Number> *__restrict shape_vals = shape_values.begin();
    const VectorizedArray<Number> *__restrict shape_grads = shape_gradients.begin();

    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const Number *__restrict input_ptr =
#ifdef READ_SINGLE_VECTOR
          input_array.begin();
#else
        input_array.begin()+vector_offsets[cell];
#endif

        // --------------------------------------------------------------------
        // apply tensor product kernels
        // z-direction
        //for (unsigned int i=0; i<dofs_per_cell; ++i)
        //  std::cout << input_ptr[i] << " ";
        //std::cout << std::endl;
        const unsigned int n_xy = Utilities::pow(degree+1,dim-1);
        const unsigned int n_chunks = (nn+VectorizedArray<Number>::n_array_elements-1)/VectorizedArray<Number>::n_array_elements;
        unsigned int offset2[VectorizedArray<Number>::n_array_elements];
        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
          offset2[v] = v*VectorizedArray<Number>::n_array_elements;
        VectorizedArray<Number> xo[n_chunks*VectorizedArray<Number>::n_array_elements];
        for (unsigned int i=nn; i<n_chunks*VectorizedArray<Number>::n_array_elements; ++i)
          xo[i] = VectorizedArray<Number>();
        const unsigned int z_offset = (n_xy + VectorizedArray<Number>::n_array_elements - 1) /
            VectorizedArray<Number>::n_array_elements * VectorizedArray<Number>::n_array_elements;
        //std::cout << "zoffset=" << z_offset << std::endl;
        for (unsigned int i1=0; i1<n_xy; i1+=VectorizedArray<Number>::n_array_elements)
          {
            VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (unsigned int i=0; i<mid; ++i)
              {
                VectorizedArray<Number> y1, y2;
                y1.load(input_ptr + i*n_xy + i1);
                y2.load(input_ptr + (nn-1-i)*n_xy + i1);
                xp[i] = y1 + y2;
                xm[i] = y1 - y2;
              }
            for (unsigned int col=0; col<mid; ++col)
              {
                VectorizedArray<Number> r0, r1;
                r0 = shape_vals[col]                 * xp[0];
                r1 = shape_vals[degree*offset + col] * xm[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_vals[ind*offset+col]          * xp[ind];
                    r1 += shape_vals[(degree-ind)*offset+col] * xm[ind];
                  }
                if (nn % 2 == 1)
                  {
                    VectorizedArray<Number> y1;
                    y1.load(input_ptr + mid*n_xy + i1);
                    r0 += shape_vals[mid*offset+col] * y1;
                  }
                xo[col] = r0 + r1;
                xo[nn-1-col] = r0 - r1;
              }
            if (nn % 2 == 1)
              xo[mid].load(input_ptr+mid*n_xy+i1);
            //for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
            //for (unsigned int i=0; i<nn; ++i)
            //  std::cout << xo[i][v] << " ";
            //std::cout << std::endl;

            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_transpose_and_store(false, VectorizedArray<Number>::n_array_elements,
                                             &xo[i*VectorizedArray<Number>::n_array_elements],
                                             offset2, &data_ptr[i1+i*z_offset][0]);
          }
        //for (unsigned int i=0; i<4; ++i)
        //  for (unsigned int v=0; v<2; ++v)
        //   std::cout << data_ptr[i][v] << " ";
        //std::cout << std::endl;

        for (unsigned int i2=0; i2<n_chunks; ++i2)
          {
            // x-direction
            VectorizedArray<Number> *__restrict in = data_ptr + i2*z_offset;
            for (unsigned int i1=0; i1<(dim>2?nn:1); ++i1)
              {
                VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    xp[i] = in[i+i1*nn] + in[nn-1-i+i1*nn];
                    xm[i] = in[i+i1*nn] - in[nn-1-i+i1*nn];
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0, r1;
                    r0 = shape_vals[col]                 * xp[0];
                    r1 = shape_vals[degree*offset + col] * xm[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      {
                        r0 += shape_vals[ind*offset+col]          * xp[ind];
                        r1 += shape_vals[(degree-ind)*offset+col] * xm[ind];
                      }
                    if (nn % 2 == 1)
                      r0 += shape_vals[mid*offset+col] * in[i1*nn+mid];

                    in[i1*nn+col]      = r0 + r1;
                    in[i1*nn+nn-1-col] = r0 - r1;
                  }
              }
            if (dim > 2)
              {
                // y-direction
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                    for (unsigned int i=0; i<mid; ++i)
                      {
                        xp[i] = in[i*nn+i1] + in[(nn-1-i)*nn+i1];
                        xm[i] = in[i*nn+i1] - in[(nn-1-i)*nn+i1];
                      }
                    for (unsigned int col=0; col<mid; ++col)
                      {
                        VectorizedArray<Number> r0, r1;
                        r0 = shape_vals[col]                 * xp[0];
                        r1 = shape_vals[degree*offset + col] * xm[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shape_vals[ind*offset+col]          * xp[ind];
                            r1 += shape_vals[(degree-ind)*offset+col] * xm[ind];
                          }
                        if (nn % 2 == 1)
                          r0 += shape_vals[mid*offset+col] * in[i1+mid*nn];

                        in[i1+col*nn]        = r0 + r1;
                        in[i1+(nn-1-col)*nn] = r0 - r1;
                      }
                  }
              }
          }
        //for (unsigned int v=0; v<nn; ++v)
        //  for (unsigned int i=0; i<n_xy; ++i)
        //    std::cout << data_ptr[i][v] << " ";
        //std::cout << std::endl;
        const VectorizedArray<Number> *__restrict jxw_ptr =
          jxw_data.begin() + data_offsets[cell];
        const VectorizedArray<Number> *__restrict jacobian_ptr =
          jacobian_data.begin() + data_offsets[cell]*dim*dim;
        if (is_cartesian)
          for (unsigned int d=0; d<dim; ++d)
            merged_array[d] = jxw_ptr[0] * jacobian_ptr[d] *
              jacobian_ptr[d];

        // z-derivative
        for (unsigned int i1=0; i1<n_xy; i1+=VectorizedArray<Number>::n_array_elements)
          {
            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_load_and_transpose(VectorizedArray<Number>::n_array_elements,
                                            &data_ptr[i1+i*z_offset][0], offset2,
                                            &xo[i*VectorizedArray<Number>::n_array_elements]);
            VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (unsigned int i=0; i<mid; ++i)
              {
                xp[i] = xo[i] - xo[nn-1-i];
                xm[i] = xo[i] + xo[nn-1-i];
              }
            for (unsigned int col=0; col<mid; ++col)
              {
                VectorizedArray<Number> r0, r1;
                r0 = shape_grads[col]                 * xp[0];
                r1 = shape_grads[degree*offset + col] * xm[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_grads[ind*offset+col]          * xp[ind];
                    r1 += shape_grads[(degree-ind)*offset+col] * xm[ind];
                  }
                if (nn % 2 == 1)
                  r1 += shape_grads[mid*offset+col] * xo[mid];

                xo[col] = r0 + r1;
                xo[nn-1-col] = r0 - r1;
              }
            if (nn%2==1)
              {
                VectorizedArray<Number> r0 = shape_grads[mid] * xp[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  r0 += shape_grads[ind*offset+mid] * xp[ind];
                xo[mid] = r0;
              }
            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_transpose_and_store(false, VectorizedArray<Number>::n_array_elements,
                                             &xo[i*VectorizedArray<Number>::n_array_elements],
                                             offset2,
                                             &data_ptr[dofs_per_cell+i1+i*z_offset][0]);
          }
          //for (unsigned int v=0; v<nn; ++v)
          //  for (unsigned int i=0; i<n_xy; ++i)
          // std::cout << data_ptr[dofs_per_cell+i][v] << " ";
        //std::cout << std::endl;

        // --------------------------------------------------------------------
        for (unsigned int i2=0; i2<n_chunks; ++i2)
          {
            VectorizedArray<Number> *__restrict in = data_ptr + i2*z_offset;
            VectorizedArray<Number> *__restrict outz = data_ptr + dofs_per_cell;
            VectorizedArray<Number> outy[nn*nn];
            if (dim>2)
              {
                // y-derivative
                for (unsigned int i1=0; i1<nn; ++i1) // loop over x layers
                  {
                    VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                    for (unsigned int i=0; i<mid; ++i)
                      {
                        xp[i] = in[i*nn+i1] - in[(nn-1-i)*nn+i1];
                        xm[i] = in[i*nn+i1] + in[(nn-1-i)*nn+i1];
                      }
                    for (unsigned int col=0; col<mid; ++col)
                      {
                        VectorizedArray<Number> r0, r1;
                        r0 = shape_grads[col]                 * xp[0];
                        r1 = shape_grads[degree*offset + col] * xm[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shape_grads[ind*offset+col]          * xp[ind];
                            r1 += shape_grads[(degree-ind)*offset+col] * xm[ind];
                          }
                        if (nn % 2 == 1)
                          r1 += shape_grads[mid*offset+col] * in[i1+mid*nn];

                        outy[i1+col*nn]        = r0 + r1;
                        outy[i1+(nn-1-col)*nn] = r0 - r1;
                      }
                    if (nn%2 == 1)
                      {
                        VectorizedArray<Number> r0 = shape_grads[mid] * xp[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          r0 += shape_grads[ind*offset+mid] * xp[ind];
                        outy[i1+nn*mid] = r0;
                      }
                  }
              }

            // x-derivative
            for (unsigned int i1=0; i1<(dim>2?nn:1); ++i1) // loop over y layers
              {
                VectorizedArray<Number> outx[nn];
                VectorizedArray<Number> inx[nn];
                VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    xp[i] = in[i+i1*nn] - in[nn-1-i+i1*nn];
                    xm[i] = in[i+i1*nn] + in[nn-1-i+i1*nn];
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0, r1;
                    r0 = shape_grads[col]                 * xp[0];
                    r1 = shape_grads[degree*offset + col] * xm[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      {
                        r0 += shape_grads[ind*offset+col]          * xp[ind];
                        r1 += shape_grads[(degree-ind)*offset+col] * xm[ind];
                      }
                    if (nn % 2 == 1)
                      r1 += shape_grads[mid*offset+col] * in[i1*nn+mid];

                    inx[col]      = r0 + r1;
                    inx[nn-1-col] = r0 - r1;
                  }
                if (nn%2 == 1)
                  {
                    VectorizedArray<Number> r0 = shape_grads[mid] * xp[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grads[ind*offset+mid] * xp[ind];
                    inx[mid] = r0;
                  }

                // operations on quadrature points
                // Cartesian cell case
                if (is_cartesian)
                  for (unsigned int i=0; i<nn; ++i)
                    {
                      //std::cout << "d ";
                      //std::cout << inx[i][0] << " " << inx[i][1] << " "
                      //          << outy[i1*nn+i][0] << " " << outy[i1*nn+i][1] << " "
                      //          << outz[i2*z_offset+i1*nn+i][0] << " " << outz[i2*z_offset+i1*nn+i][1] << std::endl;
                      const unsigned int q = (i2*nn+i1)*(dim>2?nn:1)+i;
                      const VectorizedArray<Number> weight = quadrature_weights_vect[q];
                      if (dim==2)
                        {
                          VectorizedArray<Number> t0 = inx[i]*merged_array[0] + outz[i2*z_offset+i]*merged_array[2];
                          VectorizedArray<Number> t1 = inx[i]*merged_array[2] + outz[i2*z_offset+i]*merged_array[1];
                          outx[i] = t0 * weight;
                          outz[i2*z_offset+i] = t1 * weight;
                        }
                      else if (dim==3)
                        {
                          VectorizedArray<Number> t0 = inx[i]*merged_array[0] + outy[i1*nn+i]*merged_array[3]+outz[i2*z_offset+i1*nn+i]*merged_array[4];
                          VectorizedArray<Number> t1 = inx[i]*merged_array[3] + outy[i1*nn+i]*merged_array[1]+outz[i2*z_offset+i1*nn+i]*merged_array[5];
                          VectorizedArray<Number> t2 = inx[i]*merged_array[4] + outy[i1*nn+i]*merged_array[5]+outz[i2*z_offset+i1*nn+i]*merged_array[2];
                          outx[i] = t0 * weight;
                          outy[i1*nn+i] = t1 * weight;
                          outz[i2*z_offset+i1*nn+i] = t2 * weight;
                        }
                    }
                else
                  throw;

                // x-derivative
                for (unsigned int i=0; i<mid; ++i)
                  {
                    xp[i] = outx[i] + outx[nn-1-i];
                    xm[i] = outx[i] - outx[nn-1-i];
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0, r1;
                    r0 = shape_grads[col*offset]          * xp[0];
                    r1 = shape_grads[(degree-col)*offset] * xm[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      {
                        r0 += shape_grads[col*offset+ind]          * xp[ind];
                        r1 += shape_grads[(degree-col)*offset+ind] * xm[ind];
                      }
                    if (nn % 2 == 1)
                      r0 += shape_grads[col*offset+mid] * outx[mid];

                    in[i1*nn+col]      = r0 + r1;
                    in[i1*nn+nn-1-col] = r1 - r0;
                  }
                if (nn%2 == 1)
                  {
                    VectorizedArray<Number> r0 = shape_grads[mid*offset] * xm[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      r0 += shape_grads[ind+mid*offset] * xm[ind];
                    in[i1*nn+mid] = r0;
                  }
              } // end of loop over y layers

            //std::cout << "y ";
            //for (unsigned int v=0; v<2; ++v)
            //for (unsigned i=0; i<4; ++i)
            //  std::cout << in[i][v] << " ";
            //std::cout << std::endl;

            if (dim>2)
              {
                // y-derivative
                for (unsigned int i1=0; i1<nn; ++i1) // loop over x layers
                  {
                    VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                    for (unsigned int i=0; i<mid; ++i)
                      {
                        xp[i] = outy[i*nn+i1] + outy[(nn-1-i)*nn+i1];
                        xm[i] = outy[i*nn+i1] - outy[(nn-1-i)*nn+i1];
                      }
                    for (unsigned int col=0; col<mid; ++col)
                      {
                        VectorizedArray<Number> r0, r1;
                        r0 = shape_grads[col*offset]          * xp[0];
                        r1 = shape_grads[(degree-col)*offset] * xm[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shape_grads[col*offset+ind]          * xp[ind];
                            r1 += shape_grads[(degree-col)*offset+ind] * xm[ind];
                          }
                        if (nn % 2 == 1)
                          r0 += shape_grads[col*offset+mid] * outy[i1+mid*nn];

                        in[i1+col*nn]        += r0 + r1;
                        in[i1+(nn-1-col)*nn] += r1 - r0;
                      }
                    if (nn%2 == 1)
                      {
                        VectorizedArray<Number> r0 = in[i1+nn*mid];
                        for (unsigned int ind=0; ind<mid; ++ind)
                          r0 += shape_grads[ind+mid*offset] * xm[ind];
                        in[i1+nn*mid] = r0;
                      }
                  }
              } // end of loop over z layers
          }
        //std::cout << "before z der: ";
        //for (unsigned int v=0; v<nn; ++v)
        //  for (unsigned int i=0; i<n_xy; ++i)
        //    std::cout << data_ptr[i][v] << " ";
        //std::cout << std::endl;

        // z-derivative
        for (unsigned int i1=0; i1<n_xy; i1+=VectorizedArray<Number>::n_array_elements)
          {
            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_load_and_transpose(VectorizedArray<Number>::n_array_elements,
                                            &data_ptr[dofs_per_cell+i1+i*z_offset][0], offset2,
                                            &xo[i*VectorizedArray<Number>::n_array_elements]);
            VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (unsigned int i=0; i<mid; ++i)
              {
                xp[i] = xo[i] + xo[nn-1-i];
                xm[i] = xo[i] - xo[nn-1-i];
              }
            for (unsigned int col=0; col<mid; ++col)
              {
                VectorizedArray<Number> r0, r1;
                r0 = shape_grads[col*offset]          * xp[0];
                r1 = shape_grads[(degree-col)*offset] * xm[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_grads[col*offset+ind]          * xp[ind];
                    r1 += shape_grads[(degree-col)*offset+ind] * xm[ind];
                  }
                if (nn % 2 == 1)
                  r0 += shape_grads[col*offset+mid] * xo[mid];

                xo[col] = r0 + r1;
                xo[nn-1-col] = r1 - r0;
              }
            if (nn%2==1)
              {
                VectorizedArray<Number> r0 = shape_grads[mid*offset] * xm[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  r0 += shape_grads[ind+mid*offset] * xm[ind];
                xo[mid] = r0;
              }
            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_transpose_and_store(true, VectorizedArray<Number>::n_array_elements,
                                             &xo[i*VectorizedArray<Number>::n_array_elements],
                                             offset2,
                                             &data_ptr[i1+i*z_offset][0]);
          }
        //std::cout << "after z der:  ";
        //for (unsigned int v=0; v<nn; ++v)
        //  for (unsigned int i=0; i<n_xy; ++i)
        //    std::cout << data_ptr[i][v] << " ";
        //std::cout << std::endl;

        for (unsigned int i2=0; i2<n_chunks; ++i2)
          {
            VectorizedArray<Number> *__restrict in = data_ptr + i2*z_offset;
            if (dim>2)
              {
                // y-direction
                for (unsigned int i1=0; i1<nn; ++i1)
                  {
                    VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                    for (unsigned int i=0; i<mid; ++i)
                      {
                        xp[i] = in[i*nn+i1] + in[(nn-1-i)*nn+i1];
                        xm[i] = in[i*nn+i1] - in[(nn-1-i)*nn+i1];
                      }
                    for (unsigned int col=0; col<mid; ++col)
                      {
                        VectorizedArray<Number> r0, r1;
                        r0 = shape_vals[col*offset]          * xp[0];
                        r1 = shape_vals[(degree-col)*offset] * xm[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shape_vals[col*offset+ind]          * xp[ind];
                            r1 += shape_vals[(degree-col)*offset+ind] * xm[ind];
                          }

                        in[i1+col*nn]        = r0 + r1;
                        in[i1+(nn-1-col)*nn] = r0 - r1;
                      }
                    if (nn%2==1)
                      {
                        // sum into because shape value is one in the middle
                        VectorizedArray<Number> r0 = in[i1+mid*nn];
                        for (unsigned int ind=0; ind<mid; ++ind)
                          r0 += shape_vals[ind+mid*offset] * xp[ind];
                        in[i1+mid*nn] = r0;
                      }
                  }
              }
            // x-direction
            for (unsigned int i1=0; i1<(dim>2?nn:1); ++i1)
              {
                VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                for (unsigned int i=0; i<mid; ++i)
                  {
                    xp[i] = in[i+i1*nn] + in[nn-1-i+i1*nn];
                    xm[i] = in[i+i1*nn] - in[nn-1-i+i1*nn];
                  }
                for (unsigned int col=0; col<mid; ++col)
                  {
                    VectorizedArray<Number> r0, r1;
                    r0 = shape_vals[col*offset]          * xp[0];
                    r1 = shape_vals[(degree-col)*offset] * xm[0];
                    for (unsigned int ind=1; ind<mid; ++ind)
                      {
                        r0 += shape_vals[col*offset+ind]          * xp[ind];
                        r1 += shape_vals[(degree-col)*offset+ind] * xm[ind];
                      }

                    in[col+i1*nn]        = r0 + r1;
                    in[(nn-1-col)+i1*nn] = r0 - r1;
                  }
                if (nn%2==1)
                  {
                    // sum into because shape value is one in the middle
                    VectorizedArray<Number> r0 = in[i1*nn+mid];
                    for (unsigned int ind=0; ind<mid; ++ind)
                      r0 += shape_vals[ind+mid*offset] * xp[ind];
                    in[i1*nn+mid] = r0;
                  }
              }
          }

        // z-direction
        Number *__restrict output_ptr =
#ifdef READ_SINGLE_VECTOR
          output_array.begin();
#else
        output_array.begin()+vector_offsets[cell];
#endif
        const unsigned int n_regular = (n_xy/VectorizedArray<Number>::n_array_elements)*VectorizedArray<Number>::n_array_elements;
        const unsigned int remainder = n_xy - n_regular;
        for (unsigned int i1=0; i1<n_xy; i1+=VectorizedArray<Number>::n_array_elements)
          {
            for (unsigned int i=0; i<n_chunks; ++i)
              vectorized_load_and_transpose(VectorizedArray<Number>::n_array_elements,
                                            &data_ptr[i1+i*z_offset][0], offset2,
                                            &xo[i*VectorizedArray<Number>::n_array_elements]);
            VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (unsigned int i=0; i<mid; ++i)
              {
                xp[i] = xo[i] + xo[nn-1-i];
                xm[i] = xo[i] - xo[nn-1-i];
              }
            for (unsigned int col=0; col<mid; ++col)
              {
                VectorizedArray<Number> r0, r1;
                r0 = shape_vals[col*offset]          * xp[0];
                r1 = shape_vals[(degree-col)*offset] * xm[0];
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_vals[col*offset+ind]          * xp[ind];
                    r1 += shape_vals[(degree-col)*offset+ind] * xm[ind];
                  }
                VectorizedArray<Number> t0 = r0+r1, t1=r0-r1;
                if (i1 < n_regular)
                  {
                    t0.store(output_ptr+i1+col*n_xy);
                    t1.store(output_ptr+i1+(nn-1-col)*n_xy);
                  }
                else
                  {
                    for (unsigned int i=0; i<remainder; ++i)
                      {
                        output_ptr[i1+i+col*n_xy] = t0[i];
                        output_ptr[i1+i+(nn-1-col)*n_xy] = t1[i];
                      }
                  }
              }
            if (nn%2==1)
              {
                // sum into because shape value is one in the middle
                VectorizedArray<Number> r0 = xo[mid];
                for (unsigned int ind=0; ind<mid; ++ind)
                  r0 += shape_vals[ind+mid*offset] * xp[ind];
                if (i1 < n_regular)
                  r0.store(output_ptr+i1+mid*n_xy);
                else
                  for (unsigned int i=0; i<remainder; ++i)
                    output_ptr[i1+i+mid*n_xy] = r0[i];
              }
          }
        //for (unsigned i=0; i<dofs_per_cell; ++i)
        //  std::cout << output_ptr[i] << " ";
        //std::cout << std::endl;
      }
#endif
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

    LagrangePolynomialBasis basis_gll(get_gauss_lobatto_points(degree+1));
    std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gll.value(i, gauss_points[q]);
          const double p2 = basis_gll.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values[i*stride+q] = 0.5 * (p1 + p2);
          shape_values[(degree-i)*stride+q] = 0.5 * (p1 - p2);
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
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_gradients[degree/2*stride+q] =
          basis_gauss.derivative(degree/2, gauss_points[q]);

    // get quadrature weights
    std::vector<double> gauss_weight_1d = get_gauss_weights(n_q_points_1d);
    quadrature_weights.resize(Utilities::pow(n_q_points_1d,dim));
    const unsigned int n_chunks = (n_q_points_1d+VectorizedArray<Number>::n_array_elements-1)/VectorizedArray<Number>::n_array_elements;
    quadrature_weights_vect.resize(n_chunks*Utilities::pow(n_q_points_1d,dim-1));
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
    for (unsigned int z=0; z<n_chunks; ++z)
      for (unsigned int i=0; i<Utilities::pow(n_q_points_1d,dim-1); ++i)
        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements &&
             z*VectorizedArray<Number>::n_array_elements+v<n_q_points_1d; ++v)
          quadrature_weights_vect[z*Utilities::pow(n_q_points_1d,dim-1)+i][v]
              = quadrature_weights[(z*VectorizedArray<Number>::n_array_elements+v)*
                                   Utilities::pow(n_q_points_1d,dim-1)+i];
  }

  AlignedVector<VectorizedArray<Number> > shape_values;
  AlignedVector<VectorizedArray<Number> > shape_gradients;

  AlignedVector<Number> quadrature_weights;
  AlignedVector<VectorizedArray<Number> > quadrature_weights_vect;

  AlignedVector<unsigned int> vector_offsets;
  AlignedVector<Number> input_array;
  AlignedVector<Number> output_array;

  AlignedVector<VectorizedArray<Number> > convection;

  AlignedVector<unsigned int> data_offsets;
  AlignedVector<VectorizedArray<Number> > jxw_data;
  AlignedVector<VectorizedArray<Number> > jacobian_data;
};


#endif
