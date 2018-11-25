// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// manual implementation of cell terms for Laplacian assuming interleaved
// storage, hardcoded vectorization
//
// Author: Martin Kronbichler, May 2017

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
class EvaluationCellLaplacian
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
        jxw_data.resize(n_element_batches*n_q_points);
        jacobian_data.resize(n_element_batches*n_q_points*dim*dim);
        for (unsigned int i=0; i<n_element_batches; ++i)
          {
            data_offsets[i] = i*n_q_points;
            for (unsigned int q=0; q<n_q_points; ++q)
              jxw_data[data_offsets[i]+q] =
                jacobian_determinant * quadrature_weights[q];
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int d=0; d<dim; ++d)
                jacobian_data[data_offsets[i]*dim*dim+
                              q*dim*dim+d*dim+d] = jacobian[d];
          }
      }

    for (unsigned int d=0; d<dim; ++d)
      {
        for (unsigned int e=0; e<dim; ++e)
          {
            normal_jac1[d][e] = VectorizedArray<Number>();
            normal_jac2[d][e] = VectorizedArray<Number>();
          }
        normal_jac1[d][dim-1] = jacobian[d];
        normal_jac2[d][dim-1] = jacobian[d];
        Number determinant = 1;
        for (unsigned int e=0; e<dim; ++e)
          if (d!=e)
            determinant *= jacobian[e];
        face_jxw[d] = 1./determinant;
      }

    convection.resize(dim*dofs_per_cell*(is_cartesian ? 1 : n_element_batches));
  }

  std::size_t n_elements() const
  {
    return VectorizedArray<Number>::n_array_elements * vector_offsets.size();
  }

  void do_verification()
  {
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
    if (degree < 1)
      return;

    // do not exchange data or zero out, assume DG operator does not need to
    // exchange and that the loop below takes care of the zeroing

    // global data structures
#ifdef COPY_ONLY_BENCHMARK
    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const VectorizedArray<Number> *__restrict input_ptr =
          input_array.begin()+vector_offsets[cell];
        VectorizedArray<Number> *__restrict output_ptr =
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
    const unsigned int nn_3d = dim==3 ? degree+1 : 1;
    const unsigned int mid = nn/2;
    const unsigned int offset = (nn+1)/2;

    constexpr unsigned int dofs_per_face = Utilities::pow(degree+1,dim-1);
    VectorizedArray<Number> tmp1[2*dofs_per_face], tmp2[dim*dofs_per_face], tmp3[dim*dofs_per_face], tmp4[dofs_per_face], tmp5[dofs_per_face];

    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        const VectorizedArray<Number> *__restrict input_ptr =
#ifdef READ_SINGLE_VECTOR
          input_array.begin();
#else
          input_array.begin()+vector_offsets[cell];
#endif

        // --------------------------------------------------------------------
        // apply tensor product kernels
        for (unsigned int i2=0; i2<nn_3d; ++i2)
          {
            // x-direction
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values, input_ptr+i1*nn, in+i1*nn);
              }
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values, in+i1, in+i1);
              }
            input_ptr += nn*nn;
          }

        // z direction
        if (dim == 3)
          for (unsigned int i1 = 0; i1<nn*nn_3d; ++i1)
            {
              apply_1d_matvec_kernel<nn, nn*nn, 0, true, false, VectorizedArray<Number>>
                (shape_values, data_ptr+i1, data_ptr+i1);
#ifdef DO_MASS_MATRIX
            }
        const VectorizedArray<Number> *__restrict jxw_ptr =
          jxw_data.begin() + data_offsets[cell];
        if (is_cartesian)
          for (unsigned int q=0; q<dofs_per_cell; ++q)
            data_ptr[q] *= jxw_ptr[0] * quadrature_weights[q];
        else
          for (unsigned int q=0; q<dofs_per_cell; ++q)
            data_ptr[q] *= jxw_ptr[q];
        if (dim == 3)
          for (unsigned int i1 = 0; i1<nn*nn_3d; ++i1)
            {
              VectorizedArray<Number> *__restrict out = data_ptr + i1;
              const unsigned int stride = nn*nn;
              VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
        /*
            const VectorizedArray<Number> *__restrict jxw_ptr =
              jxw_data.begin() + data_offsets[cell];
            if (is_cartesian)
              for (unsigned int q=0; q<degree+1; ++q)
                in[q*stride] *= jxw_ptr[0] * quadrature_weights[q*stride+i1];
            else
              for (unsigned int q=0; q<degree+1; ++q)
                in[q*stride] *= jxw_ptr[q*stride+i1];
            VectorizedArray<Number> *__restrict out = in;
        */
#else

#ifdef DO_CONVECTION
            }
        const VectorizedArray<Number> *__restrict convection_ptr =
          convection.begin() + (is_cartesian ? 0 : data_offsets[cell]*dim);
        for (unsigned int i2=0; i2<nn_3d; ++i2)  // loop over z layers
          {
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            VectorizedArray<Number> *__restrict outz = data_ptr + i2*nn*nn + dofs_per_cell;
            const VectorizedArray<Number> *__restrict conv = convection_ptr + i2*nn*nn*dim;
            VectorizedArray<Number> outy[nn*nn];
            for (unsigned int i1=0; i1<nn; ++i1) // loop over y layers
              {
                VectorizedArray<Number> outx[nn];
                VectorizedArray<Number> xp[mid>0?mid:1], xm[mid>0?mid:1];
                for (unsigned int i=0; i<nn; ++i)
                  {
                    outx[i] = in[i1*nn+i] * conv[dim*(i1*nn+i)];
                    outy[i1*nn+i] = in[i1*nn+i] * conv[dim*(i1*nn+i)+1];
                    if (dim == 3)
                      outz[i1*nn+i] = in[i1*nn+i] * conv[dim*(i1*nn+i)+2];
                  }
#else
              VectorizedArray<Number> *__restrict outz = data_ptr + i1 + dofs_per_cell;
              // z-derivative
              apply_1d_matvec_kernel<nn, nn*nn, 1, true, false, VectorizedArray<Number>>
                (shape_gradients, data_ptr+i1, outz);
            }

        // --------------------------------------------------------------------
        // mix with loop over quadrature points. depends on the data layout in
        // MappingInfo
        const VectorizedArray<Number> *__restrict jxw_ptr =
          jxw_data.begin() + data_offsets[cell];
        const VectorizedArray<Number> *__restrict jacobian_ptr =
          jacobian_data.begin() + data_offsets[cell]*dim*dim;
        if (is_cartesian)
          for (unsigned int d=0; d<dim; ++d)
            merged_array[d] = jxw_ptr[0] * jacobian_ptr[d] *
              jacobian_ptr[d];

        for (unsigned int i2=0; i2<nn_3d; ++i2)  // loop over z layers
          {
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            VectorizedArray<Number> *__restrict outz = data_ptr + dofs_per_cell;
            VectorizedArray<Number> outy[nn*nn];
            // y-derivative
            for (unsigned int i1=0; i1<nn; ++i1) // loop over x layers
              {
                apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients, in+i1, outy+i1);
              }

            // x-derivative
            for (unsigned int i1=0; i1<nn; ++i1) // loop over y layers
              {
                VectorizedArray<Number> outx[nn];
                apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients, in+i1*nn, outx);

                // operations on quadrature points
                // Cartesian cell case
                if (true || is_cartesian)
                  for (unsigned int i=0; i<nn; ++i)
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
                          VectorizedArray<Number> t0 = outy[i1*nn+i]*merged_array[3]+outz[i2*nn*nn+i1*nn+i]*merged_array[4] + outx[i]*merged_array[0];
                          VectorizedArray<Number> t1 = outy[i1*nn+i]*merged_array[1]+outz[i2*nn*nn+i1*nn+i]*merged_array[5] + outx[i]*merged_array[3];
                          VectorizedArray<Number> t2 = outy[i1*nn+i]*merged_array[5]+outz[i2*nn*nn+i1*nn+i]*merged_array[2] + outx[i]*merged_array[4];
                          outx[i] = t0 * weight;
                          outy[i1*nn+i] = t1 * weight;
                          outz[i2*nn*nn+i1*nn+i] = t2 * weight;
                        }
                    }
                else
                  for (unsigned int i=0; i<nn; ++i)
                    {
                      const unsigned int q=i2*nn*nn+i1*nn+i;
                      const VectorizedArray<Number> *geo_ptr = jacobian_ptr+q*dim*dim;
                      VectorizedArray<Number> grad[dim];
                      // get gradient
                      for (unsigned int d=0; d<dim; ++d)
                        {
                          grad[d] = geo_ptr[d*dim+1] * outy[i1*nn+i];
                          if (dim == 3)
                            grad[d] += geo_ptr[d*dim+2] * outz[q];
                          grad[d] += geo_ptr[d*dim] * outx[i];
                        }

                      // apply gradient of test function
                      if (dim == 3)
                        {
                          outx[i] = jxw_ptr[q] * (geo_ptr[0] * grad[0] + geo_ptr[dim] * grad[1]
                                                  + geo_ptr[2*dim] * grad[2]);
                          outy[i1*nn+i] = jxw_ptr[q] * (geo_ptr[1] * grad[0] +
                                                        geo_ptr[1+dim] * grad[1] +
                                                        geo_ptr[1+2*dim] * grad[2]);
                          outz[q] = jxw_ptr[q] * (geo_ptr[2] * grad[0] + geo_ptr[2+dim] * grad[1]
                                                  + geo_ptr[2+2*dim] * grad[2]);
                        }
                      else
                        {
                          outx[i] = jxw_ptr[q] * (geo_ptr[0] * grad[0] + geo_ptr[dim] * grad[1]);
                          outy[i1*nn+i] = jxw_ptr[q] * (geo_ptr[1] * grad[0] +
                                                        geo_ptr[1+dim] * grad[1]);
                        }
                    }
#endif // ifdef DO_CONVECTION else case


                // x-derivative
                apply_1d_matvec_kernel<nn, 1, 1, false, false, VectorizedArray<Number>>
                  (shape_gradients, outx, in + i1*nn);
              } // end of loop over y layers

            // y-derivative
            for (unsigned int i1=0; i1<nn; ++i1) // loop over x layers
              {
                apply_1d_matvec_kernel<nn, nn, 1, false, true, VectorizedArray<Number>>
                  (shape_gradients, outy+i1, in + i1, in + i1);
              }
          } // end of loop over z layers

        // z direction
        if (dim == 3)
          for (unsigned int i1 = 0; i1<nn*nn; ++i1)
            {
              // z-derivative
              VectorizedArray<Number> *__restrict inz = data_ptr + i1 + dofs_per_cell;
              VectorizedArray<Number> *__restrict out = data_ptr + i1;
              apply_1d_matvec_kernel<nn, nn*nn, 1, false, true, VectorizedArray<Number>>
                (shape_gradients, inz, out, out);
#endif // ifdef DO_MASS, else case

              // z-values
              apply_1d_matvec_kernel<nn, nn*nn, 0, false, false, VectorizedArray<Number>>
                  (shape_values, out, out);
            }

        VectorizedArray<Number> *__restrict output_ptr =
#ifdef READ_SINGLE_VECTOR
          output_array.begin();
#else
          output_array.begin()+vector_offsets[cell];
#endif
        for (unsigned int i2=0; i2<nn_3d; ++i2)
          {
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, nn, 0, false, false, VectorizedArray<Number>>
                  (shape_values, in+i1, in+i1);
              }
            // x-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                  (shape_values, in+i1*nn, output_ptr+i1*nn);
              }
            output_ptr += nn*nn;
          }

        input_ptr =
#ifdef READ_SINGLE_VECTOR
          input_array.begin();
#else
          input_array.begin()+vector_offsets[cell];
#endif
        output_ptr =
#ifdef READ_SINGLE_VECTOR
          output_array.begin();
#else
          output_array.begin()+vector_offsets[cell];
#endif
        for (unsigned int f=0; f<dim; ++f)
          {
            const unsigned int stride1 = Utilities::pow(degree+1,(f+1)%dim);
            const unsigned int stride2 = Utilities::pow(degree+1,(f+2)%dim);
            const unsigned int offset1 = Utilities::pow(degree+1,f%dim);
            if (dim == 2)
              {
                // right
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  {
                    tmp1[i1] = input_ptr[degree*offset1 + i1*stride1];
                    tmp1[i1+dofs_per_face] = hermite_derivative_on_face *
                      (input_ptr[(degree-1)*offset1+i1*stride1]-tmp1[i1]);
                  }
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values, tmp1+dofs_per_face, tmp2+dofs_per_face);
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values, tmp1, tmp4);
                apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients, tmp4, tmp2);

                // left
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  {
                    tmp1[i1] = input_ptr[i1*stride1];
                    tmp1[i1+dofs_per_face] = hermite_derivative_on_face *
                      (input_ptr[offset1+i1*stride1]-tmp1[i1]);
                  }
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values, tmp1+dofs_per_face, tmp3+dofs_per_face);
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values, tmp1, tmp5);
                apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                  (shape_gradients, tmp5, tmp3);
              }
            else
              {
                // right
                for (unsigned int i2=0; i2<degree+1; ++i2)
                  for (unsigned int i1=0; i1<degree+1; ++i1)
                    {
                      tmp1[i2*(degree+1)+i1] = input_ptr[degree*offset1 + i2*stride2 + i1*stride1];
                      tmp1[i2*(degree+1)+i1+dofs_per_face] = hermite_derivative_on_face *
                        (input_ptr[(degree-1)*offset1+i2*stride2+i1*stride1]-tmp1[i2*(degree+1)+i1]);
                    }
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+dofs_per_face+i1*(degree+1), tmp2+2*dofs_per_face+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp2+2*dofs_per_face+i1, tmp2+2*dofs_per_face+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+i1*(degree+1), tmp4+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp4+i1, tmp4+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients, tmp4+i1*(degree+1), tmp2+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients, tmp4+i1, tmp2+dofs_per_face+i1);

                // left
                for (unsigned int i2=0; i2<degree+1; ++i2)
                  for (unsigned int i1=0; i1<degree+1; ++i1)
                    {
                      tmp1[i2*(degree+1)+i1] = input_ptr[i2*stride2 + i1*stride1];
                      tmp1[i2*(degree+1)+i1+dofs_per_face] = hermite_derivative_on_face *
                        (input_ptr[offset1+i2*stride2+i1*stride1]-tmp1[i2*(degree+1)+i1]);
                    }
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+i1*(degree+1), tmp5+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp5+i1, tmp5+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients, tmp5+i1*(degree+1), tmp3+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, true, false, VectorizedArray<Number>>
                    (shape_gradients, tmp5+i1, tmp3+dofs_per_face+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+dofs_per_face+i1*(degree+1), tmp3+2*dofs_per_face+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp3+2*dofs_per_face+i1, tmp3+2*dofs_per_face+i1);
              }

            // quadrature point operation
            const VectorizedArray<Number> sigmaF = Number((degree+1)*(degree+1))*normal_jac1[f][dim-1];
            for (unsigned int q=0; q<Utilities::pow(degree+1,dim-1); ++q)
              {
                VectorizedArray<Number> average_valgrad = tmp2[q] * normal_jac1[f][0];
                for (unsigned int d=1; d<dim; ++d)
                  average_valgrad += tmp2[q+d*dofs_per_face] * normal_jac1[f][d];
                for (unsigned int d=0; d<dim; ++d)
                  average_valgrad += tmp3[q+d*dofs_per_face] * normal_jac2[f][d];
                VectorizedArray<Number> average_value = 0.5 * (tmp4[q] + tmp5[q]);
                const VectorizedArray<Number> weight = -face_quadrature_weight[q] * face_jxw[f];
                for (unsigned int d=0; d<dim; ++d)
                  tmp2[q+d*dofs_per_face] = weight * normal_jac1[f][d] * average_value;
                for (unsigned int d=0; d<dim; ++d)
                  tmp3[q+d*dofs_per_face] = weight * normal_jac2[f][d] * average_value;
                average_valgrad = weight * (average_value * 2. * sigmaF -
                                            average_valgrad * 0.5);
                tmp4[q] = -average_valgrad;
                tmp5[q] = average_valgrad;
              }

            if (dim==2)
              {
                // right
                apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                  (shape_gradients, tmp2, tmp4, tmp4);
                apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                  (shape_values, tmp4, tmp1);
                apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                  (shape_values, tmp2+dofs_per_face, tmp1+dofs_per_face);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  {
                    output_ptr[degree*offset1 + i1*stride1] += tmp1[i1] -
                      tmp1[i1+dofs_per_face] * hermite_derivative_on_face;
                    output_ptr[(degree-1)*offset1+i1*stride1] += hermite_derivative_on_face * tmp1[i1+dofs_per_face];
                  }

                // left
                apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                  (shape_gradients, tmp3, tmp5, tmp5);
                apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                  (shape_values, tmp5, tmp1);
                apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                  (shape_values, tmp3+dofs_per_face, tmp1+dofs_per_face);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  {
                    output_ptr[i1*stride1] += tmp1[i1] -
                      tmp1[i1+dofs_per_face] * hermite_derivative_on_face;
                    output_ptr[offset1+i1*stride1] += hermite_derivative_on_face * tmp1[i1+dofs_per_face];
                  }
              }
            else
              {
                // right
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients, tmp2+i1*(degree+1), tmp4+i1*(degree+1), tmp4+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients, tmp2+dofs_per_face+i1, tmp4+i1, tmp4+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                    (shape_values, tmp4+i1*(degree+1), tmp1+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, false, false, VectorizedArray<Number>>
                    (shape_values, tmp1+i1, tmp1+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp2+2*dofs_per_face+i1*(degree+1), tmp1+dofs_per_face+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+dofs_per_face+i1, tmp1+dofs_per_face+i1);
                for (unsigned int i2=0; i2<degree+1; ++i2)
                  for (unsigned int i1=0; i1<degree+1; ++i1)
                    {
                      output_ptr[degree*offset1 + i2*stride2+i1*stride1] += tmp1[i2*(degree+1)+i1] -
                      tmp1[i2*(degree+1)+i1+dofs_per_face] * hermite_derivative_on_face;
                      output_ptr[(degree-1)*offset1+i2*stride2+i1*stride1] += hermite_derivative_on_face * tmp1[i2*(degree+1)+i1+dofs_per_face];
                    }

                // left
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients, tmp3+i1*(degree+1), tmp5+i1*(degree+1), tmp5+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 1, false, true, VectorizedArray<Number>>
                    (shape_gradients, tmp3+dofs_per_face+i1, tmp5+i1, tmp5+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, false, false, VectorizedArray<Number>>
                    (shape_values, tmp5+i1*(degree+1), tmp1+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, false, false, VectorizedArray<Number>>
                    (shape_values, tmp1+i1, tmp1+i1);
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp3+2*dofs_per_face+i1*(degree+1), tmp1+dofs_per_face+i1*(degree+1));
                for (unsigned int i1=0; i1<degree+1; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values, tmp1+dofs_per_face+i1, tmp1+dofs_per_face+i1);
                for (unsigned int i2=0; i2<degree+1; ++i2)
                  for (unsigned int i1=0; i1<degree+1; ++i1)
                    {
                      output_ptr[i2*stride2+i1*stride1] += tmp1[i2*(degree+1)+i1] -
                      tmp1[i2*(degree+1)+i1+dofs_per_face] * hermite_derivative_on_face;
                      output_ptr[offset1+i2*stride2+i1*stride1] += hermite_derivative_on_face * tmp1[i2*(degree+1)+i1+dofs_per_face];
                    }
              }
          }
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

    HermiteLikePolynomialBasis basis(degree);
    std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis.value(i, gauss_points[q]);
          const double p2 = basis.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values[i*stride+q] = 0.5 * (p1 + p2);
          shape_values[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_values[degree/2*stride+q] =
          basis.value(degree/2, gauss_points[q]);

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

    hermite_derivative_on_face = basis.derivative(0, 0);
    if (std::abs(hermite_derivative_on_face[0] + basis.derivative(1, 0)) > 1e-12)
      std::cout << "Error, unexpected value of Hermite shape function derivative: "
                << hermite_derivative_on_face[0] << " vs "
                << basis.derivative(1, 0) << std::endl;

    // get quadrature weights
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

  AlignedVector<VectorizedArray<Number> > shape_values;
  AlignedVector<VectorizedArray<Number> > shape_gradients;

  std::array<std::array<VectorizedArray<Number>,dim>,dim> normal_jac1, normal_jac2;
  std::array<VectorizedArray<Number>,dim> face_jxw;

  AlignedVector<Number> quadrature_weights;
  AlignedVector<Number> face_quadrature_weight;

  VectorizedArray<Number> hermite_derivative_on_face;

  AlignedVector<unsigned int> vector_offsets;
  AlignedVector<VectorizedArray<Number> > input_array;
  AlignedVector<VectorizedArray<Number> > output_array;

  AlignedVector<VectorizedArray<Number> > convection;

  AlignedVector<unsigned int> data_offsets;
  AlignedVector<VectorizedArray<Number> > jxw_data;
  AlignedVector<VectorizedArray<Number> > jacobian_data;
};


#endif
