#include "gauss_formula.h"
#include "lagrange_polynomials.h"
#include "vectorization.h"
#include "aligned_vector.h"
#include "utilities.h"
#include "matrix_vector_kernel.h"

#include <iostream>
#include <iomanip>

#include <chrono>
#include <omp.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

const unsigned int min_degree = 1;
const unsigned int max_degree = 25;

const std::size_t  vector_size_guess = 10000;
const bool         cartesian         = true;

typedef double value_type;


namespace kernels
{
  /**
   * A namespace that contains optimized matrix-vector multiplication kernels
   * with some custom loop blocking to be used in combination with a second
   * dimension to yield a matrix-matrix multiplication.
   */
  namespace MatrixMultiplications
  {
    /**
     * A loop kernel that unrolls the inner loop 5 times and assumes two outer
     * streams to be given. This implementation needs 13 registers, 10 for the
     * accumulators, 2 for the temporary values f0 and f1 from the first input
     * (in0/in1) and 1 for the data from the second input (shape_data).
     */
    template <int dim, int mm, int nn, typename Number,
              int stride, bool dof_to_quad, bool add, bool do_both>
    inline
    void run_kernel_2x5(const Number *__restrict shape_data,
                        const Number *__restrict in0,
                        const Number *__restrict in1,
                        Number *__restrict out0,
                        Number *__restrict out1)
    {
      const int nn_5 = (nn/5)*5;
      for (int col=0; col<nn_5; col+=5)
        {
          Number val0, val1, val2, val3, val4, val5, val6, val7, val8, val9;
          if (dof_to_quad)
            {
              const Number f1 = in0[0];
              Number f2;
              if (do_both) f2 = in1[0];
              Number t = shape_data[col+0];
              val0 = t * f1;
              if (do_both) val1 = t * f2;
              t = shape_data[col+1];
              val2 = t * f1;
              if (do_both) val3 = t * f2;
              t = shape_data[col+2];
              val4 = t * f1;
              if (do_both) val5 = t * f2;
              t = shape_data[col+3];
              val6 = t * f1;
              if (do_both) val7 = t * f2;
              t = shape_data[col+4];
              val8 = t * f1;
              if (do_both) val9 = t * f2;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f1 = in0[stride*ind];
                  Number f2;
                  if (do_both) f2 = in1[stride*ind];
                  Number t = shape_data[ind*nn+col+0];
                  val0 += t * f1;
                  if (do_both) val1 += t * f2;
                  t = shape_data[ind*nn+col+1];
                  val2 += t * f1;
                  if (do_both) val3 += t * f2;
                  t = shape_data[ind*nn+col+2];
                  val4 += t * f1;
                  if (do_both) val5 += t * f2;
                  t = shape_data[ind*nn+col+3];
                  val6 += t * f1;
                  if (do_both) val7 += t * f2;
                  t = shape_data[ind*nn+col+4];
                  val8 += t * f1;
                  if (do_both) val9 += t * f2;
                }
            }
          else
            {
              const Number f1 = in0[0];
              Number f2;
              if (do_both) f2 = in1[0];
              Number t = shape_data[(col+0)*mm];
              val0 = t * f1;
              if (do_both) val1 = t * f2;
              t = shape_data[(col+1)*mm];
              val2 = t * f1;
              if (do_both) val3 = t * f2;
              t = shape_data[(col+2)*mm];
              val4 = t * f1;
              if (do_both) val5 = t * f2;
              t = shape_data[(col+3)*mm];
              val6 = t * f1;
              if (do_both) val7 = t * f2;
              t = shape_data[(col+4)*mm];
              val8 = t * f1;
              if (do_both) val9 = t * f2;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f1 = in0[stride*ind];
                  Number f2;
                  if (do_both) f2 = in1[stride*ind];
                  Number t = shape_data[(col+0)*mm+ind];
                  val0 += t * f1;
                  if (do_both) val1 += t * f2;
                  t = shape_data[(col+1)*mm+ind];
                  val2 += t * f1;
                  if (do_both) val3 += t * f2;
                  t = shape_data[(col+2)*mm+ind];
                  val4 += t * f1;
                  if (do_both) val5 += t * f2;
                  t = shape_data[(col+3)*mm+ind];
                  val6 += t * f1;
                  if (do_both) val7 += t * f2;
                  t = shape_data[(col+4)*mm+ind];
                  val8 += t * f1;
                  if (do_both) val9 += t * f2;
                }
            }
          if (add == false)
            {
              out0[stride*(col+0)]  = val0;
              out0[stride*(col+1)]  = val2;
              out0[stride*(col+2)]  = val4;
              out0[stride*(col+3)]  = val6;
              out0[stride*(col+4)]  = val8;
              if (do_both)
                {
                  out1[stride*(col+0)]  = val1;
                  out1[stride*(col+1)]  = val3;
                  out1[stride*(col+2)]  = val5;
                  out1[stride*(col+3)]  = val7;
                  out1[stride*(col+4)]  = val9;
                }
            }
          else
            {
              out0[stride*(col+0)] += val0;
              out0[stride*(col+1)] += val2;
              out0[stride*(col+2)] += val4;
              out0[stride*(col+3)] += val6;
              out0[stride*(col+4)] += val8;
              if (do_both)
                {
                  out1[stride*(col+0)] += val1;
                  out1[stride*(col+1)] += val3;
                  out1[stride*(col+2)] += val5;
                  out1[stride*(col+3)] += val7;
                  out1[stride*(col+4)] += val9;
                }
            }
        }
      Number val[8];
      const unsigned int remainder = nn - nn_5;
      if (remainder > 0)
        {
          if (dof_to_quad)
            {
              const Number f1 = in0[0];
              Number f2;
              if (do_both) f2 = in1[0];
              for (unsigned int i=0; i<remainder; ++i)
                {
                  const Number t = shape_data[nn_5+i];
                  val[2*i]   = t*f1;
                  if (do_both) val[2*i+1] = t*f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f1 = in0[stride*ind];
                  Number f2;
                  if (do_both) f2 = in1[stride*ind];
                  for (unsigned int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[ind*nn+nn_5+i];
                      val[2*i]   += t*f1;
                      if (do_both) val[2*i+1] += t*f2;
                    }
                }
            }
          else
            {
              const Number f1 = in0[0];
              Number f2;
              if (do_both) f2 = in1[0];
              for (unsigned int i=0; i<remainder; ++i)
                {
                  const Number t = shape_data[(nn_5+i)*mm];
                  val[2*i]   = t*f1;
                  if (do_both) val[2*i+1] = t*f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f1 = in0[stride*ind];
                  Number f2;
                  if (do_both) f2 = in1[stride*ind];
                  for (unsigned int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[(nn_5+i)*mm+ind];
                      val[2*i]   += t*f1;
                      if (do_both) val[2*i+1] += t*f2;
                    }
                }
            }
          if (add == false)
            for (unsigned int i=0; i<remainder; ++i)
              {
                out0[stride*(nn_5+i)]  = val[2*i];
                if (do_both) out1[stride*(nn_5+i)] = val[2*i+1];
              }
          else
            for (unsigned int i=0; i<remainder; ++i)
              {
                out0[stride*(nn_5+i)]  += val[2*i];
                if (do_both) out1[stride*(nn_5+i)] += val[2*i+1];
              }
        }
    }



    /**
     * A loop kernel that unrolls the inner loop 4 times and assumes three
     * outer streams to be given. Needs 16 registers, 12 for accumulators, 3
     * for the temporary values f0, f1, f2 from the first input (in0, in1, in2)
     * and 1 for the second input. This should fit into the 16 registers
     * available in x86-64.
     */
    template <int dim, int mm, int nn, typename Number,
              int stride, bool dof_to_quad, bool add, int n_streams>
    inline
    void run_kernel_3x4(const Number *__restrict shape_data,
                        const Number *__restrict in0,
                        const Number *__restrict in1,
                        const Number *__restrict in2,
                        Number *__restrict out0,
                        Number *__restrict out1,
                        Number *__restrict out2)
    {
      constexpr int nn_4 = (nn/4)*4;
      for (int col=0; col<nn_4; col+=4)
        {
          Number val[12];
          if (dof_to_quad)
            {
              const Number f0 = in0[0];
              Number f1, f2;
              if (n_streams>1) f1 = in1[0];
              if (n_streams>2) f2 = in2[0];
              for (int i=0; i<4; ++i)
                {
                  const Number t = shape_data[col+i];
                  val[3*i] = t * f0;
                  if (n_streams>1) val[3*i+1] = t * f1;
                  if (n_streams>2) val[3*i+2] = t * f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f0 = in0[stride*ind];
                  Number f1, f2;
                  if (n_streams>1) f1 = in1[stride*ind];
                  if (n_streams>2) f2 = in2[stride*ind];
                  for (int i=0; i<4; ++i)
                    {
                      const Number t = shape_data[ind*nn+col+i];
                      val[3*i] += t * f0;
                      if (n_streams>1) val[3*i+1] += t * f1;
                      if (n_streams>2) val[3*i+2] += t * f2;
                    }
                }
            }
          else
            {
              const Number f0 = in0[0];
              Number f1, f2;
              if (n_streams>1) f1 = in1[0];
              if (n_streams>2) f2 = in2[0];
              for (int i=0; i<4; ++i)
                {
                  const Number t = shape_data[(col+i)*mm];
                  val[3*i] = t * f0;
                  if (n_streams>1) val[3*i+1] = t * f1;
                  if (n_streams>2) val[3*i+2] = t * f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f0 = in0[stride*ind];
                  Number f1, f2;
                  if (n_streams>1) f1 = in1[stride*ind];
                  if (n_streams>2) f2 = in2[stride*ind];
                  for (int i=0; i<4; ++i)
                    {
                      const Number t = shape_data[(col+i)*mm+ind];
                      val[3*i] += t * f0;
                      if (n_streams>1) val[3*i+1] += t * f1;
                      if (n_streams>2) val[3*i+2] += t * f2;
                    }
                }
            }
          if (add == false)
            {
              for (unsigned int i=0; i<4; ++i)
                out0[stride*(col+i)]  = val[3*i];
              if (n_streams > 1)
                for (unsigned int i=0; i<4; ++i)
                  out1[stride*(col+i)]  = val[3*i+1];
              if (n_streams > 2)
                for (unsigned int i=0; i<4; ++i)
                  out2[stride*(col+i)]  = val[3*i+2];
            }
          else
            {
              for (unsigned int i=0; i<4; ++i)
                out0[stride*(col+i)] += val[3*i];
              if (n_streams > 1)
                for (unsigned int i=0; i<4; ++i)
                  out1[stride*(col+i)] += val[3*i+1];
              if (n_streams > 2)
                for (unsigned int i=0; i<4; ++i)
                  out2[stride*(col+i)] += val[3*i+2];
            }
        }
      Number val[9];
      constexpr int remainder = nn - nn_4;
      if (remainder > 0)
        {
          if (dof_to_quad)
            {
              const Number f0 = in0[0];
              Number f1, f2;
              if (n_streams>1) f1 = in1[0];
              if (n_streams>2) f2 = in2[0];
              for (int i=0; i<remainder; ++i)
                {
                  const Number t = shape_data[nn_4+i];
                  val[3*i] = t * f0;
                  if (n_streams>1) val[3*i+1] = t * f1;
                  if (n_streams>2) val[3*i+2] = t * f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f0 = in0[stride*ind];
                  Number f1, f2;
                  if (n_streams>1) f1 = in1[stride*ind];
                  if (n_streams>2) f2 = in2[stride*ind];
                  for (int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[ind*nn+nn_4+i];
                      val[3*i] += t * f0;
                      if (n_streams>1) val[3*i+1] += t * f1;
                      if (n_streams>2) val[3*i+2] += t * f2;
                    }
                }
            }
          else
            {
              const Number f0 = in0[0];
              Number f1, f2;
              if (n_streams>1) f1 = in1[0];
              if (n_streams>2) f2 = in2[0];
              for (int i=0; i<remainder; ++i)
                {
                  const Number t = shape_data[(nn_4+i)*mm];
                  val[3*i] = t * f0;
                  if (n_streams>1) val[3*i+1] = t * f1;
                  if (n_streams>2) val[3*i+2] = t * f2;
                }
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f0 = in0[stride*ind];
                  Number f1, f2;
                  if (n_streams>1) f1 = in1[stride*ind];
                  if (n_streams>2) f2 = in2[stride*ind];
                  for (int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[(nn_4+i)*mm+ind];
                      val[3*i] += t * f0;
                      if (n_streams>1) val[3*i+1] += t * f1;
                      if (n_streams>2) val[3*i+2] += t * f2;
                    }
                }
            }
          if (add == false)
            {
              for (int i=0; i<remainder; ++i)
                out0[stride*(nn_4+i)]  = val[3*i];
              if (n_streams > 1)
                for (int i=0; i<remainder; ++i)
                  out1[stride*(nn_4+i)]  = val[3*i+1];
              if (n_streams > 2)
                for (int i=0; i<remainder; ++i)
                  out2[stride*(nn_4+i)]  = val[3*i+2];
            }
          else
            {
              for (int i=0; i<remainder; ++i)
                out0[stride*(nn_4+i)] += val[3*i];
              if (n_streams > 1)
                for (int i=0; i<remainder; ++i)
                  out1[stride*(nn_4+i)] += val[3*i+1];
              if (n_streams > 2)
                for (int i=0; i<remainder; ++i)
                  out2[stride*(nn_4+i)] += val[3*i+2];
            }
        }
    }



    /**
     * A loop kernel that unrolls the inner loop 8 times and assumes only a
     * single input stream to be given. Needs 10 registers, 8 for the
     * accumulators, 1 for the first input data `in` and 1 for the second input
     * data shape_data.
     */
    template <int dim, int mm, int nn, typename Number,
              int stride, bool dof_to_quad, bool add>
    inline
    void run_kernel_1x8(const Number *__restrict shape_data,
                        const Number *__restrict in,
                        Number *__restrict out)
    {
      constexpr int nn_8 = (nn / 8) * 8;
      for (int col=0; col<nn_8; col+=8)
        {
          Number val0, val1, val2, val3, val4, val5, val6, val7;
          if (dof_to_quad)
            {
              const Number f = in[0];
              val0 = shape_data[col+0] * f;
              val1 = shape_data[col+1] * f;
              val2 = shape_data[col+2] * f;
              val3 = shape_data[col+3] * f;
              val4 = shape_data[col+4] * f;
              val5 = shape_data[col+5] * f;
              val6 = shape_data[col+6] * f;
              val7 = shape_data[col+7] * f;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f = in[stride*ind];
                  val0 += shape_data[ind*nn+col+0] * f;
                  val1 += shape_data[ind*nn+col+1] * f;
                  val2 += shape_data[ind*nn+col+2] * f;
                  val3 += shape_data[ind*nn+col+3] * f;
                  val4 += shape_data[ind*nn+col+4] * f;
                  val5 += shape_data[ind*nn+col+5] * f;
                  val6 += shape_data[ind*nn+col+6] * f;
                  val7 += shape_data[ind*nn+col+7] * f;
                }
            }
          else
            {
              const Number f = in[0];
              val0 = shape_data[(col+0)*mm] * f;
              val1 = shape_data[(col+1)*mm] * f;
              val2 = shape_data[(col+2)*mm] * f;
              val3 = shape_data[(col+3)*mm] * f;
              val4 = shape_data[(col+4)*mm] * f;
              val5 = shape_data[(col+5)*mm] * f;
              val6 = shape_data[(col+6)*mm] * f;
              val7 = shape_data[(col+7)*mm] * f;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f = in[stride*ind];
                  val0 += shape_data[(col+0)*mm+ind] * f;
                  val1 += shape_data[(col+1)*mm+ind] * f;
                  val2 += shape_data[(col+2)*mm+ind] * f;
                  val3 += shape_data[(col+3)*mm+ind] * f;
                  val4 += shape_data[(col+4)*mm+ind] * f;
                  val5 += shape_data[(col+5)*mm+ind] * f;
                  val6 += shape_data[(col+6)*mm+ind] * f;
                  val7 += shape_data[(col+7)*mm+ind] * f;
                }
            }
          if (add == false)
            {
              out[stride*(col+0)]  = val0;
              out[stride*(col+1)]  = val1;
              out[stride*(col+2)]  = val2;
              out[stride*(col+3)]  = val3;
              out[stride*(col+4)]  = val4;
              out[stride*(col+5)]  = val5;
              out[stride*(col+6)]  = val6;
              out[stride*(col+7)]  = val7;
            }
          else
            {
              out[stride*(col+0)] += val0;
              out[stride*(col+1)] += val1;
              out[stride*(col+2)] += val2;
              out[stride*(col+3)] += val3;
              out[stride*(col+4)] += val4;
              out[stride*(col+5)] += val5;
              out[stride*(col+6)] += val6;
              out[stride*(col+7)] += val7;
            }
        }
      Number val0, val1, val2, val3, val4, val5, val6;
      constexpr unsigned int remainder = nn - nn_8;
      if (remainder > 0)
        {
          if (dof_to_quad)
            {
              const Number f = in[0];
              val0 = shape_data[nn_8+0] * f;
              if (remainder > 1) val1 = shape_data[nn_8+1] * f;
              if (remainder > 2) val2 = shape_data[nn_8+2] * f;
              if (remainder > 3) val3 = shape_data[nn_8+3] * f;
              if (remainder > 4) val4 = shape_data[nn_8+4] * f;
              if (remainder > 5) val5 = shape_data[nn_8+5] * f;
              if (remainder > 6) val6 = shape_data[nn_8+6] * f;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f = in[stride*ind];
                  val0 += shape_data[ind*nn+nn_8+0] * f;
                  if (remainder > 1) val1 += shape_data[ind*nn+nn_8+1] * f;
                  if (remainder > 2) val2 += shape_data[ind*nn+nn_8+2] * f;
                  if (remainder > 3) val3 += shape_data[ind*nn+nn_8+3] * f;
                  if (remainder > 4) val4 += shape_data[ind*nn+nn_8+4] * f;
                  if (remainder > 5) val5 += shape_data[ind*nn+nn_8+5] * f;
                  if (remainder > 6) val6 += shape_data[ind*nn+nn_8+6] * f;
                }
            }
          else
            {
              const Number f = in[0];
              val0 = shape_data[(nn_8+0)*mm] * f;
              if (remainder > 1) val1 = shape_data[(nn_8+1)*mm] * f;
              if (remainder > 2) val2 = shape_data[(nn_8+2)*mm] * f;
              if (remainder > 3) val3 = shape_data[(nn_8+3)*mm] * f;
              if (remainder > 4) val4 = shape_data[(nn_8+4)*mm] * f;
              if (remainder > 5) val5 = shape_data[(nn_8+5)*mm] * f;
              if (remainder > 6) val6 = shape_data[(nn_8+6)*mm] * f;
              for (int ind=1; ind<mm; ++ind)
                {
                  const Number f = in[stride*ind];
                  val0 += shape_data[(nn_8+0)*mm+ind] * f;
                  if (remainder > 1) val1 += shape_data[(nn_8+1)*mm+ind] * f;
                  if (remainder > 2) val2 += shape_data[(nn_8+2)*mm+ind] * f;
                  if (remainder > 3) val3 += shape_data[(nn_8+3)*mm+ind] * f;
                  if (remainder > 4) val4 += shape_data[(nn_8+4)*mm+ind] * f;
                  if (remainder > 5) val5 += shape_data[(nn_8+5)*mm+ind] * f;
                  if (remainder > 6) val6 += shape_data[(nn_8+6)*mm+ind] * f;
                }
            }
          if (add == false)
            {
              out[stride*(nn_8+0)]  = val0;
              if (remainder > 1) out[stride*(nn_8+1)]  = val1;
              if (remainder > 2) out[stride*(nn_8+2)]  = val2;
              if (remainder > 3) out[stride*(nn_8+3)]  = val3;
              if (remainder > 4) out[stride*(nn_8+4)]  = val4;
              if (remainder > 5) out[stride*(nn_8+5)]  = val5;
              if (remainder > 6) out[stride*(nn_8+6)]  = val6;
            }
          else
            {
              out[stride*(nn_8+0)] += val0;
              if (remainder > 1) out[stride*(nn_8+1)] += val1;
              if (remainder > 2) out[stride*(nn_8+2)] += val2;
              if (remainder > 3) out[stride*(nn_8+3)] += val3;
              if (remainder > 4) out[stride*(nn_8+4)] += val4;
              if (remainder > 5) out[stride*(nn_8+5)] += val5;
              if (remainder > 6) out[stride*(nn_8+6)] += val6;
            }
        }
    }

  } // end of namespace MatrixMultiplications



  template <int dim, int n_functions, int n_points, typename Number, int kernel_no>
  struct EvaluatorTensorProduct
  {
    static const unsigned int dofs_per_cell = Utilities::pow(n_functions,dim);
    static const unsigned int n_q_points = Utilities::pow(n_points,dim);

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const unsigned int           dummy1 = 0,
                            const unsigned int           dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    /**
     * @tparam direction Direction that is evaluated
     * @tparam dof_to_quad If true, the input vector is on the dofs and the
     *                    output on the quadrature points, if false viceversa
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     *
     * @param shape_data Unit cell information
     * @param in address of the input data vector
     * @param out address of the output data vector
     */
    template <int direction, bool dof_to_quad, bool add>
    static void apply (const Number *__restrict shape_data,
                       const Number *__restrict in,
                       Number       *__restrict out);

    template <int face_direction, bool dof_to_quad, bool add, int max_derivative>
    void apply_face (const Number *__restrict in,
                     Number       *__restrict out) const;

    const Number *shape_values;
    const Number *shape_gradients;
  };




  template <int dim, int n_functions, int n_points, typename Number, int kernel_no>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<dim,n_functions,n_points,Number,kernel_no>
  ::apply (const Number *__restrict shape_data,
           const Number *__restrict in,
           Number       *__restrict out)
  {
    const int mm     = dof_to_quad ? n_functions : n_points,
              nn     = dof_to_quad ? n_points : n_functions;

    const int stride    = Utilities::pow(n_points,direction);
    const int jump_in   = stride * (mm-1);
    const int jump_out  = stride * (nn-1);
    const int n_blocks1 = stride;
    const int n_blocks2 = Utilities::pow(n_functions,(direction>=dim)?0:(dim-direction-1));

    if (kernel_no == 0)
      {
        for (int i2=0; i2<n_blocks2; ++i2)
          {
            for (int i1=0; i1<n_blocks1; ++i1)
              {
                for (int col=0; col<nn; ++col)
                  {
                    Number val0;
                    if (dof_to_quad == true)
                      val0 = shape_data[col];
                    else
                      val0 = shape_data[col*n_points];
                    Number res0 = val0 * in[0];
                    for (int ind=1; ind<mm; ++ind)
                      {
                        if (dof_to_quad == true)
                          val0 = shape_data[ind*n_points+col];
                        else
                          val0 = shape_data[col*n_points+ind];
                        res0 += val0 * in[stride*ind];
                      }
                    if (add == false)
                      out[stride*col]  = res0;
                    else
                      out[stride*col] += res0;
                  }

                ++in;
                ++out;
              }
            in += stride*(mm-1);
            out += stride*(nn-1);
          }
      }
    else if (kernel_no == 1)
      for (int i2=0; i2<n_blocks2; ++i2)
        {
          for (int i1=0; i1<n_blocks1; ++i1)
            {
              MatrixMultiplications::run_kernel_1x8<dim,mm,nn,Number,stride,dof_to_quad,add>
                (shape_data, in, out);
              ++in;
              ++out;
            }
          in += jump_in;
          out+= jump_out;
        }
    else if (kernel_no == 2)
      {
        const int n_blocks1a = direction > 0 ? n_blocks1/2 : n_blocks1;
        const int n_blocks2a = direction > 0 ? n_blocks2   : n_blocks2/2;

        for (int i2=0; i2<n_blocks2a; ++i2)
          {
            for (int i1=0; i1<n_blocks1a; ++i1)
              {
                MatrixMultiplications::run_kernel_2x5<dim,mm,nn,Number,stride,dof_to_quad,add,true>
                  (shape_data, in, in +(direction>0?1:mm),
                   out, out+(direction>0?1:nn));
                in  += (direction>0 ? 2 : 1);
                out += (direction>0 ? 2 : 1);
              }
            if (direction > 0 && n_blocks1%2 == 1)
              {
                MatrixMultiplications::run_kernel_2x5<dim,mm,nn,Number,stride,dof_to_quad,add,false>
                  (shape_data, in, nullptr, out, nullptr);
                in += 1;
                out += 1;
              }
            in += (direction>0 ? jump_in : (jump_in+mm));
            out+= (direction>0 ? jump_out : (jump_out+nn));
          }
        if (direction == 0 && n_blocks2%2 == 1)
          {
            MatrixMultiplications::run_kernel_2x5<dim,mm,nn,Number,stride,dof_to_quad,add,false>
              (shape_data, in, nullptr, out, nullptr);
          }
      }
    else if (kernel_no == 3)
      {
        const int n_blocks1a = direction > 0 ? n_blocks1/3 : n_blocks1;
        const int n_blocks2a = direction > 0 ? n_blocks2   : n_blocks2/3;

        for (int i2=0; i2<n_blocks2a; ++i2)
          {
            for (int i1=0; i1<n_blocks1a; ++i1)
              {
                MatrixMultiplications::run_kernel_3x4<dim,mm,nn,Number,stride,dof_to_quad,add,3>
                  (shape_data, in, in +(direction>0?1:mm),in +2*(direction>0?1:mm),
                   out, out+(direction>0?1:nn), out+2*(direction>0?1:nn));
                in  += (direction>0 ? 3 : 1);
                out += (direction>0 ? 3 : 1);
              }
            if (direction > 0 && n_blocks1%3 != 0)
              {
                MatrixMultiplications::run_kernel_3x4<dim,mm,nn,Number,stride,dof_to_quad,add,n_blocks1%3>
                  (shape_data, in, in +(direction>0?1:mm),in +2*(direction>0?1:mm),
                   out, out+(direction>0?1:nn), out+2*(direction>0?1:nn));
                in += n_blocks1%3;
                out += n_blocks1%3;
              }
            in += (direction>0 ? jump_in : (jump_in+2*mm));
            out+= (direction>0 ? jump_out : (jump_out+2*nn));
          }
        if (direction == 0 && n_blocks2%3 != 0)
          {
            MatrixMultiplications::run_kernel_3x4<dim,mm,nn,Number,stride,dof_to_quad,add,n_blocks2%3>
              (shape_data, in, in+mm, nullptr, out, out+nn, nullptr);
          }
      }
    else
      throw;
  }



  template <int dim, typename Number,int kernel>
  struct EvaluatorTensorProduct<dim,0,0,Number,kernel>
  {
    static const unsigned int dofs_per_cell = static_cast<unsigned int>(-1);
    static const unsigned int n_q_points = static_cast<unsigned int>(-1);

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other constructor
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      n_functions (-1),
      n_points (-1)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const unsigned int           n_functions,
                            const unsigned int           n_points)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      n_functions (n_functions),
      n_points (n_points)
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number *in,
            Number       *out) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number *in,
               Number       *out) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void apply (const Number *shape_data,
                const Number *in,
                Number       *out) const;

    template <int face_direction, bool dof_to_quad, bool add, int max_derivative>
    void apply_face (const Number *__restrict in,
                     Number       *__restrict out) const;

    const Number *shape_values;
    const Number *shape_gradients;
    const unsigned int n_functions;
    const unsigned int n_points;
  };

  template <int dim, typename Number, int kernel>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<dim,0,0,Number,kernel>
  ::apply (const Number *shape_data,
           const Number *in,
           Number       *out) const
  {
    const int mm     = dof_to_quad ? n_functions : n_points,
              nn     = dof_to_quad ? n_points : n_functions;

    const int stride    = direction==0 ? 1 : Utilities::pow(n_points, direction);
    const int n_blocks1 = stride;
    const int n_blocks2 = direction >= dim-1 ? 1 : Utilities::pow(n_functions,dim-direction-1);

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<nn; ++col)
              {
                Number val0;
                if (dof_to_quad == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_points];
                Number res0 = val0 * in[0];
                for (int ind=1; ind<mm; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_data[ind*n_points+col];
                    else
                      val0 = shape_data[col*n_points+ind];
                    res0 += val0 * in[stride*ind];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            ++in;
            ++out;
          }
        in += stride*(mm-1);
        out += stride*(nn-1);
      }
  }


  template <int dim, int n_functions, int n_points, typename Number>
  struct EvaluatorTensorProduct<dim,n_functions,n_points,Number,4>
  {
    static const unsigned int dofs_per_cell = Utilities::pow(n_functions,dim);
    static const unsigned int n_q_points = Utilities::pow(n_points,dim);

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other constructor
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const unsigned int           dummy1 = 0,
                            const unsigned int           dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number *in,
            Number       *out) const
    {
      apply<direction,dof_to_quad,add,0,false>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number *in,
               Number       *out) const
    {
      apply<direction,dof_to_quad,add,1,false>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add, int type, bool one_line>
    void apply (const Number *shape_data,
                const Number *in,
                Number       *out) const;

    const Number *shape_values;
    const Number *shape_gradients;
  };

  template <int dim, int n_rows, int n_columns, typename Number>
  template <int direction, bool contract_over_rows, bool add, int type, bool one_line>
  inline
  void
  EvaluatorTensorProduct<dim,n_rows,n_columns,Number,4>
  ::apply (const Number *shapes,
           const Number *in,
           Number       *out) const
  {
    static_assert(type < 3, "Only three variants type=0,1,2 implemented");
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");

    // We cannot statically assert that direction is less than dim, so must do
    // an additional dynamic check

    constexpr int nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int mm     = contract_over_rows ? n_rows : n_columns;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    constexpr int offset = (n_columns + 1) / 2;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            Number xp[mid > 0 ? mid : 1], xm[mid > 0 ? mid : 1];
            for (int i = 0; i < mid; ++i)
              {
                if (contract_over_rows == true && type == 1)
                  {
                    xp[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                    xm[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                  }
                else
                  {
                    xp[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                    xm[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                  }
              }
            Number xmid = in[stride * mid];
            for (int col = 0; col < n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (contract_over_rows == true)
                      {
                        r0 = shapes[col] * xp[0];
                        r1 = shapes[(n_rows - 1) * offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col * offset] * xp[0];
                        r1 = shapes[(n_rows - 1 - col) * offset] * xm[0];
                      }
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            r0 += shapes[ind * offset + col] * xp[ind];
                            r1 += shapes[(n_rows - 1 - ind) * offset + col] *
                                  xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col * offset + ind] * xp[ind];
                            r1 += shapes[(n_rows - 1 - col) * offset + ind] *
                                  xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && contract_over_rows == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid * offset + col] * xmid;
                    else
                      r0 += shapes[mid * offset + col] * xmid;
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0))
                  r0 += shapes[col * offset + mid] * xmid;

                if (add == false)
                  {
                    out[stride * col] = r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride * (nn - 1 - col)] = r1 - r0;
                    else
                      out[stride * (nn - 1 - col)] = r0 - r1;
                  }
                else
                  {
                    out[stride * col] += r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride * (nn - 1 - col)] += r1 - r0;
                    else
                      out[stride * (nn - 1 - col)] += r0 - r1;
                  }
              }
            if (type == 0 && contract_over_rows == true && nn % 2 == 1 &&
                mm % 2 == 1)
              {
                if (add == false)
                  out[stride * n_cols] = shapes[mid * offset + n_cols] * xmid;
                else
                  out[stride * n_cols] += shapes[mid * offset + n_cols] * xmid;
              }
            else if (contract_over_rows == true && nn % 2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    r0 = shapes[n_cols] * xp[0];
                    for (int ind = 1; ind < mid; ++ind)
                      r0 += shapes[ind * offset + n_cols] * xp[ind];
                  }
                else
                  r0 = Number();
                if (type != 1 && mm % 2 == 1)
                  r0 += shapes[mid * offset + n_cols] * xmid;

                if (add == false)
                  out[stride * n_cols] = r0;
                else
                  out[stride * n_cols] += r0;
              }
            else if (contract_over_rows == false && nn % 2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    if (type == 1)
                      {
                        r0 = shapes[n_cols * offset] * xm[0];
                        for (int ind = 1; ind < mid; ++ind)
                          r0 += shapes[n_cols * offset + ind] * xm[ind];
                      }
                    else
                      {
                        r0 = shapes[n_cols * offset] * xp[0];
                        for (int ind = 1; ind < mid; ++ind)
                          r0 += shapes[n_cols * offset + ind] * xp[ind];
                      }
                  }
                else
                  r0 = Number();

                if ((type == 0 || type == 2) && mm % 2 == 1)
                  r0 += shapes[n_cols * offset + mid] * xmid;

                if (add == false)
                  out[stride * n_cols] = r0;
                else
                  out[stride * n_cols] += r0;
              }
            if (one_line == false)
              {
                in += 1;
                out += 1;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm - 1);
            out += stride * (nn - 1);
          }
      }
  }
}



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

    runtime_degree = degree;
  }

  std::size_t n_elements() const
  {
    return VectorizedArray<Number>::n_array_elements * vector_offsets.size();
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
    AlignedVector<VectorizedArray<Number> > boundary_integral(Utilities::pow(degree+1,dim-1));
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
            VectorizedArray<Number> *data_ptr = &input_array[vector_offsets[i]];
            indices[2] = 0;
            for (unsigned int p=0; indices[2]<(dim>2?degree+1:1); ++indices[2])
              for (indices[1]=0; indices[1]<(dim>1?degree+1:1); ++indices[1])
                for (indices[0]=0; indices[0]<degree+1; ++indices[0], ++p)
                  data_ptr[p] = values_1d[indices[test]];
          }

        matrix_vector_product<4>();

        // remove the boundary integral from the cell integrals and check the
        // error
        double boundary_factor = 1.;
        for (unsigned int d=0; d<dim; ++d)
          if (d!=test)
            boundary_factor /= jacobian[d];
        double max_error = 0;
        unsigned int cell=0;
          {
            VectorizedArray<Number> *data_ptr = &output_array[vector_offsets[cell]];
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
              for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
                max_error = std::max(max_error, (double)data_ptr[i][v]);
          }

          std::cout << "Error of integral in direction " << test << ": "
                    << max_error << std::endl;
      }
  }

  template <int kernel>
  void matrix_vector_product()
  {
    if (degree < 1)
      return;

    VectorizedArray<Number> *scratch = scratch_array.begin();
    const unsigned int n_q_points_1d = degree+1;
    const unsigned int dofs_per_cell = Utilities::pow(degree+1,dim);
    const unsigned int n_q_points = Utilities::pow(n_q_points_1d,dim);

    kernels::EvaluatorTensorProduct<dim,(kernel>-1?degree+1:0),(kernel>-1?n_q_points_1d:0),
      VectorizedArray<Number>,kernel>
      tp(kernel==4 ? shape_values_eo : shape_values,
         kernel==4 ? shape_gradients_eo : shape_gradients,
         runtime_degree+1, runtime_degree+1);

    for (unsigned int cell=0; cell<vector_offsets.size(); ++cell)
      {
        if (dim == 3)
          {
            tp.template values<0,true,false>(input_array.begin(), scratch);
            tp.template values<1,true,false>(scratch, scratch+dofs_per_cell);
            tp.template values<2,true,false>(scratch+dofs_per_cell, scratch);
            tp.template gradients<2,true,false>(scratch,scratch+3*dofs_per_cell);
          }
        else
          {
            tp.template values<0,true,false>(input_array.begin(), scratch+dofs_per_cell);
            tp.template values<1,true,false>(scratch+dofs_per_cell, scratch);
          }
        tp.template gradients<1,true,false>(scratch, scratch+2*dofs_per_cell);
        tp.template gradients<0,true,false>(scratch, scratch+dofs_per_cell);
        VectorizedArray<Number>* phi_grads = scratch + dofs_per_cell;
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            if (dim==2)
              {
                VectorizedArray<Number> t0 = phi_grads[q]*jacobian_data[0] + phi_grads[q+n_q_points]*jacobian_data[2];
                VectorizedArray<Number> t1 = phi_grads[q]*jacobian_data[2] + phi_grads[q+n_q_points]*jacobian_data[1];
                phi_grads[q] = t0 * quadrature_weights[q];
                phi_grads[q+n_q_points] = t1 * quadrature_weights[q];
              }
            else if (dim==3)
              {
                VectorizedArray<Number> t0 = phi_grads[q]*jacobian_data[0] + phi_grads[q+n_q_points]*jacobian_data[3]+phi_grads[q+2*n_q_points]*jacobian_data[4];
                VectorizedArray<Number> t1 = phi_grads[q]*jacobian_data[3] + phi_grads[q+n_q_points]*jacobian_data[1]+phi_grads[q+2*n_q_points]*jacobian_data[5];
                VectorizedArray<Number> t2 = phi_grads[q]*jacobian_data[4] + phi_grads[q+n_q_points]*jacobian_data[5]+phi_grads[q+2*n_q_points]*jacobian_data[2];
                phi_grads[q] = t0 * quadrature_weights[q];
                phi_grads[q+n_q_points] = t1 * quadrature_weights[q];
                phi_grads[q+2*n_q_points] = t2 * quadrature_weights[q];
              }
          }

        tp.template gradients<0,false,false>(scratch+dofs_per_cell, scratch);
        tp.template gradients<1,false,true>(scratch+2*dofs_per_cell, scratch);
        if (dim == 3)
          tp.template gradients<2,false,true>(scratch+3*dofs_per_cell, scratch);
        if (dim == 3)
          {
            tp.template values<2,false,false>(scratch, scratch+dofs_per_cell);
            tp.template values<1,false,false>(scratch+dofs_per_cell, scratch);
            tp.template values<0,false,false>(scratch, output_array.begin());
          }
        else
          {
            tp.template values<1,false,false>(scratch, scratch+dofs_per_cell);
            tp.template values<0,false,false>(scratch+dofs_per_cell, output_array.begin());
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
    shape_values.resize((degree+1)*n_q_points_1d);
    shape_gradients.resize((degree+1)*n_q_points_1d);
    shape_values_eo.resize((degree+1)*stride);
    shape_gradients_eo.resize((degree+1)*stride);

    LagrangePolynomialBasis basis_gll(get_gauss_lobatto_points(degree+1));
    std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gll.value(i, gauss_points[q]);
          const double p2 = basis_gll.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values[i*n_q_points_1d+q] = p1;
          shape_values[(i+1)*n_q_points_1d-1-q] = p2;
          shape_values[(degree-i)*n_q_points_1d+q] = p2;
          shape_values[(degree+1-i)*n_q_points_1d-1-q] = p1;
          shape_values_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_values_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        {
          shape_values[degree/2*n_q_points_1d+q] =
            basis_gll.value(degree/2, gauss_points[q]);
          shape_values[degree/2*n_q_points_1d+n_q_points_1d-1-q] =
            shape_values[degree/2*n_q_points_1d+q];
          shape_values_eo[degree/2*stride+q] =
            shape_values[degree/2*n_q_points_1d+q];
        }

    LagrangePolynomialBasis basis_gauss(get_gauss_points(degree+1));
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gauss.derivative(i, gauss_points[q]);
          const double p2 = basis_gauss.derivative(i, gauss_points[n_q_points_1d-1-q]);
          shape_gradients[i*n_q_points_1d+q] = p1;
          shape_gradients[(i+1)*n_q_points_1d-1-q] = p2;
          shape_gradients[(degree-i)*n_q_points_1d+q] = -p2;
          shape_gradients[(degree+1-i)*n_q_points_1d-1-q] = -p1;
          shape_gradients_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_gradients_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        {
          shape_gradients[degree/2*n_q_points_1d+q] =
            basis_gauss.derivative(degree/2, gauss_points[q]);
          shape_gradients[degree/2*n_q_points_1d+n_q_points_1d-1-q] =
            -shape_gradients[degree/2*n_q_points_1d+q];
          shape_gradients_eo[degree/2*stride+q] =
            shape_gradients[degree/2*n_q_points_1d+q];
        }

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

    scratch_array.resize((dim+1)*quadrature_weights.size());
  }

  unsigned int runtime_degree;

  AlignedVector<VectorizedArray<Number> > shape_values;
  AlignedVector<VectorizedArray<Number> > shape_gradients;
  AlignedVector<VectorizedArray<Number> > shape_values_eo;
  AlignedVector<VectorizedArray<Number> > shape_gradients_eo;

  AlignedVector<Number> quadrature_weights;

  AlignedVector<unsigned int> vector_offsets;
  AlignedVector<VectorizedArray<Number> > input_array;
  AlignedVector<VectorizedArray<Number> > output_array;
  AlignedVector<VectorizedArray<Number> > scratch_array;

  AlignedVector<VectorizedArray<Number> > convection;

  AlignedVector<unsigned int> data_offsets;
  AlignedVector<VectorizedArray<Number> > jacobian_data;
};



template <int dim, int degree, typename Number, int kernel>
double run_kernel(const unsigned int n_tests)
{
  const unsigned int n_cell_batches = std::max(vector_size_guess / Utilities::pow(degree+1,dim),
                                               1UL);
  double best_avg = std::numeric_limits<double>::max();
#ifdef _OPENMP
  const unsigned int nthreads = omp_get_max_threads();
#else
  const unsigned int nthreads = 1;
#endif

  for (unsigned int i=0; i<3; ++i)
    {
      std::vector<double> min_time(nthreads), max_time(nthreads), avg_time(nthreads), std_dev(nthreads);
#pragma omp parallel shared(min_time, max_time, avg_time, std_dev)
      {
        EvaluationCellLaplacian<dim,degree,Number> evaluator;
        const unsigned int n_cell_batches = std::max(vector_size_guess / Utilities::pow(degree+1,dim),
                                                     1UL);
        evaluator.initialize(n_cell_batches, cartesian);

#ifdef LIKWID_PERFMON
        if (kernel == 4)
          LIKWID_MARKER_START(("cell_laplacian_deg_" + std::to_string(degree)).c_str());
#endif
        double tmin = 1e10, tmax = 0, tavg = 0, variance = 0;

#pragma omp barrier

#pragma omp for schedule(static)
        for (unsigned int thr=0; thr<nthreads; ++thr)
          for (unsigned int t=0; t<(50/degree); ++t)
            evaluator.template matrix_vector_product<kernel>();

#pragma omp barrier

#pragma omp for schedule(static)
        for (unsigned int thr=0; thr<nthreads; ++thr)
          {
            for (unsigned int t=0; t<n_tests; ++t)
              {
                auto t1 = std::chrono::system_clock::now();
                evaluator.template matrix_vector_product<kernel>();
                double time = std::chrono::duration<double>(std::chrono::system_clock::now()-t1).count();
                tmin = std::min(tmin, time);
                tmax = std::max(tmax, time);
                variance = (t > 0 ? (double)(t-1)/(t)*variance : 0) + (time - tavg) * (time - tavg) / (t+1);
                tavg = tavg + (time-tavg)/(t+1);
              }
            min_time[thr] = tmin;
            max_time[thr] = tmax;
            avg_time[thr] = tavg;
            std_dev[thr] = std::sqrt(variance);
          }

#ifdef LIKWID_PERFMON
        if (kernel == 4)
          LIKWID_MARKER_STOP(("cell_laplacian_deg_" + std::to_string(degree)).c_str());
#endif
      }
      double tmin = min_time[0], tmax = max_time[0], tavg = 0, stddev = 0;
      for (unsigned int i=0; i<nthreads; ++i)
        {
          std::cout << "p" << degree << " ker" << kernel << " statistics t" << std::setw(3) << std::right << i << "  ";
          std::cout << std::setw(12) << avg_time[i]
                    << "  dev " << std::setw(12) << std_dev[i]
                    << "  min " << std::setw(12) << min_time[i]
                    << "  max " << std::setw(12) << max_time[i]
                    << std::endl;
          tmin = std::min(min_time[i], tmin);
          tmax = std::max(max_time[i], tmax);
          tavg += avg_time[i] / nthreads;
          stddev = std::max(stddev, std_dev[i]);
        }
      const std::size_t n_dofs = (std::size_t)n_cell_batches * Utilities::pow(degree+1,dim) * VectorizedArray<Number>::n_array_elements * omp_get_max_threads();
      std::cout << "p" << degree << " ker" << kernel << " statistics tall  ";
      const std::size_t ops_interpolate = kernel < 4 ?
                                                   (degree+1 +
                                                    2*degree*(degree+1))
                                                   :(/*add*/2*((degree+1)/2)*2 +
                                                     /*mult*/degree+1 +
                                                     /*fma*/2*((degree-1)*(degree+1)/2));
      const std::size_t ops_approx = (std::size_t)n_cell_batches * VectorizedArray<Number>::n_array_elements * omp_get_max_threads()
        * (4 * dim * ops_interpolate * Utilities::pow(degree+1,dim-1)
           + dim * 2 * dim * Utilities::pow(degree+1,dim));
      std::cout << std::setw(12) << tavg
                << "  dev " << std::setw(12) << stddev
                << "  min " << std::setw(12) << tmin
                << "  max " << std::setw(12) << tmax
                << "  DoFs/s " << n_dofs / tavg
                << "  GFlops/s " << 1e-9*ops_approx / tavg
                << std::endl << std::endl;
      best_avg = std::min(tavg, best_avg);
    }
  return best_avg;
}


template <int dim, int degree, typename Number>
void run_program(const unsigned int n_tests)
{
  EvaluationCellLaplacian<dim,degree,Number> evaluator;
  evaluator.initialize(3, cartesian);
  evaluator.do_verification();

  run_kernel<dim,degree,Number,-1>(n_tests);
  run_kernel<dim,degree,Number,0>(n_tests);
  run_kernel<dim,degree,Number,1>(n_tests);
  run_kernel<dim,degree,Number,2>(n_tests);
  run_kernel<dim,degree,Number,3>(n_tests);
  run_kernel<dim,degree,Number,4>(n_tests);
}


template<int dim, int degree, int max_degree, typename Number>
class RunTime
{
public:
  static void run(const unsigned int degree_select,
                  const unsigned int n_tests)
  {
    if (degree==degree_select)
      run_program<dim,degree,Number>(n_tests);
    if (degree < max_degree)
      RunTime<dim,(degree<max_degree?degree+1:degree),max_degree,Number>::run(degree_select, n_tests);
  }
};

int main(int argc, char** argv)
{
#ifdef LIKWID_PERFMON
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
  }
#endif

  unsigned int n_tests           = 100;
  if (argc > 2)
    n_tests = std::atoi(argv[2]);
  unsigned int degree = 4;
  if (argc > 1)
    degree = std::atoi(argv[1]);

  //RunTime<2,min_degree,max_degree,value_type>::run();
  RunTime<3,min_degree,max_degree,value_type>::run(degree, n_tests);

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif

  return 0;
}
