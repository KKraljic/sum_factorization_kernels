// this file is inspired by the file
// deal.II/include/deal.II/matrix_free/tensor_product_kernels.h see
// www.dealii.org for information about licenses.

#include "utilities.h"
#include "aligned_vector.h"
#include "vectorization.h"
#include "gauss_formula.h"
#include "lagrange_polynomials.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>

template <int dim, int fe_degree, int n_q_points_1d, typename Number, int direction, bool dof_to_quad, bool add>
inline
void
evaluate_tensor_product (const Number *__restrict shape_data,
                         const Number *__restrict in,
                         Number       *__restrict out)
{
  const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

  const int stride    = Utilities::pow(nn,direction);
  const int jump_in   = stride * (mm-1);
  const int jump_out  = stride * (nn-1);

  const int n_blocks1 = stride;
  const int n_blocks2 = Utilities::pow(mm,(direction>=dim)?0:(dim-direction-1));

  const int nn_8 = (nn / 8) * 8;
  for (int i2=0; i2<n_blocks2; ++i2)
    {
      for (int i1=0; i1<n_blocks1; ++i1)
        {
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
                      val0 += shape_data[ind*n_q_points_1d+col+0] * f;
                      val1 += shape_data[ind*n_q_points_1d+col+1] * f;
                      val2 += shape_data[ind*n_q_points_1d+col+2] * f;
                      val3 += shape_data[ind*n_q_points_1d+col+3] * f;
                      val4 += shape_data[ind*n_q_points_1d+col+4] * f;
                      val5 += shape_data[ind*n_q_points_1d+col+5] * f;
                      val6 += shape_data[ind*n_q_points_1d+col+6] * f;
                      val7 += shape_data[ind*n_q_points_1d+col+7] * f;
                    }
                }
              else
                {
                  const Number f = in[0];
                  val0 = shape_data[(col+0)*n_q_points_1d] * f;
                  val1 = shape_data[(col+1)*n_q_points_1d] * f;
                  val2 = shape_data[(col+2)*n_q_points_1d] * f;
                  val3 = shape_data[(col+3)*n_q_points_1d] * f;
                  val4 = shape_data[(col+4)*n_q_points_1d] * f;
                  val5 = shape_data[(col+5)*n_q_points_1d] * f;
                  val6 = shape_data[(col+6)*n_q_points_1d] * f;
                  val7 = shape_data[(col+7)*n_q_points_1d] * f;
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f = in[stride*ind];
                      val0 += shape_data[(col+0)*n_q_points_1d+ind] * f;
                      val1 += shape_data[(col+1)*n_q_points_1d+ind] * f;
                      val2 += shape_data[(col+2)*n_q_points_1d+ind] * f;
                      val3 += shape_data[(col+3)*n_q_points_1d+ind] * f;
                      val4 += shape_data[(col+4)*n_q_points_1d+ind] * f;
                      val5 += shape_data[(col+5)*n_q_points_1d+ind] * f;
                      val6 += shape_data[(col+6)*n_q_points_1d+ind] * f;
                      val7 += shape_data[(col+7)*n_q_points_1d+ind] * f;
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
          const unsigned int remainder = nn - nn_8;
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
                      val0 += shape_data[ind*n_q_points_1d+nn_8+0] * f;
                      if (remainder > 1) val1 += shape_data[ind*n_q_points_1d+nn_8+1] * f;
                      if (remainder > 2) val2 += shape_data[ind*n_q_points_1d+nn_8+2] * f;
                      if (remainder > 3) val3 += shape_data[ind*n_q_points_1d+nn_8+3] * f;
                      if (remainder > 4) val4 += shape_data[ind*n_q_points_1d+nn_8+4] * f;
                      if (remainder > 5) val5 += shape_data[ind*n_q_points_1d+nn_8+5] * f;
                      if (remainder > 6) val6 += shape_data[ind*n_q_points_1d+nn_8+6] * f;
                    }
                }
              else
                {
                  const Number f = in[0];
                  val0 = shape_data[(nn_8+0)*n_q_points_1d] * f;
                  if (remainder > 1) val1 = shape_data[(nn_8+1)*n_q_points_1d] * f;
                  if (remainder > 2) val2 = shape_data[(nn_8+2)*n_q_points_1d] * f;
                  if (remainder > 3) val3 = shape_data[(nn_8+3)*n_q_points_1d] * f;
                  if (remainder > 4) val4 = shape_data[(nn_8+4)*n_q_points_1d] * f;
                  if (remainder > 5) val5 = shape_data[(nn_8+5)*n_q_points_1d] * f;
                  if (remainder > 6) val6 = shape_data[(nn_8+6)*n_q_points_1d] * f;
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f = in[stride*ind];
                      val0 += shape_data[(nn_8+0)*n_q_points_1d+ind] * f;
                      if (remainder > 1) val1 += shape_data[(nn_8+1)*n_q_points_1d+ind] * f;
                      if (remainder > 2) val2 += shape_data[(nn_8+2)*n_q_points_1d+ind] * f;
                      if (remainder > 3) val3 += shape_data[(nn_8+3)*n_q_points_1d+ind] * f;
                      if (remainder > 4) val4 += shape_data[(nn_8+4)*n_q_points_1d+ind] * f;
                      if (remainder > 5) val5 += shape_data[(nn_8+5)*n_q_points_1d+ind] * f;
                      if (remainder > 6) val6 += shape_data[(nn_8+6)*n_q_points_1d+ind] * f;
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
          ++in;
          ++out;
        }
      in += jump_in;
      out+= jump_out;
    }
}



template <int dim, int fe_degree, int n_q_points_1d, typename Number, int direction, bool dof_to_quad, bool add>
inline
void
evaluate_tensor_product_plain (const Number *__restrict shape_data,
                               const Number *__restrict in,
                               Number       *__restrict out)
{
  const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

  const int stride    = Utilities::pow(nn,direction);
  const int jump_in   = stride * (mm-1);
  const int jump_out  = stride * (nn-1);

  const int n_blocks1 = stride;
  const int n_blocks2 = Utilities::pow(mm,(direction>=dim)?0:(dim-direction-1));

  for (int i2=0; i2<n_blocks2; ++i2)
    {
      for (int i1=0; i1<n_blocks1; ++i1)
        {
          for (int col=0; col<nn; ++col)
            {
              Number val0;
              if (dof_to_quad)
                {
                  val0 = shape_data[col+0] * in[0];
                  for (int ind=1; ind<mm; ++ind)
                    val0 += shape_data[ind*n_q_points_1d+col+0] * in[stride*ind];
                }
              else
                {
                  val0 = shape_data[(col+0)*n_q_points_1d] * in[0];
                  for (int ind=1; ind<mm; ++ind)
                    val0 += shape_data[(col+0)*n_q_points_1d+ind] * in[stride*ind];
                }
              if (add == false)
                out[stride*(col+0)]  = val0;
              else
                out[stride*(col+0)] += val0;
            }
          ++in;
          ++out;
        }
      in += jump_in;
      out+= jump_out;
    }
}



template <int dim, int fe_degree, int n_q_points_1d, typename Number, int direction, bool dof_to_quad, bool add>
inline
void
evaluate_tensor_product_2x5 (const Number *__restrict shape_data,
                             const Number *__restrict in,
                             Number       *__restrict out)
{
  if (fe_degree % 2 == 0)
    std::abort();
  const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

  const int stride    = Utilities::pow(nn,direction);
  const int jump_in   = stride * (mm-1);
  const int jump_out  = stride * (nn-1);

  const int n_blocks1 = stride;
  const int n_blocks2 = Utilities::pow(mm,(direction>=dim)?0:(dim-direction-1));

  const int n_blocks1a = direction > 0 ? n_blocks1/2 : n_blocks1;
  const int n_blocks2a = direction > 0 ? n_blocks2   : n_blocks2/2;

  const unsigned nn_5 = (nn/5)*5;
  for (int i2=0; i2<n_blocks2a; ++i2)
    {
      for (int i1=0; i1<n_blocks1a; ++i1)
        {
          const Number *in2 = in +(direction>0?1:mm);
          Number      *out2 = out+(direction>0?1:nn);
          for (int col=0; col<nn_5; col+=5)
            {
              Number val0, val1, val2, val3, val4, val5, val6, val7, val8, val9;
              if (dof_to_quad)
                {
                  const Number f1 = in[0];
                  const Number f2 = in2[0];
                  Number t = shape_data[col+0];
                  val0 = t * f1;
                  val1 = t * f2;
                  t = shape_data[col+1];
                  val2 = t * f1;
                  val3 = t * f2;
                  t = shape_data[col+2];
                  val4 = t * f1;
                  val5 = t * f2;
                  t = shape_data[col+3];
                  val6 = t * f1;
                  val7 = t * f2;
                  t = shape_data[col+4];
                  val8 = t * f1;
                  val9 = t * f2;
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f1 = in[stride*ind];
                      const Number f2 = in2[stride*ind];
                      Number t = shape_data[ind*n_q_points_1d+col+0];
                      val0 += t * f1;
                      val1 += t * f2;
                      t = shape_data[ind*n_q_points_1d+col+1];
                      val2 += t * f1;
                      val3 += t * f2;
                      t = shape_data[ind*n_q_points_1d+col+2];
                      val4 += t * f1;
                      val5 += t * f2;
                      t = shape_data[ind*n_q_points_1d+col+3];
                      val6 += t * f1;
                      val7 += t * f2;
                      t = shape_data[ind*n_q_points_1d+col+4];
                      val8 += t * f1;
                      val9 += t * f2;
                    }
                }
              else
                {
                  const Number f1 = in[0];
                  const Number f2 = in2[0];
                  Number t = shape_data[(col+0)*n_q_points_1d];
                  val0 = t * f1;
                  val1 = t * f2;
                  t = shape_data[(col+1)*n_q_points_1d];
                  val2 = t * f1;
                  val3 = t * f2;
                  t = shape_data[(col+2)*n_q_points_1d];
                  val4 = t * f1;
                  val5 = t * f2;
                  t = shape_data[(col+3)*n_q_points_1d];
                  val6 = t * f1;
                  val7 = t * f2;
                  t = shape_data[(col+4)*n_q_points_1d];
                  val8 = t * f1;
                  val9 = t * f2;
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f1 = in[stride*ind];
                      const Number f2 = in2[stride*ind];
                      Number t = shape_data[(col+0)*n_q_points_1d+ind];
                      val0 += t * f1;
                      val1 += t * f2;
                      t = shape_data[(col+1)*n_q_points_1d+ind];
                      val2 += t * f1;
                      val3 += t * f2;
                      t = shape_data[(col+2)*n_q_points_1d+ind];
                      val4 += t * f1;
                      val5 += t * f2;
                      t = shape_data[(col+3)*n_q_points_1d+ind];
                      val6 += t * f1;
                      val7 += t * f2;
                      t = shape_data[(col+4)*n_q_points_1d+ind];
                      val8 += t * f1;
                      val9 += t * f2;
                    }
                }
              if (add == false)
                {
                  out[stride*(col+0)]   = val0;
                  out[stride*(col+1)]   = val2;
                  out[stride*(col+2)]   = val4;
                  out[stride*(col+3)]   = val6;
                  out[stride*(col+4)]   = val8;
                  out2[stride*(col+0)]  = val1;
                  out2[stride*(col+1)]  = val3;
                  out2[stride*(col+2)]  = val5;
                  out2[stride*(col+3)]  = val7;
                  out2[stride*(col+4)]  = val9;
                }
              else
                {
                  out[stride*(col+0)]  += val0;
                  out[stride*(col+1)]  += val2;
                  out[stride*(col+2)]  += val4;
                  out[stride*(col+3)]  += val6;
                  out[stride*(col+4)]  += val8;
                  out2[stride*(col+0)] += val1;
                  out2[stride*(col+1)] += val3;
                  out2[stride*(col+2)] += val5;
                  out2[stride*(col+3)] += val7;
                  out2[stride*(col+4)] += val9;
                }
            }
          Number val[8];
          const unsigned int remainder = nn - nn_5;
          if (remainder > 0)
            {
              if (dof_to_quad)
                {
                  const Number f1 = in[0];
                  const Number f2 = in2[0];
                  for (unsigned int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[nn_5+i];
                      val[2*i]   = t*f1;
                      val[2*i+1] = t*f2;
                    }
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f1 = in[stride*ind];
                      const Number f2 = in2[stride*ind];
                      for (unsigned int i=0; i<remainder; ++i)
                        {
                          const Number t = shape_data[ind*n_q_points_1d+nn_5+i];
                          val[2*i]   += t*f1;
                          val[2*i+1] += t*f2;
                        }
                    }
                }
              else
                {
                  const Number f1 = in[0];
                  const Number f2 = in2[0];
                  for (unsigned int i=0; i<remainder; ++i)
                    {
                      const Number t = shape_data[(nn_5+i)*n_q_points_1d];
                      val[2*i]   = t*f1;
                      val[2*i+1] = t*f2;
                    }
                  for (int ind=1; ind<mm; ++ind)
                    {
                      const Number f1 = in[stride*ind];
                      const Number f2 = in2[stride*ind];
                      for (unsigned int i=0; i<remainder; ++i)
                        {
                          const Number t = shape_data[(nn_5+i)*n_q_points_1d+ind];
                          val[2*i]   += t*f1;
                          val[2*i+1] += t*f2;
                        }
                    }
                }
              if (add == false)
                for (unsigned int i=0; i<remainder; ++i)
                  {
                    out[stride*(nn_5+i)]  = val[2*i];
                    out2[stride*(nn_5+i)] = val[2*i+1];
                  }
              else
                for (unsigned int i=0; i<remainder; ++i)
                  {
                    out[stride*(nn_5+i)]  += val[2*i];
                    out2[stride*(nn_5+i)] += val[2*i+1];
                  }
            }
          in  += (direction>0 ? 2 : 1);
          out += (direction>0 ? 2 : 1);
        }
      in += (direction>0 ? jump_in : (jump_in+mm));
      out+= (direction>0 ? jump_out : (jump_out+nn));
    }
}


template <int dim, int fe_degree, int n_q_points_1d, typename Number, int direction, bool dof_to_quad, bool add>
inline
void
evaluate_tensor_product_eo (const Number *__restrict shapes,
                            const Number *__restrict in,
                            Number       *__restrict out)
{
  const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = Utilities::pow(nn,direction);

  const int offset = (n_q_points_1d+1)/2;

  for (int i2=0; i2<n_blocks2; ++i2)
    {
      for (int i1=0; i1<n_blocks1; ++i1)
        {
          Number xp[mid>0?mid:1], xm[mid>0?mid:1];
          for (int i=0; i<mid; ++i)
            {
              xp[i] = in[stride*i] + in[stride*(mm-1-i)];
              xm[i] = in[stride*i] - in[stride*(mm-1-i)];
            }
          for (int col=0; col<n_cols; ++col)
            {
              Number r0, r1;
              if (mid > 0)
                {
                  if (dof_to_quad == true)
                    {
                      r0 = shapes[col]                    * xp[0];
                      r1 = shapes[fe_degree*offset + col] * xm[0];
                      for (int ind=1; ind<mid; ++ind)
                        {
                          r0 += shapes[ind*offset+col]             * xp[ind];
                          r1 += shapes[(fe_degree-ind)*offset+col] * xm[ind];
                        }
                    }
                  else
                    {
                      r0 = shapes[col*offset]             * xp[0];
                      r1 = shapes[(fe_degree-col)*offset] * xm[0];
                      for (int ind=1; ind<mid; ++ind)
                        {
                          r0 += shapes[col*offset+ind]             * xp[ind];
                          r1 += shapes[(fe_degree-col)*offset+ind] * xm[ind];
                        }
                    }
                }
              else
                r0 = r1 = Number();
              if (mm % 2 == 1 && dof_to_quad == true)
                r0 += shapes[mid*offset+col] * in[stride*mid];
              else if (mm % 2 == 1 && nn % 2 == 0)
                r0 += shapes[col*offset+mid] * in[stride*mid];

              if (add == false)
                {
                  out[stride*col]         = r0 + r1;
                  out[stride*(nn-1-col)]  = r0 - r1;
                }
              else
                {
                  out[stride*col]        += r0 + r1;
                  out[stride*(nn-1-col)] += r0 - r1;
                }
            }
          if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
            {
              if (add==false)
                out[stride*n_cols]  = in[stride*mid];
              else
                out[stride*n_cols] += in[stride*mid];
            }
          else if (dof_to_quad == true && nn%2==1)
            {
              Number r0;
              if (mid > 0)
                {
                  r0  = shapes[n_cols] * xp[0];
                  for (int ind=1; ind<mid; ++ind)
                    r0 += shapes[ind*offset+n_cols] * xp[ind];
                }
              else
                r0 = Number();
              if (mm % 2 == 1)
                r0 += shapes[mid*offset+n_cols] * in[stride*mid];

              if (add == false)
                out[stride*n_cols]  = r0;
              else
                out[stride*n_cols] += r0;
            }
          else if (dof_to_quad == false && nn%2 == 1)
            {
              Number r0;
              if (mid > 0)
                {
                  r0 = shapes[n_cols*offset] * xp[0];
                  for (int ind=1; ind<mid; ++ind)
                    r0 += shapes[n_cols*offset+ind] * xp[ind];
                }
              else
                r0 = Number();

              if (mm % 2 == 1)
                r0 += in[stride*mid];

              if (add == false)
                out[stride*n_cols]  = r0;
              else
                out[stride*n_cols] += r0;
            }

          // increment: in regular case, just go to the next point in
          // x-direction. If we are at the end of one chunk in x-dir, need to
          // jump over to the next layer in z-direction
          switch (direction)
            {
            case 0:
              in += mm;
              out += nn;
              break;
            case 1:
            case 2:
              ++in;
              ++out;
              break;
            default:
              std::abort();
            }
        }
      if (direction == 1)
        {
          in += nn*(mm-1);
          out += nn*(nn-1);
        }
    }
}


template <int dim, int fe_degree, int n_q_points_1d, typename Number, int direction, bool dof_to_quad, bool add>
inline
void
evaluate_tensor_product_eo_2 (const Number *__restrict shapes,
                              const Number *__restrict in,
                              Number       *__restrict out)
{
  const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = Utilities::pow(nn,direction);

  const int offset = (n_q_points_1d+1)/2;
  const int n_cols_2 = (n_cols/2)*2;

  for (int i2=0; i2<n_blocks2; ++i2)
    {
      for (int i1=0; i1<n_blocks1; ++i1)
        {
          Number xp[mid>0?mid:1], xm[mid>0?mid:1];
          for (int i=0; i<mid; ++i)
            {
              xp[i] = in[stride*i] + in[stride*(mm-1-i)];
              xm[i] = in[stride*i] - in[stride*(mm-1-i)];
            }
          for (int col=0; col<n_cols_2; col+=2)
            {
              Number r0, r1, r2, r3;
              if (dof_to_quad == true)
                {
                  r0 = shapes[col]                      * xp[0];
                  r1 = shapes[fe_degree*offset + col]   * xm[0];
                  r2 = shapes[col+1]                    * xp[0];
                  r3 = shapes[fe_degree*offset + col+1] * xm[0];
                  for (int ind=1; ind<mid; ++ind)
                    {
                      r0 += shapes[ind*offset+col]               * xp[ind];
                      r1 += shapes[(fe_degree-ind)*offset+col]   * xm[ind];
                      r2 += shapes[ind*offset+col+1]             * xp[ind];
                      r3 += shapes[(fe_degree-ind)*offset+col+1] * xm[ind];
                    }
                }
              else
                {
                  r0 = shapes[col*offset]               * xp[0];
                  r1 = shapes[(fe_degree-col)*offset]   * xm[0];
                  r2 = shapes[(col+1)*offset]           * xp[0];
                  r3 = shapes[(fe_degree-col-1)*offset] * xm[0];
                  for (int ind=1; ind<mid; ++ind)
                    {
                      r0 += shapes[col*offset+ind]               * xp[ind];
                      r1 += shapes[(fe_degree-col)*offset+ind]   * xm[ind];
                      r2 += shapes[(col+1)*offset+ind]           * xp[ind];
                      r3 += shapes[(fe_degree-col-1)*offset+ind] * xm[ind];
                    }
                }
              if (mm % 2 == 1 && dof_to_quad == true)
                {
                  r0 += shapes[mid*offset+col]   * in[stride*mid];
                  r2 += shapes[mid*offset+col+1] * in[stride*mid];
                }
              else if (mm % 2 == 1 && nn % 2 == 0)
                {
                  r0 += shapes[col*offset+mid]     * in[stride*mid];
                  r2 += shapes[(col+1)*offset+mid] * in[stride*mid];
                }

              if (add == false)
                {
                  out[stride*col]         = r0 + r1;
                  out[stride*(nn-1-col)]  = r0 - r1;
                  out[stride*col+1]       = r2 + r3;
                  out[stride*(nn-2-col)]  = r2 - r3;
                }
              else
                {
                  out[stride*col]        += r0 + r1;
                  out[stride*(nn-1-col)] += r0 - r1;
                  out[stride*col+1]      += r2 + r3;
                  out[stride*(nn-2-col)] += r2 - r3;
                }
            }
          for (int col=n_cols_2; col<n_cols; ++col)
            {
              Number r0, r1;
              if (dof_to_quad == true)
                {
                  r0 = shapes[col]                      * xp[0];
                  r1 = shapes[fe_degree*offset + col]   * xm[0];
                  for (int ind=1; ind<mid; ++ind)
                    {
                      r0 += shapes[ind*offset+col]               * xp[ind];
                      r1 += shapes[(fe_degree-ind)*offset+col]   * xm[ind];
                    }
                }
              else
                {
                  r0 = shapes[col*offset]               * xp[0];
                  r1 = shapes[(fe_degree-col)*offset]   * xm[0];
                  for (int ind=1; ind<mid; ++ind)
                    {
                      r0 += shapes[col*offset+ind]               * xp[ind];
                      r1 += shapes[(fe_degree-col)*offset+ind]   * xm[ind];
                    }
                }
              if (mm % 2 == 1 && dof_to_quad == true)
                r0 += shapes[mid*offset+col]   * in[stride*mid];
              else if (mm % 2 == 1 && nn % 2 == 0)
                r0 += shapes[col*offset+mid]     * in[stride*mid];

              if (add == false)
                {
                  out[stride*col]         = r0 + r1;
                  out[stride*(nn-1-col)]  = r0 - r1;
                }
              else
                {
                  out[stride*col]        += r0 + r1;
                  out[stride*(nn-1-col)] += r0 - r1;
                }
            }
          if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
            {
              if (add==false)
                out[stride*n_cols]  = in[stride*mid];
              else
                out[stride*n_cols] += in[stride*mid];
            }
          else if (dof_to_quad == true && nn%2==1)
            {
              Number r0;
              if (mid > 0)
                {
                  r0  = shapes[n_cols] * xp[0];
                  for (int ind=1; ind<mid; ++ind)
                    r0 += shapes[ind*offset+n_cols] * xp[ind];
                }
              else
                r0 = Number();
              if (mm % 2 == 1)
                r0 += shapes[mid*offset+n_cols] * in[stride*mid];

              if (add == false)
                out[stride*n_cols]  = r0;
              else
                out[stride*n_cols] += r0;
            }
          else if (dof_to_quad == false && nn%2 == 1)
            {
              Number r0;
              if (mid > 0)
                {
                  r0 = shapes[n_cols*offset] * xp[0];
                  for (int ind=1; ind<mid; ++ind)
                    r0 += shapes[n_cols*offset+ind] * xp[ind];
                }
              else
                r0 = Number();

              if (mm % 2 == 1)
                r0 += in[stride*mid];

              if (add == false)
                out[stride*n_cols]  = r0;
              else
                out[stride*n_cols] += r0;
            }

          // increment: in regular case, just go to the next point in
          // x-direction. If we are at the end of one chunk in x-dir, need to
          // jump over to the next layer in z-direction
          switch (direction)
            {
            case 0:
              in += mm;
              out += nn;
              break;
            case 1:
            case 2:
              ++in;
              ++out;
              break;
            default:
              std::abort();
            }
        }
      if (direction == 1)
        {
          in += nn*(mm-1);
          out += nn*(nn-1);
        }
    }
}

template <int degree, typename Number>
void do_test()
{
  const unsigned int n_q_points_1d = degree+1;
  AlignedVector<VectorizedArray<Number> > shape_values((degree+1)*n_q_points_1d);
  AlignedVector<VectorizedArray<Number> > shape_gradients((degree+1)*n_q_points_1d);

  LagrangePolynomialBasis basis_gll(get_gauss_lobatto_points(degree+1));
  LagrangePolynomialBasis basis_gauss(get_gauss_points(degree+1));
  std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
  for (unsigned int i=0; i<degree+1; ++i)
    for (unsigned int q=0; q<n_q_points_1d; ++q)
      {
        shape_values[i*n_q_points_1d+q] = basis_gll.value(i,gauss_points[q]);
        shape_gradients[i*n_q_points_1d+q] = basis_gauss.derivative(i,gauss_points[q]);
      }

  const unsigned int stride = (n_q_points_1d+1)/2;
  AlignedVector<VectorizedArray<Number> > shape_values_eo((degree+1)*stride);
  for (unsigned int i=0; i<(degree+1)/2; ++i)
    for (unsigned int q=0; q<stride; ++q)
      {
        const VectorizedArray<Number> p1 = shape_values[i*n_q_points_1d+q];
        const VectorizedArray<Number> p2 = shape_values[(i+1)*n_q_points_1d-1-q];
        shape_values_eo[i*stride+q] = 0.5 * (p1 + p2);
        shape_values_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
      }
  if (degree%2 == 0)
    for (unsigned int q=0; q<stride; ++q)
      shape_values_eo[degree/2*stride+q] = shape_values[degree/2*n_q_points_1d+q];

  AlignedVector<VectorizedArray<Number> > in(Utilities::pow(degree+1,3)),
    out(Utilities::pow(degree+1,3)), tmp(Utilities::pow(degree+1,3));
  for (unsigned int i=0; i<in.size(); ++i)
    in[i] = (Number)rand()/RAND_MAX;
  AlignedVector<VectorizedArray<Number> > in1(in.size());
  for (unsigned int i=0; i<in.size(); ++i)
    in1[i] = in[i];
  AlignedVector<VectorizedArray<Number> > in2(in.size());
  for (unsigned int i=0; i<in.size(); ++i)
    in2[i] = in[i];
  AlignedVector<VectorizedArray<Number> > in3(in.size());
  for (unsigned int i=0; i<in.size(); ++i)
    in3[i] = in[i];
  AlignedVector<VectorizedArray<Number> > in4(in.size());
  for (unsigned int i=0; i<in.size(); ++i)
    in4[i] = in[i];
  AlignedVector<VectorizedArray<Number> > in5(in.size());
  for (unsigned int i=0; i<in.size(); ++i)
    in5[i] = in[i];

  const unsigned int n_tests = 500000000 / Utilities::pow(degree+1,4);

  AlignedVector<Number> quadrature_weights(out.size());
  std::vector<double> gauss_weight_1d = get_gauss_weights(n_q_points_1d);
  for (unsigned int q=0, z=0; z<n_q_points_1d; ++z)
    for (unsigned int y=0; y<n_q_points_1d; ++y)
      for (unsigned int x=0; x<n_q_points_1d; ++x, ++q)
        quadrature_weights[q] = (gauss_weight_1d[z] * gauss_weight_1d[y]) *
          gauss_weight_1d[x] * 0.000001 * degree * degree * degree;

  struct timeval wall_timer;
  gettimeofday(&wall_timer, NULL);
  double start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              0,true,false>(shape_values.begin(), in1.begin(), out.begin());
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              1,true,false>(shape_values.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              2,true,false>(shape_values.begin(), tmp.begin(), out.begin());
      for (unsigned int i=0; i<Utilities::pow(degree+1,3); ++i)
        out[i] = out[i] * quadrature_weights[i];
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              2,false,false>(shape_values.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              1,false,false>(shape_values.begin(), tmp.begin(), out.begin());
      evaluate_tensor_product<3,degree,n_q_points_1d,VectorizedArray<Number>,
                              0,false,true> (shape_values.begin(), out.begin(), in1.begin());
    }

  gettimeofday(&wall_timer, NULL);
  double compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  const double ops = (double)n_tests * Utilities::pow(degree+1,3) * (2*degree+1) * 6 * VectorizedArray<Number>::n_array_elements;
  std::cout << "Performance unroll 8   for n=" << degree+1 << ": "
            << ops/compute_time*1e-6 << " MFLOP/s" << std::endl;


  gettimeofday(&wall_timer, NULL);
  start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    0,true,false>(shape_values.begin(), in2.begin(), out.begin());
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    1,true,false>(shape_values.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    2,true,false>(shape_values.begin(), tmp.begin(), out.begin());
      for (unsigned int i=0; i<Utilities::pow(degree+1,3); ++i)
        out[i] = out[i] * quadrature_weights[i];
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    2,false,false>(shape_values.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    1,false,false>(shape_values.begin(), tmp.begin(), out.begin());
      evaluate_tensor_product_plain<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                    0,false,true> (shape_values.begin(), out.begin(), in2.begin());
    }

  gettimeofday(&wall_timer, NULL);
  compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  std::cout << "Performance plain      for n=" << degree+1 << ": "
            << (double)ops/compute_time*1e-6 << " MFLOP/s" << std::endl;

  Number error2 = 0, ref = 0, error3 = 0, error4 = 0, error5 = 0;
  for (unsigned int i=0; i<in.size(); ++i)
    {
      ref = std::max(std::abs(in1[i][0]), ref);
      error2 = std::max(std::abs(in1[i][0]-in2[i][0]), error2);
    }

  if (degree%2 == 1)
    {
      gettimeofday(&wall_timer, NULL);
      start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

      for (unsigned int t=0; t<n_tests; ++t)
        {
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      0,true,false>(shape_values.begin(), in3.begin(), out.begin());
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      1,true,false>(shape_values.begin(), out.begin(), tmp.begin());
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      2,true,false>(shape_values.begin(), tmp.begin(), out.begin());
          for (unsigned int i=0; i<Utilities::pow(degree+1,3); ++i)
            out[i] = out[i] * quadrature_weights[i];
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      2,false,false>(shape_values.begin(), out.begin(), tmp.begin());
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      1,false,false>(shape_values.begin(), tmp.begin(), out.begin());
          evaluate_tensor_product_2x5<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                      0,false,true> (shape_values.begin(), out.begin(), in3.begin());
        }

      gettimeofday(&wall_timer, NULL);
      compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

      std::cout << "Performance 2x5 loops  for n=" << degree+1 << ": "
                << (double)ops/compute_time*1e-6 << " MFLOP/s" << std::endl;

      for (unsigned int i=0; i<in.size(); ++i)
        error3 = std::max(std::abs(in1[i][0]-in3[i][0]), error3);
    }

  gettimeofday(&wall_timer, NULL);
  start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 0,true,false>(shape_values_eo.begin(), in4.begin(), out.begin());
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 1,true,false>(shape_values_eo.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 2,true,false>(shape_values_eo.begin(), tmp.begin(), out.begin());
      for (unsigned int i=0; i<Utilities::pow(degree+1,3); ++i)
        out[i] = out[i] * quadrature_weights[i];
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 2,false,false>(shape_values_eo.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 1,false,false>(shape_values_eo.begin(), tmp.begin(), out.begin());
      evaluate_tensor_product_eo<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                 0,false,true> (shape_values_eo.begin(), out.begin(), in4.begin());
    }

  gettimeofday(&wall_timer, NULL);
  compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  const double ops_eo = (double)n_tests
    * 6 * (/*add*/2*((degree+1)/2)*2 +
           /*mult*/degree+1 +
           /*fma*/2*((degree-1)*(degree+1)/2))*(degree+1)*(degree+1)
    * VectorizedArray<Number>::n_array_elements;

  std::cout << "Performance even-odd   for n=" << degree+1 << ": "
            << ops/compute_time*1e-6 << " MFLOP/s ("
            << ops_eo/compute_time*1e-6 << " MFLOP/s effective)" << std::endl;

  for (unsigned int i=0; i<in.size(); ++i)
    error4 = std::max(std::abs(in2[i][0]-in4[i][0]), error4);

  gettimeofday(&wall_timer, NULL);
  start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   0,true,false>(shape_values_eo.begin(), in5.begin(), out.begin());
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   1,true,false>(shape_values_eo.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   2,true,false>(shape_values_eo.begin(), tmp.begin(), out.begin());
      for (unsigned int i=0; i<Utilities::pow(degree+1,3); ++i)
        out[i] = out[i] * quadrature_weights[i];
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   2,false,false>(shape_values_eo.begin(), out.begin(), tmp.begin());
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   1,false,false>(shape_values_eo.begin(), tmp.begin(), out.begin());
      evaluate_tensor_product_eo_2<3,degree,n_q_points_1d,VectorizedArray<Number>,
                                   0,false,true> (shape_values_eo.begin(), out.begin(), in5.begin());
    }

  gettimeofday(&wall_timer, NULL);
  compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  std::cout << "Performance even-odd 2 for n=" << degree+1 << ": "
            << ops/compute_time*1e-6 << " MFLOP/s ("
            << ops_eo/compute_time*1e-6 << " MFLOP/s effective)" << std::endl;

  for (unsigned int i=0; i<in.size(); ++i)
    error5 = std::max(std::abs(in2[i][0]-in4[i][0]), error5);

  std::cout << "Array verification: " << error2/ref << " " << error3/ref << " "
             << error4/ref << " "<< error5/ref << " for ref " << ref << std::endl;
  std::cout << std::endl;
}


int main()
{
  do_test<1,double>();
  do_test<2,double>();
  do_test<3,double>();
  do_test<4,double>();
  do_test<5,double>();
  do_test<6,double>();
  do_test<7,double>();
  do_test<8,double>();
  do_test<9,double>();
  do_test<11,double>();
  do_test<15,double>();
  do_test<19,double>();
  do_test<23,double>();
  do_test<29,double>();
  do_test<39,double>();

  do_test<6,float>();
  do_test<7,float>();
  do_test<8,float>();
  do_test<9,float>();


  return 0;
}
