// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// This file is inspired by the file
// deal.II/include/deal.II/base/vectorization.h, see www.dealii.org for
// information about licenses. Here is the original deal.II license statement:
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

#ifndef vectorization_h
#define vectorization_h

#include <x86intrin.h>
#include <cstdlib>

#ifndef USE_VECTOR_ARITHMETICS
#define USE_VECTOR_ARITHMETICS 1
#endif

template <typename Number>
class VectorizedArray
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 1;

  /**
   * This function assigns a scalar to this class.
   */
  VectorizedArray &
  operator = (const Number scalar)
  {
    data = scalar;
    return *this;
  }

  /**
   * Access operator (only valid with component 0)
   */
  Number &
  operator [] (const unsigned int)
  {
    return data;
  }

  /**
   * Constant access operator (only valid with component 0)
   */
  const Number &
  operator [] (const unsigned int) const
  {
    return data;
  }

  /**
   * Addition
   */
  VectorizedArray &
  operator += (const VectorizedArray<Number> &vec)
  {
    data+=vec.data;
    return *this;
  }

  /**
   * Subtraction
   */
  VectorizedArray &
  operator -= (const VectorizedArray<Number> &vec)
  {
    data-=vec.data;
    return *this;
  }

  /**
   * Multiplication
   */
  VectorizedArray &
  operator *= (const VectorizedArray<Number> &vec)
  {
    data*=vec.data;
    return *this;
  }

  /**
   * Division
   */
  VectorizedArray &
  operator /= (const VectorizedArray<Number> &vec)
  {
    data/=vec.data;
    return *this;
  }

  void load (const Number *ptr)
  {
    data = *ptr;
  }

  void store (Number *ptr) const
  {
    *ptr = data;
  }

  void streaming_store (Number *ptr) const
  {
    *ptr = data;
  }

  void gather (const Number       *base_ptr,
               const unsigned int *offsets)
  {
    data = base_ptr[offsets[0]];
  }

  void scatter (const unsigned int *offsets,
                Number             *base_ptr) const
  {
    base_ptr[offsets[0]] = data;
  }

  VectorizedArray get_abs () const
  {
    VectorizedArray res;
    res.data = std::abs(data);
    return res;
  }

  Number data;
};



template <typename Number>
VectorizedArray<Number>
make_vectorized_array (const Number &u)
{
  VectorizedArray<Number> result;
  result = u;
  return result;
}

template <typename Number>
inline
void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const Number            *in,
                              const unsigned int      *offsets,
                              VectorizedArray<Number> *out)
{
  for (unsigned int i=0; i<n_entries; ++i)
    for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
      out[i][v] = in[offsets[v]+i];
}
template <typename Number>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<Number> *in,
                               const unsigned int            *offsets,
                               Number                        *out)
{
  if (add_into)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        out[offsets[v]+i] = in[i][v];
}


#if defined(__AVX512F__)

/**
 * Specialization of VectorizedArray class for double and AVX-512.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 8;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const double x)
  {
    data = _mm512_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm512_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm512_sub_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm512_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm512_div_pd(data,vec.data);
#endif
    return *this;
  }

  void load (const double *ptr)
  {
    data = _mm512_loadu_pd (ptr);
  }

  void store (double *ptr) const
  {
    _mm512_storeu_pd (ptr, data);
  }

  void streaming_store (double *ptr) const
  {
    _mm512_stream_pd(ptr,data);
  }

  void gather (const double       *base_ptr,
               const unsigned int *offsets)
  {
    const __m256 index_val = _mm256_loadu_ps((const float *)offsets);
    const __m256i index = *((__m256i *)(&index_val));
    data = _mm512_i32gather_pd(index, base_ptr, 8);
  }

  void scatter (const unsigned int *offsets,
                double             *base_ptr) const
  {
    const __m256 index_val = _mm256_loadu_ps((const float *)offsets);
    const __m256i index = *((__m256i *)(&index_val));
    _mm512_i32scatter_pd(base_ptr, index, data, 8);
  }

  VectorizedArray get_abs () const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +
    __m512d mask = _mm512_set1_pd (-0.);
    VectorizedArray res;
    res.data = (__m512d)_mm512_andnot_epi64 ((__m512i)mask, (__m512i)data);
    return res;
  }

  __m512d data;
};

template <>
inline
void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double            *in,
                              const unsigned int      *offsets,
                              VectorizedArray<double> *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int outer=0; outer<8; outer += 4)
    {
      const double *in0 = in + offsets[0+outer];
      const double *in1 = in + offsets[1+outer];
      const double *in2 = in + offsets[2+outer];
      const double *in3 = in + offsets[3+outer];

      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m256d u0 = _mm256_loadu_pd(in0+4*i);
          __m256d u1 = _mm256_loadu_pd(in1+4*i);
          __m256d u2 = _mm256_loadu_pd(in2+4*i);
          __m256d u3 = _mm256_loadu_pd(in3+4*i);
          __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
          __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
          __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
          __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
          *(__m256d *)((double *)(&out[4*i+0].data)+outer) = _mm256_unpacklo_pd (t0, t1);
          *(__m256d *)((double *)(&out[4*i+1].data)+outer) = _mm256_unpackhi_pd (t0, t1);
          *(__m256d *)((double *)(&out[4*i+2].data)+outer) = _mm256_unpacklo_pd (t2, t3);
          *(__m256d *)((double *)(&out[4*i+3].data)+outer) = _mm256_unpackhi_pd (t2, t3);
        }
      for (unsigned int i=4*n_chunks; i<n_entries; ++i)
        for (unsigned int v=0; v<4; ++v)
          out[i][outer+v] = in[offsets[v+outer]+i];
    }
}



/**
 * Specialization for double and AVX-512.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double> *in,
                               const unsigned int            *offsets,
                               double                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  // do not do full transpose because the code is too long and will most
  // likely not pay off. rather do the transposition on the vectorized array
  // on size smaller, mm256d
  for (unsigned int outer=0; outer<8; outer += 4)
    {
      double *out0 = out + offsets[0+outer];
      double *out1 = out + offsets[1+outer];
      double *out2 = out + offsets[2+outer];
      double *out3 = out + offsets[3+outer];
      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m256d u0 = *(const __m256d *)((const double *)(&in[4*i+0].data)+outer);
          __m256d u1 = *(const __m256d *)((const double *)(&in[4*i+1].data)+outer);
          __m256d u2 = *(const __m256d *)((const double *)(&in[4*i+2].data)+outer);
          __m256d u3 = *(const __m256d *)((const double *)(&in[4*i+3].data)+outer);
          __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
          __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
          __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
          __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
          __m256d res0 = _mm256_unpacklo_pd (t0, t1);
          __m256d res1 = _mm256_unpackhi_pd (t0, t1);
          __m256d res2 = _mm256_unpacklo_pd (t2, t3);
          __m256d res3 = _mm256_unpackhi_pd (t2, t3);

          // Cannot use the same store instructions in both paths of the 'if'
          // because the compiler cannot know that there is no aliasing between
          // pointers
          if (add_into)
            {
              res0 = _mm256_add_pd(_mm256_loadu_pd(out0+4*i), res0);
              _mm256_storeu_pd(out0+4*i, res0);
              res1 = _mm256_add_pd(_mm256_loadu_pd(out1+4*i), res1);
              _mm256_storeu_pd(out1+4*i, res1);
              res2 = _mm256_add_pd(_mm256_loadu_pd(out2+4*i), res2);
              _mm256_storeu_pd(out2+4*i, res2);
              res3 = _mm256_add_pd(_mm256_loadu_pd(out3+4*i), res3);
              _mm256_storeu_pd(out3+4*i, res3);
            }
          else
            {
              _mm256_storeu_pd(out0+4*i, res0);
              _mm256_storeu_pd(out1+4*i, res1);
              _mm256_storeu_pd(out2+4*i, res2);
              _mm256_storeu_pd(out3+4*i, res3);
            }
        }
      if (add_into)
        for (unsigned int i=4*n_chunks; i<n_entries; ++i)
          for (unsigned int v=0; v<4; ++v)
            out[offsets[v+outer]+i] += in[i][v+outer];
      else
        for (unsigned int i=4*n_chunks; i<n_entries; ++i)
          for (unsigned int v=0; v<4; ++v)
            out[offsets[v+outer]+i] = in[i][v+outer];
    }
}



/**
 * Specialization for float and AVX512.
 */
template<>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 16;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const float x)
  {
    data = _mm512_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm512_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm512_sub_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm512_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm512_div_ps(data,vec.data);
#endif
    return *this;
  }

  void load (const float *ptr)
  {
    data = _mm512_loadu_ps (ptr);
  }

  void store (float *ptr) const
  {
    _mm512_storeu_ps (ptr, data);
  }

  void streaming_store (float *ptr) const
  {
    _mm512_stream_ps(ptr,data);
  }

  void gather (const float        *base_ptr,
               const unsigned int *offsets)
  {
    const __m512 index_val = _mm512_loadu_ps((const float *)offsets);
    const __m512i index = *((__m512i *)(&index_val));
    data = _mm512_i32gather_ps(index, base_ptr, 4);
  }

  void scatter (const unsigned int *offsets,
                float              *base_ptr) const
  {
    const __m512 index_val = _mm512_loadu_ps((const float *)offsets);
    const __m512i index = *((__m512i *)(&index_val));
    _mm512_i32scatter_ps(base_ptr, index, data, 4);
  }

  VectorizedArray get_abs () const
  {
    __m512 mask = _mm512_set1_ps (-0.f);
    VectorizedArray res;
    res.data = (__m512)_mm512_andnot_epi32 ((__m512i)mask, (__m512i)data);
    return res;
  }

  __m512 data;
};

template <>
inline
void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float            *in,
                              const unsigned int     *offsets,
                              VectorizedArray<float> *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int outer = 0; outer<16; outer += 8)
    {
      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m128 u0 = _mm_loadu_ps(in+4*i+offsets[0+outer]);
          __m128 u1 = _mm_loadu_ps(in+4*i+offsets[1+outer]);
          __m128 u2 = _mm_loadu_ps(in+4*i+offsets[2+outer]);
          __m128 u3 = _mm_loadu_ps(in+4*i+offsets[3+outer]);
          __m128 u4 = _mm_loadu_ps(in+4*i+offsets[4+outer]);
          __m128 u5 = _mm_loadu_ps(in+4*i+offsets[5+outer]);
          __m128 u6 = _mm_loadu_ps(in+4*i+offsets[6+outer]);
          __m128 u7 = _mm_loadu_ps(in+4*i+offsets[7+outer]);
          // To avoid warnings about uninitialized variables, need to initialize
          // one variable with zero before using it.
          __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
          t0 = _mm256_insertf128_ps (t3, u0, 0);
          t0 = _mm256_insertf128_ps (t0, u4, 1);
          t1 = _mm256_insertf128_ps (t3, u1, 0);
          t1 = _mm256_insertf128_ps (t1, u5, 1);
          t2 = _mm256_insertf128_ps (t3, u2, 0);
          t2 = _mm256_insertf128_ps (t2, u6, 1);
          t3 = _mm256_insertf128_ps (t3, u3, 0);
          t3 = _mm256_insertf128_ps (t3, u7, 1);
          __m256 v0 = _mm256_shuffle_ps (t0, t1, 0x44);
          __m256 v1 = _mm256_shuffle_ps (t0, t1, 0xee);
          __m256 v2 = _mm256_shuffle_ps (t2, t3, 0x44);
          __m256 v3 = _mm256_shuffle_ps (t2, t3, 0xee);
          *(__m256 *)((float *)(&out[4*i+0].data)+outer) = _mm256_shuffle_ps (v0, v2, 0x88);
          *(__m256 *)((float *)(&out[4*i+1].data)+outer) = _mm256_shuffle_ps (v0, v2, 0xdd);
          *(__m256 *)((float *)(&out[4*i+2].data)+outer) = _mm256_shuffle_ps (v1, v3, 0x88);
          *(__m256 *)((float *)(&out[4*i+3].data)+outer) = _mm256_shuffle_ps (v1, v3, 0xdd);
        }
      for (unsigned int i=4*n_chunks; i<n_entries; ++i)
        for (unsigned int v=0; v<8; ++v)
          out[i][v+outer] = in[offsets[v+outer]+i];
    }
}



/**
 * Specialization for float and AVX-512.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float> *in,
                               const unsigned int           *offsets,
                               float                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int outer = 0; outer<16; outer += 8)
    {
      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m256 u0 = *(const __m256 *)((const float *)(&in[4*i+0].data)+outer);
          __m256 u1 = *(const __m256 *)((const float *)(&in[4*i+1].data)+outer);
          __m256 u2 = *(const __m256 *)((const float *)(&in[4*i+2].data)+outer);
          __m256 u3 = *(const __m256 *)((const float *)(&in[4*i+3].data)+outer);
          __m256 t0 = _mm256_shuffle_ps (u0, u1, 0x44);
          __m256 t1 = _mm256_shuffle_ps (u0, u1, 0xee);
          __m256 t2 = _mm256_shuffle_ps (u2, u3, 0x44);
          __m256 t3 = _mm256_shuffle_ps (u2, u3, 0xee);
          u0 = _mm256_shuffle_ps (t0, t2, 0x88);
          u1 = _mm256_shuffle_ps (t0, t2, 0xdd);
          u2 = _mm256_shuffle_ps (t1, t3, 0x88);
          u3 = _mm256_shuffle_ps (t1, t3, 0xdd);
          __m128 res0 = _mm256_extractf128_ps (u0, 0);
          __m128 res4 = _mm256_extractf128_ps (u0, 1);
          __m128 res1 = _mm256_extractf128_ps (u1, 0);
          __m128 res5 = _mm256_extractf128_ps (u1, 1);
          __m128 res2 = _mm256_extractf128_ps (u2, 0);
          __m128 res6 = _mm256_extractf128_ps (u2, 1);
          __m128 res3 = _mm256_extractf128_ps (u3, 0);
          __m128 res7 = _mm256_extractf128_ps (u3, 1);

          // Cannot use the same store instructions in both paths of the 'if'
          // because the compiler cannot know that there is no aliasing between
          // pointers
          if (add_into)
            {
              res0 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[0+outer]), res0);
              _mm_storeu_ps(out+4*i+offsets[0+outer], res0);
              res1 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[1+outer]), res1);
              _mm_storeu_ps(out+4*i+offsets[1+outer], res1);
              res2 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[2+outer]), res2);
              _mm_storeu_ps(out+4*i+offsets[2+outer], res2);
              res3 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[3+outer]), res3);
              _mm_storeu_ps(out+4*i+offsets[3+outer], res3);
              res4 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[4+outer]), res4);
              _mm_storeu_ps(out+4*i+offsets[4+outer], res4);
              res5 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[5+outer]), res5);
              _mm_storeu_ps(out+4*i+offsets[5+outer], res5);
              res6 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[6+outer]), res6);
              _mm_storeu_ps(out+4*i+offsets[6+outer], res6);
              res7 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[7+outer]), res7);
              _mm_storeu_ps(out+4*i+offsets[7+outer], res7);
            }
          else
            {
              _mm_storeu_ps(out+4*i+offsets[0+outer], res0);
              _mm_storeu_ps(out+4*i+offsets[1+outer], res1);
              _mm_storeu_ps(out+4*i+offsets[2+outer], res2);
              _mm_storeu_ps(out+4*i+offsets[3+outer], res3);
              _mm_storeu_ps(out+4*i+offsets[4+outer], res4);
              _mm_storeu_ps(out+4*i+offsets[5+outer], res5);
              _mm_storeu_ps(out+4*i+offsets[6+outer], res6);
              _mm_storeu_ps(out+4*i+offsets[7+outer], res7);
            }
        }
      if (add_into)
        for (unsigned int i=4*n_chunks; i<n_entries; ++i)
          for (unsigned int v=0; v<8; ++v)
            out[offsets[v+outer]+i] += in[i][v+outer];
      else
        for (unsigned int i=4*n_chunks; i<n_entries; ++i)
          for (unsigned int v=0; v<8; ++v)
            out[offsets[v+outer]+i] = in[i][v+outer];
    }
}


#elif defined(__AVX__)

/**
 * Specialization of VectorizedArray class for double and AVX.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const double x)
  {
    data = _mm256_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm256_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm256_sub_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm256_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm256_div_pd(data,vec.data);
#endif
    return *this;
  }

  void load (const double *ptr)
  {
    data = _mm256_loadu_pd (ptr);
  }

  void store (double *ptr) const
  {
    _mm256_storeu_pd (ptr, data);
  }

  void streaming_store (double *ptr) const
  {
    _mm256_stream_pd(ptr,data);
  }

  void gather (const double       *base_ptr,
               const unsigned int *offsets)
  {
#if __AVX2__
    const __m128 index_val = _mm_loadu_ps((const float *)offsets);
    const __m128i index = *((__m128i *)(&index_val));
    data = _mm256_i32gather_pd(base_ptr, index, 8);
#else
    for (unsigned int i=0; i<4; ++i)
      *(reinterpret_cast<double *>(&data)+i) = base_ptr[offsets[i]];
#endif
  }

  void scatter (const unsigned int *offsets,
                double             *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i=0; i<4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data)+i);
  }

  VectorizedArray get_abs () const
  {
    __m256d mask = _mm256_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm256_andnot_pd(mask, data);
    return res;
  }

  __m256d data;
};

template <>
inline
void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double            *in,
                              const unsigned int      *offsets,
                              VectorizedArray<double> *out)
{
  const unsigned int n_chunks = n_entries/4;
  const double *in0 = in + offsets[0];
  const double *in1 = in + offsets[1];
  const double *in2 = in + offsets[2];
  const double *in3 = in + offsets[3];

  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256d u0 = _mm256_loadu_pd(in0+4*i);
      __m256d u1 = _mm256_loadu_pd(in1+4*i);
      __m256d u2 = _mm256_loadu_pd(in2+4*i);
      __m256d u3 = _mm256_loadu_pd(in3+4*i);
      __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
      __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
      __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
      __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
      out[4*i+0].data = _mm256_unpacklo_pd (t0, t1);
      out[4*i+1].data = _mm256_unpackhi_pd (t0, t1);
      out[4*i+2].data = _mm256_unpacklo_pd (t2, t3);
      out[4*i+3].data = _mm256_unpackhi_pd (t2, t3);
    }
  for (unsigned int i=4*n_chunks; i<n_entries; ++i)
    for (unsigned int v=0; v<4; ++v)
      out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double> *in,
                               const unsigned int            *offsets,
                               double                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  double *out0 = out + offsets[0];
  double *out1 = out + offsets[1];
  double *out2 = out + offsets[2];
  double *out3 = out + offsets[3];
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256d u0 = in[4*i+0].data;
      __m256d u1 = in[4*i+1].data;
      __m256d u2 = in[4*i+2].data;
      __m256d u3 = in[4*i+3].data;
      __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
      __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
      __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
      __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
      __m256d res0 = _mm256_unpacklo_pd (t0, t1);
      __m256d res1 = _mm256_unpackhi_pd (t0, t1);
      __m256d res2 = _mm256_unpacklo_pd (t2, t3);
      __m256d res3 = _mm256_unpackhi_pd (t2, t3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out0+4*i), res0);
          _mm256_storeu_pd(out0+4*i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out1+4*i), res1);
          _mm256_storeu_pd(out1+4*i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out2+4*i), res2);
          _mm256_storeu_pd(out2+4*i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out3+4*i), res3);
          _mm256_storeu_pd(out3+4*i, res3);
        }
      else
        {
          _mm256_storeu_pd(out0+4*i, res0);
          _mm256_storeu_pd(out1+4*i, res1);
          _mm256_storeu_pd(out2+4*i, res2);
          _mm256_storeu_pd(out3+4*i, res3);
        }
    }
  if (add_into)
    for (unsigned int i=4*n_chunks; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=4*n_chunks; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] = in[i][v];
}


/**
 * Specialization for float and AVX.
 */
template<>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 8;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const float x)
  {
    data = _mm256_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm256_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm256_sub_ps(data,vec.data);
#endif
    return *this;
  }

  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm256_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm256_div_ps(data,vec.data);
#endif
    return *this;
  }

  void load (const float *ptr)
  {
    data = _mm256_loadu_ps (ptr);
  }

  void store (float *ptr) const
  {
    _mm256_storeu_ps (ptr, data);
  }

  void streaming_store (float *ptr) const
  {
    _mm256_stream_ps(ptr,data);
  }

  void gather (const float        *base_ptr,
               const unsigned int *offsets)
  {
#if __AVX2__
    const __m256 index_val = _mm256_loadu_ps((const float *)offsets);
    const __m256i index = *((__m256i *)(&index_val));
    data = _mm256_i32gather_ps(base_ptr, index, 4);
#else
    for (unsigned int i=0; i<8; ++i)
      *(reinterpret_cast<float *>(&data)+i) = base_ptr[offsets[i]];
#endif
  }

  void scatter (const unsigned int *offsets,
                float              *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i=0; i<8; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data)+i);
  }

  VectorizedArray get_abs () const
  {
    __m256 mask = _mm256_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm256_andnot_ps(mask, data);
    return res;
  }

  __m256 data;
};

template <>
inline
void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float            *in,
                              const unsigned int     *offsets,
                              VectorizedArray<float> *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m128 u0 = _mm_loadu_ps(in+4*i+offsets[0]);
      __m128 u1 = _mm_loadu_ps(in+4*i+offsets[1]);
      __m128 u2 = _mm_loadu_ps(in+4*i+offsets[2]);
      __m128 u3 = _mm_loadu_ps(in+4*i+offsets[3]);
      __m128 u4 = _mm_loadu_ps(in+4*i+offsets[4]);
      __m128 u5 = _mm_loadu_ps(in+4*i+offsets[5]);
      __m128 u6 = _mm_loadu_ps(in+4*i+offsets[6]);
      __m128 u7 = _mm_loadu_ps(in+4*i+offsets[7]);
      // To avoid warnings about uninitialized variables, need to initialize
      // one variable with zero before using it.
      __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
      t0 = _mm256_insertf128_ps (t3, u0, 0);
      t0 = _mm256_insertf128_ps (t0, u4, 1);
      t1 = _mm256_insertf128_ps (t3, u1, 0);
      t1 = _mm256_insertf128_ps (t1, u5, 1);
      t2 = _mm256_insertf128_ps (t3, u2, 0);
      t2 = _mm256_insertf128_ps (t2, u6, 1);
      t3 = _mm256_insertf128_ps (t3, u3, 0);
      t3 = _mm256_insertf128_ps (t3, u7, 1);
      __m256 v0 = _mm256_shuffle_ps (t0, t1, 0x44);
      __m256 v1 = _mm256_shuffle_ps (t0, t1, 0xee);
      __m256 v2 = _mm256_shuffle_ps (t2, t3, 0x44);
      __m256 v3 = _mm256_shuffle_ps (t2, t3, 0xee);
      out[4*i+0].data = _mm256_shuffle_ps (v0, v2, 0x88);
      out[4*i+1].data = _mm256_shuffle_ps (v0, v2, 0xdd);
      out[4*i+2].data = _mm256_shuffle_ps (v1, v3, 0x88);
      out[4*i+3].data = _mm256_shuffle_ps (v1, v3, 0xdd);
    }
  for (unsigned int i=4*n_chunks; i<n_entries; ++i)
    for (unsigned int v=0; v<8; ++v)
      out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for float and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float> *in,
                               const unsigned int           *offsets,
                               float                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256 u0 = in[4*i+0].data;
      __m256 u1 = in[4*i+1].data;
      __m256 u2 = in[4*i+2].data;
      __m256 u3 = in[4*i+3].data;
      __m256 t0 = _mm256_shuffle_ps (u0, u1, 0x44);
      __m256 t1 = _mm256_shuffle_ps (u0, u1, 0xee);
      __m256 t2 = _mm256_shuffle_ps (u2, u3, 0x44);
      __m256 t3 = _mm256_shuffle_ps (u2, u3, 0xee);
      u0 = _mm256_shuffle_ps (t0, t2, 0x88);
      u1 = _mm256_shuffle_ps (t0, t2, 0xdd);
      u2 = _mm256_shuffle_ps (t1, t3, 0x88);
      u3 = _mm256_shuffle_ps (t1, t3, 0xdd);
      __m128 res0 = _mm256_extractf128_ps (u0, 0);
      __m128 res4 = _mm256_extractf128_ps (u0, 1);
      __m128 res1 = _mm256_extractf128_ps (u1, 0);
      __m128 res5 = _mm256_extractf128_ps (u1, 1);
      __m128 res2 = _mm256_extractf128_ps (u2, 0);
      __m128 res6 = _mm256_extractf128_ps (u2, 1);
      __m128 res3 = _mm256_extractf128_ps (u3, 0);
      __m128 res7 = _mm256_extractf128_ps (u3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[0]), res0);
          _mm_storeu_ps(out+4*i+offsets[0], res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[1]), res1);
          _mm_storeu_ps(out+4*i+offsets[1], res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[2]), res2);
          _mm_storeu_ps(out+4*i+offsets[2], res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[3]), res3);
          _mm_storeu_ps(out+4*i+offsets[3], res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[4]), res4);
          _mm_storeu_ps(out+4*i+offsets[4], res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[5]), res5);
          _mm_storeu_ps(out+4*i+offsets[5], res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[6]), res6);
          _mm_storeu_ps(out+4*i+offsets[6], res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[7]), res7);
          _mm_storeu_ps(out+4*i+offsets[7], res7);
        }
      else
        {
          _mm_storeu_ps(out+4*i+offsets[0], res0);
          _mm_storeu_ps(out+4*i+offsets[1], res1);
          _mm_storeu_ps(out+4*i+offsets[2], res2);
          _mm_storeu_ps(out+4*i+offsets[3], res3);
          _mm_storeu_ps(out+4*i+offsets[4], res4);
          _mm_storeu_ps(out+4*i+offsets[5], res5);
          _mm_storeu_ps(out+4*i+offsets[6], res6);
          _mm_storeu_ps(out+4*i+offsets[7], res7);
        }
    }
  if (add_into)
    for (unsigned int i=4*n_chunks; i<n_entries; ++i)
      for (unsigned int v=0; v<8; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=4*n_chunks; i<n_entries; ++i)
      for (unsigned int v=0; v<8; ++v)
        out[offsets[v]+i] = in[i][v];
}


#elif defined(__SSE2__)

/**
 * Specialization for double and SSE2.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 2;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const double x)
  {
    data = _mm_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm_sub_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm_div_pd(data,vec.data);
#endif
    return *this;
  }

  void load (const double *ptr)
  {
    data = _mm_loadu_pd (ptr);
  }

  void store (double *ptr) const
  {
    _mm_storeu_pd (ptr, data);
  }

  void streaming_store (double *ptr) const
  {
    _mm_stream_pd(ptr,data);
  }

  void gather (const double       *base_ptr,
               const unsigned int *offsets)
  {
    for (unsigned int i=0; i<2; ++i)
      *(reinterpret_cast<double *>(&data)+i) = base_ptr[offsets[i]];
  }

  void scatter (const unsigned int *offsets,
                double             *base_ptr) const
  {
    for (unsigned int i=0; i<2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data)+i);
  }

  VectorizedArray get_abs () const
  {
    __m128d mask = _mm_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm_andnot_pd(mask, data);
    return res;
  }

  __m128d data;
};



/**
 * Specialization for float and SSE2.
 */
template <>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const float x)
  {
    data = _mm_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm_sub_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#if USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm_div_ps(data,vec.data);
#endif
    return *this;
  }

  void load (const float *ptr)
  {
    data = _mm_loadu_ps (ptr);
  }

  void store (float *ptr) const
  {
    _mm_storeu_ps (ptr, data);
  }

  void streaming_store (float *ptr) const
  {
    _mm_stream_ps(ptr,data);
  }

  void gather (const float        *base_ptr,
               const unsigned int *offsets)
  {
    for (unsigned int i=0; i<4; ++i)
      *(reinterpret_cast<float *>(&data)+i) = base_ptr[offsets[i]];
  }

  void scatter (const unsigned int *offsets,
                float              *base_ptr) const
  {
    for (unsigned int i=0; i<4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data)+i);
  }

  VectorizedArray get_abs () const
  {
    __m128 mask = _mm_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm_andnot_ps(mask, data);
    return res;
  }

  __m128 data;
};


#endif // __SSE2__


template <typename Number>
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp+=v;
}

template <typename Number>
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp-=v;
}

template <typename Number>
VectorizedArray<Number>
operator * (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp*=v;
}

template <typename Number>
VectorizedArray<Number>
operator / (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp/=v;
}

template <typename Number>
VectorizedArray<Number>
operator + (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp+=v;
}

VectorizedArray<float>
operator + (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = u;
  return tmp+=v;
}

template <typename Number>
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  return u + v;
}

VectorizedArray<float>
operator + (const VectorizedArray<float> &v,
            const double                 &u)
{
  return u + v;
}

template <typename Number>
VectorizedArray<Number>
operator - (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp-=v;
}

VectorizedArray<float>
operator - (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp-=v;
}

template <typename Number>
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v-tmp;
}

VectorizedArray<float>
operator - (const VectorizedArray<float> &v,
            const double                 &u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v-tmp;
}

template <typename Number>
VectorizedArray<Number>
operator * (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp*=v;
}

VectorizedArray<float>
operator * (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp*=v;
}

template <typename Number>
VectorizedArray<Number>
operator * (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  return u * v;
}

VectorizedArray<float>
operator * (const VectorizedArray<float> &v,
            const double                 &u)
{
  return u * v;
}

template <typename Number>
VectorizedArray<Number>
operator / (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp/=v;
}

VectorizedArray<float>
operator / (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp/=v;
}

template <typename Number>
VectorizedArray<Number>
operator / (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v/tmp;
}

VectorizedArray<float>
operator / (const VectorizedArray<float> &v,
            const double                 &u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v/tmp;
}

template <typename Number>
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &u)
{
  return u;
}

template <typename Number>
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &u)
{
  // to get a negative sign, subtract the input from zero (could also
  // multiply by -1, but this one is slightly simpler)
  return VectorizedArray<Number>()-u;
}


namespace std
{
  template <typename Number>
  inline
  VectorizedArray<Number>
  abs (const VectorizedArray<Number> &x)
  {
    return x.get_abs();
  }
}

#endif
