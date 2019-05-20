// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// This file is inspired by the file
// deal.II/include/deal.II/matrix_free/tensor_product_kernels.h, see
// www.dealii.org for information about licenses. Here is the original deal.II
// license statement:
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


#ifndef matrix_vector_kernel_h
#define matrix_vector_kernel_h


template<int nn, int stride, int type, bool contract_over_rows, bool add_into,
        typename Number, typename Number2=Number, bool nontemporal_store = false, int do_dg = 0>
inline ALWAYS_INLINE
void
apply_1d_matvec_kernel(const Number2 *__restrict coefficients_eo,
                       const Number *in,
                       Number *out,
                       const Number *array_for_add = nullptr,
                       const Number2 *__restrict dg_coefficients = nullptr,
                       Number *array_face = nullptr) {
    const unsigned int mid = nn / 2;
    const unsigned int offset = (nn + 1) / 2;
    Number xp[mid > 0 ? mid : 1], xm[mid > 0 ? mid : 1];
    for (unsigned int i = 0; i < mid; ++i) {
        xp[i] = in[i * stride] + in[(nn - 1 - i) * stride];
        xm[i] = in[i * stride] - in[(nn - 1 - i) * stride];
    }


    Number xmid = in[stride * mid];

    for (unsigned int col = 0; col < mid; ++col) {
        Number r0, r1;
        r0 = coefficients_eo[col] * xp[0];
        r1 = coefficients_eo[(nn - 1) * offset + col] * xm[0];

        for (unsigned int ind = 1; ind < mid; ++ind) {
            r0 += coefficients_eo[ind * offset + col] * xp[ind];
            r1 += coefficients_eo[(nn - 1 - ind) * offset + col] * xm[ind];
        }
        if (nn % 2 == 1) {
            r0 += coefficients_eo[mid * offset + col] * xmid;
        }

        Number t = r0;
        r0 = (add_into ? array_for_add[col * stride] + t : t) + r1;
        r1 = (add_into ? array_for_add[(nn - 1 - col) * stride] + t : t) - r1;

        out[col * stride] = r0;
        out[(nn - 1 - col) * stride] = r1;

    }
    if (nn % 2 == 1) {
        Number r0 = (add_into ?
                     array_for_add[mid] + coefficients_eo[mid * offset + mid] * xmid :
                     coefficients_eo[mid * offset + mid] * xmid);
        for (unsigned int ind = 0; ind < mid; ++ind)
            r0 += coefficients_eo[ind * offset + mid] * xp[ind];

        out[mid * stride] = r0;
    }
}


#endif
