// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// This file is inspired by the file deal.II/include/deal.II/base/utilities.h,
// see www.dealii.org for information about licenses. Here is the original
// deal.II license statement:
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

#ifndef utilities_h
#define utilities_h


#if defined (__GNUC__)
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define ALWAYS_INLINE
#endif


/**
 * A replacement for <code>std::pow</code> that allows compile-time
 * calculations for constant expression arguments.
 */
namespace Utilities
{
  constexpr
  unsigned int
  pow(const unsigned int base, const unsigned int iexp)
  {
    return iexp == 0 ? 1 : base*::Utilities::pow(base, iexp - 1);
  }
}


#endif
