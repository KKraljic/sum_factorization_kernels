// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// This file is inspired by the file
// deal.II/include/deal.II/base/aligned_vector.h, see www.dealii.org for
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

#ifndef aligned_vector_h
#define aligned_vector_h

#include <limits>
#include <vector>
#include <iostream>

/**
 * Reduced aligned vector class. Very basic logic for resize, i.e., entries
 * are deleted upon resize, no re-use of memory allocations. No logic for
 * avoiding out-of-bound access.
 */
template <typename Number>
class AlignedVector
{
public:
  AlignedVector()
    :
    val (nullptr),
    my_size (0)
  {}

  AlignedVector(const std::size_t size, const Number init = Number())
    :
    val (nullptr),
    my_size (0)
  {
    resize(size);
  }

  ~AlignedVector()
  {
    if (val != nullptr)
      free(val);
  }

  void resize(const std::size_t size, const Number init = Number())
  {
    if (val != nullptr)
      {
        free(val);
        val = nullptr;
      }
    resize_fast(size);
    for (unsigned int i=0; i<size; ++i)
      val[i] = init;
  }

  void resize_fast(const std::size_t size)
  {
    if (val != nullptr)
      throw;
    if (size > 0)
      {
        int err = posix_memalign((void **)&val, 64, size * sizeof(Number));
        if (err != 0 || val == nullptr)
          {
            std::cout << "Out of memory" << std::endl;
            throw;
          }
       }
    my_size = size;
  }

  Number* begin()
  {
    return val;
  }

  const Number * begin() const
  {
    return val;
  }

  const Number* end() const
  {
    return val+my_size;
  }

  const Number operator[] (const std::size_t index) const
  {
    return val[index];
  }

  Number& operator[] (const std::size_t index)
  {
    return val[index];
  }

  std::size_t size() const
  {
    return my_size;
  }

private:
  Number *val;
  std::size_t my_size;
};

#endif
