// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// implementation of cell and face terms for DG Laplacian (interior penalty
// method) using integration on Cartesian cell geometries with integration
//
// Author: Martin Kronbichler, April 2018

#ifndef evaluation_dg_laplacian_h
#define evaluation_dg_laplacian_h

#include <mpi.h>

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


template<int dim, int degree, typename Number>
class EvaluationDGLaplacian {
public:
    static constexpr unsigned int dimension = dim;
    static constexpr unsigned int n_q_points = Utilities::pow(degree + 1, dim);
    static constexpr unsigned int dofs_per_cell = Utilities::pow(degree + 1, dim);
    unsigned int blx;
    unsigned int bly;
    unsigned int blz;

    void initialize(const unsigned int *n_cells_in) {
        n_cells[0] = n_cells_in[0] / VectorizedArray<Number>::n_array_elements;
        for (unsigned int d = 1; d < dim; ++d)
            n_cells[d] = n_cells_in[d];
        for (unsigned int d = dim; d < 3; ++d)
            n_cells[d] = 1;

        n_blocks[2] = (n_cells[2] + blz - 1) / blz;
        n_blocks[1] = (n_cells[1] + bly - 1) / bly;
        n_blocks[0] = (n_cells[0] + blx - 1) / blx;

        sol_old.resize(0);
        sol_new.resize(0);
        mat_diagonal.resize(0);
        sol_tmp.resize(0);
        sol_rhs.resize(0);

        sol_old.resize_fast(n_elements() * dofs_per_cell);
        sol_new.resize_fast(n_elements() * dofs_per_cell);
        mat_diagonal.resize_fast(n_elements() * dofs_per_cell);
        sol_tmp.resize_fast(n_elements() * dofs_per_cell);
        sol_rhs.resize_fast(n_elements() * dofs_per_cell);

#pragma omp parallel firstprivate() reduction(+: )
        {
#pragma omp for schedule (static) collapse(2)
            for (unsigned int ib = 0; ib < n_blocks[2]; ++ib)
                for (unsigned int jb = 0; jb < n_blocks[1]; ++jb)
                    for (unsigned int kb = 0; kb < n_blocks[0]; ++kb)
                        for (unsigned int i = ib * blz; i < std::min(n_cells[2], (ib + 1) * blz); ++i)
                            for (unsigned int j = jb * bly; j < std::min(n_cells[1], (jb + 1) * bly); ++j) {
                                const unsigned int ii = (i * n_cells[1] + j) * n_cells[0];
                                for (std::size_t ix =
                                        dofs_per_cell * VectorizedArray<Number>::n_array_elements * (kb * blx + ii);
                                     ix < (std::min(n_cells[0], (kb + 1) * blx) + ii) * dofs_per_cell *
                                          VectorizedArray<Number>::n_array_elements; ++ix) {
                                    sol_old[ix] = 1;
                                    sol_new[ix] = 0.;
                                    mat_diagonal[ix] = 1.;
                                    sol_tmp[ix] = 0.;
                                    sol_rhs[ix] = 1;
                                }
                            }
        }

        std::vector<double> jacobian = get_diagonal_jacobian();
        Number jacobian_determinant = 1.;
        for (unsigned int d = 0; d < dim; ++d)
            jacobian_determinant *= jacobian[d];
        jacobian_determinant = 1. / jacobian_determinant;

        jxw_data.resize(1);
        jxw_data[0] = jacobian_determinant;
        jacobian_data.resize(dim);
        for (unsigned int d = 0; d < dim; ++d)
            jacobian_data[d] = jacobian[d];

        fill_shape_values();
    }

    std::size_t n_elements() const {
        std::size_t n_element = VectorizedArray<Number>::n_array_elements;
        for (unsigned int d = 0; d < dim; ++d)
            n_element *= n_cells[d];
        return n_element;
    }


private:
    std::vector<double> get_diagonal_jacobian() const {
        std::vector<double> jacobian(dim);
        jacobian[0] = 4.;
        for (unsigned int d = 1; d < dim; ++d) {
            double entry = (double) (d + 1) / 4.;
            jacobian[d] = 1. / entry;
        }
        return jacobian;
    }

    unsigned int n_cells[3];
    unsigned int n_blocks[3];
    AlignedVector<VectorizedArray<Number> > shape_values_eo;
    AlignedVector<VectorizedArray<Number> > shape_gradients_eo;
    AlignedVector<VectorizedArray<Number> > shape_values_on_face_eo;

    AlignedVector<Number> quadrature_weights;
    AlignedVector<Number> face_quadrature_weight;

    VectorizedArray<Number> hermite_derivative_on_face;

    AlignedVector<Number> sol_old;
    AlignedVector<Number> sol_new;
    AlignedVector<Number> sol_rhs;
    AlignedVector<Number> sol_tmp;
    AlignedVector<Number> mat_diagonal;

    AlignedVector<VectorizedArray<Number> > jxw_data;
    AlignedVector<VectorizedArray<Number> > jacobian_data;
};


#endif
