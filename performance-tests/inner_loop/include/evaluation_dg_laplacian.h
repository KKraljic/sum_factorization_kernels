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

#pragma omp parallel
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


    //TODO: "Inline"
    void do_chebyshev() {
        do_cheby_iter<true>(sol_rhs, sol_old, mat_diagonal, sol_new, sol_tmp, 0.5, 0.5);
    }

    void do_matvec() {
        do_cheby_iter<false>(sol_new, sol_new, sol_new, sol_tmp, sol_tmp, 0., 0.);
    }

    //TODO: Mock subsequent calls
    template<bool evaluate_chebyshev = true>
    void do_inner_loop(const unsigned int start_x,
                       const unsigned int end_x,
                       const unsigned int iy,
                       const unsigned int iz,
                       const AlignedVector<Number> &src,
                       const AlignedVector<Number> &sol_old,
                       const AlignedVector<Number> &mat_diagonal,
                       AlignedVector<Number> &sol_new,
                       AlignedVector<Number> &vec_tm,
                       const Number coefficient_np,
                       const Number coefficient_tm) {

        constexpr unsigned int nn = degree + 1;
        constexpr unsigned int n_lanes = VectorizedArray<Number>::n_array_elements;
        constexpr unsigned int dofs_per_face = Utilities::pow(degree + 1, dim - 1);
        constexpr unsigned int dofs_per_plane = Utilities::pow(degree + 1, 2);
        const VectorizedArray<Number> *__restrict shape_values_eo = this->shape_values_eo.begin();
        const VectorizedArray<Number> *__restrict shape_gradients_eo = this->shape_gradients_eo.begin();
        AlignedVector<VectorizedArray<Number> > scratch_data_array;
        VectorizedArray<Number> my_array[degree < 13 ? 2 * dofs_per_cell : 1];
        VectorizedArray<Number> *__restrict data_ptr;
        VectorizedArray<Number> array_f[6][dofs_per_face], array_fd[6][dofs_per_face];

        if (degree < 13) {
            data_ptr = my_array;
        } else {
            scratch_data_array.resize_fast(2 * dofs_per_cell);
            data_ptr = scratch_data_array.begin();
        }

        //Outer Loop
        for (unsigned int ix = start_x; ix < end_x; ++ix) {
            const unsigned int ii = ((iz * n_cells[1] + iy) * n_cells[0] + ix) * n_lanes;
            const VectorizedArray<Number> *src_array =
                    reinterpret_cast<const VectorizedArray<Number> *>(sol_old.begin() + ii * dofs_per_cell);
            VectorizedArray<Number> *dst_array =
                    reinterpret_cast<VectorizedArray<Number> *>(sol_new.begin() + ii * dofs_per_cell);

            const VectorizedArray<Number> *inv_jac = jacobian_data.begin();
            const VectorizedArray<Number> my_jxw = jxw_data[0];

            for (unsigned int i2 = 0; i2 < (dim > 2 ? degree + 1 : 1); ++i2) {
                // x-direction
                VectorizedArray<Number> *__restrict in = data_ptr + i2 * nn * nn;
                for (unsigned int i1 = 0; i1 < nn; ++i1) {
                    apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                            (shape_values_eo, src_array + i2 * nn * nn + i1 * nn, in + i1 * nn);
                }
                // y-direction
                for (unsigned int i1 = 0; i1 < nn; ++i1) {
                    apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                            (shape_values_eo, in + i1, in + i1);
                }
            }

            const double penalty_factor = 1.;
            if (dim == 3) {
                const unsigned int index[2] = {(iz > 0 ?
                                                (ii - n_cells[1] * n_cells[0] * n_lanes) :
                                                (ii + (n_cells[2] - 1) * n_cells[1] * n_cells[0] * n_lanes)
                                               ) * dofs_per_cell,
                                               (iz < n_cells[2] - 1 ?
                                                (ii + n_cells[1] * n_cells[0] * n_lanes) :
                                                (ii - (n_cells[2] - 1) * n_cells[1] * n_cells[0] * n_lanes)
                                               ) * dofs_per_cell};
                //Exactly 2 passes: f = 4 and f = 5
                for (unsigned int f = 4; f < 6; ++f) {
                    const VectorizedArray<Number> w0 = (f == 4 ? 1. : -1.) * hermite_derivative_on_face;
                    const unsigned int offset1 = (f == 4 ? dofs_per_face * degree : 0);
                    const unsigned int offset2 = dofs_per_face * (f == 4 ? degree - 1 : 1);
                    for (unsigned int i = 0; i < dofs_per_face; ++i) {
                        array_f[f][i].load(sol_old.begin() + index[f % 2] + (offset1 + i) * n_lanes);
                        array_fd[f][i].load(sol_old.begin() + index[f % 2] + (offset2 + i) * n_lanes);
                        array_fd[f][i] = w0 * (array_fd[f][i] - array_f[f][i]);
                    }

                    // interpolate values onto quadrature points
                    for (unsigned int i1 = 0; i1 < nn; ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1 * nn, array_f[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < nn; ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1, array_f[f] + i1);

                    // interpolate derivatives onto quadrature points
                    for (unsigned int i1 = 0; i1 < nn; ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1 * nn, array_fd[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < nn; ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1, array_fd[f] + i1);
                }

                const VectorizedArray<Number> tau = (degree + 1) * (degree + 1) * penalty_factor * inv_jac[2];
                const VectorizedArray<Number> JxW_face = my_jxw * inv_jac[2];
                for (unsigned int i2 = 0; i2 < dofs_per_face; ++i2) {
                    //TODO: Whut?
                    apply_1d_matvec_kernel<degree + 1, dofs_per_face, 0, true, false, VectorizedArray<Number>>
                            (shape_values_eo, data_ptr + i2, data_ptr + i2);

                    // include evaluation from this cell onto face
                    VectorizedArray<Number> array_face[4];
                    apply_1d_matvec_kernel<degree +
                                           1, dofs_per_face, 1, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
                            (shape_gradients_eo, data_ptr + i2, data_ptr + dofs_per_cell + i2,
                             nullptr, shape_values_on_face_eo.begin(), array_face);

                    // face integrals in z direction
                    {
                        const VectorizedArray<Number> outval = array_f[4][i2];
                        const VectorizedArray<Number> outgrad = array_fd[4][i2];
                        const VectorizedArray<Number> avg_grad = 0.5 * inv_jac[2] * (array_face[2] + outgrad);
                        const VectorizedArray<Number> jump = array_face[0] - outval;
                        array_f[4][i2] = (avg_grad + jump * tau) * JxW_face * face_quadrature_weight[i2];
                        array_fd[4][i2] = (0.5 * inv_jac[2] * JxW_face) * jump * face_quadrature_weight[i2];
                    }
                    {
                        const VectorizedArray<Number> outval = array_f[5][i2];
                        const VectorizedArray<Number> outgrad = array_fd[5][i2];
                        const VectorizedArray<Number> avg_grad = -0.5 * inv_jac[2] * (array_face[3] + outgrad);
                        const VectorizedArray<Number> jump = array_face[1] - outval;
                        array_f[5][i2] = (avg_grad + jump * tau) * JxW_face * face_quadrature_weight[i2];
                        array_fd[5][i2] = (-0.5 * inv_jac[2] * JxW_face) * jump * face_quadrature_weight[i2];
                    }
                }
            } //end of dim = 3 if statement
            // interpolate external x values for faces
            {
                unsigned int indices[2 * n_lanes];
                for (unsigned int v = 1; v < n_lanes; ++v)
                    indices[v] = ii * dofs_per_cell + v - 1;
                indices[0] = (ii - n_lanes) * dofs_per_cell + n_lanes - 1;
                if (ix == 0) {
                    // assume periodic boundary conditions
                    indices[0] = (ii + (n_cells[0] - 1) * n_lanes) * dofs_per_cell + n_lanes - 1;
                }
                for (unsigned int v = 0; v < n_lanes - 1; ++v)
                    indices[n_lanes + v] = ii * dofs_per_cell + v + 1;
                indices[2 * n_lanes - 1] = (ii + n_lanes) * dofs_per_cell;
                if (ix == n_cells[0] - 1) {
                    // assume periodic boundary conditions
                    indices[2 * n_lanes - 1] = (ii - (n_cells[0] - 1) * n_lanes) + dofs_per_cell;
                }
                for (unsigned int f = 0; f < 2; ++f) {
                    const VectorizedArray<Number> w0 = (f == 0 ? 1. : -1.) * hermite_derivative_on_face;

                    const unsigned int offset1 = (f == 0 ? degree : 0);
                    const unsigned int offset2 = (f == 0 ? degree - 1 : 1);
                    for (unsigned int i = 0; i < dofs_per_face; ++i) {
                        array_f[f][i].gather(sol_old.begin() + (offset1 + i * (degree + 1)) * n_lanes,
                                             indices + f * n_lanes);
                        array_fd[f][i].gather(sol_old.begin() + (offset2 + i * (degree + 1)) * n_lanes,
                                              indices + f * n_lanes);
                        array_fd[f][i] = w0 * (array_fd[f][i] - array_f[f][i]);
                    }

                    // interpolate values onto quadrature points
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 1); ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1 * nn, array_f[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 0); ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1, array_f[f] + i1);

                    // interpolate derivatives onto quadrature points
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 1); ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1 * nn, array_fd[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 0); ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1, array_fd[f] + i1);
                }
            }
            // interpolate external y values for faces
            {
                const unsigned int index[2] = {(iy > 0 ?
                                                (ii - n_cells[0] * n_lanes) :
                                                (ii + (n_cells[1] - 1) * n_cells[0] * n_lanes)
                                               ) * dofs_per_cell,
                                               (iy < n_cells[1] - 1 ?
                                                (ii + n_cells[0] * n_lanes) :
                                                (ii - (n_cells[1] - 1) * n_cells[0] * n_lanes)
                                               ) * dofs_per_cell};
                for (unsigned int f = 2; f < 4; ++f) {
                    const VectorizedArray<Number> w0 = (f == 2 ? 1. : -1.) * hermite_derivative_on_face;

                    for (unsigned int i1 = 0; i1 < (dim > 2 ? (degree + 1) : 1); ++i1) {
                        const unsigned int base_offset1 =
                                i1 * (degree + 1) * (degree + 1) + (f == 2 ? degree : 0) * (degree + 1);
                        const unsigned int base_offset2 =
                                i1 * (degree + 1) * (degree + 1) + (f == 2 ? degree - 1 : 1) * (degree + 1);
                        for (unsigned int i2 = 0; i2 < degree + 1; ++i2) {
                            const unsigned int i = i1 * (degree + 1) + i2;
                            array_f[f][i].load(sol_old.begin() + index[f % 2] + (base_offset1 + i2) * n_lanes);
                            array_fd[f][i].load(sol_old.begin() + index[f % 2] + (base_offset2 + i2) * n_lanes);
                            array_fd[f][i] = w0 * (array_fd[f][i] - array_f[f][i]);
                        }
                    }

                    // interpolate values onto quadrature points
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 1); ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1 * nn, array_f[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 0); ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_f[f] + i1, array_f[f] + i1);

                    // interpolate derivatives onto quadrature points
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 1); ++i1)
                        apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1 * nn, array_fd[f] + i1 * nn);
                    for (unsigned int i1 = 0; i1 < (dim == 3 ? nn : 0); ++i1)
                        apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                                (shape_values_eo, array_fd[f] + i1, array_fd[f] + i1);
                }
            }

            const VectorizedArray<Number> tauy = (degree + 1) * (degree + 1) * penalty_factor * inv_jac[1];
            const VectorizedArray<Number> JxW_facey = my_jxw * inv_jac[1];
            const VectorizedArray<Number> taux = (degree + 1) * (degree + 1) * penalty_factor * inv_jac[0];
            const VectorizedArray<Number> JxW_facex = my_jxw * inv_jac[0];
            for (unsigned int i2 = 0; i2 < (dim == 3 ? degree + 1 : 1); ++i2) {
                const unsigned int offset = i2 * dofs_per_plane;
                VectorizedArray<Number> *array_ptr = data_ptr + offset;
                VectorizedArray<Number> *array_2_ptr = data_ptr + dofs_per_cell + offset;
                const Number *quadrature_ptr = quadrature_weights.begin() + offset;

                VectorizedArray<Number> array_0[dofs_per_plane], array_1[dofs_per_plane];

                for (unsigned int i = 0; i < degree + 1; ++i) {
                    const unsigned int i1 = i2 * (degree + 1) + i;
                    VectorizedArray<Number> array_face[4];
                    apply_1d_matvec_kernel<
                            degree + 1,
                            degree + 1, 1, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
                            (shape_gradients_eo, array_ptr + i, array_1 + i,
                             nullptr, shape_values_on_face_eo.begin(), array_face);

                    // face integrals in y direction
                    const VectorizedArray<Number> weight = make_vectorized_array(face_quadrature_weight[i1]);
                    {
                        const VectorizedArray<Number> outval = array_f[2][i1];
                        const VectorizedArray<Number> outgrad = array_fd[2][i1];
                        const VectorizedArray<Number> avg_grad = 0.5 * inv_jac[1] * (array_face[2] + outgrad);
                        const VectorizedArray<Number> jump = array_face[0] - outval;
                        array_f[2][i1] = (avg_grad + jump * tauy) * JxW_facey * weight;
                        array_fd[2][i1] = (0.5 * inv_jac[1] * JxW_facey) * jump * weight;
                    }
                    {
                        const VectorizedArray<Number> outval = array_f[3][i1];
                        const VectorizedArray<Number> outgrad = array_fd[3][i1];
                        const VectorizedArray<Number> avg_grad = -0.5 * inv_jac[1] * (array_face[3] + outgrad);
                        const VectorizedArray<Number> jump = array_face[1] - outval;
                        array_f[3][i1] = (avg_grad + jump * tauy) * JxW_facey * weight;
                        array_fd[3][i1] = (-0.5 * inv_jac[1] * JxW_facey) * jump * weight;
                    }

                    apply_1d_matvec_kernel<
                            degree + 1, 1, 1, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
                            (shape_gradients_eo, array_ptr + i * (degree + 1), array_0 + i * (degree + 1),
                             nullptr, shape_values_on_face_eo.begin(), array_face);

                    // face integrals in x direction
                    {
                        const VectorizedArray<Number> outval = array_f[0][i1];
                        const VectorizedArray<Number> outgrad = array_fd[0][i1];
                        const VectorizedArray<Number> avg_grad = 0.5 * inv_jac[0] * (array_face[2] + outgrad);
                        const VectorizedArray<Number> jump = array_face[0] - outval;
                        array_f[0][i1] = (avg_grad + jump * taux) * JxW_facex * weight;
                        array_fd[0][i1] = (0.5 * inv_jac[0] * JxW_facex) * jump * weight;
                    }
                    {
                        const VectorizedArray<Number> outval = array_f[1][i1];
                        const VectorizedArray<Number> outgrad = array_fd[1][i1];
                        const VectorizedArray<Number> avg_grad = -0.5 * inv_jac[0] * (array_face[3] + outgrad);
                        const VectorizedArray<Number> jump = array_face[1] - outval;
                        array_f[1][i1] = (avg_grad + jump * taux) * JxW_facex * weight;
                        array_fd[1][i1] = (-0.5 * inv_jac[0] * JxW_facex) * jump * weight;
                    }
                }

                // cell integral on quadrature points
                for (unsigned int q = 0; q < dofs_per_plane; ++q) {
                    array_0[q] *= (inv_jac[0] * inv_jac[0] * my_jxw) * quadrature_ptr[q];
                    array_1[q] *= (inv_jac[1] * inv_jac[1] * my_jxw) * quadrature_ptr[q];
                    if (dim > 2)
                        array_2_ptr[q] *= (inv_jac[2] * inv_jac[2] * my_jxw) * quadrature_ptr[q];
                }
                for (unsigned int i = 0; i < degree + 1; ++i) {
                    const unsigned int i1 = i2 * (degree + 1) + i;
                    VectorizedArray<Number> array_face[4];
                    array_face[0] = array_f[0][i1] + array_f[1][i1];
                    array_face[1] = array_f[0][i1] - array_f[1][i1];
                    array_face[2] = array_fd[0][i1] + array_fd[1][i1];
                    array_face[3] = array_fd[0][i1] - array_fd[1][i1];
#ifdef ONLY_CELL_TERMS
                    apply_1d_matvec_kernel<degree+1,1,1,false,false,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                    apply_1d_matvec_kernel<
                            degree + 1, 1, 1, false, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
#endif
                            (shape_gradients_eo, array_0 + i * (degree + 1), array_ptr + i * (degree + 1),
                             nullptr, shape_values_on_face_eo.begin(), array_face);
                }
                for (unsigned int i = 0; i < degree + 1; ++i) {
                    const unsigned int i1 = i2 * (degree + 1) + i;
                    VectorizedArray<Number> array_face[4];
                    array_face[0] = array_f[2][i1] + array_f[3][i1];
                    array_face[1] = array_f[2][i1] - array_f[3][i1];
                    array_face[2] = array_fd[2][i1] + array_fd[3][i1];
                    array_face[3] = array_fd[2][i1] - array_fd[3][i1];
#ifdef ONLY_CELL_TERMS
                    apply_1d_matvec_kernel<degree+1,degree+1,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                    apply_1d_matvec_kernel<
                            degree + 1,
                            degree + 1, 1, false, true, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
#endif
                            (shape_gradients_eo, array_1 + i, array_ptr + i,
                             array_ptr + i, shape_values_on_face_eo.begin(), array_face);
                }
            }
            if (dim == 3) {
                for (unsigned int i2 = 0; i2 < dofs_per_face; ++i2) {
                    VectorizedArray<Number> array_face[4];
                    array_face[0] = array_f[4][i2] + array_f[5][i2];
                    array_face[1] = array_f[4][i2] - array_f[5][i2];
                    array_face[2] = array_fd[4][i2] + array_fd[5][i2];
                    array_face[3] = array_fd[4][i2] - array_fd[5][i2];
#ifdef ONLY_CELL_TERMS
                    apply_1d_matvec_kernel<degree+1,dofs_per_face,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                    apply_1d_matvec_kernel<degree +
                                           1, dofs_per_face, 1, false, true, VectorizedArray<Number>, VectorizedArray<Number>, false, 2>
#endif
                            (shape_gradients_eo, data_ptr + dofs_per_cell + i2,
                             data_ptr + i2, data_ptr + i2, shape_values_on_face_eo.begin(), array_face);

                    apply_1d_matvec_kernel<degree + 1, dofs_per_face, 0, false, false, VectorizedArray<Number>>
                            (shape_values_eo, data_ptr + i2, data_ptr + i2);
                }
            }

            for (unsigned int i2 = 0; i2 < (dim > 2 ? degree + 1 : 1); ++i2) {
                const unsigned int offset = i2 * dofs_per_plane;
                // y-direction
                for (unsigned int i1 = 0; i1 < nn; ++i1) {
                    apply_1d_matvec_kernel<nn, nn, 0, false, false, VectorizedArray<Number>>
                            (shape_values_eo, data_ptr + offset + i1, data_ptr + offset + i1);
                }
                // x-direction
                VectorizedArray<Number> *__restrict in = data_ptr + i2 * nn * nn;
                for (unsigned int i1 = 0; i1 < nn; ++i1) {
                    apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                            (shape_values_eo, data_ptr + offset + i1 * nn, data_ptr + offset + i1 * nn);

                    if (evaluate_chebyshev) {
                        VectorizedArray<Number> *tmp_array =
                                reinterpret_cast<VectorizedArray<Number> *>
                                (vec_tm.begin() + ii * dofs_per_cell) + offset + i1 * nn;
                        const VectorizedArray<Number> *rhs_array =
                                reinterpret_cast<const VectorizedArray<Number> *>
                                (sol_rhs.begin() + ii * dofs_per_cell) + offset + i1 * nn;
                        const VectorizedArray<Number> *diag_array =
                                reinterpret_cast<const VectorizedArray<Number> *>
                                (mat_diagonal.begin() + ii * dofs_per_cell) + offset + i1 * nn;
                        for (unsigned int i = 0; i < degree + 1; ++i) {
                            const VectorizedArray<Number> res = data_ptr[offset + i1 * nn + i] - rhs_array[i];
                            const VectorizedArray<Number> tmp =
                                    coefficient_tm * tmp_array[i] + coefficient_np * diag_array[i] * res;
                            tmp_array[i] = tmp;
                            (src_array[offset + i1 * nn + i] - tmp).streaming_store(
                                    &dst_array[offset + i1 * nn + i][0]);
                        }
                    } else {
                        for (unsigned int i = 0; i < degree + 1; ++i)
                            data_ptr[offset + i1 * nn + i].streaming_store(&dst_array[offset + i1 * nn + i][0]);
                    }
                }
            }
        }
    }

    template<bool evaluate_chebyshev = true>
    void do_cheby_iter(const AlignedVector<Number> &src,
                       const AlignedVector<Number> &sol_old,
                       const AlignedVector<Number> &mat_diagonal,
                       AlignedVector<Number> &sol_new,
                       AlignedVector<Number> &vec_tm,
                       const Number coefficient_np,
                       const Number coefficient_tm) {
        if (degree < 1)
            return;

#pragma omp parallel
        {
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_START(("real_dg_laplacian_" + std::to_string(dim) +
                                 "d_deg_" + std::to_string(degree) +
                                 (evaluate_chebyshev ? "ch" : "mv")).c_str());
#endif

#pragma omp for schedule (static) collapse(2)
            for (unsigned int ib = 0; ib < n_blocks[2]; ++ib)
                for (unsigned int jb = 0; jb < n_blocks[1]; ++jb)
                    for (unsigned int kb = 0; kb < n_blocks[0]; ++kb)
                        for (unsigned int i = ib * blz; i < std::min(n_cells[2], (ib + 1) * blz); ++i)
                            for (unsigned int j = jb * bly; j < std::min(n_cells[1], (jb + 1) * bly); ++j)
                                do_inner_loop<evaluate_chebyshev>
                                        (kb * blx, std::min(n_cells[0], (kb + 1) * blx), j, i,
                                         src, sol_old, mat_diagonal, sol_new, vec_tm,
                                         coefficient_np, coefficient_tm);

#ifdef LIKWID_PERFMON
            LIKWID_MARKER_STOP(("real_dg_laplacian_" + std::to_string(dim) +
                                "d_deg_" + std::to_string(degree) +
                                 (evaluate_chebyshev ? "ch" : "mv")).c_str());
#endif
        }
    }

    void emulate_cheby_vector_updates() {
        const Number *ptr_old = sol_old.begin();
        Number *ptr_new = sol_new.begin();
        const Number *ptr_diag = mat_diagonal.begin();
        Number *ptr_tmp = sol_tmp.begin();
        const Number *ptr_rhs = sol_rhs.begin();
        const Number coeff1 = 0.6;
        const Number coeff2 = 0.45;

#pragma omp parallel
        {
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_START(("emu_dg_laplacian_" + std::to_string(dim) +
                                 "d_deg_" + std::to_string(degree) + "vu").c_str());
#endif


            std::size_t start_value;
            std::size_t min_value;
            std::size_t end_value;
#pragma omp for schedule (static) collapse(2) private(start_value, min_value, end_value)
            for (unsigned int ib = 0; ib < n_blocks[2]; ++ib)
                for (unsigned int jb = 0; jb < n_blocks[1]; ++jb)
                    for (unsigned int kb = 0; kb < n_blocks[0]; ++kb)
                        for (unsigned int i = ib * blz; i < std::min(n_cells[2], (ib + 1) * blz); ++i)
                            for (unsigned int j = jb * bly; j < std::min(n_cells[1], (jb + 1) * bly); ++j) {
                                const unsigned int ii = (i * n_cells[1] + j) * n_cells[0];

                                start_value = dofs_per_cell * VectorizedArray<Number>::n_array_elements * (kb * blx + ii);
                                min_value = (std::min(n_cells[0], (kb + 1) * blx) + ii);
                                end_value = min_value * dofs_per_cell * VectorizedArray<Number>::n_array_elements;
#pragma omp simd
                                for (std::size_t ix = start_value; ix < end_value; ++ix) {

                                    ptr_tmp[ix] = coeff1 * ptr_tmp[ix] + coeff2 * ptr_diag[ix] *
                                                                         (ptr_new[ix] - ptr_rhs[ix]);
                                    ptr_new[ix] = ptr_old[ix] - ptr_tmp[ix];

                                }
                            }
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_STOP(("emu_dg_laplacian_" + std::to_string(dim) +
                                "d_deg_" + std::to_string(degree) + "vu").c_str());
#endif
        }
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

    void fill_shape_values() {
        constexpr unsigned int n_q_points_1d = degree + 1;
        constexpr unsigned int stride = (n_q_points_1d + 1) / 2;
        shape_values_eo.resize((degree + 1) * stride);
        shape_gradients_eo.resize((degree + 1) * stride);

        HermiteLikePolynomialBasis basis(degree);
        std::vector<double> gauss_points(get_gauss_points(n_q_points_1d));
        for (unsigned int i = 0; i < (degree + 1) / 2; ++i)
            for (unsigned int q = 0; q < stride; ++q) {
                const double p1 = basis.value(i, gauss_points[q]);
                const double p2 = basis.value(i, gauss_points[n_q_points_1d - 1 - q]);
                shape_values_eo[i * stride + q] = 0.5 * (p1 + p2);
                shape_values_eo[(degree - i) * stride + q] = 0.5 * (p1 - p2);
            }
        if (degree % 2 == 0)
            for (unsigned int q = 0; q < stride; ++q)
                shape_values_eo[degree / 2 * stride + q] =
                        basis.value(degree / 2, gauss_points[q]);

        LagrangePolynomialBasis basis_gauss(get_gauss_points(degree + 1));
        for (unsigned int i = 0; i < (degree + 1) / 2; ++i)
            for (unsigned int q = 0; q < stride; ++q) {
                const double p1 = basis_gauss.derivative(i, gauss_points[q]);
                const double p2 = basis_gauss.derivative(i, gauss_points[n_q_points_1d - 1 - q]);
                shape_gradients_eo[i * stride + q] = 0.5 * (p1 + p2);
                shape_gradients_eo[(degree - i) * stride + q] = 0.5 * (p1 - p2);
            }
        if (degree % 2 == 0)
            for (unsigned int q = 0; q < stride; ++q)
                shape_gradients_eo[degree / 2 * stride + q] =
                        basis_gauss.derivative(degree / 2, gauss_points[q]);

        shape_values_on_face_eo.resize(2 * (degree + 1));
        for (unsigned int i = 0; i < degree / 2 + 1; ++i) {
            const double v0 = basis_gauss.value(i, 0);
            const double v1 = basis_gauss.value(i, 1);
            shape_values_on_face_eo[degree - i] = 0.5 * (v0 - v1);
            shape_values_on_face_eo[i] = 0.5 * (v0 + v1);

            const double d0 = basis_gauss.derivative(i, 0);
            const double d1 = basis_gauss.derivative(i, 1);
            shape_values_on_face_eo[degree + 1 + i] = 0.5 * (d0 + d1);
            shape_values_on_face_eo[degree + 1 + degree - i] = 0.5 * (d0 - d1);
        }

        hermite_derivative_on_face = basis.derivative(0, 0);
        if (std::abs(hermite_derivative_on_face[0] + basis.derivative(1, 0)) > 1e-12)
            std::cout << "Error, unexpected value of Hermite shape function derivative: "
                      << hermite_derivative_on_face[0] << " vs "
                      << basis.derivative(1, 0) << std::endl;

        std::vector<double> gauss_weight_1d = get_gauss_weights(n_q_points_1d);
        quadrature_weights.resize(Utilities::pow(n_q_points_1d, dim));
        if (dim == 3)
            for (unsigned int q = 0, z = 0; z < n_q_points_1d; ++z)
                for (unsigned int y = 0; y < n_q_points_1d; ++y)
                    for (unsigned int x = 0; x < n_q_points_1d; ++x, ++q)
                        quadrature_weights[q] = (gauss_weight_1d[z] * gauss_weight_1d[y]) *
                                                gauss_weight_1d[x];
        else if (dim == 2)
            for (unsigned int q = 0, y = 0; y < n_q_points_1d; ++y)
                for (unsigned int x = 0; x < n_q_points_1d; ++x, ++q)
                    quadrature_weights[q] = gauss_weight_1d[y] * gauss_weight_1d[x];
        else if (dim == 1)
            for (unsigned int q = 0; q < n_q_points_1d; ++q)
                quadrature_weights[q] = gauss_weight_1d[q];
        else
            throw;

        face_quadrature_weight.resize(Utilities::pow(n_q_points_1d, dim - 1));
        if (dim == 3)
            for (unsigned int q = 0, y = 0; y < n_q_points_1d; ++y)
                for (unsigned int x = 0; x < n_q_points_1d; ++x, ++q)
                    face_quadrature_weight[q] = gauss_weight_1d[y] * gauss_weight_1d[x];
        else if (dim == 2)
            for (unsigned int q = 0; q < n_q_points_1d; ++q)
                face_quadrature_weight[q] = gauss_weight_1d[q];
        else
            face_quadrature_weight[0] = 1.;
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
