// This file is free software. You can use it, redistribute it, and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.
//
// implementation of cell and face terms for advection using integration on
// Cartesian cell geometries with integration
//
// Author: Martin Kronbichler, April 2018

#ifndef evaluation_dg_advect_h
#define evaluation_dg_advect_h


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


template <int dim>
struct Point
{
  double data[dim];
};

template <int dim>
class ExactSolution
{
public:
  ExactSolution (const double time = 0.)
    :
    time(time),
    wave_number(4.)
  {
    advection.data[0] = 1.;
    if (dim>1)
      advection.data[1] = 0.15;
    if (dim>2)
      advection.data[2] = -0.05;
  }

  double value (const Point<dim> &p) const
  {
    const double PI = 3.14159265358979323846;
    double position[dim];
    for (unsigned int d=0; d<dim; ++d)
      position[d] = p.data[d] - time*advection.data[d];
    double result = std::sin(wave_number*position[0]*PI);
    for (unsigned int d=1; d<dim; ++d)
      result *= std::cos(wave_number*position[d]*PI);
    return result;
  }

  double time_derivative (const Point<dim> &p) const
  {
    const double PI = 3.14159265358979323846;
    double position[dim];
    for (unsigned int d=0; d<dim; ++d)
      position[d] = p.data[d] - time*advection.data[d];
    double result = -advection.data[0] * wave_number*PI*std::cos(wave_number*position[0]*PI);
    for (unsigned int d=1; d<dim; ++d)
      result *= std::cos(wave_number*position[d]*PI);
    double add = wave_number*PI*std::sin(wave_number*PI*position[0]);
    if (dim == 2)
      add *= advection.data[1] * std::sin(wave_number*PI*position[1]);
    else if (dim==3)
      add *= (advection.data[1] * std::sin(wave_number*PI*position[1])*std::cos(wave_number*position[2]*PI)
              +advection.data[2] * std::cos(wave_number*PI*position[1])*std::sin(wave_number*position[2]*PI));
    return result+add;
  }

  Point<dim> get_transport_direction() const
  {
    return advection;
  }

private:
  Point<dim> advection;
  const double time;
  const double wave_number;
};



template <int dim, int degree, typename Number>
class EvaluationDGAdvection
{
public:
  static constexpr unsigned int dimension = dim;
  static constexpr unsigned int n_q_points = Utilities::pow(degree+1,dim);
  static constexpr unsigned int dofs_per_cell = Utilities::pow(degree+1,dimension);
  unsigned int blx;
  unsigned int bly;
  unsigned int blz;

  void initialize(const unsigned int *n_cells_in)
  {
    n_cells[0] = n_cells_in[0]/VectorizedArray<Number>::n_array_elements;
    for (unsigned int d=1; d<dim; ++d)
      n_cells[d] = n_cells_in[d];
    for (unsigned int d=dim; d<3; ++d)
      n_cells[d] = 1;

    n_blocks[2] = (n_cells[2] + blz - 1)/blz;
    n_blocks[1] = (n_cells[1] + bly - 1)/bly;
    n_blocks[0] = (n_cells[0] + blx - 1)/blx;

    sol_old.resize(0);
    sol_new.resize(0);
    sol_tmp.resize(0);
    sol_old.resize_fast(n_elements() * dofs_per_cell);
    sol_new.resize_fast(n_elements() * dofs_per_cell);
    sol_tmp.resize_fast(n_elements() * dofs_per_cell);

    std::vector<double> jacobian  = get_diagonal_jacobian();
    std::vector<double> gl_points = get_gauss_lobatto_points(degree+1);

    ExactSolution<dim> exact(0.);

#pragma omp parallel
    {
#pragma omp for schedule (static) collapse(2)
      for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
        for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
          for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
            for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
              for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
                {
                  const unsigned int ii=(i*n_cells[1]+j)*n_cells[0];
                  for (std::size_t ix=kb*blx; ix<std::min(n_cells[0], (kb+1)*blx); ++ix)
                    for (unsigned int lz=0,c=0; lz<(dim>2 ? degree+1 : 1); ++lz)
                      for (unsigned int ly=0; ly<degree+1; ++ly)
                        for (unsigned int lx=0; lx<degree+1; ++lx)
                          for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v, ++c)
                            {
                              Point<dim> p;
                              p.data[0] = (ix*VectorizedArray<Number>::n_array_elements+v + gl_points[lx]) / jacobian[0];
                              p.data[1] = (j + gl_points[ly]) / jacobian[1];
                              if (dim>2)
                                p.data[2] = (i + gl_points[lz]) / jacobian[2];
                              const std::size_t idx = (ix+ii) * dofs_per_cell * VectorizedArray<Number>::n_array_elements + c;
                              sol_new[idx] = exact.value(p);
                              sol_old[idx] = 0.;
                              sol_tmp[idx] = 0.;
                            }
                }
    }

    Point<dim> transport = exact.get_transport_direction();
    double trans_norm = 0;
    for (unsigned int d=0; d<dim; ++d)
      trans_norm += transport.data[d] * transport.data[d];
    trans_norm = std::sqrt(trans_norm);
    time_step = 0.2 / jacobian[0] / trans_norm / degree;


    Number jacobian_determinant = 1.;
    for (unsigned int d=0; d<dim; ++d)
      jacobian_determinant *= jacobian[d];
    jacobian_determinant = 1./jacobian_determinant;

    jxw_data.resize(1);
    jxw_data[0] = jacobian_determinant;
    jacobian_data.resize(dim);
    for (unsigned int d=0; d<dim; ++d)
      jacobian_data[d] = jacobian[d];

    fill_shape_values();
  }

  std::size_t n_elements() const
  {
    std::size_t n_element = VectorizedArray<Number>::n_array_elements;
    for (unsigned int d=0; d<dim; ++d)
      n_element *= n_cells[d];
    return n_element;
  }

  void do_time_step(double &current_time)
  {
    const double a21 = 970286171893./4311952581923.;
    const double a32 = 6584761158862./12103376702013.;
    const double a43 = 2251764453980./15575788980749.;
    const double a54 = 26877169314380./34165994151039.;

    const double b1 =  1153189308089./22510343858157.;
    const double b2 = 1772645290293./ 4653164025191.;
    const double b3 = -1672844663538./ 4480602732383.;
    const double b4 =  2114624349019./3568978502595.;
    const double b5 =  5198255086312./14908931495163.;

    do_rk_stage(sol_new, sol_new, sol_old, sol_tmp, current_time,
                a21*time_step, (b1-a21)*time_step);
    do_rk_stage(sol_old, sol_tmp, sol_new, sol_tmp, current_time+a21*time_step,
                a32*time_step, (b2-a32)*time_step);
    do_rk_stage(sol_new, sol_tmp, sol_old, sol_tmp, current_time+(b1+a32)*time_step,
                a43*time_step, (b3-a43)*time_step);
    do_rk_stage(sol_old, sol_tmp, sol_tmp, sol_new, current_time+(b1+b2+a43)*time_step,
                a54*time_step, (b4-a54)*time_step);
    do_rk_stage(sol_tmp, sol_new, sol_new, sol_tmp, current_time+(b1+b2+b3+a54)*time_step,
                b5*time_step, 0.0);
    current_time += time_step;
  }

  void do_matvec()
  {
    do_rk_stage<false>(sol_new, sol_new, sol_tmp, sol_tmp, 0., 0., 0.);
  }

  void verify_derivative()
  {
    do_rk_stage<false>(sol_new, sol_new, sol_tmp, sol_tmp, 0., 0., 0.);
    std::vector<double> gl_points = get_gauss_lobatto_points(degree+1);
    double error = 0.;
    ExactSolution<dim> exact(0.);
    for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
      for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
        for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
          for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
            for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
              {
                const unsigned int ii=(i*n_cells[1]+j)*n_cells[0];
                for (std::size_t ix=kb*blx; ix<std::min(n_cells[0], (kb+1)*blx); ++ix)
                  for (unsigned int lz=0,c=0; lz<(dim>2 ? degree+1 : 1); ++lz)
                    for (unsigned int ly=0; ly<degree+1; ++ly)
                      for (unsigned int lx=0; lx<degree+1; ++lx)
                        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v, ++c)
                          {
                            Point<dim> p;
                            p.data[0] = (ix*VectorizedArray<Number>::n_array_elements+v + gl_points[lx]) / jacobian_data[0][0];
                            p.data[1] = (j + gl_points[ly]) / jacobian_data[1][0];
                            if (dim>2)
                              p.data[2] = (i + gl_points[lz]) / jacobian_data[2][0];
                            const std::size_t idx = (ix+ii) * dofs_per_cell * VectorizedArray<Number>::n_array_elements + c;
                            error = std::max (std::abs(exact.time_derivative(p) - sol_tmp[idx]),
                                              error);
                          }
              }
    std::cout << "Error against analytic time derivative: " << error << std::endl;
  }

  void compute_max_error(const double time)
  {
    std::vector<double> gl_points = get_gauss_lobatto_points(degree+1);
    double error = 0.;
    ExactSolution<dim> exact(time);
    for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
      for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
        for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
          for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
            for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
              {
                const unsigned int ii=(i*n_cells[1]+j)*n_cells[0];
                for (std::size_t ix=kb*blx; ix<std::min(n_cells[0], (kb+1)*blx); ++ix)
                  for (unsigned int lz=0,c=0; lz<(dim>2 ? degree+1 : 1); ++lz)
                    for (unsigned int ly=0; ly<degree+1; ++ly)
                      for (unsigned int lx=0; lx<degree+1; ++lx)
                        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v, ++c)
                          {
                            Point<dim> p;
                            p.data[0] = (ix*VectorizedArray<Number>::n_array_elements+v + gl_points[lx]) / jacobian_data[0][0];
                            p.data[1] = (j + gl_points[ly]) / jacobian_data[1][0];
                            if (dim>2)
                              p.data[2] = (i + gl_points[lz]) / jacobian_data[2][0];
                            const std::size_t idx = (ix+ii) * dofs_per_cell * VectorizedArray<Number>::n_array_elements + c;
                            error = std::max (std::abs(exact.value(p) - sol_new[idx]),
                                              error);
                          }
              }
    std::cout << "Error against analytic solution at t=" << time
              << ": " << error << std::endl;
  }

  template <bool evaluate_stage=true>
  void do_inner_loop (const unsigned int start_x,
                      const unsigned int end_x,
                      const unsigned int iy,
                      const unsigned int iz,
                      const AlignedVector<Number> &src,
                      const AlignedVector<Number> &update_vec,
                      AlignedVector<Number>       &vec_np,
                      AlignedVector<Number>       &vec_tm,
                      const double                 current_time,
                      const Number                 coefficient_np,
                      const Number                 coefficient_tm)
  {
    ExactSolution<dim> solution(current_time);
    Point<dim> advection = solution.get_transport_direction();

    constexpr unsigned int nn = degree+1;
    constexpr unsigned int mid = nn/2;
    constexpr unsigned int n_lanes = VectorizedArray<Number>::n_array_elements;
    constexpr unsigned int dofs_per_face = Utilities::pow(degree+1,dim-1);
    constexpr unsigned int dofs_per_plane = Utilities::pow(degree+1,2);
    const VectorizedArray<Number> *__restrict shape_values_eo = this->shape_values_eo.begin();
    const VectorizedArray<Number> *__restrict inv_shape_values_eo = this->inv_shape_values_eo.begin();
    const VectorizedArray<Number> *__restrict shape_gradients_eo = this->shape_gradients_eo.begin();
    AlignedVector<VectorizedArray<Number> > scratch_data_array;
    VectorizedArray<Number> my_array[degree < 13 ? 2*dofs_per_cell : 1];
    VectorizedArray<Number> *__restrict data_ptr;
    VectorizedArray<Number> array_f[6][dofs_per_face];
    if (degree < 13)
      data_ptr = my_array;
    else
      {
        scratch_data_array.resize_fast(2*dofs_per_cell);
        data_ptr = scratch_data_array.begin();
      }

    for (unsigned int ix=start_x; ix<end_x; ++ix)
      {
        const unsigned int ii=((iz*n_cells[1]+iy)*n_cells[0]+ix)*n_lanes;
        const VectorizedArray<Number>* src_array =
          reinterpret_cast<const VectorizedArray<Number>*>(src.begin()+ii*dofs_per_cell);
        VectorizedArray<Number>* dst_array =
          reinterpret_cast<VectorizedArray<Number>*>(vec_np.begin()+ii*dofs_per_cell);

        const VectorizedArray<Number> * inv_jac = jacobian_data.begin();
        const VectorizedArray<Number> my_jxw = jxw_data[0];

        for (unsigned int i2=0; i2<(dim>2 ? degree+1 : 1); ++i2)
          {
            // x-direction
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, src_array+i2*nn*nn+i1*nn, in+i1*nn);
              }
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, in+i1, in+i1);
              }
          }

        if (dim == 3)
          {
            const unsigned int index[2] = {(iz > 0 ?
                                            (ii-n_cells[1]*n_cells[0]*n_lanes) :
                                            (ii+(n_cells[2]-1)*n_cells[1]*n_cells[0]*n_lanes)
                                            )*dofs_per_cell,
                                           (iz < n_cells[2]-1 ?
                                            (ii+n_cells[1]*n_cells[0]*n_lanes) :
                                            (ii-(n_cells[2]-1)*n_cells[1]*n_cells[0]*n_lanes)
                                            )*dofs_per_cell};

            for (unsigned int f=4; f<6; ++f)
              {
                const unsigned int offset1 = (f==4 ? dofs_per_face*degree : 0);
                for (unsigned int i=0; i<dofs_per_face; ++i)
                  array_f[f][i].load(src.begin()+index[f%2]+(offset1+i)*n_lanes);

                // interpolate values onto quadrature points
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
                for (unsigned int i1=0; i1<nn; ++i1)
                  apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                    (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
              }

            const VectorizedArray<Number> JxW_face = my_jxw * inv_jac[2];
            for (unsigned int i2=0; i2<dofs_per_face; ++i2)
              {
                apply_1d_matvec_kernel<degree+1,dofs_per_face,0,true,false,VectorizedArray<Number>>
                  (shape_values_eo, data_ptr+i2, data_ptr+i2);

                // evaluation from cell onto face
                VectorizedArray<Number> r0 = shape_values_on_face_eo[0] * (data_ptr[i2]+
                                                                           data_ptr[i2+nn*nn*(nn-1)]);
                VectorizedArray<Number> r1 = shape_values_on_face_eo[nn-1] * (data_ptr[i2]-
                                                                              data_ptr[i2+nn*nn*(nn-1)]);
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_values_on_face_eo[ind] * (data_ptr[i2+nn*nn*ind]+
                                                          data_ptr[i2+nn*nn*(nn-1-ind)]);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (data_ptr[i2+nn*nn*ind]-
                                                               data_ptr[i2+nn*nn*(nn-1-ind)]);
                  }
                if (nn%2 == 1)
                  r0 += shape_values_on_face_eo[mid] * data_ptr[i2+nn*nn*mid];

                // face integrals in z direction
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(-advection.data[2]);

                  const VectorizedArray<Number> u_minus = r0 + r1;
                  const VectorizedArray<Number> u_plus  = array_f[4][i2];
                  array_f[4][i2] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_face * face_quadrature_weight[i2];
                }
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(advection.data[2]);

                  const VectorizedArray<Number> u_minus = r0 - r1;
                  const VectorizedArray<Number> u_plus  = array_f[5][i2];
                  array_f[5][i2] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_face * face_quadrature_weight[i2];
                }
              }
          }

        // interpolate external x values for faces
        {
          unsigned int indices[2*n_lanes];
          for (unsigned int v=1; v<n_lanes; ++v)
            indices[v] = ii*dofs_per_cell+v-1;
          indices[0] = (ii-n_lanes)*dofs_per_cell+n_lanes-1;
          if (ix==0)
            {
              // assume periodic boundary conditions
              indices[0] = (ii+(n_cells[0]-1)*n_lanes)*dofs_per_cell+n_lanes-1;
            }
          for (unsigned int v=0; v<n_lanes-1; ++v)
          indices[n_lanes+v] = ii*dofs_per_cell+v+1;
          indices[2*n_lanes-1] = (ii+n_lanes)*dofs_per_cell;
          if (ix==n_cells[0]-1)
            {
              // assume periodic boundary conditions
              indices[2*n_lanes-1] = (ii-(n_cells[0]-1)*n_lanes)+dofs_per_cell;
            }
          for (unsigned int f=0; f<2; ++f)
            {
              const unsigned int offset1 = (f==0 ? degree : 0);
              for (unsigned int i=0; i<dofs_per_face; ++i)
                array_f[f][i].gather(src.begin()+(offset1+i*(degree+1))*n_lanes, indices+f*n_lanes);

              // interpolate values onto quadrature points
              for (unsigned int i1=0; i1<(dim==3 ? nn : 1); ++i1)
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
              for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
            }
        }
        // interpolate external y values for faces
        {
          const unsigned int index[2] = {(iy > 0 ?
                                          (ii-n_cells[0]*n_lanes) :
                                          (ii+(n_cells[1]-1)*n_cells[0]*n_lanes)
                                          ) * dofs_per_cell,
                                         (iy < n_cells[1]-1 ?
                                          (ii+n_cells[0]*n_lanes) :
                                          (ii-(n_cells[1]-1)*n_cells[0]*n_lanes)
                                          ) * dofs_per_cell};
          for (unsigned int f=2; f<4; ++f)
            {
              for (unsigned int i1=0; i1<(dim>2 ? (degree+1) : 1); ++i1)
                {
                  const unsigned int base_offset1 = i1*(degree+1)*(degree+1)+(f==2 ? degree : 0)*(degree+1);
                  for (unsigned int i2=0; i2<degree+1; ++i2)
                    array_f[f][i1*(degree+1)+i2].load(src.begin()+index[f%2]+(base_offset1+i2)*n_lanes);
                }

              // interpolate values onto quadrature points
              for (unsigned int i1=0; i1<(dim==3 ? nn : 1); ++i1)
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1*nn, array_f[f]+i1*nn);
              for (unsigned int i1=0; i1<(dim==3 ? nn : 0); ++i1)
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (shape_values_eo, array_f[f]+i1, array_f[f]+i1);
            }
        }

        const VectorizedArray<Number> JxW_facey = my_jxw * inv_jac[1];
        const VectorizedArray<Number> JxW_facex = my_jxw * inv_jac[0];
        for (unsigned int i2=0; i2<(dim==3 ? degree+1 : 1); ++i2)
          {
            const unsigned int offset = i2*dofs_per_plane;
            VectorizedArray<Number> *array_ptr = data_ptr + offset;
            VectorizedArray<Number> *array_2_ptr = data_ptr + dofs_per_cell + offset;
            const Number *quadrature_ptr = quadrature_weights.begin() + offset;

            VectorizedArray<Number> array_0[dofs_per_plane], array_1[dofs_per_plane];

            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                const VectorizedArray<Number> weight = make_vectorized_array(face_quadrature_weight[i1]);
                // evaluation from cell onto face
                VectorizedArray<Number> r0 = shape_values_on_face_eo[0] * (array_ptr[i]+
                                                                           array_ptr[i+nn*(nn-1)]);
                VectorizedArray<Number> r1 = shape_values_on_face_eo[nn-1] * (array_ptr[i]-
                                                                              array_ptr[i+nn*(nn-1)]);
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_values_on_face_eo[ind] * (array_ptr[i+nn*ind]+
                                                          array_ptr[i+nn*(nn-1-ind)]);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (array_ptr[i+nn*ind]-
                                                               array_ptr[i+nn*(nn-1-ind)]);
                  }
                if (nn%2 == 1)
                  r0 += shape_values_on_face_eo[mid] * array_ptr[i+nn*mid];

                // face integrals in y direction
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(-advection.data[1]);

                  const VectorizedArray<Number> u_minus = r0 + r1;
                  const VectorizedArray<Number> u_plus  = array_f[2][i1];
                  array_f[2][i1] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_facey * weight;
                }
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(advection.data[1]);

                  const VectorizedArray<Number> u_minus = r0 - r1;
                  const VectorizedArray<Number> u_plus  = array_f[3][i1];
                  array_f[3][i1] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_facey * weight;
                }

                // evaluation from cell onto face
                r0 = shape_values_on_face_eo[0] * (array_ptr[i*nn]+
                                                   array_ptr[i*nn+(nn-1)]);
                r1 = shape_values_on_face_eo[nn-1] * (array_ptr[i*nn]-
                                                      array_ptr[i*nn+(nn-1)]);
                for (unsigned int ind=1; ind<mid; ++ind)
                  {
                    r0 += shape_values_on_face_eo[ind] * (array_ptr[i*nn+ind]+
                                                          array_ptr[i*nn+(nn-1-ind)]);
                    r1 += shape_values_on_face_eo[nn-1-ind] * (array_ptr[i*nn+ind]-
                                                               array_ptr[i*nn+(nn-1-ind)]);
                  }
                if (nn%2 == 1)
                  r0 += shape_values_on_face_eo[mid] * array_ptr[i*nn+mid];

                // face integrals in x direction
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(-advection.data[0]);

                  const VectorizedArray<Number> u_minus = r0 + r1;
                  const VectorizedArray<Number> u_plus  = array_f[0][i1];
                  array_f[0][i1] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_facex * weight;
                }
                {
                  const VectorizedArray<Number> normal_times_advection =
                    make_vectorized_array<Number>(advection.data[0]);

                  const VectorizedArray<Number> u_minus = r0 - r1;
                  const VectorizedArray<Number> u_plus  = array_f[1][i1];
                  array_f[1][i1] = -0.5 * ((u_minus+u_plus) * normal_times_advection +
                                           std::abs(normal_times_advection) * (u_minus-u_plus)) * JxW_facex * weight;
                }
              }

            // cell integral on quadrature points
            for (unsigned int q=0; q<dofs_per_plane; ++q)
              {
                array_0[q] = advection.data[0] * inv_jac[0]*(my_jxw*quadrature_ptr[q] * array_ptr[q]);
                array_1[q] = advection.data[1] * inv_jac[1]*(my_jxw*quadrature_ptr[q] * array_ptr[q]);
                if (dim>2)
                  array_2_ptr[q] = advection.data[2] * inv_jac[2]*(my_jxw*quadrature_ptr[q] * array_ptr[q]);
              }
            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                VectorizedArray<Number> array_face[2];
                array_face[0] = array_f[0][i1]+array_f[1][i1];
                array_face[1] = array_f[0][i1]-array_f[1][i1];
#ifdef ONLY_CELL_TERMS
                apply_1d_matvec_kernel<degree+1,1,1,false,false,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                apply_1d_matvec_kernel<degree+1,1,1,false,false,VectorizedArray<Number>,VectorizedArray<Number>,false,1>
#endif
                  (shape_gradients_eo, array_0+i*(degree+1), array_ptr+i*(degree+1),
                   nullptr, shape_values_on_face_eo.begin(), array_face);
              }
            for (unsigned int i=0; i<degree+1; ++i)
              {
                const unsigned int i1 = i2*(degree+1)+i;
                VectorizedArray<Number> array_face[2];
                array_face[0] = array_f[2][i1]+array_f[3][i1];
                array_face[1] = array_f[2][i1]-array_f[3][i1];
#ifdef ONLY_CELL_TERMS
                apply_1d_matvec_kernel<degree+1,degree+1,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                apply_1d_matvec_kernel<degree+1,degree+1,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,1>
#endif
                  (shape_gradients_eo, array_1+i, array_ptr+i,
                   array_ptr+i, shape_values_on_face_eo.begin(), array_face);
              }
          }
        if (dim == 3)
          {
            for (unsigned int i2=0; i2<dofs_per_face; ++i2)
              {
                VectorizedArray<Number> array_face[2];
                array_face[0] = array_f[4][i2]+array_f[5][i2];
                array_face[1] = array_f[4][i2]-array_f[5][i2];
#ifdef ONLY_CELL_TERMS
                apply_1d_matvec_kernel<degree+1,dofs_per_face,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,0>
#else
                apply_1d_matvec_kernel<degree+1,dofs_per_face,1,false,true,VectorizedArray<Number>,VectorizedArray<Number>,false,1>
#endif
                  (shape_gradients_eo, data_ptr+dofs_per_cell+i2,
                   data_ptr+i2, data_ptr+i2, shape_values_on_face_eo.begin(), array_face);

                // apply inverse mass matrix
                for (unsigned int i=0; i<degree+1; ++i)
                  data_ptr[i2+i*nn*nn] /= my_jxw*quadrature_weights[i2+i*nn*nn];

                apply_1d_matvec_kernel<degree+1,dofs_per_face,0,true,false,VectorizedArray<Number>>
                  (inv_shape_values_eo, data_ptr+i2, data_ptr+i2);
              }
          }
        else
          {
            // apply inverse mass matrix
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              data_ptr[i] /= my_jxw*quadrature_weights[i];
          }

        for (unsigned int i2=0; i2< (dim>2 ? degree+1 : 1); ++i2)
          {
            const unsigned int offset = i2*dofs_per_plane;
            // y-direction
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, nn, 0, true, false, VectorizedArray<Number>>
                  (inv_shape_values_eo, data_ptr+offset+i1, data_ptr+offset+i1);
              }
            // x-direction
            VectorizedArray<Number> *__restrict in = data_ptr + i2*nn*nn;
            for (unsigned int i1=0; i1<nn; ++i1)
              {
                apply_1d_matvec_kernel<nn, 1, 0, true, false, VectorizedArray<Number>>
                  (inv_shape_values_eo, data_ptr+offset+i1*nn, data_ptr+offset+i1*nn);


                if (evaluate_stage)
                  {
                    VectorizedArray<Number>* tmp_array =
                      reinterpret_cast<VectorizedArray<Number>*>
                      (vec_tm.begin()+ii*dofs_per_cell) + offset + i1*nn;
                    const VectorizedArray<Number>* upd_array =
                      reinterpret_cast<const VectorizedArray<Number>*>
                      (update_vec.begin()+ii*dofs_per_cell) + offset + i1*nn;
                    for (unsigned int i=0; i<degree+1; ++i)
                      {
                        const VectorizedArray<Number> vecnp = upd_array[i] + coefficient_np*data_ptr[offset+i1*nn+i];
                        vecnp.streaming_store(&dst_array[offset+i1*nn+i][0]);
                        if (coefficient_tm != 0.)
                          (vecnp + coefficient_tm*data_ptr[offset+i1*nn+i]).
                            streaming_store(&tmp_array[i][0]);
                      }
                  }
                else
                  {
                    for (unsigned int i=0; i<degree+1; ++i)
                      data_ptr[offset+i1*nn+i].streaming_store(&dst_array[offset+i1*nn+i][0]);
                  }
              }
          }
      }
  }

  template <bool evaluate_stage=true>
  void do_rk_stage(const AlignedVector<Number> &src,
                   const AlignedVector<Number> &update_vec,
                   AlignedVector<Number>       &vec_np,
                   AlignedVector<Number>       &vec_tm,
                   const double                 current_time,
                   const Number                 coefficient_np,
                   const Number                 coefficient_tm)
  {
    if (degree < 1)
      return;

#pragma omp parallel
    {
#ifdef LIKWID_PERFMON
      LIKWID_MARKER_START(("dg_advect_" + std::to_string(dim) +
                           "d_deg_" + std::to_string(degree) +
                           (evaluate_stage ? "rk" : "mv")).c_str());
#endif

#pragma omp for schedule (static) collapse(2)
      for (unsigned int ib=0; ib<n_blocks[2]; ++ib)
        for (unsigned int jb=0; jb<n_blocks[1]; ++jb)
          for (unsigned int kb=0; kb<n_blocks[0]; ++kb)
            for (unsigned int i=ib*blz; i<std::min(n_cells[2], (ib+1)*blz); ++i)
              for (unsigned int j=jb*bly; j<std::min(n_cells[1], (jb+1)*bly); ++j)
                do_inner_loop<evaluate_stage>
                  (kb*blx, std::min(n_cells[0], (kb+1)*blx), j, i,
                   src, update_vec, vec_np, vec_tm, current_time,
                   coefficient_np, coefficient_tm);

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_STOP(("dg_advect_" + std::to_string(dim) +
                          "d_deg_" + std::to_string(degree) +
                          (evaluate_stage ? "rk" : "mv")).c_str());
#endif
    }
  }

private:

  std::vector<double> get_diagonal_jacobian() const
  {
    const double domain_size = 2.;
    std::vector<double> jacobian(dim);
    for (unsigned int d=0; d<dim; ++d)
      jacobian[d] = n_cells[d]/domain_size;
    jacobian[0] *= VectorizedArray<Number>::n_array_elements;
    return jacobian;
  }

  void fill_shape_values()
  {
    constexpr unsigned int n_q_points_1d = degree+1;
    constexpr unsigned int stride = (n_q_points_1d+1)/2;
    inv_shape_values_eo.resize((degree+1)*stride);
    shape_values_eo.resize((degree+1)*stride);
    shape_gradients_eo.resize((degree+1)*stride);

    std::vector<double> points = get_gauss_lobatto_points(n_q_points_1d);
    std::vector<double> gauss_points = get_gauss_points(n_q_points_1d);
    std::vector<double> gauss_weights = get_gauss_weights(n_q_points_1d);
    LagrangePolynomialBasis basis(points);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis.value(i, gauss_points[q]);
          const double p2 = basis.value(i, gauss_points[n_q_points_1d-1-q]);
          shape_values_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_values_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_values_eo[degree/2*stride+q] =
          basis.value(degree/2, gauss_points[q]);

    LagrangePolynomialBasis basis_gauss(gauss_points);
    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          const double p1 = basis_gauss.derivative(i, gauss_points[q]);
          const double p2 = basis_gauss.derivative(i, gauss_points[n_q_points_1d-1-q]);
          shape_gradients_eo[i*stride+q] = 0.5 * (p1 + p2);
          shape_gradients_eo[(degree-i)*stride+q] = 0.5 * (p1 - p2);
        }
    if (degree%2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        shape_gradients_eo[degree/2*stride+q] =
          basis_gauss.derivative(degree/2, gauss_points[q]);

    shape_values_on_face_eo.resize(2*(degree+1));
    for (unsigned int i=0; i<degree/2+1; ++i)
      {
        const double v0 = basis_gauss.value(i, 0);
        const double v1 = basis_gauss.value(i, 1);
        shape_values_on_face_eo[degree-i] = 0.5 * (v0 - v1);
        shape_values_on_face_eo[i] = 0.5 * (v0 + v1);

        const double d0 = basis_gauss.derivative(i, 0);
        const double d1 = basis_gauss.derivative(i, 1);
        shape_values_on_face_eo[degree+1+i] = 0.5 * (d0 + d1);
        shape_values_on_face_eo[degree+1+degree-i] = 0.5 * (d0 - d1);
      }

    for (unsigned int i=0; i<(degree+1)/2; ++i)
      for (unsigned int q=0; q<stride; ++q)
        {
          inv_shape_values_eo[i*stride+q] =
            0.5 * (basis_gauss.value(i, points[q]) +
                   basis_gauss.value(i, points[degree-q]));
          inv_shape_values_eo[(degree-i)*stride+q] =
            0.5 * (basis_gauss.value(i, points[q]) -
                   basis_gauss.value(i, points[degree-q]));
        }
    if (degree % 2 == 0)
      for (unsigned int q=0; q<stride; ++q)
        inv_shape_values_eo[degree/2*stride+q] =
          basis_gauss.value(degree/2, points[q]);

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

  unsigned int n_cells[3];
  unsigned int n_blocks[3];
  AlignedVector<VectorizedArray<Number> > shape_values_eo;
  AlignedVector<VectorizedArray<Number> > shape_gradients_eo;
  AlignedVector<VectorizedArray<Number> > shape_values_on_face_eo;
  AlignedVector<VectorizedArray<Number> > inv_shape_values_eo;

  AlignedVector<Number> quadrature_weights;
  AlignedVector<Number> face_quadrature_weight;

  AlignedVector<Number> sol_old;
  AlignedVector<Number> sol_new;
  AlignedVector<Number> sol_tmp;

  AlignedVector<VectorizedArray<Number> > jxw_data;
  AlignedVector<VectorizedArray<Number> > jacobian_data;

  double time_step;
};


#endif
