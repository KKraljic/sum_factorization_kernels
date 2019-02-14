#include <iostream>
#include <iomanip>

#include <omp.h>
#include <chrono>

//#define ELE_BASED_LAPLACE
//#define READ_VECTOR

#ifdef ELE_BASED_ADVECT
#define DO_CONVECTION
#include "evaluation_dgele_advect.h"
#elif defined(ELE_BASED_LAPLACE)
#include "evaluation_dgele_laplacian.h"
#else
#include "evaluation_cell_laplacian.h"
#endif

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

const unsigned int min_degree = 1;
const unsigned int max_degree = 25;
const unsigned int dimension  = 3;

#ifdef READ_VECTOR
const std::size_t  vector_size_guess = 50000000;
#else
const std::size_t  vector_size_guess = 50000;
#endif
const bool         cartesian         = true;

typedef double value_type;


template <int dim, int degree, typename Number>
void run_program(const unsigned int n_tests)
{
#ifdef _OPENMP
  const unsigned int nthreads = omp_get_max_threads();
#else
  const unsigned int nthreads = 1;
#endif

  const unsigned int n_cell_batches = (vector_size_guess/
#ifdef READ_VECTOR
                                       nthreads
#else
                                       1
#endif
                                       +Utilities::pow(degree+1,dim+1)-1) / Utilities::pow(degree+1,dim+1);
  double best_avg = std::numeric_limits<double>::max();

  {
    EvaluationCellLaplacian<dim,degree,Number> evaluator;
    evaluator.initialize(n_cell_batches, cartesian);
    evaluator.do_verification();
  }

  const std::size_t n_dofs = (std::size_t)n_cell_batches * Utilities::pow(degree+1,dim) * VectorizedArray<Number>::n_array_elements * omp_get_max_threads();
  const std::size_t ops_interpolate = (/*add*/2*((degree+1)/2)*2 +
                                       /*mult*/degree+1 +
                                       /*fma*/2*((degree-1)*(degree+1)/2));
  const std::size_t ops_approx = (std::size_t)n_cell_batches * VectorizedArray<Number>::n_array_elements * omp_get_max_threads()
    * (
#ifdef DO_CONVECTION
       3 * dim * ops_interpolate * Utilities::pow(degree+1,dim-1)
       + (dim+1) * Utilities::pow(degree+1,dim)
#if defined(DO_FACES)
       + dim * ((dim-1)*2*2*ops_interpolate * Utilities::pow(degree+1,dim-2)
                + (6+2+1)*Utilities::pow(degree+1,dim-1)
                )
#elif defined (ELE_BASED_ADVECT)
       + 2*dim * ((dim-1)*ops_interpolate * Utilities::pow(degree+1,dim-2)
                  + (7)*Utilities::pow(degree+1,dim-1)
                )
       + 2*2* (degree-1 + 2*(degree+1) + 4) * Utilities::pow(degree+1,dim-1)
#endif
#else
       4 * dim * ops_interpolate * Utilities::pow(degree+1,dim-1)
       + dim * 2 * dim * Utilities::pow(degree+1,dim)
#ifdef DO_FACES
       + dim * (3*(dim-1)*2*2*ops_interpolate * Utilities::pow(degree+1,dim-2)
                + (2*2+2*3)*Utilities::pow(degree+1,dim-1)
                + (4*dim-1+2+4+1+2+3*2*dim+3)*Utilities::pow(degree+1,dim-1)
                )
#elif defined (ELE_BASED_LAPLACE)
       + 2*dim * (5*(dim-1)*ops_interpolate * Utilities::pow(degree+1,dim-2)
                  + (4*dim-1+2+2+3+2*dim)*Utilities::pow(degree+1,dim-1)
                )
       + 4*2* (degree+1 + 2*(degree-1) + 2) * Utilities::pow(degree+1,dim-1)
       + 2*(dim-2)* (2*2) * Utilities::pow(degree+1,dim-1)
#endif
#endif
       );

  for (unsigned int i=0; i<3; ++i)
    {
      std::vector<double> min_time(nthreads), max_time(nthreads), avg_time(nthreads), std_dev(nthreads);
#pragma omp parallel shared(min_time, max_time, avg_time, std_dev)
      {
        EvaluationCellLaplacian<dim,degree,Number> evaluator;
        evaluator.initialize(n_cell_batches, cartesian);

#ifdef LIKWID_PERFMON
        LIKWID_MARKER_START(("cell_laplacian_deg_" + std::to_string(degree)).c_str());
#endif
        double tmin = 1e10, tmax = 0, tavg = 0, variance = 0;

#pragma omp barrier

#pragma omp for schedule(static)
        for (unsigned int thr=0; thr<nthreads; ++thr)
#ifdef READ_VECTOR
          for (unsigned int t=0; t<5; ++t)
#else
          for (unsigned int t=0; t<500; ++t)
#endif
            evaluator.matrix_vector_product();

#pragma omp barrier

#pragma omp for schedule(static)
        for (unsigned int thr=0; thr<nthreads; ++thr)
          {
            for (unsigned int t=0; t<n_tests; ++t)
              {
                auto t1 = std::chrono::system_clock::now();
                evaluator.matrix_vector_product();
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
        LIKWID_MARKER_STOP(("cell_laplacian_deg_" + std::to_string(degree)).c_str());
#endif
      }
      double tmin = min_time[0], tmax = max_time[0], tavg = 0, stddev = 0;
      for (unsigned int i=0; i<nthreads; ++i)
        {
          std::cout << "p" << degree << " statistics t" << std::setw(3) << std::right << i << "  ";
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
      std::cout << "p" << degree << " statistics tall  ";
      std::cout << std::setw(12) << tavg
                << "  dev " << std::setw(12) << stddev
                << "  min " << std::setw(12) << tmin
                << "  max " << std::setw(12) << tmax
                << "  DoFs/s " << n_dofs / tavg
                << "  GFlops/s " << 1e-9*ops_approx / tavg
                << std::endl << std::endl;
      best_avg = std::min(tavg, best_avg);
    }
  std::cout << "Best result p" << degree
            << " --- DoFs/s " << std::setw(12) << n_dofs/best_avg
            << "  GFlops/s " << std::setw(12) << 1e-9*ops_approx/best_avg
#ifdef READ_VECTOR
            << "  GB/s " << std::setw(12) << 1e-9*3*sizeof(Number)*n_dofs/best_avg
#endif
            << "  ops/dof  " << (double)ops_approx/n_dofs
            << std::endl;
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

  RunTime<dimension,min_degree,max_degree,value_type>::run(degree, n_tests);

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif

  return 0;
}
