
#include <iostream>
#include <iomanip>

#include <omp.h>
#include <chrono>

#include "evaluation_cell_laplacian_vec_ele.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

const unsigned int min_degree = 1;
const unsigned int max_degree = 25;

const std::size_t  vector_size_guess = 200000;
const bool         cartesian         = true;

typedef double value_type;


template <int dim, int degree, typename Number>
void run_program(const unsigned int n_tests)
{
  const unsigned int n_cell_batches = (vector_size_guess+Utilities::pow(degree+1,dim+1)-1) / Utilities::pow(degree+1,dim+1);
  double best_avg = std::numeric_limits<double>::max();
#ifdef _OPENMP
  const unsigned int nthreads = omp_get_max_threads();
#else
  const unsigned int nthreads = 1;
#endif

  {
    EvaluationCellLaplacianVecEle<dim,degree,Number> evaluator;
    evaluator.initialize(n_cell_batches, cartesian);
    evaluator.do_verification();
  }
  for (unsigned int i=0; i<3; ++i)
    {
      std::vector<double> min_time(nthreads), max_time(nthreads), avg_time(nthreads), std_dev(nthreads);
#pragma omp parallel shared(min_time, max_time, avg_time, std_dev)
      {
        EvaluationCellLaplacianVecEle<dim,degree,Number> evaluator;
        evaluator.initialize(n_cell_batches, cartesian);

#ifdef LIKWID_PERFMON
        LIKWID_MARKER_START(("cell_laplacian_deg_" + std::to_string(degree)).c_str());
#endif
        double tmin = 1e10, tmax = 0, tavg = 0, variance = 0;

#pragma omp barrier

#pragma omp for schedule(static)
        for (unsigned int thr=0; thr<nthreads; ++thr)
          for (unsigned int t=0; t<500; ++t)
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
      const std::size_t n_dofs = (std::size_t)n_cell_batches * Utilities::pow(degree+1,dim) * omp_get_max_threads();
      std::cout << "p" << degree << " statistics tall  ";
      const std::size_t ops_interpolate = (/*add*/2*((degree+1)/2)*2 +
                                           /*mult*/degree+1 +
                                           /*fma*/2*((degree-1)*(degree+1)/2));
      const std::size_t ops_approx = (std::size_t)n_cell_batches * omp_get_max_threads()
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
  RunTime<2,min_degree,max_degree,value_type>::run(degree, n_tests);

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif

  return 0;
}
