
#include <iostream>
#include <iomanip>

#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "precomputed_dg_laplacian.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

const unsigned int min_degree = 3;
const unsigned int max_degree = 3;
const unsigned int dimension = 3;
typedef double value_type;
//#define DO_BLOCK_SIZE_TEST

template <int dim, int degree, typename Number>
void run_program(const unsigned int vector_size_guess,
                 const unsigned int n_tests)
{
  int rank = -1;
  int n_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

  PrecomputedDGLaplacian<dim,degree,Number> evaluator;
  const unsigned int n_cells_tot = std::max(vector_size_guess / Utilities::pow(degree+1,dim),
                                            1U);
  unsigned int n_cells[3];
  n_cells[0] = std::max(static_cast<unsigned int>(1.00001*std::pow((double)n_cells_tot, 1./dim))
                        /VectorizedArray<Number>::n_array_elements,
                        1U)*VectorizedArray<Number>::n_array_elements;
  if (dim > 2)
    {
      n_cells[1] = n_cells[0];
      n_cells[2] = std::max(n_cells_tot/(n_cells[0]*n_cells[1]), 1U);
    }
  else
    {
      n_cells[1] = std::max(n_cells_tot/n_cells[0], 1U);
      n_cells[2] = 1;
    }

  evaluator.blx = std::max(20/(degree+1), 3);
  evaluator.bly = std::max(2, 20/(degree+1));
  evaluator.blz = std::max(2, 20/(degree+1));

  // fit cells to multiple of block size
  if (dim == 3)
    n_cells[2] = std::max(1U, n_cells[2] / evaluator.blz) * evaluator.blz;
  n_cells[1] = std::max(1U, n_cells[1] / evaluator.bly) * evaluator.bly;
  n_cells[0] = std::max(n_cells_tot/(n_cells[1] * n_cells[2])/VectorizedArray<Number>::n_array_elements, 1U) * VectorizedArray<Number>::n_array_elements;

  evaluator.initialize(n_cells);

  std::size_t local_size = evaluator.n_elements()*evaluator.dofs_per_cell;
  std::size_t global_size = -1;
  MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (false && rank == 0)
    {
      std::cout << std::endl;
      std::cout << "Polynomial degree: " << degree << std::endl;
      std::cout << "Vector size:       " << global_size << " [";
      for (unsigned int d=0; d<dim; ++d)
        std::cout << n_cells[d] << (d<dim-1 ? " x " : "");
      std::cout << " times " << evaluator.dofs_per_cell << "]" << std::endl;
    }

#ifdef DO_BLOCK_SIZE_TEST
  for (unsigned int i=1; i<8192/evaluator.dofs_per_cell; ++i)
    for (unsigned int j=1; j<40; ++j)
      {
        evaluator.blx = i;
        evaluator.bly = j;
        evaluator.initialize(n_cells);
#endif

  double best_avg = std::numeric_limits<double>::max();

  for (unsigned int i=0; i<3; ++i)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      struct timeval wall_timer;
      gettimeofday(&wall_timer, NULL);
      double start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

      for (unsigned int t=0; t<n_tests; ++t)
        evaluator.matrix_vector_product();

      gettimeofday(&wall_timer, NULL);
      double compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

      double min_time = -1, max_time = -1, avg_time = -1;
      MPI_Allreduce(&compute_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&compute_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&compute_time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      best_avg = std::min(best_avg, avg_time/n_procs);
      if (false && rank == 0)
        {
          std::cout << "Time for operation (min/avg/max): "
                    << min_time/n_tests << " "
                    << avg_time/n_procs/n_tests << " "
                    << max_time/n_tests << " "
                    << std::endl;
        }
    }
  if (rank == 0)
    {
      const std::size_t mem_transfer = global_size * sizeof(Number) *
        2 * n_tests;
      const std::size_t ops_approx = global_size / evaluator.dofs_per_cell
        * ((1<<dim) * (/*add*/2*((degree+1)/2)*2 +
                       /*mult*/degree+1 +
                       /*fma*/2*((degree-1)*(degree+1)/2))*Utilities::pow(degree+1,dim-1)
           + 2*dim * 6 * Utilities::pow(degree+1,dim-1)) * n_tests;
      std::cout << "Degree " << std::setw(2) << degree << "  ";
      for (unsigned int d=0; d<dim; ++d)
        std::cout << n_cells[d] << (d<dim-1 ? " x " : "");
      std::cout << " elem " << evaluator.dofs_per_cell << ", block sizes: "
                << evaluator.blx*VectorizedArray<Number>::n_array_elements
                << " x " << evaluator.bly;
      if (dim==3)
        std::cout << " x " << evaluator.blz;
      std::cout  << ", MDoFs/s: "
                 << global_size * n_tests / best_avg/1e6 << ", GB/s: "
                 << (double)mem_transfer/best_avg*1e-9 << " GFLOP/s: "
                 << (double)ops_approx/best_avg*1e-9
                 << std::endl;
    }
#ifdef DO_BLOCK_SIZE_TEST
      }
#endif
}


template<int dim, int degree, int max_degree, typename Number>
class RunTime
{
public:
  static void run(const unsigned int vector_size_guess,
                  const unsigned int n_tests)
  {
    run_program<dim,degree,Number>(vector_size_guess, n_tests);
    if (degree<max_degree)
      RunTime<dim,(degree<max_degree?degree+1:degree),max_degree,Number>
              ::run(vector_size_guess, n_tests);
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

  MPI_Init(&argc, &argv);

#ifdef _OPENMP
  const unsigned int nthreads = omp_get_max_threads();
#else
  const unsigned int nthreads = 1;
#endif
  std::cout << "Number of threads: " << nthreads << std::endl;
  std::cout << "SIMD width:        "
            << sizeof(VectorizedArray<value_type>)/sizeof(value_type) << std::endl;
  std::size_t  vector_size_guess = 10000000;
  unsigned int n_tests           = 100;
  if (argc > 1)
    vector_size_guess = std::atoi(argv[1]);
  if (argc > 2)
    n_tests = std::atoi(argv[2]);

  RunTime<dimension,min_degree,max_degree,value_type>::run(vector_size_guess, n_tests);
  //run_program<dimension,3,value_type>(vector_size_guess, n_tests);
  //run_program<dimension,6,value_type>(vector_size_guess, n_tests);

  MPI_Finalize();

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif

  return 0;
}
