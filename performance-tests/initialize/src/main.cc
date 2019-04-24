
#include <iostream>
#include <iomanip>

#include <sys/time.h>
#include <sys/resource.h>

#include "evaluation_dg_laplacian.h"


#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#ifdef _OPENMP

#include <omp.h>

#endif

const unsigned int min_degree = 3;
const unsigned int max_degree = 12;
const unsigned int dimension = 3;
typedef double value_type;
char region_tag[100];
//#define DO_BLOCK_SIZE_TEST

template<int dim, int degree, typename Number>
void run_program(const unsigned int vector_size_guess,
                 const unsigned int n_tests,
                 const unsigned int variants) {
    // currently only degrees 3 and higher implemented
    if (degree < 3) return;

    int rank = 0, n_procs = 1;

    EvaluationDGLaplacian<dim, degree, Number> evaluator;
    const unsigned int n_cells_tot = std::max(vector_size_guess / Utilities::pow(degree + 1, dim),
                                              1U);
    unsigned int n_cells[3];
    n_cells[0] = std::max(static_cast<unsigned int>(1.00001 * std::pow((double) n_cells_tot, 1. / dim))
                          / VectorizedArray<Number>::n_array_elements,
                          1U) * VectorizedArray<Number>::n_array_elements;
    if (dim > 2) {
        n_cells[1] = n_cells[0];
        n_cells[2] = std::max(n_cells_tot / (n_cells[0] * n_cells[1]), 1U);
    } else {
        n_cells[1] = std::max(n_cells_tot / n_cells[0], 1U);
        n_cells[2] = 1;
    }

    evaluator.blx = std::max(20 / (degree + 1), 3);
    evaluator.bly = std::max(2, 20 / (degree + 1));
    evaluator.blz = std::max(2, 20 / (degree + 1));

    // fit cells to multiple of block size
    if (dim == 3) n_cells[2] = std::max(1U, n_cells[2] / evaluator.blz) * evaluator.blz;
    n_cells[1] = std::max(1U, n_cells[1] / evaluator.bly) * evaluator.bly;
    n_cells[0] = std::max(n_cells_tot / (n_cells[1] * n_cells[2]) / VectorizedArray<Number>::n_array_elements, 1U) *
                 VectorizedArray<Number>::n_array_elements;

    evaluator.initialize(n_cells);

    std::size_t local_size = evaluator.n_elements() * evaluator.dofs_per_cell;
    std::size_t global_size = local_size;


#ifdef DO_BLOCK_SIZE_TEST
    for (unsigned int i=1; i<8192/evaluator.dofs_per_cell; ++i)
      for (unsigned int j=1; j<40; ++j)
        {
          evaluator.blx = i;
          evaluator.bly = j;
          evaluator.initialize(n_cells);
#endif

    for (unsigned int i = 0; i < 5; ++i) {
        for (unsigned int t = 0; t < n_tests; ++t) evaluator.do_matvec();
    }

#ifdef DO_BLOCK_SIZE_TEST
    }
#endif

    if (variants > 1) {
        best_avg = std::numeric_limits<double>::max();

        for (unsigned int i = 0; i < 5; ++i) {
            for (unsigned int t = 0; t < n_tests; ++t)
                evaluator.do_chebyshev();
        }
    }

    if (variants > 2) {
        for (unsigned int i = 0; i < 5; ++i) {
            for (unsigned int t = 0; t < n_tests; ++t)
                evaluator.emulate_cheby_vector_updates();
        }

    }
}


template<int dim, int degree, int max_degree, typename Number>
class RunTime {
public:
    static void run(const int target_degree,
                    const unsigned int vector_size_guess,
                    const unsigned int n_tests,
                    const unsigned int variants) {
        if (degree == target_degree || target_degree == -1)
            run_program<dim, degree, Number>(vector_size_guess, n_tests, variants);
        if (degree < max_degree)
            RunTime<dim, (degree < max_degree ? degree + 1 : degree), max_degree, Number>
            ::run(target_degree, vector_size_guess, n_tests, variants);
    }
};

int main(int argc, char **argv) {
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
#pragma omp parallel
    {
      LIKWID_MARKER_THREADINIT;
    }
#endif

    //MPI_Init(&argc, &argv);

#ifdef _OPENMP
    const unsigned int nthreads = omp_get_max_threads();
#else
    const unsigned int nthreads = 1;
#endif
    std::cout << "Number of threads: " << nthreads << std::endl;
    std::cout << "SIMD width:        "
              << sizeof(VectorizedArray<value_type>) / sizeof(value_type) << std::endl;
    int degree = 4;
    std::size_t vector_size_guess = 10000000;
    unsigned int n_tests = 100;
    unsigned int variants = 2;
    if (argc > 1)
        degree = std::atoi(argv[1]);
    if (argc > 2)
        vector_size_guess = std::atoi(argv[2]);
    if (argc > 3)
        n_tests = std::atoi(argv[3]);
    if (argc > 4)
        variants = std::atoi(argv[4]);

    RunTime<dimension, min_degree, max_degree, value_type>::run(degree, vector_size_guess,
                                                                n_tests, variants);

    //MPI_Finalize();

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif

    return 0;
}
