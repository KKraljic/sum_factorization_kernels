#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include "evaluation_dg_laplacian.h"

#include <climits>
#include "../include/utilities.h"
#include <likwid.h>

//#include <omp.h>


const unsigned int min_degree = 3;
const unsigned int max_degree = 12;
const unsigned int dimension = 3;
typedef double value_type;
char region_tag[100];



template<int dim, int degree, typename Number>
void run(const unsigned int n_tests) {
        if (degree < 3){
            std::cout << "Only rank 3 and higher implemented, aborting."<< std::endl;
            return;
        }

        unsigned long long vector_size = 100000;

        EvaluationDGLaplacian<dim, degree, Number> evaluator;
        const unsigned int n_cells_tot = std::max(vector_size / Utilities::pow(degree + 1, dim),
                                                  1ULL);
        unsigned int n_cells[3];
        n_cells[0] = std::max(static_cast<unsigned int>(1.00001 * Utilities::pow((double) n_cells_tot, 1. / dim))
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

        //const unsigned int start_x,
        //                       const unsigned int end_x,
        //                       const unsigned int iy,
        //                       const unsigned int iz,
        //                       const AlignedVector<Number> &src,
        //                       const AlignedVector<Number> &sol_old,
        //                       const AlignedVector<Number> &mat_diagonal,
        //                       AlignedVector<Number> &sol_new,
        //                       AlignedVector<Number> &vec_tm,
        //                       const Number coefficient_np,
        //                       const Number coefficient_tm

        /*
         *  sol_old.resize_fast(n_elements() * dofs_per_cell);
        sol_new.resize_fast(n_elements() * dofs_per_cell);
        mat_diagonal.resize_fast(n_elements() * dofs_per_cell);
          static constexpr unsigned int dofs_per_cell = Utilities::pow(degree + 1, dim);
          sol_tmp.resize_fast(n_elements() * dofs_per_cell);
         * */

        for (unsigned int ib = 0; ib < evaluator.n_blocks[2]; ++ib)
            for (unsigned int jb = 0; jb < evaluator.n_blocks[1]; ++jb)
                for (unsigned int kb = 0; kb < evaluator.n_blocks[0]; ++kb)
                    for (unsigned int i = ib * evaluator.blz; i < std::min(n_cells[2], (ib + 1) * evaluator.blz); ++i)
                        for (unsigned int j = jb * evaluator.bly; j < std::min(n_cells[1], (jb + 1) * evaluator.bly); ++j)
                            evaluator.do_inner_loop
                                    (kb * evaluator.blx, std::min(n_cells[0], (kb + 1) * evaluator.blx), j, i,
                                     evaluator.sol_old, evaluator.sol_old, evaluator.mat_diagonal, evaluator.sol_new, evaluator.sol_tmp,
                                     0., 0.);
    }

int main(int argc, char **argv) {

    //LIKWID_MARKER_INIT;
/*#pragma omp parallel
    {
      LIKWID_MARKER_THREADINIT;
    }
#endif*/
    unsigned long long int vector_size, n_tests;
    std::string::size_type sz = 0;

    /*if (argc > 1)
        vector_size = std::stoull(argv[1], &sz,0);
    if (argc > 2)
        n_tests = std::stoull(argv[2], &sz,0);*/

    vector_size = 100000;
    n_tests = 10;

    run<2, 3, double>(n_tests);
    run<2, 4, double>(n_tests);
    run<2, 5, double>(n_tests);
    run<2, 6, double>(n_tests);
    run<2, 7, double>(n_tests);
    run<2, 8, double>(n_tests);
    run<2, 9, double>(n_tests);
    run<2, 10, double>(n_tests);
    run<2, 11, double>(n_tests);
    run<2, 12, double>(n_tests);

    run<3, 3, double>(n_tests);
    run<3, 4, double>(n_tests);
    run<3, 5, double>(n_tests);
    run<3, 6, double>(n_tests);
    run<3, 7, double>(n_tests);
    run<3, 8, double>(n_tests);
    run<3, 9, double>(n_tests);
    run<3, 10, double>(n_tests);
    run<3, 11, double>(n_tests);
    run<3, 12, double>(n_tests);

    //LIKWID_MARKER_CLOSE;

    return 0;
}
