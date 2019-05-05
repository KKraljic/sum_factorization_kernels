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
int max_num_threads = 28;


/// \brief Calculates the total number of cells.
/// Formula: n = vector_size / (degree+1)^dim
/// \param vector_size
/// \param degree
/// \param dim
/// \return 1 or total number of cells in the "grid"
template<int dim, int degree>
const unsigned int get_total_number_of_cells(unsigned long long vector_size){
    //TODO: Unit Tests
    return std::max(vector_size / Utilities::pow(degree + 1, dim),
                    1ULL);
}

template<int dim, typename Number>
void init_number_of_cells_per_dim(unsigned int output_arr[], const unsigned int n_cells_tot){
    //TODO: Unit Tests
    output_arr[0] = std::max(static_cast<unsigned int>(1.00001 * std::pow((double) n_cells_tot, 1. / dim))
                          / VectorizedArray<Number>::n_array_elements,
                          1U) * VectorizedArray<Number>::n_array_elements;
    if (dim > 2) {
        output_arr[1] = output_arr[0];
        output_arr[2] = std::max(n_cells_tot / (output_arr[0] * output_arr[1]), 1U);
    } else {
        output_arr[1] = std::max(n_cells_tot / output_arr[0], 1U);
        output_arr[2] = 1;
    }

}

template<int dim, int degree, typename Number>
void set_block_sizes(EvaluationDGLaplacian<dim, degree, Number> &evaluator){
    evaluator.blx = std::max(20 / (degree + 1), 3);
    evaluator.bly = std::max(2, 20 / (degree + 1));
    evaluator.blz = std::max(2, 20 / (degree + 1));
}

template<int dim, int degree, typename Number>
void adjust_total_number_of_cells_to_block_sizes(unsigned int in_n_cells[], unsigned int out_n_cells[], EvaluationDGLaplacian<dim, degree, Number> &evaluator, const unsigned int n_cells_tot){
    // fit cells to multiple of block size
    if (dim == 3) out_n_cells[2] = std::max(1U, in_n_cells[2] / evaluator.blz) * evaluator.blz;
    out_n_cells[1] = std::max(1U, in_n_cells[1] / evaluator.bly) * evaluator.bly;
    out_n_cells[0] = std::max(n_cells_tot / (out_n_cells[1] * out_n_cells[2]) / VectorizedArray<Number>::n_array_elements, 1U) *
                 VectorizedArray<Number>::n_array_elements;
}

template<int dim, int degree, typename Number>
void run(const unsigned int n_tests) {
        if (degree < 3 || degree > 12){
            std::cout << "Only degree 3-12 implemented, aborting."<< std::endl;
            return;
        }
        for(unsigned long long vector_size = 1; vector_size < pow(2, 63)-1; vector_size*=2){
            for(int num_threads = 1; num_threads <= max_num_threads; num_threads++) {
                omp_set_num_threads(num_threads);
                #pragma omp parallel
                {
                    LIKWID_MARKER_THREADINIT;
                }

                EvaluationDGLaplacian<dim, degree, Number> evaluator;
                const unsigned int n_cells_tot = get_total_number_of_cells<dim, degree>(vector_size);
                unsigned int initial_n_cells[3];
                unsigned int adjusted_n_cells[3];

                init_number_of_cells_per_dim<dim, Number>(initial_n_cells, n_cells_tot);
                set_block_sizes<dim, degree, Number>(evaluator);
                adjust_total_number_of_cells_to_block_sizes<dim, degree, Number>(initial_n_cells, adjusted_n_cells,
                                                                                 evaluator, n_cells_tot);

                evaluator.initialize(adjusted_n_cells);

                sprintf(region_tag, "Inner-Loop-t-%i-n-%u", omp_get_num_threads(), n_cells_tot);
                LIKWID_MARKER_START(region_tag);
                #pragma omp parallel for schedule (static) collapse(2)
                for (unsigned int ib = 0; ib < evaluator.n_blocks[2]; ++ib)
                    for (unsigned int jb = 0; jb < evaluator.n_blocks[1]; ++jb)
                        for (unsigned int kb = 0; kb < evaluator.n_blocks[0]; ++kb)
                            for (unsigned int i = ib * evaluator.blz;
                                 i < std::min(evaluator.n_cells[2], (ib + 1) * evaluator.blz); ++i)
                                for (unsigned int j = jb * evaluator.bly;
                                     j < std::min(evaluator.n_cells[1], (jb + 1) * evaluator.bly); ++j) {
                                    evaluator.do_inner_loop
                                            (kb * evaluator.blx,
                                             std::min(evaluator.n_cells[0], (kb + 1) * evaluator.blx), j, i,
                                             evaluator.sol_new, evaluator.sol_new, evaluator.sol_new, evaluator.sol_new,
                                             evaluator.sol_tmp,
                                             0., 0.);
                                }
                LIKWID_MARKER_STOP(region_tag);
            }
        }
    }

int main(int argc, char **argv) {

    LIKWID_MARKER_INIT;
    unsigned long long int vector_size, n_tests;
    std::string::size_type sz = 0;

    /*if (argc > 1)
        vector_size = std::stoull(argv[1], &sz,0);
    if (argc > 2)
        n_tests = std::stoull(argv[2], &sz,0);*/

    n_tests = 10;

    run<2, 4, double>(n_tests);


    LIKWID_MARKER_CLOSE;

    return 0;
}
