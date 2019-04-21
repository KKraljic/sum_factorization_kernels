#include "../include/vectorization.h"
#include "../include/aligned_vector.h"
#include "../include/utilities.h"
#include "../include/matrix_vector_kernel.h"


//#include <likwid.h>
char region_tag[80];//80 is an arbitratry number; large enough to keep the tags...


template<size_t number_columns,typename Number>
void run_matvec(){
    //Matrix dimensions: 1xN; M = number of rows, N = entries per row
    //in_vector dimension: Nx1
    //out_vector dimension: 1x1
    AlignedVector<VectorizedArray<Number>> in_vector, in_matrix, out_vector;
    constexpr unsigned int entries_in_row = number_columns;


    //Set sizes of vectors and matrix

    /*
     * =====================================================
     * =====Seg faults right here===========================
     * =====================================================
     *
     * */
    //in_vector.resize_fast(number_columns*2-1);
    in_vector.resize_fast(number_columns);

    /*
   * =====================================================
   * =====End of critical section...======================
   * =====================================================
   *
   * */
    out_vector.resize_fast(number_columns);
    in_matrix.resize_fast(number_columns*number_columns/2);

    //Initialize vectors and matrix with some initial values
    for(size_t i = 0; i < in_vector.size(); ++i) in_vector[i] = 1.0;
    for(size_t i = 0; i < out_vector.size(); ++i) out_vector[i] = 1.0;
    for(size_t i = 0; i < in_matrix.size(); ++i) in_matrix[i] = 2.0;


    //Typecasting in order to use the kernel...
    const VectorizedArray<Number> *__restrict in_vector_ptr = in_vector.begin();
    const VectorizedArray<Number> *in_matrix_ptr = in_matrix.begin();
    VectorizedArray<Number> *out_vector_ptr = out_vector.begin();

    //Do the maths

    //sprintf(region_tag, "1D-MATVEC_Kernel-Degree-%i", degree, degree+2);
    // LIKWID_MARKER_START(region_tag);
    apply_1d_matvec_kernel<entries_in_row, 1, 0, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 0>
                (in_vector_ptr, in_matrix_ptr, out_vector_ptr);

    //LIKWID_MARKER_STOP(region_tag);

    //Print result
    std::cout << "Resulting vector (out_vector): " << std::endl;
    std::cout << out_vector[0][0] << std::endl;
    std::cout << "Finished method" << std::endl;
}


int main(int argc, char **argv) {
     //LIKWID_MARKER_INIT;

    run_matvec<10, double>();

    //run_matvec<59100, double>();

    /*run_matvec<3, double>();
    run_matvec<4, double>();
    run_matvec<5, double>();
    run_matvec<6, double>();
    run_matvec<7, double>();
    run_matvec<8, double>();
    run_matvec<9, double>();
    run_matvec<10, double>();
    run_matvec<11, double>();*/

    //LIKWID_MARKER_CLOSE;

    return 0;
}
