#include "../include/vectorization.h"
#include "../include/aligned_vector.h"
#include "../include/utilities.h"
#include "../include/matrix_vector_kernel.h"


#include <likwid.h>
char region_tag[80];//80 is an arbitratry number; large enough to keep the tags...


template<int degree, typename Number>
void run_matvec(){
    AlignedVector<VectorizedArray<Number>> vector, matrix_in_array, matrix_out_array;

    //Initialize the vector with some initial values
    vector.resize_fast(degree+1);
    for(ssize_t i = 0; i < degree + 1; ++i) vector[i] = 1.0;

    //Initialize the input and output array with some initial values
    matrix_in_array.resize_fast((degree+1) * (degree+1));
    matrix_out_array.resize_fast((degree+1) * (degree+1));
    for(ssize_t i = 0; i < (degree+1) * (degree+1); ++i){
        matrix_in_array[i] = 2.0;
        matrix_out_array[i] = 1.0;
    }

    //Typecasting in order to use the kernel...
    const VectorizedArray<Number> *__restrict vector_begin = vector.begin();
    const VectorizedArray<Number> *in_matrix = matrix_in_array.begin();
    VectorizedArray<Number> *out_matrix = matrix_out_array.begin();
    constexpr unsigned int entries_per_row = degree + 1;

    //Do the maths

    sprintf(region_tag, "1D-MATVEC_Kernel-%i", degree);
    LIKWID_MARKER_START(region_tag);
    for (unsigned int i = 0; i < degree; ++i) {
        apply_1d_matvec_kernel<entries_per_row, 1, 0, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 0>
                (vector_begin, in_matrix + i*entries_per_row, out_matrix + i*entries_per_row);
    }
    LIKWID_MARKER_STOP(region_tag);

    //Print result
    /*std::cout << "Resulting array (matrix_out_array): " << std::endl;
    for(ssize_t j = 0; j < (degree+1); ++j) {
        for (ssize_t i = 0; i < (degree + 1); ++i) {
            std::cout << matrix_out_array[i][0] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;*/
}


int main(int argc, char **argv) {
    LIKWID_MARKER_INIT;

    run_matvec<2, double>();
    run_matvec<3, double>();
    run_matvec<4, double>();
    run_matvec<5, double>();
    run_matvec<6, double>();
    run_matvec<7, double>();
    run_matvec<8, double>();
    run_matvec<9, double>();
    run_matvec<10, double>();
    run_matvec<11, double>();

    LIKWID_MARKER_CLOSE;

    return 0;
}
