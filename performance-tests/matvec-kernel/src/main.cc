#include "../include/vectorization.h"
#include "../include/aligned_vector.h"
#include "../include/utilities.h"
#include "../include/matrix_vector_kernel.h"
#include <omp.h>
#include <likwid.h>

#define STRIDE_FACTOR 0
#define AMOUNT_ITERATIONS 100000000
char region_tag[80];//80 is an arbitratry number; large enough to keep the tags...


void setup_result_table(){
    std::cout << "vector_size,striding_factor,amount_flops,amount_stores,GFLOP/s,GSTORE/s,duration_[s],throughput_[e/s],throughput_[GB/s]" << std::endl;
}



template<int number_columns,typename Number>
void run_matvec(int stride_factor = 0){
    //Matrix dimensions: NxN/2; N = entries per row, due do symmetry of matrix divided by 2
    //in_vector dimension: Nx1
    //out_vector dimension: 1xN
    AlignedVector<VectorizedArray<Number>> in_vector, in_matrix, out_vector;
    constexpr unsigned int entries_in_row = number_columns;
    double start_point, end_point, duration;

    //Set sizes of vectors and matrix
    in_vector.resize_fast(number_columns);
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

    reset_flop_ctr();
    //Do the maths
    sprintf(region_tag, "1D-MATVEC_Kernel-Degree-%i-stride-%i", number_columns, stride_factor);

    //Start performance measurements
    start_point = omp_get_wtime();
    LIKWID_MARKER_START(region_tag);
    for (size_t i = 1; i < AMOUNT_ITERATIONS; ++i) {
        apply_1d_matvec_kernel<entries_in_row, 1, 0, true, false, VectorizedArray<Number>, VectorizedArray<Number>, false, 0>
        (in_vector_ptr, in_matrix_ptr, out_vector_ptr);
    }
    LIKWID_MARKER_STOP(region_tag);
    end_point = omp_get_wtime();
    duration = end_point - start_point;

    //End of performance measurements
    // FLOPs
    std::cout << number_columns << ",";
    std::cout << stride_factor << ",";
    print_flop_summary();
    print_gflop_per_second(duration);
    std::cout << duration << std::endl;

    //Elements per second


    //Total number of loaded Doubles

}


int main() {
    LIKWID_MARKER_INIT;

        setup_result_table();
        run_matvec<2, double>(STRIDE_FACTOR);
        run_matvec<3, double>(STRIDE_FACTOR);
        run_matvec<4, double>(STRIDE_FACTOR);
        run_matvec<5, double>(STRIDE_FACTOR);
        run_matvec<6, double>(STRIDE_FACTOR);
        run_matvec<7, double>(STRIDE_FACTOR);
        run_matvec<8, double>(STRIDE_FACTOR);
        run_matvec<9, double>(STRIDE_FACTOR);
        run_matvec<10, double>(STRIDE_FACTOR);
        run_matvec<11, double>(STRIDE_FACTOR);
        run_matvec<12, double>(STRIDE_FACTOR);
        run_matvec<13, double>(STRIDE_FACTOR);
        run_matvec<14, double>(STRIDE_FACTOR);
        run_matvec<15, double>(STRIDE_FACTOR);
        run_matvec<16, double>(STRIDE_FACTOR);
        run_matvec<17, double>(STRIDE_FACTOR);
        run_matvec<18, double>(STRIDE_FACTOR);

















        /*run_matvec<19, double>();
        run_matvec<20, double>();
        run_matvec<21, double>();
        run_matvec<22, double>();
        run_matvec<23, double>();
        run_matvec<24, double>();
        run_matvec<25, double>();
        run_matvec<26, double>();
        run_matvec<27, double>();
        run_matvec<28, double>();
        run_matvec<29, double>();
        run_matvec<30, double>();
        run_matvec<31, double>();
        run_matvec<32, double>();
        run_matvec<33, double>();
        run_matvec<34, double>();
        run_matvec<35, double>();
        run_matvec<36, double>();
        run_matvec<37, double>();
        run_matvec<38, double>();
        run_matvec<39, double>();
        run_matvec<40, double>();
        run_matvec<41, double>();
        run_matvec<42, double>();
        run_matvec<43, double>();
        run_matvec<44, double>();
        run_matvec<45, double>();
        run_matvec<46, double>();
        run_matvec<47, double>();
        run_matvec<48, double>();
        run_matvec<49, double>();
        run_matvec<50, double>();
        run_matvec<51, double>();
        run_matvec<52, double>();
        run_matvec<53, double>();
        run_matvec<54, double>();
        run_matvec<55, double>();
        run_matvec<56, double>();
        run_matvec<57, double>();
        run_matvec<58, double>();
        run_matvec<59, double>();
        run_matvec<60, double>();
        run_matvec<61, double>();
        run_matvec<62, double>();
        run_matvec<63, double>();
        run_matvec<64, double>();
        run_matvec<65, double>();
        run_matvec<66, double>();
        run_matvec<67, double>();
        run_matvec<68, double>();
        run_matvec<69, double>();
        run_matvec<70, double>();
        run_matvec<71, double>();
        run_matvec<72, double>();
        run_matvec<73, double>();
        run_matvec<74, double>();
        run_matvec<75, double>();
        run_matvec<76, double>();
        run_matvec<77, double>();
        run_matvec<78, double>();
        run_matvec<79, double>();
        run_matvec<80, double>();
        run_matvec<81, double>();
        run_matvec<82, double>();
        run_matvec<83, double>();
        run_matvec<84, double>();
        run_matvec<85, double>();
        run_matvec<86, double>();
        run_matvec<87, double>();
        run_matvec<88, double>();
        run_matvec<89, double>();
        run_matvec<90, double>();
        run_matvec<91, double>();
        run_matvec<92, double>();
        run_matvec<93, double>();
        run_matvec<94, double>();
        run_matvec<95, double>();
        run_matvec<96, double>();
        run_matvec<97, double>();
        run_matvec<98, double>();
        run_matvec<99, double>();
        run_matvec<100, double>();
        run_matvec<101, double>();
        run_matvec<102, double>();
        run_matvec<103, double>();
        run_matvec<104, double>();
        run_matvec<105, double>();
        run_matvec<106, double>();
        run_matvec<107, double>();
        run_matvec<108, double>();
        run_matvec<109, double>();
        run_matvec<110, double>();
        run_matvec<111, double>();
        run_matvec<112, double>();
        run_matvec<113, double>();
        run_matvec<114, double>();
        run_matvec<115, double>();
        run_matvec<116, double>();
        run_matvec<117, double>();
        run_matvec<118, double>();
        run_matvec<119, double>();
        run_matvec<120, double>();
        run_matvec<121, double>();
        run_matvec<122, double>();
        run_matvec<123, double>();
        run_matvec<124, double>();
        run_matvec<125, double>();
        run_matvec<126, double>();
        run_matvec<127, double>();
        run_matvec<128, double>();
        run_matvec<129, double>();
        run_matvec<130, double>();
        run_matvec<131, double>();
        run_matvec<132, double>();
        run_matvec<133, double>();
        run_matvec<134, double>();
        run_matvec<135, double>();
        run_matvec<136, double>();
        run_matvec<137, double>();
        run_matvec<138, double>();
        run_matvec<139, double>();
        run_matvec<140, double>();
        run_matvec<141, double>();
        run_matvec<142, double>();
        run_matvec<143, double>();
        run_matvec<144, double>();
        run_matvec<145, double>();
        run_matvec<146, double>();
        run_matvec<147, double>();
        run_matvec<148, double>();
        run_matvec<149, double>();
        run_matvec<150, double>();
        run_matvec<151, double>();
        run_matvec<152, double>();
        run_matvec<153, double>();
        run_matvec<154, double>();
        run_matvec<155, double>();
        run_matvec<156, double>();
        run_matvec<157, double>();
        run_matvec<158, double>();
        run_matvec<159, double>();
        run_matvec<160, double>();
        run_matvec<161, double>();
        run_matvec<162, double>();
        run_matvec<163, double>();
        run_matvec<164, double>();
        run_matvec<165, double>();
        run_matvec<166, double>();
        run_matvec<167, double>();
        run_matvec<168, double>();
        run_matvec<169, double>();
        run_matvec<170, double>();
        run_matvec<171, double>();
        run_matvec<172, double>();
        run_matvec<173, double>();
        run_matvec<174, double>();
        run_matvec<175, double>();
        run_matvec<176, double>();
        run_matvec<177, double>();
        run_matvec<178, double>();
        run_matvec<179, double>();
        run_matvec<180, double>();
        run_matvec<181, double>();
        run_matvec<182, double>();
        run_matvec<183, double>();
        run_matvec<184, double>();
        run_matvec<185, double>();
        run_matvec<186, double>();
        run_matvec<187, double>();
        run_matvec<188, double>();
        run_matvec<189, double>();
        run_matvec<190, double>();
        run_matvec<191, double>();
        run_matvec<192, double>();
        run_matvec<193, double>();
        run_matvec<194, double>();
        run_matvec<195, double>();
        run_matvec<196, double>();
        run_matvec<197, double>();
        run_matvec<198, double>();
        run_matvec<199, double>();
        run_matvec<200, double>();
        run_matvec<201, double>();
        run_matvec<202, double>();
        run_matvec<203, double>();
        run_matvec<204, double>();
        run_matvec<205, double>();
        run_matvec<206, double>();
        run_matvec<207, double>();
        run_matvec<208, double>();
        run_matvec<209, double>();
        run_matvec<210, double>();
        run_matvec<211, double>();
        run_matvec<212, double>();
        run_matvec<213, double>();
        run_matvec<214, double>();
        run_matvec<215, double>();
        run_matvec<216, double>();
        run_matvec<217, double>();
        run_matvec<218, double>();
        run_matvec<219, double>();
        run_matvec<220, double>();
        run_matvec<221, double>();
        run_matvec<222, double>();
        run_matvec<223, double>();
        run_matvec<224, double>();
        run_matvec<225, double>();
        run_matvec<226, double>();
        run_matvec<227, double>();
        run_matvec<228, double>();
        run_matvec<229, double>();
        run_matvec<230, double>();
        run_matvec<231, double>();
        run_matvec<232, double>();
        run_matvec<233, double>();
        run_matvec<234, double>();
        run_matvec<235, double>();
        run_matvec<236, double>();
        run_matvec<237, double>();
        run_matvec<238, double>();
        run_matvec<239, double>();
        run_matvec<240, double>();
        run_matvec<241, double>();
        run_matvec<242, double>();
        run_matvec<243, double>();
        run_matvec<244, double>();
        run_matvec<245, double>();
        run_matvec<246, double>();
        run_matvec<247, double>();
        run_matvec<248, double>();
        run_matvec<249, double>();
        run_matvec<250, double>();
        run_matvec<251, double>();
        run_matvec<252, double>();
        run_matvec<253, double>();
        run_matvec<254, double>();
        run_matvec<255, double>();
        run_matvec<256, double>();
         */

    LIKWID_MARKER_CLOSE;

    return 0;
}
