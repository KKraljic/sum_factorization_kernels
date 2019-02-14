cmake_minimum_required(VERSION 2.8)
project(evaluator CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")
set(CMAKE_BUILD_TYPE "Release")

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message("Using Intel compiler flags")
  set(CMAKE_CXX_FLAGS "-std=c++11 -xhost -qopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "-qopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
else()
  message("Using gcc/clang compiler flags")
  set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -fopenmp -ffp-contract=fast")
  set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
endif()

include_directories(include)
find_package(Likwid)

if (LIKWID_FOUND)
  message("${LIKWID_INCLUDE_PATH}")
  include_directories(${LIKWID_INCLUDE_PATH})

  add_executable(test_program_likwid src/test_program.cc)
  target_link_libraries(test_program_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(test_program_likwid PRIVATE LIKWID_PERFMON)


  add_executable(benchmark_loops_likwid src/benchmark_loops.cc)
  target_link_libraries(benchmark_loops_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(benchmark_loops_likwid PRIVATE LIKWID_PERFMON)

  add_executable(test_dg_precomputed_likwid src/test_dg_precomputed.cc)
  target_link_libraries(test_dg_precomputed_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(test_dg_precomputed_likwid PRIVATE LIKWID_PERFMON)

  add_executable(test_dg_integrate_likwid src/test_dg_integrate.cc)
  target_link_libraries(test_dg_integrate_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(test_dg_integrate_likwid PRIVATE LIKWID_PERFMON)

  add_executable(test_dg_advect_likwid src/test_dg_advect.cc)
  target_link_libraries(test_dg_advect_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(test_dg_advect_likwid PRIVATE LIKWID_PERFMON)

  add_executable(test_program_vect_likwid src/test_program_vect.cc)
  target_link_libraries(test_program_vect_likwid ${LIKWID_LIBRARIES})
  target_compile_definitions(test_program_vect_likwid PRIVATE LIKWID_PERFMON)
endif()
