cmake_minimum_required(VERSION 2.8)
project(matvec-measurement CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")
set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

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

  add_executable(matvec-measurement src/main.cc)
  target_link_libraries(matvec-measurement ${LIKWID_LIBRARIES})
  target_compile_definitions(matvec-measurement PRIVATE LIKWID_PERFMON)
  target_compile_options(matvec-measurement PRIVATE -Werror -Wall -Wextra)

endif()
