#! /bin/bash
module load likwid/4.3
rm -rf build/
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=icpc  ../
make -j 32
cd ..

