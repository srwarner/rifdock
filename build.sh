#!/bin/bash

#rm -r build/*
export CC=/home/bcov/linuxbrew/bin/gcc
export CXX=/home/bcov/linuxbrew/bin/g++

cd 112919_build
CMAKE_ROSETTA_PATH=/home/bcov/rifdock/rosetta_builds/bcov_fordas_structure_store2_grid_silent_19-04-03/main CMAKE_FINAL_ROSETTA_PATH=/home/bcov/rifdock/rosetta_builds/bcov_fordas_structure_store2_grid_silent_19-04-03/main/source/cmake/build_cxx11_omp_hdf5 cmake .. -DCMAKE_BUILD_TYPE=Release -DUSEHDF5=1 -DUSEGRIDSCORE=1
make -j20 rifgen rif_dock_test

