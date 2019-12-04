#!/bin/bash

#rm -r build/*
cd build_pm
CMAKE_ROSETTA_PATH=/home/bcov/rifdock/rosetta_builds/bcov_fordas_structure_store2_grid_silent_19-04-03/main CMAKE_FINAL_ROSETTA_PATH=/home/bcov/rifdock/rosetta_builds/bcov_fordas_structure_store2_grid_silent_19-04-03/main/source/cmake/build_cxx11_omp_hdf5_debug cmake .. -DCMAKE_BUILD_TYPE=Debug -DUSEHDF5=1
make -j20 scorePM

