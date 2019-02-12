# Minigia: Synergia mini apps

* Synergia
	* [https://cdcvs.fnal.gov/redmine/projects/synergia2](https://cdcvs.fnal.gov/redmine/projects/synergia2)

* vectorclass: This distribution includes the vector class library.
	* [http://www.agner.org/optimize/#vectorclass](http://www.agner.org/optimize/#vectorclass)
	* see vectorclass/license.txt

* Eigen: This distribution includes the Eigen class library.
    * [http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/)

## Building
    mkdir build 
    cd build
    CXX=/usr/local/bin/g++-7 cmake .. -DCMAKE_BUILD_TYPE=Release \
        -DDEFINES='-DGSV_AVX -ffast-math'
    make


* Kokkos with CUDA:

    cmake -DCMAKE_BUILD_TYPE=Release -DFFTW_DIR=${FFTW_ROOT} -DKOKKOS_ENABLE_OPENMP=off -DKOKKOS_ENABLE_CUDA=on -DCMAKE_CXX_COMPILER=nvcc_wrapper -DCMAKE_CXX_FLAGS="-arch=sm_61" ..

* Kokkos with OpenMP:

    cmake -DCMAKE_BUILD_TYPE=Release -DFFTW_DIR=${FFTW_ROOT} -DKOKKOS_ENABLE_OPENMP=on ..
