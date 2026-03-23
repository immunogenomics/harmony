# BLAS vs. OPENBLAS

R distributions can be bundled with different scientific computing libraries. This can drastically impact harmony's performance. Rstudio comes by default with BLAS. In contrast, conda distributions of R are bundled with OPENBLAS. Overall, our benchmarks show that **harmony+OPENBLAS is substantially faster compared harmony+BLAS**. Therefore users with large datasets will benefit using OPENBLAS.


## Install OpenBLAS in non-conda environment

Install [ropenblas](https://prdm0.github.io/ropenblas/) from CRAN:
```r
install.packages("ropenblas")
```

## Conda environment

Conda environments come by default with openblas so no action is required.

## Validate openblas
Using `sessionInfo()` from R it should report something similar to this: `BLAS/LAPACK: /usr/lib/libopenblas.so.0.3`


# Multithreading

By default harmony turns uses only one core. However, large datasets (>1M cells) may benefit from parallelization. This behavior can be controlled by the `ncores` parameter which expects a number threads which harmony will use for its math operation. Users are advised to increase gradually `ncores` and assess potential performance benefits.

# OpenMP support

[Armadillo can leverage the OpenMP backend](https://arma.sourceforge.net/faq.html#speed) to provide multithread support for its operations. Overall, performance benefits are relatively small. 

If users want to leverage OpenMP with large datasets, OpenBLAS needs to openmp and not pthread support. There is no reliable way to determine whether the loaded OpenBLAS has support for OpenMP hence this section guiding users how to enable this support.

By default environments such as conda will not install openmp supported versions.

## Building harmony with OpenMP support

0. Ensure your environment supports OpenMP (conda)

### Conda 

Ensure your environment uses OpenMP in the used compilers or install OpenMP `conda install conda-forge::openmp`. 

Then download an (openmp compatible package)[https://anaconda.org/channels/conda-forge/packages/libopenblas/files?name=openmp] for your architecture. For example for linux64 install conda using:

```sh
conda activate myenv
cd <the downloaded conda file>
conda install /home/main/Downloads/libopenblas-0.3.30-openmp_hd680484_4.conda
conda install <package_filename>.conda
```

### Other environments

Build openblas from source using the `USE_OPENMP=1` option.


1. Clone the master branch of the harmony github repository.

```sh
git clone https://github.com/immunogenomics/harmony/
```

2. Add OpenMP compilation flag

Uncomment/Include `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)` in the ~src/Makevars~ file.

This line will include OpenMP during the compilation uses an incorrect version of OpenBLAS you will get a flood of messages.

3. Install harmony in R from the local directory

```r
devtools::install_local("harmony", force=T, upgrade=F)
```


4. Ensure that harmony runs properly

You should *not* get messages as the following:

`OpenBLAS Warning : Detect OpenMP Loop and this application may hang. Please rebuild the library with USE_OPENMP=1 option.`

While the algorithm does not fail, its correctness has not been tested and can not be guaranteed.
