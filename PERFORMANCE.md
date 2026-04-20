
# BLAS vs. OPENBLAS

R distributions can be bundled with different scientific computing libraries. This can drastically impact harmony's performance. Rstudio comes by default with BLAS. In contrast, conda distributions of R are bundled with OPENBLAS. Overall, our benchmarks show that **harmony+OPENBLAS is substantially faster compared harmony+BLAS**. Therefore users with large datasets will benefit using OPENBLAS.

# Using OpenBLAS on Windows

Harmony, was pretuned to just install in windows R installations. By default, R for Windows uses BLAS and LAPACK libraries which can not be utilized by armadillo.

See this [tutorial](https://github.com/david-cortes/R-openblas-in-windows) for instructions on getting R for Windows to use OpenBLAS. 


# Install OpenBLAS in non-conda environment (linux)

Install [ropenblas](https://prdm0.github.io/ropenblas/) from CRAN:
```r
install.packages("ropenblas")
```

## Conda environment

Conda environments ships with openblas so no action is required.

## Validate openblas

Using `sessionInfo()` from R it should report something similar to this: `BLAS/LAPACK: /usr/lib/libopenblas.so.0.3`


# Multithreading

By default harmony turns uses only one core. However, large datasets (>1M cells) may benefit from parallelization. This behavior can be controlled by the `ncores` parameter which expects a number threads which harmony will use for its math operation. Users are advised to increase gradually `ncores` and assess potential performance benefits.

# OpenMP support 

[Armadillo can leverage the OpenMP backend](https://arma.sourceforge.net/faq.html#speed) to provide multithread support for its operations. Overall, performance benefits are negligible for small datasets. For >10M cells, the algorithm can benefit from the OpenMP backend.


## Step 1 - Install OpenMP compatible OpenBLAS

To utilize OpenMP, we first need need to install an OpenBLAS library with OpenMP support, instead of pthreads. This version of OpenBLAS which will allow for proper multithreading control and support when used in conjuction with Armadillo.

### Option 1 - Conda (easy)

By default, OpenBLAS on conda ships with the 'pthreads' backend for multi-threading which is incompatible with OpenMP.  Therefore, users seeking to leverage OpenMP with large datasets, OpenBLAS needs to openmp and not pthread support.

By default environments such as conda will not install openmp supported versions.


```sh
conda activate myenv
conda install -c conda-forge --override-channels r-base blas=*=*openblas* openblas=*=*openmp*
```

### Option 2 - build OpenBLAS from source (hard)

Build openblas from source using the `USE_OPENMP=1` option. After building, you will need to link the library binaries to your R install. You also need to have OpenMP installed in your environment.

## Step 2 - Enable OpenMP support in Harmony

```sh
git clone https://github.com/immunogenomics/harmony/
```

2. Add OpenMP compilation flag

Uncomment `PKG_CXXFLAGS="${PKG_CXXFLAGS} $(SHLIB_OPENMP_CXXFLAGS)"` in the ~src/configure~ file.

This line will include OpenMP during the compilation uses an incorrect version of OpenBLAS you will get a flood of messages.

3. Install harmony in R from the local directory

```r
devtools::install_local("<local-harmoyny-directory>", force=T, upgrade=F)
```

## Step 3 - Ensure that harmony runs properly with OpenMP support

You should *not* get messages as the following:

`OpenBLAS Warning : Detect OpenMP Loop and this application may hang. Please rebuild the library with USE_OPENMP=1 option.`

While the algorithm does not fail, its correctness has not been tested and can not be guaranteed.
