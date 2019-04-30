#   -*- mode: Autoconf; tab-width: 4; -*-
# Process this file with autoconf to produce a configure script.      

# require at least autoconf 2.61
AC_PREREQ(2.61)

# Every other copy of the package version number gets its value from here
# The m4_esyscmd_s invocation is due to an excellent StackOverflow answer.
#AC_INIT([harmonyFast], m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))
AC_INIT([harmonyFast], [0.1.0])

PKG_CXXFLAGS="-DUSE_FLOAT_TYPES=0"
AC_ARG_WITH([openblas],
    [AC_HELP_STRING([--with-openblas],
        [Use high performance OpenBLAS library with Harmony])],
              [PKG_CXXFLAGS="-DUSE_FLOAT_TYPES=1"],
)

AC_SUBST(PKG_CXXFLAGS)
AC_OUTPUT([src/Makevars])
