# Process this file with autoconf to produce a configure script.
AC_INIT(relax, 1.1.1, john.shumway@asu.edu)
AM_INIT_AUTOMAKE(relax, 1.1.1)
AC_CONFIG_SRCDIR([src])
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE 

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_AWK
AC_F77_WRAPPERS
AC_LANG(C++)

# Checks for libraries.
ACX_BLAS([
  LIBS="$BLAS_LIBS $LIBS $FLIBS"
  AC_DEFINE(HAVE_BLAS,[],[Define if you have a BLAS library.])])

AM_PATH_XML2(2.4.0,,AC_MSG_ERROR([
*** Libxml is required to build pi; Available at http://www.libxml.org/.]))
LIBS="$XML_LIBS $LIBS"
CPPFLAGS="$XML_CPPFLAGS $CPPFLAGS"
ACX_BLITZ(required)
AC_CHECK_LIB(hdf5,H5Fopen)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.

# Check for float option (saves memory).
AC_ARG_ENABLE(float,[  --enable-float          Uses float instead of double to save memory],
[
AC_DEFINE(ENABLE_FLOAT,[],[Flag to store coordinates as type float.]) 
])


AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT
