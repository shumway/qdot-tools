# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(geom, 1.1, john.shumway@asu.edu)
AC_CONFIG_SRCDIR([config.h.in])
AM_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_AWK
AC_LANG(C++)

# Checks for libraries.
AM_PATH_XML2(2.4.0,,AC_MSG_ERROR([
*** Libxml is required to build geom; Available at 
http://www.libxml.org/.]))
LIBS="$XML_LIBS $LIBS"
CPPFLAGS="$XML_CFLAGS $CPPFLAGS"

ACX_BLITZ(required)
AC_CHECK_LIB(hdf5,H5Fopen)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile])
AC_OUTPUT
