dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl     $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id: acinclude.m4,v 1.2 2005/01/07 18:35:55 jshumwa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
        [AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
        yes | "") ;;
        no) acx_blas_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
        *) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
        AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
        AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
        AC_MSG_RESULT($acx_blas_ok)
        LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS"
        AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
        LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(atlas, ATL_xerbla,
                [AC_CHECK_LIB(f77blas, $sgemm,
                [AC_CHECK_LIB(cblas, cblas_dgemm,
                        [acx_blas_ok=yes
                         BLAS_LIBS="-lcblas -lf77blas -latlas"],
                        [], [-lf77blas -latlas])],
                        [], [-latlas])])
fi

# BLAS in Apple vecLib framework? (Mac OS X)
if test $acx_blas_ok = no; then
        vlib_flags="-framework vecLib"
        save_LIBS="$LIBS"; LIBS="$vlib_flags $LIBS"
        AC_MSG_CHECKING([for $sgemm in $vlib_flags])
        AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes; BLAS_LIBS="$vlib_flags"],
                         [BLAS_LIBS=""])
        AC_MSG_RESULT($acx_blas_ok)
        LIBS="$save_LIBS"
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm,
                [AC_CHECK_LIB(dgemm, $dgemm,
                [AC_CHECK_LIB(sgemm, $sgemm,
                        [acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
                        [], [-lblas])],
                        [], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
        if test "x$GCC" != xyes; then # only works with Sun CC
                AC_CHECK_LIB(sunmath, acosp,
                        [AC_CHECK_LIB(sunperf, $sgemm,
                                [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
        fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(complib.sgimath, $sgemm,
                     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm,
                [AC_CHECK_LIB(essl, $sgemm,
                        [acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
                        [], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS
dnl @synopsis ACX_BLITZ([optional-string "required"])
dnl
dnl Check whether Blitz++ is installed.
dnl Blitz++ is available at http://oonumerics.org/blitz.
dnl
dnl   Set the path for Blitz++  with the option
dnl      --with-blitz[=DIR]
dnl   Blitz headers should be under DIR/includes
dnl   Blitz library should be under DIR/lib
dnl   Then try to compile and run a simple program with a Blitz Array
dnl   Optional argument `required' triggers an error if Blitz++ not installed
dnl  dnl @version $Id: acinclude.m4,v 1.2 2005/01/07 18:35:55 jshumwa Exp $
dnl @author Patrick Guio <address@bogus.example.com>
dnl
AC_DEFUN([AC_MSG_ERROR_BLITZ],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Blitz++ template library
available at http://oonumerics.org/blitz
When installed give the directory of installation with the option
  --with-blitz@<:@=DIR@:>@
])])


AC_DEFUN([ACX_BLITZ],[
 
AC_ARG_WITH(blitz,
AC_HELP_STRING([--with-blitz@<:@=DIR@:>@],[Set the path for Blitz++]),
[],[withval='yes'])
 
if test "$1" = required -a "$withval" = no ; then
        AC_MSG_ERROR_BLITZ
fi

if test "$withval" != no ; then
 
        saveCPPFLAGS=$CPPFLAGS
        saveLDFLAGS=$LDFLAGS
        saveLIBS=$LIBS
 
        if test "$withval" != 'yes'; then
                CPPFLAGS="-I$withval/include"
                LDFLAGS="-L$withval/lib"
        fi
        LIBS="-lblitz $LIBS"

        AC_CACHE_CHECK([whether Blitz++ is installed],acx_blitz,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include <blitz/array.h>
]],[[
blitz::Array<int,1> x(10);
x = blitz::tensor::i;
        ]])],[acx_blitz=yes],[acx_blitz=no])
        AC_LANG_RESTORE
        ])

        CPPFLAGS=$saveCPPFLAGS
        LDFLAGS=$saveLDFLAGS
        LIBS=$saveLIBS

 
        if test "$acx_blitz" = yes ; then
                if test "$withval" != yes ; then
                        CPPFLAGS="-I$withval/include $CPPFLAGS"
                        LDFLAGS="-L$withval/lib $LDFLAGS"
                fi
                LIBS="-lblitz $LIBS"
        else
                if test "$1" = required ; then
                        AC_MSG_ERROR_BLITZ
                fi
        fi

 
fi

 
])
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
AC_MSG_CHECKING(for lapack libraries)
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
      
        acx_lapack_ok=no
        $2
fi
AC_MSG_RESULT($acx_lapack_ok)
])dnl ACX_LAPACK
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_mpi.html
dnl
AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc, $CC)
	acx_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpiCC mpCC, mpic++, $CXX)
	acx_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r, $F77)
	acx_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_TRY_LINK([],[      call MPI_Init], [MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
	[C++], [CXX="$acx_mpi_save_CXX"],
	[Fortran 77], [F77="$acx_mpi_save_F77"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI
dnl @synopsis ACX_TVMET([optional-string "required"])
dnl
dnl Check whether tvmet is installed.
dnl tvmet is available at http://tvmet.sourceforge.net/.
dnl
dnl   Set the path for TVMET  with the option
dnl      --with-tvmet[=DIR]
dnl   tvmet headers should be under DIR/includes
dnl   tvmet library should be under DIR/lib
dnl   Then try to compile and run a simple program with a tvmet Array
dnl   Optional argument `required' triggers an error if tvmet not installed
dnl  dnl @version $Id: acinclude.m4,v 1.2 2005/01/07 18:35:55 jshumwa Exp $
dnl @author Patrick Guio <address@bogus.example.com> 
dnl (modified by John Shumway <john.shumway@asu.edu> for TVMET)
dnl
AC_DEFUN([AC_MSG_ERROR_TVMET],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the tvmet template library
available at http://tvmet.sourceforge.net
When installed give the directory of installation with the option
  --with-tvmet@<:@=DIR@:>@
])])


AC_DEFUN([ACX_TVMET],[
 
AC_ARG_WITH(tvmet,
AC_HELP_STRING([--with-tvmet@<:@=DIR@:>@],[Set the path for tvmet]),
[],[withval='yes'])
 
if test "$1" = required -a "$withval" = no ; then
        AC_MSG_ERROR_TVMET
fi

if test "$withval" != no ; then
 
        saveCPPFLAGS=$CPPFLAGS
        saveLDFLAGS=$LDFLAGS
 
        if test "$withval" != 'yes'; then
                CPPFLAGS="-I$withval/include"
                LDFLAGS="-L$withval/lib"
        fi

        AC_CACHE_CHECK([whether tvmet is installed],acx_tvmet,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
]],[[
tvmet::Vector<double,3> x(1.0);
tvmet::Matrix<double,3,3> a(1.0);
        ]])],[acx_tvmet=yes],[acx_tvmet=no])
        AC_LANG_RESTORE
        ])

        CPPFLAGS=$saveCPPFLAGS
        LDFLAGS=$saveLDFLAGS

 
        if test "$acx_tvmet" = yes ; then
                if test "$withval" != yes ; then
                        CPPFLAGS="-I$withval/include $CPPFLAGS"
                fi
        else
                if test "$1" = required ; then
                        AC_MSG_ERROR_TVMET
                fi
        fi

fi
 
])
