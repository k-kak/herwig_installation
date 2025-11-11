dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([HJets],[1.3],[simon.plaetzer@gmail.com],[HJets])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([Amplitudes/AmplitudeBase.h])
AC_CANONICAL_HOST

if test -n "$F77" ; then
  AC_MSG_ERROR([F77 is set internally. Please use FC to specify a non-default Fortran compiler.])
fi

AC_LANG([C++])

AM_INIT_AUTOMAKE([1.9 gnu subdir-objects dist-bzip2 -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
dnl Checks for programs.

AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])

AC_PROG_FC([gfortran])
AC_SUBST([F77],[$FC])

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_PREREQ([2.2.6])
LT_INIT([disable-static dlopen pic-only])

HJETS_CHECK_HERWIG

HJETS_FORTRAN_FLAGS

AC_CONFIG_FILES([Amplitudes/Makefile
		 Amplitudes/AmplitudeBase.cc
		 Virtuals/Amplitudes/Makefile
		 Virtuals/Integrals/Makefile
		 Virtuals/MatrixElement/Makefile
		 Virtuals/TenRed/Makefile
		 Virtuals/Utilities/Makefile
		 Virtuals/OneLOop/Makefile
		 Virtuals/Makefile
      		 include/Makefile
      		 src/Makefile
		 Makefile])

AC_OUTPUT
