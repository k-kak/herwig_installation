dnl --- check for Herwig 7 --
AC_DEFUN([HJETS_CHECK_HERWIG],
[
defaultlocation="${prefix}"
test "x$defaultlocation" = xNONE && defaultlocation="${ac_default_prefix}"
AC_MSG_CHECKING([for Herwig 7 in])
AC_ARG_WITH(herwig,
        AC_HELP_STRING([--with-herwig=DIR],[location of Herwig 7 installation]),
        [],
	[with_herwig="${defaultlocation}"])
AC_MSG_RESULT([$with_herwig])

AS_IF([test "x$with_herwig" != "xno"],
      [AC_CHECK_FILES(
      ${with_herwig}/bin/herwig-config,
      [have_herwig=yes], [have_herwig=no])],
      [have_herwig=no])

AS_IF([test "x$have_herwig" = "xyes"],
      [HERWIGCPPFLAGS=`${with_herwig}/bin/herwig-config --cppflags`
       HERWIGLDFLAGS=`${with_herwig}/bin/herwig-config --ldflags`
       HERWIGLDLIBS=`${with_herwig}/bin/herwig-config --ldlibs`
      ],
      [AS_IF([test "x$with_herwig" != "xno"],
             [AC_MSG_ERROR([Cannot build HJets without Herwig. Please set --with-herwig.])
      ])
      ])

AM_CPPFLAGS="-I\$(top_builddir)/include $HERWIGCPPFLAGS"
AM_LDFLAGS="$HERWIGLDFLAGS $HERWIGLDLIBS"

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

])

dnl --- set fortran flags ---
AC_DEFUN([HJETS_FORTRAN_FLAGS],
[

AM_FCFLAGS="$AM_FCFLAGS -fno-automatic -ffixed-line-length-none"
AM_FFLAGS="$AM_FFLAGS -fno-automatic -ffixed-line-length-none"

AC_SUBST(AM_FCFLAGS)
AC_SUBST(AM_FFLAGS)

])

