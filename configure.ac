AC_INIT([riotng],[0.1],[kiel@inf.uni-due.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_PROG_AR
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
	src/Makefile
	libtm.pc
])
AC_CONFIG_MACRO_DIR([m4])

dnl Compiler
USER_CXXFLAGS=$CXXFLAGS
USER_CPPFLAGS=$CPPFLAGS
AC_PROG_CXX
CXXFLAGS="$USER_CXXFLAGS $CXX0X_FLAGS"
user_libdir="$libdir"

dnl System Header
AC_CHECK_HEADERS([fenv.h float.h])
AC_CHECK_FUNCS([fesetround])
AX_CXX_HEADER_STDCXX_98
AC_TYPE_LONG_DOUBLE_WIDER

dnl IA Libraries
m4_define([ivlibs],[[filib], [cxsc]])
m4_foreach([var], [ivlibs],[
AC_MSG_CHECKING([whether to use ] m4_translit(var, [a-z], [A-Z]))
AS_IF([test x"$with_[]var" != x -o x"$m4_translit(var, [a-z], [A-Z])[]_ROOT" != x],
[AS_IF([test x"$use_iv_lib" = x],[use_iv_lib=var
AC_MSG_RESULT([yes])], [
AC_MSG_ERROR([More than one interval library specified!])])],
[AC_MSG_RESULT([no])])])
AS_CASE([$use_iv_lib],[filib],[
dnl filib++
SEARCH_HEADER([filib],[interval/interval.hpp], [], [filib])
SEARCH_LIBRARY([filib],[filib::fp_traits<double>::setup();],[],[],[#include <interval/interval.hpp>
],[prim])
IV_CPPFLAGS="$FILIB_CPPFLAGS -DFILIB_VERSION"
IV_LIBS="$FILIB_LIBS $FILIB_LDFLAGS"
],[
dnl cxsc
AS_IF([test x"$use_iv_lib" = x],[AC_MSG_NOTICE([No interval library specified. Assuming --with-cxsc=yes])
	    with_cxsc=yes; use_iv_lib=cxsc
])
CXSC_REQUIRE([2],[2],[4])
CXSC_FIND_LIB
IV_CPPFLAGS="$CXSC_CPPFLAGS -DCXSC_VERSION"
IV_LIBS="$CXSC_LIBS"
])
m4_foreach([var],[ivlibs],[
 	AM_CONDITIONAL(m4_translit(var, [a-z], [A-Z]),[test x"$use_iv_lib" = x"var"])
])
AC_SUBST([IV_CPPFLAGS],[$IV_CPPFLAGS])
AC_SUBST([IV_LIBS],[$IV_LIBS])


dnl DEBUG
AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[enables debugging built]), [enable_debug="$enableval"], [enable_debug=no])
AC_MSG_RESULT([$enable_debug])
AS_IF([test x"$enable_debug" = "xyes"], [AC_DEFINE([DEBUG], [], [Debugging built])
 	    DEBUG_FLAGS="-g -O0"],[AC_DEFINE([NDEBUG], [], [No debugging built])
	    DEBUG_FLAGS=-O3])
AC_SUBST([DEBUG_FLAGS])

dnl DGL Loeser
AC_MSG_CHECKING([whether to build with the ODE solver])
AC_ARG_ENABLE([ode],
AS_HELP_STRING([--disable-ode], [Compile without ODE solver code]), [enable_ode="$enableval"], [enable_ode=yes])
AC_MSG_RESULT([$enable_ode])
AM_CONDITIONAL([ODE_SUPPORT], [test x"$enable_ode" != xno])

dnl Libtool
LT_INIT

libdir=$user_libdir

CXXFLAGS="$USER_CXXFLAGS"
CPPFLAGS="$USER_CPPFLAGS"

dnl Ausgeben
AC_OUTPUT

