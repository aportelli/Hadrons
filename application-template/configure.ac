AC_PREREQ([2.69])
AC_INIT([my-hadrons-app], [0.1], [name@domain.suf])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [try this for a non-standard install prefix of Grid])],
    [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_grid/include"]
    [LDFLAGS="$LDFLAGS -L$with_grid/lib"])
AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [try this for a non-standard install prefix of Hadrons])],
    [PATH="$with_hadrons/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_hadrons/include"]
    [LDFLAGS="$LDFLAGS -L$with_hadrons/lib"])
AC_CHECK_PROG([HADCONF],[hadrons-config],[yes])
if test x"$HADCONF" != x"yes" ; then
    AC_MSG_ERROR([hadrons-config not found])
fi
if test x"$CXX" == x ; then
    CXX="`hadrons-config --cxx`"
elif test "$CXX" != "`hadrons-config --cxx`" ; then
    AC_MSG_WARN([CXX differs from that reported by hadrons-config])
fi
if test x"$CXXLD" == x ; then
    CXXLD="`hadrons-config --cxxld`"
elif test "$CXXLD" != "`hadrons-config --cxxld`" ; then
    AC_MSG_WARN([CXXLD differs from that reported by hadrons-config])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
AC_LANG([C++])

CXXFLAGS="$CXXFLAGS `hadrons-config --cxxflags`"
LDFLAGS="$LDFLAGS `hadrons-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS=" -lHadrons $LIBS `hadrons-config --libs`"

AC_MSG_CHECKING([that a minimal Grid/Hadrons program compiles]);
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    #include <Hadrons/Application.hpp>

    using namespace Grid;
    using namespace Hadrons;

    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Grid_finalize();
        Application application;
        application.run();

        return 0;
    }

    ]])],
        [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Could not compile a minimal Hadrons program])])

AC_SUBST([CXXLD])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
