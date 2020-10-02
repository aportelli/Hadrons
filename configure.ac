AC_PREREQ([2.63])
AC_INIT(Hadrons,[v1.0],[antonin.portelli@ed.ac.uk],[hadrons])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([Hadrons/Global.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_EXTRA_RECURSIVE_TARGETS([tests])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [try this for a non-standard install prefix of Grid])],
    [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_grid/include"]
    [LDFLAGS="$LDFLAGS -L$with_grid/lib"])

AC_CHECK_PROG([GRIDCONF],[grid-config],[yes])
if test x"$GRIDCONF" != x"yes" ; then
    AC_MSG_ERROR([grid-config not found])
fi
if test x"$CXX" == x ; then
    CXX="`grid-config --cxx`"
elif test "$CXX" != "`grid-config --cxx`" ; then
    AC_MSG_WARN([CXX differs from that reported by grid-config])
fi
if test x"$CXXLD" == x ; then
    CXXLD="`grid-config --cxxld`"
elif test "$CXXLD" != "`grid-config --cxxld`" ; then
    AC_MSG_WARN([CXXLD differs from that reported by grid-config])
fi
CXXFLAGS="$CXXFLAGS `grid-config --cxxflags`"
LDFLAGS="$LDFLAGS `grid-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS=" -ldl -lGrid $LIBS `grid-config --libs`"

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR
AC_LANG([C++])

AC_MSG_CHECKING([that a minimal Grid program compiles]);
AC_LINK_IFELSE(
	[AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    
    using namespace Grid;
    
    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Grid_finalize();
        
        return 0;
    }
    
    ]])],
	[AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([impossible to compile a minimal Grid program])])

HADRONS_CXX="$CXX"
HADRONS_CXXLD="$CXXLD"
HADRONS_CXXFLAGS="$CXXFLAGS"
HADRONS_LDFLAGS="$LDFLAGS"
HADRONS_LIBS="$LIBS"
HADRONS_SHA=`git rev-parse HEAD`
HADRONS_BRANCH=`git rev-parse --abbrev-ref HEAD`

AC_SUBST([CXXLD])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([HADRONS_CXX])
AC_SUBST([HADRONS_CXXLD])
AC_SUBST([HADRONS_CXXFLAGS])
AC_SUBST([HADRONS_LDFLAGS])
AC_SUBST([HADRONS_LIBS])
AC_SUBST([HADRONS_SHA])
AC_SUBST([HADRONS_BRANCH])

AC_CONFIG_FILES([hadrons-config], [chmod +x hadrons-config])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([Hadrons/Makefile])
AC_CONFIG_FILES([utilities/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_OUTPUT
