#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)

AC_INIT(BEETL, 1.0, acox@illumina.com)
AM_INIT_AUTOMAKE([subdir-objects])

AC_CONFIG_SRCDIR([src/frontends/BeetlBwt.cpp])
AC_CONFIG_HEADER([config.h])
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CC
AC_PROG_RANLIB

# Checks for C++11 compatibility
m4_include([m4/ax_cxx_compile_stdcxx_11.m4])
AX_CXX_COMPILE_STDCXX_11

# Checks for openmp
m4_include([m4/openmp.m4])
AC_OPENMP

# Checks for libraries.
# Seqan
AC_ARG_WITH([seqan],
    AS_HELP_STRING([--with-seqan=PATH], [Path to Seqan package]))

SEQAN_PATH="not found (Metagenomics tool 'makeBWTSkew' will be unavailable - download at www.seqan.de)"
SEQAN_CXXFLAGS=
AS_IF([test "x$with_seqan" != "x" -a "x$with_seqan" != "xno"],
      [AC_CHECK_FILE([$with_seqan/include/seqan/index.h],
        [AC_SUBST([SEQAN_PATH], ["$with_seqan"])
         AC_SUBST([SEQAN_CXXFLAGS], ["-I$with_seqan/include -DHAVE_SEQAN"])
        ],
        [AC_MSG_ERROR([Seqan header file not found in $with_seqan])]
       )
      ]
     )

# Boost
AC_ARG_WITH([boost],
    AS_HELP_STRING([--with-boost=PATH], [Path to Boost package]))

BOOST_PATH="not found (but only used for experimental stuff)"
BOOST_CXXFLAGS=
BOOST_LDADD=
BOOST_LDFLAGS=
AS_IF([test "x$with_boost" != "x" -a "x$with_boost" != "xno"],
      [AC_CHECK_FILE([$with_boost/include/boost/version.hpp],
        [AC_SUBST([BOOST_PATH], ["$with_boost"])
         AC_SUBST([BOOST_CXXFLAGS], ["-I$with_boost/include -DHAVE_BOOST"])
         AC_SUBST([BOOST_LDFLAGS], ["-L$with_boost/lib"])
         AC_SUBST([BOOST_LDADD], ["-lboost_system -lboost_thread"])
        ],
        [AC_MSG_ERROR([Boost header file not found in $with_boost])]
       )
      ]
     )

# zlib
#AC_CHECK_LIB(zlib, zlib,
#    [],
#    [AC_MSG_ERROR(["Zlib not detected! Maybe try yum/apt-get install zlib1g-dev"])]
#)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_STAT
AC_CHECK_FUNCS([floor gettimeofday memset mkdir strchr strtoul strcasecmp posix_fadvise sysinfo sysctlbyname proc_pidpath])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/frontends/Makefile
                 src/frontends/beetl
                 src/metagenomics/Makefile
                 scripts/Makefile
                 scripts/lcp/Makefile
                 scripts/tumourNormalFiltering/Makefile
                 scripts/BEETL-fastq/Makefile
                 scripts/perSubsetDistributedLauncher/Makefile
                 share/Makefile
                 tests/Makefile
                 doc/Makefile
                ])
AC_OUTPUT

echo "
Configure summary:
    Packages for metagenomics:
        Seqan ...........: $SEQAN_PATH

    Optional packages:
        Boost ...........: $BOOST_PATH

    OpenMP:
        Compiler flags...: $OPENMP_CXXFLAGS
"
