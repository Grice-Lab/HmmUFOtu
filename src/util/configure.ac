#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([util],[1.0],[e00011027@gmail.com])
AC_CONFIG_HEADERS([config.h])
LT_INIT

# Checks for programs.
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_FUNCS([strerror])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

# Checks for autoconf-archive macros
AX_CXX_COMPILE_STDCXX(11)
AX_BOOST_BASE([1.50.0])

# Checks for library functions.

AC_CONFIG_FILES(
  [Makefile
  UCSC/Makefile
  test/Makefile]
)
AC_OUTPUT
