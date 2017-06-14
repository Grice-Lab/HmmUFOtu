dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(HmmUFOtu, 1.0)
LT_INIT([shared static])

AC_CANONICAL_SYSTEM
AC_CHECK_HEADERS([src/HmmUFOtu.h])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_INT32_T
AC_FUNC_MALLOC
AC_PROG_LIBTOOL
AC_PROG_RANLIB
AC_PROG_LN_S
AM_INIT_AUTOMAKE()

#AC_LANG(C++)
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])

AC_CHECK_FUNCS([sqrt isnan floor])
AC_CHECK_FUNCS([strerror])

AC_CONFIG_FILES([Makefile
								 data/Makefile
								 src/Makefile
								 src/math/Makefile
								 test/Makefile])

AC_CONFIG_SUBDIRS([src/libcds
                   src/libdivsufsort])
AC_OPENMP

AC_OUTPUT
