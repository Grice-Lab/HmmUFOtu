dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(HmmUFOtu, 1.0)
LT_INIT([shared static])

AC_CANONICAL_SYSTEM
AC_CANONICAL_HOST

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_INT32_T
AC_PROG_LIBTOOL
AC_PROG_RANLIB
#LT_INIT
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADER_STDBOOL

AM_INIT_AUTOMAKE()

#AC_LANG(C++)
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])

AC_CHECK_FUNCS([sqrt isnan floor])
AC_CHECK_FUNCS([strerror])

#AC_FUNC_MALLOC

AX_CHECK_ZLIB()
AX_BOOST_BASE([1.50.0])
AX_BOOST_IOSTREAMS()

AC_CONFIG_FILES([Makefile
				 data/Makefile
				 src/Makefile
				 src/math/Makefile
				 test/Makefile])

AC_CONFIG_SUBDIRS([src/libcds
                   src/libdivsufsort])
AC_OPENMP

AC_OUTPUT
