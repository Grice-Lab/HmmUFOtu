AC_PREREQ([2.69])
AC_INIT([libdivsufsort], [2.0.0], [e00011027@gmail.com])
AC_CONFIG_SRCDIR([lib/utils.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_INIT_AUTOMAKE()

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h memory.h stddef.h stdint.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_PROG_RANLIB

AC_CONFIG_FILES([Makefile lib/Makefile examples/Makefile])
AC_OUTPUT
