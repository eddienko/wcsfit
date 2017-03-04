AC_INIT([CASUtools], [1.0.30], [Jim Lewis <jrl@ast.cam.ac.uk>], [casutools], [CASU http://casu.ast.cam.ac.uk/])
AC_DEFINE([PACKAGE_DATE],["$Date: 2014/11/03 08:24:01 $"],["Package release date"])
AM_INIT_AUTOMAKE
AC_PROG_CC

# Checks for header files.

AC_CHECK_HEADERS([fcntl.h float.h limits.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/socket.h unistd.h values.h])

# Checks for typedefs, structures, and compiler characteristics.

AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T

# Checks for library functions.

AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero floor gethostbyname memmove memset mkdir pow socket sqrt strdup strerror strrchr strstr strtol])

AC_PROG_RANLIB
CASU_CHECK_CFITSIO
CASU_CHECK_WCS
CASU_SET_PREFIX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile src/tools/Makefile src/dribble/Makefile src/grout/Makefile src/imcore/Makefile src/imstack/Makefile src/wcsfit/catsearch/Makefile src/wcsfit/Makefile src/nebuliser/Makefile src/mosaic/Makefile src/vstpickup/Makefile src/vstfix/Makefile])
AC_OUTPUT
