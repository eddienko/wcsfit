# CASU_SET_PREFIX
#---------------
# Sets and the directory prefix for the package installation. If no
# directory prefix was given on the command line the default prefix
# is appended to the configure argument list and thus passed to
# the subdirs configure.
AC_DEFUN([CASU_SET_PREFIX],
[
  unset CDPATH
  # make cwd the default for the installation
  AC_PREFIX_DEFAULT(${CASUDIR:-$PWD})

  if test "x$prefix" = "xNONE"; then
    prefix=$ac_default_prefix
    ac_configure_args="$ac_configure_args --prefix $prefix"
    if test -d $PWD/bin; then
        rm -rf $PWD/bin
    fi
    mkdir $PWD/bin
  fi
])

# CASU_CHECK_CFITSIO(version)
#---------------------------
# Checks for the cfitsio library and header files.
AC_DEFUN([CASU_CHECK_CFITSIO],
[

    AC_REQUIRE([CASU_CHECK_THREADS_POSIX])
        
    casu_cfitsio_check_version="$1"
    casu_cfitsio_check_header="fitsio.h"
    casu_cfitsio_check_lib="libcfitsio.a"

    casu_cfitsio_incdirs=""
    casu_cfitsio_libdirs=""
    casu_cfitsio_includes=""
    casu_cfitsio_libraries=""

    AC_ARG_WITH(cfitsio,
                AC_HELP_STRING([--with-cfitsio],
                               [location where cfitsio is installed]),
                [
                    casu_with_cfitsio=$withval
                ])

    AC_ARG_WITH(cfitsio-includes,
                AC_HELP_STRING([--with-cfitsio-includes],
                               [location of the cfitsio header files]),
                casu_with_cfitsio_includes=$withval)

    AC_ARG_WITH(cfitsio-libs,
                AC_HELP_STRING([--with-cfitsio-libs],
                               [location of the cfitsio library]),
                casu_with_cfitsio_libs=$withval)

    AC_ARG_ENABLE(cfitsio-test,
                  AC_HELP_STRING([--disable-cfitsio-test],
                                 [disables checks for the cfitsio library and headers]),
                  casu_enable_cfitsio_test=$enableval,
                  casu_enable_cfitsio_test=yes)


    # We need libpthread for the folloing tests

    if test -z "$LIBPTHREAD"; then
        AC_MSG_ERROR([POSIX thread library was not found on your system! Please check!])
    fi

    
    AC_MSG_CHECKING([for cfitsio])

    if test "x$casu_enable_cfitsio_test" = xyes; then

        # Check for the cfitsio includes

        if test -z "$casu_with_cfitsio_includes"; then
        
            if test -z "$casu_with_cfitsio"; then
            
                # Try some known system locations
                
                casu_cfitsio_incdirs="/opt/cfitsio/include"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/local/include/libcfitsio0"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/local/include/cfitsio"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/local/include"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/include/libcfitsio0"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/include/cfitsio"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs /usr/include"

                test -n "$CFITSIODIR" && \ 
                    casu_cfitsio_incdirs="$CFITSIODIR/include/libcfitsio0 \
                                         $CFITSIODIR/include/cfitsio \
                                         $CFITSIODIR/include \
                                         $casu_cfitsio_incdirs"

                test -n "$CASUDIR" && \
                    casu_cfitsio_incdirs="$CASUDIR/include/libcfitsio0 \
                                         $CASUDIR/include/cfitsio \
                                         $CASUDIR/include \
                                         $casu_cfitsio_incdirs"

            else

                casu_cfitsio_incdirs="$casu_with_cfitsio/include/libcfitsio0"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs $casu_with_cfitsio/include/cfitsio"
                casu_cfitsio_incdirs="$casu_cfitsio_incdirs $casu_with_cfitsio/include"

            fi
            
        else
            casu_cfitsio_incdirs="$casu_with_cfitsio_includes"
        fi

        CASU_FIND_FILE($casu_cfitsio_check_header, $casu_cfitsio_incdirs,
                      casu_cfitsio_includes)


        # Check for the cfitsio library

        if test -z "$casu_with_cfitsio_libs"; then

            if test -z "$casu_with_cfitsio"; then
            
                # Try some known system locations

                casu_cfitsio_libdirs="/opt/cfitsio/lib"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/local/lib64"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/local/lib"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/local/lib32"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/lib64"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/lib"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs /usr/lib32"

                test -n "$CFITSIODIR" && \
                    casu_cfitsio_libdirs="$CFITSIODIR/lib64 $CFITSIODIR/lib \
                                         $CFITSIODIR/lib32 $casu_cfitsio_libdirs"

                test -n "$CASUDIR" && \
                    casu_cfitsio_libdirs="$CASUDIR/lib64 $CASUDIR/lib $CASUDIR/lib32 \
                                         $casu_cfitsio_libdirs"

            else

                casu_cfitsio_libdirs="$casu_with_cfitsio/lib64"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs $casu_with_cfitsio/lib"
                casu_cfitsio_libdirs="$casu_cfitsio_libdirs $casu_with_cfitsio/lib32"

            fi
            
        else
            casu_cfitsio_libdirs="$casu_with_cfitsio_libs"
        fi

        CASU_FIND_FILE($casu_cfitsio_check_lib, $casu_cfitsio_libdirs,
                      casu_cfitsio_libraries)


        if test x"$casu_cfitsio_includes" = xno || \
            test x"$casu_cfitsio_libraries" = xno; then
            casu_cfitsio_notfound=""

            if test x"$casu_cfitsio_includes" = xno; then
                if test x"$casu_cfitsio_libraries" = xno; then
                    casu_cfitsio_notfound="(headers and libraries)"
                else
                    casu_cfitsio_notfound="(headers)"
                fi
            else
                casu_cfitsio_notfound="(libraries)"
            fi

            AC_MSG_ERROR([cfitsio $casu_cfitsio_notfound was not found on your system. Please check!])
        else
            AC_MSG_RESULT([libraries $casu_cfitsio_libraries, headers $casu_cfitsio_includes])
        fi

        # Set up the symbols

        # Add '-lz' to the static library symbol, as distributors apparently
        # remove the libz code from the cfitsio sources.
        
        CFITSIO_INCLUDES="-I$casu_cfitsio_includes"
        CFITSIO_LDFLAGS="-L$casu_cfitsio_libraries"
        LIBCFITSIO="-lcfitsio"
        LIBCFITSIO_STATIC="$casu_cfitsio_libraries/$casu_cfitsio_check_lib"
        
        # Do not add redundant libraries        
        echo $LIBS | grep -q -e '-lm' || LIBS="-lm $LIBS" 
        

        # Check whether cfitsio can be used

        AC_MSG_CHECKING([whether cfitsio can be used])
        AC_LANG_PUSH(C)
        
        casu_cfitsio_cflags_save="$CFLAGS"
        casu_cfitsio_ldflags_save="$LDFLAGS"
        casu_cfitsio_libs_save="$LIBS"

        CFLAGS="$CFITSIO_INCLUDES $CFLAGS"
        LDFLAGS="$CFITSIO_LDFLAGS $LDFLAGS"
        LIBS="$LIBCFITSIO_STATIC $LIBPTHREAD -lm"
        
        AC_LINK_IFELSE([AC_LANG_PROGRAM(
                       [[
                       #include <fitsio.h>
                       ]],
                       [
                       float v;
                       fits_get_version(&v);
                       ])],
                       [casu_cfitsio_is_usable="yes"],
                       [casu_cfitsio_is_usable="no"])

        AC_MSG_RESULT([$casu_cfitsio_is_usable])
        
        AC_LANG_POP(C)
        
        CFLAGS="$casu_cfitsio_cflags_save"
        LDFLAGS="$casu_cfitsio_ldflags_save"
        LIBS="$casu_cfitsio_libs_save"

        if test x"$casu_cfitsio_is_usable" = xno; then
            AC_MSG_ERROR([Linking with cfitsio failed! Please check architecture!])
        fi
 
        
        # Check cfitsio version

        AC_MSG_CHECKING([for a cfitsio version >= $casu_cfitsio_check_version])
        
        AC_LANG_PUSH(C)
        
        casu_cfitsio_cflags_save="$CFLAGS"
        casu_cfitsio_ldflags_save="$LDFLAGS"
        casu_cfitsio_libs_save="$LIBS"

        CFLAGS="$CFITSIO_INCLUDES $CFLAGS"
        LDFLAGS="$CFITSIO_LDFLAGS $LDFLAGS"
        LIBS="$LIBCFITSIO_STATIC $LIBPTHREAD -lm"
        
        AC_RUN_IFELSE([AC_LANG_PROGRAM(
                      [[
                      #include <stdio.h>
                      #include <fitsio.h>
                      ]],
                      [
                      int vlib = 0;
                      int vmin = (int)(1000. * $casu_cfitsio_check_version + 0.5);
                       
                      float v;

                      fits_get_version(&v);
                      vlib = (int)(v * 1000 + 0.5);
                                            
                      FILE* f = fopen("conftest.out", "w");
                      fprintf(f, "%5.3f\n", v);
                      fclose(f);
                      
                      if (vlib < vmin) {
                          return 1;
                      }

                      return 0;
                      ])],
                      [casu_cfitsio_version="`cat conftest.out`"],
                      [
                       casu_cfitsio_version="no";
                       casu_cfitsio_version_found="`cat conftest.out`"
                      ])

        AC_MSG_RESULT([$casu_cfitsio_version])
        
        AC_LANG_POP(C)
        
        CFLAGS="$casu_cfitsio_cflags_save"
        LDFLAGS="$casu_cfitsio_ldflags_save"
        LIBS="$casu_cfitsio_libs_save"

        if test x"$casu_cfitsio_version" = xno; then
            AC_MSG_ERROR([Installed cfitsio ($casu_cfitsio_version_found) is too old. Please update to version $casu_cfitsio_check_version or newer.])
        fi
 
        
        # Check whether cfitsio has large file support
        
        AC_LANG_PUSH(C)
        
        casu_cfitsio_cflags_save="$CFLAGS"
        casu_cfitsio_ldflags_save="$LDFLAGS"
        casu_cfitsio_libs_save="$LIBS"

        CFLAGS="$CFITSIO_INCLUDES $CFLAGS"
        LDFLAGS="$CFITSIO_LDFLAGS $LDFLAGS"
        LIBS="$LIBCFITSIO_STATIC -lm $LIBPTHREAD"
        
        AC_MSG_CHECKING([whether cfitsio provides fits_hdu_getoff()])

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                          [[
                          #include <fitsio.h>
                          ]],
                          [
                          fitsfile f;
                          int sts;
                          fits_get_hduoff(&f, NULL, NULL, NULL, &sts);
                          ])],
                          [casu_cfitsio_have_fits_get_hduoff="yes"],
                          [casu_cfitsio_have_fits_get_hduoff="no"])

        AC_MSG_RESULT([$casu_cfitsio_have_fits_get_hduoff])
        
        AC_MSG_CHECKING([whether cfitsio provides fits_get_hduaddrll()])

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                          [[
                          #include <fitsio.h>
                          ]],
                          [
                          fitsfile f;
                          int sts;
                          fits_get_hduaddrll(&f, NULL, NULL, NULL, &sts);
                          ])],
                          [casu_cfitsio_have_fits_get_hduaddrll="yes"],
                          [casu_cfitsio_have_fits_get_hduaddrll="no"])

        AC_MSG_RESULT([$casu_cfitsio_have_fits_get_hduaddrll])

        AC_LANG_POP(C)
        
        CFLAGS="$casu_cfitsio_cflags_save"
        LDFLAGS="$casu_cfitsio_ldflags_save"
        LIBS="$casu_cfitsio_libs_save"


        # Check whether cfitsio is thread-safe
        
        AC_MSG_CHECKING([whether cfitsio requires libpthread])
        AC_LANG_PUSH(C)
        
        casu_cfitsio_cflags_save="$CFLAGS"
        casu_cfitsio_ldflags_save="$LDFLAGS"
        casu_cfitsio_libs_save="$LIBS"

        CFLAGS="$CFITSIO_INCLUDES $CFLAGS"
        LDFLAGS="$CFITSIO_LDFLAGS $LDFLAGS"
        LIBS="$LIBCFITSIO_STATIC -lm"
        
        
        casu_cfitsio_provides_pthread="no"
        
        AC_LINK_IFELSE([AC_LANG_PROGRAM(
                       [[
                       #include <fitsio.h>
                       ]],
                       [
                       float v;
                       fits_get_version(&v);
                       ])],
                       [casu_cfitsio_requires_pthread="no"],
                       [casu_cfitsio_requires_pthread="undefined"])

		if test x"$casu_cfitsio_requires_pthread" = xno; then
		
	    	# If libpthread is not required this means either cfitsio is
    		# not compiled with thread support, or the library dependencies
    		# are compiled into cfitsio.

			# Check whether shared library dependencies are present
			# Don't use pthread_mutex_init/destroy here, since glibc
			# provides these symbols too! Sigh!
		
	        AC_LINK_IFELSE([AC_LANG_PROGRAM(
		                   [[
	    	               #include <pthread.h>
    	    	           #include <fitsio.h>
        	    	       ]],
            	    	   [                	       
	                	   float v;
    	            	   pthread_mutexattr_t attrb;
    	            	   pthread_mutexattr_init(&attrb);
            	           fits_get_version(&v);
                		   pthread_mutexattr_destroy(&attrb);
	                	   ])],
	    	               [casu_cfitsio_provides_pthread="yes"],
    	    	           [casu_cfitsio_provides_pthread="no"])
            	           
        else
        
            LIBS="$LIBCFITSIO_STATIC -lm $LIBPTHREAD"
        
            AC_LINK_IFELSE([AC_LANG_PROGRAM(
                           [[
                           #include <fitsio.h>
                           ]],
                           [
                           float v;
                           fits_get_version(&v);
                           ])],
                           [casu_cfitsio_requires_pthread="yes"],
                           AC_MSG_FAILURE([Cannot link with cfitsio! Please check!]))
          
        fi
                       
        AC_MSG_RESULT([$casu_cfitsio_requires_pthread])
        
        AC_LANG_POP(C)
        
        CFLAGS="$casu_cfitsio_cflags_save"
        LDFLAGS="$casu_cfitsio_ldflags_save"
        LIBS="$casu_cfitsio_libs_save"


		AC_MSG_CHECKING([whether cfitsio was compiled with thread support])
		
		if test x"$casu_cfitsio_requires_pthread" = xyes || \
			test x"$casu_cfitsio_provides_pthread" = xyes; then
			casu_cfitsio_is_thread_safe=yes
		else
			casu_cfitsio_is_thread_safe=no
		fi
		
		AC_MSG_RESULT([$casu_cfitsio_is_thread_safe])
		
		
        # Set compiler flags and libraries
        
        if test x"$casu_cfitsio_have_fits_get_hduoff" = xyes || \
          test x"$casu_cfitsio_have_fits_get_hduaddrll" = xyes; then

            CFLAGS="-D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 $CFLAGS"
            
            if test x"$casu_cfitsio_have_fits_get_hduoff"; then
                AC_DEFINE([HAVE_FITS_GET_HDUOFF], [1],
                          [Define if you have the `fits_get_hduoff' function])
            else
                AC_DEFINE([HAVE_FITS_GET_HDUADDRLL],  [1],
                          [Define if you have the `fits_get_hduaddrll' function])
            fi
                    
        fi
                
        if test x"$casu_cfitsio_requires_pthread" = xyes; then
            echo $LIBS | grep -q -e "$LIBPTHREAD" || LIBS="$LIBPTHREAD $LIBS"
        fi
        
    else
        AC_MSG_RESULT([disabled])
        AC_MSG_WARN([cfitsio checks have been disabled! This package may not build!])
        CFITSIO_INCLUDES=""
        CFITSIO_LDFLAGS=""
        LIBCFITSIO=""
        
        casu_cfitsio_is_thread_safe="undefined"
        casu_cfitsio_requires_pthread="undefined"
    fi

	AC_CACHE_VAL(casu_cv_cfitsio_requires_pthread,
	             casu_cv_cfitsio_requires_pthread=$casu_cfitsio_requires_pthread)
	AC_CACHE_VAL(casu_cv_cfitsio_is_thread_safe,
	             casu_cv_cfitsio_is_thread_safe=$casu_cfitsio_is_thread_safe)
	             
    AC_SUBST(CFITSIO_INCLUDES)
    AC_SUBST(CFITSIO_LDFLAGS)
    AC_SUBST(LIBCFITSIO)

])

# CASU_CHECK_WCS
#--------------
# Checks for the wcs library and header files.
AC_DEFUN([CASU_CHECK_WCS],
[
    AC_MSG_CHECKING([for wcs])

    casu_wcs_check_header="wcslib/wcslib.h"
    casu_wcs_check_lib="libwcs.a"

    casu_wcs_includes=""
    casu_wcs_libraries=""

    AC_ARG_WITH(wcs,
                AC_HELP_STRING([--with-wcs],
                               [location where wcs is installed]),
                [
                    casu_with_wcs_includes=$withval/include
                    casu_with_wcs_libs=$withval/lib
                ])

    # Check for the wcs includes
    if test -z "$casu_with_wcs_includes"; then
        test -n "$WCSDIR" && casu_wcs_incdirs="$WCSDIR/include"
    else
        casu_wcs_incdirs="$casu_with_wcs_includes"
    fi
    CASU_FIND_FILE($casu_wcs_check_header, $casu_wcs_incdirs, casu_wcs_includes)

    # Check for the wcs library
    if test -z "$casu_with_wcs_libs"; then
        test -n "$WCSDIR" && casu_wcs_libdirs="$WCSDIR/lib"
    else
        casu_wcs_libdirs="$casu_with_wcs_libs"
    fi
    CASU_FIND_FILE($casu_wcs_check_lib, $casu_wcs_libdirs, casu_wcs_libraries)

    if test x"$casu_wcs_includes" = xno || test x"$casu_wcs_libraries" = xno; then
        AC_MSG_WARN([wcs was not found on your system.])
    else
        AC_MSG_RESULT([libraries $casu_wcs_libraries, headers $casu_wcs_includes])
        # Attempt to check the version by checking the include files
#        casu_wcs_check_vers4=`grep "WCSLIB 4.[3-9] - an implementation of the FITS WCS standard" $casu_wcs_includes/wcslib/wcslib.h`
#        if test -z "$casu_wcs_check_vers4" ; then
#            casu_wcs_check_vers5=`grep "WCSLIB [5-9].[0-9] - an implementation of the FITS WCS standard" $casu_wcs_includes/wcslib/wcslib.h`
#            if test -z "$casu_wcs_check_vers5" ; then
#                AC_MSG_WARN([wcs version seems to be older than 4.3])
#            fi
#        fi
        AC_DEFINE_UNQUOTED(CASU_WCS_INSTALLED, 1, [Defined if WCS is available])
        # Set up the symbols
        WCS_INCLUDES="-I$casu_wcs_includes/wcslib"
        WCS_LDFLAGS="-L$casu_wcs_libraries"
        LIBWCS="-lwcs"

        AC_SUBST(WCS_INCLUDES)
        AC_SUBST(WCS_LDFLAGS)
        AC_SUBST(LIBWCS)
    fi

])

# CASU_FIND_FILE(file, directories, variable)
#------------------------------------------
# Search for file in directories. Set variable to the first location
# where file was found, if file is not found at all variable is set to NO.
AC_DEFUN([CASU_FIND_FILE],
[
    $3=no

    for i in $2; do
        for j in $1; do

            echo "configure: __oline__: $i/$j" >&AC_FD_CC

            if test -r "$i/$j"; then
                echo "taking that" >&AC_FD_CC
                $3=$i
                break 2
            fi
        done
    done
])

# CASU_CHECK_THREADS_POSIX
#------------------------
# Check whether the POSIX threads are available. The cached result is
# set to 'yes' if either the compiler supports the '-pthread' flag, or linking
# with the POSIX thread library works, and the header file defining the POSIX
# threads symbols is present. If POSIX threads are not supported, the
# result is set to 'no'. Whether the compiler supports POSIX threads,
# or whether the library, and the header are available is stored in cache
# variables.  
AC_DEFUN([CASU_CHECK_THREADS_POSIX],
[
    AC_REQUIRE([AC_PROG_CC])

    CASU_PROG_CC_FLAG([pthread], [], [])
    
    AC_CHECK_LIB([pthread], [pthread_create],
                 [casu_threads_have_libpthread=yes],
                 [casu_threads_have_libpthread=no])

    AC_CHECK_HEADER([pthread.h],
                    [casu_threads_have_pthread_h=yes],
                    [casu_threads_have_pthread_h=no])

    if test x"$casu_threads_have_pthread_h" != xyes; then
        casu_threads_posix=no
    else
        if test x"$casu_threads_have_libpthread" != xyes && \
          test x"$casu_cv_prog_cc_pthread" != xyes; then
            casu_threads_posix=no
        else
            casu_threads_posix=yes
        fi
    fi
    
    
    # Setup the POSIX thread symbols

    if test x"$casu_threads_have_pthread_h" = xyes; then
        AC_DEFINE([HAVE_PTHREAD_H], [1],
                  [Define to 1 if you have <pthread.h>.])
    fi
    
    if test x"$casu_threads_posix" = xyes; then
    
        if test x"$casu_cv_prog_cc_pthread" = xyes; then
            PTHREAD_CFLAGS="-pthread"
        else
            PTHREAD_CFLAGS=""
        fi
        
        if test x"$casu_threads_have_libpthread" = xyes; then
            LIBPTHREAD="-lpthread"
        else
            LIBPTHREAD=""
        fi
        
    fi  

    AC_CACHE_VAL(casu_cv_threads_posix_header,
                 casu_cv_threads_posix_header=$casu_threads_have_pthread_h)          
    AC_CACHE_VAL(casu_cv_threads_posix_lib,
                 casu_cv_threads_posix_lib=$casu_threads_have_libpthread)          
    AC_CACHE_VAL(casu_cv_threads_posix_flags,
                 casu_cv_threads_posix_flags=$casu_cv_prog_cc_pthread)          
    AC_CACHE_VAL(casu_cv_threads_posix,
                 casu_cv_threads_posix=$casu_threads_posix)

    AC_SUBST(PTHREAD_CFLAGS)
    AC_SUBST(LIBPTHREAD)
    
])

# CASU_PROG_CC_FLAG(FLAG, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#-----------------------------------------------------------------
AC_DEFUN([CASU_PROG_CC_FLAG],
[
    AC_REQUIRE([AC_PROG_CC])

    flag=`echo $1 | sed 'y%.=/+-%___p_%'`
    AC_CACHE_CHECK([whether $CC supports -$1],
                   [casu_cv_prog_cc_$flag],
                   [
                       eval "casu_cv_prog_cc_$flag=no"
                       AC_LANG_PUSH(C)

                       echo 'int main() { return 0; }' >conftest.$ac_ext

                       try_compile="`$CC -$1 -c conftest.$ac_ext 2>&1`"
                       if test -z "$try_compile"; then
                           try_link="`$CC -$1 -o conftest$ac_exeext \
                                    conftest.$ac_ext 2>&1`"
                           if test -z "$try_link"; then
                               eval "casu_cv_prog_cc_$flag=yes"
                           fi
                       fi
                       rm -f conftest*

                       AC_LANG_POP(C)
                   ])

    if eval "test \"`echo '$casu_cv_prog_cc_'$flag`\" = yes"; then
        :
        $2
    else
        :
        $3
    fi
])
