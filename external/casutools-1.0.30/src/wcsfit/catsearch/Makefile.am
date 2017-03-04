AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

noinst_HEADERS = common.h cvtunit.h extract.h floatmath.h psc.h usnob.h util.h

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES)

noinst_LIBRARIES = libcatsearch.a

libcatsearch_a_SOURCES = common.c \
			 cvtunit.c \
			 extract.c \
			 psc.c \
			 usnob.c \
			 util.c

