AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES)

noinst_HEADERS =    casu_wcs.h \
                    cir_memory.h \
		    misc.h \
                    tools.h 


noinst_LIBRARIES = libcasutools.a

libcasutools_a_SOURCES = casu_wcs.c \
			 cir_catcoord.c \
			 cir_memory.c \
			 misc.c 

