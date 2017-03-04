AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = mosaic.h

mosaic_SOURCES = mosaic_main.c \
		  cir_mosaic.c

bin_PROGRAMS = mosaic

mosaic_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
mosaic_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

