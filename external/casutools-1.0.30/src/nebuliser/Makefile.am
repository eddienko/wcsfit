AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = nebuliser.h

nebuliser_SOURCES = nebuliser_main.c \
		    cir_2dfilt.c \
                    filtsubs.c

bin_PROGRAMS = nebuliser

nebuliser_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
nebuliser_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

