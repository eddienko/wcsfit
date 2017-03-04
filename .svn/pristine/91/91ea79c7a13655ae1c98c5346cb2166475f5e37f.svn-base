AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools

grout_SOURCES = grout_main.c

ungrout_SOURCES = ungrout_main.c

bin_PROGRAMS = grout ungrout

grout_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
grout_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm
ungrout_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
ungrout_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

