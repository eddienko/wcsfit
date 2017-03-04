AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools

dribble_SOURCES = dribble_main.c

bin_PROGRAMS = dribble

dribble_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
dribble_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

