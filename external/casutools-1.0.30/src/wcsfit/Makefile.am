AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

SUBDIRS = catsearch

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools -I catsearch
noinst_HEADERS = wcsfit.h

wcsfit_SOURCES = wcsfit_main.c \
		 cir_getstds.c \
		 cir_matchstds.c \
		 cir_platesol.c \
	 	 wcsfit.c

bin_PROGRAMS = wcsfit

wcsfit_LDFLAGS = -L../tools -Lcatsearch $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
wcsfit_LDADD = -lcasutools -lcatsearch $(LIBCFITSIO) $(LIBWCS) -lm

