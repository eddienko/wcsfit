AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = imstack.h

imstack_SOURCES = imstack_main.c \
		  cir_imstack_cat.c

bin_PROGRAMS = imstack

imstack_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
imstack_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

