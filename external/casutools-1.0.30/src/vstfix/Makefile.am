AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = vstfix.h

vstfix_SOURCES = vstfix_main.c \
	         cir_vstfix.c

vsttab_a2f_SOURCES = vsttab_a2f_main.c

bin_PROGRAMS = vstfix vsttab_a2f

vstfix_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
vstfix_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

vsttab_a2f_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
vsttab_a2f_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

