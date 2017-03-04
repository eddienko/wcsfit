AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = vstpickup.h

vstpickup_SOURCES = vstpickup_main.c \
		    cir_filt1d.c \
                    cir_vstpickup.c

bin_PROGRAMS = vstpickup

vstpickup_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
vstpickup_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

