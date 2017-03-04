AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

INCLUDES = $(CFITSIO_INCLUDES) $(WCS_INCLUDES) -I ../tools
noinst_HEADERS = ap.h cir_classify.h errcodes.h floatmath.h imcore.h \
	 	 imcore_list.h imcore_radii.h imcore_version.h util.h

imcore_SOURCES = imcore_main.c \
		 apclust.c \
		 apinit.c \
		 apline.c \
		 areals.c \
		 cir_classify.c \
		 create_table_1.c \
		 create_table_2.c \
		 create_table_3.c \
		 create_table_4.c \
		 create_table_6.c \
		 create_table.c \
		 imcore_background.c \
		 imcore_conf.c \
		 imcore_extend.c \
		 imcore_filter.c \
		 imcore_list.c \
		 imcore_opm.c \
		 imcore_overlp.c \
		 imcore_phopt.c \
		 imcore_radii.c \
		 imcore_rdbuf_mef.c \
		 moments.c \
		 polynm.c \
		 seeing.c \
		 solve.c \
		 terminate.c 

imcore_list_SOURCES = imcore_main_list.c \
		 apclust.c \
		 apinit.c \
		 apline.c \
		 areals.c \
		 cir_classify.c \
		 create_table_1.c \
		 create_table_2.c \
		 create_table_3.c \
		 create_table_4.c \
		 create_table_6.c \
		 create_table.c \
		 imcore_background.c \
		 imcore_conf.c \
		 imcore_extend.c \
		 imcore_filter.c \
		 imcore_list.c \
		 imcore_opm.c \
		 imcore_overlp.c \
		 imcore_phopt.c \
		 imcore_radii.c \
		 imcore_rdbuf_mef.c \
		 moments.c \
		 polynm.c \
		 seeing.c \
		 solve.c \
		 terminate.c 

classify_SOURCES = classify_main.c \
                   cir_classify.c

bin_PROGRAMS = imcore imcore_list classify

imcore_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
imcore_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm
imcore_list_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
imcore_list_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm
classify_LDFLAGS = -L../tools $(CFITSIO_LDFLAGS) $(WCS_LDFLAGS) 
classify_LDADD = -lcasutools $(LIBCFITSIO) $(LIBWCS) -lm

