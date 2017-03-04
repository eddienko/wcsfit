/*

$Id: cir_getstds.c,v 1.3 2012/08/13 10:48:12 jim Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <fcntl.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <wcsfit.h>
#include <tools.h>

#include <extract.h>

#define BUFMAX 32768
#define CASU  "archive.ast.cam.ac.uk"
#define CDS   "webviz.u-strasbg.fr"
#define UKIRT "www.ukirt.jach.hawaii.edu"
#define PORT 80

#define TWOMASS "II/246"
#define USNOB   "I/284"
#define LANDOLT "II/183A"
#define PPMXL "I/317"

#define LOCALSEARCH    1
#define INTERNETSEARCH 2
#define CATSEARCH      3

#define CACHEDIR "catcache"
#define CACHEIND "catcache/index"

static char *ttype[2] = {"X_coordinate","Y_coordinate"};
static char *tform1[2] = {"F10.2","F10.2"};
static char *tform2[2] = {"1E","1E"};

static char *form_request(float, float, float, float, float, char *);
static char *url_encode(char *);
static int connect_site(char *, int *, char *);
static int send_request(int, char *, char *);
static int get_response(int, char *, char *);
static int get_local(char *, float, float, float, float, char *, char *);
static int check_cache(char *, float, float, float, float, char *);
static void addto_cache(char *, char *, float, float, float, float);
static void tidy();

static int cir_coverage(char *infile, int fudge, float *ra1, float *ra2, 
                        float *dec1, float *dec2, char *errmsg);
/*+
 *  Name:
 *      cir_getstds
 *
 *  Purpose:
 *      Get standards stars from a specified source.
 *
 *  Description:
 *      Standard star catalogues are searched for objects that may
 *      appear on a input frame.  The WCS coverage for the frame is
 *      calculated and the catalogues are searched for objects in that range.
 *      The catalogue sources are specified in the input arguments.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          The name of the input FITS image. A valid WCS must be present
 *          in the header.
 *      outfile = char * (Given)
 *          The name for the output fits table.  Each row of the output will
 *          contain the x,y,ra,dec positions for a standard star
 *          from the catalogue searched as well as any other information the
 *          input catalogue had to offer.  The cartesian coordinates are 
 *          transformed world coordinates using the image WCS and as such are 
 *          an indication of where the standard _should_ be on the image if 
 *          the WCS and the input equitorial coordinates are right.
 *      catsrc = char * (Given)
 *          This can be one of several values.  An internet connection is
 *          required for any starting with 'viz':
 *              localfits:
 *                  A locally held fits table will be searched. The full
 *                  path to this table must be given in 'catpath'.  The only
 *                  columns that the table must have are 'ra' and 'dec'.  Both
 *                  must be in degrees.
 *              viz2mass:
 *                  The 2MASS point source catalogue searched with VizieR
 *              vizlandolt:
 *                  The Landolt SA region standards catalogue search with VizieR
 *              vizusnob:
 *                  The USNO-B catalogue searched with Vizier.
 *              localusnob:
 *                  A local copy of the USNOB catalogue is searched in its 
 *                  native format. The path to the directory holding all the 
 *                  .cat and .acc files should be given in 'catpath'.
 *              local2mass:
 *                  A local copy of the 2MASS psc catalogue is searched in its 
 *                  native format. The path to the directory holding all the 
 *                  fits table files should be given in 'catpath'.
 *      site = char * (Given)
 *          The site for the VizieR search (e.g. catsrc = "vizusnob"). The copy
 *          of VizieR that is searched can be one of the following:
 *              casu:
 *                  CASU, Cambridge, UK
 *              cds:
 *                  CDS, Strasbourg, France
 *              ukirt:
 *                  UKIRT, Hawaii, USA
 *      catpath = char * (Given)
 *          The full path to a FITS table with a catalogue (catsrc=localfits)
 *          or the full path to a directory with local catalogue copies 
 *          (e.g. catsrc=vizusnob). 
 *      equinox = float (Given)
 *          The equinox that you want the standards transformed to.  This will
 *          only apply to catalogues that are not 'local'.
 *      fudge = int (Given)
 *          An integer value which specifies a percentage 'fudge' factor. This
 *          will expand the size of the field by the given percentage and is
 *          useful if you don't have much faith in your input WCS.
 *      cache = int (Given)
 *          If set, then this routine will check a cache area for a local
 *          table.  This may cut down the time required to do a large number
 *          of calls to the routine.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      None
 *              
 *  Dependencies:
 *      cfitsio, wcslib, cir_coverage.c, cir_wcssubs.c
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_getstds(char *infile, char *outfile, char *catsrc, char *site,
		       char *catpath, float equinox, int fudge, int cache,
		       char *errmsg) {
    int status,search,sock,racol,deccol,xcol,ycol,i,anynul,hdutype;
    int cache_status,retval;
    long nrows,naxis[2];
    double x,y;
    fitsfile *iptr,*tptr;
    char msg[BUFSIZ],siteaddr[BUFSIZ],catid[BUFSIZ],*req_string,key[16];
    char cache_cat[BUFSIZ],catpath2[BUFSIZ],catsrc2[BUFSIZ],cc[BUFSIZ];
    float ramin,ramax,decmin,decmax,ra,dec,dra,ddec,width,x1,x2,y1,y2;
    struct wcsprm *wcs;

    /* Find the RA and Dec limits of the image */

    if (cir_coverage(infile,fudge,&ramin,&ramax,&decmin,&decmax,msg) != CIR_OK) {
        sprintf(errmsg,"GETSTDS: Error getting WCS coverage for %s: %s\n",
                infile,msg);
        return(CIR_FATAL);
    }
    dra = ramax - ramin;
    ddec = decmax - decmin;
    ra = 0.5*(ramin + ramax);
    dec = 0.5*(decmin + decmax);

    /* Make copies of these into local variables. This is so that we can 
       change them if using the caching */

    strcpy(catpath2,catpath);
    strcpy(catsrc2,catsrc);
    if (!strcmp(catsrc,"localfits"))
        strcpy(cc,catpath);
    else 
	strcpy(cc,catsrc);

    /* If using a cache, then check to see if you have a file in the
       cache that will do for this call.  If you do, then try and
       open it.  If that works, then copy that name into the local
       catalogue string and do a local search. */

    cache_status = 1;
    if (cache) {
        cache_status = check_cache(cc,ramin,ramax,decmin,decmax,cache_cat); 
        status = 0;
        if (cache_status == 0) {
            (void)fits_open_file(&iptr,cache_cat,READONLY,&status);
            (void)fits_close_file(iptr,&status);
	}
        if (cache_status == 0 && status == 0) {
	    strcpy(catpath2,cache_cat);
	    strcpy(catsrc2,"localfits");
	}
    }


    /* Do I know about the sites and the sources? First check and see
       if a local catalogue (fits table) has been requested */

    if (! strcmp(catsrc2,"localfits")) {
        status = 0;
        (void)fits_open_file(&iptr,catpath2,READONLY,&status);
        if (status != 0) {
            fits_get_errstatus(status,msg);
            (void)sprintf(errmsg,"GETSTDS: Can't open local catalogue %s -- %s",
			  catpath2,msg);
            tidy();
            return(CIR_FATAL);
        }
        (void)fits_close_file(iptr,&status);
        search = LOCALSEARCH;
    
    /* Check to see if I know about the site that you've requested */

    } else if (! strcmp(catsrc2,"vizusnob") || ! strcmp(catsrc2,"viz2mass") ||
	       ! strcmp(catsrc2,"vizlandolt") || ! strcmp(catsrc2,"vizppmxl")) {
        if (! strcmp(site,"casu"))
            strcpy(siteaddr,CASU);
        else if (! strcmp(site,"ukirt"))
            strcpy(siteaddr,UKIRT);
        else if (! strcmp(site,"cds"))
            strcpy(siteaddr,CDS);
        else {
            (void)sprintf(errmsg,"GETSTDS: Don't recognise site: %s\n",site);
            tidy();
            return(CIR_FATAL);
        }
	if (! strcmp(catsrc2,"vizusnob"))
	    strcpy(catid,USNOB);
        else if (! strcmp(catsrc2,"viz2mass"))
	    strcpy(catid,TWOMASS);
        else if (! strcmp(catsrc2,"vizlandolt"))
	    strcpy(catid,LANDOLT);
        else if (! strcmp(catsrc2,"vizppmxl"))
	    strcpy(catid,PPMXL);
        search = INTERNETSEARCH;
    } else if (! strcmp(catsrc2,"localusnob")) {
	strcpy(catid,"usnob");
        search = CATSEARCH;
    } else if (! strcmp(catsrc2,"local2mass")) {
	strcpy(catid,"2mass");
	search = CATSEARCH;
    } else {
	(void)sprintf(errmsg,"GETSTDS: Unrecognised catalogue source: %s\n",
		      catsrc2);
	tidy();
	return(CIR_FATAL);
    }

    /* Ok, if this is an internet query, then do that now.  Start by forming
       the query parameters into a string and then encoding them */

    if (search == INTERNETSEARCH) {
        req_string = form_request(ra,dec,dra,ddec,equinox,catid);

        /* Now make the connection to the site */

        if (connect_site(siteaddr,&sock,msg) != CIR_OK) {
            (void)sprintf(errmsg,"GETSTDS: Couldn't connect to site: %s -- %s\n",
			  site,msg);
            tidy();
	    return(CIR_FATAL);
	}

        /* Send the request to the site */

        if (send_request(sock,req_string,msg) != CIR_OK) {
            (void)sprintf(errmsg,"GETSTDS: Couldn't send request to site: %s -- %s\n",
			  site,msg);
            tidy();
	    return(CIR_FATAL);
	}

        /* Get the response from the the server and write it to the output file */

        if (get_response(sock,outfile,msg) != CIR_OK) {       
            (void)sprintf(errmsg,"GETSTDS: Error receiving info from  site: %s -- %s\n",
			  site,msg);
            tidy();
	    return(CIR_FATAL);
	}

        /* Close up the connection */

        close(sock);

        /* Do a little surgery on the output file to change the RA, Dec column
           names */

        status = 0;
        (void)fits_open_file(&iptr,outfile,READWRITE,&status);
        (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
        if (status != 0) {
	    (void)sprintf(errmsg,"GETSTDS: No standards written to catalogue file\n");
	    closefits(iptr);
	    remove(outfile);
	    tidy();
	    return(CIR_WARN);
	}
        (void)fits_get_num_rows(iptr,&nrows,&status);
        if (status != 0 || nrows <= 0) {
	    (void)sprintf(errmsg,"GETSTDS: No rows written to catalogue file\n");
	    closefits(iptr);
	    remove(outfile);
	    tidy();
	    return(CIR_WARN);
	}     
	(void)fits_get_colnum(iptr,CASEINSEN,"RA*",&racol,&status);
	(void)sprintf(key,"TTYPE%d",racol);
	(void)fits_update_key(iptr,TSTRING,key,"ra",NULL,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"DE*",&deccol,&status);
	(void)sprintf(key,"TTYPE%d",deccol);
	(void)fits_update_key(iptr,TSTRING,key,"dec",NULL,&status);
        (void)fits_close_file(iptr,&status);
	if (status != 0) {
	    fits_get_errstatus(status,msg);
	    (void)sprintf(errmsg,"GETSTDS: Unable to update extracted catalogue %s -- %s\n",outfile,msg);
	    tidy();
	    return(CIR_FATAL);
	}
        if (cache)
            addto_cache(outfile,cc,ramin,ramax,decmin,decmax);
    
    /* If this is a local file search then do that now */

    } else if (search == LOCALSEARCH) {
        if (get_local(catpath2,ramin,ramax,decmin,decmax,outfile,msg) != CIR_OK) {
            (void)sprintf(errmsg,"GETSTDS: Unable to search local file %s -- %s\n",
			  catpath2,msg);
	    tidy();
	    return(CIR_FATAL);
	}
        (void)fits_open_file(&iptr,outfile,READWRITE,&status);
        (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
        if (status != 0) {
	    (void)sprintf(errmsg,"GETSTDS: No standards written to catalogue file\n");
	    closefits(iptr);
	    remove(outfile);
	    tidy();
	    return(CIR_WARN);
	}
        (void)fits_get_num_rows(iptr,&nrows,&status);
        if (status != 0 || nrows <= 0) {
	    (void)sprintf(errmsg,"GETSTDS: No rows written to catalogue file\n");
	    closefits(iptr);
	    remove(outfile);
	    tidy();
	    return(CIR_WARN);
	}     
        closefits(iptr);
        if (cache && cache_status)
            addto_cache(outfile,cc,ramin,ramax,decmin,decmax);
    } else if (search == CATSEARCH) {

        /* Extract the data from a local catalogue */

	width = 0.5*max(dra,ddec);
	status = extract_cat(catid,catpath2,ra,dec,width,&nrows,outfile,
			     msg);
	if (status != 0) {
	    (void)sprintf(errmsg,"GETSTDS: Error extracting local data: %s\n",
			  msg);
	    remove(outfile);
	    tidy();
	    return(CIR_FATAL);
	} else if (nrows <= 0) {
	    (void)sprintf(errmsg,"GETSTDS: No rows extracted from local catalogue\n");
	    remove(outfile);
	    tidy();
	    return(CIR_WARN);
	}
	if (cache)
	    addto_cache(outfile,cc,ramin,ramax,decmin,decmax);
    }

    /* Ok, we have a FITS table with the objects from the catalogue source.  Now
       we need to create columns with the XY positions of the objects according
       to the WCS in the input file. First open the WCS descriptor */

    retval = cir_wcsopen(infile,&wcs,msg);
    if (retval != CIR_OK) {
	(void)sprintf(errmsg,"GETSTDS: Unable to read WCS from %s -- %s\n",
		      infile,msg);
	tidy();
	return(CIR_FATAL);
    }

    /* Now open the extracted catalogue file */

    status = 0;
    (void)fits_open_file(&iptr,outfile,READWRITE,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"GETSTDS: Unable to open extracted catalogue %s -- %s\n",outfile,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Create two new columns for X and Y coordinates if they don't already
       exist */

    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"X_coordinate",&xcol,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"Y_coordinate*",&ycol,&status);
    if (status != 0) {
        status = 0;
        if (hdutype == ASCII_TBL) 
            (void)fits_insert_cols(iptr,1,2,ttype,tform1,&status);    
        else
            (void)fits_insert_cols(iptr,1,2,ttype,tform2,&status);    
        if (status != 0) {
	    fits_get_errstatus(status,msg);
    	    (void)sprintf(errmsg,"GETSTDS: Unable to add XY columns %s -- %s\n",outfile,msg);
            tidy();
            return(CIR_FATAL);
        }
        (void)fits_get_colnum(iptr,CASEINSEN,"X_coordinate",&xcol,&status);
        (void)fits_get_colnum(iptr,CASEINSEN,"Y_coordinate*",&ycol,&status);
    }

    /* Find the RA and Dec columns and the number of objects */

    (void)fits_get_colnum(iptr,CASEINSEN,"ra",&racol,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"dec",&deccol,&status);
    (void)fits_get_num_rows(iptr,&nrows,&status);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"GETSTDS: Unable to get column information %s -- %s\n",outfile,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Right, loop for each row. Calculate the XY position and write it out */

    for (i = 1; i <= nrows; i++) {
        (void)fits_read_col(iptr,TFLOAT,racol,i,1,1,NULL,&ra,&anynul,&status);
        (void)fits_read_col(iptr,TFLOAT,deccol,i,1,1,NULL,&dec,&anynul,&status);
        cir_radectoxy(wcs,ra,dec,&x,&y);
	ra = (float)x;
        dec = (float)y;
	(void)fits_write_col(iptr,TFLOAT,xcol,i,1,1,&ra,&status);
	(void)fits_write_col(iptr,TFLOAT,ycol,i,1,1,&dec,&status);
    }
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"GETSTDS: Unable to add XY data to %s -- %s\n",outfile,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Now select out any rows that are too far away from the edges. This helps
       correct searches that are near the poles */

    status = 0;
    (void)fits_open_file(&tptr,infile,READONLY,&status);
    (void)fits_get_img_size(tptr,2,naxis,&status);
    closefits(tptr);
    x1 = -0.01*(float)(fudge*naxis[0]);
    y1 = -0.01*(float)(fudge*naxis[1]);
    x2 = (float)naxis[0] - x1;
    y2 = (float)naxis[1] - y1;
    (void)sprintf(msg,
		  "X_coordinate > %g && X_coordinate < %g && Y_coordinate > %g && Y_coordinate < %g",
		  x1,x2,y1,y2);
    (void)fits_select_rows(iptr,iptr,msg,&status);

    /* OK, close everything up and get out of here */

    (void)fits_close_file(iptr,&status);
    cir_wcsclose(wcs);
    return(CIR_OK);
}

static int get_local(char *catpath, float ramin, float ramax, float decmin, 
		     float decmax, char *outfile, char *errmsg) {
    fitsfile *cptr,*optr;
    int status,i,wrap;
    float ra1[2],ra2[2];
    char buf[BUFMAX],msg[BUFSIZ];

    /* Start by creating a clone in the output file */

    status = 0;
    (void)fits_open_file(&cptr,catpath,READONLY,&status);
    (void)fits_create_file(&optr,outfile,&status);
    (void)fits_copy_hdu(cptr,optr,0,&status);
    (void)fits_movabs_hdu(cptr,2,NULL,&status);
    (void)fits_copy_header(cptr,optr,&status);
    (void)fits_modify_key_lng(optr,"NAXIS2",0,NULL,&status);
    (void)fits_close_file(optr,&status);
    (void)fits_open_file(&optr,outfile,READWRITE,&status);
    (void)fits_movabs_hdu(optr,2,NULL,&status);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"Unable to create new local file %s -- %s\n",outfile,msg);
        return(CIR_FATAL);
    }

    /* See if you have a wrap around problem at the equinox */

    if (ramin < 0.0 && ramax > 0.0) {
        wrap = 1;
        ra1[0] = ramin + 360;
	ra2[0] = 360.0;
        ra1[1] = 0.000001;
	ra2[1] = ramax;
    } else {
        wrap = 0;
        ra1[0] = ramin;
        ra2[0] = ramax;
    }

    /* Now do a select and pump the results into the output file */

    for (i = 0; i <= wrap; i++) {
        sprintf(buf,"ra >= %g && ra <= %g && dec >= %g && dec <= %g",
            ra1[i],ra2[i],decmin,decmax);
        (void)fits_select_rows(cptr,optr,buf,&status);
        if (status != 0) {
	    fits_get_errstatus(status,msg);
	    (void)sprintf(errmsg,"Unable to select from  file %s -- %s\n",outfile,msg);
            return(CIR_FATAL);
        }
    }

    /* Close up now */

    (void)fits_close_file(cptr,&status);
    (void)fits_close_file(optr,&status);
    return(CIR_OK);
}

static int get_response(int sock, char *outfile, char *errmsg) {
    int nnewline,i,nv,irem,fd,hdutype,status;
    long nrows;
    char buf[BUFMAX];
    fitsfile *iptr;

    /* Read from the socket until we find the double newline */

    nnewline = 0;
    while (1) {
        nv = recv(sock,buf,sizeof(buf),0);
        if (nv <= 0) {
	    (void)sprintf(errmsg,"Unable to find double newline\n");
	    return(CIR_FATAL);
	}
	for (i = 0; i < nv; i++) {
            if (buf[i] == '\n') {
		nnewline++;
                if (nnewline == 2) 
		    break;
	    } else {
		if (nnewline && buf[i] != '\r')
		    nnewline = 0;
	    }
	}
	if (nnewline == 2) {
	    irem = i + 1;
	    break;
	}
    }

    /* Shift the buffer */

    for (i = irem; i < nv; i++)
        buf[i-irem] = buf[i];
    nv -= irem;

    /* Open the output file */

    fd = open(outfile,O_WRONLY|O_CREAT,0644);

    /* Now write the data to the output file */

    while (1) {
        write(fd,buf,nv);
	nv = recv(sock,buf,sizeof(buf),0);
        if (nv == 0) 
	    break;
	else if (nv < 0) {
	    (void)sprintf(errmsg,"Read from socket failed\n");
	    close(fd);
	    remove(outfile);
	    return(CIR_FATAL);
	}
    }

    /* Close the output file. Check and see if a table was created and whether
       or not it contains any rows... */

    close(fd);
    status = 0;
    (void)fits_open_file(&iptr,outfile,READWRITE,&status);
    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    if (status != 0) {
	sprintf(errmsg,"GETSTDS: No standards were found\n");
	closefits(iptr);
	remove(outfile);
        return(CIR_FATAL);
    }
    (void)fits_get_num_rows(iptr,&nrows,&status);
    if (status != 0 || nrows <= 0) {
        sprintf(errmsg,"GETSTDS: Standards table contained no rows\n");
	closefits(iptr);
	remove(outfile);
	return(CIR_FATAL);
    }
    (void)fits_close_file(iptr,&status);
    return(CIR_OK);
}
	
static int send_request(int sock, char *req_string, char *errmsg) {
    char buf[BUFMAX];

    /* Form the full request contents and send it. */

    sprintf(buf,"GET /viz-bin/asu-fits?%s HTTP/1.0\r\n\r\n",req_string);
    if (send(sock,buf,strlen(buf),0) < 0) {
        (void)sprintf(errmsg,"Attemp to send message failed, error: %d\n",errno);
        return(CIR_FATAL);
    }
    return(CIR_OK);
}

static int connect_site(char *site, int *sock, char *errmsg) {
    struct hostent *hp;
    struct sockaddr_in sin;

    /* Get host pointer */

    hp = gethostbyname(site);
    if (hp == NULL) {
        sprintf(errmsg,"Unable to get host information for %s\n",site);
        return(CIR_FATAL);
    }
    sin.sin_family = hp->h_addrtype;
    memcpy(&sin.sin_addr,hp->h_addr_list[0],hp->h_length);
    sin.sin_port = htons(PORT);

    /*  Create an IP-family socket on which to make the connection */

    *sock = socket(hp->h_addrtype,SOCK_STREAM,0);
    if (*sock < 0) {
        sprintf(errmsg,"Unable to create socket descriptor for %s\n",site);
        return (CIR_FATAL);
    }

    /* Connect to that address...*/

    if (connect(*sock,(struct sockaddr*)&sin,sizeof(sin)) < 0) {
        sprintf(errmsg,"Unable to connect to site: %s\n",site);
        return(CIR_FATAL);
    }

    /* Otherwise send back a good status */

    return(CIR_OK);
}

static char *form_request(float ra, float dec, float dra, float ddec,
			  float equinox, char *catid) {
    static char buf[BUFSIZ];
    char buf2[BUFSIZ],equi[1];

    /* Format each of the fields and then add them to the buffer that will
       become the query string.  URLencode them before appending. First
       the coordinates */

    sprintf(buf2,"-c=%8.3f %8.3f",ra,dec);
    strcpy(buf,url_encode(buf2));

    /* Now the box size */

    sprintf(buf2,"-c.bd=%g/%g",dra,ddec);
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));

    /* The catalogue name */

    sprintf(buf2,"-source=%s",catid);
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));

    /* The equinox */

    equi[0] = (equinox == 1950.0 ? 'B' : 'J');
    sprintf(buf2,"-c.eq=%c%g",equi[0],equinox);
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));

    /* Finally make sure that we don't truncate the result and that we
       sort by RA */

    sprintf(buf2,"-out.max=unlimited");
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));
    sprintf(buf2,"-sort=_RA*-c.eq");
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));
    sprintf(buf2,"-out.add=_RA*-c.eq");
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));
    sprintf(buf2,"-out.add=_DEC*-c.eq");
    strcat(buf,"&");
    strcat(buf,url_encode(buf2));

    /* Right, get out of here */

    return(buf);
}


static char *url_encode(char *instring) {
    static char buf[BUFSIZ];
    int i,j,k,len;

    /* First copy over everything before the equal sign */

    k = 0;
    do {
        buf[k] = instring[k];
    } while (instring[k++] != '=');
    len = strlen(instring);

    /* Copy over the string, substituting the things that we don't want in a URL */

    j = k;
    for (i = k; i < len; i++) {
        assert(j < sizeof(buf));
        if (instring[i] == ' ') 
            buf[j++] = '+';
        else if (isalnum((unsigned char)instring[i]))
            buf[j++] = instring[i];
        else {
            sprintf(buf+j,"%%%2x",(unsigned char)instring[i]);
            j += 3;
        }
    }
    buf[j] = '\0';
    return(buf);
}

static int check_cache(char *src, float ra1_im, float ra2_im, float dec1_im, 
		       float dec2_im, char *cache_cat) {
    int wrap1,wrap2;
    FILE *fd;
    char catname[BUFSIZ],csrc[BUFSIZ];
    float best,ra1_cat,ra2_cat,dec1_cat,dec2_cat,d1,d2,fra,fdec,ftot;

    /* Check to see if there is wrap around in the coordinates */

    wrap1 = (ra1_im < 0.0 ? 1 : 0);

    /* Open the index file.  NB the path and file name are hardcoded */

    fd = fopen(CACHEIND,"r");
    if (fd == NULL)
        return(1);

    /* Now see if you have any matching entries */

    best = 0.0;
    while (fscanf(fd,"%s %s %g %g %g %g",catname,csrc,&ra1_cat,&ra2_cat,
		  &dec1_cat,&dec2_cat) != EOF) {
        if (strcmp(csrc,src))
            continue;
        wrap2 = (ra1_cat < 0.0 ? 1 : 0);
        if (wrap1 != wrap2)
	    continue;
	
        /* Check to see if there is at least some overlap */

        if (!(((ra1_im >= ra1_cat && ra1_im <= ra2_cat) ||
	     (ra2_im >= ra1_cat && ra2_im <= ra2_cat)) && 
	    ((dec1_im >= dec1_cat && dec1_im <= dec2_cat) ||
	     (dec2_im >= dec1_cat && dec2_im <= dec2_cat))))
	    continue;

	/* Work out exactly how much there is in each coordinate */

        d1 = max(0.0,ra1_cat-ra1_im);
        d2 = max(0.0,ra2_im-ra2_cat);
	fra = 1.0 - (d1 + d2)/(ra2_im - ra1_im);
        d1 = max(0.0,dec1_cat-dec1_im);
        d2 = max(0.0,dec2_im-dec2_cat);
	fdec = 1.0 - (d1 + d2)/(dec2_im - dec1_im);
	ftot = fra*fdec;

        /* Keep track of which is the best one */

        if (ftot > best) {
	    sprintf(cache_cat,"%s/%s",CACHEDIR,catname);
	    best = ftot;
        }
    }
    fclose(fd);

    /* Return good status if there is sufficient overlap */

    if (best > 0.9) 
	return(0);
    else
	return(1);
}
    
static void addto_cache(char *fname, char *src, float ramin, float ramax, 
			float decmin, float decmax) {
    FILE *fd,*fd2;
    char newname[BUFSIZ],new2[BUFSIZ];
    int ch,i;

    /* Check to see if the cache directory exists.  If it doesn't, then create
       it. */

    if (access(CACHEDIR,0) != 0)
        mkdir(CACHEDIR,0775);

    /* Open the index file with 'append' access */

    fd = fopen(CACHEIND,"a");

    /* Check the files in the directory to get a number that isn't already
       being used */

    i = 0;
    while (1) {
        i++;
        sprintf(newname,"%s/cch_%08d",CACHEDIR,i);
        if (access(newname,F_OK) != 0) 
	    break;
    }

    /* Now write the current entry and make a copy of the file into the 
       directory */

    sprintf(newname,"cch_%08d",i);
    fprintf(fd,"%s %s %g %g %g %g\n",newname,src,ramin-0.0005,ramax+0.0005,
	    decmin-0.0005,decmax+0.0005);
    fclose(fd);
    sprintf(new2,"%s/%s",CACHEDIR,newname);
    fd = fopen(fname,"r");
    fd2 = fopen(new2,"wb");
    while ((ch = getc(fd)) != EOF)
	putc(ch,fd2);
    fclose(fd);
    fclose(fd2);
}


static int cir_coverage(char *infile, int fudge, float *ra1, float *ra2, 
                        float *dec1, float *dec2, char *errmsg) {
    struct wcsprm *wcs;
    float dra,ddec,boxfudge,min_4q,max_1q;
    double ra,dec,x,y;
    int nx,ny,i,j,first_quad,fourth_quad,retval,status;
    char msg[BUFSIZ];
    fitsfile *iptr = NULL;
    long naxes[2];

    /* Open the file and grab its WCS info */

    retval =  cir_wcsopen(infile,&wcs,msg);
    if (retval != CIR_OK) {
        sprintf(errmsg,"COVERAGE: Can't get wcs info from: %s -- %s\n",infile,
		msg);
        return(CIR_FATAL);
    }

    /* Get the size of the data array */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_img_size(iptr,2,naxes,&status);
    (void)fits_close_file(iptr,&status);
    if (status != 0) {
	(void)fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"COVERAGE: Can't open image %s -- %s\n",infile,
		      msg);
	freewcs(wcs);
	return(CIR_FATAL);
    }

    /* Find the RA and Dec limits of the image */

    nx = (int)naxes[0];
    ny = (int)naxes[1];
    *ra1 = 370.0;
    *ra2 = -370.0;
    *dec1 = 95.0;
    *dec2 = -95.0;
    first_quad = 0;
    fourth_quad = 0;
    min_4q = 370.0;
    max_1q = 0.0;
    for (j = 1; j < ny; j += 10) {
        j = min(j,ny);
        y = (double)j;
        for (i = 1; i < nx; i += 10) {
            i = min(i,nx);
            x = (double)i;
            cir_xytoradec(wcs,x,y,&ra,&dec);
            if (ra >= 0.0 && ra <= 90.0) {
                first_quad = 1;
                max_1q = max((float)ra,max_1q);
            } else if (ra >= 270.0 && ra <= 360.0) {
                fourth_quad = 1;
                min_4q = min((float)(ra-360.0),min_4q);
            }
            *ra1 = min(*ra1,(float)ra);
            *ra2 = max(*ra2,(float)ra);
            *dec1 = min(*dec1,(float)dec);
            *dec2 = max(*dec2,(float)dec);
        }
    }
    cir_wcsclose(wcs);

    /* Now have a look to see if you had RA values in both the first and
       fourth quadrants.  If you have, then make the minimum RA a negative
       value.  This will be the signal to the caller that you have the
       wraparound... */
 
    if (first_quad && fourth_quad) {
        *ra1 = min_4q;
        *ra2 = max_1q;
    }

    /* Pad out search a bit */

    if (fudge) {
        boxfudge = 0.01*(float)fudge;
        dra = 0.5*boxfudge*(*ra2 - *ra1);
        *ra1 -= dra;
        *ra2 += dra;
        ddec = 0.5*boxfudge*(*dec2 - *dec1);
        *dec1 -= ddec;
        *dec2 += ddec;
    }

    /* Exit */

    return(CIR_OK);
}



static void tidy () {
    return;
}

/*

$Log: cir_getstds.c,v $
Revision 1.3  2012/08/13 10:48:12  jim
Added PPMXL to list of astrometric catalogues

Revision 1.2  2011-09-16 12:05:00  jim
Added code to weed out entries that fall outside the fudged box. This is
important when looking near the poles

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.25  2008/09/15 08:06:22  jim
Updated url for CDS

Revision 1.24  2008/03/07 11:50:53  jim
Vizier lookups now target columns starting with RA and DE rather than _RA and
_DE

Revision 1.23  2005/04/13 12:34:06  jim
Changed naming scheme for cached files

Revision 1.22  2004/10/22 10:10:34  jim
yet another bug

Revision 1.21  2004/10/22 09:36:31  jim
small bug fix

Revision 1.20  2004/10/22 08:22:57  jim
fixed little bug where vizier catalogues weren't being trapped correctly

Revision 1.19  2004/10/07 08:59:38  jim
Fixed small bug in main routine

Revision 1.18  2004/09/07 14:18:52  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.17  2004/08/31 09:14:28  jim
Added support for local catalogues

Revision 1.16  2004/08/02 11:50:01  jim
Modified to use new version of cir_wcssubs routines

Revision 1.15  2004/03/04 09:54:26  jim
Changed format for XY columns to cope with larger numbers

Revision 1.14  2004/02/05 09:48:29  jim
Fixed bug where CFITSIO status variable wasn't initialised before the first
call in the cache testing section.  Caused a seg. violation with some
compilers

Revision 1.13  2004/02/03 11:33:53  jim
Modifed so that caching now logs the source of the query. This is done in
case several different sources of standards are used during a run

Revision 1.12  2004/01/30 14:57:59  jim
Slight change to get around stupid tempnam behaviour

Revision 1.11  2004/01/30 14:26:28  jim
Added cache parameter

Revision 1.10  2003/11/13 11:29:15  jim
Now checks to see if any rows were written.  Also fixed a small bug in
the querying part whereby if the ra or the dec were specified without a
decimal point, the VizieR parser gagged.

Revision 1.9  2003/09/30 20:32:01  jim
Modified so that if the X,Y columns already exist in the output table, then
it isn't an error and it just overwrites the values in that column

Revision 1.8  2003/09/11 09:55:56  jim
Minor modification to some comments

Revision 1.7  2003/09/01 10:35:03  jim
teensy bug fix

Revision 1.6  2003/09/01 08:54:44  jim
Nearly complete rewrite.  Now uses various VizieR sites as the main internet
source.

Revision 1.5  2003/06/03 07:06:22  jim
Changed the way that ptrans detects the beginning of the info in response
to new version of apache.  This whole routine needs rewriting...

Revision 1.4  2003/01/29 15:43:54  jim
Modified to deal with the wraparound problem at the equinox

Revision 1.3  2003/01/07 10:32:16  jim
Fixed a few typos

Revision 1.2  2002/12/16 10:23:34  jim
Plugged in routine cir_coverage in the relevant location.  Added prologue

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/
