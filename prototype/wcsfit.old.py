import os
import sys
import subprocess
import shutil
import argparse
import logging
from time import sleep
from distutils import spawn

import numpy
from astropy.io import fits
from astropy.wcs import WCS
import pywcs

# Find 'imcore' and 'wcsfit' 
IMCORE=os.getenv('IMCORE', spawn.find_executable('imcore'))
WCSFIT=os.getenv('WCSFIT', spawn.find_executable('wcsfit'))

logging.basicConfig(level=logging.INFO)

def checkExeInPath():
    logging.info('IMCORE: %s', IMCORE)
    logging.info('WCSFIT: %s', WCSFIT)
    if (IMCORE is None) or (not os.access(IMCORE, os.X_OK)):
        logging.error('imcore not found in path')
        sys.exit(1)
    if (WCSFIT is None) or (not os.access(WCSFIT, os.X_OK)):
        logging.error('wcsfit not found in path')
        sys.exit(1)

def getNumberOfExtensions(input):
    """Return number of extensions in a FITS file"""
    fh = fits.open(input)
    n = len(fh) - 1
    fh.close()
    return n

def mergeQuadrants(input):
    """Merge quadrants into CCD image."""
    import pyfits
    
    output = input.replace('.fit','_ccd.fit')
    
    if os.access(output, os.F_OK):
        return
    
    logging.info('merging quadrants into CCD image (%s)' % input)
    
    fh = pyfits.open(input)
    hdulist=[fh[0]]

    for i in range(36):
            logging.debug('CCD #%d' % (i+1))
            ccd=pyfits.ImageHDU(data=numpy.zeros([4132,4096]).astype(numpy.float32), header=fh[i*4+1].header)
            ccd.data[0:2066,0:2048]=fh[4*i+1].data.astype(numpy.float32)
            ccd.data[0:2066,2048:4096]=fh[4*i+2].data.astype(numpy.float32)
            ccd.data[2066:4132,0:2048]=fh[4*i+3].data.astype(numpy.float32)
            ccd.data[2066:4132,2048:4096]=fh[4*i+4].data.astype(numpy.float32)
            hdulist+=[ccd]
        
    nfh=pyfits.HDUList(hdulist)
    
    logging.info('writing CCD image to disk (%s)' % output)
    
    nfh.writeto(output)

def initAstrometry(input, moreKeys=False):
    """Initialise the astrometric keywords in the headers and allocate space in headers for more keywords."""
    
    logging.info('writing initial astrometric solution')
    
    fh = fits.open(input, mode='update')

    for i in range(len(fh)-1):
        w = WCS(fh[i+1].header)
        crpix = numpy.array([fh[i+1].header['NAXIS1']/2.0, fh[i+1].header['NAXIS2']/2.0], numpy.float_)
        crval = w.wcs_pix2world([crpix],1)[0]
        fh[i+1].header['CRVAL1']=crval[0]
        fh[i+1].header['CRVAL2']=crval[1]
        fh[i+1].header['CRPIX1']=crpix[0]
        fh[i+1].header['CRPIX2']=crpix[1]
    
        if 'PC1_1' in fh[i+1].header: del fh[i+1].header['PC1_1']
        if 'PC2_2' in fh[i+1].header: del fh[i+1].header['PC2_2']
    
        # Allocate space in header
        if moreKeys:
            for j in range(10):
                fh[i+1].header.add_blank()
            
    fh.flush() 
    fh.close()


def run_imcore(input, confmap='noconf', cattype=6):
    """Extract objects from image using imcore"""
    
    logging.info('detecting objects on images')
    
    catalogue = input.replace('.fit', '_cat.fit')
    command = [IMCORE, input, confmap, catalogue, '3', '1.5', '--noell', '--nocrowd', '--cattype=%s' % cattype]
    #print ' '.join(command)
    res = subprocess.call(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    # Clean up
    if not res: shutil.rmtree('catcache', True)
    return res

def weightToConf(input, confmap):
    """Convert weight map to 100 median and integer."""
    import pyfits
    
    logging.info('making confidence map')
    
    output = input.replace('.fit','_conf.fit')
    fh = pyfits.open(confmap)
    hdulist=[fh[0]]

    for i in range(36):
        weight = fh[i+1].data * 100.0
        ccd=pyfits.ImageHDU(data=weight.astype(numpy.int16))
        hdulist+=[ccd]
        
    nfh=pyfits.HDUList(hdulist)
    
    logging.info('writing confidence image to disk (%s)' % output)
    
    nfh.writeto(output, clobber=True)
    
    return output
    
def run_wcsfit(input, refcat, extno=0, async=False):
    
    logging.info('fitting astrometry')
    
    # Read position from header
    fh = fits.open(input)
    # Workaround until RA, DEC are in primary header
    crval1, crval2 = fh[1].header['CRVAL1'], fh[1].header['CRVAL2']
    fh.close()
    catalogue = input.replace('.fit', '_cat.fit')
    command = [WCSFIT, input, catalogue, '--catsrc=localfits', '--catpath=%s' % refcat]
    if extno>0:
        command += ['--extno=%d' % extno]
    if async:
        res = subprocess.Popen(command)
    else:
        res = subprocess.call(command)
        if not res: shutil.rmtree('catcache', True)
    return res

def print_summary(input):
    fh = fits.open(input)
    for i in range(1, len(fh)):
        logging.info('%s[%d] %d %f' % (input, i, fh[i].header.get('NUMBRMS', 0), fh[i].header.get('STDCRMS', -1)))
    fh.close()

def usage():
    print "Usage: wcsfit.py [-q] image.fits"
    print ""
    print " -q : given an input image with quadrants do not merge into CCDs and fit the quadrants individually"
    print ""
    sys.exit(0)

if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Perform astrometric calibration.')
    parser.add_argument('-q', action='store_true', help='process the image without merging quadrants')
    parser.add_argument('--async', action='store_true', help='perform astrometric calibration asynchronously')
    parser.add_argument('--copy', action='store_true', help='copy input file before modification')
    parser.add_argument('--nocat', action='store_true', help='skip catalogue generation')
    parser.add_argument('--cattype', default='6', help='catalogue type to generate (3=basic, 6=long)')
    parser.add_argument('--weight', default='noconf', help='weight map')
    parser.add_argument('image', help='FITS image')
    parser.add_argument('refcat', help='Reference catalogue')
    args = parser.parse_args()
    
    if args.copy:
        shutil.copy(args.image, 'w%s' % args.image)
        input = 'w%s' % args.image
    else:
        input = args.image

    refcat = args.refcat
        
    if args.q:
        fitQuadrants = True
    else:
        fitQuadrants = False

    # Check that executables are in path
    checkExeInPath()

    # Check number of extensions and merge quadrants into ccd image if necessary
    nExtensions = getNumberOfExtensions(input)
    if nExtensions==144: # Quadrants
        if not fitQuadrants:
            mergeQuadrants(input)   
            input = input.replace('.fit', '_ccd.fit')
    elif nExtensions==36:
        pass
    else:
        logging.error('fits file does not contain correct number of extensions')
        sys.exit(1)

    if not args.nocat:
        initAstrometry(input, moreKeys=True)
        logging.info('running astrometric calibration on %s' % input)
        if args.weight.find('noconf')==-1:
            confmap = weightToConf(input, args.weight)
        else:
            confmap = args.weight
        if run_imcore(input, confmap=confmap, cattype = args.cattype):
            logging.error('imcore failed')
            sys.exit(1)
        else:
            logging.info('imcore success')
        
    # We may run the fitting individually for each chip [BETA TEST]
    # Better to move this block to another function e.g. runWcsfitAsync
    if args.async:
        a = []
        for i in range(nExtensions):
            a.append(run_wcsfit(input, refcat, extno=i+1, async=True))
        p = map(lambda x: x.poll(), a)
        while p.count(None)>0:
            j = p.index(None)
            a[j].wait()
            p = map(lambda x: x.poll(), a)
        if p.count(0)!=nExtensions:
            logging.error('wcsfit failed')
        else:
            logging.info('wcsfit success')
    else:
        if run_wcsfit(input, refcat):
            logging.error('wcsfit failed')
        else:
            logging.info('wcsfit success')

    print_summary(input)
