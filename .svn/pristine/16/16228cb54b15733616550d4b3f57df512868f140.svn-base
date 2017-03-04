#!/usr/bin/python

"""
Implement astrometric fit to Euclid OU-VIS Images.

Given an image and a reference catalogue, extract objects from the image
and solve the astrometric distortion.

Author(s)
----------

2016 Eduardo Gonzalez-Solares

Changelog
----------

TODO

 - Better handling of cases with no stars or not enough stars
 - Star/Galaxy separation
 - We may want to let the software decide if it should fit the PV coefficients
   base e.g. on the number of stars

201611

 - Star/Galaxy separation (partially done for imcore)
 - Compute offset between catalogues now using FFT
 - When fitting affine transformation provide estimation on success before 
   updating a solution
 - When fitting PV provide estimation on success before updating a solution
 - Added experimental Python based source extractor.

201610

 - Nearest neighbour distribution used to compute optimum match radius.
 - Polynomial of distortion defined to have 12 elements per axis (3 degrees). 
   Added option to change it when initializing.
 - Added check of astrometry RMS in order to end iteration.
 - Removed need for wcsfit from casutools. Dependency on imcore still present.
   Probably needs some bug squishing for cases with no stars in detector, etc.
 - Reworked most of the run_wcsfit function.
 - Allow to fit full astrometric model for dense fields removing the need of
   using SCAMP.
 - Added CCDImage class to drive the fitting
 - Dropped support for running in quadrants only.
 - Dropped multithreading support.
 
"""

__version__ = 201611

# Module imports
import os
import sys
import subprocess
import shutil
import argparse
import logging
from datetime import datetime
from distutils import spawn
from configobj import ConfigObj

# Numeric Python and AstroPy modules
import random
import numpy as np
from scipy import optimize
from scipy import interpolate
import scipy.ndimage as sn
import scipy.fftpack as fft

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

# Median Absolute Deviation
mad = lambda x: 1.4*np.median(np.absolute(x - np.median(x)))

# Define local status variables
OK = 0
ERROR =1
STATUS = OK

# Find 'imcore'
# It first checks whether a IMCORE enviroment variables has been set
# and defaults to whatever version is in path if not
IMCORE=os.getenv('IMCORE', spawn.find_executable('imcore'))

# Set logging
logger = logging.getLogger('wcsfit')
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(name)s:%(levelname)s:%(message)s'))
logger.addHandler(console)
logger.setLevel(logging.INFO)

def check_wcslib():
    """
    Try to determine wcslib version.
    """
    try:
        import astropy.wcs
        version = astropy.wcs._wcs.__version__
        version = tuple(map(int, version.split('.')))
    except:
        logger.warn('Cannot determine wcslib version.')
        return OK

    if version<(5,9):
        logger.error('wcslib version too old {}'.format(version))
        return ERROR

    return OK

def check_exe_in_path():
    """
    Check that imcore executable is available. It must be given 
    either in the IMCORE environment variable or present in PATH.
    """
    logger.info('IMCORE: %s', IMCORE)
    if (IMCORE is None) or (not os.access(IMCORE, os.X_OK)):
        logger.error('imcore not found in path')
        return ERROR
        
def get_number_of_extensions(input):
    """
    Return number of extensions in a FITS file
    """
    fh = fits.open(input)
    n = len(fh) - 1
    fh.close()
    return n
    
def extract_extension(input, output, nhdu, **kwargs):
    """
    Extract a extension from a MEF a write to output.
    """

    fh = fits.open(input)
    fh[nhdu].writeto(output, **kwargs)
    fh.close()


def merge_quadrants(input, append='_ccd', overwrite=False):
    """
    Given an input image containing quadrants (image.fit), produce
    an image containing CCDs (image_ccd.fit)
    """
    
    logger.info('merging quadrants into CCD image (%s)' % input)
    
    output = input.replace('.fit','{}.fit'.format(append))
    
    if os.access(output, os.F_OK):
        logger.info('output image already exists {}'.format(output))
        if not overwrite: return
    
    fh = fits.open(input)
    hdulist=[fh[0]]

    logger.debug('input image has {} extensions'.format(len(fh)-1))

    for i in range(int((len(fh)-1)/4)):
            logger.debug('merging CCD #%d' % (i+1))
            ccd=fits.ImageHDU(data=np.zeros([4132,4096]).astype(np.float32), 
                              header=fh[i*4+1].header)
            ccd.data[0:2066,0:2048]=fh[4*i+1].data.astype(np.float32)
            ccd.data[0:2066,2048:4096]=fh[4*i+2].data.astype(np.float32)
            ccd.data[2066:4132,0:2048]=fh[4*i+3].data.astype(np.float32)
            ccd.data[2066:4132,2048:4096]=fh[4*i+4].data.astype(np.float32)
            hdulist+=[ccd]
        
    nfh=fits.HDUList(hdulist)
    
    logger.info('writing CCD image to disk (%s)' % output)
    
    nfh.writeto(output, clobber=True)
    
def weight_to_conf(input, confmap):
    """
    Convert weight map to 100 median and integer as required by IMCORE.
    """
    
    logger.info('making confidence map')
    
    output = input.replace('.fit','_conf.fit')
    fh = fits.open(confmap)
    hdulist=[fh[0]]

    for i in range(36):
        weight = fh[i+1].data * 100.0
        ccd=fits.ImageHDU(data=weight.astype(np.int16))
        hdulist+=[ccd]
        
    nfh=fits.HDUList(hdulist)
    
    logger.info('writing confidence image to disk (%s)' % output)
    
    nfh.writeto(output, clobber=True)
    
    return output
    
def imcore(image_name, confmap='noconf'):
    """
    Extract objects from image using imcore.

    TODO: Possibly replace this by some other method (?). Need to take
    into account PSF shape for centroiding.

    Parameters
    ----------

    image_name: name of input image
    
    confmap: name of confidence map, either 'noconf' to not use a map,
             'auto' to look for one called image_conf.fit or name.

    cattype: type of catalogue to generate. Just leave the default.
    """

    logger.info('detecting objects on images')

    check_exe_in_path()
    
    catalogue = image_name.replace('.fit', '_cat.fit')
    if confmap=='auto':
        confmap = image_name.replace('.fit', '_conf.fit')
        if not os.access(confmap, os.R_OK):
            confmap='noconf'
        
    command = [IMCORE, image_name, confmap, catalogue, 
                  '3', '1.5', '--noell', '--nocrowd', '--cattype=6']
    res = subprocess.call(command, stderr = subprocess.PIPE, 
                                   stdout = subprocess.PIPE)
    # Clean up
    if not res: shutil.rmtree('catcache', True)
    return res
        
class CCDImage(object):
    """
    Class to deal with astrometry of Euclid CCD images.
    """
    def __init__(self, image_name, weight_map=None, mode='readonly', save_backup=False, npv=12):
        """
        Initialize the class. Note that errors are not fully handled, i.e., 
        it is the user responsability to feed a valid FITS file.

        Parameters
        -----------

        image_name: string
            name of MEF image

        weight_map: string, optional
            name of MEF weight map. Set to None for no weight map

        mode: string, optional
            readonly or update

        save_backup: bool, optional
            if True and mode='update' save a backup before updating

        npv: integer, optional
            number of distortion elements to fit 
        """
        self.image_name = image_name
        self.weight_map = weight_map

        # Number of PV coefficients in each axis
        self.npv = npv

        self.mode = mode
        
        self.cattype = 'imcore'
        
        self.use_leastsq = True
        
        self.fh = fits.open(self.image_name, mode=self.mode, save_backup=save_backup)

    def fit_astrometry(self, refcat=None,
                             wcsconf='wcsfit.conf',
                             fitpv=False):
        """
        Covenience function to run all commands to fit the astrometry to an 
        image.

        Parameters
        ----------

        refcat: string
             name of the catalogue containing reference stars

        wcsconf: string
            file containing astrometric model parameters

        fitpv: bool
            fit the model distortion coefficients
        """
        self.init_astrometry()
        self.add_astrometric_distortion(wcsconf)
        self.wcsfit(refcat=refcat, fitpv=fitpv)

    def init_astrometry(self, centre='ccd', more_keys=False):
        """
        Initialise the astrometric keywords in the headers and optionally allocate space for more keywords.
        This function only moves the reference point to the centre of each CCD. See function
        add_astrometric_distortion to add the astrometric model.
        
        Parameters
        ----------
        
        centre: location of the reference point. Either in the centre of each CCD or in the centre of the FPA.
        """
    
        logger.info('writing initial astrometric solution ({})'.format(centre))
    
        for i in range(len(self.fh)-1):
            w = WCS(self.fh[i+1].header)
            if centre == 'fpa':
                # The centre of the FPA is in the primary HDU
                cenra, cendec = self.fh[0].header['RA'], self.fh[0].header['DEC']
                crval = np.array([cenra, cendec], np.float_)
                crpix = w.wcs_world2pix([crval],1)[0]
            elif centre == 'ccd':
                # We set the reference point to the centre of the CCD
                crpix = np.array([self.fh[i+1].header['NAXIS1']/2.0, self.fh[i+1].header['NAXIS2']/2.0], np.float_)
                crval = w.wcs_pix2world([crpix],1)[0]
            
            self.fh[i+1].header['CRVAL1']=crval[0]
            self.fh[i+1].header['CRVAL2']=crval[1]
            self.fh[i+1].header['CRPIX1']=crpix[0]
            self.fh[i+1].header['CRPIX2']=crpix[1]

            self.fh[i+1].header['CTYPE1']='RA---TPV'
            self.fh[i+1].header['CTYPE2']='DEC--TPV'

            # Initialize PV coefficients
            # Two loops so that they are ordered in the header
            for k in range(self.npv):
                self.fh[i+1].header['PV1_{}'.format(k)] = 0.0
            for k in range(self.npv):
                self.fh[i+1].header['PV2_{}'.format(k)] = 0.0

            self.fh[i+1].header['PV1_1'] = 1.0
            self.fh[i+1].header['PV2_1'] = 1.0
        
            # Delete PC keywords since we are using the CD notation
            if 'PC1_1' in self.fh[i+1].header: del self.fh[i+1].header['PC1_1']
            if 'PC2_2' in self.fh[i+1].header: del self.fh[i+1].header['PC2_2']
            if 'PC2_1' in self.fh[i+1].header: del self.fh[i+1].header['PC2_1']
            if 'PC1_2' in self.fh[i+1].header: del self.fh[i+1].header['PC1_2']
    
            # Allocate space in header
            if more_keys:
                for j in range(20):
                    self.fh[i+1].header.add_blank()
            
    def add_astrometric_distortion(self, config):
        """
        Add astrometric distortion model to the header. This function reads the 
        astrometric model (i.e. PV keywords) from a configuration file and updates
        the header accordingly.
        """
    
        logger.info('adding distortion coefficients')
    
        # Open image to update
        #fh = fits.open(self.image_name, mode='update')

        # Open the file containing coefficients
        logger.info('Distortion model {}'.format(config))
        conf = ConfigObj(config)

        # Add coefficents for each extension
        for i in range(len(self.fh)-1):
            self.fh[i+1].header['CTYPE1'] = 'RA---TPV'
            self.fh[i+1].header['CTYPE2'] = 'DEC--TPV'
    
            self.fh[i+1].header['CUNIT1'] = 'deg'
            self.fh[i+1].header['CUNIT2'] = 'deg'

            pvs = conf['CCD{}'.format(self.fh[i+1].header['CCDID'].strip())]
            for k,v in pvs.items():
                self.fh[i+1].header[k] = float(v)
        
    def find_stars(self, catalogue_name=None, cattype='imcore', extension=1):
        """
        Find stars on image.

        Parameters
        -----------

        catalogue_name: string, optional
            file name of extracted catalogue

        cattype: string, optional
            type of catalogue (imcore|sextractor|python)

        extension: integer, optional
            catalogue extension to read
        """

        self.cattype=cattype

        if self.cattype=='imcore':
            logger.info('Using imcore object extractor')
            if not catalogue_name:
                catalogue_name = self.image_name.replace('.fit', '_cat.fit')
            if not os.access(catalogue_name, os.R_OK):
                imcore(self.image_name, confmap='auto')
                
            cat = fits.open(catalogue_name)
            xcat = cat[extension].data.field('x_coordinate')
            ycat = cat[extension].data.field('y_coordinate')
            
            # We read in classification and ellipticity 
            # and select only point-like objects
            classification = cat[extension].data.field('classification')
            ellipticity = cat[extension].data.field('ellipticity')
            mask = (classification == -1) & (ellipticity < 0.1)
            
            xcat, ycat = xcat[mask], ycat[mask]
            cat.close()
        elif self.cattype=='sextractor':
            logger.error('Sextractor catalogue type not yet implemented')
            xcat = ycat = 0.0
        elif self.cattype=='python':
            logger.info('Using Python object extractor')
            import extractor
            t = extractor.extract(self.fh[extension].data)
            xcat, ycat = t['x'], t['y']

        return (xcat, ycat)
        
    def __fill_coords(self):
        """
        Calculates RA and DEC coordinates from X and Y and stores in catalogue.
        This step is necessary if the RA, DEC coordinates are needed after a run
        of imcore.

        Currently this is a convenience function if one needs to work with the
        catalogue generated by imcore elsewhere. It is not used in the fitting
        procedure.
        """
    
        logger.info('filling RA & Dec coordinates')
    
        # Open image and catalogue
        img = fits.open(self.image_name)
        cat = fits.open(self.image_name.replace('.fit', '_cat.fit'), 
                                                           mode='update')
    
        # Use the header on the image to compute coordinates and update catalogue
        for i in range(len(img)-1):
            w=WCS(img[i+1].header)
            x=cat[i+1].data.field('x_coordinate')
            y=cat[i+1].data.field('y_coordinate')
            a,d = w.all_pix2world(x,y,1)
            cat[i+1].data['ra']=a
            cat[i+1].data['dec']=d
            cat.flush()
        
        cat.close()
        img.close()
        
    def _errfunc(self, p, x, y, xi, eta):
        """
        Function to be minimized during least squares fit to calculate the 
        affine transformation between two sets of points.

        Parameters
        ----------

        p: vector to be mimimized containing rotation, scale and translation
        x, y: image coordinates of objects detected in image
        xi, eta: standard coordinates of reference objects
        """
        r, s, tx, ty   = p
    
        delta_xi = xi - ( -r*x + s*y + tx )
        delta_eta = eta - (s*x + r*y + ty )
        res = np.concatenate((delta_xi*3600.0, delta_eta*3600.0))
    
        return res.flatten()
        
    def _tpv(self, pv, x, y, hdr, refra, refdec):
        """
        Function to be minimized when computing the full distortion model.
        """
        pv = pv.reshape(2,self.npv)
        
        tmphdr = hdr.copy()
        for k in range(self.npv):
            tmphdr['PV1_{}'.format(k)] = pv[0][k]
            tmphdr['PV2_{}'.format(k)] = pv[1][k]
        w=WCS(tmphdr)

        xip, etap = w.wcs_pix2world(x, y, 1)
        xip = xip - hdr['CRVAL1']
        etap = etap - hdr['CRVAL2']
        
        delta_xi = xip - (refra - hdr['CRVAL1'])
        delta_eta = etap - (refdec - hdr['CRVAL2'])
        res = np.concatenate((delta_xi*3600., delta_eta*3600.0))
        
        return res.flatten()
        
        
    def _xmatch(self, refra, refdec, ra, dec):
        """
        Compute cross match between two sets of coordinates.
        """
        c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        refc = SkyCoord(ra=refra*u.degree, dec=refdec*u.degree)
        idx, d2d, d3d = refc.match_to_catalog_sky(c)
        return idx, d2d.arcsec

    def _find_shift_crosscorr(self, refx, refy, x, y, scl=6):
        """
        Find the shift between two coordinates by FFT crosscorrelation.
        """
        dimx = int(max(refx.max(), x.max())/scl)
        dimy = int(max(refy.max(), y.max())/scl)
        im1 = np.zeros([dimx+1, dimy+1])
        im2 = np.zeros([dimx+1, dimy+1])
        im1[map(int, refx/scl), map(int, refy/scl)]=1
        im2[map(int, x/scl), map(int, y/scl)]=1
        imag1 = sn.filters.gaussian_filter(im1,1.0)
        imag2 = sn.filters.gaussian_filter(im2,1.0)
        im1 = fft.fft2(imag1)
        im2 = fft.fft2(imag2)
        crosscorr = fft.ifft2(im1*np.ma.conjugate(im2))
        modulus = fft.fftshift(np.abs(crosscorr))
        dx, dy = np.where(modulus==np.max(modulus.flatten()))
        delta_x = imag1.shape[0]/2.-dx[0]
        delta_y = imag1.shape[0]/2.-dy[0]
        return (delta_x*scl, delta_y*scl)

    def _nearest_neighbours(self, refra, refdec, ra, dec):
        """
        Nearest neighbour distribution.
        """
        
        def smooth(y, box_pts):
            box = np.ones(box_pts)/box_pts
            y_smooth = np.convolve(y, box, mode='same')
            return y_smooth
            
        c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        refc = SkyCoord(ra=refra*u.degree, dec=refdec*u.degree)
        _, _, d2d, _ = refc.search_around_sky(c, 30*u.arcsec)
        hist, bins = np.histogram(d2d.arcsec, bins=np.arange(0,30,1))
        local_min = (np.diff(np.sign(np.diff(hist))) > 0).nonzero()[0] + 1
        # a = diff(sign(diff(data))).nonzero()[0] + 1 # local min+max
        # c = (diff(sign(diff(data))) < 0).nonzero()[0] + 1 # local max

        rad = bins[0]+1
        for n in local_min:
            # ~50% of the reference stars should be fine
            if hist[:n].sum()>len(refra)/2.1:
                rad = bins[n] + 1
                break

        return rad 
        
    def _platesol(self, A, B):
        """
        Least square solution of affine transformation.
        """

        # Allow arrays of different length but
        # ignore the extra points.
        N = min(len(A), len(B))

        a1 = b1 = c1 = d1 = 0.0
        a2 = b2 = 0.0
        ad = bc = ac = bd = 0.0
        for i in range(N):
            a = A[i][0]
            b = A[i][1]
            c = B[i][0]
            d = B[i][1]
            a1 += a
            b1 += b
            c1 += c
            d1 += d
            a2 += a * a
            b2 += b * b
            ad += a * d
            bc += b * c
            ac += a * c
            bd += b * d

        # Denominator.
        # It is zero iff X[i] = X[j] for every i and j in [0, n).
        # In other words, iff all the domain points are the same.
        den = N * a2 + N * b2 - a1 * a1 - b1 * b1

        if abs(den) < 1e-8:
            # The domain points are the same.
            # We assume the translation to the mean of the range
            # to be the best guess. However if N=0, assume identity.
            if N == 0:
                return [1.0, 0.0, 0.0, 0.0]
            else:
                return [1.0, 0.0, (c1 / N) - a, (d1 / N) - b]

        # Estimators
        s = (N * (ac + bd) - a1 * c1 - b1 * d1) / den
        r = (N * (ad - bc) + b1 * c1 - a1 * d1) / den
        tx = (-a1 * (ac + bd) + b1 * (ad - bc) + a2 * c1 + b2 * c1) / den
        ty = (-b1 * (ac + bd) - a1 * (ad - bc) + a2 * d1 + b2 * d1) / den

        return [s, r, tx, ty]
       
    @property
    def reference_catalogue(self):
        return self.refcat 

    @reference_catalogue.setter
    def reference_catalogue(self, val):
        self.refcat = val

    def wcsfit(self, refcat=None, niter=100, fitpv=False, maxstars=None, tol=1e-4):
        """
        Fit astrometry model.

        Parameters
        ----------

        refcat: reference astrometric catalogue

        niter: number of iterations to perform

        fitpv: if True refine the astrometric model.
        
        maxstars: maximum number of reference stars used for the match 

        tol: tolerance value to stop iteration loop
        """
        logger.info('fitting astrometry')

        # Open reference catalogue
        if refcat: self.refcat = refcat
        refcat = fits.open(self.refcat)
    
        # Read reference positions
        # TODO: We will want to do some filtering here
        refra = refcat[1].data.field('ra')
        refdec = refcat[1].data.field('dec')
            
        for i in range(1, len(self.fh)):
            # Read position of objects detected in image
            # TODO: We will want to do some filtering here

            if self.cattype=='python':
                # Experimental!! Python based source extractor
                xcat, ycat = self.find_stars(cattype='python', extension=i)
            else:
                xcat, ycat = self.find_stars(cattype='imcore', extension=i)

            # Read header
            hdr = self.fh[i].header
            w = WCS(hdr)
            # Compute x, y of reference stars
            refx, refy = w.wcs_world2pix(refra, refdec, 1)
            # Select reference stars inside the CCD
            ref_mask = (refx > 1) & (refx < hdr['NAXIS1']) & \
                       (refy > 1) & (refy < hdr['NAXIS2'])
            trefra = refra[ref_mask]
            trefdec = refdec[ref_mask]
            
            if maxstars and len(trefra)>maxstars:
                idx = random.sample(range(len(trefra)), maxstars)
                trefra = np.take(trefra, idx)
                trefdec = np.take(trefdec, idx)

            # Match to remove offset and update CRVAL
            ra, dec = w.wcs_pix2world(xcat,ycat,1)

            # Find optimum match radius
            trefx, trefy = w.wcs_world2pix(trefra, trefdec, 1)
            dx, dy = self._find_shift_crosscorr(trefx, trefy, xcat, ycat)
            xx_rad = np.sqrt(dx*dx+dy*dy)*0.101
            #nn_rad = self._nearest_neighbours(trefra, trefdec, ra, dec)
            xx_rad = max(1.0, xx_rad)
            
            idx, d2d = self._xmatch(trefra, trefdec, ra, dec)
            
            rad = min(xx_rad, np.median(d2d[d2d<xx_rad])+mad(d2d[d2d<xx_rad]))
            logger.debug('Setting match radius {}'.format(rad))
            mask = d2d<rad
            if sum(mask)<3: 
                logger.error('No matches in reference catalogue.')
                continue
            offset_ra =  np.median((ra[idx]-trefra)[mask])
            offset_dec = np.median((dec[idx]-trefdec)[mask])
            hdr['CRVAL1'] -= offset_ra
            hdr['CRVAL2'] -= offset_dec

            hdr['NUMBRMS'] = sum(mask)
            hdr['STDCRMS'] = mad(d2d[mask])
            
            logger.info('First pass offset dRA: {}, dDEC: {}, RMS: {}, NUM: {}'.format(offset_ra*3600.0, offset_dec*3600.0, mad(d2d[mask]), sum(mask)))

            # Start the iteration loop
            for iter in range(niter):
                # First compute affine transformation fitting only
                # rotation, scale (shear) and translation
                w=WCS(hdr)
                # ... ra,dec of detected objects using the current WCS
                ra, dec = w.wcs_pix2world(xcat,ycat,1)
                # ... x,y of reference objects using current WCS
                refx,refy = w.wcs_world2pix(trefra,trefdec,1)
                # ... standard coordinates of reference stars using current WCS
                refxi = hdr['CD1_1']*(refx - hdr['CRPIX1']) + \
                        hdr['CD1_2']*(refy - hdr['CRPIX2'])  
                refeta = hdr['CD2_1']*(refx - hdr['CRPIX1']) + \
                         hdr['CD2_2']*(refy - hdr['CRPIX2'])

                # ... match and select only close matches
                trefx, trefy = w.wcs_world2pix(trefra, trefdec, 1)
                dx, dy = self._find_shift_crosscorr(trefx, trefy, xcat, ycat)
                xx_rad = np.sqrt(dx*dx+dy*dy)*0.101
                xx_rad = max(1.0, xx_rad)
                
                idx, d2d = self._xmatch(trefra, trefdec, ra, dec)
                rad = min(xx_rad, np.median(d2d[d2d<xx_rad])+mad(d2d[d2d<xx_rad]))
                logger.debug('Setting match radius {}'.format(rad))
                mask = d2d<rad

                # ... perform the optimization
                if not self.use_leastsq:
                    # Using 'analytic' solution to the least square problem
                    inc = np.array(zip(xcat[idx] , ycat[idx]))
                    outc = np.array(zip(-refxi, refeta))

                    trans = self._platesol(inc[mask], outc[mask])
                    r, s, tx, ty = trans
                else:
                    # Using scipy and not worrying about all those matrix manipulations
                    # Benefit is that in the future we can define here our own minimizing function
                    # So instead of sum of squares we may want to explore sum of absolute values
                    # We also set lower and upper bounds for the scale
                    scale = hdr['CD2_2']
                    rotation = hdr['CD1_2']
                    offset = 0.0
                    p0 = (scale, rotation, offset, offset)
                    p1 = optimize.least_squares(self._errfunc, p0, 
                                                bounds=((0.01/3600.0, -1, -1, -1), (0.11/3600.0, 1, 1, 1)),
                                                args=(xcat[idx][mask] - hdr['CRPIX1'], ycat[idx][mask] - hdr['CRPIX2'], refxi[mask], refeta[mask]))
                    r, s, tx, ty = p1.x
                    logger.debug('Iter {} - cost: {}, optimality: {}'.format(iter, p1.cost, p1.optimality))

                    if not p1.success:
                        logger.ERROR('Fitting unsuccessful ({})'.format(p.message))

                logger.debug('Iter {} - r: {} s: {}, tx, ty: {} , {}'.format(iter, r, s, tx, ty))
                
                # ... update header with results
                hdr['CD1_1']=-r
                hdr['CD2_2']=r
                hdr['CD1_2']=s
                hdr['CD2_1']=s
                

                # Match again to remove offset and update CRVAL
                w = WCS(hdr)
                ra, dec = w.wcs_pix2world(xcat,ycat,1)
        
                idx, d2d = self._xmatch(trefra, trefdec, ra, dec)
                rad = min(xx_rad, np.median(d2d[d2d<xx_rad])+mad(d2d[d2d<xx_rad]))
                mask = d2d<rad
                offset_ra =  np.median((ra[idx]-trefra)[mask])
                offset_dec = np.median((dec[idx]-trefdec)[mask])
                hdr['CRVAL1'] -= offset_ra
                hdr['CRVAL2'] -= offset_dec
                
                # We enter here only if we want to fit the PV coefficients
                # By default we fit a 3-degree polynomial (22 coefficients)
                # It is solved using leastsq optimization and we need a good number of
                # stars to make it work
                if fitpv:
                    # Check for optimization of PVs
                    pv = np.zeros(self.npv*2).reshape(2,self.npv)
                    for k in range(self.npv):
                        if 'PV1_{}'.format(k) in hdr: pv[0][k] = hdr['PV1_{}'.format(k)]
                        if 'PV2_{}'.format(k) in hdr: pv[1][k] = hdr['PV2_{}'.format(k)]
                    pv = pv.flatten()
                
                    # Run least squares optimization
                    # Need to look again at the docs for better configuration
                    # Set lower and upper bounds
                    bounds_lo = np.zeros(self.npv*2).reshape(2,self.npv)
                    bounds_hi = np.zeros(self.npv*2).reshape(2,self.npv)
                    bounds_lo -= 1
                    bounds_hi += 1
                    bounds_lo[0][1] = bounds_lo[1][1] = 0.8
                    bounds_hi[0][1] = bounds_hi[1][1] = 1.2
                    p1 = optimize.least_squares(self._tpv, pv, bounds = (bounds_lo.flatten(), bounds_hi.flatten()),
                                                  args=(xcat[idx][mask], ycat[idx][mask], hdr, trefra[mask], trefdec[mask]))

                    logger.debug('Iter {} - PV fit cost: {}, optimality: {}'.format(iter, p1.cost, p1.optimality))

                    # Write to header
                    if p1.success:
                        pv = p1.x.reshape(2,self.npv)
                        for k in range(self.npv):
                            hdr['PV1_{}'.format(k)] = pv[0][k]
                            hdr['PV2_{}'.format(k)] = pv[1][k]
                    else:
                        logger.error('Least squares fitting to PV model failed ({})'.format(p1.message))
                        
                    # Match again to remove offset and update CRVAL
                    w = WCS(hdr)
                    ra, dec = w.wcs_pix2world(xcat,ycat,1)
        
                    idx, d2d = self._xmatch(trefra, trefdec, ra, dec)
                    rad = min(xx_rad, np.median(d2d[d2d<xx_rad])+mad(d2d[d2d<xx_rad]))
                    mask = d2d<rad
                    #offset_ra =  np.median((ra[idx]-trefra)[mask])
                    #offset_dec = np.median((dec[idx]-trefdec)[mask])
                    #hdr['CRVAL1'] -= offset_ra
                    #hdr['CRVAL2'] -= offset_dec
                
                logger.debug('Iter {} - RMS: {} Difference: {}'.format(iter, mad(d2d[mask]), np.abs(hdr['STDCRMS'] - mad(d2d[mask]))))
                    
                current_diff = abs(hdr['STDCRMS'] - mad(d2d[mask]))
                hdr['NUMBRMS'] = sum(mask)
                hdr['STDCRMS'] = mad(d2d[mask])
                
                if current_diff<tol:
                    break
        
            logger.info('%s CCDID: %s NUMB: %s Offset: %s, RMS: %s', i, self.fh[i].header['CCDID'],sum(mask), np.median(d2d[mask]), mad(d2d[mask]))

            # This is how to compute the residual image
            # off for now...
            #w = WCS(hdr)
            #ra, dec = w.wcs_pix2world(xcat,ycat,1)
            #trefx, trefy = w.wcs_world2pix(trefra[mask], trefdec[mask], 1)
            #val = (ra[idx]-trefra)[mask]
            #xs, ys = np.meshgrid(np.arange(1, w.wcs.crpix[0]*2, 100), np.arange(1, w.wcs.crpix[1]*2, 100))
            #res = interpolate.griddata(np.array([trefx, trefy]).transpose(), val*3600.0, np.array([xs, ys]).transpose(), method='linear')
            #out = fits.PrimaryHDU(res.T)
            #out.writeto('test.fits', clobber=True)
        
            self.fh[i].header = hdr # is this needed?
            
            
    def update_model(self, config='wcsfit.conf'):
        conf = ConfigObj(os.path.join(os.path.dirname(__file__), config))

        # Write backup
        backup = config+'.'+datetime.now().strftime('%Y%m%d')
        for i in range(100):
            b = backup+'.{}'.format(i)
            if not os.access(b, os.F_OK):
               conf.write(open(b, 'w'))
               break

        fh = fits.open(self.image_name)

        for i in range(len(fh)-1):
            pvs={}
            for k in range(self.npv):
                pvs['PV1_{}'.format(k)] = fh[i+1].header['PV1_{}'.format(k)]
                pvs['PV2_{}'.format(k)] = fh[i+1].header['PV2_{}'.format(k)]
        
            conf['CCD{}'.format(fh[i+1].header['CCDID'].strip())] = pvs

        conf.write()
        
    def writeto(self, output, **kwargs):
        self.fh.writeto(output, **kwargs)
        
    def close(self):
        if self.mode == 'update':
            self.fh.flush()
            
        self.fh.close()
        
        #if self.cattype == 'imcore':
        #    self.__fill_coords()
        
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        self.close()
    
def print_summary(input):
    """
    Print a summary of results from astrometric fit.
    """
    fh = fits.open(input)
    for i in range(1, len(fh)):
        logger.info('%s[%d] CCDID: %s NUMBRMS: %d STDCRMS: %f' % (input, i, fh[i].header.get('CCDID'), fh[i].header.get('NUMBRMS', 0), fh[i].header.get('STDCRMS', -1)))
    fh.close()
    
if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Perform astrometric calibration.')
    parser.add_argument('--copy', action='store_true', help='copy input file before modification')
    parser.add_argument('--nocat', action='store_true', help='skip catalogue generation')
    parser.add_argument('--weight', default='noconf', help='weight map')
    parser.add_argument('--debug', action='store_true', help='set logger to debug')
    parser.add_argument('--mode', default='update', help='readonly|update')
    parser.add_argument('--fitpv', action='store_true', help='fit PV coefficients')
    parser.add_argument('--wcsconf', default=os.path.join(os.path.dirname(__file__), 'wcsfit.conf'), help='location of astrometric model')
    parser.add_argument('--cattype', default="imcore", help='catalogue type to use')
    parser.add_argument('image', help='MEF image')
    parser.add_argument('refcat', help='Reference catalogue')
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)

    image_name = args.image
        
    if check_wcslib() == ERROR:
        sys.exit(1)

    if check_exe_in_path() == ERROR:
        sys.exit(1)

    # Check number of extensions and merge quadrants into ccd image if necessary
    n_extensions = get_number_of_extensions(image_name)
    if n_extensions in (4,144): # Quadrants
        merge_quadrants(image_name)   
        image_name = image_name.replace('.fit', '_ccd.fit')
    elif n_extensions in (1, 36): # CCD
        pass
    else:
        logger.error('fits file does not contain correct number of extensions ({})'.format(n_extensions))
        sys.exit(1)
    
    # Prepare the weight map so that imcore can understand it
    if args.weight.find('noconf')==-1:
        weight_name = args.weight
        n_extensions = get_number_of_extensions(image_name)
        if n_extensions in (4,144): # Quadrants
            merge_quadrants(weight_name)   
            weight_name = weight_name.replace('.fit', '_ccd.fit')
        weight_to_conf(image_name, weight_name)
    else:
        weight_name = None
        
    # Extracto objects from the image
    if not args.nocat:
        imcore(image_name, confmap='auto')

    # Fit astrometry
    with CCDImage(image_name, weight_name, mode=args.mode, save_backup=args.copy) as img:
        img.fit_astrometry(refcat=args.refcat, 
                           wcsconf=args.wcsconf,
                           fitpv=args.fitpv)
        
    # Primt summary
    print_summary(image_name)
    
