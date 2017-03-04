"""
Detects objets in an image.
"""

import numpy as np
import scipy.ndimage as sn
from scipy.optimize import leastsq

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

# MAD
mad = lambda x: 1.4*np.median(np.absolute(x - np.median(x)))

def detect_peaks(image, detection_area = 2):
    """
    Takes an image and detects the peaks using a local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when the pixel's value is a 
    local maximum, 0 otherwise)

    Parameters
    ----------
    
    image:  2D marray
        the data of one FITS image.
    
    detection_area: integer, optional
         optional argument. The local maxima are found in asquare neighborhood 
         with size  (2*detection_area+1)*(2*detection_area+1). 
         Default value is detection_area = 2.
    
    Returns
    -------
    
    detected_peaks:  2D array
        Numpy array with same size as the input image, with 1 at the image 
        local maxima, 0 otherwise.
    """
    # define an 8-connected neighborhood --> image dimensionality is 2 + there are 8 neighbors for a single pixel
    neighborhood = np.ones((1+2*detection_area,1+2*detection_area))
    # apply the local maximum filter; all pixel of maximal value in their
    # neighborhood are set to 1
    local_max = sn.maximum_filter(image, footprint=neighborhood)==image
    # local_max is a mask that contains the peaks we are looking for, but also 
    # the background. In order to isolate the peaks we must remove the 
    # background from the mask. We create the mask of the background
    background = (image==0)
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_max, otherwise a line will appear 
    # along the background border (artifact of the local maximum filter)
    eroded_background = sn.binary_erosion(background, structure=neighborhood, border_value=1)
    # we obtain the final mask, containing only peaks, by removing the 
    # background from the local_max mask
    detected_peaks = local_max - eroded_background
    return detected_peaks

    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710     
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html

def detect_stars(image, weight, threshold, detection_area = 2, saturation = [False,0], margin = 5):
    """
    Takes an image and detects all objects above a certain threshold.
    Returns a boolean mask giving the objects positions (i.e. 1 when the pixel's 
    value is at an object's center, 0 otherwise)

    Parameters
    ----------
    image: array
        the data of one FITS image
    
    threshold: integer
        the objects are only detected if their maximum signal is above threshold.
    
    detection_area: integer, optional
        The object's local maxima are found in a square neighborhood with size 
        (2*detection_area+1)*(2*detection_area+1). Default value is 
        detection_area = 2.
    
    saturation: list, optional
        Two elements lists: if saturation[0] is True, the objects with maximum 
        signal >= saturation[1] are disregarded.
    
    margin: integer, optional
        The default value is 5. The objects having their center closer than 
        margin pixels from the image border are disregarded. Margin can also be 
        a list of 4 integers, in which case it represents the top, bottom, 
        right, left margins respectively.

    Returns
    -------
    
    Output: array
        Numpy array with same size as the input image, with 1 at the objects 
        centers, 0 otherwise.
    """
    # 
    sigmafilter = 1

    # select only good pixels
    imageT = (image>=threshold) & (weight>0)

    # smooth threshold, mask and filter initial image
    neighborhood = sn.generate_binary_structure(2,1)
    imageT = sn.binary_erosion(imageT, structure=neighborhood, border_value=1)
    imageT = sn.binary_dilation(imageT, structure=neighborhood, border_value=1)
    imageF = sn.filters.gaussian_filter(imageT*image,sigmafilter)

    # find local max
    peaks = detect_peaks(imageF, detection_area = detection_area)
    if isinstance(margin,list):
        peaks[-margin[0]:,:]=0
        peaks[:margin[1],:]=0
        peaks[:,-margin[2]:]=0
        peaks[:,:margin[3]]=0
    else:
        peaks[-margin:,:]=0
        peaks[:margin,:]=0
        peaks[:,-margin:]=0
        peaks[:,:margin]=0

    # removes too bright stars (possibly saturated)
    if saturation[0]:
        peaks = peaks * (image<saturation[1])
    
    return peaks

    
def find_stars(image, weight, excluded = [], included = [], detection_area = 2, saturation = [False,0], margin = 5, text_display = True):
    """
    ---------------------
    Purpose
    Takes an image and returns the approximate positions of the N brightest stars, N being specified by the user.
    ---------------------
    Inputs
    * image (2D Numpy array) = the data of one FITS image, as obtained for instance by reading one FITS image file with the function astro.get_imagedata().
    * nb_stars (integer) = the number of stars requested.
    * excluded (list) = optional argument. List of 4 elements lists, with the form [[x11,x12,y11,y12], ..., [xn1,xn2,yn1,yn2]]. Each element is a small part of the image that will be excluded for the star detection.
    This optional parameter allows to reject problematic areas in the image (e.g. saturated star(s), artefact(s), moving object(s)..)
    * included (list) = optional argument. List of 4 integers, with the form [x1,x2,y1,y2]. It defines a subpart of the image where all the stars will be detected.
    Both excluded and included can be used simultaneously.
    * detection_area (integer) = optional argument. The object's local maxima are found in a square neighborhood with size (2*detection_area+1)*(2*detection_area+1). Default value is detection_area = 2.
    * saturation (list) = optional argument. Two elements lists: if saturation[0] is True, the objects with maximum signal >= saturation[1] are disregarded.
    * margin (integer) = optional argument. The objects having their center closer than margin pixels from the image border are disregarded.
    * text_display (boolean) = optional argument. If True then some information about the star detection is printed in the console.
    ---------------------
    Output (2D Numpy array) = Numpy array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the approximate x and y positions for each star.
    ---------------------
    """

    #if text_display:
    #    print()
    #    print("--- detecting stars ---")
    #    print ("threshold, nb of detected stars:")

    #--- included ---
    if len(included)>0:
        x0, x1, y0, y1 = included
        im = image[x0:x1,y0:y1]
        wg = weight[x0:x1,y0:y1]
    else:
        im = image
        wg = weight
        x0, y0 = [0,0]
    
    #--- get subimage information ---
    ima = np.ma.array(im, mask=weight<0.5)
    #maxi = np.max(ima)
    median = np.median(ima)
    #ima=im
    #median = np.median(ima)

    #--- searching threshold ---
    #th = median + 0.1*maxi
    th = median + 3*mad(ima)
    peaks = detect_stars(im, wg, th, detection_area = detection_area, saturation = saturation, margin = margin)
    for area in excluded:
        x2,x3,y2,y3 = area
        peaks[max(x2-x0,0):max(x3-x0,0),max(y2-y0,0):max(y3-y0,0)] = 0
    x, y = np.where(peaks==1)
    x = x + x0
    y = y + y0
    nb = len(x)
    #if text_display:
    #    print( th, nb)

    xy = np.vstack([x,y]).transpose()
    maxi = image[list(x),list(y)]
    xy = xy[np.flipud(np.argsort(maxi)),:]

    #if text_display:
    #    print ("number of detected stars =",nb)

    return xy[:nb,:]

def fit_stars(image, stars, detection_area, maxshift = [False,5], error_messages = 1, full_info=False):
    """
    ---------------------
    Purpose
    This function calculates the accurate positions of some stars on one image, starting from their approximate positions.
    The calculation is done by fitting an elliptical gaussian PSF on each star.
    It can be controlled that the PSF centroid that is found is not too far from the initial position (parameter maxshift).
    ---------------------
    Inputs
    * image (2D Numpy array) = the data of one FITS image, as obtained for instance by reading one FITS image file with the function astro.get_imagedata().
    * stars (2D Numpy array) = array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the approximate positions of the stars.
    * detection_area (integer, must be >0) = defines the +/- size of the square window used for star fitting (the square window has a size 2*detection_area+1).
    The window size must be large enough to include the stars for the PSF fit, but small enough not to include contributions from other stars.
    * maxshift (list) = optional argument. 2 elements list with the form [boolean, integer]. The accurate positions are rejected if their distance to the initial approximate position is larger than a given limit.
        - if (maxshift[0]==True) the limit is maxshift[1]
        - if (maxshift[0]==False) the limit is detection_area
    If the shift is too large then the initial star position is returned.
    * error_messages (integer) = optional argument. Defines the warning messages which are displayed.
    If error_messages = 0 no warnings are displayed.
    If error_messages = 1 short warnings are displayed, giving the number of occurrence for each problem.
    If error_messages = 2 detailed warnings are displayed, listing all stars for which a particular problem occurred.
    * full_info (boolean) = optional argument. Depending on its value, the function will return different results (see explanations below).
    ---------------------
    Possible warnings
    - star can be too close from image border for psf fit
    - accurate position can be too much distant from initial position (= poor psf fit)
    ---------------------
    Output (2D Numpy array)
    If full_info = False (default), then the output is a Numpy array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains respectively, for each star:
    - the star centroid x and y positions (fit result)
    If full_info = True, then the output is a Numpy array with the form [[max1,floor1,height1,x1,y1,fhwmS1,fhwmL1,angle1],...,[maxn,floorn,heightn,xn,yn,fhwmSn,fhwmLn,anglen]], that contains respectively, for each star:
    - the value of the star maximum signal,
    - the sky background (fit result),
    - the PSF height (fit result),
    - the star centroid x and y positions (fit results), 
    - the smallest and largest full width half maximum (fit result) 
    - the angle, measured clockwise starting from the vertical direction, of the direction with the largest fwhm (fit result), and expressed in degrees. The direction of the smallest fwhm is obtained by adding 90 deg to this angle.
    ---------------------
    """
    #things to be controlled:
    #1- input position can be np.nan
    #2- input position can be too close to border to make psf fit
    #3- output position can be too far away from input position, with shift >= maxshift
    #4- output position can be too far away from input position, with shift >= detection_area and then the risk of being outside the image
    #output position is not controlled if (full_info==True)
    
    #image sizes
    size_x = image.shape[0] - detection_area - 1
    size_y = image.shape[1] - detection_area - 1
                
    #finding accurate star centers
    stars2 = []
    error = []

    if full_info:
        for i,xy in enumerate(stars):
            if np.isnan(xy[0]*xy[1]):
                stars2.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
            elif (xy[0]<detection_area)or(xy[0]>size_x)or(xy[1]<detection_area)or(xy[1]>size_y):
                stars2.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                error.append(i+1)
            else:
                p = fit_gauss_elliptical([xy[0]-detection_area,xy[1]-detection_area],image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1])
                max, floor, height, mean_x, mean_y, fwhm_1, fwhm_2, angle = p
                stars2.append([max, floor, height, mean_x, mean_y, fwhm_1, fwhm_2, angle])
        stars2 = np.array(stars2)

    else:
        for i,xy in enumerate(stars):
            if np.isnan(xy[0]*xy[1]):
                stars2.append([np.nan, np.nan])
            elif (xy[0]<detection_area)or(xy[0]>size_x)or(xy[1]<detection_area)or(xy[1]>size_y):
                stars2.append([np.nan, np.nan])
                error.append(i+1)
            else:
                p = fit_gauss_elliptical([xy[0]-detection_area,xy[1]-detection_area],image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1])
                max, floor, height, mean_x, mean_y, fwhm_1, fwhm_2, angle = p
                stars2.append([mean_x, mean_y]) 
        stars2 = np.array(stars2)

        if maxshift[0]:
            max_shift = min(maxshift[1],detection_area)
        else:
            max_shift = detection_area

        test = (np.abs(stars2[:,0]-stars[:,0]) < max_shift)*(np.abs(stars2[:,1]-stars[:,1]) < max_shift) + (np.isnan(stars2[:,0]))  #the last term is needed to have correct error messages
        test2 = np.vstack([test,test]).transpose()
        stars2 = stars2*test2 + stars*(~test2)

    #displaying error messages
    #if error_messages==1:
    #    if len(error)>0:
    #        print ('astro.fit_stars():',len(error),'/',len(stars),'object(s) too close from image border')
    #    if not(full_info):
    #        if np.sum(~test)>0:
    #            print ('astro.fit_stars():',np.sum(~test),'/',len(stars),'object(s) with poor psf fit')

    #elif error_messages==2:
    #    if len(error)>0:
    #        print ('astro.fit_stars():',len(error),'/',len(stars),'object(s) too close from image border:'+str(error)[1:-1])
    #    if not(full_info):
    #        if np.sum(~test)>0:
    #            print ('astro.fit_stars():',np.sum(~test),'/',len(stars),'object(s) with poor psf fit:'+str(list(np.arange(1,len(~test)+1)[~test]))[1:-1])

    return stars2

def centroid_stars(image, stars, detection_area, saturation = [False,0], maxshift = 2, error_messages = 1):
    """
    ---------------------
    Purpose
    This function calculates the accurate positions of some stars on one image, starting from their approximate positions.
    The star centroid is calculated with an approximate formula proposed by N.Kaiser (A photometric study of the supercluster MS0302 with the UH8K CCD camera - image processing and object catalogs, arXiv:astro-ph/9907229v1 16 Jul 1999).
    The centroid calculation involves finding the star local maximum. It can be checked that the star local maximum is not above saturation (parameter saturation).
    It can be also checked that the star centroid that is found is not too far from the star local maximum (parameter maxshift).
    ---------------------
    Inputs
    * image (2D Numpy array) = the data of one FITS image, as obtained for instance by reading one FITS image file with the function astro.get_imagedata().
    * stars (2D Numpy array) = array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the approximate positions of the stars.
    * detection_area (integer, must be >0) = defines the +/- size of the square window used to find the star maximum (the square window has a size 2*detection_area+1).
    The window size must be large enough to find the stars maxima, but small enough not to include contributions from other stars.
    * saturation (list) = optional argument. Two elements lists: if saturation[0] is True, a warning is displayed if the star maximum signal is >= than saturation[1].
    * maxshift (integer) = optional argument. The centroid positions are rejected if their distance to the star maximum position is larger than maxshift. Then the position of the star maximum is returned.
    * error_messages (integer) = optional argument. Defines the warning messages which are displayed.
    If error_messages = 0 no warnings are displayed.
    If error_messages = 1 short warnings are displayed, giving the number of occurrence for each problem.
    If error_messages = 2 detailed warnings are displayed, listing all stars for which a particular problem occurred.
    ---------------------
    Possible warnings
    - star can be too close from image border to find maximum or calculate centroid
    - star maximum signal can be above saturation
    - the centroid position can be too much distant from initial position (= poor centroid position)
    ---------------------
    Output (2D Numpy array) = Numpy array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the centroid x and y positions for each star.
    ---------------------
    """

    #--- other possibility ---
    #stars2 = np.array(map(lambda xy: np.array(np.where(image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1]== \
    #np.max(np.max(image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1])))).transpose()[0] \
    #+ xy - detection_area, stars))
    #-------------------------

    #things to be controlled:
    #1- input position can be too close to border to make max detection
    #2- star can be saturated
    #3- calculated maximum X and Y positions +/-1 must be within image boundaries to allow to calculate derivatives
    #4- final centroid position must be within image boundaries

    #image sizes
    size_0 = detection_area + max(1,maxshift)
    size_x = image.shape[0] - detection_area - 1 - max(1,maxshift)
    size_y = image.shape[1] - detection_area - 1 - max(1,maxshift)

    #finding centroids
    stars2 = []
    error_position = []
    error_saturate = []
    for i,xy in enumerate(stars):
        if (xy[0]<size_0)or(xy[0]>size_x)or(xy[1]<size_0)or(xy[1]>size_y):
            stars2.append([np.nan, np.nan])
            error_position.append(i+1)
        else:
            subimage = image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1]
            x,y = divmod(subimage.argmax(),subimage.shape[1])
            if saturation[0]:
                if subimage[x,y]>=saturation[1]:
                    stars2.append([np.nan, np.nan])
                    error_saturate.append(i+1)
                else:
                    stars2.append([x + xy[0] - detection_area, y + xy[1] - detection_area])
            else:
                stars2.append([x + xy[0] - detection_area, y + xy[1] - detection_area])

    stars2 = np.array(stars2)
    ind = np.array(np.nan_to_num(stars2),dtype=np.int)

    maxi = image[ind[:,0],ind[:,1]]
    Fx = (image[ind[:,0]+1,ind[:,1]] - image[ind[:,0]-1,ind[:,1]])/2.
    Fy = (image[ind[:,0],ind[:,1]+1] - image[ind[:,0],ind[:,1]-1])/2.
    Fxx = image[ind[:,0]+1,ind[:,1]] - 2.*maxi + image[ind[:,0]-1,ind[:,1]] + 1e-20
    Fyy = image[ind[:,0],ind[:,1]+1] - 2.*maxi + image[ind[:,0],ind[:,1]-1] + 1e-20
    
    #checks for shift smaller than maxshift on both directions
    Sx = Fx/Fxx
    Sy = Fy/Fyy
    test_ok = (np.abs(Sx)<maxshift)*(np.abs(Sy)<maxshift) + (np.isnan(stars2[:,0]))     #the second term is needed to have correct error messages

    #stars3 is automatically set to np.nan where stars2 is equal to np.nan
    #so no need for this equation: stars3 = stars3 + 0.*stars2
    stars3 = np.vstack([stars2[:,0]-Sx*test_ok,stars2[:,1]-Sy*test_ok]).transpose()
    
    #warnings
    #if error_messages == 1:
    #    if len(error_position)>0:
    #        print ('astro.centroid_stars():',len(error_position),'/',len(stars),'object(s) too close from image border')
    #    if len(error_saturate)>0:
    #        print ('astro.centroid_stars():',len(error_saturate),'/',len(stars),'object(s) above saturation')
    #    if np.sum(~test_ok)>0:
    #        print ('astro.centroid_stars():',np.sum(~test_ok),'/',len(stars),'object(s) with poor centroid calculation')

    #elif error_messages == 2:
    #    if len(error_position)>0:
    #        print ('astro.centroid_stars():',len(error_position),'/',len(stars),'object(s) too close from image border:'+str(error_position)[1:-1])
    #    if len(error_saturate)>0:
    #        print ('astro.centroid_stars():',len(error_saturate),'/',len(stars),'object(s) above saturation:'+str(error_saturate)[1:-1])
    #    if np.sum(~test_ok)>0:
    #        print ('astro.centroid_stars():',np.sum(~test_ok),'/',len(stars),'object(s) with poor centroid calculation:'+str(list(np.arange(1,len(~test_ok)+1)[~test_ok]))[1:-1])
        
    return stars3

def max_stars(image, stars, detection_area, saturation = [False,0], error_messages = 1):
    """
    ---------------------
    Purpose
    This function calculates the positions of the star local maxima, starting from an estimate.
    It can be checked that the local maxima are not above saturation (parameter saturation).
    ---------------------
    Inputs
    * image (2D Numpy array) = the data of one FITS image, as obtained for instance by reading one FITS image file with the function astro.get_imagedata().
    * stars (2D Numpy array) = array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the approximate positions of the stars.
    * detection_area (integer, must be >0) = defines the +/- size of the square window used to find the star maximum (the square window has a size 2*detection_area+1).
    The window size must be large enough to find the stars maxima, but small enough not to include contributions from other stars.
    * saturation (list) = optional argument. Two elements lists: if saturation[0] is True, a warning is displayed if the star maximum signal is >= than saturation[1].
    * error_messages (integer) = optional argument. Defines the warning messages which are displayed.
    If error_messages = 0 no warnings are displayed.
    If error_messages = 1 short warnings are displayed, giving the number of occurrence for each problem.
    If error_messages = 2 detailed warnings are displayed, listing all stars for which a particular problem occurred.
    ---------------------
    Possible warnings
    - star can be too close from image border to find maximum
    - star maximum signal can be above saturation
    ---------------------
    Output (2D Numpy array) = Numpy array with the form [[x1,y1],[x2,y2],...,[xn,yn]], that contains the local max x and y positions for each star.
    ---------------------
    """

    #things to be controlled:
    #1- input position can be too close to border to make max detection
    #2- star can be saturated

    #image sizes
    size_x = image.shape[0] - detection_area - 1
    size_y = image.shape[1] - detection_area - 1

    #finding local max
    stars2 = []
    error_position = []
    error_saturate = []
    for i,xy in enumerate(stars):
        if (xy[0]<detection_area)or(xy[0]>size_x)or(xy[1]<detection_area)or(xy[1]>size_y):
            stars2.append([np.nan, np.nan])
            error_position.append(i+1)
        else:
            subimage = image[xy[0]-detection_area:xy[0]+detection_area+1, xy[1]-detection_area:xy[1]+detection_area+1]
            x,y = divmod(subimage.argmax(),subimage.shape[1])
            if saturation[0]:
                if subimage[x,y]>=saturation[1]:
                    stars2.append([np.nan, np.nan])
                    error_saturate.append(i+1)
                else:
                    stars2.append([x + xy[0] - detection_area, y + xy[1] - detection_area])
            else:
                stars2.append([x + xy[0] - detection_area, y + xy[1] - detection_area])

    stars2 = np.array(stars2)
            
    #warnings
    #if error_messages == 1:
    #    if len(error_position)>0:
    #        print ('astro.max_stars():',len(error_position),'/',len(stars),'object(s) too close from image border')
    #    if len(error_saturate)>0:
    #        print( 'astro.max_stars():',len(error_saturate),'/',len(stars),'object(s) above saturation')
    #elif error_messages == 2:
    #    if len(error_position)>0:
    #        print ('astro.max_stars():',len(error_position),'/',len(stars),'object(s) too close from image border:'+str(error_position)[1:-1])
    #    if len(error_saturate)>0:
    #        print ('astro.max_stars():',len(error_saturate),'/',len(stars),'object(s) above saturation:'+str(error_saturate)[1:-1])

    return stars2

def fit_gauss_circular(xy, data):
    """
    ---------------------
    Purpose
    Fitting a star with a 2D circular gaussian PSF.
    ---------------------
    Inputs
    * xy (list) = list with the form [x,y] where x and y are the integer positions in the complete image of the first pixel (the one with x=0 and y=0) of the small subimage that is used for fitting.
    * data (2D Numpy array) = small subimage, obtained from the full FITS image by slicing. It must contain a single object : the star to be fitted, placed approximately at the center.
    ---------------------
    Output (list) = list with 6 elements, in the form [maxi, floor, height, mean_x, mean_y, fwhm]. The list elements are respectively:
    - maxi is the value of the star maximum signal,
    - floor is the level of the sky background (fit result),
    - height is the PSF amplitude (fit result),
    - mean_x and mean_y are the star centroid x and y positions, on the full image (fit results), 
    - fwhm is the gaussian PSF full width half maximum (fit result) in pixels
    ---------------------
    """
    
    #find starting values
    maxi = data.max()
    floor = np.ma.median(data.flatten())
    height = maxi - floor
    if height==0.0:             #if star is saturated it could be that median value is 32767 or 65535 --> height=0
        floor = np.mean(data.flatten())
        height = maxi - floor

    mean_x = (np.shape(data)[0]-1)/2
    mean_y = (np.shape(data)[1]-1)/2

    fwhm = np.sqrt(np.sum((data>floor+height/2.).flatten()))
    
    #---------------------------------------------------------------------------------
    sig = fwhm / (2.*np.sqrt(2.*np.log(2.)))
    width = 0.5/np.square(sig)
    
    p0 = floor, height, mean_x, mean_y, width

    #---------------------------------------------------------------------------------
    #fitting gaussian
    def gauss(floor, height, mean_x, mean_y, width):        
        return lambda x,y: floor + height*np.exp(-np.abs(width)*((x-mean_x)**2+(y-mean_y)**2))

    def err(p,data):
        return np.ravel(gauss(*p)(*np.indices(data.shape))-data)
    
    p = leastsq(err, p0, args=(data), maxfev=1000)
    p = p[0]
    
    #---------------------------------------------------------------------------------
    #formatting results
    floor = p[0]
    height = p[1]
    mean_x = p[2] + xy[0]
    mean_y = p[3] + xy[1]

    sig = np.sqrt(0.5/np.abs(p[4]))
    fwhm = sig * (2.*np.sqrt(2.*np.log(2.)))    
    
    output = [maxi, floor, height, mean_x, mean_y, fwhm]
    return output

def fit_gauss_elliptical(xy, data):
    """
    ---------------------
    Purpose
    Fitting a star with a 2D elliptical gaussian PSF.
    ---------------------
    Inputs
    * xy (list) = list with the form [x,y] where x and y are the integer positions in the complete image of the first pixel (the one with x=0 and y=0) of the small subimage that is used for fitting.
    * data (2D Numpy array) = small subimage, obtained from the full FITS image by slicing. It must contain a single object : the star to be fitted, placed approximately at the center.
    ---------------------
    Output (list) = list with 8 elements, in the form [maxi, floor, height, mean_x, mean_y, fwhm_small, fwhm_large, angle]. The list elements are respectively:
    - maxi is the value of the star maximum signal,
    - floor is the level of the sky background (fit result),
    - height is the PSF amplitude (fit result),
    - mean_x and mean_y are the star centroid x and y positions, on the full image (fit results), 
    - fwhm_small is the smallest full width half maximum of the elliptical gaussian PSF (fit result) in pixels
    - fwhm_large is the largest full width half maximum of the elliptical gaussian PSF (fit result) in pixels
    - angle is the angular direction of the largest fwhm, measured clockwise starting from the vertical direction (fit result) and expressed in degrees. The direction of the smallest fwhm is obtained by adding 90 deg to angle.
    ---------------------
    """

    #find starting values
    maxi = data.max()
    floor = np.ma.median(data.flatten())
    height = maxi - floor
    if height==0.0:             #if star is saturated it could be that median value is 32767 or 65535 --> height=0
        floor = np.mean(data.flatten())
        height = maxi - floor
    
    mean_x = (np.shape(data)[0]-1)/2
    mean_y = (np.shape(data)[1]-1)/2

    fwhm = np.sqrt(np.sum((data>floor+height/2.).flatten()))
    fwhm_1 = fwhm
    fwhm_2 = fwhm
    sig_1 = fwhm_1 / (2.*np.sqrt(2.*np.log(2.)))
    sig_2 = fwhm_2 / (2.*np.sqrt(2.*np.log(2.)))    

    angle = 0.

    p0 = floor, height, mean_x, mean_y, sig_1, sig_2, angle

    #---------------------------------------------------------------------------------
    #fitting gaussian
    def gauss(floor, height, mean_x, mean_y, sig_1, sig_2, angle):
    
        A = (np.cos(angle)/sig_1)**2. + (np.sin(angle)/sig_2)**2.
        B = (np.sin(angle)/sig_1)**2. + (np.cos(angle)/sig_2)**2.
        C = 2.0*np.sin(angle)*np.cos(angle)*(1./(sig_1**2.)-1./(sig_2**2.))

        #do not forget factor 0.5 in exp(-0.5*r**2./sig**2.)    
        return lambda x,y: floor + height*np.exp(-0.5*(A*((x-mean_x)**2)+B*((y-mean_y)**2)+C*(x-mean_x)*(y-mean_y)))

    def err(p,data):
        return np.ravel(gauss(*p)(*np.indices(data.shape))-data)
    
    p = leastsq(err, p0, args=(data), maxfev=1000)
    p = p[0]
    
    #---------------------------------------------------------------------------------
    #formatting results
    floor = p[0]
    height = p[1]
    mean_x = p[2] + xy[0]
    mean_y = p[3] + xy[1]
    
    #angle gives the direction of the p[4]=sig_1 axis, starting from x (vertical) axis, clockwise in direction of y (horizontal) axis
    if np.abs(p[4])>np.abs(p[5]):

        fwhm_large = np.abs(p[4]) * (2.*np.sqrt(2.*np.log(2.)))
        fwhm_small = np.abs(p[5]) * (2.*np.sqrt(2.*np.log(2.))) 
        angle = np.arctan(np.tan(p[6]))
            
    else:   #then sig_1 is the smallest : we want angle to point to sig_y, the largest
    
        fwhm_large = np.abs(p[5]) * (2.*np.sqrt(2.*np.log(2.)))
        fwhm_small = np.abs(p[4]) * (2.*np.sqrt(2.*np.log(2.))) 
        angle = np.arctan(np.tan(p[6]+np.pi/2.))
    
    output = [maxi, floor, height, mean_x, mean_y, fwhm_small, fwhm_large, angle]
    return output

def fit_moffat_circular(xy, data):
    """
    ---------------------
    Purpose
    Fitting a star with a 2D circular moffat PSF.
    ---------------------
    Inputs
    * xy (list) = list with the form [x,y] where x and y are the integer positions in the complete image of the first pixel (the one with x=0 and y=0) of the small subimage that is used for fitting.
    * data (2D Numpy array) = small subimage, obtained from the full FITS image by slicing. It must contain a single object : the star to be fitted, placed approximately at the center.
    ---------------------
    Output (list) = list with 7 elements, in the form [maxi, floor, height, mean_x, mean_y, fwhm, beta]. The list elements are respectively:
    - maxi is the value of the star maximum signal,
    - floor is the level of the sky background (fit result),
    - height is the PSF amplitude (fit result),
    - mean_x and mean_y are the star centroid x and y positions, on the full image (fit results), 
    - fwhm is the gaussian PSF full width half maximum (fit result) in pixels
    - beta is the "beta" parameter of the moffat function
    ---------------------
    """
    
    #---------------------------------------------------------------------------------
    #find starting values
    maxi = data.max()
    floor = np.ma.median(data.flatten())
    height = maxi - floor
    if height==0.0:             #if star is saturated it could be that median value is 32767 or 65535 --> height=0
        floor = np.mean(data.flatten())
        height = maxi - floor

    mean_x = (np.shape(data)[0]-1)/2
    mean_y = (np.shape(data)[1]-1)/2

    fwhm = np.sqrt(np.sum((data>floor+height/2.).flatten()))

    beta = 4
    
    p0 = floor, height, mean_x, mean_y, fwhm, beta

    #---------------------------------------------------------------------------------
    #fitting gaussian
    def moffat(floor, height, mean_x, mean_y, fwhm, beta):
        alpha = 0.5*fwhm/np.sqrt(2.**(1./beta)-1.)  
        return lambda x,y: floor + height/((1.+(((x-mean_x)**2+(y-mean_y)**2)/alpha**2.))**beta)

    def err(p,data):
        return np.ravel(moffat(*p)(*np.indices(data.shape))-data)
    
    p = leastsq(err, p0, args=(data), maxfev=1000)
    p = p[0]
    
    #---------------------------------------------------------------------------------
    #formatting results
    floor = p[0]
    height = p[1]
    mean_x = p[2] + xy[0]
    mean_y = p[3] + xy[1]
    fwhm = np.abs(p[4])
    beta = p[5]
    
    output = [maxi, floor, height, mean_x, mean_y, fwhm, beta]
    return output

def fit_moffat_elliptical(xy, data):
    """
    ---------------------
    Purpose
    Fitting a star with a 2D elliptical moffat PSF.
    ---------------------
    Inputs
    * xy (list) = list with the form [x,y] where x and y are the integer positions in the complete image of the first pixel (the one with x=0 and y=0) of the small subimage that is used for fitting.
    * data (2D Numpy array) = small subimage, obtained from the full FITS image by slicing. It must contain a single object : the star to be fitted, placed approximately at the center.
    ---------------------
    Output (list) = list with 9 elements, in the form [maxi, floor, height, mean_x, mean_y, fwhm_small, fwhm_large, angle, beta]. The list elements are respectively:
    - maxi is the value of the star maximum signal,
    - floor is the level of the sky background (fit result),
    - height is the PSF amplitude (fit result),
    - mean_x and mean_y are the star centroid x and y positions, on the full image (fit results), 
    - fwhm_small is the smallest full width half maximum of the elliptical gaussian PSF (fit result) in pixels
    - fwhm_large is the largest full width half maximum of the elliptical gaussian PSF (fit result) in pixels
    - angle is the angular direction of the largest fwhm, measured clockwise starting from the vertical direction (fit result) and expressed in degrees. The direction of the smallest fwhm is obtained by adding 90 deg to angle.
    - beta is the "beta" parameter of the moffat function   
    ---------------------
    """
    
    #---------------------------------------------------------------------------------
    #find starting values
    maxi = data.max()
    floor = np.ma.median(data.flatten())
    height = maxi - floor
    if height==0.0:             #if star is saturated it could be that median value is 32767 or 65535 --> height=0
        floor = np.mean(data.flatten())
        height = maxi - floor

    mean_x = (np.shape(data)[0]-1)/2
    mean_y = (np.shape(data)[1]-1)/2

    fwhm = np.sqrt(np.sum((data>floor+height/2.).flatten()))
    fwhm_1 = fwhm
    fwhm_2 = fwhm

    angle = 0.
    beta = 4
    
    p0 = floor, height, mean_x, mean_y, fwhm_1, fwhm_2, angle, beta

    #---------------------------------------------------------------------------------
    #fitting gaussian
    def moffat(floor, height, mean_x, mean_y, fwhm_1, fwhm_2, angle, beta):
        
        alpha_1 = 0.5*fwhm_1/np.sqrt(2.**(1./beta)-1.)
        alpha_2 = 0.5*fwhm_2/np.sqrt(2.**(1./beta)-1.)
    
        A = (np.cos(angle)/alpha_1)**2. + (np.sin(angle)/alpha_2)**2.
        B = (np.sin(angle)/alpha_1)**2. + (np.cos(angle)/alpha_2)**2.
        C = 2.0*np.sin(angle)*np.cos(angle)*(1./alpha_1**2. - 1./alpha_2**2.)
        
        return lambda x,y: floor + height/((1.+ A*((x-mean_x)**2) + B*((y-mean_y)**2) + C*(x-mean_x)*(y-mean_y))**beta)

    def err(p,data):
        return np.ravel(moffat(*p)(*np.indices(data.shape))-data)
    
    p = leastsq(err, p0, args=(data), maxfev=1000)
    p = p[0]
    
    #---------------------------------------------------------------------------------
    #formatting results
    floor = p[0]
    height = p[1]
    mean_x = p[2] + xy[0]
    mean_y = p[3] + xy[1]
    beta = p[7]
    
    #angle gives the direction of the p[4]=fwhm_1 axis, starting from x (vertical) axis, clockwise in direction of y (horizontal) axis
    if np.abs(p[4])>np.abs(p[5]):

        fwhm_large = np.abs(p[4])
        fwhm_small = np.abs(p[5])
        angle = np.arctan(np.tan(p[6]))
            
    else:   #then fwhm_1 is the smallest : we want angle to point to sig_y, the largest
    
        fwhm_large = np.abs(p[5])
        fwhm_small = np.abs(p[4])
        angle = np.arctan(np.tan(p[6]+np.pi/2.))

    output = [maxi, floor, height, mean_x, mean_y, fwhm_small, fwhm_large, angle, beta]
    return output

def extract(image, weight=None):
    if not weight:
        weight = image*0 + 1
        
    stars1 = find_stars(image, weight)
    stars2 = fit_stars(image, stars1, 3, full_info=True)
    
    t = Table()
    t['x'] = stars2[:, 4]+1
    t['y'] = stars2[:, 3]+1
    t['max'] = stars2[:, 0]
    t['floor'] = stars2[:, 1]
    t['fwhm_x'] = stars2[:, 5]
    t['fwhm_y'] = stars2[:, 6]
    t['angle'] = stars2[:, 7]
    
    return t
    

if __name__ == '__main__':
    fh = fits.open('simone_ccd_2.fits')
    
    img = fh[1].data
    
    t = extract(img)
    
    wcs = WCS(fh[1].header)
    a, d = wcs.wcs_pix2world(t['x'], t['y'], 0)

    t['ra']=a
    t['dec']=d

    t.write('test_cat.fits', overwrite=True)
