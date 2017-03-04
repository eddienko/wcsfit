"""
@file: python/wcsfit/runwcsfit.py
@author: user
@date: 03/11/16
"""
#
# History
# Oct 28. 2016 (EGS): Updated to new version of wcsfit.
# Jun 17, 2016 (OH,KG): Now compatible with CCD and Quadrants handling
# Mar 18, 2016 (CG): replaced --binpath with --imcore_path and --wcsfits_path

import argparse
import ElementsKernel.Logging as log

import wcsfit

def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program
    
    @details
        See the Elements documentation for more details.
    @return
        An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('-q', action='store_true', 
                        help='process the image without merging quadrants [DELETED]')
    parser.add_argument('--async', action='store_true', 
                        help='perform astrometric calibration asynchronously [DELETED]')
    parser.add_argument('--copy', action='store_true', 
                        help='copy input file before modification')
    parser.add_argument('--nocat', action='store_true', 
                        help='skip catalogue generation')
    parser.add_argument('--refcat',
                        help='reference catalogue')
    parser.add_argument('--cattype', default='6', 
                        help='catalogue type to generate (3=basic, 6=long) [DELETED]')
    parser.add_argument('--weight', default='noconf',
                        help='weight map')
    parser.add_argument('--imcore_path', 
                        help='CASUtools imcore path')
    parser.add_argument('--wcsfit_path', 
                        help='CASUtools wcsfit path [DELETED]')
    parser.add_argument('--wcsconf', default=os.path.join(os.path.dirname(__file__), 'wcsfit.conf'),
                        help='location of astrometric model')
    parser.add_argument('--fitpv', action='store_true', help='fit PV coefficients')
    parser.add_argument('--mode', default='update', help='readonly|update')
    parser.add_argument('--image', 
                        help='FITS image')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is 
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('runwcsfit')

    logger.info('#')
    logger.info('# Entering runwcsfit mainMethod()')
    logger.info('#')
    
    wcsfit.IMCORE=args.imcore_path
    
    #if args.debug:
    #    logger.setLevel(logging.DEBUG)

    image_name = args.image
       
    if wcsfit.check_wcslib() == ERROR:
        sys.exit(1)

    if wcsfit.check_exe_in_path() == ERROR:
        sys.exit(1)

    # Check number of extensions and merge quadrants into ccd image if necessary
    n_extensions = wcsfit.get_number_of_extensions(image_name)
    if n_extensions in (4,144): # Quadrants
        wcsfit.merge_quadrants(image_name)   
        image_name = image_name.replace('.fit', '_ccd.fit')
    elif n_extensions in (1, 36): # CCD
        pass
    else:
        logger.error('fits file does not contain correct number of extensions ({})'.format(n_extensions))
        sys.exit(1)
    
    # Deal with the weight map if needed
    if args.weight.find('noconf')==-1:
        weight_name = args.weight
        n_extensions = wcsfit.get_number_of_extensions(image_name)
        if n_extensions in (4,144): # Quadrants
            wcsfit.merge_quadrants(weight_name)   
            weight_name = weight_name.replace('.fit', '_ccd.fit')
        wcsfit.weight_to_conf(image_name, weight_name)
        
    # Fit astrometry
    with wcsfit.CCDImage(image_name, mode=args.mode, save_backup=args.copy) as img:
        img.fit_astrometry(refcat=args.refcat, makecat=args.nocat, 
                           wcsconf=args.wcsconf,
                           fitpv=args.fitpv)
        
    # Primt summary
    print_summary(image_name)

    logger.info('#')
    logger.info('# Exiting runwcsfit mainMethod()')
    logger.info('#')
