import pytest
import os
import sys
import astropy.io.fits as pf

from AstrometricCalibration.wcsfit import WcsFit
import glob as g
import ElementsKernel.Logging as log      # for Elements logging support

def testAstrometricCalibration():
    w = WcsFit("image.fits", "catalogue.fit")
    w.run()