import os,sys
from astropy.io import fits
from sys import *
import math
import numpy as np
import scipy
from scipy.ndimage import gaussian_filter

def emissionpeak(fitsimg,kern_px=20):# The search is performed within the searching radius (searchdis in pixels) within xc,yc (also in image pixels). kern_px is the smoothing kernel
       img=fits.open(fitsimg)[0].data
       head=fits.open(fitsimg)[0].header
       imgsize=int(head['NAXIS1'])
       smooth_img=scipy.ndimage.gaussian_filter(img, kern_px/(2*math.sqrt(2*math.log(2))))
       highest=np.max(img)
       yc=int(img.argmax()/int(imgsize))
       xc=img.argmax()-yc*900
       return xc,yc


def EWC(fitsimg): #emission-weighted centre, equation from Mohr 1994(Cosmological Constraints from Cluster X{ray Morphologies)
    img=fits.open(fitsimg)[0].data
    head=fits.open(fitsimg)[0].header
    imgsize=int(head['NAXIS1'])
    N=np.sum(img)
    xc,yc=0,0
    for i in range(0,imgsize):
        for j in range(0,imgsize):
           xc=xc+j*img[i][j]
           yc=yc+i*img[i][j]
    return xc/N,yc/N
   


infile=["beta_cluster.fits","linear_cluster.fits"]

for fitsfile in infile:
    cx,cy=emissionpeak(fitsfile)
    cx2,cy2=EWC(fitsfile)
    print("emission peak for "+fitsfile+": ",cx,cy)
    print("emission-weighted centre for "+fitsfile+": ",cx2,cy2)

