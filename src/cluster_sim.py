import os
import sys,string
import numpy as np;
from astropy.io import fits


pixsize=2.5/3600. #pix size in degree for an 900*900 XMM Newton image
Nx=900   #This is a 900*900 pixels image
Ny=900

xc=Nx/2  #x and y centre of the simulated images
yc=Ny/2

xgrid=np.zeros([Ny,Nx]);
ygrid=np.zeros([Ny,Nx]);

for ix in np.arange(Nx):
    for iy in np.arange(Ny):
        xgrid[iy,ix]=1.*ix;
        ygrid[iy,ix]=1.*iy;

d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)


def betamodel(n0,rc,beta):
    return n0*(1+d*d/rc/rc)**(0.5-3*beta);
    
def linear(b):
     return d+b
   
#beta and linear model parameters
n0=20
rc=500.0
beta=5
b=5


output_fitsname1="beta_cluster.fits"
output_fitsname2="linear_cluster.fits" 

output_img1=betamodel(n0,rc,beta)
output_img2=linear(b)

fits.writeto(output_fitsname1, output_img1,overwrite=True)
fits.writeto(output_fitsname2, output_img2,overwrite=True)

#write the header
DATA1,HEADER1 = fits.getdata(output_fitsname1, header=True)   
DATA2,HEADER2 = fits.getdata(output_fitsname2, header=True)  

HEADER1['CDELT2'] =  (pixsize)
HEADER1['CDELT1'] =  (-pixsize)
HEADER2['CDELT2'] =  (pixsize)
HEADER2['CDELT1'] =  (-pixsize)

fits.writeto(output_fitsname1,DATA1,HEADER1, overwrite=True) 
fits.writeto(output_fitsname2,DATA2,HEADER2, overwrite=True) 
   


