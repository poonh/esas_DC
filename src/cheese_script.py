from sys import *
import numpy as np
#import pyfits
import commands
import os,sys,string
from cheese_commands import cheese_commands
#produce rate image

a=cheese_commands()


#=============================================================
#write extracted parameters from emllist.fits to emllist.txt using fdump
#a.fitstotxt("emllist.fits","emllist.txt","ID_INST ID_BAND ID_CLUSTER FLUX DET_ML DIST_NN") 
#=============================================================

#=============================================================
#make rate image using farith #counts img,exposure img,outputfile,DIV/ADD/SUB/MUL
#a.makeImg("mos2S002-obj-im.fits","mos2S002-obj-imexp.fits","mos2S002-rate.fits","DIV")
#a.makeImg("mos1S001-obj-im.fits","mos1S001-obj-imexp.fits","mos1S001-rate.fits","DIV")
#a.makeImg("pnS003-obj-im.fits","pnS003-obj-imexp.fits","pnS003-rate.fits","DIV")

#=============================================================
#convert mos2S002-bkg_region-sky.fits to text format  to see the coordinates of the mask point sources. mos2_sky.reg can be loaded in ds9.
"""
a.regionextract("mos1S001-bkg_region-sky.fits","mos1_sky.reg")#infile,outfile
a.regionextract("mos1S001-bkg_region-det.fits","mos1_det.reg")
a.regionextract("mos2S002-bkg_region-sky.fits","mos2_sky.reg")
a.regionextract("mos2S002-bkg_region-det.fits","mos2_det.reg")
a.regionextract("pnS003-bkg_region-sky.fits","pn_sky.reg")
a.regionextract("pnS003-bkg_region-det.fits","pn_det.reg")
"""
#=============================================================

#=============================================================
#update the flux of the summary band of emllist.fits by taking the average so that there is no INDEF value. Make a copy first in case this step goes wrong
#os.system("cp emllist.fits emllist.fits.ori")
#a.emllistUpdate()
#run region after updating emllist.fits so that all ccds have the same masks
#=============================================================
"""
a.remakeRegion("mos2S002-clean.fits","mos2S002-bkg_region-det.fits","detxy")#infile,outfile,coor
a.remakeRegion("mos2S002-clean.fits","mos2S002-bkg_region-sky.fits","xy")
a.remakeRegion("mos1S001-clean.fits","mos1S001-bkg_region-det.fits","detxy")
a.remakeRegion("mos1S001-clean.fits","mos1S001-bkg_region-sky.fits","xy")
a.remakeRegion("pnS003-clean.fits","pnS003-bkg_region-det.fits","detxy")
a.remakeRegion("pnS003-clean.fits","pnS003-bkg_region-sky.fits","xy")
"""
#=============================================================

#output the updated fits file to a text file to check which row you need to delete
"""
a.fitstotxt("mos2S002-bkg_region-sky.fits","mos2_skynew.txt","X Y R") 
a.fitstotxt("mos2S002-bkg_region-det.fits","mos2_detnew.txt","DETX DETY R") 
a.fitstotxt("mos1S001-bkg_region-sky.fits","mos1_skynew.txt","X Y R") 
a.fitstotxt("mos1S001-bkg_region-det.fits","mos1_detnew.txt","DETX DETY R") 
a.fitstotxt("pnS003-bkg_region-sky.fits","pn_skynew.txt","X Y R") 
a.fitstotxt("pnS003-bkg_region-det.fits","pn_detnew.txt","DETX DETY R") 
"""
#=============================================================
#convert the updated fits file to ds9 readable text file to be loaded with a fits image
#a.regionextract("mos2S002-bkg_region-sky.fits","mos2_sky_v2.reg")

#=============================================================

#set the radius of the mask masking the cluster in the centre to 0. enter either mos2S002-bkg_region-sky.fits or mos2S002-bkg_region-det.fits. A new mos2S002-cheese.fits is also produced
#rowchange() uses fpartab, rowchange_py() uses python. They are the same.
"""
a.rowchange_py("mos2S002-bkg_region-sky.fits","R",0,[1,12],[1,2]) #filename,column_name,new_value,list of row_num,list of element_num
a.rowchange_py("mos1S001-bkg_region-sky.fits","R",0,[1,12],[1,2])
a.rowchange_py("pnS003-bkg_region-sky.fits","R",0,[1,12],[1,2])
"""
"""
a.rowchange("mos2S002-bkg_region-sky.fits","R",0,[1,12],[1,2])
a.rowchange("mos1S001-bkg_region-sky.fits","R",0,[1,12],[1,2])
a.rowchange("pnS003-bkg_region-sky.fits","R",0,[1,12],[1,2])
"""


#==============================================
"""
#add an extra mask by merging two files
#create an empty fits file to be merged to another file. The input file is the one you want to merge to and the empty file name.  
a.emptyFits("mos2S002-bkg_region-sky.fits","empty.fits") #create empty.fits in the same format as mos2S002-bkg_region-sky.fits
a.emptyFits("mos2S002-bkg_region-det.fits","emptydet.fits")

a.NewMask("mos2S002","empty.fits",1,17515.893,19017.939,500,500,0) #ccd, fits file to be updated,row number to be inserted in empty.fits,x_coor,y_coor,ellipse semi major and minor axis,rotation angle
a.NewMask("mos2S002","emptydet.fits",1,17515.893,19017.939,500,500,0) #The detector coordinate. Coordinate conversion is performed using ecoordconv.
a.NewMask("mos2S002","empty.fits",2,28728.593,19418.211,1006.1019,553.08489,50)
a.NewMask("mos2S002","emptydet.fits",2,28728.593,19418.211,1006.1019,553.08489,50)

a.Mergefits("empty.fits","mos2S002-bkg_region-sky.fits") #merge the two files. The file in the latter position is updated.
a.Mergefits("emptydet.fits","mos2S002-bkg_region-det.fits") 
#make a new mask image(mos2S002-cheese-new.fits) with a new region list(reglist)
a.MakeMask("mos2S002","mos2S002-cheese-new.fits","mos2S002-bkg_region-sky.fits")#make a new cheese image,mos2S002-cheese-new.fits, using an updated region file, mos2S002-bkg_region-sky.fits.

"""
