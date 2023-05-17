import os
from os import *
from sys import *
#import shlex,subprocess
import numpy as np
from filter_sp_commands import filter_sp_commands


#========================================================
#write the results of mos-filter to InputFiles/goodccdlist_mos1/mos2.txt, which will be used for mos-spectra.
#a.write_ccd_anol("mosproblemccd.txt") #mosproblemccd.txt is the command.csh file after you run mos-filter
#========================================================
if len(argv)!=2:
    print("Run like this:")
    print("python3 src/filter_sp_script.py mos1S001")
    print("or")
    print("python3 src/filter_sp_script.py mos2S002")
    exit()

CCD=argv[1]

    
a=filter_sp_commands()
elow=2.5 #in keV
ehigh=12.0

#========================================================
#filter the unfiltered event list mos1S001-ori.fits for soft flare protons by making a gaussian fit to the photon counts of the FOV lightcurves
inputfile="%s-ori.fits"%CCD   #product of emchain/epchain, which is an unfiltered event list
LC=a.Lightcurves("FOV",inputfile,elow,ehigh)   #produce the FOV LCs from the inputfile in the energy range:[elow,ehigh]. Enter either "FOV" or "corner"
LCcorn=a.Lightcurves("corner",inputfile,elow,ehigh) #produce the corner LCs
a.make_gti(CCD,LC,LCcorn,elow,ehigh,binsize=60,sig=2,histogram_bin=0.05,fit_limit=1.4,pltshow="False")#the gti text file is created here. It also creates a plot showing the histogram with a gaussian fit called mos1S001_2.5_12.0_gti_ED.png. binsize in seconds, signifance = 2sigma. See the code for details
a.txt2fits("%s_gti_ED.txt"%CCD,"%s_gti_ED.fits"%CCD)  #convert gti text file to fits using ftcreate

#========================================================

#========================================================
#make the final product - mos1S001-clean.fits
a.make_Clean(CCD,"%s_gti_ED.fits"%CCD)#name of the final product is "%s-clean_ED.fits"%CCD, created using evselect with the gti file input as the 2nd parameter
#========================================================


#========================================================
#check for the hardness ratio of the ccd chips of mos1 and mos2
"""
elow1,ehigh1=0.5,0.8 #keV
elow2,ehigh2=2.5,5.0
a.problemccdcheck(CCD,elow1,ehigh1,elow2,ehigh2)
cornImage1=a.CornImage(CCD,elow1,ehigh1) #produce corner in the energy range [elow1,ehigh2]keV
cornImage2=a.CornImage(CCD,elow2,ehigh2)
print("The corner images %s and %s are produced. You can check for ccd anomaly using ds9."%(cornImage1,cornImage2))
"""
#========================================================

