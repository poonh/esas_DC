#!/usr/bin/env python
import os
from os import *
from sys import *
import shlex,subprocess
import numpy as np
from astropy.io import fits
from spectrum_commands import spectrum_commands
from spectrum_commands import *
#write in notes: back obj have the same area and exposure


mos1=subprocess.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=subprocess.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=subprocess.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")

ccd="mos2S002"
a=spectrum_commands()

xc,yc=420,430  #cluster emission centre from emssionpeak() in src/spectrum_commands.csh


#produce a rate image with unit of your choice.input 1.counts image,2.background img, 3.epoxsure,4.cheese,5.new_file_name,6,unit=deg/arcmin/pixel (default is deg).The sigma image is name sigma_new_file_name. The newly produced files are found in the folder you are working on.
##################################################################
"""
a.rateimg("full_spectrum_mos2/mos2S002-obj-im-400-2300.fits","full_spectrum_mos2/mos2S002-back-im-sky-400-2300.fits","full_spectrum_mos2/mos2S002-exp-im-400-2300.fits","full_spectrum_mos2/mos2S002-cheese.fits","rate_mos2.fits")
a.rateimg("full_spectrum_mos1/mos1S001-obj-im-400-2300.fits","full_spectrum_mos1/mos1S001-back-im-sky-400-2300.fits","full_spectrum_mos1/mos1S001-exp-im-400-2300.fits","full_spectrum_mos1/mos1S001-cheese.fits","rate_mos1.fits")
a.rateimgpn("full_spectrum_pn/pnS003-obj-im-400-2300.fits","full_spectrum_pn/pnS003-obj-im-400-2300-oot.fits","full_spectrum_pn/pnS003-back-im-sky-400-2300.fits","full_spectrum_pn/pnS003-exp-im-400-2300.fits","full_spectrum_pn/pnS003-cheese.fits",0.023,"rate_pn.fits") #0.023 is the oot scale factor. When you run the command comb, you can see it. It is in log/comb_pn.log

a.rateimg("full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits","full_spectrum_mos2/comb-back-im-sky-400-2300-mos2.fits","full_spectrum_mos2/comb-exp-im-400-2300-mos2.fits","full_spectrum_mos2/mos2S002-cheese.fits","comb_rate_mos2.fits")
a.rateimg("full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits","full_spectrum_mos1/comb-back-im-sky-400-2300-mos1.fits","full_spectrum_mos1/comb-exp-im-400-2300-mos1.fits","full_spectrum_mos1/mos1S001-cheese.fits","comb_rate_mos1.fits")
a.rateimg("full_spectrum_pn/comb-obj-im-400-2300-pn.fits","full_spectrum_pn/comb-back-im-sky-400-2300-pn.fits","full_spectrum_pn/comb-exp-im-400-2300-pn.fits","full_spectrum_pn/pnS003-cheese.fits","comb_rate_pn.fits") #comb-obj-im-400-2300-pn.fits has already been corrected for oot, so you don't use rateimgpn here
"""
##################################################################
"""
#combine mos1,mos2 and pn, then search for the emission peak. The final products are combined counts, combined background img, combined exposure and combined rate image, namely,counts_comb.fits,back_comb.fits,expo_comb.fits,rate_comb.fits respectively
#input counts,bkg and exposure images as a list

a.imgcomb(["full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits","full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits","full_spectrum_pn/comb-obj-im-400-2300-pn.fits"],
          ["full_spectrum_mos1/comb-back-im-sky-400-2300-mos1.fits","full_spectrum_mos2/comb-back-im-sky-400-2300-mos2.fits","full_spectrum_pn/comb-back-im-sky-400-2300-pn.fits"],
          ["full_spectrum_mos1/comb-exp-im-400-2300-mos1.fits","full_spectrum_mos2/comb-exp-im-400-2300-mos2.fits","full_spectrum_pn/comb-exp-im-400-2300-pn.fits"],"comb.fits")
a.emissionpeak("rate_comb.fits",xc=450,yc=450,searchdis=100,kern_px=20,outfile="Centre.txt") #use the rate image produced from above as the input to search for the emission peak. The image is first smoothed. The search is performed within the searching radius (searchdis in pixels) within xc,yc (also in image pixels).kern_px is the smoothing kernel.outfile is the result.
"""
###################################################################
"""
#surface brightness profile(rate=(source-bkg)/exposure). Input the following in list: source_counts_img,particle_bkg_img,cheese,exposure.xc,yc=cluster centre in img pixels. rmax=maximum radius in arcmins, Nbin=number_of_bins,constant=unit in y axis, choose one of them: "deg","arcmin","pixel"
a.plotradial(["full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits","full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits","full_spectrum_pn/comb-obj-im-400-2300-pn.fits"],
          ["full_spectrum_mos1/comb-back-im-sky-400-2300-mos1.fits","full_spectrum_mos2/comb-back-im-sky-400-2300-mos2.fits","full_spectrum_pn/comb-back-im-sky-400-2300-pn.fits"],
          ["full_spectrum_mos1/mos1S001-cheese.fits","full_spectrum_mos2/mos2S002-cheese.fits","full_spectrum_pn/pnS003-cheese.fits"],
          ["full_spectrum_mos1/mos1S001-exp-im-400-2300.fits","full_spectrum_mos2/mos2S002-exp-im-400-2300.fits","full_spectrum_pn/pnS003-exp-im-400-2300.fits"],
             xc=420,yc=430,rmax=12.,Nbin=18,constant="deg",outfile="spectrum_exp_sb.png")

a.plotradial(["full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits","full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits","full_spectrum_pn/comb-obj-im-400-2300-pn.fits"],
               ["full_spectrum_mos1/comb-back-im-sky-400-2300-mos1.fits","full_spectrum_mos2/comb-back-im-sky-400-2300-mos2.fits","full_spectrum_pn/comb-back-im-sky-400-2300-pn.fits"],
                ["full_spectrum_mos1/mos1S001-cheese.fits","full_spectrum_mos2/mos2S002-cheese.fits","full_spectrum_pn/pnS003-cheese.fits"],
                ["full_spectrum_mos1/comb-exp-im-400-2300-mos1.fits","full_spectrum_mos2/comb-exp-im-400-2300-mos2.fits","full_spectrum_pn/comb-exp-im-400-2300-pn.fits"],xc=420,yc=430,rmax=12.,Nbin=18,constant="deg",outfile="spectrum_combexp_sb.png")
"""
###################################################################

#This plots the radial profile of mos2S002-back-im-sky-400-2300.fits, input the files as a list(1.particle bkg img, 2.exposure image, 3.cheese ;default binsze=1arcmin,maxrad=12arcmin. 
#a.bkgrate(xc,yc,["full_spectrum_mos1/mos1S001-back-im-sky-400-2300.fits","full_spectrum_mos2/mos2S002-back-im-sky-400-2300.fits","full_spectrum_pn/pnS003-back-im-sky-400-2300.fits"],["full_spectrum_mos1/mos1S001-exp-im-400-2300.fits","full_spectrum_mos2/mos2S002-exp-im-400-2300.fits","full_spectrum_pn/pnS003-exp-im-400-2300.fits"],["full_spectrum_mos1/mos1S001-cheese.fits","full_spectrum_mos2/mos2S002-cheese.fits","full_spectrum_pn/pnS003-cheese.fits"])#xc,yc and file list


#################################################################
#compare the total counts of mos2S002-back-im-det-400-2300.fits with the bkg img from fwc data,
#work in the folder full_spectrum_* or put mos2S002-*obj.pi,mos2S002-im*-400-2300.fits and mos2S002-back-im-det-400-2300.fits in the folder you are running this
#check log/mos2_back.log

#a.checkbkimg("mos1S001",[1,2,4,5,7]) 
#a.checkbkimg("mos2S002",[1,2,3,4,6,7]) 
#a.checkbkimg("pnS003",[1,2,3,4])  #input either mos1 or mos2, and the ccd numbers you want in a list



#######

#This plots spectra.  
#work in the folder full_spectrum_* 
"""
ccd="mos2S002"
pitype="fc"   #fw,oc,ff,obj
infilelist=["%s-2%s.pi"%(ccd,pitype),"%s-3%s.pi"%(ccd,pitype),"%s-4%s.pi"%(ccd,pitype),"%s-6%s.pi"%(ccd,pitype),"%s-7%s.pi"%(ccd,pitype)]
inrmflist=["%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd]
a.pltspec(infilelist,inrmflist,binning=20,outfile="mos2_fc_spec.png")


ccd="pnS003"
pitype="obj-oot"
infilelist=["%s-2%s.pi"%(ccd,pitype),"%s-3%s.pi"%(ccd,pitype),"%s-4%s.pi"%(ccd,pitype),"%s-1%s.pi"%(ccd,pitype)]
inrmflist=["%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd]
a.pltspec(infilelist,inrmflist,binning=20,outfile="pn_obj_oot.png")

ccd="mos1S001"
pitype="ff"   #fw,oc,ff,obj
infilelist=["%s-1%s.pi"%(ccd,pitype),"%s-2%s.pi"%(ccd,pitype),"%s-4%s.pi"%(ccd,pitype),"%s-5%s.pi"%(ccd,pitype),"%s-7%s.pi"%(ccd,pitype)]
inrmflist=["%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd,"%s.rmf"%ccd]
a.pltspec(infilelist,inrmflist,binning=20,outfile="mos1_ff_spec.png")
"""

##########


##############compare the counts from mos2S002-back.pi in a certain energy range to that of the corresponding image;0.4=400eV (lower limit),2.3=2300eV(upper limit)
"""
a.backcounts("full_spectrum_mos1/mos1S001-back.pi","full_spectrum_mos1/mos1S001.rmf",0.4,2.3)
a.totalcounts("full_spectrum_mos1/mos1S001-back-im-det-400-2300.fits")
a.backcounts("full_spectrum_mos2/mos2S002-back.pi","full_spectrum_mos2/mos2S002.rmf",0.4,2.3)
a.totalcounts("full_spectrum_mos2/mos2S002-back-im-det-400-2300.fits")
a.backcountspn("full_spectrum_pn/pnS003-back.pi","full_spectrum_pn/pnS003.rmf",0.4,2.3)
a.totalcounts("full_spectrum_pn/pnS003-back-im-det-400-2300.fits")

"""
################


#a.checkfullspec("full_spectrum_mos2/mos2S002-obj.pi","full_spectrum_mos2/mos2S002-back.pi") #output the rate and count and compare with full-pi

#a.qdpplot("mos2-qpb.fits.gz",["2","3","4","6","7"]) #replot mos2S002-aug.qdp, input the qdp file of the calibration data and a list of the ccd chips you want, here ccd2,3,4,6,7



