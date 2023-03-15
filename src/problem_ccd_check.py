from sys import *
import numpy as np
import pyfits
import commands
import os,sys,string
from esas_commands import esas_commands

#create two filtered event list from the corner event list in the energy range [0.5-0.8]keV and [2.5-5.0]keV and compare the total count rates to check for ccd anomaly.
ccd="mos2" #change this to mos2 if you want to do mos2

elow1=0.5         #units in keV 
elow2=2.5
ehigh1=0.8
ehigh2=5.0

mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")

if ccd=="mos1":
   CCD=mos1
elif ccd=="mos2":
   CCD=mos2



a=esas_commands()

#a.make_corn_evtlist(CCD,"%s_corn.fits"%CCD)# mos1S001-corn.fits is supposed to have been created by mos-filter. If not, you can run this step.
hd1=" ccd "+"||"+" total      "+"||"+" total      "+"||"+" hardness"
hd2=" num "+"||"+" events(LE) "+"||"+" events(HE) "+"||"+" ratio"
sp="============================================="
ccdFile=open("%s_problemccd.txt"%CCD,"w")
print(CCD)
print(hd1)
print(hd2)
print(sp)

ccdFile.write(CCD+"\n")
ccdFile.write(hd1+"\n")
ccdFile.write(hd2+"\n")
ccdFile.write(sp+"\n")

for i in range(2,8):
   lowEoutput="%s_%sccd_%s_%s.fits"%(CCD,i,elow1,ehigh1)
   highEoutput="%s_%sccd_%s_%s.fits"%(CCD,i,elow2,ehigh2)
   if os.path.isfile(lowEoutput) == True:
      os.system("rm %s"%lowEoutput)
   if os.path.isfile(highEoutput) == True:
      os.system("rm %s"%highEoutput)
   expr=a.expr_ccdCheck(CCD,i,elow1,ehigh1)
   a.filteredCornEvts(CCD,"%s"%lowEoutput,expr)   #make filtered event list in the energy range[elow1,ehigh2]
   expr=a.expr_ccdCheck(CCD,i,elow2,ehigh2)
   a.filteredCornEvts(CCD,"%s"%highEoutput,expr)
   lowE_list = pyfits.open(lowEoutput)[1].data  #open the total number of events in the file
   highE_list = pyfits.open(highEoutput)[1].data
   lowEvts_num=len(lowE_list['TIME'])  #retrieve the total number of events in the file
   highEvts_num=len(highE_list['TIME'])
   if lowEvts_num!=highEvts_num: 
       line=str(i)+" "+str(lowEvts_num)+" "+str(highEvts_num)+" "+str("{:.3}".format(highEvts_num/float(lowEvts_num)))
       print(line) 
       ccdFile.write(line+"\n")
   else:
       print(i,"missing ccd")
       ccdFile.write(str(i)+" missing ccd"+"\n")

ccdFile.close()

#make corner imag in the low energy range
cornImage=a.CornImage(CCD,elow1,ehigh1)
print "The corner image %s is produced. You can check for ccd anomaly using ds9."%cornImage
