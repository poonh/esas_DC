import os,sys
from astropy.io import fits
from sys import *
import numpy as np
import shlex,subprocess

class residual_sp_commands():
  def makeImage(self,ccd,elow,ehigh): #make a FOV image in the energy of choice
    eloweV=elow*1000
    ehigheV=ehigh*1000
    outputfile="%s-obj-im-%s-%skeV.fits"%(ccd,elow,ehigh)
    if os.path.isfile(outputfile) == True:
       os.system("rm %s"%outputfile)
    if ccd.find("mos1")>=0:       
       os.system("evselect table=%s-clean.fits:EVENTS filtertype=expression expression='(PATTERN<=12)&&(FLAG == 0) &&region(%s-bkg_region-sky.fits)&&(PI in [%s:%s]) ' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=%s squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes  >log/tmp.log"%(ccd,ccd,eloweV,ehigheV,outputfile))
    elif ccd.find("mos2")>=0:
       os.system("evselect table=%s-clean.fits:EVENTS filtertype=expression expression='(PATTERN<=12)&&(FLAG == 0) &&region(%s-bkg_region-sky.fits)&&(PI in [%s:%s]) ' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=%s squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes  >log/tmp.log"%(ccd,ccd,eloweV,ehigheV,outputfile))
    elif ccd.find("pn")>=0:
       os.system("evselect table=%s-clean.fits:EVENTS expression='(PATTERN <= 4)&&(FLAG == 0)&&region(%s-bkg_region-sky.fits)&&(PI in [%s:%s])' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=%s squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes  >log/tmp.log"%(ccd,ccd,eloweV,ehigheV,outputfile))
    return outputfile
  
  def SB(self,ccd,inner,outer,Countsimg,Exp): #calculate the surface brightness in the area of choice; the expousre map is only for area calculation; coordinate conversion is done here. Make sure "ecoordconv" works
       inner=inner*60/2.5
       outer=outer*60/2.5
       head=fits.open(Countsimg)[0].header
       imgsize=int(head['NAXIS1'])
       exposure=float(head['EXPOSURE'])
       if Countsimg.find("mos1")>=0:
          xc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(112,199,0)' imageset=%s-obj-image-sky.fits| grep IM_X | awk '{print $3}'"%(ccd)))
          yc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(112,199,0)' imageset=%s-obj-image-sky.fits| grep IM_X | awk '{print $4}'"%(ccd)))
       elif Countsimg.find("mos2")>=0:
          xc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(155.661,285.833,0)' imageset=%s-obj-image-sky.fits | grep IM_X | awk '{print $3}'"%(ccd)))
          yc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(155.661,285.833,0)' imageset=%s-obj-image-sky.fits | grep IM_X | awk '{print $4}'"%(ccd)))
       elif Countsimg.find("pn")>=0:
          xc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(285.8,719.722,0)' imageset=%s-obj-image-sky.fits | grep IM_X | awk '{print $3}'"%(ccd)))
          yc=float(subprocess.getoutput("ecoordconv srcexp='(DETX,DETY) IN circle(285.8,719.722,0)' imageset=%s-obj-image-sky.fits | grep IM_X | awk '{print $4}'"%(ccd)))
       pixsq2amsq=(2.5/60)**2
       exp=fits.open(Exp)[0].data  
       img=fits.open(Countsimg)[0].data  
       xgrid,ygrid=np.meshgrid(1.*np.arange(imgsize),1.*np.arange(imgsize)) 
       d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)
       ii = np.logical_and(np.logical_and(d >=inner,d<=outer),exp>0)
       rate=img[ii]/exposure
       return (np.sum(rate)/np.sum(ii))/pixsq2amsq

  def cornerspec(self,ccd,elow,ehigh): #create a corner spectrum in the energy of choice. The file mos1S001-corn.fits is needed
       if ccd.find("S")>=0:
          ccdshort=ccd.split("S")[0]
       elif ccd.find("U")>=0:
          ccdshort=ccd.split("U")[0]
       eloweV=elow*1000
       ehigheV=ehigh*1000
       outfile='%s-corn-%s-%skeV.pi'%(ccd,elow,ehigh)
       if os.path.isfile(outfile) == True:
          os.system("rm %s"%outfile)
       expr=""
       backscal=0
       if ccd.find("mos1")>=0 or ccd.find("mos2")>=0:
          ccd2=subprocess.getoutput("grep 'ccd2' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd2=="1":
             expr=expr+"(CCDNR==2)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-2oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          ccd3=subprocess.getoutput("grep 'ccd3' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd3=="1":
             expr=expr+"(CCDNR==3)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-3oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          ccd4=subprocess.getoutput("grep 'ccd4' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd4=="1":
             expr=expr+"(CCDNR==4)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-4oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          ccd5=subprocess.getoutput("grep 'ccd5' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd5=="1":
             expr=expr+"(CCDNR==5)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-5oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          ccd6=subprocess.getoutput("grep 'ccd6' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd6=="1":
             expr=expr+"(CCDNR==6)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-6oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          ccd7=subprocess.getoutput("grep 'ccd7' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccdshort)
          if ccd7=="1":
             expr=expr+"(CCDNR==7)" 
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-7oc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd)))
          expr= expr.replace(")(",")||(")
          os.system("evselect table=%s-corn.fits filtertype=expression expression='(%s)&&(PI in [%s:%s])' spectrumset='%s' energycolumn='PI' spectralbinsize=5  >log/tmp.log"%(ccd,expr,eloweV,ehigheV,outfile))
       elif ccd.find("pn")>=0: 
          for i in range(1,5):
             backscal=backscal+float(subprocess.getoutput("fkeyprint full_spectrum_%s/%s-%soc.pi+1 BACKSCAL | grep = | awk '{print (int($2))}'"%(ccdshort,ccd,str(i))))         
          os.system("evselect table=%s-corn.fits filtertype=expression expression='(PI in [%s:%s])' spectrumset='%s' energycolumn='PI' spectralbinsize=5  >log/tmp.log"%(ccd,eloweV,ehigheV,outfile))
       exposure=float(subprocess.getoutput("fkeyprint %s+1 EXPOSURE | grep = | awk '{print (int($2))}'"%outfile))
       area=float(backscal*((0.05/60.)**2))  #convert total solid angle from pixels**2 to arcmin**2
       counts=float(subprocess.getoutput("fstatistic %s COUNTS rows='-' | grep sum | awk '{print ($8)}'"%outfile))
       return counts/area/exposure

 
