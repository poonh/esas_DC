import os,sys
import shlex,subprocess
from astropy.io import fits
from sys import *
import math
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
from scipy.ndimage import gaussian_filter


def openspec(infile,inrmf,binsize):#
      specfile = fits.open(infile)[1].data
      rmffile = fits.open(inrmf)[2].data
      channel=specfile["CHANNEL"]
      counts=specfile["COUNTS"]
      elow=rmffile["E_MIN"]
      ehigh=rmffile["E_MAX"]
      emid=(elow+ehigh)/2.
      subCounts = [counts[n:n+binsize] for n in range(0, len(counts), binsize)] 
      subCounts.append(np.concatenate((subCounts[-2],subCounts[-1]))) #combine the last two sublists together to make a new one, then remove these two sublists
      subEmid = [emid[n:n+binsize] for n in range(0, len(emid), binsize)]
      mean_counts,mean_ene=[],[]
      if binsize>1:
        for i in range(0,len(subCounts)-1):   #calculate the rate of the smoothed LC in units of counts per second
           tmp_counts=np.cumsum(subCounts[i])
           mean_counts.append(tmp_counts[-1]/float(binsize))
           tmp_ene=np.cumsum(subEmid[i])
           mean_ene.append(tmp_ene[-1]/float(binsize))
      else:        
        for i in range(0,len(emid)):
           mean_ene.append(emid[i])
           mean_counts.append(counts[i])
      return mean_ene,mean_counts



class spectrum_commands():

#combine mos1,mos2 and pn, then search for the emission peak. The final products are combined counts, combined background img, combined exposure and combined rate image, namely,counts_comb.fits,back_comb.fits,expo_comb.fits,rate_comb.fits respectively
   def imgcomb(self,Countsimg,Backimg,Expo,newimgname): #input counts,bkg and exposure images as a list
       total_num=len(Countsimg)
       head=fits.open(Countsimg[0])[0].header
       imgsize=int(head['NAXIS1'])
       pixsize=float(head['CDELT2']) #pixel size in degree
       if total_num==3:
          countsimg1=fits.open(Countsimg[0])[0].data
          countsimg2=fits.open(Countsimg[1])[0].data
          countsimg3=fits.open(Countsimg[2])[0].data
          backimg1=fits.open(Backimg[0])[0].data
          backimg2=fits.open(Backimg[1])[0].data
          backimg3=fits.open(Backimg[2])[0].data
          expo1=fits.open(Expo[0])[0].data
          expo2=fits.open(Expo[1])[0].data
          expo3=fits.open(Expo[2])[0].data
          counts_tmp=np.add(countsimg1,countsimg2)
          counts_total=np.add(counts_tmp,countsimg3)
          back_tmp=np.add(backimg1,backimg2)
          back_total=np.add(back_tmp,backimg3)
          expo_tmp=np.add(expo1,expo2)
          expo_total=np.add(expo_tmp,expo3)
          tmp=np.subtract(counts_total,back_total)          
       elif total_num==2:
          countsimg1=fits.open(Countsimg[0])[0].data
          countsimg2=fits.open(Countsimg[1])[0].data
          backimg1=fits.open(Backimg[0])[0].data
          backimg2=fits.open(Backimg[1])[0].data
          expo1=fits.open(Expo[0])[0].data
          expo2=fits.open(Expo[1])[0].data
          counts_total=np.add(countsimg1,countsimg2)
          back_total=np.add(backimg1,backimg2)
          expo_total=np.add(expo1,expo2)
          tmp=np.subtract(counts_total,back_total)
       rate=np.zeros([imgsize,imgsize]);
       ii = np.where(expo_total >0)
       rate[ii]=(counts_total[ii]-back_total[ii])/expo_total[ii]
       fits.writeto("counts_%s"%newimgname,counts_total,overwrite=True)
       fits.writeto("back_%s"%newimgname,back_total,overwrite=True)
       fits.writeto("expo_%s"%newimgname,expo_total,overwrite=True)
       fits.writeto("rate_%s"%newimgname,rate,overwrite=True)
       DATAc,HEADERc = fits.getdata("counts_%s"%newimgname, header=True)   
       DATAb,HEADERb = fits.getdata("back_%s"%newimgname, header=True)  
       DATAe,HEADERe = fits.getdata("expo_%s"%newimgname, header=True) 
       DATAr,HEADERr = fits.getdata("rate_%s"%newimgname, header=True)  
       HEADERc['CDELT2'] =  (pixsize)
       HEADERc['CDELT1'] =  (-pixsize)
       HEADERb['CDELT2'] =  (pixsize)
       HEADERb['CDELT1'] =  (-pixsize)
       HEADERe['CDELT2'] =  (pixsize)
       HEADERe['CDELT1'] =  (-pixsize)
       HEADERr['CDELT2'] =  (pixsize)
       HEADERr['CDELT1'] =  (-pixsize)
       fits.writeto("counts_%s"%newimgname,DATAc,HEADERc, overwrite=True) #combined counts img
       fits.writeto("back_%s"%newimgname,DATAb,HEADERb, overwrite=True) #combined background img
       fits.writeto("expo_%s"%newimgname,DATAe,HEADERe, overwrite=True) #combined exposure img
       fits.writeto("rate_%s"%newimgname,DATAr,HEADERr, overwrite=True) #combined rate img

   def emissionpeak(self,Rate,xc=450,yc=450,searchdis=100,kern_px=20,outfile="Centre.txt"):# The search is performed within the searching radius (searchdis in pixels) within xc,yc (also in image pixels). kern_px is the smoothing kernel
       rateimg=fits.open(Rate)[0].data
       head=fits.open(Rate)[0].header
       imgsize=int(head['NAXIS1'])
       pixsize=float(head['CDELT2']) #pixel size in degree
#       smooth_img=np.zeros([imgsize,imgsize]);
       xgrid,ygrid=np.meshgrid(1.*np.arange(imgsize),1.*np.arange(imgsize)) 
       d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)
       smooth_img=scipy.ndimage.gaussian_filter(rateimg, kern_px/(2*math.sqrt(2*math.log(2))))
       sel=[0<d,d<=searchdis,rateimg!=0.]
       flags=np.all(sel,axis=0); #number of pixels in that annulus
       highest=np.max(rateimg[flags])
       highest2=np.max(rateimg)
       yc=rateimg.argmax()/int(imgsize)
       xc=rateimg.argmax()-yc*imgsize
       os.system("ecoordconv imageset=full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits x=%s y=%s coordtype=impix > tmpcentermos1.txt"%(xc,yc))
       os.system("ecoordconv imageset=full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits x=%s y=%s coordtype=impix> tmpcentermos2.txt"%(xc,yc))
       os.system("ecoordconv imageset=full_spectrum_pn/comb-obj-im-400-2300-pn.fits x=%s y=%s coordtype=impix> tmpcenterpn.txt"%(xc,yc))
       mos1detx=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $3}' tmpcentermos1.txt"))
       mos1dety=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $4}' tmpcentermos1.txt"))
       mos2detx=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $3}' tmpcentermos2.txt"))
       mos2dety=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $4}' tmpcentermos2.txt"))
       pndetx=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $3}' tmpcenterpn.txt"))
       pndety=str(subprocess.getoutput("gawk '/DETX: DETY: / {print $4}' tmpcenterpn.txt"))
       ra=str(subprocess.getoutput("gawk '/RA: DEC: / {print $3}' tmpcenterpn.txt"))
       dec=str(subprocess.getoutput("gawk '/RA: DEC: / {print $4}' tmpcenterpn.txt"))
       os.system('echo "ra,dec = %s %s" > %s'%(ra,dec,outfile))
       os.system('echo "xc,yc(pix in python) = %s %s" >> Centre.txt'%(xc,yc))
       os.system('echo "detx,dety(mos1) = %s %s" >> Centre.txt'%(mos1detx,mos1dety))
       os.system('echo "detx,dety(mos2) = %s %s" >> Centre.txt'%(mos2detx,mos2dety))
       os.system('echo "detx,dety(pn) = %s %s" >> Centre.txt'%(pndetx,pndety))
       os.system("cat Centre.txt")

      
         
   def plotradial(self,Countsimg,Backimg,Cheese,Expo,xc,yc,rmax=12.,Nbin=20,constant="deg",outfile=None):   
       color=["b","k","r","c","m","g","b","k","r","c","m","g","b","k","r","c","m","g"]
       marker=["+","o","v","D","*","+","o","v","D","*","+","o","v","D","*"]
       total_num=len(Countsimg)
       if constant=="deg" or const=="degrees" or const=="degree":
          const=2073650.  #2073650=(3600/2.5)**2 for degree**2; size of image pixel= 2.5 arcsec
       if constant=="arcmin" or const=="arcmins" or const=="am":
          const=(60/2.5)**2 
       if constant=="pix" or const=="pixels" or const=="pixel":
          const=1
       for l in range(0,total_num):
          countsimg=fits.open(Countsimg[l])[0].data 
          backimg=fits.open(Backimg[l])[0].data 
          expo=fits.open(Expo[l])[0].data 
          cheese=fits.open(Cheese[l])[0].data 
          head=fits.open(Countsimg[l])[0].header
          imgsize=int(head['NAXIS1'])
          pixsize=float(head['CDELT2']) #pixel size in degree
          pixtoam=pixsize*60.;
          grid_img = np.zeros([imgsize,imgsize]);
          grid_sigma = np.zeros([imgsize,imgsize]);
          ii = np.where(expo*cheese >0)
          grid_img[ii]=((const*(countsimg[ii]-backimg[ii]))/expo[ii])*cheese[ii]
          ii = np.where(np.logical_and(expo*cheese >0,grid_img>0))
          grid_sigma[ii]=((const*countsimg[ii]**0.5)/expo[ii])*cheese[ii]
          xgrid,ygrid=np.meshgrid(1.*np.arange(imgsize),1.*np.arange(imgsize)) 
          d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)*pixtoam 
          binsize=np.log10(rmax)/Nbin;
          r1=10**(binsize*np.linspace(0,Nbin-1,Nbin))
          r2=10**(binsize*np.linspace(1,Nbin,Nbin))
          rmid=0.5*(r1+r2)
          sb,sberr=np.zeros(Nbin),np.zeros(Nbin)
          for ibin in range(0,Nbin):
             sel=[r1[ibin]<d,d<=r2[ibin],expo*cheese!=0.]
             flags=np.all(sel,axis=0); 
             sb[ibin]=np.sum(const*(countsimg[flags]-backimg[flags])/expo[flags])/np.sum(flags)
             sberr[ibin] = np.sqrt(np.sum(const*countsimg[flags] / expo[flags] ** 2)) / np.sum(flags) #formular for error propagation 
          cm=color[l]+marker[l]
          label_tmp=Countsimg[l].split("comb-obj-im-400-2300-")[1]
          labelname=label_tmp.split(".fits")[0]
          plt.errorbar(rmid,sb,yerr=sberr,fmt=cm,label=labelname)
       plt.xlabel("radius(arcmin)",fontsize=14)
       plt.ylabel("counts/s/$%s^{2}$"%constant,fontsize=14)
       plt.yticks(size=14)
       plt.xticks(size=14)
       tickx1=[1,2,3,4,5,6,7,8,9,10]
       tickx2=['1','2','3','4','5','6','7','8','9','10']
       ticky1=[2,4,6,10,30,50,70,100]
       ticky2=['2','4','6','10','30','50','70','100']
       plt.xscale("log")
       plt.yscale("log")
       plt.text(3,50 ,'*-exp-im-400-2300.fits',fontsize=14,color="g")
#      xlab_phrase="$M_{tot}$ ($10^{14}$$M_{\odot}$)"
#      ylab_phrase="$M_{gas}$ ($10^{14}$$M_{\odot}$)"
       plt.tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
       plt.tick_params(length=5, width=1, which='major', direction='in', right=True, top=True)  
       plt.tick_params(right="on",which='minor',color='k') #"true" for "on" for python3
       plt.tick_params(right="on",which='major', color='k')
       plt.tick_params(top="on",which='minor', color='k')
       plt.tick_params(top="on",which='major', color='k')
       plt.yticks(ticky1,ticky2,size=14)
       plt.xticks(tickx1,tickx2,size=14)
       plt.legend(loc='upper right',fontsize=14, frameon=False)
       if outfile is not None:
          plt.savefig(outfile)
       plt.show()


       

####produce a rate image with unit of your choice.input 1.counts image,2.background img, 3.epoxsure,4.cheese,5.new_file_name,6,unit=deg/arcmin/pixel (default is deg).The sigma image is name sigma_new_file_name 
#   def rateimg(self,Countsimg,Bkgimg,Exposure,Exposurenovig,Cheese,newfilename): 
   def rateimg(self,Countsimg,Bkgimg,Exposure,Cheese,newfilename,unit="deg"): 
       if unit=="deg" or unit=="degree" or unit=="degrees":
          const=(60*60/2.5)**2
       elif unit=="pixel" or unit=="pix" or unit=="pixels":
          const=1.
       elif unit=="arcmin" or unit=="am" or unit=="arcmins":
          const=(60/2.5)**2
       countsimg=fits.open(Countsimg)[0].data  
       countshead=fits.open(Countsimg)[0].header
       bkgimg=fits.open(Bkgimg)[0].data    
       exposure=fits.open(Exposure)[0].data      
       cheese=fits.open(Cheese)[0].data    
       imgsize=int(countshead['NAXIS1'])
       pixsize=countshead["CDELT2"]
       grid_img = np.zeros([imgsize,imgsize]);
       grid_sigma = np.zeros([imgsize,imgsize]);
       for i in range(imgsize):
           for j in range(imgsize):
               if exposure[i][j]!=0:
                  grid_img[i][j]=((const*(countsimg[i][j]-bkgimg[i][j]))/exposure[i][j])*cheese[i][j]
                  if grid_img[i][j]>=0:
                      grid_sigma[i][j]=((const*countsimg[i][j]**0.5)/exposure[i][j])*cheese[i][j]
                  else:
                      grid_sigma[i][j]=0
           else:
                  grid_img[i][j] = 0
                  grid_sigma[i][j] = 0
       fits.writeto("%s"%newfilename,grid_img,overwrite=True)
       fits.writeto("sigma_%s"%newfilename,grid_sigma,overwrite=True)
       DATA,HEADER = fits.getdata("%s"%newfilename, header=True)   
       DATA_sig,HEADER_sig = fits.getdata("sigma_%s"%newfilename, header=True)   
       HEADER['CDELT2'] =  (pixsize)
       HEADER['CDELT1'] =  (-pixsize)
       HEADER_sig['CDELT2'] =  (pixsize)
       HEADER_sig['CDELT1'] =  (-pixsize)
       fits.writeto("%s"%newfilename,DATA,HEADER, overwrite=True)
       fits.writeto("sigma_%s"%newfilename,DATA_sig,HEADER_sig, overwrite=True)

   def rateimgpn(self,Countsimg,Ootimg,Bkgimg,Exposure,Cheese,ootfactor,newfilename,unit="deg"): 
       if unit=="pixel":
          const=1.
       elif unit=="arcmin":
          const=(60/2.5)**2  
       elif unit=="deg":
          const=(3600/2.5)**2  
       countsimg=fits.open(Countsimg)[0].data  
       countshead=fits.open(Countsimg)[0].header
       ootimg=fits.open(Ootimg)[0].data  
       bkgimg=fits.open(Bkgimg)[0].data    
       exposure=fits.open(Exposure)[0].data    
#       exposurenovig=fits.open(Exposurenovig)[0].data    
       cheese=fits.open(Cheese)[0].data    
       imgsize=int(countshead['NAXIS1'])
       pixsize=countshead["CDELT2"]
       grid_img = np.zeros([imgsize,imgsize]);
       grid_sigma = np.zeros([imgsize,imgsize]);
       for i in range(imgsize):
           for j in range(imgsize):
               if exposure[i][j]!=0:
                  grid_img[i][j]=((const*(countsimg[i][j]-bkgimg[i][j]-ootimg[i][j]*float(ootfactor)))/exposure[i][j])*cheese[i][j]#(countsimg[i][j]/exposure[i][j])*cheese[i][j]-(bkgimg[i][j]/exposurenovig[i][j])*cheese[i][j]
                  if (countsimg[i][j]-ootimg[i][j]*float(ootfactor)) >=0:
                      grid_sigma[i][j]=((const*(((countsimg[i][j]-ootimg[i][j]*float(ootfactor)))**0.5))/exposure[i][j])*cheese[i][j]
                  else:
                      grid_sigma[i][j]=0
           else:
                  grid_img[i][j] = 0
                  grid_sigma[i][j] = 0
       fits.writeto("%s"%newfilename,grid_img,overwrite=True)
       fits.writeto("sigma_%s"%newfilename,grid_sigma,overwrite=True)
       DATA,HEADER = fits.getdata("%s"%newfilename, header=True)   
       DATA_sig,HEADER_sig = fits.getdata("sigma_%s"%newfilename, header=True)   
       HEADER['CDELT2'] =  (pixsize)
       HEADER['CDELT1'] =  (-pixsize)
       HEADER_sig['CDELT2'] =  (pixsize)
       HEADER_sig['CDELT1'] =  (-pixsize)
       fits.writeto("%s"%newfilename,DATA,HEADER, overwrite=True)
       fits.writeto("sigma_%s"%newfilename,DATA_sig,HEADER_sig, overwrite=True)

   #plot the radial profile of mos2S002-back-im-sky-400-2300.fits with binsize and up to maxrad;xc and yc are the centroid;
   def bkgrate(self,xc,yc,infile,Expo,Cheese,binsize=1.,maxrad=12.,outfile=None): #input file name as a list;binsize and max radius in arcmin
       for l in range(0,len(infile)):
           color=["b","k","r","c","m","g","b","k","r","c","m","g","b","k","r","c","m","g"]
           marker=["+","o","v","D","*","+","o","v","D","*","+","o","v","D","*"]
           bkgimg=fits.open(infile[l])[0].data
           hdr = fits.open(infile[l])[0].header
           cheese=fits.open(Cheese[l])[0].data
           expo=fits.open(Expo[l])[0].data
           imgsize=hdr['NAXIS1']
           pixsize=hdr["CDELT2"] #in degree
           xgrid,ygrid=np.meshgrid(1.*np.arange(imgsize),1.*np.arange(imgsize)) 
           pixtoarcmin=pixsize*60.
           d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)*pixtoarcmin 
           r_width=binsize*pixtoarcmin
           Nbin=int(maxrad/float(binsize))
           r=np.linspace(0,maxrad-binsize,Nbin)
           rmid=[]
           averate=np.zeros(Nbin-1);
           for i in range(0,Nbin-1): 
               selpix=[r[i]<d,d<=r[i+1],expo*cheese!=0.] #selection criteria
               flags=np.all(selpix,axis=0);  #select the pixels you want
               averate[i]=np.sum(bkgimg[flags])/np.sum(flags) 
               rmid.append(0.5*(r[i]+r[i+1]))
           cm=str(color[l])+str(marker[l])
#           plt.xlim(4,11)
#           plt.ylim(0,5e-6)
           plt.plot(rmid,averate,"%s"%cm,label="%s"%infile[l]) 
           plt.tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
           plt.tick_params(length=5, width=1, which='major', direction='in', right=True, top=True)  
           plt.tick_params(right="on",which='minor',color='k') #true for on for python3
           plt.tick_params(right="on",which='major', color='k')
           plt.tick_params(top="on",which='minor', color='k')
           plt.tick_params(top="on",which='major', color='k')
           plt.xlabel("radius (arcmin)")
           plt.ylabel("counts/pixel")
           plt.legend(loc='upper right',fontsize=11, frameon=False)
       if outfile is not None:
          plt.savefig(outfile)       
       plt.show()

   def backcounts(self,infile,inrmf,Elow,Ehigh):#input mos2S002-back.pi(infile) and the rmf to check the number of counts, Elow and Ehigh are the lower and upper energy limits in keV
      specfile = fits.open(infile)[1].data
      rmffile = fits.open(inrmf)[2].data
      channel=specfile["CHANNEL"]
      if infile.find("pn")<0:
         counts=specfile["COUNTS"]
      else:
         counts=specfile["RATE"]
      elow=rmffile["E_MIN"]
      ehigh=rmffile["E_MAX"]
      total=0
      channel_low,channel_high=0,0
      for i in range(0,len(ehigh)):
          if np.abs(Elow-elow[i])<1e-5:
             channel_low=i
          if np.abs(Ehigh-ehigh[i])<1e-5:
             channel_high=i
      for i in range(0,len(counts)):
          if i>=channel_low and i<=channel_high:
               total=total+counts[i] #int(counts[i]) for xspec
      print("total counts in %s-%s: "%(Elow,Ehigh),total)

   def backcountspn(self,infile,inrmf,Elow,Ehigh):#input pnS003-back.pi(infile) and the rmf to check the number of counts, Elow and Ehigh are the lower and upper energy limits in eV
      specfile = fits.open(infile)[1].data
      rmffile = fits.open(inrmf)[2].data
      channel=specfile["CHANNEL"]
      head = fits.open(infile)[1].header
      exposure=head["EXPOSURE"]
      counts=specfile["RATE"]
      elow=rmffile["E_MIN"]
      ehigh=rmffile["E_MAX"]
      total=0
      channel_low,channel_high=0,0
      for i in range(0,len(ehigh)-1):
          if Elow<ehigh[i] and Elow>elow[i] or Elow==elow[i]:
             channel_low=i
          if (Ehigh<ehigh[i]) and (Ehigh>elow[i]) or Ehigh==ehigh[i]:
             channel_high=i+1
      for i in range(0,len(counts)):
          if i>=channel_low and i<=channel_high:
               total=total+counts[i] #int(counts[i]) for xspec
      print("total counts in %s-%s: "%(Elow,Ehigh),total*exposure)
   
   def checkbkimg(self,ccd,ccdnum):#compare the total counts of mos2S002-back-im-det-400-2300.fits with the bkg img from fwc data,input mos1S001/mos2S002/pnS003 and the ccd chip numbers as a list("mos2S002",[1,2,3,4,5,6,7]). Put mos2S002-*obj.pi,mos2S002-im*-400-2300.fits and mos2S002-back-im-det-400-2300.fits in the folder you are running this. Check the file log/mos2_back.log
      back="%s-back-im-det-400-2300.fits"%ccd
      backfile=fits.open(back)[0].data
      print("Total counts of %s"%back,": ",np.sum(backfile))
      for k in ccdnum:
           infile="%s-im%s-400-2300.fits"%(ccd,k)
           img=fits.open(infile)[0].data
           head = fits.open(infile)[0].header
           imgsize=head['NAXIS1']
           exposure=head["EXPOSURE"]
           pifile="%s-%sobj.pi"%(ccd,k)
           pihead = fits.open(pifile)[1].header
           piexposure=pihead["EXPOSURE"]
           pix,count=0,0
           for i in range(0,imgsize):
              for j in range(0,imgsize):
                 if img[i][j]!=0:
                    pix=pix+1
                    count=count+img[i][j]        
           print("Average non-zero value and total counts of %s"%infile,": ", count/float(pix),count)
           print("exposure ratio for observation/fwc of ccd%s"%k,": ", float(piexposure)/float(exposure))  

   def totalcounts(self,infile):
           img=fits.open(infile)[0].data
           head = fits.open(infile)[0].header
           imgsize=head['NAXIS1']
           pix,count=0,0
           for i in range(0,imgsize):
              for j in range(0,imgsize):
                 if img[i][j]!=0:
                    pix=pix+1
                    count=count+img[i][j]        
           print("Total counts of image:",count)

   def pltspec(self,infilelist,inrmflist,binning,outfile=None):
       color=["b","k","r","c","m","g","b","k","r","c","m","g","b","k","r","c","m","g"]
       marker=["+","o","v","D","*","+","o","v","D","*","+","o","v","D","*"]
       ticksize=14
       plt.xlim(0.3,10)
       plt.xscale("log")
       plt.yscale("log")
       for i in range(0,len(infilelist)):
          emid,counts=openspec(infilelist[i],inrmflist[i],binsize=binning)
          plt.plot(emid,counts,color[i]+marker[i],label="%s"%infilelist[i])
       plt.legend(loc='upper right',fontsize=12, frameon=False)
       plt.xlabel("Energy(keV)",fontsize=14)
       plt.ylabel("photon counts",fontsize=14)
       plt.tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
       plt.tick_params(length=5, width=1, which='major', direction='in', right=True, top=True)  
       plt.tick_params(right="on",which='minor',color='k') #"true" for "on" for python3
       plt.tick_params(right="on",which='major', color='k')
       plt.tick_params(top="on",which='minor', color='k')
       plt.tick_params(top="on",which='major', color='k')
       tickx1=[0.3,0.5,1,5,10]
       tickx2=['0.3','0.5','1','5','10']
#       ticky1=[0.1,1,10]
#       ticky2=['0.1','1','10']
       plt.yticks(size=ticksize)
       plt.xticks(size=ticksize)
#       plt.yticks(ticky1,ticky2,size=ticksize)
       plt.xticks(tickx1,tickx2,size=ticksize)
       if outfile is not None:
          plt.savefig(outfile)
       plt.show()

   def checkfullspec(self,spec,back): #check the count rates by inputting mos2S002-obj.pi and mos2S002-back.pi
       specfile = fits.open(spec)[1].data
       backfile = fits.open(back)[1].data
       channel=specfile["CHANNEL"]
       counts=specfile["COUNTS"]
       countsbg=backfile["COUNTS"]
       head = fits.open(spec)[1].header
       expo=float(head['EXPOSURE'])
       countRate=(np.cumsum(counts)[-1]-np.cumsum(countsbg)[-1])/float(expo)
       err=(np.sqrt(np.cumsum(counts)[-1])/float(expo))
       print("exposure:", "{:.2e}".format(expo))
       print("count rates: ", "{:.5f}".format(countRate)+" +/- "+ "{:.4f}".format(err))


   def qdpplot(self,infile,ccdnum): #plot mos2S002-aug.qdp,infile=mos2-qpb.fits.gz,ccdnum=["2","3","4","5","6","7"] - a list of ccd numbers
       fig, axarr = plt.subplots(2, 3,figsize=(16,13))
       for i in ccdnum:
           qdpfile = fits.open(infile)[int(i)-1].data
           rate=qdpfile["RATE"]
           rate_err=qdpfile["RATE_ERR"]
           hard=qdpfile["HARD"]
           hard_err=qdpfile["HARD_ERR"]
           Filter=qdpfile["FILTER"]
           expo=qdpfile["EXPO"]
           tmp_evts_counts=fits.open("temp_events_mos2ccd%s_0.3-10.fits"%i)[1].header['NAXIS2']
           expo_ccd=fits.open("mos2S002-%sobj.pi"%i)[1].header['EXPOSURE']
           tmp_evts_soft=fits.open("temp_events_mos2ccd%s_0.4-0.8.fits"%i)[1].header['NAXIS2']
           tmp_evts_hard=fits.open("temp_events_mos2ccd%s_2.5-5.0.fits"%i)[1].header['NAXIS2']
           hardness=tmp_evts_hard/float(tmp_evts_soft)
           hardness_err=hardness*((np.sqrt(tmp_evts_hard)/tmp_evts_hard)**2+(np.sqrt(tmp_evts_soft)/tmp_evts_soft)**2)
           rate_fullE=tmp_evts_counts/float(expo_ccd)
           rate_fullE_err=np.sqrt(tmp_evts_counts)/float(expo_ccd)
           totalrate,numrate=0,0       
           totalhard,numhard=0,0   
           newrate,newrateerr,newhard,newharderr=[],[],[],[]
           for j in range(0,len(rate)): 
              if i=="2":  a,b=0,0
              elif i=="3": a,b=0,1
              elif i=="4": a,b=0,2
              elif i=="6": a,b=1,0
              elif i=="7": a,b=1,1
              if math.isinf(hard[j])==False and math.isnan(hard[j])==False and  math.isinf(rate[j])==False and math.isnan(rate[j])==False and Filter[j].find("Medium")>=0 or Filter[j].find("Thin")>=0 or Filter[j].find("Thick")>=0:              
                 newrate.append(rate[j])
                 newrateerr.append(rate_err[j])
                 newhard.append(hard[j])
                 newharderr.append(hard_err[j])
           axarr[a,b].errorbar(newrate,newhard,xerr=newrateerr,yerr=newharderr,fmt="r+",label="ccd"+i)
           axarr[a,b].errorbar(rate_fullE,hardness,xerr=rate_fullE_err,yerr=hardness_err,fmt="g+",label="corner data")
           axarr[a,b].set_xlim(0,0.1)
           axarr[a,b].set_ylim(0,10)
           axarr[a,b].set_xlabel("rate[0.3-10 keV]",fontsize=14)
           axarr[a,b].set_ylabel("hardness[2.5-5]/[0.4-0.8] keV",fontsize=14)
           axarr[a,b].legend(loc='upper left',fontsize=14, frameon=False)
       plt.savefig("spec_mos2_qdp.png")
       plt.show()
