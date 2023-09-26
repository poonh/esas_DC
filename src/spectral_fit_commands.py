import os,sys
from astropy.io import fits
from sys import *
import numpy as np
import commands

class spectral_fit_commands():
#produce relevent images in [0.7-10]keV for the production of background images in the same energy range
#commands are written in full_spectrum_mos2/backgd_command_mos2.csh. It is run automatically
   def makeImage(self,CCD,elow_ori,ehigh_ori,elow_new,ehigh_new,caldb):#input elow,ehigh(original and new) in keV
       elow_ori,ehigh_ori=str(int(elow_ori*1000)),str(int(ehigh_ori*1000))
       elow_new,ehigh_new=str(int(elow_new*1000)),str(int(ehigh_new*1000))
       if CCD.find("mos1")>=0:
          ccd="mos1"
          ccdshort=CCD.split("mos")[1]
       if CCD.find("mos2")>=0:
          ccd="mos2"
          ccdshort=CCD.split("mos")[1]
       if CCD.find("pn")>=0:
          ccd="pn"
          ccdshort=CCD.split("pn")[1]
       outputname="full_spectrum_%s/backgd_command_%s.csh"%(ccd,ccd)    #the commands are written in full_spectrum_mos2/backgd_command_mos2.csh 
       if os.path.isfile(outputname) == True:
          os.system("rm %s"%outputname)
       outputfile=open(outputname,"w")        
       commandextract=file("command/FOV_imgspe_%s.csh"%ccd).readlines() #the commands are extracted from command/FOV_imgspe_mos2.csh to produce the files in the filelist
    
       if ccd=="mos1" or ccd=="mos2":
          #list of files to be created
          filelist=["%s-obj-im-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-obj-im-det-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im1-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im2-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im3-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im4-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im5-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im6-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im7-%s-%s.fits"%(CCD,elow_ori,ehigh_ori)]
       elif ccd=="pn":
          filelist=["%s-obj-im-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-obj-im-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori),"%s-obj-im-det-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-obj-im-det-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori),"%s-im1-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im1-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori),"%s-im2-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im2-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori),"%s-im3-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im3-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori),"%s-im4-%s-%s.fits"%(CCD,elow_ori,ehigh_ori),"%s-im4-%s-%s-oot.fits"%(CCD,elow_ori,ehigh_ori)]
       def action(ccd): #more to do for mos1 and mos2
           if ccd=="mos1" or ccd=="mos2":
              tmp="farith %s-im1-%s-%s.fits '%s-mask-im-%s-%s-ccd1.fits[MASK]' %s-im1-%s-%s-mask.fits MUL copyprime=yes"%(CCD,elow_new,ehigh_new,CCD,elow_ori,ehigh_ori,CCD,elow_new,ehigh_new)
              outputfile.write(tmp+"\n")
           elif ccd=="pn":
              tmp="farith %s-im1-%s-%s.fits '%s-mask-im-det-%s-%s.fits[MASK]' %s-im1-%s-%s-mask.fits MUL copyprime=yes"%(CCD,elow_new,ehigh_new,CCD,elow_ori,ehigh_ori,CCD,elow_new,ehigh_new)
              outputfile.write(tmp+"\n")
           tmp="mv %s-im1-%s-%s-mask.fits %s-im1-%s-%s.fits"%(CCD,elow_new,ehigh_new,CCD,elow_new,ehigh_new)
           outputfile.write(tmp+"\n")
 
       for name in filelist: #extract commands to create files and write the commands in full_spectrum_mos2/backgd_command_mos2.csh 
           for i in commandextract:
             if i.find(name)>=0 and i.find("#")<0:
                if i.find("evselect")>=0:
                   tmp=i.replace("%s:%s"%(elow_ori,ehigh_ori),"%s:%s"%(elow_new,ehigh_new))
                   tmp=tmp.replace("%s-%s"%(elow_ori,ehigh_ori),"%s-%s"%(elow_new,ehigh_new))
                   commandextract.remove(i)
                   newname=name.replace("%s-%s"%(elow_ori,ehigh_ori),"%s-%s"%(elow_new,ehigh_new))
                   outputfile.write("#%s\n"%newname)           
                   outputfile.write(tmp+"\n")
                   if name=="%s-im1-%s-%s.fits"%(CCD,elow_ori,ehigh_ori):
                      action(ccd)

       outputfile.write("cp ../ccf.cif . \n")
       outputfile.write("setenv SAS_CCF ../ccf.cif\n")
       if ccd=="mos1" or ccd=="mos2": #check for good ccd chips for mos1 and mos2
          ccd2=commands.getoutput("grep 'ccd2' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd)
          ccd3=commands.getoutput("grep 'ccd3' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd)
          ccd4=commands.getoutput("grep 'ccd4' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd) 
          ccd5=commands.getoutput("grep 'ccd5' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd) 
          ccd6=commands.getoutput("grep 'ccd6' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd)
          ccd7=commands.getoutput("grep 'ccd7' InputFiles/goodccdlist_%s.dat | awk '{print $2}'"%ccd)
          outputfile.write("mos_back prefix=%s caldb=%s diag=0 elow=%s ehigh=%s ccd1=1 ccd2=%s ccd3=%s ccd4=%s ccd5=%s ccd6=%s ccd7=%s >tmpback.log \n"%(ccdshort,caldb,elow_new,ehigh_new,ccd2,ccd3,ccd4,ccd5,ccd6,ccd7))

       elif ccd=="pn":
          outputfile.write("pn_back prefix=%s caldb=%s diag=0 elow=%s ehigh=%s quad1=1 quad2=1 quad3=1 quad4=1\n "%(ccdshort,caldb,elow_new,ehigh_new))
       
       #produce exposure map in the new energy range, but it is practically the same as the old one 
       outputfile.write("eexpmap attitudeset=atthk.fits eventset=%s-clean.fits:EVENTS expimageset=%s-exp-im-%s-%s.fits imageset=%s-obj-im-%s-%s.fits pimax=%s pimin=%s withdetcoords=no \n"%(CCD,CCD,elow_new,ehigh_new,CCD,elow_new,ehigh_new,ehigh_new,elow_new))
       outputfile.write("emask detmaskset=%s-mask-im-%s-%s.fits expimageset=%s-exp-im-%s-%s.fits threshold1=0.1 threshold2=0.5 \n"%(CCD,elow_new,ehigh_new,CCD,elow_new,ehigh_new))
       #produce background images in the sky coordinate
       outputfile.write("rot-im-det-sky.old prefix=%s mask=1 elow=%s ehigh=%s mode=1 \n"%(ccdshort,elow_new,ehigh_new))
       #produce the combined source image, background image and exposure map
       outputfile.write('comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=%s ehighlist=%s mask=1 prefixlist="%s" \n'%(elow_new,ehigh_new,ccdshort))
       #rename the newly produced files
       outputfile.write("mv comb-obj-im-%s-%s.fits comb-obj-im-%s-%s-%s.fits \n"%(elow_new,ehigh_new,elow_new,ehigh_new,ccd))
       outputfile.write("mv comb-back-im-sky-%s-%s.fits comb-back-im-sky-%s-%s-%s.fits \n"%(elow_new,ehigh_new,elow_new,ehigh_new,ccd))
       outputfile.write("mv comb-exp-im-%s-%s.fits comb-exp-im-%s-%s-%s.fits \n"%(elow_new,ehigh_new,elow_new,ehigh_new,ccd))
       outputfile.close()
       #run the command file,full_spectrum_mos2/backgd_command_mos2.csh,that has just been written 
       os.chdir("full_spectrum_%s"%ccd)
       os.system('tcsh -c "source %s"'%outputname.split("/")[1])
       os.chdir("..")


   def coordinateconversion(self,image,x,y):
       x=commands.getoutput("ecoordconv imageset=%s x=%s y=%s coordtype=det | grep IM_X | awk '{print $3}'"%(image,x,y))
       y=commands.getoutput("ecoordconv imageset=%s x=%s y=%s coordtype=det | grep IM_X | awk '{print $4}'"%(image,x,y))
       return float(x),float(y)
   
   #change inner here so that this can be used for CXB spectral region
   def DefineRegionCounts(self,inner,limit,width,sel,xc,yc,Img,Backimg,checkregion="yes",Region=None,ootimg=None,ootscalefactor=None):
       limit=limit#in arcsecs
       hdul = fits.open(Img)
       img=hdul[0].data;
       hdr = hdul[0].header
       pixtodeg=2.5/3600.#float(hdr['CDELT1']). If you don't bin the image, this is the value. If you read the value from the image directly, pn has a wrong value.
       pixtoam=pixtodeg*60.
       pixtoas=pixtoam*60
       hdulb = fits.open(Backimg)
       backimg=hdulb[0].data;
       if ootimg is not None:
          hduloot = fits.open(ootimg)
          pnootimg=hduloot[0].data;

       Nx=int(hdul[0].header['NAXIS1']);
       Ny=int(hdul[0].header['NAXIS2']);
       xgrid,ygrid=np.meshgrid(1.*np.arange(Nx),1.*np.arange(Ny)) 
       d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)*pixtoas
       def checkcount(in_rad,out_rad,limit,checkperform="yes"):    
           ii = np.logical_and(d >=in_rad,d<=out_rad)

           if ootimg is None:
              counts=np.sum(img[ii])-np.sum(backimg[ii])
              sn=(np.sum(img[ii])-np.sum(backimg[ii]))/np.sqrt(np.sum(img[ii]))
           else:
              counts=np.sum(img[ii])-np.sum(pnootimg[ii])*ootscalefactor-np.sum(backimg[ii])
              sn=(np.sum(img[ii])-np.sum(backimg[ii])-np.sum(pnootimg[ii])*ootscalefactor)/np.sqrt(np.sum(img[ii]))
           bkgcounts=np.sum(backimg[ii])
           if checkperform == "yes":
              if counts > sel:       
#              if sn > sel:          #uncomment this if you want to use s/n
                 return in_rad,out_rad,counts,sn,counts/bkgcounts
              else:
                 out_rad=out_rad+10
                 if out_rad <= limit: 
                    return checkcount(in_rad,out_rad,limit)
                 else: 
                    return in_rad,out_rad,counts,sn,counts/bkgcounts
           elif checkperform == "no":
              return in_rad,out_rad,counts,sn,counts/bkgcounts

       if checkregion=="yes":
#          inner=0
          outer=inner+width
          region,counts,sn,sbratio=[],[],[],[]
          region.append(inner)
          while outer< limit:
              Inner,Outer,Counts,Sn,SBratio=checkcount(inner,outer,limit)
              inner=Outer
              outer=inner+width
              region.append(Outer)
              counts.append(Counts)
              sn.append(Sn)
              sbratio.append(SBratio)
       
          if region[-1] > limit:   #if the last element in the region list is > limit, remove it 
             region.remove(region[-1])
             counts.remove(counts[-1])
             sbratio.remove(sbratio[-1])
             region[-1]=limit
             Inner,Outer,Counts,Sn,SBratio=checkcount(region[-2],region[-1],limit,checkperform="no")
             counts[-1],sn[-1],sbratio[-1]=Counts,Sn,SBratio
          elif region[-1] < limit:
             region[-1]=limit
             Inner,Outer,Counts,Sn,SBratio=checkcount(region[-2],region[-1],limit,checkperform="no")
             counts[-1],sn[-1],sbratio[-1]=Counts,Sn,SBratio

          return region,counts,sn,sbratio
        
       elif checkregion=="no":
            region=Region
            counts,sn,sbratio=[],[],[]
            for i in range(0,len(Region)-1):
                Inner,Outer,Counts,Sn,SBratio=checkcount(Region[i],Region[i+1],limit,checkperform="no")
                counts.append(Counts)
                sn.append(Sn)    
                sbratio.append(SBratio)          
            return counts,sn,sbratio                 

   def writeannulus(self,outfile,region):
       if os.path.isfile(outfile) == True:
          os.system("rm %s"%outfile)
       output=file(outfile,"w") #open
       if region.ndim==1:
           for i in region:
               output.write(str(i)+" \n")
       elif region.ndim==2:
           for i in range(0,len(region)):
               output.write(str(region[i][0])+"   "+str(region[i][1])+" \n")             
       output.close()

   def writeregion(self,ccd,xc,yc,region):
       if ccd.find("mos1")>=0:
          prefix="mos1"
       if ccd.find("mos2")>=0:
          prefix="mos2"
       if ccd.find("pn")>=0:
          prefix="pn"
       for i in range(0,len(region)-1):
           in_rad,out_rad=region[i]/0.05,region[i+1]/0.05
           outfilename="%s-%s-%s.reg"%(prefix,int(region[i]),int(region[i+1]))
           if os.path.isfile(outfilename) == True:
              os.system("rm %s"%outfilename)
           outfile=file(outfilename,"w")
           outfile.write("&&((DETX,DETY) IN circle(%s,%s,%s)&&!((DETX,DETY) IN circle(%s,%s,%s)) \n"%(str(xc),str(yc),str(out_rad),str(xc),str(yc),str(in_rad)))
           outfile.close()
           


   def howmanyCounts(self,xc,yc,Img,Backimg,inner,outer,ootimg=None,ootscalefactor=None):
       inner=inner*60
       outer=outer*60
       hdul = fits.open(Img)
       img=hdul[0].data;
       hdr = hdul[0].header
       pixtodeg=2.5/3600.#float(hdr['CDELT1']). If you don't bin the image, this is the value. If you read the value from the image directly, pn has a wrong value.
       pixtoam=pixtodeg*60.
       pixtoas=pixtoam*60
       hdulb = fits.open(Backimg)
       backimg=hdulb[0].data;
       if ootimg is not None:
          hduloot = fits.open(ootimg)
          pnootimg=hduloot[0].data;

       Nx=int(hdul[0].header['NAXIS1']);
       Ny=int(hdul[0].header['NAXIS2']);
       xgrid,ygrid=np.meshgrid(1.*np.arange(Nx),1.*np.arange(Ny)) 
       d=np.sqrt((xgrid-xc)**2+(ygrid-yc)**2)*pixtoas  
       ii = np.logical_and(d >=inner,d<=outer)

       if ootimg is None:
              counts=np.sum(img[ii])-np.sum(backimg[ii])
              sn=(np.sum(img[ii])-np.sum(backimg[ii]))/np.sqrt(np.sum(img[ii]))
       else:
              counts=np.sum(img[ii])-np.sum(pnootimg[ii])*ootscalefactor-np.sum(backimg[ii])
              sn=(np.sum(img[ii])-np.sum(backimg[ii])-np.sum(pnootimg[ii])*ootscalefactor)/np.sqrt(np.sum(img[ii]))
       sbratio=counts/np.sum(backimg[ii])
       return counts,sn,sbratio

