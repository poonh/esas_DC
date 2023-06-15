import os,sys
import commands
from astropy.io import fits
from sys import *
import math
import numpy as np


def extractcoor(ccd,x,y): 
        img="%s-obj-image-det.fits"%ccd
        if os.path.isfile(img) == False:
           print("%s is missing. Cannot perform coordinate conversion. %s is a product of mos-filter. If it is not present, do the following. Just copy and paste at the terminal and do this part again."%(img,img))
           print("evselect table=%s-clean.fits withimageset=yes imageset=%s-obj-image-det.fits xcolumn='DETX' ximagesize=780 ximagemax=19500 ximagemin=-19499 ycolumn='DETY' yimagesize=780 yimagemax=19500 yimagemin=-19499 imagebinning=imageSize ignorelegallimits=yes"%(ccd,ccd))
           exit()
        elif os.path.isfile(img) == True:
           os.system("ecoordconv imageset=%s-obj-image-det.fits x=%s y=%s coordtype=POS > tmp.dat" %(ccd,x,y))#The image is a product of mos-filter. If it is not present, just do the following         
           DETX=commands.getoutput("head -6 tmp.dat | tail -1 | gawk '{print $3}'")
           DETY=commands.getoutput("head -6 tmp.dat | tail -1 | gawk '{print $4}'")       
           os.system("cat tmp.dat")
           os.system("rm tmp.dat")
        return DETX,DETY

class cheese_commands():
   def regionextract(self,infile,outfile):#convert mos2S002-bkg_region-sky.fits (product of cheese) into ds9 readable region. 
      outputfile=open(outfile,"w")
      fitsfile = fits.open(infile)[1].data
      head = fits.open(infile)[1].header
      x=str(head['TTYPE2'])
      y=str(head['TTYPE3'])
      X=fitsfile[x][:,0].astype(str)
      Y=fitsfile[y][:,0].astype(str)
      R1=fitsfile["R"][:,0].astype(str)
      R2=fitsfile["R"][:,1].astype(str)
      Rotang=fitsfile["ROTANG"][:,0].astype(str)
      firstline='global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n'
      outputfile.write(firstline)
      for i in range(0,len(X)):
          outputfile.write("ellipse("+X[i]+","+Y[i]+","+R1[i]+","+R2[i]+","+Rotang[i]+")\n")
      outputfile.close();
   def fitstotxt(self,infile,outfile,parameters):#convert a fits text file into a text file using fdump 
      os.system('fdump %s+1 %s "%s" - pagewidth=256 prhead=no clobber=yes'%(infile,outfile,parameters))
   def makeImg(self,infileimg,infileexp,outfile,choice):#produce an image using farith,choice=MUL/DIV/SUB/ADD
      os.system('farith %s %s %s %s clobber=yes'%(infileimg,infileexp,outfile,choice))
   def emllistUpdate(self):  #update the flux of emllist.fits so that there is no INDEF for the summary band
      filename="emllist.fits"
      fitsfile = fits.open(filename,mode="update")
      hdu=fitsfile[1].data
      n=hdu.shape[0]
      flux=hdu["FLUX"]
      for i in xrange(0,n,4):
          allflux=[]
          ccdnum=0
          if math.isnan(flux[i])==True:
             if math.isnan(flux[i+1])==False:
                allflux.append(flux[i+1])
             if math.isnan(flux[i+2])==False:
                allflux.append(flux[i+2])
             if math.isnan(flux[i+3])==False:
                allflux.append(flux[i+3]    )       
             flux[i]=np.mean(allflux)
      fitsfile.close()
   def remakeRegion(self,infile,outfile,coor):#infile=mos2S002-clean.fits,outfile=your choice,coor=detxy/xy
          os.system("region eventset=%s operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(FLUX >= 1e-18)&&(DET_ML >= 15)&&(ID_INST == 0)&&(DIST_NN >= 0)&&(ID_BAND == 0)' bkgregionset=%s  bkgfraction=0.25 radiusstyle=contour nosrcellipse=no outunit=%s verbosity=1"%(infile,outfile,coor))
          if coor=="xy":
            ccd=infile.split("-clean")[0]
            os.system("make_mask inimage=%s-obj-im.fits inmask=%s-mask-im.fits outmask=%s-cheese.fits reglist=%s"%(ccd,ccd,ccd,outfile))
   def rowchange(self,infile,parameter,value,row,element,ccd=None):#change the value of an element in a row of a parameter in a filename
          if infile.find("sky")>=0:
             skyfile=infile
             detfile=skyfile.replace("sky","det")
          elif infile.find("det")>=0:
             detfile=infile
             skyfile=skyfile.replace("det","sky")
          for Files in (detfile,skyfile):
            for i in range(0,len(row)):
               for j in range(0,len(element)):
                  os.system("fpartab %s %s+1 %s row=%s element=%s"%(str(value),str(Files),str(parameter),str(row[i]),str(element[j])))  
          os.system("make_mask inimage=%s-obj-im.fits inmask=%s-mask-im.fits outmask=%s-cheese.fits reglist=%s"%(ccd,ccd,ccd,skyfile))        
   def rowchange_py(self,infile,parameter,value,row,element):
       ccd=infile.split("-bkg_region")[0]
       fitsFile = fits.open(infile,mode="update")
       if infile.find("sky")>=0:
          skyfile=infile
          detfile=skyfile.replace("sky","det")
       elif infile.find("det")>=0:
          detfile=infile
          skyfile=detfile.replace("det","sky")
       for Files in (skyfile,detfile):
          fitsFile = fits.open(Files,mode="update")
          for i in range(0,len(row)):
              for j in range(0,len(element)):
                  fitsFile[1].data[parameter][row[i]-1][element[j]-1]=value
          fitsFile.close()
       os.system("make_mask inimage=%s-obj-im.fits inmask=%s-mask-im.fits outmask=%s-cheese.fits reglist=%s"%(ccd,ccd,ccd,skyfile))  
   def DeleteMask(self,oldfile,newfile):
       Old=open(oldfile,"r").readlines()
       New=open(newfile,"r").readlines()
       old_x,new_x,remove_index=[],[],[]
       for i in Old:
           if i.find("ellipse")>=0:
              old_x.append(i.split(",")[1])
       for i in New:
           if i.find("ellipse")>=0:  
              tmp=i.split(",")[1]
              if tmp.find(".")<0:
                 tmp=str(str(tmp)+(".0"))
              new_x.append(tmp)    
       for i in old_x: 
           if i not in new_x:
              remove_index.append(old_x.index(i)+1) 
       return remove_index
   def emptyFits(self,infile,output):#create an empty fits file to be merged to another file. The input file is the one you want to merge to.Set everything to 0 in output.
       os.system("cp %s %s"%(infile,output))
       fitsfile = fits.open("%s"%output,mode="update")
       hdu=fitsfile[1].data
       n=hdu.shape[0]
       x=fitsfile[1].columns[1].name
       y=fitsfile[1].columns[2].name
       R=hdu["R"]
       X=hdu[x]
       Y=hdu[y]
       Rotang=hdu["ROTANG"]
       for i in xrange(0,n):
           R[i][0]=0
           R[i][1]=0
           X[i][0]=0
           Y[i][0]=0
           Rotang[i][0]=0
       fitsfile.close()
   def NewMask(self,ccd,infile,row,x,y,r1,r2,rotang): #input the values of the new mask in empty.fits
       fitsfile = fits.open("%s"%infile)
       hdu=fitsfile[1].data
       xcoor=fitsfile[1].columns[1].name
       if xcoor=="X":
          os.system("fpartab %s %s+1 X row=%s element=1"%(str(x),infile,row))
          os.system("fpartab %s %s+1 Y row=%s element=1"%(str(y),infile,row))
          os.system("fpartab %s %s+1 R row=%s element=1"%(str(r1),infile,row))
          os.system("fpartab %s %s+1 R row=%s element=2"%(str(r2),infile,row))
          os.system("fpartab %s %s+1 ROTANG row=%s element=1"%(str(rotang),infile,row))
       elif xcoor=="DETX":
          detx,dety=extractcoor(ccd,x,y)
          os.system("fpartab %s %s+1 DETX row=%s element=1"%(detx,infile,row))
          os.system("fpartab %s %s+1 DETY row=%s element=1"%(dety,infile,row))
          os.system("fpartab %s %s+1 R row=%s element=1"%(r1,infile,row))
          os.system("fpartab %s %s+1 R row=%s element=2"%(r2,infile,row))
          os.system("fpartab %s %s+1 ROTANG row=%s element=1"%(rotang,infile,row))
       else:
          print("You input the wrong file")
   def ChangeMask(self,ccd,infile,row,x,y,r1,r2,rotang=None): #enlarge a mask which is already in mos2S002-bkg_region-sky.fits, the det file is automatically updated.
       fitsfile = fits.open("%s"%infile)
       hdu=fitsfile[1].data
       xcoor=fitsfile[1].columns[1].name
       detfile=infile.replace("sky","det")
       detx,dety=extractcoor(ccd,x,y)
       os.system("fpartab %s %s+1 X row=%s element=1"%(str(x),infile,row))
       os.system("fpartab %s %s+1 Y row=%s element=1"%(str(y),infile,row))
       os.system("fpartab %s %s+1 R row=%s element=1"%(str(r1),infile,row))
       os.system("fpartab %s %s+1 R row=%s element=2"%(str(r2),infile,row))
       if rotang is not None:
          os.system("fpartab %s %s+1 ROTANG row=%s element=1"%(str(rotang),infile,row))     
       os.system("fpartab %s %s+1 DETX row=%s element=1"%(detx,detfile,row))
       os.system("fpartab %s %s+1 DETY row=%s element=1"%(dety,detfile,row))
       os.system("fpartab %s %s+1 R row=%s element=1"%(r1,detfile,row))
       os.system("fpartab %s %s+1 R row=%s element=2"%(r2,detfile,row))
       if rotang is not None:
          os.system("fpartab %s %s+1 ROTANG row=%s element=1"%(rotang,detfile,row))
       fitsfile.close()
       os.system("make_mask inimage=%s-obj-im.fits inmask=%s-mask-im.fits outmask=%s-cheese.fits reglist=%s"%(ccd,ccd,ccd,infile)) 
   def Mergefits(self,newfile,oldfile):#merge two fits tables together
       os.system('fmerge "%s[1] %s[1]" %s - clobber=yes'%(oldfile,newfile,oldfile))
   def MakeMask(self,ccd,newcheesefile,reglist): #make a new mask image(newcheesefile) with a new region list(reglist)
       os.system("make_mask inimage=%s-obj-im.fits inmask=%s-mask-im.fits outmask=%s reglist=%s"%(ccd,ccd,newcheesefile,reglist))  

