import os,sys
from sys import *
import numpy as np
import commands
from spectral_fit_commands import spectral_fit_commands
from spectrum_commands import spectrum_commands

whichstep=argv[1]
#limit=4.6 #limit for defining cluster area in arcmin. You can take r500 from the MCXC catalogue as reference.

if whichstep=="4":
   inradius=input("input the inner radius in arcsecs: ")
   outradius=input("input the outer radius in arcsecs: ")
   if outradius<inradius:
      print("Your outer radius,%s, is smaller than inner radius,%s"%(outradius,inradius))
      exit()

if whichstep=="3":
   limit=input("input the limiting radius in arcsecs: ")
   

mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")

xc=float(commands.getoutput("grep 'PixX,PixY' log/Center.log | gawk '{print $3}'"))
yc=float(commands.getoutput("grep 'PixX,PixY' log/Center.log | gawk '{print $4}'"))

#r500=float(commands.getoutput("grep -i 'r500_mcxc_arcmin' r500_z.txt | gawk '{print $3}'"))

#print "r500: ",r500
elow_ori,ehigh_ori=0.4,2.3 #old energy range in keV, the files of which have already been produced
#elow_new,ehigh_new=0.7,10
#eloweV_new,ehigheV_new=int(0.7*1000),int(10*1000)

elow_new,ehigh_new=0.7,10.0 #the new energy range in keV
eloweV_new,ehigheV_new=int(elow_new*1000),int(ehigh_new*1000) 
a=spectral_fit_commands()
b=spectrum_commands()



mos1img,mos1back="%s-obj-im-det-700-10000.fits"%mos1,"%s-back-im-det-700-10000.fits"%mos1
mos2img,mos2back="%s-obj-im-det-700-10000.fits"%mos2,"%s-back-im-det-700-10000.fits"%mos2
pnimg,pnootimg,pnback="%s-obj-im-det-700-10000.fits"%pn,"%s-obj-im-det-700-10000-oot.fits"%pn,"%s-back-im-det-700-10000.fits"%pn

mos1x=commands.getoutput("grep 'mos1' log/Center.log | gawk '{print $5}'")
mos1y=commands.getoutput("grep 'mos1' log/Center.log | gawk '{print $6}'")

mos2x=commands.getoutput("grep 'mos2' log/Center.log | gawk '{print $5}'")
mos2y=commands.getoutput("grep 'mos2' log/Center.log | gawk '{print $6}'")

pnx=commands.getoutput("grep 'pn' log/Center.log | gawk '{print $5}'")
pny=commands.getoutput("grep 'pn' log/Center.log | gawk '{print $6}'")

ootfactor=float(commands.getoutput("grep -i 'OOT scale factor: ' log/17comb-pn.log | awk '{print $5}'"))
xcmos1,ycmos1=a.coordinateconversion("full_spectrum_mos1/%s"%mos1img,mos1x,mos1y)
xcmos2,ycmos2=a.coordinateconversion("full_spectrum_mos2/%s"%mos2img,mos2x,mos2y)
xcpn,ycpn=a.coordinateconversion("full_spectrum_pn/%s"%pnimg,pnx,pny)


#produce the files necessary for counts estimation in different annuli. The commands are extracted from command/FOV_imgspe_mos2.csh and written into full_spectrum_mos2/backgd_command_mos2.csh. This csh file is run automatically. Check this file and run the commands line by line if things go wrong.
if whichstep=="1":
   path_to_caldb="/net/cluster491/software/newton/caldb/esas"
   a.makeImage(mos1,elow_ori,ehigh_ori,elow_new,ehigh_new,path_to_caldb)
   a.makeImage(mos2,elow_ori,ehigh_ori,elow_new,ehigh_new,path_to_caldb)
   a.makeImage(pn,elow_ori,ehigh_ori,elow_new,ehigh_new,path_to_caldb)



if whichstep=="2":
   b.plotradial(["full_spectrum_mos1/comb-obj-im-400-2300-mos1.fits","full_spectrum_mos2/comb-obj-im-400-2300-mos2.fits","full_spectrum_pn/comb-obj-im-400-2300-pn.fits"],
               ["full_spectrum_mos1/comb-back-im-sky-400-2300-mos1.fits","full_spectrum_mos2/comb-back-im-sky-400-2300-mos2.fits","full_spectrum_pn/comb-back-im-sky-400-2300-pn.fits"],
                ["full_spectrum_mos1/%s-cheese.fits"%mos1,"full_spectrum_mos2/%s-cheese.fits"%mos2,"full_spectrum_pn/%s-cheese.fits"%pn],
                ["full_spectrum_mos1/comb-exp-im-400-2300-mos1.fits","full_spectrum_mos2/comb-exp-im-400-2300-mos2.fits","full_spectrum_pn/comb-exp-im-400-2300-pn.fits"],xc=xc,yc=yc,rmax=14.,Nbin=40,constant="deg",outfile="sb_0.4-2.3eV.png")


if whichstep=="2a":
   b.plotradial(["full_spectrum_mos1/comb-obj-im-%s-%s-mos1.fits"%(eloweV_new,ehigheV_new),"full_spectrum_mos2/comb-obj-im-%s-%s-mos2.fits"%(eloweV_new,ehigheV_new),"full_spectrum_pn/comb-obj-im-%s-%s-pn.fits"%(eloweV_new,ehigheV_new)],
                ["full_spectrum_mos1/comb-back-im-sky-%s-%s-mos1.fits"%(eloweV_new,ehigheV_new),"full_spectrum_mos2/comb-back-im-sky-%s-%s-mos2.fits"%(eloweV_new,ehigheV_new),"full_spectrum_pn/comb-back-im-sky-%s-%s-pn.fits"%(eloweV_new,ehigheV_new)],
                ["full_spectrum_mos1/%s-cheese.fits"%mos1,"full_spectrum_mos2/%s-cheese.fits"%mos2,"full_spectrum_pn/%s-cheese.fits"%pn],
                ["full_spectrum_mos1/comb-exp-im-%s-%s-mos1.fits"%(eloweV_new,ehigheV_new),"full_spectrum_mos2/comb-exp-im-%s-%s-mos2.fits"%(eloweV_new,ehigheV_new),"full_spectrum_pn/comb-exp-im-%s-%s-pn.fits"%(eloweV_new,ehigheV_new)],xc=xc,yc=yc,rmax=14.,Nbin=40,constant="deg",outfile="sb_0.7-10eV.png")

###

if whichstep=="3":
   inner=0
   region,counts,sn,sbratiom2=a.DefineRegionCounts(inner,limit,30,1500,xcmos2,ycmos2,"full_spectrum_mos2/%s"%mos2img,"full_spectrum_mos2/%s"%mos2back)
   b,c,sbratiom1=a.DefineRegionCounts(inner,limit,30,1500,xcmos1,ycmos1,"full_spectrum_mos1/%s"%mos1img,"full_spectrum_mos1/%s"%mos1back,checkregion="no",Region=region)
   d,e,sbratiopn=a.DefineRegionCounts(inner,limit,30,1500,xcpn,ycpn,"full_spectrum_pn/%s"%pnimg,"full_spectrum_pn/%s"%pnback,checkregion="no",Region=region,ootimg="full_spectrum_pn/%s"%pnootimg,ootscalefactor=ootfactor)



   print("source region(arcsec) counts(bkg subtracted)")
   print("                   mos1        mos2      pn")
   for i in range(0,len(region)-1):
       print("{:.2f}".format(region[i])+" - "+"{:.2f}".format(region[i+1])+": "+ " || "+"{:.2f}".format(b[i])+" || "+"{:.2f}".format(counts[i])+" || "+"{:.2f}".format(d[i]))
   print("S/N ratio")
   for i in range(0,len(region)-1):
       print("{:.2f}".format(region[i])+" - "+"{:.2f}".format(region[i+1])+": "+" || "+"{:.2f}".format(sn[i])+" || "+"{:.2f}".format(c[i])+" || "+"{:.2f}".format(e[i]))
   print("source/bkg ratio")
   for i in range(0,len(region)-1):
       print("{:.2f}".format(region[i])+" - "+"{:.2f}".format(region[i+1])+": "+" || "+"{:.2f}".format(sbratiom1[i])+" || "+"{:.2f}".format(sbratiom2[i])+" || "+"{:.2f}".format(sbratiopn[i]))



   cxbinner=(limit+60.)
   cxbouter=15.*60. #outer limit = 15 arcmin
   cxbregion,cxbcounts,cxbsn,cxbsbm2=a.DefineRegionCounts(cxbinner,cxbouter,30,1500,xcmos2,ycmos2,"full_spectrum_mos2/%s"%mos2img,"full_spectrum_mos2/%s"%mos2back)

   cxbcounts1,cxbsn1,cxbsbm1=a.DefineRegionCounts(cxbinner,cxbouter,30,1500,xcmos1,ycmos1,"full_spectrum_mos1/%s"%mos1img,"full_spectrum_mos1/%s"%mos1back,checkregion="no",Region=cxbregion)

   cxbcounts3,cxbsn3,cxbsbpn=a.DefineRegionCounts(cxbinner,cxbouter,30,1500,xcpn,ycpn,"full_spectrum_pn/%s"%pnimg,"full_spectrum_pn/%s"%pnback,checkregion="no",Region=cxbregion,ootimg="full_spectrum_pn/%s"%pnootimg,ootscalefactor=ootfactor)


   print("cxb region(arcsec) counts(bkg subtracted)")
   print("                   mos1        mos2      pn")
   for i in range(0,len(cxbregion)-1):
       print("{:.2f}".format(cxbregion[i])+" - "+"{:.2f}".format(cxbregion[i+1])+": "+"{:.2f}".format(cxbcounts1[i])+" || "+"{:.2f}".format(cxbcounts[i])+" || "+"{:.2f}".format(cxbcounts3[i]))
   print("S/N ratio")
   for i in range(0,len(cxbregion)-1):
       print("{:.2f}".format(cxbregion[i])+" - "+"{:.2f}".format(cxbregion[i+1])+": "+"{:.2f}".format(cxbsn1[i])+" || "+"{:.2f}".format(cxbsn[i])+" || "+"{:.2f}".format(cxbsn3[i]))
   print("source/bkg ratio")
   for i in range(0,len(cxbregion)-1):
       print("{:.2f}".format(cxbregion[i])+" - "+"{:.2f}".format(cxbregion[i+1])+": "+"{:.2f}".format(cxbsbm1[i])+" || "+"{:.2f}".format(cxbsbm2[i])+" || "+"{:.2f}".format(cxbsbpn[i]))


   source,cxb=[],[]

   for i in range(0,len(region)):
       source.append("source")

   for i in range(0,len(cxbregion)):
       cxb.append("cxb")

   TotReg=region+cxbregion
   TotSrc=source+cxb

   AllRegion=np.vstack((TotReg,TotSrc)).T

   lasttolimitratio=(region[-1]+region[-2])*0.5/float(limit)
   print("last region(mid point) to limiting radius ratio: ","{:.2f}".format(lasttolimitratio))
   print("The result is written in spectral_region.dat")

   a.writeannulus("spectral_region_tmp.dat",AllRegion)



if whichstep=="4":
   m1cts,m1sn,sbratiom1 = a.howmanyCounts(xcmos1,ycmos1,"full_spectrum_mos1/%s"%mos1img,"full_spectrum_mos1/%s"%mos1back,inradius/60.,outradius/60.)
   m2cts,m2sn,sbratiom2 = a.howmanyCounts(xcmos2,ycmos2,"full_spectrum_mos2/%s"%mos2img,"full_spectrum_mos2/%s"%mos2back,inradius/60.,outradius/60.)
   pncts,pnsn,sbratiopn = a.howmanyCounts(xcpn,ycpn,"full_spectrum_pn/%s"%pnimg,"full_spectrum_pn/%s"%pnback,inradius/60.,outradius/60.,ootimg="full_spectrum_pn/%s"%pnootimg,ootscalefactor=ootfactor)
   print("source region(arcsec) counts")
   print("                  mos1      mos2      pn")
   print(str(inradius)+" - "+str(outradius)+" :  ||  "+"{:.2f}".format(m1cts)+ " || "+"{:.2f}".format(m2cts)+" || "+"{:.2f}".format(pncts))
   print("source region(arcsec) S/N")
   print("                  mos1      mos2      pn")
   print(str(inradius)+" - "+str(outradius)+" :  ||  "+"{:.2f}".format(m1sn)+ " || "+"{:.2f}".format(m2sn)+" || "+"{:.2f}".format(pnsn))
   print("source region(arcsec) source/bkg")
   print("                  mos1      mos2      pn")
   print(str(inradius)+" - "+str(outradius)+" :  ||  "+"{:.2f}".format(sbratiom1)+ " || "+"{:.2f}".format(sbratiom2)+" || "+"{:.2f}".format(sbratiopn))




