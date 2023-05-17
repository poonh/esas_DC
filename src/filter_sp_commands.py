import os,sys,string
import shlex,subprocess
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#           expression="(PATTERN<=12)&&!((DETX,DETY) in BOX(10167,13005,3011,6575,0))&&(((FLAG & 0x766a0f63)==0)||((FLAG & 0x766a0763) == 0))&&(CCDNR==%s)&&(PI in [%s:%s])"%(ccdNum,int(elow*1000),int(ehigh*1000))

def func(x, a, x0, sigma):  #gaussian function
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def filteredCornEvts(CCD,outputfile,ccdNum,elow,ehigh): #make filtered event list using mos1S001-corn.fits of ccd chip 2-7
    expression="(PATTERN<=12)&&((FLAG & 0x766a0f63)==0)&&(CCDNR==%s)&&(PI in [%s:%s])"%(ccdNum,int(elow*1000),int(ehigh*1000))
    os.system("evselect table=%s-corn.fits withfilteredset=yes expression='%s' filtertype=expression filteredset=%s &> tmp_command.csh"%(CCD,expression,outputfile))

class filter_sp_commands():  #only for SAS19.0.0
   def write_ccd_anol(self,infile): #write the good ccd list to InputFiles/goodccdlist_mos1/mos2.txt, which will be used for mos-spectra.
       if os.path.isdir("InputFiles")==False:
          os.system("mkdir InputFiles")
       Infile=file(infile,"r").readlines() #open for file
       Outfilemos1=file("InputFiles/goodccdlist_mos1.dat","w")
       Outfilemos2=file("InputFiles/goodccdlist_mos2.dat","w")
       Outfilemos1.write("mos1ccd1: 1\n")
       Outfilemos2.write("mos2ccd1: 1\n")
       for i in range(0,6):
           if Infile[i].find("****")>=0 or Infile[i].find("####")>=0:
              Outfilemos1.write("mos1ccd%s: 0\n"%(i+2))
           else:
              Outfilemos1.write("mos1ccd%s: 1\n"%(i+2))

       for i in range(6,12):
           if Infile[i].find("****")>=0 or Infile[i].find("####")>=0:
              Outfilemos2.write("mos2ccd%s: 0\n"%(i-4))
           else:
              Outfilemos2.write("mos2ccd%s: 1\n"%(i-4))
       Outfilemos1.close()
       Outfilemos2.close()
       if os.path.isfile("InputFiles/goodccdlist_mos1.dat")==True:
           os.system("cat InputFiles/goodccdlist_mos1.dat")
       else:
           print("Inputfiles/goodccdlist_mos1.dat is not produced")
       if os.path.isfile("InputFiles/goodccdlist_mos2.dat")==True:
           os.system("cat InputFiles/goodccdlist_mos2.dat")
       else:
           print("Inputfiles/goodccdlist_mos2.dat is not produced")

   def Lightcurves(self,choice,inputfile,elow,ehigh): #make FOV or corner LCs,choice="FOV"or "corner",log file written in tmp_command.csh
      CCD=inputfile.split("-ori.fits")[0]
      if choice!="FOV":
         if choice!="corner":
             print("Please enter 'FOV' or 'corner' as your choice")
             exit()
      if choice=="FOV":         
         outputfile="%s-LC-%s-%s.fits"%(CCD,str(elow),str(ehigh))
         if os.path.isfile(outputfile) == True:
            os.system("rm %s"%outputfile)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.txt"%CCD)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.fits"%CCD)
         os.system("evselect table=%s expression='(PATTERN<=12)&&(PI in [%s:%s])&&((FLAG & 0xfb0000) == 0)&&!((DETX,DETY) in BOX(10167,13005,3011,6575,0))' filtertype=expression rateset=%s timecolumn=TIME timebinsize=1 maketimecolumn=yes makeratecolumn=yes withrateset=yes &> tmp_command.csh"%(inputfile,int(elow*1000),int(ehigh*1000),outputfile))
#         return str(outputfile)
      elif choice=="corner": 
         outputfile="%s-LC-corn-%s-%s.fits"%(CCD,str(elow),str(ehigh))
         if os.path.isfile(outputfile) == True:
            os.system("rm %s"%outputfile)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.txt"%CCD)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.fits"%CCD)
         os.system("evselect table=%s expression='(PATTERN<=12)&&(PI in [%s:%s])&&(((FLAG & 0x766a0f63) == 0)||((FLAG & 0x766a0763) == 0))&&!((DETX,DETY) in BOX(13280,-306,6610,6599,0))&&!((DETX,DETY) in BOX(-13169,-105,6599,6599,0))&&((FLAG & 0x766a0f63) == 0)&&!(((DETX,DETY) in CIRCLE(100,-200,17700))||((DETX,DETY) in CIRCLE(834,135,17100))||((DETX,DETY) in CIRCLE(770,-803,17100))||((DETX,DETY) in BOX(-20,-17000,6500,500,0))||((DETX,DETY) in BOX(5880,-20500,7500,1500,10))||((DETX,DETY) in BOX(-5920,-20500,7500,1500,350))||((DETX,DETY) in BOX(-20,-20000,5500,500,0)))' filtertype=expression rateset=%s timecolumn=TIME timebinsize=1 maketimecolumn=yes makeratecolumn=yes withrateset=yes &> tmp_command.csh"%(inputfile,int(elow*1000),int(ehigh*1000),outputfile))
      if os.path.isfile(outputfile) == True:
           print(outputfile+" is produdced")
      else:
           print(outputfile+" is not produdced. Check tmp_command.csh.")
      return str(outputfile)


   def make_gti(self,CCD,LC,LCcorn,elow=2.5,ehigh=12.0,binsize=60,sig=2,histogram_bin=0.05,fit_limit=1.4,pltshow="True"):
#From the FOV and corner LCs, you make a gaussian fit and select the GTI for the production of mos1S001-clean.fits
#LC and LCcorn are the names of the FOV and corner LCs extracted in the same elow and ehigh
#elow and ehigh in keV
#binsize in seconds
#sig=significance=2sigma . It is the gti extraction range, which is 2 sigma within the gaussian fit
#histogram_bin=0.05,the histogram is in uniform binsize in units of counts per second
#fit_limit: the range for fitting the histogram, which is the (hightest count rates from the histogram)*fit_limit and the same value divided by fit_limit 
      legendsize=10
      FOVLC = fits.open(LC)[1].data
      rate=FOVLC['RATE'];
      err=FOVLC['ERROR'];
      time=FOVLC['TIME']

      CORNLC= fits.open(LCcorn)[1].data
      ratec=CORNLC['RATE'];
      errc=CORNLC['ERROR'];
      timec=CORNLC['TIME']

      subListRate = [rate[n:n+binsize] for n in range(0, len(rate), binsize)]  #divide the rate, error and time into different sublist of time bin=binsize
      subListErr = [err[n:n+binsize] for n in range(0, len(err), binsize)]
      subListTime = [time[n:n+binsize] for n in range(0, len(time), binsize)]

      subListRateC = [ratec[n:n+binsize] for n in range(0, len(ratec), binsize)]  #the same for the corner LC
      subListErrC = [errc[n:n+binsize] for n in range(0, len(errc), binsize)]
      """
      subListRate.append(np.concatenate((subListRate[-2],subListRate[-1]))) #combine the last two sublists together to make a new one, then remove these two sublists
      subListRate=np.delete(subListRate,-3)
      subListRate=np.delete(subListRate,-2)

      subListRateC.append(np.concatenate((subListRateC[-2],subListRateC[-1]))) #combine the last two sublists together to make a new one, then remove these two sublists
      subListRateC=np.delete(subListRateC,-3)
      subListRateC=np.delete(subListRateC,-2)

      subListErr.append(np.concatenate((subListErr[-2],subListErr[-1])))
      subListErr=np.delete(subListErr,-3)
      subListErr=np.delete(subListErr,-2)

      subListErrC.append(np.concatenate((subListErrC[-2],subListErrC[-1])))
      subListErrC=np.delete(subListErrC,-3)
      subListErrC=np.delete(subListErrC,-2)

      subListTime.append(np.concatenate((subListTime[-2],subListTime[-1])))
      subListTime=np.delete(subListTime,-3)
      subListTime=np.delete(subListTime,-2)
      """
      mean_cr,mean_err,mid_time=[],[],[]
      mean_crC,mean_errC=[],[]

      for i in range(0,len(subListRate)-1):   #calculate the rate of the smoothed LC in units of counts per second
          tmp_rate=np.cumsum(subListRate[i])
          tmp_ratec=np.cumsum(subListRateC[i])
          tmp_err=np.cumsum(subListErr[i]**2) 
          tmp_errc=np.cumsum(subListErrC[i]**2)       
          mean_cr.append(tmp_rate[-1]/float(binsize))
          mean_crC.append(tmp_ratec[-1]/float(binsize))
          mean_err.append(tmp_err[-1]**0.5/float(binsize))
          mean_errC.append(tmp_errc[-1]**0.5/float(binsize))
          mid_time.append(i*binsize+binsize/2.)

      tmp_rate=np.cumsum(subListRate[-1])     #the mean count rate of the last sub list, which has a duation larger than the binsize
      tmp_rateC=np.cumsum(subListRateC[-1])
      tmp_err=np.cumsum(subListErr[-1]**2)
      tmp_errC=np.cumsum(subListErrC[-1]**2)
      duration= subListTime[-1][-1]-subListTime[-1][0]

      mean_cr.append(tmp_rate[-1]/float(duration))
      mean_crC.append(tmp_rateC[-1]/float(duration))
      mean_err.append(tmp_err[-1]**0.5/float(duration))
      mean_errC.append(tmp_errC[-1]**0.5/float(duration))
      mid_time.append((len(subListRate)-1)*binsize+(subListTime[-1][-1]-subListTime[-1][-0])/2.)

      fig, axarr = plt.subplots(3, 1,figsize=(16,13))

      axarr[0].hist(mean_cr,bins=np.arange(min(mean_cr), max(mean_cr) + histogram_bin, histogram_bin))
      axarr[1].errorbar(mid_time,mean_cr,yerr=mean_err,fmt="o",markersize=1,label="extracted events")
      axarr[2].errorbar(mid_time,mean_crC,yerr=mean_errC,fmt="o",markersize=1,label="extracted events")
      axarr[0].tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
      axarr[0].tick_params(length=5, width=1, which='major', direction='in', right=True, top=True) 
      axarr[1].tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
      axarr[1].tick_params(length=5, width=1, which='major', direction='in', right=True, top=True) 
      axarr[2].tick_params(length=5, width=1, which='minor', direction='in', right=True, top=True)
      axarr[2].tick_params(length=5, width=1, which='major', direction='in', right=True, top=True) 
      axarr[0].set_title ("Count Rate Histogram")
      axarr[1].set_title ("FOV LC")
      axarr[2].set_title ("corner LC")
      axarr[0].set_ylabel("N",fontsize=legendsize) 
      axarr[1].set_xlabel("Time(s)",fontsize=legendsize)
      axarr[1].set_ylabel("Count Rate (counts/s)",fontsize=legendsize)
      axarr[2].set_xlabel("Time(s)",fontsize=legendsize)
      axarr[2].set_ylabel("Count Rate (counts/s)",fontsize=legendsize)

#The following retrieves the x, y values from the histogram and then make a gaussian fit
      ax = plt.gca() 
      p = axarr[0].patches  
      heights = [patch.get_height() for patch in p] #y value of the histogram
      x_value = [patch.get_xy() for patch in p]   #x value of the histogram
      x=[row[0]+histogram_bin/2. for row in x_value] #the mid value
      y_peak=np.max(heights)  
      index= heights.index(y_peak)
      max_x=x[index]

      upper_lim,lower_lim=x[-1],x[0] #max and min x range

#find out the exact fitting range by cutting off the highest and lowest tail
      for i in range(0,len(x)):
          if x[i]>fit_limit*max_x:
             upper_lim=x[i-1]
             break
      for i in range(0,len(x)):
          if x[i]<max_x/fit_limit and i <len(x)-1:
             lower_lim=x[i+1]
             break
          else:
             lower_lim=x[i]
             break
      
#extract the data fulfulling the criteria for gaussian fit
      x_range_tmp = np.asarray(x)[np.asarray(x)<= upper_lim]
      x_range = x_range_tmp[x_range_tmp>=lower_lim]

      y_range_tmp = np.asarray(heights)[np.asarray(x)<= upper_lim]
      y_range = y_range_tmp[x_range_tmp>=lower_lim]

      popt, pcov = curve_fit(func, x_range, y_range)
      x0,sigma = popt[1],np.abs(popt[2])

#      print("best fit parameters (x0 and sigma): ","{:.2f}".format(x0),"{:.2f}".format(sigma))

      evt_Lbound = x0-sig*sigma
      evt_Rbound = x0+sig*sigma

#produce the best fit on the plot
      x_data = np.arange(int(np.min(x)*100/1.5),int(np.max(x)*100)*1.1)/100.
      ym = func(x_data, popt[0], popt[1], popt[2])

      axarr[0].plot(x_data, ym, c='r', label='Best fit')
      axarr[0].set_xlim(0,x0+np.max(x)*1.1)
      axarr[0].axvline(upper_lim+histogram_bin/2., ymin=0., ymax = max(ym)*1.2, linewidth=2, color='b',label='Fitting range')
      axarr[0].axvline(lower_lim-histogram_bin/2., ymin=0., ymax = max(ym)*1.2, linewidth=2, color='b')
      axarr[0].axvline(evt_Lbound, ymin=0., ymax = max(ym)*1.2, linewidth=2, color='g',label='evts extraction range')
      axarr[0].axvline(evt_Rbound, ymin=0., ymax = max(ym)*1.2, linewidth=2, color='g')
      axarr[0].tick_params(length=5, width=1, direction='in', right=True, top=True)
      axarr[1].tick_params(length=5, width=1, direction='in', right=True, top=True)
      axarr[2].tick_params(length=5, width=1, direction='in', right=True, top=True)

#Remove the unwanted events that do not pass the filtering criteria
      unwanted_evts,unwanted_evts_C,unwanted_time,unwanted_err,unwanted_errC=[],[],[],[],[]

      for i in range(0,len(mean_cr)):
          if mean_cr[i]<evt_Lbound or mean_cr[i]>evt_Rbound:
             unwanted_evts.append(mean_cr[i])
             unwanted_evts_C.append(mean_crC[i])
             unwanted_err.append(mean_err[i])
             unwanted_errC.append(mean_errC[i])
             unwanted_time.append(mid_time[i])

      gtifile=open("%s_gti_ED.txt"%CCD,"w") 

      if unwanted_time[0]!=mid_time[0]: 
         gtifile.write(str(time[0]-0.5)+"       "+str(unwanted_time[0]-binsize/2.-0.5+time[0])+"\n")
      for i in range(0,len(unwanted_time)-1):
         gtifile.write(str(time[0]+unwanted_time[i]+binsize/2.-0.5)+"       "+str(time[0]+unwanted_time[i+1]-binsize/2.-0.5)+"\n")
      if unwanted_time[-1]!=mid_time[-1]: 
         gtifile.write(str(time[0]+unwanted_time[-1]+binsize/2.-0.5)+"       "+str(time[-1]-0.5)+"\n")
      gtifile.close()
      axarr[0].text(evt_Lbound*0.3, y_peak*0.95, "sigma="+str("{:.2f}".format(sigma)),fontsize=legendsize+2,color="k")
      axarr[0].text(evt_Lbound*0.3, y_peak*0.85, "x0="+str("{:.2f}".format(x0)),fontsize=legendsize+2,color="k")
      axarr[1].errorbar(unwanted_time,unwanted_evts,yerr=unwanted_err,fmt="ko",markersize=1,label="unwanted events")
      axarr[2].errorbar(unwanted_time,unwanted_evts_C,yerr=unwanted_err,fmt="ko",markersize=1,label="unwanted events")
      axarr[0].legend(loc='upper right',fontsize=legendsize, frameon=False)
      axarr[1].legend(loc='upper right',fontsize=legendsize, frameon=False)
      axarr[2].legend(loc='upper right',fontsize=legendsize, frameon=False)
      plt.savefig("%s_%s_%s_gti_ED.png"%(CCD,elow,ehigh))
      print("The good time interval file, %s_gti_ED.txt, is produced"%CCD)
      if pltshow == "True":
         plt.show()

   def txt2fits(self,infiletxt,outfilefits):
       os.system('ftcreate colname.lis %s %s extname = "STDGTI" clobber=yes'%(infiletxt,outfilefits))
       if os.path.isfile(outfilefits) == True:
          print(outfilefits+" is produced")
       else:
          print(outfilefits+" is not produced")

   def make_Clean(self,CCD,gtifile): #make mos1S001-clean.fits
         final_product="%s-clean_ED.fits"%CCD
         if os.path.isfile("%s-clean_ED.fits"%CCD) == True:
            os.system("rm %s-clean_ED.fits"%CCD)
         os.system("evselect table=%s-ori.fits filteredset=%s expression='(PATTERN<=12)&&GTI(%s,TIME)&&(((FLAG & 0x766a0f63)==0)||((FLAG & 0x766a0763) == 0))' filtertype=expression &> tmp_command.csh"%(CCD,final_product,gtifile))
         if os.path.isfile(final_product) == True:
            print(final_product+" is produced")
         else:
            print(final_product+" is not produced")
#         return final_product
   def make_corn_evtlist(self,CCD,outputfile): #make mos1S001-corn.fits
       if CCD.find("mos1")>=0:
          os.system("evselect table=%s-clean.fits withfilteredset=yes expression='(((FLAG & 0x766a0f63) == 0)||((FLAG & 0x766a0763) == 0))&&!(((DETX,DETY) in CIRCLE(100,-200,17700))||((DETX,DETY) in CIRCLE(834,135,17100))||((DETX,DETY) in CIRCLE(770,-803,17100))||((DETX,DETY) in BOX(-20,-17000,6500,500,0))||((DETX,DETY) in BOX(5880,-20500,7500,1500,10))||((DETX,DETY) in BOX(-5920,-20500,7500,1500,350))||((DETX,DETY) in BOX(-20,-20000,5500,500,0))||((DETX,DETY) in BOX(-12900,16000,250,4000,0))||((DETX,DETY) in BOX(80,18600,150,1300,0)))||((DETX,DETY) in BOX(-10,-18800,125,1500,0))' filteredset=%s filtertype=expression &> tmp_command.csh"%(CCD,outputfile))
       elif CCD.find("mos2")>=0:
          os.system("evselect table=%s-clean.fits:EVENTS withfilteredset=yes expression='((FLAG & 0x766a0f63) == 0)&&!(CIRCLE(435,1006,17100,DETX,DETY)||CIRCLE(-34,68,17700,DETX,DETY)||BOX(-20,-17000,6500,500,0,DETX,DETY)||BOX(5880,-20500,7500,1500,10,DETX,DETY)||BOX(-5920,-20500,7500,1500,350,DETX,DETY)||BOX(-20,-20000,5500,500,0,DETX,DETY))' filteredset=%s filtertype=expression &> tmp_command.csh"%(CCD,outputfile))

   def problemccdcheck(self,CCD,elow1,ehigh1,elow2,ehigh2): #reproduce command.csh after mos-filter, which indicates the anomalous ccd chips
       hd1=" ccd "+"||"+" total      "+"||"+" total      "+"||"+" hardness"
       hd2=" chip "+"||"+" events(LE) "+"||"+" events(HE) "+"||"+" ratio"
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
          filteredCornEvts(CCD,lowEoutput,i,elow1,ehigh1)   #make filtered event list in the energy range[elow1,ehigh2]
          filteredCornEvts(CCD,highEoutput,i,elow2,ehigh2)
          lowE_list = fits.open(lowEoutput)[1].data  #open the total number of events in the file
          highE_list = fits.open(highEoutput)[1].data
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

   def CornImage(self,CCD,elow,ehigh): #make corner image in a certain energy band
        outputfile="%s-corn-image-%s-%s.fits"%(CCD,str(elow),str(ehigh))
        if os.path.isfile("%s"%outputfile) == True:
            os.system("rm %s"%outputfile)
        os.system("evselect table=%s-corn.fits withimageset=yes imageset=%s xcolumn='DETX' ximagesize=780 filtertype=expression ignorelegallimits=yes expression='PI in [%s:%s]' ximagemax=19500 ximagemin=-19499 ycolumn='DETY' yimagesize=780 yimagemax=19500 yimagemin=-19499 imagebinning=imageSize &> tmp_command.csh"%(CCD,outputfile,int(elow*1000),int(ehigh*1000)))
        return outputfile
   

        
