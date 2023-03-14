from sys import *
import numpy as np
import pyfits
import commands
import matplotlib.pyplot as plt
import os,sys,string
from scipy.optimize import curve_fit
from esas_commands import esas_commands


ccd="mos1"
elow=2.5         #units in keV 
ehigh=12.0
binsize=60 #60 seconds is the default for binning lightcurves, see https://xmm-tools.cosmos.esa.int/external/sas/current/doc/espfilt.pdf
sig=2.    #when making the gaussian fit to the histogram, data within 2 sigma are filtered
histogram_bin=0.05   #the histogram is in uniform binsize in units of counts per second
fit_limit=1.4    #the range for fitting the histogram, which is the (hightest count rates from the histogram)*fit_limit and the same value divided by fit_limit 
legendsize=10

if ccd=="mos1":
   CCD="mos1S001"
elif ccd=="mos2":
   CCD="mos2S002"


inputfile="%s-ori.fits"%CCD   #product of emchain/epchain, which is an unfiltered event list
a=esas_commands()
LC=a.Lightcurves("FOV",inputfile,elow,ehigh)   #produce FOV LCs from the inputfile in the energy range:[elow,ehigh] 
LCcorn=a.Lightcurves("corner",inputfile,elow,ehigh) #produce corner LCs


if os.path.isfile(str(LC)) == True:
   print "The LCs for the FOV and corner are produced: ",LC,"and",LCcorn
else:
   print "The LCs for the FOV and corner are NOT produced. Please check."
   exit()




def func(x, a, x0, sigma):  #gaussian function
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
  
FOVLC = pyfits.open(LC)[1].data
rate=FOVLC['RATE'];
err=FOVLC['ERROR'];
time=FOVLC['TIME']

CORNLC= pyfits.open(LCcorn)[1].data
ratec=CORNLC['RATE'];
errc=CORNLC['ERROR'];
timec=CORNLC['TIME']

subListRate = [rate[n:n+binsize] for n in range(0, len(rate), binsize)]  #divide the rate, error and time into different sublist of time bin=binsize
subListErr = [err[n:n+binsize] for n in range(0, len(err), binsize)]
subListTime = [time[n:n+binsize] for n in range(0, len(time), binsize)]

subListRateC = [ratec[n:n+binsize] for n in range(0, len(ratec), binsize)]  #the same for the corner LC
subListErrC = [errc[n:n+binsize] for n in range(0, len(errc), binsize)]


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


print "best fit parameters (x0 and sigma): ","{:.2f}".format(x0),"{:.2f}".format(sigma)

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
#   print time1[0]-0.5,"       ",unwanted_time[0]-binsize/2.-0.5+time1[0]
   gtifile.write(str(time[0]-0.5)+"       "+str(unwanted_time[0]-binsize/2.-0.5+time[0])+"\n")

for i in range(0,len(unwanted_time)-1):
#   print  time1[0]+unwanted_time[i]+binsize/2.-0.5,"       ",time1[0]+unwanted_time[i+1]-binsize/2.-0.5
   gtifile.write(str(time[0]+unwanted_time[i]+binsize/2.-0.5)+"       "+str(time[0]+unwanted_time[i+1]-binsize/2.-0.5)+"\n")

if unwanted_time[-1]!=mid_time[-1]: 
#   print time1[0]+unwanted_time[len(unwanted_time)-1]+binsize/2.-0.5,"       ",time1[len(time1)-1]-0.5
   gtifile.write(str(time[0]+unwanted_time[-1]+binsize/2.-0.5)+"       "+str(time[-1]-0.5)+"\n")

gtifile.close()

axarr[0].text(evt_Lbound*0.3, y_peak*0.95, "sigma="+str("{:.2f}".format(sigma)),fontsize=legendsize+2,color="k")
axarr[0].text(evt_Lbound*0.3, y_peak*0.85, "x0="+str("{:.2f}".format(x0)),fontsize=legendsize+2,color="k")
axarr[1].errorbar(unwanted_time,unwanted_evts,yerr=unwanted_err,fmt="ko",markersize=1,label="unwanted events")
axarr[2].errorbar(unwanted_time,unwanted_evts_C,yerr=unwanted_err,fmt="ko",markersize=1,label="unwanted events")
axarr[0].legend(loc='upper right',fontsize=legendsize, frameon=False)
axarr[1].legend(loc='upper right',fontsize=legendsize, frameon=False)
axarr[2].legend(loc='upper right',fontsize=legendsize, frameon=False)

if os.path.isfile("%s_%s_%s_gti_ED.png"%(CCD,elow,ehigh)) == True:
   os.system("rm %s_%s_%s_gti_ED.png"%(CCD,elow,ehigh))
plt.savefig("%s_%s_%s_gti_ED.png"%(CCD,elow,ehigh))
#plt.show()

#creation of the good time interval file for making the filtered event list - mos1S001_clean.fits
final_product=a.make_Clean(CCD)

if os.path.isfile(str(final_product)) == True:
   print "The filtered clean event list is created: ",final_product
else:
   print "The filtered clean event list is NOT created. Please check."
   exit()
