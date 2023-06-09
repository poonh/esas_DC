from sys import *
import numpy as np
import shlex,subprocess
import os,sys,string
from residual_sp_commands import residual_sp_commands

mos1="mos1S001"
mos2="mos2S002"
pn="pnS003"

a=residual_sp_commands()

#use the combined exposure to calculate the area. Any pixels with a value=0 is not counted
mos1exp="comb-exp-im-400-2300-mos1.fits"
mos2exp="comb-exp-im-400-2300-mos2.fits"
pnexp="comb-exp-im-400-2300-pn.fits"

elow,ehigh=8,12   #lower and upper energy range in keV 

inner,outer=10,12 #inner and outer annuli in arcmin

mos1img=a.makeImage(mos1,elow,ehigh)  #make a FOV image in the energy of choice
mos2img=a.makeImage(mos2,elow,ehigh)
pnimg=a.makeImage(pn,elow,ehigh)

mos1rate=a.SB(mos1,inner,outer,mos1img,mos1exp) #calculate the surface brightness in the area of choice; the expousre map is only for area calculation; coordinate conversion is done here. Make sure "ecoordconv" works
mos2rate=a.SB(mos2,inner,outer,mos2img,mos2exp)
pnrate=a.SB(pn,inner,outer,pnimg,pnexp)

mos1cornrate=a.cornerspec(mos1,elow,ehigh) #create a corner spectrum in the energy of choice. The file mos1S001-corn.fits is needed
mos2cornrate=a.cornerspec(mos2,elow,ehigh)
pncornrate=a.cornerspec(pn,elow,ehigh)

def checkstate(ratio):
    if ratio < 1.15:
       state="not"
    elif ratio > 1.15 and ratio < 1.3:
       state="slightly"
    elif ratio > 1.3 and ratio < 1.5:
       state="very"
    elif ratio > 1.5:
       state="extremely"
    return state

ratiomos1=mos1rate/mos1cornrate
ratiomos2=mos2rate/mos2cornrate
ratiopn=pnrate/pncornrate

statemos1=checkstate(ratiomos1)
statemos2=checkstate(ratiomos2)
statepn=checkstate(ratiopn)

print("mos1 ratio="+"{:.2f}".format(ratiomos1)+". It is "+statemos1+" contaminated")
print("mos2 ratio="+"{:.2f}".format(ratiomos2)+". It is "+statemos2+" contaminated")
print("pn ratio="+"{:.2f}".format(ratiopn)+". It is "+statepn+" contaminated")




