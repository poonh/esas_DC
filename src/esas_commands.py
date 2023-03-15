import os,sys
import commands


class esas_commands():
   def Lightcurves(self,choice,inputfile,elow,ehigh): #make FOV or corner LCs
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
         return str(outputfile)
      elif choice=="corner": 
         outputfile="%s-LC-corn-%s-%s.fits"%(CCD,str(elow),str(ehigh))
         if os.path.isfile(outputfile) == True:
            os.system("rm %s"%outputfile)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.txt"%CCD)
         if os.path.isfile("%s_gti_ED.txt"%CCD) == True:
            os.system("rm %s_gti_ED.fits"%CCD)
         os.system("evselect table=%s expression='(PATTERN<=12)&&(PI in [%s:%s])&&(((FLAG & 0x766a0f63) == 0)||((FLAG & 0x766a0763) == 0))&&!((DETX,DETY) in BOX(13280,-306,6610,6599,0))&&!((DETX,DETY) in BOX(-13169,-105,6599,6599,0))&&((FLAG & 0x766a0f63) == 0)&&!(((DETX,DETY) in CIRCLE(100,-200,17700))||((DETX,DETY) in CIRCLE(834,135,17100))||((DETX,DETY) in CIRCLE(770,-803,17100))||((DETX,DETY) in BOX(-20,-17000,6500,500,0))||((DETX,DETY) in BOX(5880,-20500,7500,1500,10))||((DETX,DETY) in BOX(-5920,-20500,7500,1500,350))||((DETX,DETY) in BOX(-20,-20000,5500,500,0)))' filtertype=expression rateset=%s timecolumn=TIME timebinsize=1 maketimecolumn=yes makeratecolumn=yes withrateset=yes &> tmp_command.csh"%(inputfile,int(elow*1000),int(ehigh*1000),outputfile))
         return str(outputfile)
   def make_Clean(self,CCD): #make mos1S001-clean.fits
         final_product="%s-clean_ED.fits"%CCD
         if os.path.isfile("%s-clean_ED.fits"%CCD) == True:
            os.system("rm %s-clean_ED.fits"%CCD)
         os.system('ftcreate colname.lis %s_gti_ED.txt %s_gti_ED.fits extname = "STDGTI" clobber=yes'%(CCD,CCD))
         os.system("evselect table=%s-ori.fits filteredset=%s expression='(PATTERN<=12)&&GTI(%s_gti_ED.fits,TIME)&&(((FLAG & 0x766a0f63)==0)||((FLAG & 0x766a0763) == 0))' filtertype=expression &> tmp_command.csh"%(CCD,final_product,CCD))
         return final_product
   def make_corn_evtlist(self,CCD,outputfile): #make mos1S001-corn.fits
       if CCD.find("mos1")>=0:
          os.system("evselect table=%s-clean.fits withfilteredset=yes expression='(((FLAG & 0x766a0f63) == 0)||((FLAG & 0x766a0763) == 0))&&!(((DETX,DETY) in CIRCLE(100,-200,17700))||((DETX,DETY) in CIRCLE(834,135,17100))||((DETX,DETY) in CIRCLE(770,-803,17100))||((DETX,DETY) in BOX(-20,-17000,6500,500,0))||((DETX,DETY) in BOX(5880,-20500,7500,1500,10))||((DETX,DETY) in BOX(-5920,-20500,7500,1500,350))||((DETX,DETY) in BOX(-20,-20000,5500,500,0))||((DETX,DETY) in BOX(-12900,16000,250,4000,0))||((DETX,DETY) in BOX(80,18600,150,1300,0)))||((DETX,DETY) in BOX(-10,-18800,125,1500,0))' filteredset=%s filtertype=expression &> tmp_command.csh"%(CCD,outputfile))
       elif CCD.find("mos2")>=0:
          os.system("evselect table=%s-clean.fits:EVENTS withfilteredset=yes expression='((FLAG & 0x766a0f63) == 0)&&!(CIRCLE(435,1006,17100,DETX,DETY)||CIRCLE(-34,68,17700,DETX,DETY)||BOX(-20,-17000,6500,500,0,DETX,DETY)||BOX(5880,-20500,7500,1500,10,DETX,DETY)||BOX(-5920,-20500,7500,1500,350,DETX,DETY)||BOX(-20,-20000,5500,500,0,DETX,DETY))' filteredset=%s filtertype=expression &> tmp_command.csh"%(CCD,outputfile))
   def expr_ccdCheck(self,CCD,ccdNum,elow,ehigh):
       if CCD.find("mos1")>=0:
           expression="(PATTERN<=12)&&!((DETX,DETY) in BOX(10167,13005,3011,6575,0))&&(((FLAG & 0x766a0f63)==0)||((FLAG & 0x766a0763) == 0))&&(CCDNR==%s)&&(PI in [%s:%s])"%(ccdNum,int(elow*1000),int(ehigh*1000))
       elif CCD.find("mos2")>=0:
           expression="(PATTERN<=12)&&((FLAG & 0x766a0f63)==0)&&(CCDNR==%s)&&(PI in [%s:%s])"%(ccdNum,int(elow*1000),int(ehigh*1000))
       return expression
   def filteredCornEvts(self,CCD,outputfile,expression): #make filtered event list using mos1S001-corn.fits of ccd 2-7
        os.system("evselect table=%s-corn.fits withfilteredset=yes expression='%s' filtertype=expression filteredset=%s &> tmp_command.csh"%(CCD,expression,outputfile))
   def CornImage(self,CCD,elow,ehigh): #make corner image in a certain energy band
        outputfile="%s-corn-image-%s-%s.fits"%(CCD,str(elow),str(ehigh))
        if os.path.isfile("%s"%outputfile) == True:
            os.system("rm %s"%outputfile)
        os.system("evselect table=%s-corn.fits withimageset=yes imageset=%s xcolumn='DETX' ximagesize=780 filtertype=expression ignorelegallimits=yes expression='PI in [%s:%s]' ximagemax=19500 ximagemin=-19499 ycolumn='DETY' yimagesize=780 yimagemax=19500 yimagemin=-19499 imagebinning=imageSize &> tmp_command.csh"%(CCD,outputfile,int(elow*1000),int(ehigh*1000)))
        return outputfile
   
