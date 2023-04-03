setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

atthkgen atthkset=atthk.fits timestep=1

evselect table=mos1S001-clean.fits:EVENTS withfilteredset=yes expression='(PATTERN<=12)&&(FLAG == 0)&&(PI in [400:2300])' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=mos1S001-obj-im.fits squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes verbosity=1

evselect table=mos2S002-clean.fits:EVENTS withfilteredset=yes expression='(PATTERN<=12)&&(FLAG == 0)&&(PI in [400:2300])' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=mos2S002-obj-im.fits squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes verbosity=1
   
evselect table=pnS003-clean.fits:EVENTS withfilteredset=yes expression='(PATTERN<=12)&&(FLAG == 0)&&(PI in [400:2300])&&(DETY in [-16510:14345])' filtertype=expression imagebinning='imageSize' imagedatatype='Int32' imageset=pnS003-obj-im.fits squarepixels=yes ignorelegallimits=yes withxranges=yes withyranges=yes xcolumn='X' ximagesize=900 ximagemax=48400 ximagemin=3401 ycolumn='Y' yimagesize=900 yimagemax=48400 yimagemin=3401 updateexposure=yes filterexposure=yes verbosity=1


#edetect_chain
eexpmap imageset=pnS003-obj-im.fits attitudeset=atthk.fits eventset=pnS003-clean.fits expimageset=pnS003-obj-imexp.fits withdetcoords=no withvignetting=yes usefastpixelization=no usedlimap=no attrebin=4 pimin=400 pimax=2300  

emask expimageset=pnS003-obj-imexp.fits detmaskset=pnS003-obj-immask.fits detmasktable=MASK threshold1=0.3 threshold2=0.5  

eexpmap imageset=mos1S001-obj-im.fits attitudeset=atthk.fits eventset=mos1S001-clean.fits expimageset=mos1S001-obj-imexp.fits withdetcoords=no withvignetting=yes usefastpixelization=no usedlimap=no attrebin=4 pimin=400 pimax=2300  

emask expimageset=mos1S001-obj-imexp.fits detmaskset=mos1S001-obj-immask.fits detmasktable=MASK threshold1=0.3 threshold2=0.5 

eexpmap imageset=mos2S002-obj-im.fits attitudeset=atthk.fits eventset=mos2S002-clean.fits expimageset=mos2S002-obj-imexp.fits withdetcoords=no withvignetting=yes usefastpixelization=no usedlimap=no attrebin=4 pimin=400 pimax=2300 

emask expimageset=mos2S002-obj-imexp.fits detmaskset=mos2S002-obj-immask.fits detmasktable=MASK threshold1=0.3 threshold2=0.5  

eboxdetect detmasksets='pnS003-obj-immask.fits mos1S001-obj-immask.fits mos2S002-obj-immask.fits' withdetmask=yes expimagesets='pnS003-obj-imexp.fits mos1S001-obj-imexp.fits mos2S002-obj-imexp.fits' withexpimage=yes nruns=3 likemin=15 boxsize=5 withimagebuffersize=no imagebuffersize=640 ecf='3.2 1.2 1.2' imagesets='pnS003-obj-im.fits mos1S001-obj-im.fits mos2S002-obj-im.fits' boxlistset=eboxlist_l.fits withoffsets=no mergedlistset=mergedlist.fits hrdef='1 2 2 3 3 4' pimin='400 400 400' pimax='2300 2300 2300' obsmode=pointing  

esplinemap scut=0.01 mlmin=1 nsplinenodes=20 excesssigma=4 nfitrun=3 idband=1 boxlistset=eboxlist_l.fits imageset=pnS003-obj-im.fits expimageset=pnS003-obj-imexp.fits withexpimage=yes expimageset2=pnS003-obj-imexpnovig.fits withexpimage2=no detmaskset=pnS003-obj-immask.fits withdetmask=yes bkgimageset=pnS003-obj-imbkg.fits cheeseimageset=pnS003-obj-imcheese.fits withcheese=no cheesemaskset=cheesemask.fits withcheesemask=no ooteventset=pnS003-clean-oot.fits withootset=yes pimin=400 pimax=2300 fitmethod=spline snrmin=30 smoothsigma=15  

esplinemap scut=0.01 mlmin=1 nsplinenodes=20 excesssigma=4 nfitrun=3 idband=1 boxlistset=eboxlist_l.fits imageset=mos1S001-obj-im.fits expimageset=mos1S001-obj-imexp.fits withexpimage=yes expimageset2=mos1S001-obj-imexpnovig.fits withexpimage2=no detmaskset=mos1S001-obj-immask.fits withdetmask=yes bkgimageset=mos1S001-obj-imbkg.fits cheeseimageset=mos1S001-obj-imcheese.fits withcheese=no cheesemaskset=cheesemask.fits withcheesemask=no ooteventset=ootevents.fits withootset=no pimin=400 pimax=2300 fitmethod=spline snrmin=30 smoothsigma=15  

esplinemap scut=0.01 mlmin=1 nsplinenodes=20 excesssigma=4 nfitrun=3 idband=1 boxlistset=eboxlist_l.fits imageset=mos2S002-obj-im.fits expimageset=mos2S002-obj-imexp.fits withexpimage=yes expimageset2=mos2S002-obj-imexpnovig.fits withexpimage2=no detmaskset=mos2S002-obj-immask.fits withdetmask=yes bkgimageset=mos2S002-obj-imbkg.fits cheeseimageset=mos2S002-obj-imcheese.fits withcheese=no cheesemaskset=cheesemask.fits withcheesemask=no ooteventset=ootevents.fits withootset=no pimin=400 pimax=2300 fitmethod=spline snrmin=30 smoothsigma=15 

eboxdetect bkgimagesets='pnS003-obj-imbkg.fits mos1S001-obj-imbkg.fits mos2S002-obj-imbkg.fits' usemap=yes usematchedfilter=no detmasksets='pnS003-obj-immask.fits mos1S001-obj-immask.fits mos2S002-obj-immask.fits' withdetmask=yes expimagesets='pnS003-obj-imexp.fits mos1S001-obj-imexp.fits mos2S002-obj-imexp.fits' withexpimage=yes nruns=3 likemin=15 boxsize=5 withimagebuffersize=no imagebuffersize=640 ecf='3.2 1.2 1.2' imagesets='pnS003-obj-im.fits mos1S001-obj-im.fits mos2S002-obj-im.fits' boxlistset=eboxlist_m.fits withoffsets=no mergedlistset=mergedlist.fits hrdef='1 2 2 3 3 4' pimin='400 400 400' pimax='2300 2300 2300' obsmode=pointing  -w 1 -V 4

emldetect boxlistset=eboxlist_m.fits mllistset=emllist.fits mlmin=10 determineerrors=yes fitposition=yes psfmodel=ellbeta nmaxfit=1 nmulsou=1 ecut=15 scut=15 fitextent=yes extentmodel=beta dmlextmin=6 minextent=1.5 maxextent=20 withthreshold=yes threshold=15 threshcolumn=LIKE withtwostage=yes withimagebuffersize=no imagebuffersize=640 withxidband=no xidfixed=no rateonly=no simulate=no pimin='400 400 400' pimax='2300 2300 2300' hrpndef='1 2 2 3 3 4 4 5' xidpndef='2 3 4' hrm1def='1 2 2 3 3 4 4 5' xidm1def='2 3 4' hrm2def='1 2 2 3 3 4 4 5' xidm2def='2 3 4' ecf='3.2 1.2 1.2' xidecf=1 fitcounts=yes fitnegative=no mergedlistset=mergedlist.fits useevents=no usecalpsf=yes withhotpixelfilter=no withoffsets=no withrawrows=no tmpwrite=0 imagesets='pnS003-obj-im.fits mos1S001-obj-im.fits mos2S002-obj-im.fits' bkgimagesets='pnS003-obj-imbkg.fits mos1S001-obj-imbkg.fits mos2S002-obj-imbkg.fits' detmasksets='pnS003-obj-immask.fits mos1S001-obj-immask.fits mos2S002-obj-immask.fits' withdetmask=yes expimagesets='pnS003-obj-imexp.fits mos1S001-obj-imexp.fits mos2S002-obj-imexp.fits' withexpimage=yes sourceimagesets='pnS003-obj-imsmap.fits mos1S001-obj-imsmap.fits mos2S002-obj-imsmap.fits' withsourcemap=yes  -w 1 -V 4

emask detmaskset=mos1S001-mask-im.fits expimageset=mos1S001-obj-imexp.fits threshold1=0.1 threshold2=0.5 verbosity=1
   
emask detmaskset=mos2S002-mask-im.fits expimageset=mos2S002-obj-imexp.fits threshold1=0.1 threshold2=0.5 verbosity=1
   
emask detmaskset=pnS003-mask-im.fits expimageset=pnS003-obj-imexp.fits threshold1=0.1 threshold2=0.5 verbosity=1
   
region eventset=mos1S001-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(FLUX >= 1e-18)&&(DET_ML >= 15)&&(ID_INST == 2)&&(DIST_NN >= 0)&&(ID_BAND == 1)' bkgregionset=mos1S001-bkg_region-det.fits  bkgfraction=0.25 radiusstyle=contour nosrcellipse=no outunit=detxy verbosity=1
   
region eventset=mos1S001-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(DIST_NN >= 0)&&(DET_ML >= 15)&&(ID_INST == 2)&&(FLUX >= 1e-18)&&(ID_BAND == 1)' bkgregionset=mos1S001-bkg_region-sky.fits radiusstyle=contour nosrcellipse=no bkgfraction=0.25 outunit=xy verbosity=1
   
make_mask inimage=mos1S001-obj-im.fits inmask=mos1S001-mask-im.fits outmask=mos1S001-cheese.fits reglist=mos1S001-bkg_region-sky.fits
   
region eventset=mos2S002-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(FLUX >= 1e-18)&&(DET_ML >= 15)&&(ID_INST == 3)&&(DIST_NN >= 0)&&(ID_BAND == 1)' bkgregionset=mos2S002-bkg_region-det.fits  bkgfraction=0.25 radiusstyle=contour nosrcellipse=no outunit=detxy verbosity=1
   
region eventset=mos2S002-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(DIST_NN >= 0)&&(DET_ML >= 15)&&(ID_INST == 3)&&(FLUX >= 1e-18)&&(ID_BAND == 1)' bkgregionset=mos2S002-bkg_region-sky.fits radiusstyle=contour nosrcellipse=no bkgfraction=0.25 outunit=xy verbosity=1
   
make_mask inimage=mos2S002-obj-im.fits inmask=mos2S002-mask-im.fits outmask=mos2S002-cheese.fits reglist=mos2S002-bkg_region-sky.fits  
   
region eventset=pnS003-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(DIST_NN >= 0)&&(DET_ML >= 15)&&(ID_INST == 1)&&(FLUX >= 1e-18)&&(ID_BAND == 1)' bkgregionset=pnS003-bkg_region-det.fits  radiusstyle=contour nosrcellipse=no bkgfraction=0.25 outunit=detxy verbosity=1
   
region eventset=pnS003-clean.fits operationstyle=global srclisttab=emllist.fits:SRCLIST expression='(DIST_NN >= 0)&&(DET_ML >= 15)&&(ID_INST == 1)&&(FLUX >= 1e-18)&&(ID_BAND == 1)' bkgregionset=pnS003-bkg_region-sky.fits radiusstyle=contour nosrcellipse=no bkgfraction=0.25 outunit=xy verbosity=1
   
make_mask inimage=pnS003-obj-im.fits inmask=pnS003-mask-im.fits outmask=pnS003-cheese.fits reglist=pnS003-bkg_region-sky.fits
   
