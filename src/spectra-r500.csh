#!/bin/csh -f
set ccd=$1
set CALDB="/net/cluster491/software/newton/caldb/esas"
set coorfile="Centre.txt"

echo $ccd 
echo "calibration path" $CALDB


set r500_arcmin=`grep 'r500_mcxc_arcmin = ' r500_z.txt | gawk '{print $3}'`
set inner=`echo "$r500_arcmin * 0.15 * 60 * 20" | bc`  #inner radius = r500_in_arcmin*0.15, (60*20) is the conversion factor from arcmin to pixel, 1 pixel = 0.05 arcsec
set outer=`echo "$r500_arcmin * 60 * 20" | bc`

echo "inner radius in pixels: $inner"
echo "outer radius in pixels: $outer"

set mos1x=`grep 'mos1' Centre.txt | gawk '{print $3}'`
set mos1y=`grep 'mos1' Centre.txt | gawk '{print $4}'`

set mos2x=`grep 'mos2' Centre.txt | gawk '{print $3}'`
set mos2y=`grep 'mos2' Centre.txt | gawk '{print $4}'`

set pnx=`grep 'pn' Centre.txt | gawk '{print $3}'`
set pny=`grep 'pn' Centre.txt | gawk '{print $4}'`


if ($ccd =~ mos1*) then
   set ccdname="mos1"
   set ccdshort = ($ccd:as/mos/ /) 
endif


if ($ccd =~ mos2*) then
   set ccdname="mos2"
   set ccdshort = ($ccd:as/mos/ /) 
endif

if ($ccd =~ pn*) then
   set ccdname="pn"
   set ccdshort = ($ccd:as/pn/ /) 
endif


if ($ccdname == "mos1" || $ccdname == "mos2") then

   if (! -f InputFiles/goodccdlist_${ccdname}.dat) then
   echo InputFiles/goodccdlist_${ccdname}.dat is missing
   exit
   endif

   set goodccdlist=InputFiles/goodccdlist_${ccdname}.dat
   set ccd1=`head -1 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd2=`head -2 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd3=`head -3 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd4=`head -4 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd5=`head -5 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd6=`head -6 $goodccdlist | tail -1 | gawk '{print $2}'`
   set ccd7=`head -7 $goodccdlist | tail -1 | gawk '{print $2}'`

   if ($ccdname == "mos1") then
      set xc=$mos1x
      set yc=$mos1y
   endif

   if ($ccdname == "mos2") then
      set xc=$mos2x
      set yc=$mos2y
   endif

   rm -r ${ccd}-r500
   mkdir ${ccd}-r500
   cd ${ccd}-r500

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outer}))&&!((DETX,DETY) IN circle(${xc},${yc},${inner}))" > ${ccdname}-r500.reg

   echo "${ccdname}-r500 in progress"
   pwd
   date
   mos-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-r500.reg mask=1 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}-spectra-r500.log
   mos_back prefix=$ccdshort caldb=$CALDB diag=2 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}_back-r500.log

   rm ../${ccd}-obj-r500.pi
   rm ../${ccd}-back-r500.pi
   rm ../${ccd}-r500.rmf
   rm ../${ccd}-r500.arf
   rm ../${ccd}-sp-r500.fits
   rm ../${ccd}-r500-spec.qdp
   rm ../${ccd}-r500-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-r500.pi
   cp ${ccd}-back.pi ../${ccd}-back-r500.pi
   cp ${ccd}.rmf ../${ccd}-r500.rmf
   cp ${ccd}.arf ../${ccd}-r500.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-r500.fits
   cp ${ccd}-spec.qdp ../${ccd}-r500-spec.qdp
   cp ${ccd}-back.qdp ../${ccd}-r500-back.fits

   mv command.csh ../command/${ccdname}-spectra-r500.csh
   cd ..
  

endif


if ($ccdname == "pn") then

   set xc=$pnx
   set yc=$pny

   rm -r ${ccd}-r500
   mkdir ${ccd}-r500
   cd ${ccd}-r500

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-clean-oot.fits . 
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outer}))&&!((DETX,DETY) IN circle(${xc},${yc},${inner}))" > ${ccdname}-r500.reg

   echo "${ccdname}-r500 in progress"
   pwd
   date
   pn-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-r500.reg mask=1 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn-spectra-r500.log

   pn_back prefix=$ccdshort caldb=$CALDB diag=0 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn_back-r500.log

   rm ../${ccd}-obj-os-r500.pi
   rm ../${ccd}-obj-oot-r500.pi
   rm ../${ccd}-obj-r500.pi
   rm ../${ccd}-back-r500.pi
   rm ../${ccd}-r500.rmf
   rm ../${ccd}-r500.arf
   rm ../${ccd}-sp-r500.fits
   rm ../${ccd}-r500-spec.qdp
   rm ../${ccd}-r500-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-r500.pi
   cp ${ccd}-back.pi ../${ccd}-back-r500.pi
   cp ${ccd}.rmf ../${ccd}-r500.rmf
   cp ${ccd}.arf ../${ccd}-r500.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-r500.fits
   cp ${ccd}-spec.qdp ../${ccd}-r500-spec.qdp
   cp ${ccd}-back.qdp ../${ccd}-r500-back.fits
   cp ${ccd}-obj-os.pi ../${ccd}-obj-os-r500.pi
   cp ${ccd}-obj-oot.pi ../${ccd}-obj-oot-r500.pi

   mv command.csh ../command/${ccdname}-spectra-r500.csh
   cd ..


endif


