#!/bin/csh -f
set ccd=$1
set CALDB=$2
set coorfile="Centre.txt"

echo $ccd 
echo "calibration path" $CALDB


set Region="spectral_region.dat"

set SRCCOUNT=`grep -i source $Region | wc -l`
set CXBCOUNT=`grep -i cxb $Region | wc -l`

@ TOTAL = $SRCCOUNT + $CXBCOUNT

set mos1x=`grep 'mos1' $coorfile | gawk '{print $3}'`
set mos1y=`grep 'mos1' $coorfile | gawk '{print $4}'`

set mos2x=`grep 'mos2' $coorfile | gawk '{print $3}'`
set mos2y=`grep 'mos2' $coorfile | gawk '{print $4}'`

set pnx=`grep 'pn' $coorfile | gawk '{print $3}'`
set pny=`grep 'pn' $coorfile | gawk '{print $4}'`



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

   @ i = 1
   while({$i} <= $SRCCOUNT - 1)
      set inner=`head -$i $Region  | tail -1 | gawk '{print int($1)}'`
      @ j = $i + 1
      set outer=`head -$j $Region  | tail -1 | gawk '{print int($1)}'`
      @ i ++
      @ innerradius = ${inner} * 20
      @ outerradius = ${outer} * 20

   rm -r ${ccd}-${inner}-${outer}
   mkdir ${ccd}-${inner}-${outer}
   cd ${ccd}-${inner}-${outer}

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outerradius}))&&!((DETX,DETY) IN circle(${xc},${yc},${innerradius}))" > ${ccdname}-${inner}-${outer}.reg

   echo "${ccdname}-${inner}-${outer} in progress"
   pwd
   date
   mos-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-${inner}-${outer}.reg mask=1 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}-spectra-${inner}-${outer}.log
   mos_back prefix=$ccdshort caldb=$CALDB diag=2 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}_back-${inner}-${outer}.log

   rm ../${ccd}-obj-${inner}-${outer}.pi
   rm ../${ccd}-back-${inner}-${outer}.pi
   rm ../${ccd}-${inner}-${outer}.rmf
   rm ../${ccd}-${inner}-${outer}.arf
   rm ../${ccd}-sp-${inner}-${outer}.fits
   rm ../${ccd}-${inner}-${outer}-spec.qdp
   rm ../${ccd}-${inner}-${outer}-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-${inner}-${outer}.pi
   cp ${ccd}-back.pi ../${ccd}-back-${inner}-${outer}.pi
   cp ${ccd}.rmf ../${ccd}-${inner}-${outer}.rmf
   cp ${ccd}.arf ../${ccd}-${inner}-${outer}.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-${inner}-${outer}.fits
   cp ${ccd}-spec.qdp ../${ccd}-${inner}-${outer}-spec.qdp
   cp ${ccd}-back.qdp ../${ccd}-${inner}-${outer}-back.fits

   mv command.csh ../command/${ccdname}-spectra-${inner}-${outer}.csh
   cd ..
   end

   @ i = $SRCCOUNT + 1
   while({$i} <= $TOTAL - 1)
      set inner=`head -$i $Region  | tail -1 | gawk '{print int($1)}'`
      @ j = $i + 1
      set outer=`head -$j $Region  | tail -1 | gawk '{print int($1)}'`
      @ i ++
      @ innerradius = ${inner} * 20
      @ outerradius = ${outer} * 20

   rm -r ${ccd}-${inner}-${outer}
   mkdir ${ccd}-${inner}-${outer}
   cd ${ccd}-${inner}-${outer}

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outerradius}))&&!((DETX,DETY) IN circle(${xc},${yc},${innerradius}))" > ${ccdname}-${inner}-${outer}.reg

   echo "${ccdname}-${inner}-${outer} in progress"
   pwd
   date
   mos-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-${inner}-${outer}.reg mask=1 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}-spectra-${inner}-${outer}.log
   mos_back prefix=$ccdshort caldb=$CALDB diag=2 elow=0 ehigh=0 ccd1=${ccd1} ccd2=${ccd2} ccd3=${ccd3} ccd4=${ccd4} ccd5=${ccd5} ccd6=${ccd6} ccd7=${ccd7} >& ../log/${ccdname}_back-${inner}-${outer}.log

   rm ../${ccd}-obj-${inner}-${outer}.pi
   rm ../${ccd}-back-${inner}-${outer}.pi
   rm ../${ccd}-${inner}-${outer}.rmf
   rm ../${ccd}-${inner}-${outer}.arf
   rm ../${ccd}-sp-${inner}-${outer}.fits
   rm ../${ccd}-${inner}-${outer}-spec.qdp
   rm ../${ccd}-${inner}-${outer}-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-${inner}-${outer}.pi
   cp ${ccd}-back.pi ../${ccd}-back-${inner}-${outer}.pi
   cp ${ccd}.rmf ../${ccd}-${inner}-${outer}.rmf
   cp ${ccd}.arf ../${ccd}-${inner}-${outer}.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-${inner}-${outer}.fits
   cp ${ccd}-spec.qdp ../${ccd}-${inner}-${outer}-spec.qdp
   cp ${ccd}-back.qdp ../${ccd}-${inner}-${outer}-back.fits

   mv command.csh ../command/${ccdname}-spectra-${inner}-${outer}.csh
   cd ..
   end

endif


if ($ccdname == "pn") then

   set xc=$pnx
   set yc=$pny

   @ i = 1
   while({$i} <= $SRCCOUNT - 1)
      set inner=`head -$i $Region  | tail -1 | awk '{print int($1)}'`
      @ j = $i + 1
      set outer=`head -$j $Region  | tail -1 | awk '{print int($1)}'`
      @ i ++
      @ innerradius = ${inner} * 20
      @ outerradius = ${outer} * 20


   rm -r ${ccd}-${inner}-${outer}
   mkdir ${ccd}-${inner}-${outer}
   cd ${ccd}-${inner}-${outer}

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-clean-oot.fits . 
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outerradius}))&&!((DETX,DETY) IN circle(${xc},${yc},${innerradius}))" > ${ccdname}-${inner}-${outer}.reg

   echo "${ccdname}-${inner}-${outer} in progress"
   date
   pn-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-${inner}-${outer}.reg mask=1 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn-spectra-${inner}-${outer}.log

   pn_back prefix=$ccdshort caldb=$CALDB diag=0 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn_back-${inner}-${outer}.log

   rm ../${ccd}-obj-os-${inner}-${outer}.pi
   rm ../${ccd}-obj-oot-${inner}-${outer}.pi
   rm ../${ccd}-obj-${inner}-${outer}.pi
   rm ../${ccd}-back-${inner}-${outer}.pi
   rm ../${ccd}-${inner}-${outer}.rmf
   rm ../${ccd}-${inner}-${outer}.arf
   rm ../${ccd}-sp-${inner}-${outer}.fits
   rm ../${ccd}-${inner}-${outer}-spec.qdp
   rm ../${ccd}-${inner}-${outer}-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-${inner}-${outer}.pi
   cp ${ccd}-back.pi ../${ccd}-back-${inner}-${outer}.pi
   cp ${ccd}.rmf ../${ccd}-${inner}-${outer}.rmf
   cp ${ccd}.arf ../${ccd}-${inner}-${outer}.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-${inner}-${outer}.fits
   cp ${ccd}-spec.qdp ../${ccd}-${inner}-${outer}-spec.qdp
   cp ${ccd}-back.qdp ../${ccd}-${inner}-${outer}-back.fits
   cp ${ccd}-obj-os.pi ../${ccd}-obj-os-${inner}-${outer}.pi
   cp ${ccd}-obj-oot.pi ../${ccd}-obj-oot-${inner}-${outer}.pi

   mv command.csh ../command/${ccdname}-spectra-${inner}-${outer}.csh
   cd ..
   end

   @ i = $SRCCOUNT + 1
   while({$i} <= $TOTAL - 1)
      set inner=`head -$i $Region  | tail -1 | awk '{print int($1)}'`
      @ j = $i + 1
      set outer=`head -$j $Region  | tail -1 | awk '{print int($1)}'`
      @ i ++
      @ innerradius = ${inner} * 20
      @ outerradius = ${outer} * 20


   rm -r ${ccd}-${inner}-${outer}
   mkdir ${ccd}-${inner}-${outer}
   cd ${ccd}-${inner}-${outer}

   setenv SAS_CCF ../ccf.cif 
   setenv SAS_ODF `ls -1 ../*SUM.SAS`

   ln -s ../${ccd}-clean.fits .
   ln -s ../${ccd}-clean-oot.fits . 
   ln -s ../${ccd}-bkg_region-sky.fits .
   ln -s ../${ccd}-bkg_region-det.fits .

   echo "&&((DETX,DETY) IN circle(${xc},${yc},${outerradius}))&&!((DETX,DETY) IN circle(${xc},${yc},${innerradius}))" > ${ccdname}-${inner}-${outer}.reg

   echo "${ccdname}-${inner}-${outer} in progress"
   pwd
   date
   pn-spectra prefix=$ccdshort caldb=$CALDB region=${ccdname}-${inner}-${outer}.reg mask=1 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn-spectra-${inner}-${outer}.log

   pn_back prefix=$ccdshort caldb=$CALDB diag=0 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn_back-${inner}-${outer}.log

   rm ../${ccd}-obj-os-${inner}-${outer}.pi
   rm ../${ccd}-obj-oot-${inner}-${outer}.pi
   rm ../${ccd}-obj-${inner}-${outer}.pi
   rm ../${ccd}-back-${inner}-${outer}.pi
   rm ../${ccd}-${inner}-${outer}.rmf
   rm ../${ccd}-${inner}-${outer}.arf
   rm ../${ccd}-sp-${inner}-${outer}.fits
   rm ../${ccd}-${inner}-${outer}-spec.qdp
   rm ../${ccd}-${inner}-${outer}-back.fits

   cp ${ccd}-obj.pi ../${ccd}-obj-${inner}-${outer}.pi
   cp ${ccd}-back.pi ../${ccd}-back-${inner}-${outer}.pi
   cp ${ccd}.rmf ../${ccd}-${inner}-${outer}.rmf
   cp ${ccd}.arf ../${ccd}-${inner}-${outer}.arf
   cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-${inner}-${outer}.fits
   cp ${ccd}-spec.qdp ../${ccd}-${inner}-${outer}-spec.qdp
   cp ${ccd}-obj-os.pi ../${ccd}-obj-os-${inner}-${outer}.pi
   cp ${ccd}-obj-oot.pi ../${ccd}-obj-oot-${inner}-${outer}.pi

   mv command.csh ../command/${ccdname}-spectra-${inner}-${outer}.csh
   cd ..
   end

endif


