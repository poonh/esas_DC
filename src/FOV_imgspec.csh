#!/bin/csh -f
set ccd=$1
set elow=$2
set ehigh=$3
set CALDB=$4

echo $ccd "in" $elow "-" $ehigh "eV"
echo "calibration path" $CALDB


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


rm -r full_spectrum_${ccdname}
mkdir full_spectrum_${ccdname}
cd full_spectrum_${ccdname}

ln -s ../${ccd}-clean.fits .
ln -s ../${ccd}-bkg_region-sky.fits .
ln -s ../${ccd}-bkg_region-det.fits .
ln -s ../${ccd}-cheese.fits .

setenv SAS_CCF ../ccf.cif 
setenv SAS_ODF `ls -1 ../*SUM.SAS`


if ($ccdname == "mos1" || $ccdname == "mos2") then

if (! -f ../InputFiles/goodccdlist_${ccdname}.dat) then
echo ../InputFiles/goodccdlist_${ccdname}.dat is missing
exit
endif

set goodccdlist=../InputFiles/goodccdlist_${ccdname}.dat
set ccd1=`head -1 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd2=`head -2 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd3=`head -3 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd4=`head -4 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd5=`head -5 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd6=`head -6 $goodccdlist | tail -1 | gawk '{print $2}'`
set ccd7=`head -7 $goodccdlist | tail -1 | gawk '{print $2}'`


echo ${ccdname}
echo "mos-spectra in progress"
mos-spectra prefix=$ccdshort caldb=$CALDB region=nonfile mask=1 elow=$elow ehigh=$ehigh ccd1=$ccd1 ccd2=$ccd2 ccd3=$ccd3 ccd4=$ccd4 ccd5=$ccd5 ccd6=$ccd6 ccd7=$ccd7 >& ../log/mos-spectra_${ccdname}.log

echo "mos_back in progress"
mos_back prefix=$ccdshort caldb=$CALDB diag=0 elow=$elow ehigh=$ehigh ccd1=$ccd1 ccd2=$ccd2 ccd3=$ccd3 ccd4=$ccd4 ccd5=$ccd5 ccd6=$ccd6 ccd7=$ccd7 >& ../log/${ccdname}_back.log


endif


if ($ccdname == "pn") then

echo "pwd"
pwd

ln -s ../${ccd}-clean-oot.fits .


echo "pn-spectra in progress"
pn-spectra prefix=$ccdshort caldb=$CALDB region=nonfile mask=1 elow=$elow ehigh=$ehigh quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn-spectra.log

echo "pn_back in progress"
pn_back prefix=$ccdshort caldb=$CALDB diag=0 elow=$elow ehigh=$ehigh quad1=1 quad2=1 quad3=1 quad4=1 >& ../log/pn_back.log


rm ../${ccd}-obj-oot-FOV.pi
cp ${ccd}-obj-oot.pi ../${ccd}-obj-oot-FOV.pi 


endif



rm ../${ccd}-obj-FOV.pi
rm ../${ccd}-FOV.rmf
rm ../${ccd}-FOV.arf
rm ../${ccd}-back-FOV.pi
rm ../${ccd}-sp-FOV.fits


cp ${ccd}-obj.pi ../${ccd}-obj-FOV.pi
cp ${ccd}.rmf ../${ccd}-FOV.rmf
cp ${ccd}.arf ../${ccd}-FOV.arf
cp ${ccd}-back.pi ../$ccd}-back-FOV.pi
cp ${ccd}-obj-im-sp-det.fits ../${ccd}-sp-FOV.fits


echo "Producing particle bkg in sky coordinates for $ccdname"
rot-im-det-sky prefix=$ccdshort mask=1 elow=$elow ehigh=$ehigh mode=1 >& ../log/rot-im-det-sky_${ccdname}.log


echo "Producing combined images for $ccdname"
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$elow ehighlist=$ehigh mask=1 prefixlist="$ccdshort" >& ../log/comb_${ccdname}.log
bin_image thresholdmasking=0.02 detector=1 binning=1 prefix="$ccdshort" elow=$elow ehigh=$ehigh withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& ../log/bin_image-${ccdname}.log

 
mv rate-$elow-$ehigh.fits rate-$elow-$ehigh-${ccdname}.fits
mv sigma-$elow-$ehigh.fits sigma-$elow-$ehigh-${ccdname}.fits
mv comb-obj-im-$elow-$ehigh.fits comb-obj-im-$elow-$ehigh-${ccdname}.fits
mv comb-back-im-sky-$elow-$ehigh.fits comb-back-im-sky-$elow-$ehigh-${ccdname}.fits
mv comb-exp-im-$elow-$ehigh.fits comb-exp-im-$elow-$ehigh-${ccdname}.fits


cp rate-$elow-$ehigh-${ccdname}.fits sigma-$elow-$ehigh-${ccdname}.fits comb-obj-im-$elow-$ehigh-${ccdname}.fits comb-back-im-sky-$elow-$ehigh-${ccdname}.fits comb-exp-im-$elow-$ehigh-${ccdname}.fits ../.

mv command.csh FOV_imgspec_${ccdname}.csh
cd ../.



