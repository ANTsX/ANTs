#!/bin/bash

NUMPARAMS=$#

MAXITERATIONS=30x90x20
if [ ${#ANTSPATH} -le 3 ] ; then
  echo we guess at your ants path
  export ANTSPATH=${ANTSPATH:="/ipldev/scratch/aghayoor/ANTS/release/bin"} # EDIT THIS
fi
if [ ! -s ${ANTSPATH}/ANTS ] ; then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

if [ $NUMPARAMS -lt 3  ]
then
echo " USAGE ::  "
echo "  sh   ants.sh  ImageDimension  fixed.ext  moving.ext  OPTIONAL-OUTPREFIX   OPTIONAL-max-iterations  OPTIONAL-Labels-In-Fixed-Image-Space-To-Deform-To-Moving-Image     Option-DoANTSQC "
echo " be sure to set ANTSPATH environment variable "
echo " Max-Iterations in form :    JxKxL where "
echo "  J = max iterations at coarsest resolution (here, reduce by power of 2^2) "
echo " K = middle resolution iterations ( here, reduce by power of 2 ) "
echo " L = fine resolution iterations ( here, full resolution ) -- this level takes much more time per iteration "
echo " an extra  Ix before JxKxL would add another level "
echo " Default ierations is  $MAXITERATIONS  -- you can often get away with fewer for many apps "
echo " Other parameters are updates of the defaults used in the A. Klein evaluation paper in Neuroimage, 2009 "
exit
fi

#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
then
echo " Please set ANTSPATH=LocationOfYourAntsBinaries "
echo " Either set this in your profile or directly, here in the script. "
echo " For example : "
echo " ANTSPATH=/home/yourname/bin/ants/ "
exit
else
echo " ANTSPATH is $ANTSPATH "
fi

#initialization, here, is unbiased
DIM=$1

if  [ ${#DIM} -gt 1 ]
then
echo " Problem with specified ImageDimension => User Specified Value = $DIM "
exit
fi

FIXED=$2

if  [ ${#FIXED} -lt 1 -o  ! -f $FIXED ]
then
echo " Problem with specified Fixed Image => User Specified Value = $FIXED "
exit
fi

MOVING=$3

if  [ ${#MOVING} -lt 1 -o  ! -f $MOVING ]
then
echo " Problem with specified Moving Image => User Specified Value = $MOVING "
exit
fi

OUTPUTNAME=` echo $MOVING | cut -d '.' -f 1 `


if [ $NUMPARAMS -gt 4  ]
then
 MAXITERATIONS=$5
fi

LABELIMAGE=0
if [ $NUMPARAMS -gt 5  ]
then
 LABELIMAGE=$6
fi

DoANTSQC=0
if [ $NUMPARAMS -gt 6  ]
then
 DoANTSQC=$7
fi

if [ $NUMPARAMS -gt 3  ]
then
  OUTPUTNAME=${4}
fi

# Mapping Parameters
  TRANSFORMATION=SyN[0.25]
  ITERATLEVEL=(`echo $MAXITERATIONS | tr 'x' ' '`)
  NUMLEVELS=${#ITERATLEVEL[@]}
  echo $NUMLEVELS
  REGULARIZATION=Gauss[3,0]
  METRIC=CC[
    METRICPARAMS=1,4]
# echo " $METRICPARAMS  &  $METRIC "
# exit


echo  " ANTSPATH  is $ANTSPATH     "
echo " Mapping Parameters  :: "
echo  " Transformation is:  $TRANSFORMATION "
echo " MaxIterations :   $MAXITERATIONS "
echo " Number Of MultiResolution Levels   $NUMLEVELS "
echo " Regularization :  $REGULARIZATION "
echo " Metric :  ${METRIC} "
echo " OutputName :  $OUTPUTNAME "
echo " "
echo " if the files and parameters are all ok then uncomment the exit call below this line  "
echo " "
#exit


if [[ ! -s ${OUTPUTNAME}repaired.nii.gz ]] ; then
${ANTSPATH}N4BiasFieldCorrection -d $DIM -i $MOVING -o ${OUTPUTNAME}repaired.nii.gz -s 2 -c [50x50x50x50,0.000001] -b [200]
# ${ANTSPATH}N3BiasFieldCorrection $DIM $MOVING   ${OUTPUTNAME}repaired.nii.gz  4
fi
exe=" ${ANTSPATH}ANTS $DIM -m  ${METRIC}${FIXED},${OUTPUTNAME}repaired.nii.gz,${METRICPARAMS}  -t $TRANSFORMATION  -r $REGULARIZATION -o ${OUTPUTNAME}   -i $MAXITERATIONS   --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000  "

 echo " $exe "

  $exe

  #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

    ${ANTSPATH}WarpImageMultiTransform $DIM  ${OUTPUTNAME}repaired.nii.gz    ${OUTPUTNAME}deformed.nii.gz  ${OUTPUTNAME}Warp.nii.gz  ${OUTPUTNAME}Affine.txt  -R ${FIXED}

if  [ ${#LABELIMAGE} -gt 3 ]
then
      ${ANTSPATH}WarpImageMultiTransform $DIM  $LABELIMAGE   ${OUTPUTNAME}labeled.nii.gz  -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz   -R ${MOVING}   --use-NN
fi

exit

if [ $DoANTSQC -eq 1 ] ;  then
#  measure image similarity
for SIM in 0 1 2 ; do
${ANTSPATH}MeasureImageSimilarity $DIM $SIM $FIXED ${OUTPUTNAME}deformed.nii.gz
done
#  measure dice overlap and mds
${ANTSPATH}ThresholdImage $DIM $FIXED ${OUTPUTNAME}fixthresh.nii.gz Otsu 4
${ANTSPATH}ThresholdImage $DIM $MOVING ${OUTPUTNAME}movthresh.nii.gz Otsu 4
 ${ANTSPATH}WarpImageMultiTransform $DIM   ${OUTPUTNAME}movthresh.nii.gz  ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}Warp.nii.gz  ${OUTPUTNAME}Affine.txt  -R ${FIXED}   --use-NN
${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}dicestats.txt DiceAndMinDistSum  ${OUTPUTNAME}fixthresh.nii.gz   ${OUTPUTNAME}movthresh.nii.gz   ${OUTPUTNAME}mindistsum.nii.gz
#  labelstats for jacobian wrt segmenation
 # below, to compose
# ${ANTSPATH}ComposeMultiTransform $DIM   ${OUTPUTNAME}CompWarp.nii.gz   -R $FIXED ${OUTPUTNAME}Warp.nii.gz  ${OUTPUTNAME}Affine.txt
# ${ANTSPATH}CreateJacobianDeterminantImage $DIM ${OUTPUTNAME}CompWarp.nii.gz  ${OUTPUTNAME}jacobian.nii.gz   0
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}movlabstat.txt LabelStats ${OUTPUTNAME}movthresh.nii.gz ${OUTPUTNAME}movthresh.nii.gz
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}jaclabstat.txt LabelStats ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}jacobian.nii.gz
# we compare the output of these last two lines:
#  the Volume of the movlabstat computation vs. the mass of the jaclabstat
fi

