#!/bin/sh

NUMPARAMS=$#

if [ $NUMPARAMS -lt 3  ]
then
echo " USAGE ::  "
echo "  sh   ants.sh  ImageDimension  fixed.ext  moving.ext  OPTIONAL-OUTPREFIX     "
echo " be sure to set ANTSPATH environment variable "
echo " affine only registration "
exit
fi

#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
then
echo " Please set ANTSPATH=LocationOfYourAntsBinaries "
echo " Either set this in your profile or directly, here in the script. "
echo " For example : "
echo $ANTSPATH
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

if [ $NUMPARAMS -gt 3  ]
then
  OUTPUTNAME=${4}
fi


echo  " ANTSPATH  is $ANTSPATH     "

# finally check the image headers
compareheaders=`${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}repaired.nii  CompareHeadersAndImages $FIXED $MOVING  | grep FailureState | cut -d ' ' -f 4  `
if [ $compareheaders -ne 0 ]
then
echo " You may have a problem with your header definition "
echo " The repaired image is in : ${OUTPUTNAME}repaired.nii "
echo " Call ImageMath's   CompareHeadersAndImages  on your Fixed and Moving image "
exit
fi

 #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

exe=" ${ANTSPATH}ANTS $DIM -m  MI[${FIXED},${MOVING},1,32] -o ${OUTPUTNAME}   -i 0   --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000 "

 echo " $exe "

  $exe

    ${ANTSPATH}WarpImageMultiTransform $DIM  ${OUTPUTNAME}repaired.nii   ${OUTPUTNAME}deformed.nii ${OUTPUTNAME}Affine.txt  -R ${FIXED}


exit
