#!/bin/bash

if ! command -v ANTS &> /dev/null
then
  echo "Cannot find the ANTS program. Please \(re\)define \$PATH in your environment."
  exit 1
fi

NUMPARAMS=$#

if [ $NUMPARAMS -lt 3 ]
then
echo " USAGE ::  "
echo "  sh   antsaffine.sh  ImageDimension  fixed.ext  moving.ext  OPTIONAL-OUTPREFIX   PURELY-RIGID  "
echo " affine only registration "
exit
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

if [ $NUMPARAMS -gt 3 ]
then
  OUTPUTNAME=${4}
fi

RIGID=" --rigid-affine false  "
if [ $NUMPARAMS -gt 4 ]
then
  RIGID=" --rigid-affine true  --affine-gradient-descent-option  0.5x0.95x1.e-4x1.e-4  "
fi
echo " Will this mapping be purely rigid?  $RIGID "

 #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

exe=" ANTS $DIM -m  MI[ ${FIXED},${MOVING},1,32 ] -o ${OUTPUTNAME}   -i 0   --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000  $RIGID   "

 echo " $exe "

  $exe

    WarpImageMultiTransform $DIM  $MOVING  ${OUTPUTNAME}deformed.nii.gz ${OUTPUTNAME}Affine.txt  -R ${FIXED}
