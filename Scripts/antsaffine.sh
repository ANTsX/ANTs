#!/bin/sh
if [ ${#ANTSPATH} -le 3 ] ; then
  echo we guess at your ants path
  export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS
fi
if [ ! -s ${ANTSPATH}/ANTS ] ; then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

NUMPARAMS=$#

if [ $NUMPARAMS -lt 3  ]
then
echo " USAGE ::  "
echo "  sh   ants.sh  ImageDimension  fixed.ext  moving.ext  OPTIONAL-OUTPREFIX   PURELY-RIGID  "
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

RIGID=" --rigid-affine false  "
if [ $NUMPARAMS -gt 4 ]
then
  RIGID=" --rigid-affine true  --affine-gradient-descent-option  0.5x0.95x1.e-4x1.e-4  "
fi
echo " Will this mapping be purely rigid?  $RIGID "

echo  " ANTSPATH  is $ANTSPATH     "
 #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

exe=" ${ANTSPATH}ANTS $DIM -m  MI[${FIXED},${MOVING},1,32] -o ${OUTPUTNAME}   -i 0   --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000  $RIGID   "

 echo " $exe "

  $exe

    ${ANTSPATH}WarpImageMultiTransform $DIM  $MOVING  ${OUTPUTNAME}deformed.nii.gz ${OUTPUTNAME}Affine.txt  -R ${FIXED}


