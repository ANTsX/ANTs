#!/bin/bash
if [ $# -lt 7  ]
then
echo " USAGE \n  sh command.sh  fixed.nii fixedhipp.nii  moving.nii movinghipp.nii outputname  iterations DIM "
echo " the template = fixed.nii   ,   the individual = moving.nii "
echo " iterations should be of the form  100x100x10 "
exit
fi
if [ ${#ANTSPATH} -le 3 ] ; then
  echo we guess at your ants path
  export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS
fi
if [ ! -s ${ANTSPATH}/ANTS ] ; then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

FIX=$1
FIXH=$2
MOV=$3
MOVH=$4
OUT=$5
ITS=$6
DIM=$7

#ANTSPATH='/Users/stnava/Code/bin/ants/'

if  [ ${#FIX} -lt 1 -o  ! -f $FIX ]
then
echo " Problem with specified Fixed Image => User Specified Value = $FIX "
exit
fi
if  [ ${#MOV} -lt 1 -o  ! -f $MOV ]
then
echo " Problem with specified Moving Image => User Specified Value = $MOV "
exit
fi
if  [ ${#FIXH} -lt 1 -o  ! -f $FIXH ]
then
echo " Problem with specified Fixed Label Image => User Specified Value = $FIXH "
exit
fi
if  [ ${#MOVH} -lt 1 -o  ! -f $MOVH ]
then
echo " Problem with specified Moving Label Image => User Specified Value = $MOVH "
exit
fi

#  == Important Parameters Begin ==

LMWT=$9  # weight on landmarks

INTWT=$8   # weight on intensity --  twice the landmarks

# PSE/point-set-expectation/PointSetExpectation[fixedImage,movingImage,fixedPoints,movingPoints,weight,pointSetPercentage,pointSetSigma,boundaryPointsOnly,kNeighborhood, PartialMatchingIterations=100000]
# the partial matching option assumes the complete labeling is in the first set of label parameters ...
# more iterations leads to more symmetry in the matching  - 0 iterations means full asymmetry
PCT=0.1 # percent of labeled voxels to use
PARZ=100 # PARZEN sigma
LM=PSE[${FIX},${MOV},$FIXH,$MOVH,${LMWT},${PCT},${PARZ},0,25,100]

INTENSITY=CC[$FIX,${MOV},${INTWT},4]

#  == Important Parameters end? ==

 ${ANTSPATH}ANTS $DIM -o $OUT  -i $ITS -t SyN[0.25]  -r Gauss[3,0] -m $INTENSITY   -m   $LM

 ${ANTSPATH}WarpImageMultiTransform $DIM $MOV ${OUT}toTemplate.nii.gz ${OUT}Warp.nii.gz ${OUT}Affine.txt  -R $FIX

 ${ANTSPATH}WarpImageMultiTransform $DIM  $FIX ${OUT}toMov.nii.gz -i ${OUT}Affine.txt  ${OUT}InverseWarp.nii.gz  -R $MOV

 ${ANTSPATH}WarpImageMultiTransform $DIM $FIXH  ${OUT}hipp.nii.gz -i ${OUT}Affine.txt  ${OUT}InverseWarp.nii.gz  -R $MOV --UseNN

