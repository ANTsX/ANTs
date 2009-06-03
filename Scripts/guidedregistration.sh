#!/bin/sh
if [ $# -lt 5  ]
then
echo " USAGE \n  sh command.sh  fixed.nii fixedhipp.nii  moving.nii movinghipp.nii outputname "
echo " the template = fixed.nii   ,   the individual = moving.nii "
exit
fi

FIX=$1
FIXH=$2
MOV=$3
MOVH=$4
OUT=$5

ANTSPATH='/Users/stnava/Code/bin/ants/'

#  == Important Parameters Begin ==

LMWT=0.5  # weight on landmarks

INTWT=1    # weight on intensity --  twice the landmarks

LM=PSE[${FIX},${MOV},$FIXH,$MOVH,${LMWT},0.1,100,0,25,10]

INTENSITY=PR[$FIX,${MOV},${INTWT},2]

#  == Important Parameters end? ==

 ${ANTSPATH}ANTS 3  -o $OUT  -i 55x50x30 -t SyN[1.5]  -r Gauss[3,0] -m $INTENSITY   -m   $LM

 ${ANTSPATH}WarpImageMultiTransform 3 $MOV ${OUT}toTemplate.nii ${OUT}Warp.nii ${OUT}Affine.txt  -R $FIX


 ${ANTSPATH}WarpImageMultiTransform 3  $FIX ${OUT}toMov.nii -i ${OUT}Affine.txt  ${OUT}InverseWarp.nii  -R $MOV

 ${ANTSPATH}WarpImageMultiTransform 3 $FIXH  ${OUT}hipp.nii -i ${OUT}Affine.txt  ${OUT}InverseWarp.nii  -R $MOV --UseNN
