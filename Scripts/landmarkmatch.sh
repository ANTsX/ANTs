#!/bin/bash
if [ $# -lt 6 ]
then
echo " USAGE \n  sh $0  fixed.nii fixedhipp.nii  moving.nii movinghipp.nii  ITERATIONS LandmarkWeight "
echo " you should set LandmarkWeight yourself --- the intensity weight is 1.  so the landmark weight should be in that range (e.g. 0.5 => 1 ) "
echo "  the usage indicates .nii but you can use whatever image type you like (as long as itk can read it).
also, the fixedhipp.nii can contain either the full hippocampus or the landmarks for the hippocampus --- this would be your template.
the moving.nii can contain either one as well.  if the moving contains the landmarks and the fixed contains the full
hippocampus, then  this is known as "partial matching" as in Pluta, et al Hippocampus 2009.  "
exit
fi

TS=$5
LMWT=$6
INTWT=1

FIX=$1
FIXH=$2
MOV=$3
MOVH=$4



OUT=${MOV%.*}
NN=${#OUT}
MM=${#MOV}
let DD=MM-NN
if [ $DD -eq 3 ]
then
OUT=${MOV%.*.*}
fi

OUT2=${FIX%.*}
NN=${#OUT2}
MM=${#FIX}
let DD=MM-NN
if [ $DD -eq 3 ]
then
OUT2=${FIX%.*.*}
fi

dothis=0
if [ $dothis -eq 1 ]
then
# here, we select a subset of the labeled regions
FIXHMOD=${OUT2}locmod.nii.gz
MOVHMOD=${OUT}locmod.nii.gz
for i in  $MOVH
do
PRE=$OUT
TOUT=${PRE}locmod.nii.gz
echo " DOING $TOUT "
ThresholdImage 3 $i ${PRE}tempb.nii.gz 30.5 32.5
MultiplyImages 3 $i  ${PRE}tempb.nii.gz $TOUT
ThresholdImage 3 $i ${PRE}tempb.nii.gz 20.5 22.5
MultiplyImages 3 $i ${PRE}tempb.nii.gz ${PRE}tempb.nii
ImageMath 3 $TOUT + $TOUT ${PRE}tempb.nii.gz
ThresholdImage 3 $i ${PRE}tempb.nii.gz 8.5 10.5
MultiplyImages 3 $i ${PRE}tempb.nii.gz ${PRE}tempb.nii.gz
ImageMath 3 $TOUT + $TOUT ${PRE}tempb.nii.gz
done
else
FIXHMOD=$FIXH
MOVHMOD=$MOVH
fi

KNN=10
PERCENTtoUSE=0.5
PARZSIG=10
LM=PSE[ ${FIX},${MOV},${FIXHMOD},${MOVHMOD},${LMWT},${PERCENTtoUSE},${PARZSIG},0,${KNN}]
echo $LM

if [ $# -gt 6 ]
then
OUT=$7
fi

STEPL=0.25
INTENSITY=PR[ $FIX,${MOV},${INTWT},4]
if [ $LMWT -le 0 ]
then
exe="ANTS 3  -o $OUT  -i $ITS -t SyN[ ${STEPL}]  -r Gauss[ 3,0 ]   -m $INTENSITY   "
else
exe="ANTS 3  -o $OUT  -i $ITS -t SyN[ ${STEPL}]  -r Gauss[ 3,0 ]   -m   $LM  -m $INTENSITY    "
fi
echo " $exe "

 $exe

 WarpImageMultiTransform 3 $MOV ${OUT}deformed.nii.gz  ${OUT}Warp.nii.gz ${OUT}Affine.txt  -R $FIX
 WarpImageMultiTransform 3 $MOVH  ${OUT}label.nii.gz   ${OUT}Warp.nii.gz ${OUT}Affine.txt  -R $FIX --use-NN
 WarpImageMultiTransform 3 $FIXH ${OUT}labelinv.nii.gz  -i  ${OUT}Affine.txt  ${OUT}InverseWarp.nii.gz   -R $MOV --use-NN




