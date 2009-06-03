#!/bin/sh
if [ $# -lt 6  ]
then
echo " USAGE \n  sh $0  fixed.nii fixedhipp.nii  moving.nii movinghipp.nii  ITERATIONS LandmarkWeight "
exit
fi

ANTSPATH="/mnt/aibs1/avants/bin/ants/"
ITS=$5
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
FIXHMOD=${OUT2}locmod.nii
MOVHMOD=${OUT}locmod.nii
for i in  $MOVH
do
PRE=$OUT
TOUT=${PRE}locmod.nii
echo " DOING $TOUT "
${ANTSPATH}ThresholdImage 3 $i ${PRE}tempb.nii 30.5 32.5
${ANTSPATH}MultiplyImages 3 $i  ${PRE}tempb.nii $TOUT
${ANTSPATH}ThresholdImage 3 $i ${PRE}tempb.nii 20.5 22.5
${ANTSPATH}MultiplyImages 3 $i ${PRE}tempb.nii ${PRE}tempb.nii
${ANTSPATH}ImageMath 3 $TOUT + $TOUT ${PRE}tempb.nii
${ANTSPATH}ThresholdImage 3 $i ${PRE}tempb.nii 8.5 10.5
${ANTSPATH}MultiplyImages 3 $i ${PRE}tempb.nii ${PRE}tempb.nii
${ANTSPATH}ImageMath 3 $TOUT + $TOUT ${PRE}tempb.nii
done
else
FIXHMOD=$FIXH
MOVHMOD=$MOVH
fi

KNN=10
PERCENTtoUSE=0.5
PARZSIG=10
LM=PSE[${FIX},${MOV},${FIXHMOD},${MOVHMOD},${LMWT},${PERCENTtoUSE},${PARZSIG},0,${KNN}]
echo $LM

if [ $# -gt 6  ]
then
OUT=$7
fi

STEPL=0.25
INTENSITY=PR[$FIX,${MOV},${INTWT},2]
if [ $LMWT -le 0  ]
then
exe="${ANTSPATH}ANTS 3  -o $OUT  -i $ITS -t SyN[${STEPL}]  -r Gauss[3,0]   -m $INTENSITY   "
else
exe="${ANTSPATH}ANTS 3  -o $OUT  -i $ITS -t SyN[${STEPL}]  -r Gauss[3,0]   -m   $LM  -m $INTENSITY    "
fi
echo " $exe "

 $exe

 ${ANTSPATH}WarpImageMultiTransform 3 $MOV ${OUT}deformed.nii  ${OUT}Warp.nii ${OUT}Affine.txt  -R $FIX
 ${ANTSPATH}WarpImageMultiTransform 3 $MOVH  ${OUT}label.nii   ${OUT}Warp.nii ${OUT}Affine.txt  -R $FIX --use-NN
 ${ANTSPATH}WarpImageMultiTransform 3 $FIXH ${OUT}labelinv.nii  -i  ${OUT}Affine.txt  ${OUT}InverseWarp.nii   -R $MOV --use-NN
