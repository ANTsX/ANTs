#!/bin/bash
NUMPARAMS=$#


if [ $NUMPARAMS -lt 2  ]
then
echo " USAGE ::  "
echo "  sh   geodesicinterpolation.sh image1  image2  N-interpolation-points  N-step-size Use-Affine-Initialization?  Mask-Image Auxiliary-Template "
echo " N-interpolation-points -- if you choose 10 points and step-size of 1 then you will get 8 interpolated/morphed images -- end-points are the original images "
echo " if you choose 10 points and step-size of  2 then you will get 4 interpolated/morphed images -- end-points are the original images "
echo " if you use affine initialization - sometimes it could cause unpleasantness, if the affine mapping is large -- sometimes it may be needed, though "
echo " the mask image restricts matching to an ROI "
echo " The Auxiliary-Template is another image through which the morph images (image1 and image2) are mapped -- this often helps, if the auxiliary template is an average face image "
echo " Be sure to set the ANTSPATH path variable "
echo " "
echo " You may need to tune the ANTS parameters coded within this script for your application "
exit
fi

# should we use the generative model to find landmarks then
#  drive individual mappings by the LMs?

TEMPLATE=$1
TARGET=$2
NUMSTEPS=10
STEP=1
OUTNAME=ANTSMORPH
if [ $NUMPARAMS -gt 2  ]
then
NUMSTEPS=$3
fi
if [ $NUMPARAMS -gt 3  ]
then
STEP=$4
fi

USEAFF=0
if [ $NUMPARAMS -gt 4 ]
then
USEAFF=$5
fi

MASK=0
if [ $NUMPARAMS -gt 5 ]
then
MASK=$6
fi

TEMPLATEB=0
if [ $NUMPARAMS -gt 6 ]
then
TEMPLATEB=$7
fi

echo " Morphing in total $NUMSTEPS time-points by $STEP steps -- end-points are original images "
echo " Use-Affine ?  $USEAFF "
echo " templateB $TEMPLATEB "


    for (( n = 0 ; n <= ${NUMSTEPS}; n=n+${STEP} ))
    do
BLENDINGA=$(echo "scale=2; ${n}/${NUMSTEPS}"  | bc )
BLENDINGB=$(echo "scale=2;  1-$BLENDINGA"  | bc )
BLENDNAME=$(echo "scale=2;  100*$BLENDINGA"  | bc )
echo " Blending values:   $BLENDINGA and $BLENDINGB"
  done

ITS=100x100x10
TRAN=SyN[0.5]
REG=Gauss[3,1]
RADIUS=6

if [ $NUMPARAMS -le 5 ]
then
if [ $USEAFF -eq 0 ]
then
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATE,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o $OUTNAME   -i $ITS    --number-of-affine-iterations 0    -x $MASK
fi
if [ $USEAFF -ne 0 ]
then
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATE,$TARGET,1,${RADIUS}]   -t  $TRAN -r $REG -o $OUTNAME   -i $ITS     -x  $MASK  #--MI-option 16x2000
fi
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}Warp.nii  -R $TEMPLATE ${OUTNAME}Warp.nii ${OUTNAME}Affine.txt
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii  -R $TARGET -i ${OUTNAME}Affine.txt ${OUTNAME}InverseWarp.nii
fi

if [ $NUMPARAMS -gt 5 ]
then
if [ $USEAFF -eq 0 ]
then
#echo "  Pseudo-Morphing "
#echo " method 1 "
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATEB,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o ${OUTNAME}B  -i $ITS    --number-of-affine-iterations 0    -x $MASK
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATE,$TEMPLATEB,1,${RADIUS}]   -t  $TRAN -r $REG -o ${OUTNAME}A   -i $ITS  --number-of-affine-iterations 0    -x $MASK
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}AWarp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}AAffine.txt
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}AInverseWarp.nii  -R $TEMPLATEB -i ${OUTNAME}AAffine.txt ${OUTNAME}AInverseWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}BWarp.nii  -R $TEMPLATEB  ${OUTNAME}BWarp.nii ${OUTNAME}BAffine.txt
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}BInverseWarp.nii  -R $TARGET -i ${OUTNAME}BAffine.txt ${OUTNAME}BInverseWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}Warp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}BWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii  -R $TARGET  ${OUTNAME}BInverseWarp.nii ${OUTNAME}AInverseWarp.nii
#
#echo " method 2 "
#${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}Warp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}AAffine.txt  ${OUTNAME}BWarp.nii ${OUTNAME}BAffine.txt
#${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii  -R $TARGET -i ${OUTNAME}BAffine.txt   ${OUTNAME}BInverseWarp.nii -i ${OUTNAME}AAffine.txt  ${OUTNAME}AInverseWarp.nii
#
fi
if [ $USEAFF -ne 0 ]
then
#echo "  Pseudo-Morphing "
#echo " method 1 "
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATEB,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o ${OUTNAME}B  -i $ITS    -x $MASK
${ANTSPATH}/ANTS 2 -m PR[$TEMPLATE,$TEMPLATEB,1,${RADIUS}]   -t  $TRAN -r $REG -o ${OUTNAME}A   -i $ITS     -x $MASK
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}AWarp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}AAffine.txt
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}AInverseWarp.nii  -R $TEMPLATEB -i ${OUTNAME}AAffine.txt ${OUTNAME}AInverseWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}BWarp.nii  -R $TEMPLATEB  ${OUTNAME}BWarp.nii ${OUTNAME}BAffine.txt
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}BInverseWarp.nii  -R $TARGET -i ${OUTNAME}BAffine.txt ${OUTNAME}BInverseWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}Warp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}BWarp.nii
${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii  -R $TARGET  ${OUTNAME}BInverseWarp.nii ${OUTNAME}AInverseWarp.nii
#
#echo " method 2 "
#${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}Warp.nii  -R $TEMPLATE  ${OUTNAME}AWarp.nii ${OUTNAME}AAffine.txt  ${OUTNAME}BWarp.nii ${OUTNAME}BAffine.txt
#${ANTSPATH}/ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii  -R $TARGET -i ${OUTNAME}BAffine.txt   ${OUTNAME}BInverseWarp.nii -i ${OUTNAME}AAffine.txt  ${OUTNAME}AInverseWarp.nii
#
fi
fi

    for (( n = 0 ; n <= ${NUMSTEPS}; n=n+${STEP} ))
    do

BLENDINGA=$(echo "scale=2; ${n}/${NUMSTEPS}"  | bc )
BLENDINGB=$(echo "scale=2;  1-$BLENDINGA"  | bc )
BLENDNAME=$(echo "scale=2;  100*$BLENDINGA"  | bc )
echo " Blending values:   $BLENDINGA and $BLENDINGB"

BASEA=${TEMPLATE%.*}
BASEB=${TARGET%.*}

${ANTSPATH}/MultiplyImages 2 ${OUTNAME}InverseWarpxvec.nii  $BLENDINGA  SM${OUTNAME}InverseWarpxvec.nii
${ANTSPATH}/MultiplyImages 2 ${OUTNAME}InverseWarpyvec.nii  $BLENDINGA  SM${OUTNAME}InverseWarpyvec.nii

${ANTSPATH}/MultiplyImages 2 ${OUTNAME}Warpxvec.nii  $BLENDINGB SM${OUTNAME}Warpxvec.nii
${ANTSPATH}/MultiplyImages 2 ${OUTNAME}Warpyvec.nii  $BLENDINGB SM${OUTNAME}Warpyvec.nii

  ${ANTSPATH}/WarpImageMultiTransform 2 $TARGET temp.nii SM${OUTNAME}Warp.nii  -R $TEMPLATE
  ${ANTSPATH}/ImageMath 2 temp.nii Normalize temp.nii 1
  ${ANTSPATH}/ImageMath 2 temp.nii m temp.nii 1.   #$BLENDINGA
  ${ANTSPATH}/WarpImageMultiTransform 2 $TEMPLATE temp2.nii   SM${OUTNAME}InverseWarp.nii -R $TEMPLATE
  ${ANTSPATH}/ImageMath 2 temp2.nii Normalize temp2.nii  1
  ${ANTSPATH}/ImageMath 2 temp2.nii m temp2.nii 0  #$BLENDINGB
  echo "  ImageMath 2 ${BASEA}${BASEB}${BLENDNAME}morph.nii + temp2.nii temp.nii  "
  ${ANTSPATH}/ImageMath 2 ${BASEA}${BASEB}${BLENDNAME}morph.nii + temp2.nii temp.nii
  ${ANTSPATH}/ConvertToJpg  ${BASEA}${BASEB}${BLENDNAME}morph.nii ${BASEA}${BASEB}${BLENDNAME}morph.jpg
  rm  -f  ${BASEA}${BASEB}${BLENDNAME}morph.nii

   done

rm -f  ${OUTNAME}*  SM${OUTNAME}*

exit
