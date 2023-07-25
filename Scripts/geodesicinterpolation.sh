#!/bin/bash
NUMPARAMS=$#


if [ $NUMPARAMS -lt 2 ]
then
echo " USAGE ::  "
echo "  sh   geodesicinterpolation.sh image1  image2  N-interpolation-points  N-step-size Use-Affine-Initialization?  Mask-Image Auxiliary-Template "
echo " N-interpolation-points -- if you choose 10 points and step-size of 1 then you will get 8 interpolated/morphed images -- end-points are the original images "
echo " if you choose 10 points and step-size of  2 then you will get 4 interpolated/morphed images -- end-points are the original images "
echo " if you use affine initialization - sometimes it could cause unpleasantness, if the affine mapping is large -- sometimes it may be needed, though "
echo " the mask image restricts matching to an ROI "
echo " The Auxiliary-Template is another image through which the morph images (image1 and image2) are mapped -- this often helps, if the auxiliary template is an average face image "
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
if [ $NUMPARAMS -gt 2 ]
then
NUMSTEPS=$3
fi
if [ $NUMPARAMS -gt 3 ]
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
ANTS 2 -m PR[ $TEMPLATE,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o $OUTNAME   -i $ITS    --number-of-affine-iterations 0    -x $MASK
fi
if [ $USEAFF -ne 0 ]
then
ANTS 2 -m PR[ $TEMPLATE,$TARGET,1,${RADIUS}]   -t  $TRAN -r $REG -o $OUTNAME   -i $ITS     -x  $MASK  #--MI-option 16x2000
fi
ComposeMultiTransform 2   ${OUTNAME}Warp.nii.gz  -R $TEMPLATE ${OUTNAME}Warp.nii.gz ${OUTNAME}Affine.txt
ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii.gz  -R $TARGET -i ${OUTNAME}Affine.txt ${OUTNAME}InverseWarp.nii.gz
fi

if [ $NUMPARAMS -gt 5 ]
then
if [ $USEAFF -eq 0 ]
then
#echo "  Pseudo-Morphing "
#echo " method 1 "
ANTS 2 -m PR[ $TEMPLATEB,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o ${OUTNAME}B  -i $ITS    --number-of-affine-iterations 0    -x $MASK
ANTS 2 -m PR[ $TEMPLATE,$TEMPLATEB,1,${RADIUS}]   -t  $TRAN -r $REG -o ${OUTNAME}A   -i $ITS  --number-of-affine-iterations 0    -x $MASK
ComposeMultiTransform 2   ${OUTNAME}AWarp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}AAffine.txt
ComposeMultiTransform 2   ${OUTNAME}AInverseWarp.nii.gz  -R $TEMPLATEB -i ${OUTNAME}AAffine.txt ${OUTNAME}AInverseWarp.nii.gz
ComposeMultiTransform 2   ${OUTNAME}BWarp.nii.gz  -R $TEMPLATEB  ${OUTNAME}BWarp.nii.gz ${OUTNAME}BAffine.txt
ComposeMultiTransform 2   ${OUTNAME}BInverseWarp.nii.gz  -R $TARGET -i ${OUTNAME}BAffine.txt ${OUTNAME}BInverseWarp.nii.gz
ComposeMultiTransform 2   ${OUTNAME}Warp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}BWarp.nii.gz
ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii.gz  -R $TARGET  ${OUTNAME}BInverseWarp.nii.gz ${OUTNAME}AInverseWarp.nii.gz
#
#echo " method 2 "
#ComposeMultiTransform 2   ${OUTNAME}Warp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}AAffine.txt  ${OUTNAME}BWarp.nii.gz ${OUTNAME}BAffine.txt
#ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii.gz  -R $TARGET -i ${OUTNAME}BAffine.txt   ${OUTNAME}BInverseWarp.nii.gz -i ${OUTNAME}AAffine.txt  ${OUTNAME}AInverseWarp.nii
#
fi
if [ $USEAFF -ne 0 ]
then
#echo "  Pseudo-Morphing "
#echo " method 1 "
ANTS 2 -m PR[ $TEMPLATEB,$TARGET,1,${RADIUS}]   -t $TRAN -r $REG -o ${OUTNAME}B  -i $ITS    -x $MASK
ANTS 2 -m PR[ $TEMPLATE,$TEMPLATEB,1,${RADIUS}]   -t  $TRAN -r $REG -o ${OUTNAME}A   -i $ITS     -x $MASK
ComposeMultiTransform 2   ${OUTNAME}AWarp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}AAffine.txt
ComposeMultiTransform 2   ${OUTNAME}AInverseWarp.nii.gz  -R $TEMPLATEB -i ${OUTNAME}AAffine.txt ${OUTNAME}AInverseWarp.nii.gz
ComposeMultiTransform 2   ${OUTNAME}BWarp.nii.gz  -R $TEMPLATEB  ${OUTNAME}BWarp.nii.gz ${OUTNAME}BAffine.txt
ComposeMultiTransform 2   ${OUTNAME}BInverseWarp.nii.gz  -R $TARGET -i ${OUTNAME}BAffine.txt ${OUTNAME}BInverseWarp.nii.gz
ComposeMultiTransform 2   ${OUTNAME}Warp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}BWarp.nii
ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii.gz  -R $TARGET  ${OUTNAME}BInverseWarp.nii.gz ${OUTNAME}AInverseWarp.nii.gz
#
#echo " method 2 "
#ComposeMultiTransform 2   ${OUTNAME}Warp.nii.gz  -R $TEMPLATE  ${OUTNAME}AWarp.nii.gz ${OUTNAME}AAffine.txt  ${OUTNAME}BWarp.nii.gz ${OUTNAME}BAffine.txt
#ComposeMultiTransform 2   ${OUTNAME}InverseWarp.nii.gz  -R $TARGET -i ${OUTNAME}BAffine.txt   ${OUTNAME}BInverseWarp.nii.gz -i ${OUTNAME}AAffine.txt  ${OUTNAME}AInverseWarp.nii
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

MultiplyImages 2 ${OUTNAME}InverseWarp.nii.gz  $BLENDINGA  SM${OUTNAME}InverseWarp.nii.gz
MultiplyImages 2 ${OUTNAME}InverseWarp.nii.gz  $BLENDINGA  SM${OUTNAME}InverseWarp.nii.gz

MultiplyImages 2 ${OUTNAME}Warp.nii.gz  $BLENDINGB SM${OUTNAME}Warp.nii.gz
MultiplyImages 2 ${OUTNAME}Warp.nii.gz  $BLENDINGB SM${OUTNAME}Warp.nii.gz

  WarpImageMultiTransform 2 $TARGET temp.nii.gz SM${OUTNAME}Warp.nii.gz  -R $TEMPLATE
  ImageMath 2 temp.nii.gz Normalize temp.nii.gz 1
  ImageMath 2 temp.nii.gz m temp.nii.gz 1.   #$BLENDINGA
  WarpImageMultiTransform 2 $TEMPLATE temp2.nii.gz   SM${OUTNAME}InverseWarp.nii.gz -R $TEMPLATE
  ImageMath 2 temp2.nii.gz Normalize temp2.nii.gz  1
  ImageMath 2 temp2.nii.gz m temp2.nii.gz 0  #$BLENDINGB
  echo "  ImageMath 2 ${BASEA}${BASEB}${BLENDNAME}morph.nii.gz + temp2.nii.gz temp.nii.gz  "
  ImageMath 2 ${BASEA}${BASEB}${BLENDNAME}morph.nii.gz + temp2.nii.gz temp.nii.gz
  ConvertToJpg  ${BASEA}${BASEB}${BLENDNAME}morph.nii.gz ${BASEA}${BASEB}${BLENDNAME}morph.jpg
  rm  -f  ${BASEA}${BASEB}${BLENDNAME}morph.nii.gz

   done

rm -f  ${OUTNAME}*  SM${OUTNAME}*

exit
