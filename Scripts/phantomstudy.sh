#!/bin/bash

ANTSPATH="/Users/stnava/Code/bin/ants/"

#initialization, here, is unbiased
count=0
for x in  phantomAwmgm.jpg phantomBwmgm.jpg phantomCwmgm.jpg phantomDwmgm.jpg phantomEwmgm.jpg phantomFwmgm.jpg phantomGwmgm.jpg phantomHwmgm.jpg
do
    count=`expr $count + 1`					# Increment the counter
   echo " count is $count  and file is $x "
    ${ANTSPATH}/ANTS 2 -m PR[ phantomtemplate.jpg,$x,1,4 ]   -i 20x161x161x161  -o TEST{$count}  -t SyN[ 2 ] -r Gauss[ 3,0.]
#
#  -t SyN[ 1,2,0.05 ] -r Gauss[ 3,0.25]
    ${ANTSPATH}/WarpImageMultiTransform 2  $x  TEST{$count}registered.nii  TEST{$count}Warp.nii TEST{$count}Affine.txt
done

# to look at the stack of results, do this:

${ANTSPATH}/StackSlices volo.hdr -1 -1 0 *$1
${ANTSPATH}/StackSlices volr.hdr -1 -1 0 *registered.nii

count=0
for x in  phantomAwmgm.jpg phantomBwmgm.jpg phantomCwmgm.jpg phantomDwmgm.jpg phantomEwmgm.jpg phantomFwmgm.jpg phantomGwmgm.jpg phantomHwmgm.jpg
do
    count=`expr $count + 1`					# Increment the counter
   echo " count is $count  and file is $x "
 echo "    CreateJacobianDeterminantImage  2   TEST{$count}Warp.nii  TEST{$count}logjac.nii  1 "
   ${ANTSPATH}/CreateJacobianDeterminantImage  2   TEST{$count}Warp.nii  TEST{$count}logjac.nii  1
#    ${ANTSPATH}/SmoothImage  2 TEST{$count}logjac.nii  1  TEST{$count}logjac.nii
done


${ANTSPATH}/ThresholdImage 2 phantomtemplate.jpg mask.nii 100 200
#  use the GLM with mask
 ${ANTSPATH}/GLM 2  mask.nii designmat.txt contrast.txt Result.nii 1000 TEST*logjac.nii
