#!/bin/bash

echo "Usage: \n sh  weightedaverage.sh  \"Faces*tiff\" \n "

count=9
for x in `ls -tr  $1 `
do
    count=`expr $count + 1`					# Increment the counter
   echo " count is $count  and file is $x at it $i "
      ANTS 2 -m PR[ template.nii,$x,1,8,-0.9 ]   -t SyN[ 3 ] -r Gauss[ 3,1.5 ] -o ROBU{$count}   -i 22x21x10  --MI-option 16x8000 #-a InitAffine.txt --continue-affine 0
      WarpImageMultiTransform 2  $x    ROBU{$count}{$i}registered.nii ROBU{$count}Warp.nii ROBU{$count}Affine.txt
ComposeMultiTransform 2   ROBU{$count}Warp.nii  -R template.nii ROBU{$count}Warp.nii ROBU{$count}Affine.txt
     ComposeMaps 2 ROBU{$count}Warp.nii ROBU{$count}Affine.txt ROBU{$count}Warp.nii  2
       ConvertToJpg  ROBU{$count}{$i}registered.nii ROBU{$count}{$i}registered.jpg
       MeasureImageSimilarity 2 2 template.nii ROBU{$count}{$i}registered.nii metriclog.txt
done

# end loop
