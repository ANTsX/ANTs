#!/bin/bash

TNAME=$1
GRADSTEP=$2
shift 2

FLIST2=$*

ImageMath 2 $TNAME Normalize $TNAME
for x in $(seq 1 5 ) ; do
  for y in  $FLIST2 ; do
      ImageMath 2 turdx${y} Normalize $y
   ImageMath 2 turdx${y} CorrelationUpdate  $TNAME turdx${y} 4
#   MeasureImageSimilarity 2 1 $y  $TNAME
   SmoothImage 2 turdx${y} 0.5  turdx${y}
#  MeasureMinMaxMean 2 turdx$y
 done

  AverageImages 2 turdu.nii.gz 0 turdx**
   MultiplyImages 2  turdu.nii.gz  $GRADSTEP turdu.nii.gz
   ImageMath 2  $TNAME +  $TNAME turdu.nii.gz
done
