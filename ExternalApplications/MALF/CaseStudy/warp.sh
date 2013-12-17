#!/bin/bash
#$ -S /bin/bash
set -x -e
ref=$1
mov=$2
### change these paths to your local paths
ROOT=/home/hwang3/MA2012 # path of the folder holding the 2012 multi atlas label challenge data

$ROOT/antsApplyTransforms -d 3 -r $ROOT/Training/${ref}.nii.gz -i $ROOT/Training/${mov}.nii.gz -o $ROOT/warped/${mov}_to_${ref}_image.nii.gz -t $ROOT/displacement-fields-training-to-training/${ref}x${mov}TotalWarp.nii.gz

$ROOT/antsApplyTransforms -d 3 -r $ROOT/Training/${ref}_glm.nii.gz -i $ROOT/Training/${mov}_glm.nii.gz -o $ROOT/warped/${mov}_to_${ref}_seg_NN.nii.gz -t $ROOT/displacement-fields-training-to-training/${ref}x${mov}TotalWarp.nii.gz -n NearestNeighbor
