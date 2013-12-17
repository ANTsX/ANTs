#!/bin/bash -x
#$ -S /bin/bash

### change these paths to your local paths
ROOT=/home/hwang3/MA2012 # path of the folder holding the 2012 multi atlas label challenge data
BIN=/home/hwang3/MA2012/PICSL_MALF # path to the PICSL_MALF software
ID=$1
R=$2
S=$3
sig=$4
lamda=$5

IDs=`cat $ROOT/Training/Train`

im=""
seg=""

for id in $IDs; do
        if [[ $id == $ID ]]; then
            continue;
        fi

	im=${im}" "$ROOT/warped/${id}_to_${ID}_image.nii.gz
	seg=${seg}" "$ROOT/warped/${id}_to_${ID}_seg_NN.nii.gz
done
echo $im
echo $seg
mkdir $ROOT/malf


/home/hwang3/pkg/PICSL_MALF/jointfusion 3 1 -g $im -l $seg -m Joint[$lamda,$sig] -rp ${R}x${R}x${R} -rs ${S}x${S}x${S} -tg $ROOT/Training/$ID.nii.gz $ROOT/malf/${ID}_JointLabel.nii.gz

