#!/bin/bash -x
#$ -S /bin/bash

### change these paths to your local paths
ROOT=/home/hwang3/MA2012 # path of the folder holding the 2012 multi atlas label challenge data
BIN=/home/hwang3/MA2012/PICSL_MALF # path to the PICSL_MALF software


LABELS=`cat $ROOT/labels`
L=$1
DR=$2
rate=$3
rp=$4

IDs=`cat $ROOT/Training/Train`

feature=""
seg=""
manual=""
for id in $IDs; do
        if [[ $id == $ID ]]; then
            continue;
        fi
	feature=${feature}" "$ROOT/Training/${id}.nii.gz" "$ROOT/malf/${id}_JointLabel_posterior0$L.nii.gz
	seg=${seg}" "$ROOT/malf/${id}_JointLabel.nii.gz
        manual=${manual}" "$ROOT/Training/${id}_glm.nii.gz
done
echo $im
echo $seg
mkdir $ROOT/malf/BL
output=$ROOT/malf/BL/JointLabel_BL


$BIN/bl 3 -ms $manual -as $seg -tl $L -rd $DR -i 500 -rate $rate -rf ${rp}x${rp}x${rp} -c 2 -f $feature $output

