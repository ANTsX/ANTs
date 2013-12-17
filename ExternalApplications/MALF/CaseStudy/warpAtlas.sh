### change these paths to your local paths
ROOT=/home/hwang3/MA2012 # path of the folder holding the 2012 multi atlas label challenge data

IDs0=`cat $ROOT/Testing/Test`
IDs=`cat $ROOT/Training/Train`

for id in $IDs ;do
    for ref in $IDs0;do 
        check=`ls $ROOT/warped/${id}_to_${ref}_image.nii.gz`
        if [[ $check != '' ]]; then
            echo $check 'has been processed!';
            continue;
        fi
	qsub -p -0 -o $ROOT/running -e $ROOT/running -cwd -N "warp_${ref}_${id}" -V -pe serial 1 -l h_stack=64M $ROOT/warpTest.sh $ref $id; 
    done
done


for id in $IDs ;do
    for ref in $IDs;do
        check=`ls $ROOT/warped/${id}_to_${ref}_image.nii.gz`
        if [[ $check != '' ]]; then
            echo $check 'has been processed!';
            continue;
        fi
        qsub -p -0 -o $ROOT/running -e $ROOT/running -cwd -N "warp_${ref}_${id}" -V -pe serial 1 -l h_stack=64M $ROOT/warp.sh $ref $id;
    done
done

