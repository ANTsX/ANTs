### change this path to your local path
ROOT=/home/hwang3/MA2012

R=2
S=3
sig=2
lambda=0.1
IDs=`cat $ROOT/Testing/Test`
TrainIDs=`cat $ROOT/Training/Train`

for id in $IDs ;do
	qsub -p -1 -o $ROOT/running -e $ROOT/running -cwd -N "lf_${R}${S}${sig}_${id}" -V -pe serial 4 -l h_stack=128M $ROOT/label_fusion_qsub.sh $id $R $S $sig $lambda; 
done

for id in $TrainIDs ;do
        qsub -p -1 -o $ROOT/running -e $ROOT/running -cwd -N "lf_${R}${S}${sig}_${id}" -V -pe serial 4 -l h_stack=128M $ROOT/label_fusion_train_qsub.sh $id $R $S $sig $lambda;
done

