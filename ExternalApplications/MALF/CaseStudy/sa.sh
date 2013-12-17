### change this path to your local path
ROOT=/home/hwang3/MA2012

IDs=`cat $ROOT/Testing/Test`
TrainIDs=`cat $ROOT/Training/Train`

RP=$1
for id in $IDs; do
    qsub -p -1 -o $ROOT/running -e $ROOT/running -cwd -N "sa_${id}" -V -pe serial 1 -l h_stack=128M $ROOT/sa_qsub.sh $id $RP;
done
