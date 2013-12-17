BIN=/home/hwang3/ahead/turnkey/bin

### change this path to your local path
ROOT=/home/hwang3/MA2012 
labels=`cat $ROOT/labels`
labels1=`cat $ROOT/labels1`

rates=`cat $ROOT/samplerates1`
DR=1
IDs=`cat $ROOT/Testing/Test`
TrainIDs=`cat $ROOT/Training/Train`

c=0
for L in $labels ;do
    let c=c+1
    c1=0
    for r in $rates; do 
        let c1=c1+1
        if [ $c == $c1 ]; then
            rate=$r
        fi
    done
    c1=0
    for l in $labels1; do
        let c1=c1+1
        if [ $c == $c1 ]; then
            L2=$l
        fi
    done

    echo $c $rate $L2
    adaboostfile=$ROOT/malf/BL/JointLabel_BL-AdaBoostResults-Tlabel$L2
 
    if [ -f $adaboostfile ]; then
        continue;  
    fi
    qsub -p -1 -o $ROOT/running -e $ROOT/running -cwd -N "bl_${L}" -V -pe serial 4 -l h_stack=128M $ROOT/bl_qsub.sh $L $DR $rate 2 ; 
done

