# !/bin/sh

echo "sh $0 BrainIn=Image.nii PriorImage=Prior.nii PriorLocalWT=0 MRF=0.3  LOOPS=2 MASK=image.nii "
if [ $# -lt 5 ]
then
exit
fi

BRAIN=$1
OUT=` echo $BRAIN | cut -d '.' -f 1  `
echo $OUT
PRIOR=$2
NCLASS=` MeasureMinMaxMean 3 $2 | cut -d ':' -f 2 | cut -d ' ' -f 2 `
WT=$3
MRF=$4
BSP=1
BSPLEV=1
ITS=$5
MASK=$6
#MASK=` echo $1 | cut -d '.' -f 1 `
#MASK=${MASK}mask.nii.gz
if [ ! -s $MASK ] ; then
echo " No mask $MASK "
exit
fi
if [ ! -s $BRAIN ] ; then
echo " No brain $BRAIN "
exit
fi

if [ ! -s $PRIOR ] ; then
echo " No prior $PRIOR "
exit
fi
if [ ! -s $MASK  ] ; then
echo " No mask  "
ThresholdImage 3 $BRAIN $MASK 0.01 99999
fi
dd=3.
ddloc=1.8
csf=10
wm=${csf}
dist=" "
for x in $(seq 1 $NCLASS) ; do
dist=" $dist  -l ${x}[${dd}] "
done

#if [ ${#APOC} -lt 1 ] ; then
#echo " JERK !! "
#exit
#fi
dist2=" -l 1[${dd}] -l 2[${dd}] -l 3[${dd}] -l 4[${dd}] -l 5[${dd}] -l 6[${dd}] -l 7[${dd}] -l 8[${dd}] -l 9[${dd}] -l 10[${dd}] -l 11[${dd}] -l 12[${dd}] -l 13[${dd}] -l 14[${dd}] -l 15[${dd}] -l 16[${dd}] -l 17[${dd}] -l 18[${dd}] -l 19[${dd}] -l 20[${dd}] -l 21[${dd}] -l 22[${dd}] -l 23[${dd}] -l 24[${dd}] -l 25[${dd}] -l 26[${dd}] -l 27[${dd}] -l 28[${dd}] -l 29[${dd}] -l 30[${dd}] -l 31[${dd}] -l 32[${dd}] -l 33[${ddloc}] -l 34[${ddloc}] -l 35[${ddloc}] -l 36[${ddloc}] -l 37[${ddloc}] -l 38[${ddloc}] -l 39[${ddloc}] -l 40[${ddloc}] -l 41[${ddloc}] -l 42[${ddloc}] -l 43[${ddloc}] -l 44[${ddloc}] -l 45[${wm}] -l 46[${csf}]  "
SEG="${ANTSPATH}Apocrita 3 "
exe="$SEG -i PriorLabelImage[${BRAIN},${NCLASS},${PRIOR},${WT}] -n $ITS -c 0 -x $MASK  -m [${MRF},1,0,0]  -o [${OUT}${ITS}_seg.nii.gz,outPosterior%02d.nii.gz] -h 0 $dist2 "

BSP=" -b [${BSPLEV},${BSP}x${BSP}x${BSP},3] "
#exe="$SEG -i Kmeans[${BRAIN},${NCLASS}] -n 2 -c 0 -x $MASK  -m [$1,1,0.1,0.25] -b [4,8x8,3] -o [out.nii.gz,outPosterior%02d.nii.gz,outBSpline%02d.nii.gz] -h 0  "
echo $exe

$exe
echo $exe

exit

# To compare the final results, you can compare the slice segmentation from before %

# ~/snap -g $BRAIN -s $PRIOR  &

# and after

# ~/snap -g $BRAIN -s out.nii.gz &



