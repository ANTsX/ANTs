
NUMPARAMS=$#

echo $NUMPARAMS

exit
RADIUS=1
SIGMA=10
LAMBDA=10
if [ $NUMPARAMS -lt 3  ]
then
echo "  see http://www.insight-journal.org/browse/publication/306 for details on the Boykov algorithm implementation "
echo " USAGE ::  "
echo "  sh   antsbrainsegmentation.sh  ImageDimension imagein.nii Optional-Radius Optional-Sigma Optional-Lambda "
exit
fi

#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
then
echo " Please set ANTSPATH=LocationOfYourAntsBinaries "
echo " Either set this in your profile or directly, here in the script. "
echo " For example : "
echo " ANTSPATH=/home/yourname/bin/ants/ "
exit
else
echo " ANTSPATH is $ANTSPATH "
fi

#initialization, here, is unbiased
DIM=$1

if  [ ${#DIM} -gt 1 ]
then
echo " Problem with specified ImageDimension => User Specified Value = $DIM "
exit
fi

IMG=$1
NTISSUES=3
PRE=` echo $IMG | cut -d '.' -f 1 `
EXT=` echo $IMG | cut -d '.' -f 2 `
SEGNM=segment
PRNM=prob
SEGOUT=${PRE}${SEGNM}Labels.nii.gz
BIASIMG=${PRE}bias.nii.gz
${ANTSPATH}N3BiasFieldCorrection $DIM $IMG $BIASIMG 4
${ANTSPATH}ThresholdImage $DIM $BIASIMG $SEGOUT Otsu $NTISSUES
for x in $(seq 1 $NTISSUES ) ; do
 PRIOROUT=${PRE}${x}${PRNM}.${EXT}
  ${ANTSPATH}ThresholdImage $DIM $SEGOUT $PRIOROUT $x $x
 ${ANTSPATH}SmoothImage $DIM $PRIOROUT 1 $PRIOROUT
done
${ANTSPATH}BoykovGraphCutFilter $DIM $BIASIMG $SEGOUT ${PRE}${SEGNM}  $NTISSUES   ${PRE}1${PRNM}.${EXT} ${PRE}2${PRNM}.${EXT} ${PRE}3${PRNM}.${EXT}  $SMOOTH $SIGMA $LAMBDA
exe=" ${ANTSPATH}MultiplyImages $DIM $SEGOUT 1  turd.hdr "
echo $exe
$exe

exit
