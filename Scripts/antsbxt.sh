# !/bin/sh

NUMPARAMS=$#

MAXITERATIONS=30x90x20

if [ $NUMPARAMS -lt 3  ]
then
echo " USAGE ::  "
echo "  sh   ants.sh  ImageDimension  fixed.ext  movingtostrip.ext  OUTPREFIX  max-iterations   BrainMask  param "
echo " be sure to set ANTSPATH environment variable "
echo " Max-Iterations in form :    JxKxL where "
echo "  J = max iterations at coarsest resolution (here, reduce by power of 2^2) "
echo " K = middle resolution iterations ( here, reduce by power of 2 ) "
echo " L = fine resolution iterations ( here, full resolution ) -- this level takes much more time per iteration "
echo " an extra  Ix before JxKxL would add another level "
echo " Default ierations is  $MAXITERATIONS  -- you can often get away with fewer for many apps "
echo " Other parameters are defaults used in the A. Klein evaluation paper in Neuroimage, 2009 "
echo " "
echo " experimental skull stripping -- requires brain mask as input "
echo " param should be in range -0.05  to -0.2   --  -0.1 is a good start "
echo " "
echo " we perform mapping and then mask the brain from skull -- then segment the image and try to use the wm to extract the cerebrum "
echo " CURRENTLY BEING EVALUATED! "
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

FIXED=$2

if  [ ${#FIXED} -lt 1 -o  ! -f $FIXED ]
then
echo " Problem with specified Fixed Image => User Specified Value = $FIXED "
exit
fi

MOVING=$3

if  [ ${#MOVING} -lt 1 -o  ! -f $MOVING ]
then
echo " Problem with specified Moving Image => User Specified Value = $MOVING "
exit
fi

OUTPUTNAME=` echo $MOVING | cut -d '.' -f 1 `


if [ $NUMPARAMS -gt 4  ]
then
 MAXITERATIONS=$5
fi

LABELIMAGE=0
if [ $NUMPARAMS -gt 5  ]
then
 LABELIMAGE=$6
fi

DoANTSQC=-0.1
if [ $NUMPARAMS -gt 6  ]
then
 DoANTSQC=$7
fi

if [ $NUMPARAMS -gt 3  ]
then
  OUTPUTNAME=${4}
fi

# Mapping Parameters
  TRANSFORMATION=SyN[0.25]
  ITERATLEVEL=(`echo $MAXITERATIONS | tr 'x' ' '`)
  NUMLEVELS=${#ITERATLEVEL[@]}
  echo $NUMLEVELS
  REGULARIZATION=Gauss[3,0.]
  METRIC=PR[
    METRICPARAMS=1,2]
#echo " $METRICPARAMS  &  $METRIC "
#exit


echo  " ANTSPATH  is $ANTSPATH     "
echo " Mapping Parameters  :: "
echo  " Transformation is:  $TRANSFORMATION "
echo " MaxIterations :   $MAXITERATIONS "
echo " Number Of MultiResolution Levels   $NUMLEVELS "
echo " Regularization :  $REGULARIZATION "
echo " Metric :  ${METRIC} "
echo " OutputName :  $OUTPUTNAME "
echo " "
echo " if the files and parameters are all ok then uncomment the exit call below this line  "
echo " "
#exit


# finally check the image headers
compareheaders=`${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}repaired.nii.gz  CompareHeadersAndImages $FIXED $MOVING  | grep FailureState | cut -d ' ' -f 4  `
if [ $compareheaders -ne 0 ]
then
echo " You may have a problem with your header definition "
echo " The repaired image is in : ${OUTPUTNAME}repaired.nii.gz "
echo " Call ImageMath's   CompareHeadersAndImages  on your Fixed and Moving image "
exit
fi
BRN=${OUTPUTNAME}brain.nii.gz

dothis=1
if [ $dothis -eq 1 ] ; then
#${ANTSPATH}N3BiasFieldCorrection $DIM ${OUTPUTNAME}repaired.nii.gz  ${OUTPUTNAME}repaired.nii.gz 8
#${ANTSPATH}MultiplyImages $DIM ${OUTPUTNAME}repaired.nii.gz  1 ${OUTPUTNAME}repaired.nii.gz

exe=" ${ANTSPATH}ANTS $DIM -m  ${METRIC}${FIXED},${OUTPUTNAME}repaired.nii.gz,${METRICPARAMS}  -t $TRANSFORMATION  -r $REGULARIZATION -o ${OUTPUTNAME}.nii.gz   -i $MAXITERATIONS   --use-Histogram-Matching "

 echo " $exe "

#  $exe

  #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

    ${ANTSPATH}WarpImageMultiTransform $DIM  ${OUTPUTNAME}repaired.nii.gz   ${OUTPUTNAME}deformed.nii.gz ${OUTPUTNAME}Warp.nii.gz ${OUTPUTNAME}Affine.txt  -R ${FIXED}

BRN=${OUTPUTNAME}brain.nii.gz
if  [ ${#LABELIMAGE} -gt 3 ]
then
      ${ANTSPATH}WarpImageMultiTransform $DIM  $LABELIMAGE   ${OUTPUTNAME}labeled.nii.gz -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz  -R ${MOVING}   --use-NN
      ${ANTSPATH}MultiplyImages $DIM ${OUTPUTNAME}labeled.nii.gz ${OUTPUTNAME}repaired.nii.gz $BRN
fi


#      ${ANTSPATH}WarpImageMultiTransform $DIM  PCEtemplate_session2_brain_prob_0.nii.gz   prior0.nii.gz -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz  -R ${MOVING}
 #     ${ANTSPATH}WarpImageMultiTransform $DIM  PCEtemplate_session2_brain_prob_1.nii.gz   prior1.nii.gz -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz  -R ${MOVING}
 #     ${ANTSPATH}WarpImageMultiTransform $DIM  PCEtemplate_session2_brain_prob_2.nii.gz   prior2.nii.gz -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz  -R ${MOVING}

#${ANTSPATH}ImageMath $DIM $BRN Segment $BRN 4 0.0 # prior0.nii.gz prior1.nii.gz prior2.nii.gz


fi # dothis

SIG2=${DoANTSQC}
SUB=${OUTPUTNAME}brain
PWM=${SUB}_prob_2.nii.gz
PGM=${SUB}_prob_1.nii.gz
PCS=${SUB}_prob_0.nii.gz
BWM=${SUB}2b.nii.gz
BGM=${SUB}1b.nii.gz
BCS=${SUB}_brain0b.nii.gz

BMK=${OUTPUTNAME}labeled.nii.gz

ThresholdImage 3 $BRN $BMK 0.01 999
ThresholdImage 3 $PWM $BWM 0.5 999
ThresholdImage 3 $PGM $BGM 0.5 999
ThresholdImage 3 $PCS $BCS 0.5 999

#first, get the largest wm component
ImageMath 3 ${OUTPUTNAME}temp.nii.gz ME $BWM 1
ImageMath 3 ${OUTPUTNAME}temp.nii.gz GetLargestComponent ${OUTPUTNAME}temp.nii.gz
# now extract the csf components and dilate them
ThresholdImage 3 $PCS ${OUTPUTNAME}tempc.nii.gz 0.25 999

ImageMath 3 ${OUTPUTNAME}tempc.nii.gz MD ${OUTPUTNAME}tempc.nii.gz 1
ImageMath 3 ${OUTPUTNAME}tempc.nii.gz Neg ${OUTPUTNAME}tempc.nii.gz
ImageMath 3 ${OUTPUTNAME}tempc2.nii.gz Neg $PCS

# get the largest part of the brain mask
ImageMath 3 $BMK GetLargestComponent $BMK
MultiplyImages 3 ${OUTPUTNAME}tempc.nii.gz $BMK ${OUTPUTNAME}tempc.nii.gz
MultiplyImages 3 ${OUTPUTNAME}tempc.nii.gz ${OUTPUTNAME}tempc2.nii.gz ${OUTPUTNAME}tempc.nii.gz

# allow the wm component to propagate through the non-csf mask
ImageMath 3 ${OUTPUTNAME}temp2.nii.gz PropagateLabelsThroughMask ${OUTPUTNAME}tempc.nii.gz ${OUTPUTNAME}temp.nii.gz 30
# now get rid of big bits that are separate from the WM
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz ME ${OUTPUTNAME}temp2.nii.gz 1
 ImageMath 3 ${OUTPUTNAME}temp3.nii.gz GetLargestComponent ${OUTPUTNAME}temp3.nii.gz
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz MD ${OUTPUTNAME}temp3.nii.gz 1
ImageMath 3 ${OUTPUTNAME}temp4.nii.gz - ${OUTPUTNAME}temp2.nii.gz ${OUTPUTNAME}temp3.nii.gz
LabelClustersUniquely 3 ${OUTPUTNAME}temp4.nii.gz ${OUTPUTNAME}temp4.nii.gz 250
MeasureMinMaxMean 3 ${OUTPUTNAME}temp4.nii.gz
ThresholdImage 3 ${OUTPUTNAME}temp4.nii.gz ${OUTPUTNAME}temp4.nii.gz 1 1000
# above, the 2 largest components are taken away -- repropagate within remainder
ImageMath 3  ${OUTPUTNAME}tempc.nii.gz - ${OUTPUTNAME}tempc.nii.gz ${OUTPUTNAME}temp4.nii.gz
ImageMath 3 ${OUTPUTNAME}tempc.nii.gz abs ${OUTPUTNAME}tempc.nii.gz 1
ImageMath 3 ${OUTPUTNAME}temp2.nii.gz PropagateLabelsThroughMask ${OUTPUTNAME}tempc.nii.gz ${OUTPUTNAME}temp.nii.gz 500
ImageMath 3 ${OUTPUTNAME}temp2.nii.gz MD  ${OUTPUTNAME}temp2.nii.gz 1
# get the brain mask depth
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz Neg $BMK
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz D  ${OUTPUTNAME}temp3.nii.gz
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz exp ${OUTPUTNAME}temp3.nii.gz $SIG2
ImageMath 3 ${OUTPUTNAME}temp3.nii.gz Neg ${OUTPUTNAME}temp3.nii.gz
ThresholdImage 3 ${OUTPUTNAME}temp3.nii.gz ${OUTPUTNAME}temp3.nii.gz 0.5 9999
ImageMath 3 ${OUTPUTNAME}temp2.nii.gz addtozero ${OUTPUTNAME}temp2.nii.gz ${OUTPUTNAME}temp3.nii.gz
ImageMath 3 ${OUTPUTNAME}temp2.nii.gz MD ${OUTPUTNAME}temp2.nii.gz 1
MultiplyImages 3 $BRN ${OUTPUTNAME}temp2.nii.gz $BRN

rm -f ${OUTPUTNAME}temp*

exit
