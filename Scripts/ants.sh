# !/bin/sh

NUMPARAMS=$#

MAXITERATIONS=30x90x20

if [ $NUMPARAMS -lt 3  ]
then
echo " USAGE ::  "
echo "  sh   ants.sh  ImageDimension  fixed.ext  moving.ext  OPTIONAL-OUTPREFIX   OPTIONAL-max-iterations  OPTIONAL-Labels-In-Fixed-Image-Space-To-Deform-To-Moving-Image     Option-DoANTSQC "
echo " be sure to set ANTSPATH environment variable "
echo " Max-Iterations in form :    JxKxL where "
echo "  J = max iterations at coarsest resolution (here, reduce by power of 2^2) "
echo " K = middle resolution iterations ( here, reduce by power of 2 ) "
echo " L = fine resolution iterations ( here, full resolution ) -- this level takes much more time per iteration "
echo " an extra  Ix before JxKxL would add another level "
echo " Default ierations is  $MAXITERATIONS  -- you can often get away with fewer for many apps "
echo " Other parameters are updates of the defaults used in the A. Klein evaluation paper in Neuroimage, 2009 "
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

DoANTSQC=0
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
  REGULARIZATION=Gauss[3,0]
  METRIC=PR[
    METRICPARAMS=1,4]
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
compareheaders=`${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}repaired.nii  CompareHeadersAndImages $FIXED $MOVING  | grep FailureState | cut -d ' ' -f 4  `
if [ $compareheaders -ne 0 ]
then
echo " You may have a problem with your header definition "
echo " The repaired image is in : ${OUTPUTNAME}repaired.nii "
echo " Call ImageMath's   CompareHeadersAndImages  on your Fixed and Moving image "
exit
fi

${ANTSPATH}N3BiasFieldCorrection $DIM ${OUTPUTNAME}repaired.nii  ${OUTPUTNAME}repaired.nii 4

exe=" ${ANTSPATH}ANTS $DIM -m  ${METRIC}${FIXED},${OUTPUTNAME}repaired.nii,${METRICPARAMS}  -t $TRANSFORMATION  -r $REGULARIZATION -o ${OUTPUTNAME}   -i $MAXITERATIONS   --use-Histogram-Matching "

 echo " $exe "

  $exe

  #below, some affine options
  #--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

    ${ANTSPATH}WarpImageMultiTransform $DIM  ${OUTPUTNAME}repaired.nii   ${OUTPUTNAME}deformed.nii ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt  -R ${FIXED}

if  [ ${#LABELIMAGE} -gt 3 ]
then
      ${ANTSPATH}WarpImageMultiTransform $DIM  $LABELIMAGE   ${OUTPUTNAME}labeled.nii -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii  -R ${MOVING}   --use-NN
fi

exit

if [ $DoANTSQC -eq 1 ] ;  then
#  measure image similarity
for SIM in 0 1 2 ; do
${ANTSPATH}MeasureImageSimilarity $DIM $SIM $FIXED ${OUTPUTNAME}deformed.nii
done
#  measure dice overlap and mds
${ANTSPATH}ThresholdImage $DIM $FIXED ${OUTPUTNAME}fixthresh.nii.gz Otsu 4
${ANTSPATH}ThresholdImage $DIM $MOVING ${OUTPUTNAME}movthresh.nii.gz Otsu 4
 ${ANTSPATH}WarpImageMultiTransform $DIM   ${OUTPUTNAME}movthresh.nii.gz  ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt  -R ${FIXED}   --use-NN
${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}dicestats.txt DiceAndMinDistSum  ${OUTPUTNAME}fixthresh.nii.gz   ${OUTPUTNAME}movthresh.nii.gz   ${OUTPUTNAME}mindistsum.nii.gz
#  labelstats for jacobian wrt segmenation
 # below, to compose
# ${ANTSPATH}ComposeMultiTransform $DIM   ${OUTPUTNAME}CompWarp.nii  -R $FIXED ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt
# ${ANTSPATH}CreateJacobianDeterminantImage $DIM ${OUTPUTNAME}CompWarp.nii ${OUTPUTNAME}jacobian.nii  0
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}movlabstat.txt LabelStats ${OUTPUTNAME}movthresh.nii.gz ${OUTPUTNAME}movthresh.nii.gz
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}jaclabstat.txt LabelStats ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}jacobian.nii
# we compare the output of these last two lines:
#  the Volume of the movlabstat computation vs. the mass of the jaclabstat
fi

exit
