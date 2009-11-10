# !/bin/sh

NUMPARAMS=$#


if [ $NUMPARAMS -lt 2  ]
then
echo " USAGE ::  "
echo "  sh   removeskull.txt imagefile MorphologyRadius-Int TemplateImage TemplateBrainProbabilityImage "
echo " N-Classes - usually 3 or 4 "
echo " MorphologyRadius = Usually6to8 "
echo " If-Available-Template and TemplateBrainMask "
echo " "
echo " if the masking is too aggressive, reduce the percentofmean or change the morphology radius "
echo " Be sure to set ANTSPATH "
exit
fi

IMG=$1
MORPH=$2
USEPRIOR=0
MINSIZE=50000

if [ $NUMPARAMS -gt 3  ]
then
TEIMG=$3 #PCEtemplate.nii
TEPROB=$4   #PCEtemplatebrainprob.nii
echo " USING PRIOR From Template data:  $4 and prior for brain $5 "
USEPRIOR=1

fi

PRE=${ANTSPATH}
#/home/avants/bin/ants/



# for gzipped images use line below
BASE=${IMG%.*.*}
echo "  ${#BASE} and ${#IMG} "
if [ ${#BASE} -eq ${#IMG} ]
then
BASE=${IMG%.*}
fi

OUTPUT=${BASE}brain.nii.gz
echo "output:  ${OUTPUT}  "

#${PRE}ImageMath 3 ${OUTPUT} ThresholdAtMean ${IMG} ${PERCMEAN}
${PRE}ThresholdImage 3  ${IMG} ${OUTPUT} Otsu 4

PADVAL=15
NPADVAL=-15
${PRE}ImageMath 3 ${OUTPUT} PadImage ${OUTPUT} $PADVAL

${PRE}ThresholdImage 3  ${OUTPUT} ${OUTPUT} 2 1.e9

${PRE}ImageMath 3 ${OUTPUT} ME ${OUTPUT} $MORPH

${PRE}ImageMath 3 ${OUTPUT} MD ${OUTPUT} $MORPH

${PRE}ImageMath 3 ${OUTPUT} GetLargestComponent ${OUTPUT} $MINSIZE

#exit
#VAL=2
#MORPHB=4

echo " Closing $MORPHB "

${PRE}ImageMath 3 ${OUTPUT} MD ${OUTPUT} $MORPH

${PRE}ImageMath 3 ${OUTPUT} ME ${OUTPUT} $MORPH

${PRE}ImageMath 3 ${OUTPUT} GetLargestComponent ${OUTPUT} $MINSIZE

#echo " EXIT"
#exit

# now smooth & convert to probability ?
SIG=1
${PRE}ImageMath 3 ${OUTPUT} G ${OUTPUT} $SIG


# now use atlas prior - if available

if [ ${USEPRIOR} = 1 ]
then
echo " ATLAS is:: $TEPROB $TEIMG "

${PRE}ImageMath 3 ${BASE}temptar.nii PadImage $IMG $PADVAL

exe=" ${PRE}ANTS 3 -m MI[${BASE}temptar.nii,${TEIMG},1,32] -o ${BASE}brain  -i 11x0x0  -t SyN[0.3] "
echo " $exe "
$exe

#${PRE}WarpImageMultiTransform 3 $TEPROB ${BASE}brainprob.nii  -i PCE${BASE}Affine.txt  PCE${BASE}InverseWarp.nii  -R $IMG

${PRE}WarpImageMultiTransform 3 $TEPROB ${BASE}brainprob.nii  ${BASE}brainWarp.nii   ${BASE}brainAffine.txt  -R $OUTPUT

#${PRE}ImageMath 3 ${OUTPUT} exp ${OUTPUT} 1
${PRE}ImageMath 3 ${BASE}brainprob.nii  G ${BASE}brainprob.nii $SIG
${PRE}ImageMath 3 ${OUTPUT} m ${OUTPUT} ${BASE}brainprob.nii
echo "  ${BASE}brainprob.nii  "

rm -f ${BASE}brain*Warp* ${BASE}brain*Affine* ${BASE}temptar.nii
fi


 ${PRE}ThresholdImage 3 ${OUTPUT} ${OUTPUT}  0.25  1.0

 ${PRE}ImageMath 3 ${OUTPUT} GetLargestComponent ${OUTPUT} $MINSIZE


 ${PRE}ImageMath 3 ${OUTPUT} MD ${OUTPUT} 10

 ${PRE}ImageMath 3 ${OUTPUT} ME ${OUTPUT} 10


${PRE}ImageMath 3 ${OUTPUT} PadImage ${OUTPUT} $NPADVAL
 ${PRE}ImageMath 3 ${OUTPUT} m ${OUTPUT} ${IMG}

# ${PRE}ImageMath 3 ${OUTPUT} FlattenImage ${OUTPUT} 0.9




if [ ${USEPRIOR} = 1 ]
then
${PRE}ImageMath 3 ${BASE}brainprob.nii PadImage ${BASE}brainprob.nii $NPADVAL
fi


#/home/avants/bin/startmricro ${OUTPUT}
