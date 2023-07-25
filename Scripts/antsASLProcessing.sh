#!/bin/bash

if ! command -v antsRegistration &> /dev/null
then
  echo "Cannot find antsRegistration.  Please \(re\)define \$PATH in your environment."
fi
if ! command -v antsApplyTransforms &> /dev/null
then
  echo "Cannot find antsApplyTransforms.  Please \(re\)define \$PATH in your environment."
fi
if ! command -v antsIntermodalityIntrasubject.sh &> /dev/null
then
  echo "Cannot find antsIntermodalityIntrasubject.sh script.  Please \(re\)define \$PATH in your environemnt."
fi

function Usage {
  cat <<USAGE
`basename $0` performs ASL processing based on ANTs tools.  Preprocessing of T1 images using antsCorticalThickness.sh is assumed.  The following steps are performed:
  1) Calculation of average pCASL image.
  2) Skull stripping of average pCASL image.
  3) Registration of average pCASL image to T1 image.
  4) (Robust) calculation of mean CBF.
  5) Warping tissue priors and labels to ASL space.
  6) Warping mean CBF image to template space for VBM analysis.

Usage:
Required arguments:
`basename $0` -a anatomical image (skull stripped)
              -p brain segmentation priors (C-style, e.g. priors%d.nii.gz)
              -g hard brain segmentation
              -x t1 brain mask
              -s raw pCASL image
              -e brain template
	      -l template labels
              -t skull-stripped t1 to template transform prefix
              -o output prefix
Optional arguments:
              -b blood T1 value (defaults to 0.67 s^-1)
	      -r robustness parameter (defaults to 0.95)
	      -h print help and exit
	      -n number of bootstrap samples (defaults to 20)
	      -c percent to sample per bootstrap run (defaults to 70)
	      -k keep tmp files, including warps (defaults to false--takes lots of space to save)
	      -i use inverse warps.  Warps are assumed to go in the direction of subject to template.
	         If you are instead using template to subject warps (e.g. the brain segmentation prior warps from
	         antsCorticalThickness.sh), use -i (binary switch--no arguments).
	      -f bootstrap with replacement?  takes arguments "false" or "true"; defaults to false.
USAGE
  exit 1
}

echoParameters() {
  cat <<PARAMETERS
  Using `basename $0` with the following parameters:

   anatomical image:      ${ANATOMICAL_IMAGE}
   brain mask:            $BRAINMASK
   segmentation probs:    $SEGMENTATION_PROB
   hard segmentation:     $SEGMENTATION
   pCASL image:           $PCASL
   template:              $TEMPLATE
   template labels:       $LABELS
   transform prefix:      $TRANSFORM_PREFIX
   output prefix:         $OUTNAME
   blood relaxation:      $BLOODT1 s^-1
   robustness:            $ROBUST
   num bootstraps:        $NBOOTSTRAP
   pct per bootstrap:     $PCTBOOTSTRAP
   keep tmp files:        $KEEP_TMP_FILES
   use inverse warps:     $USE_INVERSE_WARPS
   sample w/ replacement: $SAMPLE_WITH_REPLACEMENT
PARAMETERS
}

function logCmd() {
  cmd="$*"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  $cmd
  echo "END   <<<<<<<<<<<<<<<<<<<<"
  echo
  echo
}

KEEP_TMP_FILES=false
USE_INVERSE_WARPS=false
if [[ $# -lt 3 ]]
then
  Usage >&2
  exit 1
else
  while getopts "a:s:e:p:t:o:x:l:b:r:g:n:c:f:k:i:h" OPT
  do
    case $OPT in
      a) #anatomical t1 image
    ANATOMICAL_IMAGE=$OPTARG
    ;;
      s) # raw pCASL image
    PCASL=$OPTARG
    ;;
      e) # brain template
    TEMPLATE=$OPTARG
    ;;
      p) # segmentation probabilities
    SEGMENTATION_PROB=$OPTARG
    ;;
      g) # hard seg
    SEGMENTATION=$OPTARG
    ;;
      t) # transform prefix
    TRANSFORM_PREFIX=$OPTARG
    ;;
      o) # output prefix
    OUTNAME=$OPTARG
    ;;
      x) # mask
    BRAINMASK=$OPTARG
    ;;
      l) # labels
    LABELS=$OPTARG
    ;;
      b) # blood t1
    BLOODT1=$OPTARG
    ;;
      r) # robustness
    ROBUST=$OPTARG
    ;;
      n) # number of bootstrap runs
    NBOOTSTRAP=$OPTARG
    ;;
      c) # pct to sample per bootstrap
    PCTBOOTSTRAP=$OPTARG
    ;;
      f) # sample with replacement?
    SAMPLE_WITH_REPLACEMENT=$OPTARG
    ;;
      k) # keep tmp files
    KEEP_TMP_FILES=true
    ;;
      i) # use inverse warps
    USE_INVERSE_WARPS=true
    ;;
      h) # help
    Usage >&2
    exit 0
    ;;
      *)
    echo "ERROR: unrecognized option -$OPT $OPTARG"
    exit 1
    ;;
    esac
  done
fi

if [[ -z $BLOODT1 ]]
then
  BLOODT1=0.67
fi
if [[ -z $ROBUST ]]
then
  ROBUST=0.95
fi
if [[ -z $NBOOTSTRAP ]]
then
  NBOOTSTRAP=20
fi
if [[ -z $PCTBOOTSTRAP ]]
then
  PCTBOOTSTRAP=0.70
fi
if [[ -z $SAMPLE_WITH_REPLACEMENT ]]
then
  SAMPLE_WITH_REPLACEMENT=false
fi
echoParameters >&2

# parse prior syntax
FORMAT=${SEGMENTATION_PROB}
PREFORMAT=${FORMAT%%\%*}
POSTFORMAT=${FORMAT##*d}
FORMAT=${FORMAT#*\%}
FORMAT=${FORMAT%%d*}
REPCHARACTER=''
TOTAL_LENGTH=0
if [ ${#FORMAT} -eq 2 ]
then
  REPCHARACTER=${FORMAT:0:1}
  TOTAL_LENGTH=${FORMAT:1:1}
fi
MAXNUMBER=1000
PRIOR_IMAGE_FILENAMES=()
for (( i = 1; i < $MAXNUMBER; i++ ))
do
  NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#i} ))
  ROOT='';
  for(( j=0; j < $NUMBER_OF_REPS; j++ ))
  do
    ROOT=${ROOT}${REPCHARACTER}
  done
  FILENAME=${PREFORMAT}${ROOT}${i}${POSTFORMAT}
  if [[ -f $FILENAME ]];
  then
    PRIOR_IMAGE_FILENAMES=( ${PRIOR_IMAGE_FILENAMES[@]} $FILENAME )
  else
    break 1
  fi
done

echo "Warning: antsASLProcessing.sh is being deprecated.  Use antsASLProcessing.R instead."

# check for existence of all images
if [[ ! -f ${ANATOMICAL_IMAGE} ]]
then
  echo "ERROR: Anatomical image ${ANATOMICAL_IMAGE} does not exist."
  exit 1
fi
if [[ ! -f $BRAINMASK ]]
then
  echo "ERROR: Brain mask $BRAINMASK does not exist."
  exit 1
fi
if [[ ! -f $PCASL ]]
then
  echo "ERROR: pCASL image $PCASL does not exist."
  exit 1
fi
if [[ ! -f $TEMPLATE ]]
then
  echo "ERROR: template image $TEMPLATE does not exist."
  exit 1
fi
if [[ ! -f $LABELS ]]
then
  echo "ERROR: Template label image $LABELS does not exist."
fi
if [[ ! -f $SEGMENTATION ]]
then
  echo "ERROR: Segmentation image $SEGMENTATION does not exist."
  exit 1
fi
if [[ ! -f ${TRANSFORM_PREFIX}SubjectToTemplate1Warp.nii.gz ]]
then
  echo "ERROR: Warp ${TRANSFORM_PREFIX}SubjectToTemplate1Warp.nii.gz does not exist."
  exit 1
fi
if [[ ${#PRIOR_IMAGE_FILENAMES[@]} -lt 3 ]]
then
  echo "ERROR: Fewer than 3 prior images specified."
  echo "       Check that you defined prior file names correctly."
fi

time_start=`date +%s`


# Do processing.
if [[ ! -d `dirname $OUTNAME` ]]
then
  mkdir -p `dirname $OUTNAME`
fi
if [[ ! -s ${OUTNAME}/AveragePCASL.nii.gz ]] ; then
  logCmd antsMotionCorr -d 3 -a $PCASL -o ${OUTNAME}AveragePCASL.nii.gz
fi
logCmd ThresholdImage 3 ${OUTNAME}AveragePCASL.nii.gz ${OUTNAME}tmp.nii.gz 600 999999
logCmd ImageMath 3 ${OUTNAME}tmp.nii.gz ME ${OUTNAME}tmp.nii.gz 2
logCmd ImageMath 3 ${OUTNAME}tmp.nii.gz GetLargestComponent ${OUTNAME}tmp.nii.gz
logCmd ImageMath 3 ${OUTNAME}tmp.nii.gz MD ${OUTNAME}tmp.nii.gz 3
logCmd ImageMath 3 ${OUTNAME}pCASLBrain.nii.gz m ${OUTNAME}AveragePCASL.nii.gz ${OUTNAME}tmp.nii.gz

INTERSUBJECT_PARAMS=" -d 3 -i ${OUTNAME}pCASLBrain.nii.gz -r $ANATOMICAL_IMAGE -x $BRAINMASK -w ${TRANSFORM_PREFIX}SubjectToTemplate -t 2 -o $OUTNAME "
if [[ -n $LABELS ]]
then
  INTERSUBJECT_PARAMS=" ${INTERSUBJECT_PARAMS} -l $LABELS "
fi

logCmd antsIntermodalityIntrasubject.sh $INTERSUBJECT_PARAMS
logCmd N4BiasFieldCorrection -d 3 -i ${OUTNAME}AveragePCASL.nii.gz -o  ${OUTNAME}AveragePCASL.nii.gz -r 1 -s 4
logCmd N4BiasFieldCorrection -d 3 -i ${OUTNAME}AveragePCASL.nii.gz -o  ${OUTNAME}AveragePCASL.nii.gz -r 1 -s 2
logCmd N4BiasFieldCorrection -d 3 -i ${OUTNAME}AveragePCASL.nii.gz -o  ${OUTNAME}AveragePCASL.nii.gz -r 1 -s 1
logCmd ThresholdImage 3 ${OUTNAME}AveragePCASL.nii.gz ${OUTNAME}OtsuMask.nii.gz Otsu 4
logCmd ThresholdImage 3 ${OUTNAME}OtsuMask.nii.gz ${OUTNAME}OtsuMask.nii.gz 2 4
logCmd ImageMath 3 ${OUTNAME}OtsuMask.nii.gz ME ${OUTNAME}OtsuMask.nii.gz 1
logCmd ImageMath 3 ${OUTNAME}OtsuMask.nii.gz MD ${OUTNAME}OtsuMask.nii.gz 1
logCmd ThresholdImage 3 ${OUTNAME}brainmask.nii.gz ${OUTNAME}BrainThresh.nii.gz 1 999
logCmd MultiplyImages 3 ${OUTNAME}OtsuMask.nii.gz ${OUTNAME}BrainThresh.nii.gz  ${OUTNAME}OtsuMask.nii.gz

if [ ! -f ${OUTNAME}_kcbf.nii.gz ]; then
  logCmd antsNetworkAnalysis.R \
    -o $OUTNAME \
    --freq 0.01x0.1 \
    --mask ${OUTNAME}OtsuMask.nii.gz \
    --labels ${OUTNAME}labels.nii.gz \
    --fmri $PCASL \
    --modality ASLCBF \
    --bloodt1 $BLOODT1 \
    --robust $ROBUST \
    --nboot $NBOOTSTRAP \
    --pctboot $PCTBOOTSTRAP \
    --replace $SAMPLE_WITH_REPLACEMENT
fi

logCmd antsApplyTransforms -d 3 \
    -i ${OUTNAME}_kcbf.nii.gz \
    -r $TEMPLATE \
    -o ${OUTNAME}MeanCBFWarpedToTemplate.nii.gz \
    -n Linear \
    -t ${TRANSFORM_PREFIX}SubjectToTemplate0GenericAffine.mat \
    -t ${TRANSFORM_PREFIX}SubjectToTemplate1Warp.nii.gz \
    -t ${OUTNAME}1Warp.nii.gz \
    -t ${OUTNAME}0GenericAffine.mat

logCmd antsApplyTransforms -d 3 \
    -i $LABELS \
    -r ${OUTNAME}AveragePCASL.nii.gz \
    -o ${OUTNAME}LabelsWarpedToPCASL.nii.gz \
    -n MultiLabel \
    -t ${TRANSFORM_PREFIX}TemplateToSubject0Warp.nii.gz \
    -t ${TRANSFORM_PREFIX}TemplateToSubject1GenericAffine.mat \
    -t [ ${OUTNAME}0GenericAffine.mat,1 ] \
    -t ${OUTNAME}1InverseWarp.nii.gz

logCmd antsApplyTransforms -d 3 \
  -i ${OUTNAME}_kcbf.nii.gz \
  -r $ANATOMICAL_IMAGE \
  -o ${OUTNAME}MeanCBFWarpedToT1.nii.gz \
  -n Linear \
  -t ${OUTNAME}1Warp.nii.gz \
  -t ${OUTNAME}0GenericAffine.mat \

logCmd antsApplyTransforms -d 3 \
  -i $SEGMENTATION \
  -r ${OUTNAME}AveragePCASL.nii.gz \
  -o ${OUTNAME}SegmentationWarpedToPCASL.nii.gz \
  -n MultiLabel \
  -t [ ${OUTNAME}0GenericAffine.mat,1 ] \
  -t ${OUTNAME}1InverseWarp.nii.gz \

if ! $KEEP_TMP_FILES
then
  for FILE in 0GenericAffine.mat 1Warp.nii.gz 1InverseWarp.nii.gz anatomical.nii.gz template.nii.gz tmp.nii.gz
  do
    logCmd rm ${OUTNAME}${FILE}
  done
fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with ANTs ASL processing pipeline."
echo " Script executed in $time_elapsed seconds."
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
