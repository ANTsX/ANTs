#!/bin/bash

VERSION="0.0"

if ! command -v antsRegistration &> /dev/null
  then
    echo "we cant find the antsRegistration program -- does not seem to exist.  please \(re\)define \$PATH in your environment."
    exit
  fi

if ! command -v antsApplyTransforms &> /dev/null
  then
    echo "we cant find the antsApplyTransforms program -- does not seem to exist.  please \(re\)define \$PATH in your environment."
    exit
  fi

function Usage {
    cat <<USAGE

`basename $0` performs registration between a scalar image and a T1 image:

Usage:

`basename $0` -d imageDimension
              -r anatomicalT1image (brain or whole-head, depending on modality) to align to
              -R anatomicalReference image to warp to (often higher resolution than anatomicalT1image)
              -i scalarImageToMatch  (such as avgerage bold, averge dwi, etc.)
              -x anatomicalT1brainmask (should mask out regions that do not appear in scalarImageToMatch)
              -t transformType (0=rigid, 1=affine, 2=rigid+small_def, 3=affine+small_def)
              -w prefix of T1 to template transform
              -T template space


              <OPTARGS>
              -o outputPrefix
              -l labels in template space
              -a auxiliary scalar image/s to warp to template
              -b auxiliary dt image to warp to template

Example:

  bash $0 -d 3 -i pcasl_control.nii.gz -r t1.nii.gz -x t1_mask.nii.gz -a cbf.nii.gz -l template_aal.nii.gz -w t12template_ -t 2 -o output

minimal paramters that must be passed include:

-d  -i -r -x -w -o

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using antsIntermodalityIntrasubject with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${BRAIN}
      t1 subject brain        = ${ANATOMICAL_BRAIN}
      t1 subject brain mask   = ${TEMPLATE_MASK}
      output prefix           = ${OUTPUT_PREFIX}
      template labels         = ${TEMPLATE_LABELS}
      auxiliary images        = ${AUX_IMAGES[@]}
      diffusion tensor image  = ${DTI}

    ANTs parameters:
      metric                  = ${ANTS_METRIC}[fixedImage,movingImage,${ANTS_METRIC_PARAMS}]
      regularization          = ${ANTS_REGULARIZATION}
      transformations         = ${ANTS_TRANS1} ${ANTS_TRANS2}
      max iterations          = ${ANTS_MAX_ITERATIONS}

PARAMETERS
}

# Echos a command to both stdout and stderr, then runs it
function logCmd() {
  cmd="$@"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  ( "$@" )
  echo "END   <<<<<<<<<<<<<<<<<<<<"
  echo
  echo
}

################################################################################
#
# Main routine
#
################################################################################

HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"

KEEP_TMP_IMAGES=0

DIMENSION=3

BRAIN=""
AUX_IMAGES=()
TEMPLATE_TRANSFORM=""
ANATOMICAL_BRAIN=""
ANATOMICAL_SPACE="";
TEMPLATE_MASK=""
TEMPLATE_LABELS=""
TRANSFORM_TYPE="0"
DTI=""

TEMPLATE=""

################################################################################
#
# Programs and their parameters
#
################################################################################

ANTS=antsRegistration
WARP=antsApplyTransforms

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:r:x:d:h:i:l:o:t:w:R:T:" OPT
    do
      case $OPT in
          a) # auxiliary scalar images
              AUX_IMAGES[${#AUX_IMAGES[@]}]=$OPTARG
              ;;
          b) #brain extraction registration mask
              DTI=$OPTARG
              ;;
          d) #dimensions
              DIMENSION=$OPTARG
              if [[ ${DIMENSION} -gt 3 || ${DIMENSION} -lt 2 ]];
              then
                  echo " Error:  ImageDimension must be 2 or 3 "
                  exit 1
              fi
              ;;
          r) #brain extraction anatomical image
              ANATOMICAL_BRAIN=$OPTARG
              ;;
          R) # anatomical warp space
              ANATOMICAL_SPACE=$OPTARG
              ;;
          T) # template
              TEMPLATE=$OPTARG
              ;;
          x) #brain extraction registration mask
              TEMPLATE_MASK=$OPTARG
              ;;
          l) #brain extraction registration mask
              TEMPLATE_LABELS=$OPTARG
              ;;
          h) #help
              Usage >&2
              exit 0
              ;;
          i) #max_iterations
              BRAIN=$OPTARG
              ;;
          o) #output prefix
              OUTPUT_PREFIX=$OPTARG
              ;;
          w) #template registration image
              TEMPLATE_TRANSFORM=$OPTARG
              ;;
          t) #atropos prior weight
              TRANSFORM_TYPE=$OPTARG
              ;;
          *) # getopts issues an error message
              echo "ERROR:  unrecognized option -$OPT $OPTARG"
              exit 1
              ;;
      esac
  done
fi


if [[  ${#TEMPLATE} -lt 3 ]] ; then
  TEMPLATE=$TEMPLATE_LABELS
fi
################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#  3. See if $REGISTRATION_TEMPLATE is the same as $BRAIN_TEMPLATE
#
################################################################################

# requires DIMENSION BRAIN ANATOMICAL_BRAIN TEMPLATE_MASK TEMPLATE_TRANSFORM OUTPUT_PREFIX
if [[ ${#DIMENSION} -lt 1 ]];
  then
    echo "Please specify -d DIMENSION "
    exit 1
  fi
if [[ ! -f ${BRAIN} ]];
  then
    echo "The scalar brain:"
    echo "   $BRAIN"
    exit 1
  fi
if [[ ! -f ${ANATOMICAL_BRAIN} ]];
  then
    echo "The extraction template doesn't exist:"
    echo "   $ANATOMICAL_BRAIN"
    exit 1
  fi
if [[ ! -f ${TEMPLATE_MASK} ]];
  then
    echo "The brain extraction prior doesn't exist:"
    echo "   $TEMPLATE_MASK"
    exit 1
  fi

if [[ ${#ANATOMICAL_SPACE} -lt 1 ]];
  then
    echo setting ANATOMICAL_SPACE = ANATOMICAL_BRAIN
    ANATOMICAL_SPACE=$ANATOMICAL_BRAIN
  fi

if [[ ${OUTPUT_PREFIX} == */ ]];
  then
    OUTPUT_DIR=${OUTPUT_PREFIX%/}
  else
    OUTPUT_DIR=$(dirname $OUTPUT_PREFIX)
  fi

if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi


ANTS_MAX_ITERATIONS="50x50x0"
ANTS_TRANSFORMATION="SyN[ 0.1,3,0 ]"
ANTS_LINEAR_METRIC="MI"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[ 1000x500x250x0,1e-7,5 ]"
ANTS_METRIC="mattes"
ANTS_METRIC_PARAMS="1,32"
ANTS_TRANS1=""
ANTS_TRANS2=""

if [ ${TRANSFORM_TYPE} -eq 0 ];
then
    ANTS_TRANS1="Rigid[ 0.1 ]"
elif [ ${TRANSFORM_TYPE} -eq 1 ];
then
   ANTS_TRANS1="Affine[ 0.1 ]"
elif [ ${TRANSFORM_TYPE} -eq 2 ];
then
    ANTS_TRANS1="Rigid[ 0.1 ]"
    ANTS_TRANS2="SyN[ 0.1,3,0 ]"
elif [ ${TRANSFORM_TYPE} -eq 3 ];
then
   ANTS_TRANS1="Affine[ 0.1 ]"
   ANTS_TRANS2="SyN[ 0.1,3,0 ]"
fi

echoParameters >&2

echo "---------------------  Running `basename $0` on $HOSTNAME  ---------------------"

time_start=`date +%s`

################################################################################
#
# Output images
#
################################################################################

#BRAIN_EXTRACTION_MASK=${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
#BRAIN_SEGMENTATION=${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}
#CORTICAL_THICKNESS_IMAGE=${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}

################################################################################
#
# Intermodality matching
#
################################################################################

stage1="-m ${ANTS_LINEAR_METRIC}[ ${ANATOMICAL_BRAIN},${BRAIN},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t ${ANTS_TRANS1} -f 8x4x2x1 -s 4x2x1x0 -u 1"

stage2=""
if [[ ${TRANSFORM_TYPE} -gt 1 ]] ; then
    stage2="-m ${ANTS_METRIC}[ ${ANATOMICAL_BRAIN},${BRAIN},${ANTS_METRIC_PARAMS} ] -c [ ${ANTS_MAX_ITERATIONS},1e-7,5 ] -t ${ANTS_TRANS2} -f 4x2x1 -s 2x1x0mm -u 1"
fi
globalparams=" -z 1 --winsorize-image-intensities [ 0.005, 0.995 ] "

cmd="antsRegistration -d $DIMENSION $stage1 $stage2 $globalparams -o ${OUTPUT_PREFIX}"
echo $cmd
if [[ ! -s ${OUTPUT_PREFIX}0GenericAffine.mat ]] ; then
  $cmd
else
  echo we already have ${OUTPUT_PREFIX}0GenericAffine.mat
fi

# warp input image to t1
warp=""
iwarp=""
if [[ -s  ${OUTPUT_PREFIX}1Warp.nii.gz ]];
then
  warp="-t ${OUTPUT_PREFIX}1Warp.nii.gz"
  iwarp="-t ${OUTPUT_PREFIX}1InverseWarp.nii.gz"
fi
antsApplyTransforms -d $DIMENSION -i $BRAIN -o ${OUTPUT_PREFIX}anatomical.nii.gz -r $ANATOMICAL_SPACE $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear
if  [[ ! -s ${TEMPLATE_TRANSFORM}1Warp.nii.gz ]] ; then
  echo ${TEMPLATE_TRANSFORM}1Warp.nii.gz does not exist - please specify in order to proceed to steps that map to the template
  exit
fi
cmd="antsApplyTransforms -d $DIMENSION -i $BRAIN -o ${OUTPUT_PREFIX}template.nii.gz -r ${TEMPLATE} -t ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}0GenericAffine.mat $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear"
$cmd

echo "AUX IMAGES"
# warp auxiliary images to t1
for (( i = 0; i < ${#AUX_IMAGES[@]}; i++ ))
  do
    # FIXME - how to name these reasonably
    AUXO=`basename ${AUX_IMAGES[$i]} .nii.gz`
    #AUXO=${OUTPUT_PREFIX}_aux_${i}_warped.nii.gz
    antsApplyTransforms -d $DIMENSION -i ${AUX_IMAGES[$i]} -r $ANATOMICAL_SPACE $warp -n Linear -o ${OUTPUT_DIR}/${AUXO}_anatomical.nii.gz $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat

    antsApplyTransforms -d $DIMENSION -i ${AUX_IMAGES[$i]} -r ${TEMPLATE} -n Linear -o ${OUTPUT_DIR}/${AUXO}_template.nii.gz -t ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}0GenericAffine.mat $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat

done

echo "DTI"
# warp DT image to t1
if [[ -f $DTI ]];
then
    antsApplyTransforms -d $DIMENSION -i ${DTI} -r $ANATOMICAL_SPACE $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear -o ${OUTPUT_PREFIX}dt_anatomical.nii.gz -e 2

   antsApplyTransforms -d $DIMENSION -e 2 -i ${DTI} -r ${TEMPLATE} -n Linear -o ${OUTPUT_PREFIX}dt_template.nii.gz -t ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}0GenericAffine.mat $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat

   ImageMath 3 ${OUTPUT_PREFIX}fa_template.nii.gz TensorFA ${OUTPUT_PREFIX}dt_template.nii.gz
   ImageMath 3 ${OUTPUT_PREFIX}md_template.nii.gz TensorMeanDiffusion ${OUTPUT_PREFIX}dt_template.nii.gz

fi

# warp brainmask from anatomy to subject
if [[ -f $TEMPLATE_MASK ]]; then
    antsApplyTransforms -d $DIMENSION -i $TEMPLATE_MASK -o ${OUTPUT_PREFIX}brainmask.nii.gz \
	-r $BRAIN  -n NearestNeighbor \
	-t [ ${OUTPUT_PREFIX}0GenericAffine.mat, 1 ]       \
	$iwarp
fi

# warp Labels from template to subject (if passed)
if [[ -f $TEMPLATE_LABELS ]]; then
    antsApplyTransforms -d $DIMENSION -i $TEMPLATE_LABELS -o ${OUTPUT_PREFIX}labels.nii.gz \
	-r $BRAIN  -n NearestNeighbor \
	-t [ ${OUTPUT_PREFIX}0GenericAffine.mat, 1 ]       \
	$iwarp           \
	-t   ${TEMPLATE_TRANSFORM}1Warp.nii.gz \
	-t   ${TEMPLATE_TRANSFORM}0GenericAffine.mat
fi

################################################################################
#
# End of main routine
#
################################################################################

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with ANTs intermodality intrasubject processing pipeline"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
