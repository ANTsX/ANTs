#!/bin/bash

VERSION="0.0"

if [[ ! -s ${ANTSPATH}/antsRegistration ]]; then
  echo we cant find the antsRegistration program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
if [[ ! -s ${ANTSPATH}/antsApplyTransforms ]]; then
  echo we cant find the antsApplyTransforms program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

function Usage {
    cat <<USAGE

`basename $0` performs restration between a scalar image and a T1 image:

Usage:

`basename $0` -d imageDimension
              -r anatomicalT1image (brain or whole-head, depending on modality)
              -i scalarImageToMatch
              -x anatomicalT1brainmask
              -t transformType (0=rigid, 1=affine, 2=rigid+small_def, 3=affine+small_def)
              -w prefix of T1 to template transform

           
              <OPTARGS>
              -o outputPrefix
              -l labels in template space
              -a auxiliary scalar image/s to warp to template
              -b auxiliary dt image to warp to template

Example:

  bash $0 -d 3 -i pcasl_control.nii.gz -r t1.nii.gz -x t1_mask.nii.gz -a cbf.nii.gz -l template_aal.nii.gz -w t12template_ -t 2 -o output

Required arguments:

We use *intensity* to denote the original anatomical image of the brain.

We use *probability* to denote a probability image with values in range 0 to 1.

We use *label* to denote a label image with values in range 0 to N.

     -d:  Image dimension                       2 or 3 (for 2- or 3-dimensional image)
     -i:  Anatomical image                      Structural *intensity* scalar image such as avgerage bold, averge dwi, etc. 
     -r:  T1 Anatomical image                   Structural *intensity* image, typically T1.  
     -x:  Brain extraction probability mask     Brain *probability* mask created using e.g. LPBA40 labels which
                                                have brain masks defined, and warped to anatomical template and
                                                averaged resulting in a probability image.
     -w:  T1 to template transform prefix       Prefix for transform files that map T1 (-r) to the template space
                                                -p labelsPriors%02d.nii.gz.
     -o:  Output prefix                         The following images are created:
                                                  * ${OUTPUT_PREFIX}N4Corrected.${OUTPUT_SUFFIX}


Optional arguments:

     -t:  TranformType                          0=rigid, 1=affine, 2=rigid+small_def, 3=affine+small_def [default=1]
     -l:  Brain segmentation priors             Anatomical *label* image in template space to be mapped to modality space
     -a:  Auxiliary scalar files                Additional scalar files to warp to template space
     -b:  DT image                              DTI to warp to template space
      

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using antsIntramodalityInterSubject with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${BRAIN}
      t1 subject brain        = ${ANATOMICAL_BRAIN}
      t1 subject brain mask   = ${TEMPLATE_BRAIN}
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
  cmd="$*"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  $cmd
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
TEMPLATE_MASK=""
TEMPLATE_LABELS=""
TRANSFORM_TYPE="0"
DTI=""


################################################################################
#
# Programs and their parameters
#
################################################################################

ANTS=${ANTSPATH}antsRegistration
WARP=${ANTSPATH}antsApplyTransforms

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:r:x:d:h:i:l:o:t:w:" OPT
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

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#  3. See if $REGISTRATION_TEMPLATE is the same as $BRAIN_TEMPLATE
#
################################################################################



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
if [[ ! -f ${BRAIN} ]];
  then
    echo "The scalar brain:"
    echo "   $BRAIN"
    exit 1
  fi

OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi



ANTS_MAX_ITERATIONS="50x50x0"
ANTS_TRANSFORMATION="SyN[0.1,3,0]"
ANTS_LINEAR_METRIC="MI"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[1000x500x250x0,1e-7,5]"
ANTS_METRIC="mattes"
ANTS_METRIC_PARAMS="1,32"
ANTS_TRANS1=""
ANTS_TRANS2=""

if [ ${TRANSFORM_TYPE} -eq 0 ];
then
    ANTS_TRANS1="Rigid[0.1]"
elif [ ${TRANSFORM_TYPE} -eq 1 ];
then
   ANTS_TRANS1="Affine[0.1]"
elif [ ${TRANSFORM_TYPE} -eq 2 ];
then
    ANTS_TRANS1="Rigid[0.1]"
    ANTS_TRANS2="SyN[0.1,3,0]"
elif [ ${TRANSFORM_TYPE} -eq 3 ];
then
   ANTS_TRANS1="Affine[0.1]"
   ANTS_TRANS2="SyN[0.1,3,0]"
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

stage1="-m ${ANTS_LINEAR_METRIC}[${ANATOMICAL_BRAIN},${BRAIN},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t ${ANTS_TRANS1} -f 8x4x2x1 -s 4x2x1x0 -u 1"

stage2=""
if [[ ${TRANSFORM_TYPE} -gt 1 ]] ; then 
    stage2="-m ${ANTS_METRIC}[${ANATOMICAL_BRAIN},${BRAIN},${ANTS_METRIC_PARAMS}] -c [${ANTS_MAX_ITERATIONS},1e-7,5] -t ${ANTS_TRANS2} -f 4x2x1 -s 2x1x0mm -u 1"
fi
globalparams=" -z 1 --winsorize-image-intensities [0.005, 0.995] "

cmd="${ANTSPATH}/antsRegistration -d $DIMENSION $stage1 $stage2 $globalparams -o ${OUTPUT_PREFIX}"
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
${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i $BRAIN -o ${OUTPUT_PREFIX}anatomical.nii.gz -r $ANATOMICAL_BRAIN $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear
${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i $BRAIN -o ${OUTPUT_PREFIX}template.nii.gz -r ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}0GenericAffine.mat $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear

echo "AUX IMAGES"
# warp auxiliary images to t1
for (( i = 0; i < ${#AUX_IMAGES[@]}; i++ ))
  do
    # FIXME - how to name these reasonably
    AUXO=`basename ${AUX_IMAGES[$i]} .nii.gz`
    #AUXO=${OUTPUT_PREFIX}_aux_${i}_warped.nii.gz
    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i ${AUX_IMAGES[$i]} -r $ANATOMICAL_BRAIN $warp -n Linear -o ${OUTPUT_DIR}/${AUXO}_anatomical.nii.gz $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat

    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i ${AUX_IMAGES[$i]} -r ${TEMPLATE_TRANSFORM}1Warp.nii.gz $warp -n Linear -o ${OUTPUT_DIR}/${AUXO}_template.nii.gz -t ${TEMPLATE_TRANSFORM}1Warp.nii.gz -t ${TEMPLATE_TRANSFORM}0GenericAffine.mat $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat 

done

echo "DTI"
# warp DT image to t1
if [[ -f $DTI ]]; 
then
    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i ${DTI} -r $ANATOMICAL_BRAIN $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear -o ${OUTPUT_PREFIX}dt_anatomical.nii.gz -e 2
#    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i ${DTI} -r ${TEMPLATE_MASK} $warp -t ${OUTPUT_PREFIX}0GenericAffine.mat -n Linear -o ${OUTPUT_PREFIX}dt_anatomical.nii.gz -e 2
   # FIXME - reorientation
fi

# warp brainmask from anatomy to subject
if [[ -f $TEMPLATE_MASK ]]; then
    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i $TEMPLATE_MASK -o ${OUTPUT_PREFIX}brainmask.nii.gz \
	-r $BRAIN  -n NearestNeighbor \
	-t [ ${OUTPUT_PREFIX}0GenericAffine.mat, 1]       \
	$iwarp 
fi

# warp Labels from template to subject (if passed)
if [[ -f $TEMPLATE_LABELS ]]; then
    ${ANTSPATH}/antsApplyTransforms -d $DIMENSION -i $TEMPLATE_LABELS -o ${OUTPUT_PREFIX}labels.nii.gz \
	-r $BRAIN  -n NearestNeighbor \
	-t [ ${OUTPUT_PREFIX}0GenericAffine.mat, 1]       \
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
