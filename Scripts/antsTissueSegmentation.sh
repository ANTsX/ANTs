#!/bin/bash

VERSION="0.0"

# Check dependencies

PROGRAM_DEPENDENCIES=( 'antsRegistration' 'antsApplyTransforms' 'N4BiasFieldCorrection' 'Atropos' )
SCRIPTS_DEPENDENCIES=( 'antsBrainExtraction.sh' 'antsAtroposN4.sh' )

for D in ${PROGRAM_DEPENDENCIES[@]};
  do
    if [[ ! -s ${ANTSPATH}/${D} ]];
      then
        echo "Error:  we can't find the $D program."
        echo "Perhaps you need to \(re\)define \$ANTSPATH in your environment."
        exit
      fi
  done

for D in ${SCRIPT_DEPENDENCIES[@]};
  do
    if [[ ! -s ${ANTSPATH}/${D} ]];
      then
        echo "We can't find the $D script."
        echo "Perhaps you need to \(re\)define \$ANTSPATH in your environment."
        exit
      fi
  done

function Usage {
    cat <<USAGE

`basename $0` : Run ANTs Atropos to run tissue segmentation based on input template/priors

  1. Brain extraction
  2. Brain n-tissue segmentation

Usage:

`basename $0` -d imageDimension
              -a anatomicalImage
              -e brainTemplate
              -m brainExtractionProbabilityMask
              -p brainSegmentationPriors
              <OPTARGS>
              -o outputPrefix

Example:

  bash $0 -d 3 -a t1.nii.gz -e brainWithSkullTemplate.nii.gz -m brainPrior.nii.gz -p segmentationPriors%d.nii.gz -o output

Required arguments:

We use *intensity* to denote the original anatomical image of the brain.

We use *probability* to denote a probability image with values in range 0 to 1.

We use *label* to denote a label image with values in range 0 to N.

     -d:  Image dimension                       2 or 3 (for 2- or 3-dimensional image)
     -a:  Anatomical image                      Structural *intensity* image, typically T1.  If more than one
                                                anatomical image is specified, subsequently specified
                                                images are used during the segmentation process.  However,
                                                only the first image is used in the registration of priors.
                                                Our suggestion would be to specify the T1 as the first image.
     -e:  Brain template                        Anatomical *intensity* template (possibly created using a population
                                                data set with buildtemplateparallel.sh in ANTs).  This template is
                                                *not* skull-stripped.
     -m:  Brain extraction probability mask     Brain *probability* mask created using e.g. LPBA40 labels which
                                                have brain masks defined, and warped to anatomical template and
                                                averaged resulting in a probability image.
     -p:  Brain segmentation priors             Tissue *probability* priors corresponding to the image specified
                                                with the -e option.  Specified using c-style formatting, e.g.
                                                -p labelsPriors%02d.nii.gz.  We assume that the first four priors
                                                are ordered as follows
                                                  1:  csf
                                                  2:  cortical gm
                                                  3:  wm
                                                  4:  deep gm
     -o:  Output prefix                         The following images are created:
                                                  * ${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentation*N4.${OUTPUT_SUFFIX} One for each anatomical input
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*1.${OUTPUT_SUFFIX}  CSF
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*2.${OUTPUT_SUFFIX}  GM
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*3.${OUTPUT_SUFFIX}  WM
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*4.${OUTPUT_SUFFIX}  DEEP GM
                                                  * ...
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*N.${OUTPUT_SUFFIX} where there are N priors
                                                  *                              Number formatting of posteriors matches that of the priors.

Optional arguments:

     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd
     -t:  template for t1 registration          Anatomical *intensity* template (assumed to be skull-stripped).  A common
                                                use case would be where this would be the same template as specified in the
                                                -e option which is not skull stripped.
                                                We perform the registration (fixed image = individual subject
                                                and moving image = template) to produce the files.
                                                The output from this step is
                                                  * ${OUTPUT_PREFIX}TemplateToSubject0GenericAffine.mat
                                                  * ${OUTPUT_PREFIX}TemplateToSubject1Warp.nii.gz
                                                  * ${OUTPUT_PREFIX}TemplateToSubject1InverseWarp.nii.gz
                                                  * ${OUTPUT_PREFIX}TemplateToSubjectLogJacobian.${OUTPUT_SUFFIX}
     -f:  extraction registration mask          Mask (defined in the template space) used during registration
                                                for brain extraction.
     -k:  keep temporary files                  Keep brain extraction/segmentation warps, etc (default = 0).
     -g:  denoise anatomical images             Denoise anatomical images (default = 0).
     -i:  max iterations for registration       ANTS registration max iterations (default = 100x100x70x20)
     -w:  Atropos prior segmentation weight     Atropos spatial prior *probability* weight for the segmentation (default = 0.25)
     -n:  number of segmentation iterations     N4 -> Atropos -> N4 iterations during segmentation (default = 3)
     -b:  posterior formulation                 Atropos posterior formulation and whether or not to use mixture model proportions.
                                                e.g 'Socrates[ 1 ]' (default) or 'Aristotle[ 1 ]'.  Choose the latter if you
                                                want use the distance priors (see also the -l option for label propagation
                                                control).
     -j:  use floating-point precision          Use floating point precision in registrations (default = 0)
     -u:  use random seeding                    Use random number generated from system clock in Atropos (default = 1)
     -v:  use b-spline smoothing                Use B-spline SyN for registrations and B-spline exponential mapping in DiReCT.
     -l:  label propagation                     Incorporate a distance prior one the posterior formulation.  Should be
                                                of the form 'label[ lambda,boundaryProbability ]' where label is a value
                                                of 1,2,3,... denoting label ID.  The label probability for anything
                                                outside the current label

                                                  = boundaryProbability * exp( -lambda * distanceFromBoundary )

                                                Intuitively, smaller lambda values will increase the spatial capture
                                                range of the distance prior.  To apply to all label values, simply omit
                                                specifying the label, i.e. '-l "[ lambda,boundaryProbability ]"'.
     -c                                         Add prior combination to combined gray and white matters.  For example,
                                                when calling KK for normal subjects, we combine the deep gray matter
                                                segmentation/posteriors with the white matter segmentation/posteriors.
                                                An additional example would be performing cortical thickness in the presence
                                                of white matter lesions.  We can accommodate this by specifying a lesion mask
                                                posterior as an additional posterior (suppose label '7'), and then combine
                                                this with white matter by specifying '-c "WM[ 7 ]"' or '-c "3[ 7 ]"'.
     -q:  Use quick registration parameters     If = 1, use antsRegistrationSyNQuick.sh as the basis for registration
                                                during brain extraction, brain segmentation, and (optional) normalization
                                                to a template.  Otherwise use antsRegistrationSyN.sh (default = 0).
     -x:                                        Number of iterations within Atropos (default 5).
     -y:                                        Which stage of ATS to run (default = 0, run all).  Tries to split in 2 hour chunks.
                                                Will produce OutputNameATSStageNComplete.txt for each completed stage.
                                                1: brain extraction
                                                2: template registration
                                                3: tissue segmentation
                                                4: template registration (improved, optional)
                                                5: qc, quality control
     -z:  Test / debug mode                     If > 0, runs a faster version of the script. Only for testing. Implies -u 0.
                                                Requires single thread computation for complete reproducibility.
USAGE
    exit 1
}

# Check outputs exist, runs at the end of the script
# List of outputs is taken from the usage
function checkOutputExists() {

  singleOutputs=( ${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX} )
  singleOutputs=( ${singleOutputs[@]} ${OUTPUT_PREFIX}BrainSegmentationTiledMosaic.png )

  if [[ -f ${REGISTRATION_TEMPLATE} ]];
    then
      singleOutputs=( ${singleOutputs[@]} ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}0GenericAffine.mat ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1Warp.nii.gz ${REGISTRATION_SUBJECT_OUTPUT_PREFIX}0Warp.nii.gz ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}LogJacobian.${OUTPUT_SUFFIX} )
    fi

  missingOutput=0

  for img in ${singleOutputs[@]};
    do
      if [[ ! -f $img ]];
        then
          echo "Missing output image $img"
          missingOutput=1
        fi
    done

  # Now check numbered output, numbers based on images
  for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
    do
      if [[ ! -f ${OUTPUT_PREFIX}BrainSegmentation${i}N4.${OUTPUT_SUFFIX} ]];
        then
          echo "Missing output image ${OUTPUT_PREFIX}BrainSegmentation${i}N4.${OUTPUT_SUFFIX}"
          missingOutput=1
        fi
    done

  # Segmentation output depends on the number of priors and the numbering format
  segNumWidth=${#GRAY_MATTER_LABEL_FORMAT}

  for (( j = 1; j <= ${NUMBER_OF_PRIOR_IMAGES}; j++ ));
    do
      num=$(printf "%0${segNumWidth}d" $j)

      if [[ ! -f ${OUTPUT_PREFIX}BrainSegmentationPosteriors${num}.${OUTPUT_SUFFIX} ]];
        then
          echo "Missing output image ${OUTPUT_PREFIX}BrainSegmentationPosteriors${num}.${OUTPUT_SUFFIX}"
          missingOutput=1
        fi
    done

  if [[ $missingOutput -gt 0 ]];
    then
      echo "Some of the output does not exist"
      return 1
    fi

  return 0
}

echoParameters() {
    cat <<PARAMETERS

    Using antsTissueSegmentation with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      brain template          = ${BRAIN_TEMPLATE}
      extraction prior        = ${EXTRACTION_PRIOR}
      extraction reg. mask    = ${EXTRACTION_REGISTRATION_MASK}
      segmentation prior      = ${SEGMENTATION_PRIOR}
      output prefix           = ${OUTPUT_PREFIX}
      output image suffix     = ${OUTPUT_SUFFIX}
      registration template   = ${REGISTRATION_TEMPLATE}

    ANTs parameters:
      metric                  = ${ANTS_METRIC}[ fixedImage,movingImage,${ANTS_METRIC_PARAMS} ]
      transformation          = ${ANTS_TRANSFORMATION}
      max iterations          = ${ANTS_MAX_ITERATIONS}

    Other parameters:
      run quick               = ${RUN_QUICK}
      debug mode              = ${DEBUG_MODE}
      denoise images          = ${DENOISE_ANATOMICAL_IMAGES}
      float precision         = ${USE_FLOAT_PRECISION}
      use random seeding      = ${USE_RANDOM_SEEDING}
      prior combinations      = ${PRIOR_COMBINATIONS[@]}

PARAMETERS
}

# Echos a command to stdout, then runs it
# Will immediately exit on error unless you set debug flag here
DEBUG_MODE=0

ATS_STAGE=0 # run all stages

function logCmd() {
  cmd="$@"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  # Preserve quoted parameters by running "$@" instead of $cmd
  ( "$@" )

  cmdExit=$?

  if [[ $cmdExit -gt 0 ]];
    then
      echo "ERROR: command exited with nonzero status $cmdExit"
      echo "Command: $cmd"
      echo
      if [[ ! $DEBUG_MODE -gt 0 ]];
        then
          exit 1
        fi
    fi

  echo "END   <<<<<<<<<<<<<<<<<<<<"
  echo
  echo

  return $cmdExit
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
OUTPUT_SUFFIX="nii"

KEEP_TMP_IMAGES=0

DIMENSION=3

ANATOMICAL_IMAGES=()
REGISTRATION_TEMPLATE=""
DO_REGISTRATION_TO_TEMPLATE=0

USE_RANDOM_SEEDING=1
RUN_QUICK=0

DENOISE_ANATOMICAL_IMAGES=0

BRAIN_TEMPLATE=""
EXTRACTION_PRIOR=""
EXTRACTION_REGISTRATION_MASK=""
SEGMENTATION_PRIOR=""

CSF_MATTER_LABEL=1
GRAY_MATTER_LABEL=2
WHITE_MATTER_LABEL=3
DEEP_GRAY_MATTER_LABEL=4

ATROPOS_SEGMENTATION_PRIOR_WEIGHT=0.25

################################################################################
#
# Programs and their parameters
#
################################################################################

ANTS=${ANTSPATH}/antsRegistration
ANTS_MAX_ITERATIONS="100x100x70x20"
ANTS_TRANSFORMATION="SyN[ 0.1,3,0 ]"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[ 1000x500x250x100,1e-8,10 ]"
ANTS_METRIC="CC"
ANTS_METRIC_PARAMS="1,4"

WARP=${ANTSPATH}/antsApplyTransforms

N4=${ANTSPATH}/N4BiasFieldCorrection
N4_CONVERGENCE_1="[ 50x50x50x50,0.0000001 ]"
N4_CONVERGENCE_2="[ 50x50x50x50,0.0000001 ]"
N4_SHRINK_FACTOR_1=4
N4_SHRINK_FACTOR_2=2
N4_BSPLINE_PARAMS="[ 200 ]"

ATROPOS=${ANTSPATH}/Atropos

ATROPOS_SEGMENTATION_INITIALIZATION="PriorProbabilityImages"
ATROPOS_SEGMENTATION_LIKELIHOOD="Gaussian"
ATROPOS_SEGMENTATION_CONVERGENCE="[ 5,0.0 ]"
ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION="Socrates[ 1 ]"
ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=3
ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5 # to be backward compatible but i like 25
ATROPOS_SEGMENTATION_LABEL_PROPAGATION=()

PRIOR_COMBINATIONS=( 'WM[ 4 ]' )

USE_FLOAT_PRECISION=0
USE_BSPLINE_SMOOTHING=0

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:" OPT
    do
      case $OPT in
          a) #anatomical t1 image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
       ;;
          b) # posterior formulation
       ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION=$OPTARG
       ;;
          c) # prior combinations
       PRIOR_COMBINATIONS[${#PRIOR_COMBINATIONS[@]}]=$OPTARG
       ;;
          d) #dimensions
       DIMENSION=$OPTARG
       if [[ ${DIMENSION} -gt 3 || ${DIMENSION} -lt 2 ]];
         then
           echo " Error:  ImageDimension must be 2 or 3 "
           exit 1
         fi
       ;;
          e) #brain extraction anatomical image
       BRAIN_TEMPLATE=$OPTARG
       ;;
          f) #brain extraction registration mask
       EXTRACTION_REGISTRATION_MASK=$OPTARG
       ;;
          g) # denoise anatomical images
       DENOISE_ANATOMICAL_IMAGES=$OPTARG
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          i) #max_iterations
       ANTS_MAX_ITERATIONS=$OPTARG
       ;;
          j) #use floating point precision
       USE_FLOAT_PRECISION=$OPTARG
       ;;
          k) #keep tmp images
       KEEP_TMP_IMAGES=$OPTARG
       ;;
          l)
       ATROPOS_SEGMENTATION_LABEL_PROPAGATION[${#ATROPOS_SEGMENTATION_LABEL_PROPAGATION[@]}]=$OPTARG
       ;;
          m) #brain extraction prior probability mask
       EXTRACTION_PRIOR=$OPTARG
       ;;
          n) #atropos segmentation iterations
       ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=$OPTARG
       ;;
          x) #atropos segmentation internal iterations
       ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=$OPTARG
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          p) #brain segmentation label prior image
       SEGMENTATION_PRIOR=$OPTARG
       ;;
          q) # run quick
       RUN_QUICK=$OPTARG
       ;;
          s) #output suffix
       OUTPUT_SUFFIX=$OPTARG
       ;;
          t) #template registration image
       REGISTRATION_TEMPLATE=$OPTARG
       DO_REGISTRATION_TO_TEMPLATE=1
       ;;
          u) #use random seeding
       USE_RANDOM_SEEDING=$OPTARG
       ;;
          v) #use b-spline smoothing in registration and direct
       USE_BSPLINE_SMOOTHING=$OPTARG
       ;;
          w) #atropos prior weight
       ATROPOS_SEGMENTATION_PRIOR_WEIGHT=$OPTARG
       ;;
          y) #which stage
       ATS_STAGE=$OPTARG
       ;;
          z) #debug mode
       DEBUG_MODE=$OPTARG
       ;;
          *) # getopts issues an error message
       echo "ERROR:  unrecognized option -$OPT $OPTARG"
       exit 1
       ;;
      esac
  done
fi

if [[ $USE_BSPLINE_SMOOTHING -ne 0 ]];
  then
    ANTS_TRANSFORMATION="BSplineSyN[ 0.1,26,0,3 ]"
  fi

if [[ $DEBUG_MODE -gt 0 ]];
  then

    echo "    WARNING - Running in test / debug mode. Results will be suboptimal "

    OUTPUT_PREFIX="${OUTPUT_PREFIX}testMode_"

    # Speed up by doing fewer its. Careful about changing this because
    # certain things are hard coded elsewhere, eg number of levels

    ANTS_MAX_ITERATIONS="40x40x20x0"
    ANTS_LINEAR_CONVERGENCE="[ 100x100x50x0,1e-8,10 ]"
    ANTS_METRIC_PARAMS="1,2"

    # I think this is the number of times we run the whole N4 / Atropos thing, at the cost of about 10 minutes a time
    ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=1
    ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5 # internal to atropos

    # Fix random seed to replicate exact results on each run
    USE_RANDOM_SEEDING=0

  fi

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#
################################################################################

for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do
  if [[ ! -f ${ANATOMICAL_IMAGES[$i]} ]];
    then
      echo "The specified image \"${ANATOMICAL_IMAGES[$i]}\" does not exist."
      exit 1
    fi
  done

FORMAT=${SEGMENTATION_PRIOR}
PREFORMAT=${FORMAT%%\%*}
POSTFORMAT=${FORMAT##*d}
FORMAT=${FORMAT#*\%}
FORMAT=${FORMAT%%d*}

REPCHARACTER=''
TOTAL_LENGTH=0
if [[ ${#FORMAT} -eq 2 ]]
  then
    REPCHARACTER=${FORMAT:0:1}
    TOTAL_LENGTH=${FORMAT:1:1}
  fi

# MAXNUMBER=$(( 10 ** $TOTAL_LENGTH ))
MAXNUMBER=1000

PRIOR_IMAGE_FILENAMES=()
WARPED_PRIOR_IMAGE_FILENAMES=()
BRAIN_SEGMENTATION_OUTPUT=${OUTPUT_PREFIX}BrainSegmentation
SEGMENTATION_WARP_OUTPUT_PREFIX=${BRAIN_SEGMENTATION_OUTPUT}Prior
SEGMENTATION_PRIOR_WARPED=${SEGMENTATION_WARP_OUTPUT_PREFIX}Warped
for (( i = 1; i < $MAXNUMBER; i++ ))
  do
    NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#i} ))
    ROOT='';
    for(( j=0; j < $NUMBER_OF_REPS; j++ ))
      do
        ROOT=${ROOT}${REPCHARACTER}
      done
    FILENAME=${PREFORMAT}${ROOT}${i}${POSTFORMAT}
    WARPED_FILENAME=${SEGMENTATION_PRIOR_WARPED}${ROOT}${i}.${OUTPUT_SUFFIX}
    if [[ -f $FILENAME ]];
      then
        PRIOR_IMAGE_FILENAMES=( ${PRIOR_IMAGE_FILENAMES[@]} $FILENAME )
        WARPED_PRIOR_IMAGE_FILENAMES=( ${WARPED_PRIOR_IMAGE_FILENAMES[@]} $WARPED_FILENAME )
      else
        break 1
      fi
  done

NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#DEEP_GRAY_MATTER_LABEL} ))
ROOT='';
for(( j=0; j < $NUMBER_OF_REPS; j++ ))
  do
    ROOT=${ROOT}${REPCHARACTER}
  done
DEEP_GRAY_MATTER_LABEL_FORMAT=${ROOT}${DEEP_GRAY_MATTER_LABEL}

NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#WHITE_MATTER_LABEL} ))
ROOT='';
for(( j=0; j < $NUMBER_OF_REPS; j++ ))
  do
    ROOT=${ROOT}${REPCHARACTER}
  done
WHITE_MATTER_LABEL_FORMAT=${ROOT}${WHITE_MATTER_LABEL}

NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#GRAY_MATTER_LABEL} ))
ROOT='';
for(( j=0; j < $NUMBER_OF_REPS; j++ ))
  do
    ROOT=${ROOT}${REPCHARACTER}
  done
GRAY_MATTER_LABEL_FORMAT=${ROOT}${GRAY_MATTER_LABEL}

NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#CSF_MATTER_LABEL} ))
ROOT='';
for(( j=0; j < $NUMBER_OF_REPS; j++ ))
  do
    ROOT=${ROOT}${REPCHARACTER}
  done
CSF_MATTER_LABEL_FORMAT=${ROOT}${CSF_MATTER_LABEL}

SEGMENTATION_PRIOR_WARPED=${SEGMENTATION_PRIOR_WARPED}\%${FORMAT}d.${OUTPUT_SUFFIX}
NUMBER_OF_PRIOR_IMAGES=${#WARPED_PRIOR_IMAGE_FILENAMES[*]}

if [[ ${NUMBER_OF_PRIOR_IMAGES} -lt 4 ]];
  then
    echo "Expected at least 4 prior images (${NUMBER_OF_PRIOR_IMAGES} are specified).  Check the command line specification."
    exit 1
  fi

for(( j=0; j < $NUMBER_OF_PRIOR_IMAGES; j++ ))
  do
    if [[ ! -f ${PRIOR_IMAGE_FILENAMES[$j]} ]];
      then
        echo "Prior image $j ${PRIOR_IMAGE_FILENAMES[$j]} does not exist."
        exit 1
      fi
  done

if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 1 ]];
  then
    if [[ ! -f ${REGISTRATION_TEMPLATE} ]]
      then
        echo "Template for registration, ${REGISTRATION_TEMPLATE}, does not exist."
        exit 1
      fi
  fi

OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi


echoParameters >&2

echo "---------------------  Running `basename $0` on $HOSTNAME  ---------------------"

time_start=`date +%s`


################################################################################
#
# Output images
#
################################################################################

BRAIN_EXTRACTION_MASK=${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION=${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}

################################################################################
#
# Brain extraction
#
################################################################################
EXTRACTED_SEGMENTATION_BRAIN=${OUTPUT_PREFIX}BrainExtractionBrain.${OUTPUT_SUFFIX}
EXTRACTION_GENERIC_AFFINE=${OUTPUT_PREFIX}BrainExtractionPrior0GenericAffine.mat
EXTRACTED_BRAIN_TEMPLATE=${OUTPUT_PREFIX}ExtractedTemplateBrain.${OUTPUT_SUFFIX}
if [[ ! -s ${OUTPUT_PREFIX}ATSStage1Complete.txt ]]; then
if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -eq 1 ]] ; then # BAStages bxt
if [[ ! -f ${BRAIN_EXTRACTION_MASK} ]];
  then

    if [[ ! -f ${BRAIN_TEMPLATE} ]];
      then
        echo "The extraction template doesn't exist:"
        echo "   $BRAIN_TEMPLATE"
        exit 1
      fi
    if [[ ! -f ${EXTRACTION_PRIOR} ]];
      then
        echo "The brain extraction prior doesn't exist:"
        echo "   $EXTRACTION_PRIOR"
        exit 1
      fi

    if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]]
      then
        logCmd ${ANTSPATH}/antsBrainExtraction.sh \
          -d ${DIMENSION} \
          -a ${ANATOMICAL_IMAGES[0]} \
          -e ${BRAIN_TEMPLATE} \
          -f ${EXTRACTION_REGISTRATION_MASK} \
          -m ${EXTRACTION_PRIOR} \
          -o ${OUTPUT_PREFIX} \
          -k ${KEEP_TMP_IMAGES} \
          -s ${OUTPUT_SUFFIX} \
          -q ${USE_FLOAT_PRECISION} \
          -u ${USE_RANDOM_SEEDING} \
          -z ${DEBUG_MODE}
      else
        logCmd ${ANTSPATH}/antsBrainExtraction.sh \
          -d ${DIMENSION} \
          -a ${ANATOMICAL_IMAGES[0]} \
          -e ${BRAIN_TEMPLATE} \
          -m ${EXTRACTION_PRIOR} \
          -o ${OUTPUT_PREFIX} \
          -k ${KEEP_TMP_IMAGES} \
          -s ${OUTPUT_SUFFIX} \
          -q ${USE_FLOAT_PRECISION} \
          -u ${USE_RANDOM_SEEDING} \
          -z ${DEBUG_MODE}
      fi

  fi

if [[ ! -f ${EXTRACTED_SEGMENTATION_BRAIN} ]];
  then
    logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN} m ${ANATOMICAL_IMAGES[0]} ${BRAIN_EXTRACTION_MASK}

    # Do a quick N4 on the brain before registration
    logCmd $N4 -d ${DIMENSION} -i ${EXTRACTED_SEGMENTATION_BRAIN} -s ${N4_SHRINK_FACTOR_1} -c ${N4_CONVERGENCE_1} -o ${EXTRACTED_SEGMENTATION_BRAIN} -x ${BRAIN_EXTRACTION_MASK} -b ${N4_BSPLINE_PARAMS}

  fi

if [[ -f ${BRAIN_TEMPLATE} ]] && [[ ! -f ${EXTRACTED_BRAIN_TEMPLATE} ]];
  then
    logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_PRIOR} ${EXTRACTED_BRAIN_TEMPLATE} 0.1 1.01 1 0
    logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTED_BRAIN_TEMPLATE} m ${BRAIN_TEMPLATE} ${EXTRACTED_BRAIN_TEMPLATE}
  fi
echo ${OUTPUT_PREFIX}ATSStage1Complete.txt > ${OUTPUT_PREFIX}ATSStage1Complete.txt
fi # BAStages
fi # check completion

################################################################################
#
# Brain segmentation
#
################################################################################
SEGMENTATION_WARP=${SEGMENTATION_WARP_OUTPUT_PREFIX}1Warp.nii.gz
SEGMENTATION_INVERSE_WARP=${SEGMENTATION_WARP_OUTPUT_PREFIX}1InverseWarp.nii.gz
SEGMENTATION_GENERIC_AFFINE=${SEGMENTATION_WARP_OUTPUT_PREFIX}0GenericAffine.mat
SEGMENTATION_MASK_DILATED=${BRAIN_SEGMENTATION_OUTPUT}MaskDilated.${OUTPUT_SUFFIX}
SEGMENTATION_CONVERGENCE_FILE=${BRAIN_SEGMENTATION_OUTPUT}Convergence.txt

if [[ ! -s ${OUTPUT_PREFIX}ATSStage2Complete.txt ]]  && \
   [[   -s ${OUTPUT_PREFIX}ATSStage1Complete.txt ]]; then
  if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -eq 2 ]] ; then # BAStages reg
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Brain segmentation using the following steps:"
    echo "   1) Register ${EXTRACTED_BRAIN_TEMPLATE} and ${SEGMENTATION_PRIOR} to ${ANATOMICAL_IMAGES[0]}"
    echo "   2) Warp priors to ${ANATOMICAL_IMAGES[0]}"
    echo "   3) N-tissue segmentation using Atropos and N4"
    echo "--------------------------------------------------------------------------------------"
    echo

    time_start_brain_registration=`date +%s`

    # Check to see if the warped priors exist.  If so, we don't need to warp
    # the template and priors to the individual subject.

    WARPED_PRIORS_ALREADY_EXIST=1;
    for (( i = 0; i < ${NUMBER_OF_PRIOR_IMAGES}; i++ ))
      do
        if [[ ! -f ${WARPED_PRIOR_IMAGE_FILENAMES[$i]} ]];
          then
            WARPED_PRIORS_ALREADY_EXIST=0
            break
          fi
      done

    TMP_FILES=()

    if [[ $WARPED_PRIORS_ALREADY_EXIST -eq 0 ]]
      then

      # Check inputs
      if [[ ! -f ${EXTRACTED_BRAIN_TEMPLATE} ]];
        then
          echo "The segmentation template doesn't exist:"
          echo "   ${EXTRACTED_BRAIN_TEMPLATE}"
          rm -f ${OUTPUT_PREFIX}ATSStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${EXTRACTED_SEGMENTATION_BRAIN} ]];
        then
          echo "The extracted brain doesn't exist:"
          echo "   ${EXTRACTED_SEGMENTATION_BRAIN}"
          rm -f ${OUTPUT_PREFIX}ATSStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${BRAIN_EXTRACTION_MASK} ]];
        then
          echo "The brain extraction mask does not exist:"
          echo "   ${BRAIN_EXTRACTION_MASK}"
          rm -f ${OUTPUT_PREFIX}ATSStage1Complete.txt
          exit 1
        fi

        if [[ ! -f ${SEGMENTATION_WARP} ]];
          then
            logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${SEGMENTATION_MASK_DILATED} MD ${BRAIN_EXTRACTION_MASK} 20

            basecall=''
            if [[ ${RUN_QUICK} -ne 0 ]];
              then
                TMP_FILES=( ${TMP_FILES[@]} "${SEGMENTATION_WARP_OUTPUT_PREFIX}Warped.nii.gz" "${SEGMENTATION_WARP_OUTPUT_PREFIX}InverseWarped.nii.gz" )

                basecall="${ANTSPATH}/antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${EXTRACTED_SEGMENTATION_BRAIN}"
                basecall="${basecall} -m ${EXTRACTED_BRAIN_TEMPLATE} -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} -j 1"
                if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
                  then
                    basecall="${basecall} -p f"
                  fi
            else
              basecall="${ANTS} -d ${DIMENSION} -u 1 -w [ 0.0,0.999 ] -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} --float ${USE_FLOAT_PRECISION} --verbose 1"
              IMAGES="${EXTRACTED_SEGMENTATION_BRAIN},${EXTRACTED_BRAIN_TEMPLATE}"
              if [[ -f ${EXTRACTION_GENERIC_AFFINE} ]];
                then
                  basecall="${basecall} -r [ ${EXTRACTION_GENERIC_AFFINE},1 ]"
                else
                  basecall="${basecall} -r [ ${IMAGES},1 ]"
                fi
              basecall="${basecall} -x [ ${SEGMENTATION_MASK_DILATED} ]"
              stage1="-m MI[ ${IMAGES},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[ 0.1 ] -f 8x4x2x1 -s 4x2x1x0"
              stage2="-m CC[ ${IMAGES},1,4 ] -c [ ${ANTS_MAX_ITERATIONS},1e-9,15 ] -t ${ANTS_TRANSFORMATION} -f 6x4x2x1 -s 3x2x1x0"

              basecall="${basecall} ${stage1} ${stage2}"
            fi

            exe_brain_segmentation_1=${basecall}


            # Precision errors in .nii (which stores things as float) headers can cause problems, so attempt to make everything consistent.
            # Won't be perfectly consistent because we don't change ${ANATOMICAL_IMAGES[0]} and CopyImageHeaderInformation does not make
            # a perfect copy. But hopefully close enough
            for img in ${BRAIN_EXTRACTION_MASK} ${EXTRACTED_SEGMENTATION_BRAIN} ${SEGMENTATION_MASK_DILATED};
              do
                logCmd ${ANTSPATH}/CopyImageHeaderInformation ${ANATOMICAL_IMAGES[0]} ${img} ${img} 1 1 1
              done

            logCmd $exe_brain_segmentation_1

          fi

        ## check to see if the output registration transforms exist
        if [[ ! -f ${SEGMENTATION_GENERIC_AFFINE} ]];
          then
            echo "The registration component of the segmentation step didn't complete properly."
            echo "The transform file ${SEGMENTATION_GENERIC_AFFINE} does not exist."
            exit 1
          fi

        if [[ ! -f ${SEGMENTATION_WARP} ]];
          then
            echo "The registration component of the segmentation step didn't complete properly."
            echo "The transform file ${SEGMENTATION_WARP} does not exist."
            exit 1
          fi

        ## Step 2 ##

        for (( i = 0; i < ${NUMBER_OF_PRIOR_IMAGES}; i++ ))
          do
            if [[ ! -f ${PRIOR_IMAGE_FILENAMES[$i]} ]];
              then
                echo "The prior image file name does not exist:"
                echo "   ${PRIOR_IMAGE_FILENAMES[$i]}"
                exit 1
              fi

            exe_brain_segmentation_2="${WARP} -d ${DIMENSION} -i ${PRIOR_IMAGE_FILENAMES[$i]} -o ${WARPED_PRIOR_IMAGE_FILENAMES[$i]} -r ${ANATOMICAL_IMAGES[0]} -n Gaussian  -t ${SEGMENTATION_WARP} -t ${SEGMENTATION_GENERIC_AFFINE} --float ${USE_FLOAT_PRECISION} --verbose 1"

            logCmd $exe_brain_segmentation_2
          done
      fi
      echo ${OUTPUT_PREFIX}ATSStage2Complete.txt > ${OUTPUT_PREFIX}ATSStage2Complete.txt
      time_end_brain_registration=`date +%s`
      time_elapsed_brain_registration=$((time_end_brain_registration - time_start_brain_registration))
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " Done with brain registration:  $(( time_elapsed_brain_segmentation / 3600 ))h $(( time_elapsed_brain_registration %3600 / 60 ))m $(( time_elapsed_brain_registration % 60 ))s"
      echo "--------------------------------------------------------------------------------------"
      echo
    fi # BAStages check reg
  fi # BAStages check reg

# We do two stages of antsAtroposN4.  The first stage is to get a good N4
# bias corrected image(s).  This bias corrected image(s) is used as input to the
# second stage where we only do 2 iterations.
if [[ ! -s ${OUTPUT_PREFIX}ATSStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage1Complete.txt ]] ; then
  if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -eq 3 ]] ; then # BAStages seg
    time_start_brain_segmentation=`date +%s`
    ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE='';
    for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
      do
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${ANATOMICAL_IMAGES[$j]}"
      done

    ATROPOS_LABEL_PROPAGATION_COMMAND_LINE=''
    for (( j = 0; j < ${#ATROPOS_SEGMENTATION_LABEL_PROPAGATION[@]}; j++ ))
      do
        ATROPOS_LABEL_PROPAGATION_COMMAND_LINE="${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} -l ${ATROPOS_SEGMENTATION_LABEL_PROPAGATION[$j]}";
      done


    # include everything but the csf
    N4_INCLUDE_PRIORS_COMMAND_LINE=''
    for (( j = 2; j <= ${NUMBER_OF_PRIOR_IMAGES}; j++ ))
      do
        N4_INCLUDE_PRIORS_COMMAND_LINE="${N4_INCLUDE_PRIORS_COMMAND_LINE} -y $j";
      done

    logCmd ${ANTSPATH}/antsAtroposN4.sh \
      -d ${DIMENSION} \
      -b "${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}" \
      ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} \
      ${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} \
      -x ${BRAIN_EXTRACTION_MASK} \
      -m ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS} \
      -n ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
      -c ${NUMBER_OF_PRIOR_IMAGES} \
      ${N4_INCLUDE_PRIORS_COMMAND_LINE} \
      -p ${SEGMENTATION_PRIOR_WARPED} \
      -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT} \
      -o ${OUTPUT_PREFIX}Brain \
      -u ${USE_RANDOM_SEEDING} \
      -g ${DENOISE_ANATOMICAL_IMAGES} \
      -k ${KEEP_TMP_IMAGES} \
      -s ${OUTPUT_SUFFIX} \
      -z ${DEBUG_MODE}

    ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE=''
    for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
      do
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${OUTPUT_PREFIX}BrainSegmentation${j}N4.${OUTPUT_SUFFIX}";
      done

    ## Don't de-noise a second time
    logCmd ${ANTSPATH}/antsAtroposN4.sh \
      -d ${DIMENSION} \
      -b "${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}" \
      ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} \
      ${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} \
      -x ${BRAIN_EXTRACTION_MASK} \
      -m 2 \
      -n ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
      -c ${NUMBER_OF_PRIOR_IMAGES} \
      ${N4_INCLUDE_PRIORS_COMMAND_LINE} \
      -p ${SEGMENTATION_PRIOR_WARPED} \
      -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT} \
      -o ${OUTPUT_PREFIX}Brain \
      -u ${USE_RANDOM_SEEDING} \
      -g 0 \
      -k ${KEEP_TMP_IMAGES} \
      -s ${OUTPUT_SUFFIX} \
      -z ${DEBUG_MODE}

    ## Step 3 ###
    TMP_FILES=( ${TMP_FILES[@]} $EXTRACTION_GENERIC_AFFINE $EXTRACTED_SEGMENTATION_BRAIN $SEGMENTATION_MASK_DILATED $EXTRACTED_BRAIN_TEMPLATE )
    TMP_FILES=( ${TMP_FILES[@]} ${WARPED_PRIOR_IMAGE_FILENAMES[@]}  )
    if [[ $TEMPLATES_ARE_IDENTICAL -eq 0 ]];
      then
        TMP_FILES=( ${TMP_FILES[@]} $SEGMENTATION_WARP $SEGMENTATION_INVERSE_WARP $SEGMENTATION_GENERIC_AFFINE )
      fi

    if [[ $KEEP_TMP_IMAGES -eq 0 ]];
      then
        for f in ${TMP_FILES[@]}
          do
            if [[ -e $f ]];
          then
            logCmd rm -f $f
          else
            echo "WARNING: expected temp file doesn't exist: $f"
          fi
      done
    fi


    time_end_brain_segmentation=`date +%s`
    time_elapsed_brain_segmentation=$((time_end_brain_segmentation - time_start_brain_segmentation))

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with brain segmentation:  $(( time_elapsed_brain_segmentation / 3600 ))h $(( time_elapsed_brain_segmentation %3600 / 60 ))m $(( time_elapsed_brain_segmentation % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo

   echo ${OUTPUT_PREFIX}ATSStage3Complete.txt > ${OUTPUT_PREFIX}ATSStage3Complete.txt
  fi
fi # BAStages seg

################################################################################
#
# Registration to a template
#
################################################################################

# These affect output; keep them consistent with usage and checkOutputExists function
REGISTRATION_TEMPLATE_OUTPUT_PREFIX=${OUTPUT_PREFIX}SubjectToTemplate
REGISTRATION_TEMPLATE_GENERIC_AFFINE=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}0GenericAffine.mat
REGISTRATION_TEMPLATE_WARP=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1Warp.nii.gz
REGISTRATION_TEMPLATE_INVERSE_WARP=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1InverseWarp.nii.gz
REGISTRATION_LOG_JACOBIAN=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}LogJacobian.${OUTPUT_SUFFIX}

# Want to have transforms for both directions
REGISTRATION_SUBJECT_OUTPUT_PREFIX=${OUTPUT_PREFIX}TemplateToSubject
REGISTRATION_SUBJECT_GENERIC_AFFINE=${REGISTRATION_SUBJECT_OUTPUT_PREFIX}1GenericAffine.mat
REGISTRATION_SUBJECT_WARP=${REGISTRATION_SUBJECT_OUTPUT_PREFIX}0Warp.nii.gz

# Use first N4 corrected segmentation image, which we assume to be T1
HEAD_N4_IMAGE=${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE=${OUTPUT_PREFIX}ExtractedBrain0N4.${OUTPUT_SUFFIX}

if [[ ! -s ${OUTPUT_PREFIX}ATSStage4Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage1Complete.txt ]] ; then
if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -eq 4 ]] ; then # BAStages reg

if [[ -f ${REGISTRATION_TEMPLATE} ]] && [[ ! -f $REGISTRATION_LOG_JACOBIAN ]];
  then

    TMP_FILES=()

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Registration brain masked ${HEAD_N4_IMAGE} to ${REGISTRATION_TEMPLATE} "
    echo "--------------------------------------------------------------------------------------"
    echo

    logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} m ${HEAD_N4_IMAGE} ${BRAIN_EXTRACTION_MASK}

    TMP_FILES=( ${TMP_FILES[@]} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} )

    time_start_template_registration=`date +%s`

    basecall=''
    if [[ ${RUN_QUICK} -ne 0 ]];
      then
        TMP_FILES=( ${TMP_FILES[@]} "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}Warped.nii.gz" "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}InverseWarped.nii.gz" )

        basecall="${ANTSPATH}/antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${REGISTRATION_TEMPLATE}"
        basecall="${basecall} -m ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -j 1"
        if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
          then
            basecall="${basecall} -p f"
          fi
      else
        IMAGES="${REGISTRATION_TEMPLATE},${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE}"
        basecall="${ANTS} -d ${DIMENSION} -v 1 -u 1 -w [ 0.0,0.999 ] -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -r [ ${IMAGES},1 ] --float ${USE_FLOAT_PRECISION}"
        stage1="-m MI[ ${IMAGES},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[ 0.1 ] -f 8x4x2x1 -s 3x2x1x0"
        stage2="-m MI[ ${IMAGES},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[ 0.1 ] -f 8x4x2x1 -s 3x2x1x0"
        stage3="-m CC[ ${IMAGES},1,4 ] -c [ ${ANTS_MAX_ITERATIONS},1e-9,15 ] -t ${ANTS_TRANSFORMATION} -f 6x4x2x1 -s 3x2x1x0"
        basecall="${basecall} ${stage1} ${stage2} ${stage3}"
      fi
    exe_template_registration_1="${basecall}"

    if [[ ! -f ${REGISTRATION_TEMPLATE_WARP} ]];
      then
        logCmd $exe_template_registration_1
      fi

    ## check to see if the output registration transforms exist
    if [[ ! -f ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} ]];
      then
        echo "The registration component of the segmentation step didn't complete properly."
        echo "The transform file ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} does not exist."
        exit 1
      fi

    if [[ ! -f ${REGISTRATION_TEMPLATE_WARP} ]];
      then
        echo "The registration component of the segmentation step didn't complete properly."
        echo "The transform file ${REGISTRATION_TEMPLATE_WARP} does not exist."
        exit 1
      fi

    ## Create symmetric transforms for template to subject warping
    if [[ -s ${REGISTRATION_TEMPLATE_INVERSE_WARP} ]] && [[ ! -s ${REGISTRATION_SUBJECT_WARP} ]] ; then
      logCmd mv ${REGISTRATION_TEMPLATE_INVERSE_WARP} ${REGISTRATION_SUBJECT_WARP}
    fi
    if [[ ! -s  ${REGISTRATION_SUBJECT_WARP} ]] ; then
      echo "The transform file ${REGISTRATION_SUBJECT_WARP} does not exist."
      exit 1
    fi
    logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -o Linear[ $REGISTRATION_SUBJECT_GENERIC_AFFINE,1 ] -t $REGISTRATION_TEMPLATE_GENERIC_AFFINE --verbose 1

    time_end_template_registration=`date +%s`
    time_elapsed_template_registration=$((time_end_template_registration - time_start_template_registration))

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with registration:  $(( time_elapsed_template_registration / 3600 ))h $(( time_elapsed_template_registration %3600 / 60 ))m $(( time_elapsed_template_registration % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo

    if [[ $KEEP_TMP_IMAGES -eq 0 ]];
      then
        for f in ${TMP_FILES[@]}
          do
            if [[ -e $f ]];
              then
                logCmd rm -f $f
              else
                echo "WARNING: expected temp file doesn't exist: $f"
              fi
          done
      fi
  logCmd ${ANTSPATH}/CreateJacobianDeterminantImage ${DIMENSION} ${REGISTRATION_TEMPLATE_WARP} ${REGISTRATION_LOG_JACOBIAN} 1 1
  fi # if registration template & jacobian check
  if [[ -s ${REGISTRATION_LOG_JACOBIAN} ]] ; then
    echo ${OUTPUT_PREFIX}ATSStage4Complete.txt > ${OUTPUT_PREFIX}ATSStage4Complete.txt
  fi
fi # BAStages
fi # check completion


################################################################################
#
# Create QA/QC output:
#   - tiled mosaic of ${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
################################################################################
if [[ ! -s ${OUTPUT_PREFIX}ATSStage5Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ATSStage1Complete.txt ]] ; then
if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -eq 5 ]] ; then # BAStages qc
HEAD_N4_IMAGE=${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
HEAD_N4_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentation0N4Resampled.${OUTPUT_SUFFIX}"
BRAIN_EXTRACTION_MASK_RESAMPLED="${OUTPUT_PREFIX}BrainExtractionMaskResampled.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentationResampled.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_IMAGE_RGB="${OUTPUT_PREFIX}BrainSegmentationRGB.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_MOSAIC="${OUTPUT_PREFIX}BrainSegmentationTiledMosaic.png"
ITKSNAP_COLORMAP="${OUTPUT_PREFIX}ItkSnapColormap.txt"

if [[ ! -f ${BRAIN_SEGMENTATION_MOSAIC} ]];
  then
    TMP_FILES=( $BRAIN_SEGMENTATION_IMAGE_RGB $ITKSNAP_COLORMAP )
    TMP_FILES=( ${TMP_FILES[@]} $HEAD_N4_IMAGE_RESAMPLED $BRAIN_EXTRACTION_MASK_RESAMPLED $BRAIN_SEGMENTATION_IMAGE_RESAMPLED )

    # Resample images

    resample0="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${BRAIN_EXTRACTION_MASK}"
    resample0="${resample0} ${BRAIN_EXTRACTION_MASK_RESAMPLED} 1 1"
    resample1="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${HEAD_N4_IMAGE}"
    resample1="${resample1} ${HEAD_N4_IMAGE_RESAMPLED} 1 1"
    resample2="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${BRAIN_SEGMENTATION}"
    resample2="${resample2} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED} 1 1"
    if [[ ${DIMENSION} -eq 3 ]];
      then
        resample0="${resample0} 1 0 0 1"
        resample1="${resample1} 1"
        resample2="${resample2} 1 0 0 1"
      else
        resample0="${resample0} 0 0 1"
        resample2="${resample2} 0 0 1"
      fi
    logCmd $resample0
    logCmd $resample1
    logCmd $resample2

    # Segmentation

    echo "0 1 0 0 1 0 1" > $ITKSNAP_COLORMAP
    echo "0 0 1 0 1 1 0" >> $ITKSNAP_COLORMAP
    echo "0 0 0 1 0 1 1" >> $ITKSNAP_COLORMAP

    conversion="${ANTSPATH}/ConvertScalarImageToRGB ${DIMENSION} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED}"
    conversion="${conversion} ${BRAIN_SEGMENTATION_IMAGE_RGB} none custom $ITKSNAP_COLORMAP 0 6"
    logCmd $conversion

    mosaic="${ANTSPATH}/CreateTiledMosaic -i ${HEAD_N4_IMAGE_RESAMPLED} -r ${BRAIN_SEGMENTATION_IMAGE_RGB}"
    mosaic="${mosaic} -o ${BRAIN_SEGMENTATION_MOSAIC} -a 0.3 -t -1x-1 -d 2 -p mask"
    mosaic="${mosaic} -s [ 2,mask,mask ] -x ${BRAIN_EXTRACTION_MASK_RESAMPLED}"
    logCmd $mosaic


    if [[ $KEEP_TMP_IMAGES -eq 0 ]];
      then
        for f in ${TMP_FILES[@]}
          do
            if [[ -e $f ]];
             then
              logCmd rm -f $f
            else
              echo "WARNING: expected temp file doesn't exist: $f"
            fi
        done
      fi
  fi
  echo ${OUTPUT_PREFIX}ATSStage5Complete.txt > ${OUTPUT_PREFIX}ATSStage5Complete.txt
fi # BAStages
fi # check completion
################################################################################
#
# End of main routine
#
################################################################################

if [[ ${ATS_STAGE} -eq 0 ]] || [[ ${ATS_STAGE} -ge 4 ]] ; then
  logCmd checkOutputExists
fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with ANTs processing pipeline: stage $ATS_STAGE"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
