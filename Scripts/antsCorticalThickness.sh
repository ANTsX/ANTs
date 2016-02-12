#!/bin/bash

VERSION="0.0"

# Check dependencies

PROGRAM_DEPENDENCIES=( 'antsRegistration' 'antsApplyTransforms' 'N4BiasFieldCorrection' 'Atropos' 'KellyKapowski' )
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

`basename $0` performs T1 anatomical brain processing where the following steps are currently applied:

  1. Brain extraction
  2. Brain n-tissue segmentation
  3. Cortical thickness
  4. (Optional) registration to a template

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
                                                  * ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}

Optional arguments:

     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd
     -t:  template for t1 registration          Anatomical *intensity* template (assumed to be skull-stripped).  A common
                                                use case would be where this would be the same template as specified in the
                                                -e option which is not skull stripped.
                                                We perform the registration (fixed image = individual subject
                                                and moving image = template) to produce the files.
                                                The output from this step is
                                                  * ${OUTPUT_PREFIX}TemplateToSubject0GenericAffine.mat
                                                  * ${OUTPUT_PREFIX}TemplateToSubject1Warp.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}TemplateToSubject1InverseWarp.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}TemplateToSubjectLogJacobian.${OUTPUT_SUFFIX}
     -f:  extraction registration mask          Mask (defined in the template space) used during registration
                                                for brain extraction.
     -k:  keep temporary files                  Keep brain extraction/segmentation warps, etc (default = 0).
     -g:  denoise anatomical images             Denoise anatomical images (default = 0).
     -i:  max iterations for registration       ANTS registration max iterations (default = 100x100x70x20)
     -w:  Atropos prior segmentation weight     Atropos spatial prior *probability* weight for the segmentation (default = 0.25)
     -n:  number of segmentation iterations     N4 -> Atropos -> N4 iterations during segmentation (default = 3)
     -b:  posterior formulation                 Atropos posterior formulation and whether or not to use mixture model proportions.
                                                e.g 'Socrates[1]' (default) or 'Aristotle[1]'.  Choose the latter if you
                                                want use the distance priors (see also the -l option for label propagation
                                                control).
     -j:  use floating-point precision          Use floating point precision in registrations (default = 0)
     -u:  use random seeding                    Use random number generated from system clock in Atropos (default = 1)
     -v:  use b-spline smoothing                Use B-spline SyN for registrations and B-spline exponential mapping in DiReCT.
     -r:  cortical label image                  Cortical ROI labels to use as a prior for ATITH.
     -l:  label propagation                     Incorporate a distance prior one the posterior formulation.  Should be
                                                of the form 'label[lambda,boundaryProbability]' where label is a value
                                                of 1,2,3,... denoting label ID.  The label probability for anything
                                                outside the current label

                                                  = boundaryProbability * exp( -lambda * distanceFromBoundary )

                                                Intuitively, smaller lambda values will increase the spatial capture
                                                range of the distance prior.  To apply to all label values, simply omit
                                                specifying the label, i.e. -l [lambda,boundaryProbability].
     -c                                         Add prior combination to combined gray and white matters.  For example,
                                                when calling KK for normal subjects, we combine the deep gray matter
                                                segmentation/posteriors with the white matter segmentation/posteriors.
                                                An additional example would be performing cortical thickness in the presence
                                                of white matter lesions.  We can accommodate this by specifying a lesion mask
                                                posterior as an additional posterior (suppose label '7'), and then combine
                                                this with white matter by specifying '-c WM[7]' or '-c 3[7]'.
     -q:  Use quick registration parameters     If = 1, use antsRegistrationSyNQuick.sh as the basis for registration
                                                during brain extraction, brain segmentation, and (optional) normalization
                                                to a template.  Otherwise use antsRegistrationSyN.sh (default = 0).
     -x:                                        Number of iterations within Atropos (default 5).
     -y:                                        Which stage of ACT to run (default = 0, run all).  Tries to split in 2 hour chunks.
                                                Will produce OutputNameACTStageNComplete.txt for each completed stage.
                                                1: brain extraction
                                                2: template registration
                                                3: tissue segmentation
                                                4: template registration (improved, optional)
                                                5: DiReCT cortical thickness
                                                6: qc, quality control and summary measurements
     -z:  Test / debug mode                     If > 0, runs a faster version of the script. Only for testing. Implies -u 0.
                                                Requires single thread computation for complete reproducibility.
USAGE
    exit 1
}

# Check outputs exist, runs at the end of the script
# List of outputs is taken from the usage
function checkOutputExists() {

  singleOutputs=( ${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX} )
  singleOutputs=( ${singleOutputs[@]} ${OUTPUT_PREFIX}BrainSegmentationTiledMosaic.png ${OUTPUT_PREFIX}CorticalThicknessTiledMosaic.png )

  if [[ -f ${REGISTRATION_TEMPLATE} ]];
    then
      singleOutputs=( ${singleOutputs[@]} ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}0GenericAffine.mat ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1Warp.${OUTPUT_SUFFIX} ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1InverseWarp.${OUTPUT_SUFFIX} ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}LogJacobian.${OUTPUT_SUFFIX} )
    fi

  missingOutput=0

  for img in $singleOutputs;
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

    Using antsCorticalThickness with the following arguments:
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
      metric                  = ${ANTS_METRIC}[fixedImage,movingImage,${ANTS_METRIC_PARAMS}]
      regularization          = ${ANTS_REGULARIZATION}
      transformation          = ${ANTS_TRANSFORMATION}
      max iterations          = ${ANTS_MAX_ITERATIONS}

    DiReCT parameters:
      convergence             = ${DIRECT_CONVERGENCE}
      thickness prior         = ${DIRECT_THICKNESS_PRIOR}
      gradient step size      = ${DIRECT_GRAD_STEP_SIZE}
      smoothing sigma         = ${DIRECT_SMOOTHING_PARAMETER}

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

ACT_STAGE=0 # run all stages

function logCmd() {
  cmd="$*"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  $cmd

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
OUTPUT_SUFFIX="nii.gz"

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
CORTICAL_LABEL_IMAGE=""

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
ANTS_TRANSFORMATION="SyN[0.1,3,0]"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[1000x500x250x100,1e-8,10]"
ANTS_METRIC="CC"
ANTS_METRIC_PARAMS="1,4"

WARP=${ANTSPATH}/antsApplyTransforms

N4=${ANTSPATH}/N4BiasFieldCorrection
N4_CONVERGENCE_1="[50x50x50x50,0.0000001]"
N4_CONVERGENCE_2="[50x50x50x50,0.0000001]"
N4_SHRINK_FACTOR_1=4
N4_SHRINK_FACTOR_2=2
N4_BSPLINE_PARAMS="[200]"

ATROPOS=${ANTSPATH}/Atropos

ATROPOS_SEGMENTATION_INITIALIZATION="PriorProbabilityImages"
ATROPOS_SEGMENTATION_LIKELIHOOD="Gaussian"
ATROPOS_SEGMENTATION_CONVERGENCE="[5,0.0]"
ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION="Socrates[1]"
ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=3
ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5 # to be backward compatible but i like 25
ATROPOS_SEGMENTATION_LABEL_PROPAGATION=()

DIRECT=${ANTSPATH}/KellyKapowski
DIRECT_CONVERGENCE="[45,0.0,10]"
DIRECT_THICKNESS_PRIOR="10"
DIRECT_GRAD_STEP_SIZE="0.025"
DIRECT_SMOOTHING_PARAMETER="1.5"
DIRECT_NUMBER_OF_DIFF_COMPOSITIONS="10"

PRIOR_COMBINATIONS=( 'WM[4]' )

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
          r) #cortical label image
       CORTICAL_LABEL_IMAGE=$OPTARG
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
       ACT_STAGE=$OPTARG
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
    ANTS_TRANSFORMATION="BSplineSyN[0.1,26,0,3]"
    DIRECT_SMOOTHING_PARAMETER="5.75"
  fi

if [[ $DEBUG_MODE -gt 0 ]];
  then

    echo "    WARNING - Running in test / debug mode. Results will be suboptimal "

    OUTPUT_PREFIX="${OUTPUT_PREFIX}testMode_"

    # Speed up by doing fewer its. Careful about changing this because
    # certain things are hard coded elsewhere, eg number of levels

    ANTS_MAX_ITERATIONS="40x40x20x0"
    ANTS_LINEAR_CONVERGENCE="[100x100x50x0,1e-8,10]"
    ANTS_METRIC_PARAMS="1,2"

    # I think this is the number of times we run the whole N4 / Atropos thing, at the cost of about 10 minutes a time
    ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=1
    ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5 # internal to atropos

    DIRECT_CONVERGENCE="[5,0.0,10]"

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
if [ ${#FORMAT} -eq 2 ]
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

# These arrays contain the formatted labels that will be used to combine with the white
# and gray matter posteriors for the cortical thickness section.

CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS=()
CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS=()

for(( j=0; j < ${#PRIOR_COMBINATIONS[@]}; j++ ))
  do
    echo ${PRIOR_COMBINATIONS[$j]}
    COMBINATION=( $( echo ${PRIOR_COMBINATIONS[$j]} | tr "[]," "\n" ) )

    echo ${COMBINATION[@]}

    if [[ ${COMBINATION[0]} == ${WHITE_MATTER_LABEL} || ${COMBINATION[0]} == 'WM' ]];
      then
        for(( k=1; k < ${#COMBINATION[@]}; k++ ))
          do
            OTHER_LABEL=${COMBINATION[$k]}
            CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[${#CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[@]}]=${OTHER_LABEL}
          done
      elif [[ ${COMBINATION[0]} == ${GRAY_MATTER_LABEL} || ${COMBINATION[0]} == 'GM' ]];
        then
          for(( k=1; k < ${#COMBINATION[@]}; k++ ))
            do
              OTHER_LABEL=${COMBINATION[$k]}
              CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS[${#CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS[@]}]=${OTHER_LABEL}
            done
      else
        echo "We only combine with the gray matter or the white matter."
        echo "The label ${COMBINATION[0]} does not correspond to the gray or white matters."
        exit 1
      fi
  done
CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS_FORMAT=( $( echo "${CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS_FORMAT[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ' ) )
CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS_FORMAT=( $( echo "${CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS_FORMAT[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ' ) )

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
CORTICAL_THICKNESS_IMAGE=${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}

################################################################################
#
# Brain extraction
#
################################################################################
EXTRACTED_SEGMENTATION_BRAIN=${OUTPUT_PREFIX}BrainExtractionBrain.${OUTPUT_SUFFIX}
EXTRACTION_GENERIC_AFFINE=${OUTPUT_PREFIX}BrainExtractionPrior0GenericAffine.mat
EXTRACTED_BRAIN_TEMPLATE=${OUTPUT_PREFIX}ExtractedTemplateBrain.${OUTPUT_SUFFIX}
if [[ ! -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]]; then
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 1  ]] ; then # BAStages bxt
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
  fi

if [[ -f ${BRAIN_TEMPLATE} ]] && [[ ! -f ${EXTRACTED_BRAIN_TEMPLATE} ]];
  then
    logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_PRIOR} ${EXTRACTED_BRAIN_TEMPLATE} 0.1 1.01 1 0
    logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTED_BRAIN_TEMPLATE} m ${BRAIN_TEMPLATE} ${EXTRACTED_BRAIN_TEMPLATE}
  fi
echo ${OUTPUT_PREFIX}ACTStage1Complete.txt > ${OUTPUT_PREFIX}ACTStage1Complete.txt
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
SEGMENTATION_MASK_DILATED=${BRAIN_SEGMENTATION_OUTPUT}MaskDilated.nii.gz
SEGMENTATION_CONVERGENCE_FILE=${BRAIN_SEGMENTATION_OUTPUT}Convergence.txt

if [[ ! -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]]  && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]]; then
  if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 2  ]] ; then # BAStages reg
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
          rm ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${EXTRACTED_SEGMENTATION_BRAIN} ]];
        then
          echo "The extracted brain doesn't exist:"
          echo "   ${EXTRACTED_SEGMENTATION_BRAIN}"
          rm ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${BRAIN_EXTRACTION_MASK} ]];
        then
          echo "The brain extraction mask does not exist:"
          echo "   ${BRAIN_EXTRACTION_MASK}"
          rm ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi

        if [[ ! -f ${SEGMENTATION_WARP} ]];
          then
            logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${SEGMENTATION_MASK_DILATED} MD ${BRAIN_EXTRACTION_MASK} 20

            basecall=''
            if [[ ${RUN_QUICK} -ne 0 ]];
              then
                TMP_FILES=( ${TMP_FILES[@]} "${SEGMENTATION_WARP_OUTPUT_PREFIX}Warped.nii.gz" )

                basecall="${ANTSPATH}/antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${EXTRACTED_SEGMENTATION_BRAIN}"
                basecall="${basecall} -m ${EXTRACTED_BRAIN_TEMPLATE} -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} -j 1"
                if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
                  then
                    basecall="${basecall} -p f"
                  fi
            else
              basecall="${ANTS} -d ${DIMENSION} -u 1 -w [0.01,0.99] -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} --float ${USE_FLOAT_PRECISION} --verbose 1"
              IMAGES="${EXTRACTED_SEGMENTATION_BRAIN},${EXTRACTED_BRAIN_TEMPLATE}"
              if [[ -f ${EXTRACTION_GENERIC_AFFINE} ]];
                then
                  basecall="${basecall} -r [${EXTRACTION_GENERIC_AFFINE},1]"
                else
                  basecall="${basecall} -r [${IMAGES},1]"
                fi
              basecall="${basecall} -x [${SEGMENTATION_MASK_DILATED}]"
              stage1="-m MI[${IMAGES},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0"
              stage2="-m CC[${IMAGES},1,4] -c [${ANTS_MAX_ITERATIONS},1e-9,15] -t ${ANTS_TRANSFORMATION} -f 6x4x2x1 -s 3x2x1x0"

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
      echo ${OUTPUT_PREFIX}ACTStage2Complete.txt > ${OUTPUT_PREFIX}ACTStage2Complete.txt
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
if [[ ! -s ${OUTPUT_PREFIX}ACTStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]] ; then
  if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 3  ]] ; then # BAStages seg
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
      -b ${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION} \
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

    logCmd ${ANTSPATH}/antsAtroposN4.sh \
      -d ${DIMENSION} \
      -b ${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION} \
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
      -g ${DENOISE_ANATOMICAL_IMAGES} \
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
            logCmd rm $f
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

   echo ${OUTPUT_PREFIX}ACTStage3Complete.txt > ${OUTPUT_PREFIX}ACTStage3Complete.txt
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
REGISTRATION_TEMPLATE_WARP=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1Warp.${OUTPUT_SUFFIX}
REGISTRATION_TEMPLATE_INVERSE_WARP=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1InverseWarp.${OUTPUT_SUFFIX}
REGISTRATION_LOG_JACOBIAN=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}LogJacobian.${OUTPUT_SUFFIX}

# Want to have transforms for both directions
REGISTRATION_SUBJECT_OUTPUT_PREFIX=${OUTPUT_PREFIX}TemplateToSubject
REGISTRATION_SUBJECT_GENERIC_AFFINE=${REGISTRATION_SUBJECT_OUTPUT_PREFIX}1GenericAffine.mat
REGISTRATION_SUBJECT_WARP=${REGISTRATION_SUBJECT_OUTPUT_PREFIX}0Warp.${OUTPUT_SUFFIX}

# Use first N4 corrected segmentation image, which we assume to be T1
HEAD_N4_IMAGE=${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE=${OUTPUT_PREFIX}ExtractedBrain0N4.nii.gz

if [[ ! -s ${OUTPUT_PREFIX}ACTStage4Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]] ; then
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 4 ]] ; then # BAStages reg

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
        TMP_FILES=( ${TMP_FILES[@]} "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}Warped.nii.gz" )

        basecall="${ANTSPATH}/antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${REGISTRATION_TEMPLATE}"
        basecall="${basecall} -m ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -j 1"
        if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
          then
            basecall="${basecall} -p f"
          fi
      else
        IMAGES="${REGISTRATION_TEMPLATE},${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE}"
        basecall="${ANTS} -d ${DIMENSION} -v 1 -u 1 -w [0.01,0.99] -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -r [${IMAGES},1] --float ${USE_FLOAT_PRECISION}"
        stage1="-m MI[${IMAGES},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[0.1] -f 8x4x2x1 -s 3x2x1x0"
        stage2="-m MI[${IMAGES},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0"
        stage3="-m CC[${IMAGES},1,4] -c [${ANTS_MAX_ITERATIONS},1e-9,15] -t ${ANTS_TRANSFORMATION} -f 6x4x2x1 -s 3x2x1x0"
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
    logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -o Linear[$REGISTRATION_SUBJECT_GENERIC_AFFINE,1] -t $REGISTRATION_TEMPLATE_GENERIC_AFFINE --verbose 1

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
              logCmd rm $f
            else
              echo "WARNING: expected temp file doesn't exist: $f"
            fi
        done
      fi
  logCmd ${ANTSPATH}/CreateJacobianDeterminantImage ${DIMENSION} ${REGISTRATION_TEMPLATE_WARP} ${REGISTRATION_LOG_JACOBIAN} 1 1
  fi # if registration template & jacobian check
  if [[ -s ${REGISTRATION_LOG_JACOBIAN} ]] ; then
    echo ${OUTPUT_PREFIX}ACTStage4Complete.txt > ${OUTPUT_PREFIX}ACTStage4Complete.txt
  fi
fi # BAStages
fi # check completion

################################################################################
#
# Cortical thickness
#
################################################################################

BRAIN_SEGMENTATION_GM=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${GRAY_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION_WM=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${WHITE_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION_DEEP_GM=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${DEEP_GRAY_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}

CORTICAL_THICKNESS_GM=${OUTPUT_PREFIX}CorticalThicknessPosteriors${GRAY_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_WM=${OUTPUT_PREFIX}CorticalThicknessPosteriors${WHITE_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_SEGMENTATION=${OUTPUT_PREFIX}CorticalThicknessSegmentation.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_GM_SEGMENTATION=${OUTPUT_PREFIX}CorticalThicknessGMSegmentation.${OUTPUT_SUFFIX}
CORTICAL_LABEL_THICKNESS_PRIOR=${OUTPUT_PREFIX}CorticalLabelThicknessPrior.${OUTPUT_SUFFIX}

if [[ ! -s ${OUTPUT_PREFIX}ACTStage5Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]] ; then
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 5  ]] ; then # BAStages thk
if [[ ! -f ${CORTICAL_THICKNESS_IMAGE} ]];
  then

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Cortical thickness using DiReCT (KellyKapowski)"
    echo "--------------------------------------------------------------------------------------"
    echo

    # combine posteriors

    logCmd cp ${BRAIN_SEGMENTATION_GM} ${CORTICAL_THICKNESS_GM}
    for(( i=0; i < ${#CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS[@]}; i++ ))
      do
        OTHER_LABEL=${CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS[$i]}
        NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#OTHER_LABEL} ))
        ROOT='';
        for(( l=0; l < $NUMBER_OF_REPS; l++ ))
          do
            ROOT=${ROOT}${REPCHARACTER}
          done
        OTHER_LABEL_FORMAT=${ROOT}${OTHER_LABEL}
        BRAIN_SEGMENTATION_OTHER_LABEL=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${OTHER_LABEL_FORMAT}.${OUTPUT_SUFFIX}

        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_GM} + ${CORTICAL_THICKNESS_GM} ${BRAIN_SEGMENTATION_OTHER_LABEL}
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ReplaceVoxelValue ${BRAIN_SEGMENTATION} ${OTHER_LABEL} ${OTHER_LABEL} ${GRAY_MATTER_LABEL}
      done

    logCmd cp ${BRAIN_SEGMENTATION_WM} ${CORTICAL_THICKNESS_WM}
    for(( i=0; i < ${#CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[@]}; i++ ))
      do
        OTHER_LABEL=${CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[$i]}
        NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#OTHER_LABEL} ))
        ROOT='';
        for(( l=0; l < $NUMBER_OF_REPS; l++ ))
          do
            ROOT=${ROOT}${REPCHARACTER}
          done
        OTHER_LABEL_FORMAT=${ROOT}${OTHER_LABEL}
        BRAIN_SEGMENTATION_OTHER_LABEL=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${OTHER_LABEL_FORMAT}.${OUTPUT_SUFFIX}

        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_WM} + ${CORTICAL_THICKNESS_WM} ${BRAIN_SEGMENTATION_OTHER_LABEL}
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ReplaceVoxelValue ${BRAIN_SEGMENTATION} ${OTHER_LABEL} ${OTHER_LABEL} ${WHITE_MATTER_LABEL}
      done

    # Check inputs
    if [[ ! -f ${CORTICAL_THICKNESS_SEGMENTATION} ]];
      then
        echo "The brain segmentation image doesn't exist:"
        echo "   $CORTICAL_THICKNESS_SEGMENTATION"
        exit 1
      fi
    if [[ ! -f ${CORTICAL_THICKNESS_GM} ]];
      then
        echo "The cortical gray matter probability image doesn't exist:"
        echo "   $CORTICAL_THICKNESS_GM"
        exit 1
      fi
    if [[ ! -f ${CORTICAL_THICKNESS_WM} ]]
      then
        echo "The cortical white matter probability image doesn't exist:"
        echo "   $CORTICAL_THICKNESS_WM"
        exit 1
      fi

    time_start_direct=`date +%s`

    TMP_FILES=( $CORTICAL_THICKNESS_GM $CORTICAL_THICKNESS_WM $CORTICAL_THICKNESS_SEGMENTATION )

    exe_direct="${DIRECT} -d ${DIMENSION} -s [${CORTICAL_THICKNESS_SEGMENTATION},${GRAY_MATTER_LABEL},${WHITE_MATTER_LABEL}] --verbose 1"
    exe_direct="${exe_direct} -g ${CORTICAL_THICKNESS_GM} -w ${CORTICAL_THICKNESS_WM} -o ${CORTICAL_THICKNESS_IMAGE}"
    exe_direct="${exe_direct} -c ${DIRECT_CONVERGENCE} -r ${DIRECT_GRAD_STEP_SIZE}"
    exe_direct="${exe_direct} -m ${DIRECT_SMOOTHING_PARAMETER} -n ${DIRECT_NUMBER_OF_DIFF_COMPOSITIONS} -b ${USE_BSPLINE_SMOOTHING}"
    if [[ -f ${CORTICAL_LABEL_IMAGE} ]] && [[ -f $REGISTRATION_SUBJECT_WARP ]] && [[ -f $REGISTRATION_SUBJECT_GENERIC_AFFINE ]] ;
      then
        # Calculate ATITH and multiply by a heuristically derived scalar factor
#        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} LabelThickness2 ${CORTICAL_LABEL_IMAGE}
#        logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ${CORTICAL_THICKNESS_GM_SEGMENTATION} 2 2 1 0
#        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} m ${CORTICAL_LABEL_THICKNESS_PRIOR} ${CORTICAL_THICKNESS_GM_SEGMENTATION}
#        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} m ${CORTICAL_LABEL_THICKNESS_PRIOR} 2.0

    	   logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -i ${CORTICAL_LABEL_IMAGE} -o ${CORTICAL_LABEL_THICKNESS_PRIOR} \
	         -t $REGISTRATION_SUBJECT_GENERIC_AFFINE -t $REGISTRATION_SUBJECT_WARP -r ${ANATOMICAL_IMAGES[0]} --verbose 1

        exe_direct="${exe_direct} -a ${CORTICAL_LABEL_THICKNESS_PRIOR}"

        TMP_FILES=( ${TMP_FILES[@]} $CORTICAL_LABEL_THICKNESS_PRIOR ${CORTICAL_THICKNESS_GM_SEGMENTATION} )

      else
        exe_direct="${exe_direct} -t ${DIRECT_THICKNESS_PRIOR}"
      fi
    logCmd $exe_direct

    if [[ $KEEP_TMP_IMAGES -eq 0 ]];
      then
        for f in ${TMP_FILES[@]}
          do
            if [[ -e $f ]];
              then
                logCmd rm $f
              else
                echo "WARNING: expected temp file doesn't exist: $f"
              fi
          done
      fi

    if [[ ! -f ${CORTICAL_THICKNESS_IMAGE} ]];
      then
        echo "Expected output was not produced.  The cortical thickness image doesn't exist:"
        echo "   $CORTICAL_THICKNESS_IMAGE"
        exit 1
      fi

    time_end_direct=`date +%s`
    time_elapsed_direct=$((time_end_direct - time_start_direct))

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with cortical thickness estimation:  $(( time_elapsed_direct / 3600 ))h $(( time_elapsed_direct %3600 / 60 ))m $(( time_elapsed_direct % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo

  fi
  echo ${OUTPUT_PREFIX}ACTStage5Complete.txt > ${OUTPUT_PREFIX}ACTStage5Complete.txt
fi # BAStages
fi # check completion

#### BA Edits Begin ####
if [[ ! -s ${OUTPUT_PREFIX}ACTStage6Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage5Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]] ; then
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 6  ]] ; then # BAStages qc
echo "--------------------------------------------------------------------------------------"
echo "Compute summary measurements"
echo "--------------------------------------------------------------------------------------"
if [[ ! -s ${OUTPUT_PREFIX}CorticalThickness.nii.gz ]] ; then
  echo ${OUTPUT_PREFIX}CorticalThickness.nii.gz incomplete!
  exit 1
fi
if [[ -f ${REGISTRATION_TEMPLATE_WARP} ]];
  then
    exe_template_registration_3="${WARP} -d ${DIMENSION} -i ${CORTICAL_THICKNESS_IMAGE} -o ${OUTPUT_PREFIX}CorticalThicknessNormalizedToTemplate.${OUTPUT_SUFFIX} -r ${REGISTRATION_TEMPLATE} -n Gaussian  -t ${REGISTRATION_TEMPLATE_WARP}  -t ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} --float ${USE_FLOAT_PRECISION} --verbose 1"
    logCmd $exe_template_registration_3

    EXTRACTED_SEGMENTATION_BRAIN_DEFORMED=${OUTPUT_PREFIX}BrainNormalizedToTemplate.${OUTPUT_SUFFIX}

    REGISTRATION_TEMPLATE_BRAIN_MASK=${OUTPUT_PREFIX}RegistrationTemplateBrainMask.nii.gz

    logCmd ${ANTSPATH}/ThresholdImage 3 ${REGISTRATION_TEMPLATE} ${REGISTRATION_TEMPLATE_BRAIN_MASK} 1E-6 Inf

    EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE=${OUTPUT_PREFIX}ExtractedBrain0N4.nii.gz

    logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} m ${HEAD_N4_IMAGE} ${BRAIN_EXTRACTION_MASK}

    TMP_FILES=( ${TMP_FILES[@]} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} ${REGISTRATION_TEMPLATE_BRAIN_MASK} )

    logCmd ${WARP} -d ${DIMENSION} -i ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} -o ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} -r ${REGISTRATION_TEMPLATE} -n Linear -t ${REGISTRATION_TEMPLATE_WARP}  -t ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} --float ${USE_FLOAT_PRECISION} --verbose 1
  fi

ccmetric=`${ANTSPATH}/ImageMath ${DIMENSION} a PearsonCorrelation ${REGISTRATION_TEMPLATE} ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} ${REGISTRATION_TEMPLATE_BRAIN_MASK}`
bvol=`${ANTSPATH}/ImageMath ${DIMENSION} a total ${BRAIN_EXTRACTION_MASK}  | cut -d ':' -f 2 | cut -d ' ' -f 2 `
gvol=`${ANTSPATH}/ImageMath ${DIMENSION} a total ${BRAIN_SEGMENTATION_GM}  | cut -d ':' -f 2 | cut -d ' ' -f 2 `
wvol=`${ANTSPATH}/ImageMath ${DIMENSION} a total ${BRAIN_SEGMENTATION_WM}  | cut -d ':' -f 2 | cut -d ' ' -f 2 `
thks=`${ANTSPATH}/ImageMath ${DIMENSION} a total ${CORTICAL_THICKNESS_IMAGE} | cut -d ':' -f 2 | cut -d ' ' -f 2 `
echo "PearsonCorrelation,BVOL,GVol,WVol,ThicknessSum" >   ${OUTPUT_PREFIX}brainvols.csv
echo "${ccmetric},${bvol},${gvol},${wvol},${thks}" >>  ${OUTPUT_PREFIX}brainvols.csv
if [[ -f ${ANTSPATH}/GetMeshAndTopology ]] && [[ ${DIMENSION} -eq 3 ]] ; then
  ${ANTSPATH}/ThresholdImage ${DIMENSION} ${BRAIN_SEGMENTATION} ${OUTPUT_PREFIX}temp.nii.gz 3 3
  ${ANTSPATH}/ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.nii.gz ME ${OUTPUT_PREFIX}temp.nii.gz 1
  ${ANTSPATH}/ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.nii.gz GetLargestComponent ${OUTPUT_PREFIX}temp.nii.gz 1
  ${ANTSPATH}/ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.nii.gz MD ${OUTPUT_PREFIX}temp.nii.gz 2
  ${ANTSPATH}/SmoothImage 3 ${CORTICAL_THICKNESS_IMAGE} 1 ${OUTPUT_PREFIX}temp2.nii.gz
  #          ${ANTSPATH}/GetMeshAndTopology ${OUTPUT_PREFIX}temp.nii.gz ${OUTPUT_PREFIX}temp2.nii.gz ${OUTPUT_PREFIX}.vtk thickness   0.3 0.001 ${OUTPUT_PREFIX}_Thickness.png
  rm ${OUTPUT_PREFIX}temp.nii.gz ${OUTPUT_PREFIX}temp2.nii.gz
fi
echo "--------------------------------------------------------------------------------------"
#### BA Edits End ####


################################################################################
#
# Create QA/QC output:
#   - tiled mosaic of ${OUTPUT_PREFIX}BrainSegmentation0N4.nii.gz with
#     ${OUTPUT_PREFIX}CorticalThickness.nii.gz overlay
################################################################################

HEAD_N4_IMAGE=${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
HEAD_N4_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentation0N4Resampled.nii.gz"
CORTICAL_THICKNESS_IMAGE_RESAMPLED="${OUTPUT_PREFIX}CorticalThicknessHotResampled.nii.gz"
CORTICAL_THICKNESS_IMAGE_RGB="${OUTPUT_PREFIX}CorticalThicknessHotRGB.nii.gz"
CORTICAL_THICKNESS_MOSAIC="${OUTPUT_PREFIX}CorticalThicknessTiledMosaic.png"
CORTICAL_THICKNESS_MASK="${OUTPUT_PREFIX}CorticalThicknessMask.${OUTPUT_SUFFIX}"
BRAIN_EXTRACTION_MASK_RESAMPLED="${OUTPUT_PREFIX}BrainExtractionMaskResampled.nii.gz"
BRAIN_SEGMENTATION_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentationResampled.nii.gz"
BRAIN_SEGMENTATION_IMAGE_RGB="${OUTPUT_PREFIX}BrainSegmentationRGB.nii.gz"
BRAIN_SEGMENTATION_MOSAIC="${OUTPUT_PREFIX}BrainSegmentationTiledMosaic.png"
ITKSNAP_COLORMAP="${OUTPUT_PREFIX}ItkSnapColormap.txt"

if [[ ! -f ${CORTICAL_THICKNESS_MOSAIC} || ! -f ${BRAIN_SEGMENTATION_MOSAIC} ]];
  then
    TMP_FILES=( $CORTICAL_THICKNESS_IMAGE_RGB $CORTICAL_THICKNESS_MASK $BRAIN_SEGMENTATION_IMAGE_RGB $ITKSNAP_COLORMAP )
    TMP_FILES=( ${TMP_FILES[@]} $HEAD_N4_IMAGE_RESAMPLED $CORTICAL_THICKNESS_IMAGE_RESAMPLED $BRAIN_EXTRACTION_MASK_RESAMPLED $BRAIN_SEGMENTATION_IMAGE_RESAMPLED )

    # Resample images

    resample0="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${BRAIN_EXTRACTION_MASK}"
    resample0="${resample0} ${BRAIN_EXTRACTION_MASK_RESAMPLED} 1 1"
    resample1="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${HEAD_N4_IMAGE}"
    resample1="${resample1} ${HEAD_N4_IMAGE_RESAMPLED} 1 1"
    resample2="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE}"
    resample2="${resample2} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED} 1 1"
    resample3="${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${BRAIN_SEGMENTATION}"
    resample3="${resample3} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED} 1 1"
    if [[ ${DIMENSION} -eq 3 ]];
      then
        resample0="${resample0} 1 0 0 1"
        resample1="${resample1} 1"
        resample2="${resample2} 1"
        resample3="${resample3} 1 0 0 1"
      else
        resample0="${resample0} 0 0 1"
        resample3="${resample3} 0 0 1"
      fi
    logCmd $resample0
    logCmd $resample1
    logCmd $resample2
    logCmd $resample3

    # Cortical thickness

    mask="${ANTSPATH}/ThresholdImage ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED} ${CORTICAL_THICKNESS_MASK} 0 0 0 1"
    logCmd $mask

    conversion="${ANTSPATH}/ConvertScalarImageToRGB ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED}"
    conversion="${conversion} ${CORTICAL_THICKNESS_IMAGE_RGB} none hot none 0 ${DIRECT_THICKNESS_PRIOR}"
    logCmd $conversion

    mosaic="${ANTSPATH}/CreateTiledMosaic -i ${HEAD_N4_IMAGE_RESAMPLED} -r ${CORTICAL_THICKNESS_IMAGE_RGB}"
    mosaic="${mosaic} -o ${CORTICAL_THICKNESS_MOSAIC} -a 1.0 -t -1x-1 -d z -p mask"
    mosaic="${mosaic} -s [2,mask,mask] -x ${CORTICAL_THICKNESS_MASK}"
    logCmd $mosaic

    # Segmentation

    echo "0 1 0 0 1 0 1" > $ITKSNAP_COLORMAP
    echo "0 0 1 0 1 1 0" >> $ITKSNAP_COLORMAP
    echo "0 0 0 1 0 1 1" >> $ITKSNAP_COLORMAP

    conversion="${ANTSPATH}/ConvertScalarImageToRGB ${DIMENSION} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED}"
    conversion="${conversion} ${BRAIN_SEGMENTATION_IMAGE_RGB} none custom $ITKSNAP_COLORMAP 0 6"
    logCmd $conversion

    mosaic="${ANTSPATH}/CreateTiledMosaic -i ${HEAD_N4_IMAGE_RESAMPLED} -r ${BRAIN_SEGMENTATION_IMAGE_RGB}"
    mosaic="${mosaic} -o ${BRAIN_SEGMENTATION_MOSAIC} -a 0.3 -t -1x-1 -d 2 -p mask"
    mosaic="${mosaic} -s [2,mask,mask] -x ${BRAIN_EXTRACTION_MASK_RESAMPLED}"
    logCmd $mosaic


    if [[ $KEEP_TMP_IMAGES -eq 0 ]];
      then
        for f in ${TMP_FILES[@]}
          do
            if [[ -e $f ]];
             then
              logCmd rm $f
            else
              echo "WARNING: expected temp file doesn't exist: $f"
            fi
        done
      fi
  fi
  echo ${OUTPUT_PREFIX}ACTStage6Complete.txt > ${OUTPUT_PREFIX}ACTStage6Complete.txt
fi # BAStages
fi # check completion
################################################################################
#
# End of main routine
#
################################################################################

if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -ge 5  ]] ; then
  logCmd checkOutputExists
fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with ANTs processing pipeline: stage $ACT_STAGE"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
