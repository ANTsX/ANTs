#!/bin/bash

VERSION="0.0"

# Check dependencies

PROGRAM_DEPENDENCIES=( 'antsRegistration' 'antsApplyTransforms' 'N4BiasFieldCorrection' 'Atropos' 'KellyKapowski' )
SCRIPTS_DEPENDENCIES=( 'antsBrainExtraction.sh' 'antsAtroposN4.sh' )

for D in ${PROGRAM_DEPENDENCIES[@]};
  do
    if ! command -v ${D} &> /dev/null
      then
        echo "Error:  we can't find the $D program."
        echo "Perhaps you need to \(re\)define \$PATH in your environment."
        exit
      fi
  done

for D in ${SCRIPT_DEPENDENCIES[@]};
  do
    if ! command -v ${D} &> /dev/null
      then
        echo "We can't find the $D script."
        echo "Perhaps you need to \(re\)define \$PATH in your environment."
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
                                                We recommend using the T1 as the first image.

     -e:  Brain segmentation template           Anatomical *intensity* template. This template is *not* skull-stripped.
                                                The following images must be in the same space as this template:
                                                    * Brain probability mask (-m)
                                                    * Segmentation priors (-p).
                                                If used, the following optional images must also be in the same space as
                                                this template:
                                                    * Registration metric mask (-f)
                                                    * Thickness prior image (-r).

     -m:  Brain extraction probability mask     Brain *probability* mask in the segmentation template space. A binary mask
                                                is an acceptable probability image.

     -p:  Brain segmentation priors             Tissue *probability* priors corresponding to the template image specified
                                                with the -e option.  Specified using c-style formatting, either with numeric indices
                                                e.g.
                                                  -p template/priors/Priors%02d.nii.gz
                                                or BIDS style, eg
                                                  -p tpl-templateName/tpl-templateName_res-01_label-%s_probseg.nii.gz

                                                At least four priors must exist, corresponding to CSF, Cortical GM, WM, Subcortical GM.
                                                If priors are numbered numerically, the classes must be ordered in the same way, ie
                                                  1:  CSF
                                                  2:  Cortical GM
                                                  3:  WM
                                                  4:  Subcortical GM

                                                In BIDS format, the labels must include CSF, CGM, WM, SGM, and may optionally include BS, CBM.
                                                Other labels will not be used. Templateflow labels use "SCGM" instead of "SGM", if "SGM" is not
                                                found, "SCGM" will be used instead.

     -o:  Output prefix                         A partial list of output images:
                                                  * ${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentation*N4.${OUTPUT_SUFFIX} One for each anatomical input
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*1.${OUTPUT_SUFFIX}  CSF
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*2.${OUTPUT_SUFFIX}  Cortical GM
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*3.${OUTPUT_SUFFIX}  WM
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*4.${OUTPUT_SUFFIX}  Subcortical GM
                                                  * ... and so on for additional segmentation classes
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors*N.${OUTPUT_SUFFIX} where there are N priors
                                                  *                              Number formatting of posteriors matches that of the priors.
                                                  * ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}
                                                More information on the output can be found on the ANTs Wiki
                                                https://github.com/ANTsX/ANTs/wiki.

Optional arguments:

     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd.

     -t:  template for t1 registration          Anatomical *intensity* template. This template *must* be skull-stripped.
                                                This template is used to produce a final, high-quality registration between
                                                the bias-corrected, skull-stripped subject anatomical image and the template.
                                                This template will commonly be a skull-stripped version of the template passed
                                                with -e. We perform the registration with fixed image = (this template)
                                                and moving image = (input anatomical image).
                                                The output from this step is
                                                  * Forward warps:
                                                    - ${OUTPUT_PREFIX}SubjectToTemplate1Warp.nii.gz
                                                    - ${OUTPUT_PREFIX}SubjectToTemplate0GenericAffine.mat
                                                  * Inverse warps:
                                                    - ${OUTPUT_PREFIX}TemplateToSubject1GenericAffine.mat
                                                    - ${OUTPUT_PREFIX}TemplateToSubject0Warp.nii.gz
                                                  * Jacobian:
                                                    - ${OUTPUT_PREFIX}SubjectToTemplateLogJacobian.${OUTPUT_SUFFIX}

                                                More information on the how to use these images can be found on the ANTs Wiki
                                                https://github.com/ANTsX/ANTs/wiki.

     -f:  extraction registration mask          Binary metric mask defined in the segmentation template space (-e). During the
                                                registration for brain extraction, the similarity metric is only computed within
                                                this mask.

     -k:  keep temporary files                  Keep brain extraction/segmentation warps, etc (default = 0).

     -g:  denoise anatomical images             Denoise anatomical images (default = 0).

     -i:  max iterations for registration       ANTS registration max iterations (default = 100x100x70x20).

     -w:  Atropos prior segmentation weight     Atropos spatial prior *probability* weight for the segmentation (default = 0.25).

     -n:  number of segmentation iterations     N4 -> Atropos -> N4 iterations during segmentation (default = 3).

     -b:  posterior formulation                 Atropos posterior formulation and whether or not to use mixture model proportions.
                                                e.g 'Socrates[ 1 ]' (default) or 'Aristotle[ 1 ]'.  Choose the latter if you
                                                want use the distance priors (see also the -l option for label propagation
                                                control).

     -j:  use floating-point precision          Use single float precision in registrations (default = 0).

     -u:  use random seeding                    Use random number generated from system clock (default = 1). If 0, a fixed
                                                random seed is used. To set your own seed, set this option to 0 and export the
                                                environment variable ANTS_RANDOM_SEED. To achieve exactly identical results, you must
                                                also set ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS to 1.

     -v:  use b-spline smoothing                Use B-spline SyN for registrations and B-spline exponential mapping in DiReCT (default = 0).

     -r:  cortical thickness prior image        Cortical thickness prior image in the template space, which contains an estimated
                                                upper limit to the cortical thickness at each voxel. If not specified, the prior is
                                                set to 10mm throughout the brain.

     -l:  label propagation                     Incorporate a distance prior one the posterior formulation.  Should be
                                                of the form 'label[ lambda,boundaryProbability ]' where label is a value
                                                of 1,2,3,... denoting label ID.  The label probability for anything
                                                outside the current label

                                                  = boundaryProbability * exp( -lambda * distanceFromBoundary )

                                                Intuitively, smaller lambda values will increase the spatial capture
                                                range of the distance prior.  To apply to all label values, simply omit
                                                specifying the label, i.e. '-l "[ lambda,boundaryProbability ]"'.

     -c:  Additional priors for thickness       Add segmentation classes to be treated as gray or white matter for thickness
                                                estimation. For example, when calling KellyKapowski for normal subjects, we
                                                combine the deep gray matter segmentation/posteriors (class 4) with the white
                                                matter segmentation/posteriors (class 3).
                                                Another example would be computing cortical thickness in the presence
                                                of white matter lesions. We can accommodate this by specifying a lesion mask
                                                posterior as an additional posterior (suppose label '7'), combining this with
                                                normal white matter in the thickness estimation by specifying '-c "WM[ 7 ]"'
                                                or '-c "3[ 7 ]"'.

     -q:  Use quick registration parameters     If = 1, use antsRegistrationSyNQuick.sh as the basis for registration
                                                during brain extraction, brain segmentation, and (optional) normalization
                                                to a template.  Otherwise use a slower registration comparable to
                                                antsRegistrationSyN.sh (default = 0).

     -x:  Atropos iterations                    Number of iterations within Atropos (default 5).

     -y:  Script stage to run                   Which stage of ACT to run (default = 0, run all).  Tries to split in 2 hour chunks.
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

  for (( j = 1; j <= ${NUMBER_OF_PRIOR_IMAGES}; j++ ));
    do
      num=$(printf "${WARPED_PRIOR_FORMAT}" $j)

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
      metric                  = ${ANTS_METRIC}[ fixedImage,movingImage,${ANTS_METRIC_PARAMS} ]
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

ANTS=antsRegistration
ANTS_MAX_ITERATIONS="100x100x70x20"
ANTS_TRANSFORMATION="SyN[ 0.2,3,0 ]"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[ 1000x500x250x100,1e-8,10 ]"
ANTS_METRIC="CC"
ANTS_METRIC_PARAMS="1,4"

WARP=antsApplyTransforms

N4=N4BiasFieldCorrection
N4_CONVERGENCE_1="[ 50x50x50x50,0.0000001 ]"
N4_CONVERGENCE_2="[ 50x50x50x50,0.0000001 ]"
N4_SHRINK_FACTOR_1=4
N4_SHRINK_FACTOR_2=2
N4_BSPLINE_PARAMS="[ 200 ]"

ATROPOS=Atropos

ATROPOS_SEGMENTATION_INITIALIZATION="PriorProbabilityImages"
ATROPOS_SEGMENTATION_LIKELIHOOD="Gaussian"
ATROPOS_SEGMENTATION_CONVERGENCE="[ 5,0.0 ]"
ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION="Socrates[ 1 ]"
ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=3
ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5 # to be backward compatible but i like 25
ATROPOS_SEGMENTATION_LABEL_PROPAGATION=()

DIRECT=KellyKapowski
DIRECT_CONVERGENCE="[ 45,0.0,10 ]"
DIRECT_THICKNESS_PRIOR="10"
DIRECT_GRAD_STEP_SIZE="0.025"
DIRECT_SMOOTHING_PARAMETER="1.5"
DIRECT_NUMBER_OF_DIFF_COMPOSITIONS="10"

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
    ANTS_TRANSFORMATION="BSplineSyN[ 0.2,26,0,3 ]"
    DIRECT_SMOOTHING_PARAMETER="5.75"
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

    DIRECT_CONVERGENCE="[ 5,0.0,10 ]"

    # Fix random seed to replicate exact results on each run
    USE_RANDOM_SEEDING=0

  fi

if [[ ${USE_RANDOM_SEEDING} -eq 0 ]]; then
  # Use fixed random seed unless one is already defined
  if [[ -z $ANTS_RANDOM_SEED ]] ; then
    export ANTS_RANDOM_SEED=19650218
  fi
fi

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#  3. Check template priors and define numbering for warped priors
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

# SEGMENTATION_PRIOR is a c-style string containing the path to priors, eg /path/to/priors%d.nii.gz
# Allow either numeric formatting eg priors%02d.nii.gz or BIDS string formatting using BIDS common
# derived labels eg tpl-templateName_res-01_label-%s_probseg.nii.gz or templateflow labels (same
# except for SGM).
#
# Warped priors and posteriors are always formatted numerically, eg warpedPriors%02d.nii.gz.
#
# If the template priors use a different format, eg %0d or %03d, the warped priors and posteriors
# will be formatted accordingly.


# List of file names in template space, in order: CSF, GM, WM, deep GM[, BS, CBM, others]
PRIOR_IMAGE_FILENAMES=()

# Default, overwrite with input numeric format if used
WARPED_PRIOR_FORMAT="%02d"

# Test if the pattern contains "%d" or "%s"
if [[ "$SEGMENTATION_PRIOR" =~ %([0-9]*)d ]]; then
  # Numeric formatting - normally 1 through 6, though only the first four are mandatory
  WARPED_PRIOR_FORMAT="%${BASH_REMATCH[1]}d"

  PRIOR_COUNTER=1
  FOUND_PRIOR=1
  while [[ $FOUND_PRIOR -gt 0 ]]; do
    PRIOR_FILE=$(printf "$SEGMENTATION_PRIOR" "$PRIOR_COUNTER")
    if [[ -f "$PRIOR_FILE" ]]; then
      PRIOR_IMAGE_FILENAMES+=("$PRIOR_FILE")
    else
        FOUND_PRIOR=0
    fi
    PRIOR_COUNTER=$((PRIOR_COUNTER+1))
  done
elif [[ "$SEGMENTATION_PRIOR" =~ label-%s ]]; then

  # BIDS label-%s formatting
  # Order matters here - they will be added to an array in order, required labels first and then
  # optional. Order must match ants expectations, ie CSF, GM, WM, deep GM [, BS, CBM]
  REQUIRED_LABELS=("CSF" "CGM" "WM" "SGM")
  OPTIONAL_LABELS=("BS" "CBM")

  # Check for the existence of required labels
  for BIDS_LABEL in "${REQUIRED_LABELS[@]}"; do
    PRIOR_FILE=$(printf "$SEGMENTATION_PRIOR" "$BIDS_LABEL")
    if [[ -f "$PRIOR_FILE" ]]; then
      PRIOR_IMAGE_FILENAMES+=("$PRIOR_FILE")
    elif [[ "$BIDS_LABEL" == "SGM" ]]; then
      # Hack for templateflow - try SCGM instead of SGM
      # Unfortunately templateflow does not use BIDS standard labels for SGM
      PRIOR_FILE=$(printf "$SEGMENTATION_PRIOR" "SCGM")
        if [[ -f "$PRIOR_FILE" ]]; then
          PRIOR_IMAGE_FILENAMES+=("$PRIOR_FILE")
        else
          echo "Error: Required segmentation label "SGM" not found (also tried templateflow "SCGM")"
          exit 1
        fi
    else
      echo "Error: Required segmentation label $BIDS_LABEL not found."
      exit 1
    fi
  done

  # Check for the existence of optional labels
  for BIDS_LABEL in "${OPTIONAL_LABELS[@]}"; do
    PRIOR_FILE=$(printf "$SEGMENTATION_PRIOR" "$BIDS_LABEL")
    if [[ -f "$PRIOR_FILE" ]]; then
      PRIOR_IMAGE_FILENAMES+=("$PRIOR_FILE")
    fi
  done
else
  echo "Error: Invalid prior specification $SEGMENTATION_PRIOR"
  exit 1
fi

NUMBER_OF_PRIOR_IMAGES=${#PRIOR_IMAGE_FILENAMES[*]}

if [[ ${NUMBER_OF_PRIOR_IMAGES} -lt 4 ]]; then
  echo "Expected at least 4 prior images (${NUMBER_OF_PRIOR_IMAGES} are specified).  Check the command line specification."
  exit 1
fi

BRAIN_SEGMENTATION_OUTPUT_PREFIX=${OUTPUT_PREFIX}BrainSegmentation
SEGMENTATION_WARP_OUTPUT_PREFIX=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Prior
SEGMENTATION_PRIOR_WARPED_PREFIX=${SEGMENTATION_WARP_OUTPUT_PREFIX}Warped

SEGMENTATION_PRIOR_WARPED_FORMAT=${SEGMENTATION_PRIOR_WARPED_PREFIX}${WARPED_PRIOR_FORMAT}.${OUTPUT_SUFFIX}

# Whatever the input spec for priors, we always number warped priors sequentially as required by Atropos
WARPED_PRIOR_IMAGE_FILENAMES=()

for (( i = 0; i < ${NUMBER_OF_PRIOR_IMAGES}; i++ )); do
    WARPED_PRIOR_IMAGE_FILENAMES[${#WARPED_PRIOR_IMAGE_FILENAMES[@]}]=$( printf ${SEGMENTATION_PRIOR_WARPED_FORMAT} $((i+1)) )
done

# These arrays contain the formatted labels that will be used to combine with the white
# and gray matter posteriors for the cortical thickness section.

CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS=()
CORTICAL_THICKNESS_GRAY_MATTER_OTHER_LABELS=()

for(( j=0; j < ${#PRIOR_COMBINATIONS[@]}; j++ ))
  do
    echo ${PRIOR_COMBINATIONS[$j]}
    COMBINATION=( $( echo ${PRIOR_COMBINATIONS[$j]} | tr -d ' ' | tr '[],' '\n' ) )

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

if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 1 ]];
  then
    if [[ ! -f ${REGISTRATION_TEMPLATE} ]]
      then
        echo "Template for registration, ${REGISTRATION_TEMPLATE}, does not exist."
        exit 1
      fi
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
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 1 ]] ; then # BAStages bxt
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
        logCmd antsBrainExtraction.sh \
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
        logCmd antsBrainExtraction.sh \
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
    logCmd ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN} m ${ANATOMICAL_IMAGES[0]} ${BRAIN_EXTRACTION_MASK}

    # Do a quick N4 on the brain before registration
    logCmd $N4 -d ${DIMENSION} -i ${EXTRACTED_SEGMENTATION_BRAIN} -s ${N4_SHRINK_FACTOR_1} -c ${N4_CONVERGENCE_1} -o ${EXTRACTED_SEGMENTATION_BRAIN} -x ${BRAIN_EXTRACTION_MASK} -b ${N4_BSPLINE_PARAMS}

  fi

if [[ -f ${BRAIN_TEMPLATE} ]] && [[ ! -f ${EXTRACTED_BRAIN_TEMPLATE} ]];
  then
    logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_PRIOR} ${EXTRACTED_BRAIN_TEMPLATE} 0.1 1.01 1 0
    logCmd ImageMath ${DIMENSION} ${EXTRACTED_BRAIN_TEMPLATE} m ${BRAIN_TEMPLATE} ${EXTRACTED_BRAIN_TEMPLATE}
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
SEGMENTATION_MASK_DILATED=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}MaskDilated.${OUTPUT_SUFFIX}
SEGMENTATION_CONVERGENCE_FILE=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Convergence.txt

if [[ ! -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]]  && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]]; then
  if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 2 ]] ; then # BAStages reg
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
          rm -f ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${EXTRACTED_SEGMENTATION_BRAIN} ]];
        then
          echo "The extracted brain doesn't exist:"
          echo "   ${EXTRACTED_SEGMENTATION_BRAIN}"
          rm -f ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi
      if [[ ! -f ${BRAIN_EXTRACTION_MASK} ]];
        then
          echo "The brain extraction mask does not exist:"
          echo "   ${BRAIN_EXTRACTION_MASK}"
          rm -f ${OUTPUT_PREFIX}ACTStage1Complete.txt
          exit 1
        fi

        if [[ ! -f ${SEGMENTATION_WARP} ]];
          then
            logCmd ImageMath ${DIMENSION} ${SEGMENTATION_MASK_DILATED} MD ${BRAIN_EXTRACTION_MASK} 20

            basecall=''
            if [[ ${RUN_QUICK} -ne 0 ]];
              then
                TMP_FILES=( ${TMP_FILES[@]} "${SEGMENTATION_WARP_OUTPUT_PREFIX}Warped.nii.gz" "${SEGMENTATION_WARP_OUTPUT_PREFIX}InverseWarped.nii.gz" )

                basecall="antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${EXTRACTED_SEGMENTATION_BRAIN}"
                basecall="${basecall} -m ${EXTRACTED_BRAIN_TEMPLATE} -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} -j 1"
                if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
                  then
                    basecall="${basecall} -p f"
                  fi
            else
              basecall="${ANTS} -d ${DIMENSION} -u 0 -w [ 0.0,0.999 ] -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} --float ${USE_FLOAT_PRECISION} --verbose 1"
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
                logCmd CopyImageHeaderInformation ${ANATOMICAL_IMAGES[0]} ${img} ${img} 1 1 1
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
      echo " Done with brain registration:  $(( time_elapsed_brain_registration / 3600 ))h $(( time_elapsed_brain_registration %3600 / 60 ))m $(( time_elapsed_brain_registration % 60 ))s"
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
  if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 3 ]] ; then # BAStages seg
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

    logCmd antsAtroposN4.sh \
      -d ${DIMENSION} \
      -b "${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}" \
      ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} \
      ${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} \
      -x ${BRAIN_EXTRACTION_MASK} \
      -m ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS} \
      -n ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
      -c ${NUMBER_OF_PRIOR_IMAGES} \
      ${N4_INCLUDE_PRIORS_COMMAND_LINE} \
      -p ${SEGMENTATION_PRIOR_WARPED_FORMAT} \
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
    logCmd antsAtroposN4.sh \
      -d ${DIMENSION} \
      -b "${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}" \
      ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} \
      ${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} \
      -x ${BRAIN_EXTRACTION_MASK} \
      -m 2 \
      -n ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
      -c ${NUMBER_OF_PRIOR_IMAGES} \
      ${N4_INCLUDE_PRIORS_COMMAND_LINE} \
      -p ${SEGMENTATION_PRIOR_WARPED_FORMAT} \
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

    logCmd ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} m ${HEAD_N4_IMAGE} ${BRAIN_EXTRACTION_MASK}

    TMP_FILES=( ${TMP_FILES[@]} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} )

    time_start_template_registration=`date +%s`

    basecall=''
    if [[ ${RUN_QUICK} -ne 0 ]];
      then
        TMP_FILES=( ${TMP_FILES[@]} "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}Warped.nii.gz" "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}InverseWarped.nii.gz" )

        basecall="antsRegistrationSyNQuick.sh -d ${DIMENSION} -f ${REGISTRATION_TEMPLATE}"
        basecall="${basecall} -m ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -j 1"
        if [[ ${USE_FLOAT_PRECISION} -ne 0 ]];
          then
            basecall="${basecall} -p f"
          fi
      else
        IMAGES="${REGISTRATION_TEMPLATE},${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE}"
        basecall="${ANTS} -d ${DIMENSION} -v 1 -u 0 -w [ 0.0,0.999 ] -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -r [ ${IMAGES},1 ] --float ${USE_FLOAT_PRECISION}"
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
    logCmd antsApplyTransforms -d ${DIMENSION} -o Linear[ $REGISTRATION_SUBJECT_GENERIC_AFFINE,1 ] -t $REGISTRATION_TEMPLATE_GENERIC_AFFINE --verbose 1

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
  logCmd CreateJacobianDeterminantImage ${DIMENSION} ${REGISTRATION_TEMPLATE_WARP} ${REGISTRATION_LOG_JACOBIAN} 1 1
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

WARPED_GRAY_MATTER_LABEL=$( printf "${WARPED_PRIOR_FORMAT}" ${GRAY_MATTER_LABEL} )
WARPED_WHITE_MATTER_LABEL=$( printf "${WARPED_PRIOR_FORMAT}" ${WHITE_MATTER_LABEL} )
WARPED_DEEP_GRAY_MATTER_LABEL=$( printf "${WARPED_PRIOR_FORMAT}" ${DEEP_GRAY_MATTER_LABEL} )

BRAIN_SEGMENTATION_GM=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Posteriors${WARPED_GRAY_MATTER_LABEL}.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION_WM=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Posteriors${WARPED_WHITE_MATTER_LABEL}.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION_DEEP_GM=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Posteriors${WARPED_DEEP_GRAY_MATTER_LABEL}.${OUTPUT_SUFFIX}

CORTICAL_THICKNESS_GM=${OUTPUT_PREFIX}CorticalThicknessPosteriors${WARPED_GRAY_MATTER_LABEL}.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_WM=${OUTPUT_PREFIX}CorticalThicknessPosteriors${WARPED_WHITE_MATTER_LABEL}.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_SEGMENTATION=${OUTPUT_PREFIX}CorticalThicknessSegmentation.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_GM_SEGMENTATION=${OUTPUT_PREFIX}CorticalThicknessGMSegmentation.${OUTPUT_SUFFIX}
CORTICAL_LABEL_THICKNESS_PRIOR=${OUTPUT_PREFIX}CorticalLabelThicknessPrior.${OUTPUT_SUFFIX}

if [[ ! -s ${OUTPUT_PREFIX}ACTStage5Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage3Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage2Complete.txt ]] && \
   [[   -s ${OUTPUT_PREFIX}ACTStage1Complete.txt ]] ; then
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 5 ]] ; then # BAStages thk
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
        OTHER_LABEL_FORMATTED=$( printf "${WARPED_PRIOR_FORMAT}" ${OTHER_LABEL} )
        BRAIN_SEGMENTATION_OTHER_LABEL=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Posteriors${OTHER_LABEL_FORMATTED}.${OUTPUT_SUFFIX}

        logCmd ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_GM} + ${CORTICAL_THICKNESS_GM} ${BRAIN_SEGMENTATION_OTHER_LABEL}
        logCmd ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ReplaceVoxelValue ${BRAIN_SEGMENTATION} ${OTHER_LABEL} ${OTHER_LABEL} ${GRAY_MATTER_LABEL}
      done

    logCmd cp ${BRAIN_SEGMENTATION_WM} ${CORTICAL_THICKNESS_WM}
    for(( i=0; i < ${#CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[@]}; i++ ))
      do
        OTHER_LABEL=${CORTICAL_THICKNESS_WHITE_MATTER_OTHER_LABELS[$i]}
        OTHER_LABEL_FORMATTED=$( printf "${WARPED_PRIOR_FORMAT}" ${OTHER_LABEL} )
        BRAIN_SEGMENTATION_OTHER_LABEL=${BRAIN_SEGMENTATION_OUTPUT_PREFIX}Posteriors${OTHER_LABEL_FORMATTED}.${OUTPUT_SUFFIX}

        logCmd ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_WM} + ${CORTICAL_THICKNESS_WM} ${BRAIN_SEGMENTATION_OTHER_LABEL}
        logCmd ImageMath ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ReplaceVoxelValue ${BRAIN_SEGMENTATION} ${OTHER_LABEL} ${OTHER_LABEL} ${WHITE_MATTER_LABEL}
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

    exe_direct="${DIRECT} -d ${DIMENSION} -s [ ${CORTICAL_THICKNESS_SEGMENTATION},${GRAY_MATTER_LABEL},${WHITE_MATTER_LABEL} ] --verbose 1"
    exe_direct="${exe_direct} -g ${CORTICAL_THICKNESS_GM} -w ${CORTICAL_THICKNESS_WM} -o ${CORTICAL_THICKNESS_IMAGE}"
    exe_direct="${exe_direct} -c ${DIRECT_CONVERGENCE} -r ${DIRECT_GRAD_STEP_SIZE}"
    exe_direct="${exe_direct} -m ${DIRECT_SMOOTHING_PARAMETER} -n ${DIRECT_NUMBER_OF_DIFF_COMPOSITIONS} -b ${USE_BSPLINE_SMOOTHING}"
    if [[ -f ${CORTICAL_LABEL_IMAGE} ]] && [[ -f $REGISTRATION_SUBJECT_WARP ]] && [[ -f $REGISTRATION_SUBJECT_GENERIC_AFFINE ]] ;
      then
        # Calculate ATITH and multiply by a heuristically derived scalar factor
#        logCmd ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} LabelThickness2 ${CORTICAL_LABEL_IMAGE}
#        logCmd ThresholdImage ${DIMENSION} ${CORTICAL_THICKNESS_SEGMENTATION} ${CORTICAL_THICKNESS_GM_SEGMENTATION} 2 2 1 0
#        logCmd ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} m ${CORTICAL_LABEL_THICKNESS_PRIOR} ${CORTICAL_THICKNESS_GM_SEGMENTATION}
#        logCmd ImageMath ${DIMENSION} ${CORTICAL_LABEL_THICKNESS_PRIOR} m ${CORTICAL_LABEL_THICKNESS_PRIOR} 2.0

    	   logCmd antsApplyTransforms -d ${DIMENSION} -i ${CORTICAL_LABEL_IMAGE} -o ${CORTICAL_LABEL_THICKNESS_PRIOR} \
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
                logCmd rm -f $f
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
if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -eq 6 ]] ; then # BAStages qc
echo "--------------------------------------------------------------------------------------"
echo "Compute summary measurements"
echo "--------------------------------------------------------------------------------------"
if [[ ! -s ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX} ]] ; then
  echo ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX} incomplete!
  exit 1
fi
if [[ -f ${REGISTRATION_TEMPLATE_WARP} ]];
  then
    exe_template_registration_3="${WARP} -d ${DIMENSION} -i ${CORTICAL_THICKNESS_IMAGE} -o ${OUTPUT_PREFIX}CorticalThicknessNormalizedToTemplate.${OUTPUT_SUFFIX} -r ${REGISTRATION_TEMPLATE} -n Gaussian  -t ${REGISTRATION_TEMPLATE_WARP}  -t ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} --float ${USE_FLOAT_PRECISION} --verbose 1"
    logCmd $exe_template_registration_3

    EXTRACTED_SEGMENTATION_BRAIN_DEFORMED=${OUTPUT_PREFIX}BrainNormalizedToTemplate.${OUTPUT_SUFFIX}

    REGISTRATION_TEMPLATE_BRAIN_MASK=${OUTPUT_PREFIX}RegistrationTemplateBrainMask.${OUTPUT_SUFFIX}

    logCmd ThresholdImage 3 ${REGISTRATION_TEMPLATE} ${REGISTRATION_TEMPLATE_BRAIN_MASK} 1E-6 Inf

    EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE=${OUTPUT_PREFIX}ExtractedBrain0N4.${OUTPUT_SUFFIX}

    logCmd ImageMath ${DIMENSION} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} m ${HEAD_N4_IMAGE} ${BRAIN_EXTRACTION_MASK}

    TMP_FILES=( ${TMP_FILES[@]} ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} ${REGISTRATION_TEMPLATE_BRAIN_MASK} )

    logCmd ${WARP} -d ${DIMENSION} -i ${EXTRACTED_SEGMENTATION_BRAIN_N4_IMAGE} -o ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} -r ${REGISTRATION_TEMPLATE} -n Linear -t ${REGISTRATION_TEMPLATE_WARP}  -t ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} --float ${USE_FLOAT_PRECISION} --verbose 1
  fi

ccmetric=`ImageMath ${DIMENSION} a PearsonCorrelation ${REGISTRATION_TEMPLATE} ${EXTRACTED_SEGMENTATION_BRAIN_DEFORMED} ${REGISTRATION_TEMPLATE_BRAIN_MASK}`
bvol=`ImageMath ${DIMENSION} a total ${BRAIN_EXTRACTION_MASK}  | cut -d ':' -f 3 | cut -d ' ' -f 2 `
gvol=`ImageMath ${DIMENSION} a total ${BRAIN_SEGMENTATION_GM}  | cut -d ':' -f 3 | cut -d ' ' -f 2 `
wvol=`ImageMath ${DIMENSION} a total ${BRAIN_SEGMENTATION_WM}  | cut -d ':' -f 3 | cut -d ' ' -f 2 `
thks=`ImageMath ${DIMENSION} a total ${CORTICAL_THICKNESS_IMAGE} | cut -d ':' -f 3 | cut -d ' ' -f 2 `
echo "PearsonCorrelation,BVOL,GVol,WVol,ThicknessSum" >   ${OUTPUT_PREFIX}brainvols.csv
echo "${ccmetric},${bvol},${gvol},${wvol},${thks}" >>  ${OUTPUT_PREFIX}brainvols.csv
if [[ -f GetMeshAndTopology ]] && [[ ${DIMENSION} -eq 3 ]] ; then
  ThresholdImage ${DIMENSION} ${BRAIN_SEGMENTATION} ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} 3 3
  ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} ME ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} 1
  ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} GetLargestComponent ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} 1
  ImageMath ${DIMENSION} ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} MD ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} 2
  SmoothImage 3 ${CORTICAL_THICKNESS_IMAGE} 1 ${OUTPUT_PREFIX}temp2.${OUTPUT_SUFFIX}
  #          GetMeshAndTopology ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}temp2.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}.vtk thickness   0.3 0.001 ${OUTPUT_PREFIX}_Thickness.png
  rm -f ${OUTPUT_PREFIX}temp.${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}temp2.${OUTPUT_SUFFIX}
fi
echo "--------------------------------------------------------------------------------------"
#### BA Edits End ####


################################################################################
#
# Create QA/QC output:
#   - tiled mosaic of ${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX} with
#     ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX} overlay
################################################################################

HEAD_N4_IMAGE=${OUTPUT_PREFIX}BrainSegmentation0N4.${OUTPUT_SUFFIX}
HEAD_N4_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentation0N4Resampled.${OUTPUT_SUFFIX}"
CORTICAL_THICKNESS_IMAGE_RESAMPLED="${OUTPUT_PREFIX}CorticalThicknessHotResampled.${OUTPUT_SUFFIX}"
CORTICAL_THICKNESS_IMAGE_RGB="${OUTPUT_PREFIX}CorticalThicknessHotRGB.${OUTPUT_SUFFIX}"
CORTICAL_THICKNESS_MOSAIC="${OUTPUT_PREFIX}CorticalThicknessTiledMosaic.png"
CORTICAL_THICKNESS_MASK="${OUTPUT_PREFIX}CorticalThicknessMask.${OUTPUT_SUFFIX}"
BRAIN_EXTRACTION_MASK_RESAMPLED="${OUTPUT_PREFIX}BrainExtractionMaskResampled.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_IMAGE_RESAMPLED="${OUTPUT_PREFIX}BrainSegmentationResampled.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_IMAGE_RGB="${OUTPUT_PREFIX}BrainSegmentationRGB.${OUTPUT_SUFFIX}"
BRAIN_SEGMENTATION_MOSAIC="${OUTPUT_PREFIX}BrainSegmentationTiledMosaic.png"
ITKSNAP_COLORMAP="${OUTPUT_PREFIX}ItkSnapColormap.txt"

if [[ ! -f ${CORTICAL_THICKNESS_MOSAIC} || ! -f ${BRAIN_SEGMENTATION_MOSAIC} ]];
  then
    TMP_FILES=( $CORTICAL_THICKNESS_IMAGE_RGB $CORTICAL_THICKNESS_MASK $BRAIN_SEGMENTATION_IMAGE_RGB $ITKSNAP_COLORMAP )
    TMP_FILES=( ${TMP_FILES[@]} $HEAD_N4_IMAGE_RESAMPLED $CORTICAL_THICKNESS_IMAGE_RESAMPLED $BRAIN_EXTRACTION_MASK_RESAMPLED $BRAIN_SEGMENTATION_IMAGE_RESAMPLED )

    # Resample images

    resample0="ResampleImage ${DIMENSION} ${BRAIN_EXTRACTION_MASK} ${BRAIN_EXTRACTION_MASK_RESAMPLED}"
    resample1="ResampleImage ${DIMENSION} ${HEAD_N4_IMAGE} ${HEAD_N4_IMAGE_RESAMPLED}"
    resample2="ResampleImage ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED}"
    resample3="ResampleImage ${DIMENSION} ${BRAIN_SEGMENTATION} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED}"

    if [[ ${DIMENSION} -eq 3 ]];
      then
        resample0="${resample0} 1x1x1 0 1 6"
        resample1="${resample1} 1x1x1 0 0 6"
        resample2="${resample2} 1x1x1 0 0 6"
        resample3="${resample3} 1x1x1 0 1 6"
      else
        resample0="${resample0} 1x1 0 1 6"
        resample1="${resample1} 1x1 0 0 6"
        resample2="${resample2} 1x1 0 0 6"
        resample3="${resample3} 1x1 0 1 6"
      fi
    logCmd $resample0
    logCmd $resample1
    logCmd $resample2
    logCmd $resample3

    # Cortical thickness

    mask="ThresholdImage ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED} ${CORTICAL_THICKNESS_MASK} 0 0 0 1"
    logCmd $mask

    conversion="ConvertScalarImageToRGB ${DIMENSION} ${CORTICAL_THICKNESS_IMAGE_RESAMPLED}"
    conversion="${conversion} ${CORTICAL_THICKNESS_IMAGE_RGB} none hot none 0 ${DIRECT_THICKNESS_PRIOR}"
    logCmd $conversion

    mosaic="CreateTiledMosaic -i ${HEAD_N4_IMAGE_RESAMPLED} -r ${CORTICAL_THICKNESS_IMAGE_RGB}"
    mosaic="${mosaic} -o ${CORTICAL_THICKNESS_MOSAIC} -a 1.0 -t -1x-1 -d z -p mask"
    mosaic="${mosaic} -s [ 2,mask,mask ] -x ${CORTICAL_THICKNESS_MASK}"
    logCmd $mosaic

    # Segmentation

    echo "0 1 0 0 1 0 1" > $ITKSNAP_COLORMAP
    echo "0 0 1 0 1 1 0" >> $ITKSNAP_COLORMAP
    echo "0 0 0 1 0 1 1" >> $ITKSNAP_COLORMAP

    conversion="ConvertScalarImageToRGB ${DIMENSION} ${BRAIN_SEGMENTATION_IMAGE_RESAMPLED}"
    conversion="${conversion} ${BRAIN_SEGMENTATION_IMAGE_RGB} none custom $ITKSNAP_COLORMAP 0 6"
    logCmd $conversion

    mosaic="CreateTiledMosaic -i ${HEAD_N4_IMAGE_RESAMPLED} -r ${BRAIN_SEGMENTATION_IMAGE_RGB}"
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
  echo ${OUTPUT_PREFIX}ACTStage6Complete.txt > ${OUTPUT_PREFIX}ACTStage6Complete.txt
fi # BAStages
fi # check completion
################################################################################
#
# End of main routine
#
################################################################################

if [[ ${ACT_STAGE} -eq 0 ]] || [[ ${ACT_STAGE} -ge 5 ]] ; then
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
