#!/bin/bash

VERSION="0.0"

if ! command -v N4BiasFieldCorrection &> /dev/null
then
  echo we cant find the N4 program -- does not seem to exist.  please \(re\)define \$PATH in your environment.
  exit
fi
if ! command -v Atropos &> /dev/null
then
  echo we cant find the Atropos program -- does not seem to exist.  please \(re\)define \$PATH in your environment.
  exit
fi

################################################################################
#
# Main routine
#
################################################################################

HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=""
OUTPUT_PREFIX=""
OUTPUT_SUFFIX="nii.gz"

KEEP_TMP_IMAGES=0

USE_RANDOM_SEEDING=1
DENOISE_ANATOMICAL_IMAGES=0

DIMENSION=""

ANATOMICAL_IMAGES=()

ATROPOS_SEGMENTATION_PRIORS=""

DEBUG_MODE=0

################################################################################
#
# Programs and their parameters
#
################################################################################

N4_ATROPOS_NUMBER_OF_ITERATIONS=15

N4=N4BiasFieldCorrection
N4_CONVERGENCE="[ 50x50x50x50,0.0000000001 ]"
N4_SHRINK_FACTOR=2
N4_BSPLINE_PARAMS="[ 200 ]"
N4_WEIGHT_MASK_POSTERIOR_LABELS=()

ATROPOS=Atropos
ATROPOS_SEGMENTATION_PRIOR_WEIGHT=0.0
ATROPOS_SEGMENTATION_LIKELIHOOD="Gaussian"
ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION="Socrates[ 1 ]"
ATROPOS_SEGMENTATION_MASK=""
ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=5
ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES=""
ATROPOS_SEGMENTATION_ICM=1
ATROPOS_SEGMENTATION_USE_EUCLIDEAN_DISTANCE=0
ATROPOS_SEGMENTATION_MRF=""
ATROPOS_SEGMENTATION_LABEL_PROPAGATION=()

function Usage {
    cat <<USAGE

`basename $0` iterates between N4 <-> Atropos to improve segmentation results.

Usage:

`basename $0` -d imageDimension
              -a inputImage
              -x maskImage
              -c numberOfClasses
              -o outputPrefix
              [Options]

Example:

  bash $0 -d 3 -a t1.nii.gz -x mask.nii.gz -c 4 -m 3 -p segmentationPriors%d.nii.gz -o output

Required arguments:

     -d:  image dimension                       2 or 3, for 2- or 3-dimensional image.
     -a:  input image                           Anatomical image, typically T1.  If more than one anatomical image
                                                is specified, subsequent images are also used during the segmentation process.
     -x:  mask image                            Binary mask defining the region of interest.
     -c:  number of segmentation classes        Number of classes defining the segmentation.
     -o:  output prefix                         The following images are created:
                                                  * \${OUTPUT_PREFIX}N4Corrected.\${OUTPUT_SUFFIX}
                                                  * \${OUTPUT_PREFIX}Segmentation.\${OUTPUT_SUFFIX}
                                                  * \${OUTPUT_PREFIX}SegmentationPosteriors.\${OUTPUT_SUFFIX}

Optional arguments:

     -m:  max. N4 <-> Atropos iterations        Maximum number of (outer loop) iterations between N4 <-> Atropos (default = ${N4_ATROPOS_NUMBER_OF_ITERATIONS}).
     -n:  max. Atropos iterations               Maximum number of (inner loop) iterations in Atropos (default = ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS}).
     -p:  segmentation priors                   Prior probability images initializing the segmentation.
                                                Specified using c-style formatting, e.g. -p labelsPriors%02d.nii.gz. If this
                                                is not specified, k-means initialization is used instead.
     -r:  mrf                                   Specifies MRF prior (of the form '[ weight,neighborhood ]', e.g.
                                                '[ 0.1,1x1x1 ]' which is default).
     -g:  denoise anatomical images             Denoise anatomical images (1) or not (0) (default = ${DENOISE_ANATOMICAL_IMAGES}).
     -b:  posterior formulation                 Posterior formulation and whether or not to use mixture model proportions.
                                                e.g 'Socrates[ 1 ]' (default) or 'Aristotle[ 1 ]'.  Choose the latter if you
                                                want to use the distance priors, see also the -l option for label propagation
                                                control (default = '${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}').
     -l:  label propagation                     Incorporate a distance prior into the 'Aristotle' posterior formulation.  Should be
                                                of the form 'label[ lambda,boundaryProbability ]' where label is a value
                                                of 1,2,3,... denoting label ID.  The label probability for anything
                                                outside the current label

                                                  = boundaryProbability * exp( -lambda * distanceFromBoundary )

                                                Intuitively, smaller lambda values will increase the spatial capture
                                                range of the distance prior.  To apply to all label values, simply omit
                                                specifying the label, i.e. -l '[ lambda,boundaryProbability ]'.
     -y:  posterior label for N4 weight mask    Which posterior probability image should be used to define the
                                                N4 weight mask.  Can also specify multiple posteriors in which
                                                case the chosen posteriors are combined.
     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd.
     -k:  keep temporary files                  Keep temporary files on disk (1) or delete them (0) (default = ${KEEP_TMP_IMAGES}).
     -u:  use random seeding                    Use random number generated from system clock in Atropos (default = ${USE_RANDOM_SEEDING}).
     -w:  Atropos prior segmentation weight     Atropos spatial prior probability weight for the segmentation (default = ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT}).

     -e: N4 convergence                         Convergence parameters for N4, see "-c" option in N4BiasFieldCorrection (default = ${N4_CONVERGENCE}).
     -f: N4 shrink factor                       Shrink factor for N4 (default = ${N4_SHRINK_FACTOR}).
     -q: N4 B-spline parameters                 N4 b-spline specification, see "-b" option in N4BiasFieldCorrection (default = ${N4_BSPLINE_PARAMS}).
     -i: Atropos icm                            ICM parameters for segmentation, see "-g" option in Atropos (default = ${ATROPOS_SEGMENTATION_ICM}).
     -j: Atropos use-euclidean-distance         Use euclidean distances in distance prior formulation (1) or not (0), see Atropos usage for
                                                details (default = ${ATROPOS_SEGMENTATION_USE_EUCLIDEAN_DISTANCE}).

     -z:  Test / debug mode                     If > 0, attempts to continue after errors.

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using antsAtroposN4 with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      segmentation priors     = ${ATROPOS_SEGMENTATION_PRIORS}
      output prefix           = ${OUTPUT_PREFIX}
      output image suffix     = ${OUTPUT_SUFFIX}
      denoise images          = ${DENOISE_ANATOMICAL_IMAGES}

    N4 parameters (segmentation):
      convergence             = ${N4_CONVERGENCE}
      shrink factor           = ${N4_SHRINK_FACTOR}
      B-spline parameters     = ${N4_BSPLINE_PARAMS}
      weight mask post. label = ${N4_WEIGHT_MASK_POSTERIOR_LABELS[@]}

    Atropos parameters (segmentation):
       convergence            = ${ATROPOS_SEGMENTATION_CONVERGENCE}
       likelihood             = ${ATROPOS_SEGMENTATION_LIKELIHOOD}
       prior weight           = ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT}
       posterior formulation  = ${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}
       mrf                    = ${ATROPOS_SEGMENTATION_MRF}
       Max N4->Atropos iters. = ${N4_ATROPOS_NUMBER_OF_ITERATIONS}
       Max Atropos iters.     = ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS}
       use clock random seed  = ${USE_RANDOM_SEEDING}
       icm                    = ${ATROPOS_SEGMENTATION_ICM}
       use-euclidean-distance = ${ATROPOS_SEGMENTATION_USE_EUCLIDEAN_DISTANCE}

PARAMETERS
}


# Echos a command to stdout, then runs it
# Will immediately exit on error unless you set debug flag

function logCmd() {
  cmd="$@"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd

  exec 5>&1
  logCmdOutput=$( "$@" | tee >(cat - >&5); exit ${PIPESTATUS[0]} )

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


if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:w:x:y:z:" OPT
    do
      case $OPT in
          h) #help
            Usage >&2
            exit 0
        ;;
          a) #anatomical t1 image
            ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
        ;;
          b) #atropos prior weight
            ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION=$OPTARG
        ;;
          c) #number of segmentation classes
            ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES=$OPTARG
        ;;
          d) #dimensions
            DIMENSION=$OPTARG
            if [[ ${DIMENSION} -gt 4 || ${DIMENSION} -lt 2 ]];
              then
                echo " Error:  ImageDimension must be 2, 3, or 4 "
                exit 1
              fi
        ;;
          e)
            N4_CONVERGENCE=$OPTARG
        ;;
          f)
            N4_SHRINK_FACTOR=$OPTARG
        ;;
          g) # denoise anatomical images
            DENOISE_ANATOMICAL_IMAGES=$OPTARG
        ;;
          i)
            ATROPOS_SEGMENTATION_ICM=$OPTARG
        ;;
          j)
              ATROPOS_SEGMENTATION_USE_EUCLIDEAN_DISTANCE=$OPTARG
        ;;
          k) #keep tmp images
            KEEP_TMP_IMAGES=$OPTARG
        ;;
          l)
            ATROPOS_SEGMENTATION_LABEL_PROPAGATION[${#ATROPOS_SEGMENTATION_LABEL_PROPAGATION[@]}]=$OPTARG
        ;;
          m) #atropos segmentation iterations
            N4_ATROPOS_NUMBER_OF_ITERATIONS=$OPTARG
        ;;
          n) #atropos segmentation iterations
            ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=$OPTARG
        ;;
          o) #output prefix
            OUTPUT_PREFIX=$OPTARG
        ;;
          p) # segmentation label prior image
            ATROPOS_SEGMENTATION_PRIORS=$OPTARG
        ;;
          q)
              N4_BSPLINE_PARAMS=$OPTARG
        ;;
          r) #mrf
            ATROPOS_SEGMENTATION_MRF=$OPTARG
        ;;
          s) #output suffix
            OUTPUT_SUFFIX=$OPTARG
        ;;
          t) #n4 convergence
            N4_CONVERGENCE=$OPTARG
        ;;
          u) #use random seeding
            USE_RANDOM_SEEDING=$OPTARG
        ;;
          w) #atropos prior weight
            ATROPOS_SEGMENTATION_PRIOR_WEIGHT=$OPTARG
        ;;
          x) #atropos segmentation mask
            ATROPOS_SEGMENTATION_MASK=$OPTARG
        ;;
          y) #
            N4_WEIGHT_MASK_POSTERIOR_LABELS[${#N4_WEIGHT_MASK_POSTERIOR_LABELS[@]}]=$OPTARG
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

# Check required args
if [[ $DIMENSION -lt 2 ]] || [[ $DIMENSION -gt 3 ]];
  then
    echo "Script only supports image dimension 2 or 3"
    exit 1
  fi

if [[ ${#ANATOMICAL_IMAGES[@]} -eq 0 ]];
  then
    echo "No input anatomical images to segment"
    exit 1
  fi

if [[ ! ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES} -gt 0 ]];
  then
    echo "Require number of segmentation classes"
    exit 1
  fi

if [[ -z "${OUTPUT_PREFIX}" ]];
  then
    echo "Output prefix is required"
    exit 1
  fi

if [[ -z "$ATROPOS_SEGMENTATION_MRF" ]];
  then
    ATROPOS_SEGMENTATION_MRF="[ 0.1,1x1x1 ]";
    if [[ DIMENSION -eq 2 ]];
      then
        ATROPOS_SEGMENTATION_MRF="[ 0.1,1x1 ]"
      fi
  fi

ATROPOS_SEGMENTATION_CONVERGENCE="[ ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS},0.0 ]"

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

if [[ ! -f "${ATROPOS_SEGMENTATION_MASK}" ]];
  then
    echo "Required mask image \"${ATROPOS_SEGMENTATION_MASK}\" does not exist."
    exit 1
  fi

FORMAT=${ATROPOS_SEGMENTATION_PRIORS}
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
POSTERIOR_IMAGE_FILENAMES=()
POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION=()
for (( i = 1; i <= $ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES; i++ ))
  do
    NUMBER_OF_REPS=$(( $TOTAL_LENGTH - ${#i} ))
    ROOT='';
    for(( j=0; j < $NUMBER_OF_REPS; j++ ))
      do
        ROOT=${ROOT}${REPCHARACTER}
      done
    PRIOR_FILENAME=${PREFORMAT}${ROOT}${i}${POSTFORMAT}
    POSTERIOR_FILENAME=${OUTPUT_PREFIX}SegmentationPosteriors${ROOT}${i}.${OUTPUT_SUFFIX}
    POSTERIOR_FILENAME_PREVIOUS_ITERATION=${OUTPUT_PREFIX}SegmentationPosteriorsPreviousIteration${ROOT}${i}.${OUTPUT_SUFFIX}
    POSTERIOR_IMAGE_FILENAMES=( ${POSTERIOR_IMAGE_FILENAMES[@]} $POSTERIOR_FILENAME )
    POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION=( ${POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION[@]} $POSTERIOR_FILENAME_PREVIOUS_ITERATION )
    if [[ -f $PRIOR_FILENAME ]];
      then
        PRIOR_IMAGE_FILENAMES=( ${PRIOR_IMAGE_FILENAMES[@]} $PRIOR_FILENAME )
      fi
  done

NUMBER_OF_PRIOR_IMAGES=${#PRIOR_IMAGE_FILENAMES[*]}

INITIALIZE_WITH_KMEANS=0
if [[ -z "${ATROPOS_SEGMENTATION_PRIORS}" ]];
  then
    echo "Initializing with kmeans segmentation."
    INITIALIZE_WITH_KMEANS=1
elif [[ ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES} -ne ${NUMBER_OF_PRIOR_IMAGES} ]];
  then
    echo "Expected ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES} prior images (${NUMBER_OF_PRIOR_IMAGES} are specified).  Check the command line specification."
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

ATROPOS_SEGMENTATION_OUTPUT=${OUTPUT_PREFIX}Segmentation
ATROPOS_SEGMENTATION=${ATROPOS_SEGMENTATION_OUTPUT}.${OUTPUT_SUFFIX}
ATROPOS_SEGMENTATION_POSTERIORS=${ATROPOS_SEGMENTATION_OUTPUT}Posteriors%${FORMAT}d.${OUTPUT_SUFFIX}

################################################################################
#
# Preprocess anatomical images
#    1. Truncate input intensity (needed for N4)
#    2. Denoise image (if requested)
#
################################################################################

if [[ ${DENOISE_ANATOMICAL_IMAGES} -ne 0 ]];
  then
    if ! command -v DenoiseImage &> /dev/null
      then
        echo "Error:  we can't find the DenoiseImage program."
        echo "Perhaps you need to \(re\)define \$PATH in your environment or update your repository."
        exit
      fi
  fi

# Anatomical images that have been truncated, and optionally denoised
SEGMENTATION_PREPROCESSED_IMAGES=()

for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
  do
    SEGMENTATION_PREPROCESSED_IMAGES=( ${SEGMENTATION_PREPROCESSED_IMAGES[@]} ${ATROPOS_SEGMENTATION_OUTPUT}PreprocessedAnatomical${j}.${OUTPUT_SUFFIX} )
    # Truncate on the whole head to get outliers over the whole volume, without losing contrast in the brain
    logCmd ImageMath ${DIMENSION} ${SEGMENTATION_PREPROCESSED_IMAGES[$j]} TruncateImageIntensity ${ANATOMICAL_IMAGES[$j]} 0 0.995 256
    if [[ ${DENOISE_ANATOMICAL_IMAGES} -ne 0 ]];
      then
        logCmd DenoiseImage -d ${DIMENSION} -i ${SEGMENTATION_PREPROCESSED_IMAGES[$j]} -o ${SEGMENTATION_PREPROCESSED_IMAGES[$j]} --verbose 1
      fi
  done

################################################################################
#
# Segmentation
#
################################################################################

SEGMENTATION_WEIGHT_MASK=${OUTPUT_PREFIX}SegmentationWeightMask.nii.gz
SEGMENTATION_CONVERGENCE_FILE=${OUTPUT_PREFIX}SegmentationConvergence.txt
SEGMENTATION_PREVIOUS_ITERATION=${OUTPUT_PREFIX}SegmentationPreviousIteration.${OUTPUT_SUFFIX}

N4_WEIGHT_MASK_POSTERIOR_IDXS=()
for (( i = 0; i < ${#N4_WEIGHT_MASK_POSTERIOR_LABELS[@]}; i++ ))
  do
    N4_WEIGHT_MASK_POSTERIOR_IDXS[$i]=$((N4_WEIGHT_MASK_POSTERIOR_LABELS[$i]-1))
  done

time_start_segmentation=`date +%s`

if [[ $INITIALIZE_WITH_KMEANS -eq 0 ]]
  then

    N4_WEIGHT_MASK_IMAGES=()
    for (( i = 0; i < ${#N4_WEIGHT_MASK_POSTERIOR_LABELS[@]}; i++ ))
      do
        N4_WEIGHT_MASK_IMAGES=( ${N4_WEIGHT_MASK_IMAGES[@]} ${PRIOR_IMAGE_FILENAMES[${N4_WEIGHT_MASK_POSTERIOR_IDXS[$i]}]} )
      done

    if [[ ${#N4_WEIGHT_MASK_IMAGES[@]} -gt 0 ]];
      then
        logCmd ImageMath ${DIMENSION} ${SEGMENTATION_WEIGHT_MASK} PureTissueN4WeightMask ${N4_WEIGHT_MASK_IMAGES[@]}
      fi
  fi

if [[ -f ${SEGMENTATION_CONVERGENCE_FILE} ]];
  then
    logCmd rm -f ${SEGMENTATION_CONVERGENCE_FILE}
  fi

POSTERIOR_PROBABILITY_CONVERGED=0
for (( i = 0; i < ${N4_ATROPOS_NUMBER_OF_ITERATIONS}; i++ ))
  do
    SEGMENTATION_N4_IMAGES=()
    for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
      do
        SEGMENTATION_N4_IMAGES=( ${SEGMENTATION_N4_IMAGES[@]} ${ATROPOS_SEGMENTATION_OUTPUT}${j}N4.${OUTPUT_SUFFIX} )
        exe_n4_correction="${N4} -d ${DIMENSION} -i ${SEGMENTATION_PREPROCESSED_IMAGES[$j]} -x ${ATROPOS_SEGMENTATION_MASK} -s ${N4_SHRINK_FACTOR} -c ${N4_CONVERGENCE} -b ${N4_BSPLINE_PARAMS} -o ${SEGMENTATION_N4_IMAGES[$j]} --verbose 1"
        if [[ -f ${SEGMENTATION_WEIGHT_MASK} ]];
          then
            exe_n4_correction="${exe_n4_correction} -w ${SEGMENTATION_WEIGHT_MASK}"
          fi
        logCmd $exe_n4_correction
        logCmd ImageMath ${DIMENSION} ${SEGMENTATION_N4_IMAGES[$j]} Normalize ${SEGMENTATION_N4_IMAGES[$j]}
        logCmd ImageMath ${DIMENSION} ${SEGMENTATION_N4_IMAGES[$j]} m ${SEGMENTATION_N4_IMAGES[$j]} 1000
      done

    ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE=''
    for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
      do
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${SEGMENTATION_N4_IMAGES[$j]}"
      done

    INITIALIZATION="PriorProbabilityImages[ ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES},${ATROPOS_SEGMENTATION_PRIORS},${ATROPOS_SEGMENTATION_PRIOR_WEIGHT}]"
    if [[ INITIALIZE_WITH_KMEANS -eq 1 ]];
      then
        if [[ $i -eq 0 ]];
          then
            INITIALIZATION="kmeans[ ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES} ]"
          else
            INITIALIZATION="PriorProbabilityImages[ ${ATROPOS_SEGMENTATION_NUMBER_OF_CLASSES},${ATROPOS_SEGMENTATION_POSTERIORS},${ATROPOS_SEGMENTATION_PRIOR_WEIGHT} ]"
          fi
      fi

    ATROPOS_LABEL_PROPAGATION_COMMAND_LINE=''
    for (( j = 0; j < ${#ATROPOS_SEGMENTATION_LABEL_PROPAGATION[@]}; j++ ))
      do
        ATROPOS_LABEL_PROPAGATION_COMMAND_LINE="${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} -l ${ATROPOS_SEGMENTATION_LABEL_PROPAGATION[$j]}";
      done

    exe_segmentation="${ATROPOS} -d ${DIMENSION} -x ${ATROPOS_SEGMENTATION_MASK} -c ${ATROPOS_SEGMENTATION_CONVERGENCE} ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} ${ATROPOS_LABEL_PROPAGATION_COMMAND_LINE} --verbose 1"
    exe_segmentation="${exe_segmentation} -i ${INITIALIZATION} -k ${ATROPOS_SEGMENTATION_LIKELIHOOD} -m ${ATROPOS_SEGMENTATION_MRF} -g ${ATROPOS_SEGMENTATION_ICM} -o [ ${ATROPOS_SEGMENTATION},${ATROPOS_SEGMENTATION_POSTERIORS} ] -r ${USE_RANDOM_SEEDING} -e ${ATROPOS_SEGMENTATION_USE_EUCLIDEAN_DISTANCE}"

    if [[ $i -eq 0 ]];
      then
        exe_segmentation="${exe_segmentation} -p Socrates[ 0 ]"
      else
        exe_segmentation="${exe_segmentation} -p ${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}"

        logCmd cp -f ${ATROPOS_SEGMENTATION} ${SEGMENTATION_PREVIOUS_ITERATION}

        for (( j = 0; j < ${#POSTERIOR_IMAGE_FILENAMES[@]}; j++ ))
          do
            logCmd cp -f ${POSTERIOR_IMAGE_FILENAMES[$j]} ${POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION[$j]}
          done

        for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
          do
            ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${SEGMENTATION_N4_IMAGES[$j]}";
          done
      fi

    logCmd $exe_segmentation

    if [[ $i -eq 0 ]];
      then
        if [[ ! -f ${SEGMENTATION_CONVERGENCE_FILE} ]];
          then
            echo "Iteration,Posterior" > ${SEGMENTATION_CONVERGENCE_FILE}
          fi

        POSTERIOR_PROBABILITY=0
        while read line;
          do
            tokens=( $line )
            if [[ ${tokens[0]} == "Iteration" ]];
              then
                POSTERIOR_PROBABILITY=${tokens[7]}
              fi
          done <<< "$logCmdOutput"

        echo "${i},${POSTERIOR_PROBABILITY}" >> ${SEGMENTATION_CONVERGENCE_FILE}
      fi

    if [[ $i -gt 0 && -f ${SEGMENTATION_PREVIOUS_ITERATION} ]];
      then

        POSTERIOR_PROBABILITY_PREVIOUS_ITERATION=$POSTERIOR_PROBABILITY

        POSTERIOR_PROBABILITY=0
        while read line;
          do
            tokens=( $line )
            if [[ ${tokens[0]} == "Iteration" ]];
              then
                POSTERIOR_PROBABILITY=${tokens[7]}
              fi
          done <<< "$logCmdOutput"

        if [[ $( echo "${POSTERIOR_PROBABILITY} < ${POSTERIOR_PROBABILITY_PREVIOUS_ITERATION}"|bc ) -eq 1 ]];
          then
            POSTERIOR_PROBABILITY_CONVERGED=1

            POSTERIOR_PROBABILITY=${POSTERIOR_PROBABILITY_PREVIOUS_ITERATION}
            logCmd cp -f ${SEGMENTATION_PREVIOUS_ITERATION} ${ATROPOS_SEGMENTATION}

            for (( j = 0; j < ${#POSTERIOR_IMAGE_FILENAMES[@]}; j++ ))
              do
                logCmd cp -f ${POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION[$j]} ${POSTERIOR_IMAGE_FILENAMES[$j]}
              done

            break
          else
            echo "${i},${POSTERIOR_PROBABILITY}" >> ${SEGMENTATION_CONVERGENCE_FILE}
          fi
      fi

    N4_WEIGHT_MASK_IMAGES=()
    for (( j = 0; j < ${#N4_WEIGHT_MASK_POSTERIOR_LABELS[@]}; j++ ))
      do
        N4_WEIGHT_MASK_IMAGES=( ${N4_WEIGHT_MASK_IMAGES[@]} ${POSTERIOR_IMAGE_FILENAMES[${N4_WEIGHT_MASK_POSTERIOR_IDXS[$j]}]} )
      done

    if [[ ${#N4_WEIGHT_MASK_IMAGES[@]} -gt 0 ]];
      then
        logCmd ImageMath ${DIMENSION} ${SEGMENTATION_WEIGHT_MASK} PureTissueN4WeightMask ${N4_WEIGHT_MASK_IMAGES[@]}
      fi

  done

TMP_FILES=( $SEGMENTATION_WEIGHT_MASK ${POSTERIOR_IMAGE_FILENAMES_PREVIOUS_ITERATION[@]} ${SEGMENTATION_PREVIOUS_ITERATION} ${SEGMENTATION_PREPROCESSED_IMAGES[@]} )

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

time_end_segmentation=`date +%s`
time_elapsed_segmentation=$((time_end_segmentation - time_start_segmentation))

echo
echo "--------------------------------------------------------------------------------------"
if [[ POSTERIOR_PROBABILITY_CONVERGED -eq 1 ]];
  then
    echo " Done with segmentation (posterior prob. converged):  $(( time_elapsed_segmentation / 3600 ))h $(( time_elapsed_segmentation %3600 / 60 ))m $(( time_elapsed_segmentation % 60 ))s"
  else
    echo " Done with segmentation (exceeded max. iterations):  $(( time_elapsed_segmentation / 3600 ))h $(( time_elapsed_segmentation %3600 / 60 ))m $(( time_elapsed_segmentation % 60 ))s"
  fi
echo "--------------------------------------------------------------------------------------"
echo


################################################################################
#
# End of main routine
#
################################################################################

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with N4 <-> Atropos processing"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

SEGMENTATION_CONVERGENCE_SCRIPT=${ATROPOS_SEGMENTATION_OUTPUT}Convergence.R
SEGMENTATION_CONVERGENCE_PLOT=${ATROPOS_SEGMENTATION_OUTPUT}Convergence.pdf

if [[ `type -p RScript` > /dev/null ]];
  then
    echo "library( ggplot2 )" > $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "conv <- read.csv( \"${SEGMENTATION_CONVERGENCE_FILE}\" )" >>  $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "myPlot <- ggplot( conv, aes( x = Iteration, y = Posterior ) ) +" >>  $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "  geom_point( data = conv, aes( colour = Iteration ), size = 4 ) +" >>  $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "  scale_y_continuous( breaks = seq( 0.8  , 1, by = 0.025 ), labels = seq( 0.8, 1, by = 0.025 ), limits = c( 0.8, 1 ) ) +" >>  $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "  theme( legend.position = \"none\" )" >>  $SEGMENTATION_CONVERGENCE_SCRIPT
    echo "ggsave( filename = \"$SEGMENTATION_CONVERGENCE_PLOT\", plot = myPlot, width = 4, height = 3, units = 'in' )" >>  $SEGMENTATION_CONVERGENCE_SCRIPT

    `RScript $SEGMENTATION_CONVERGENCE_SCRIPT`
    rm -f $SEGMENTATION_CONVERGENCE_SCRIPT
  fi

exit 0
