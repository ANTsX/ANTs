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
if ! command -v antsRegistration &> /dev/null
then
  echo we cant find the antsRegistration program -- does not seem to exist.  please \(re\)define \$PATH in your environment.
  exit
fi
if ! command -v antsApplyTransforms &> /dev/null
then
  echo we cant find the antsApplyTransforms program -- does not seem to exist.  please \(re\)define \$PATH in your environment.
  exit
fi
if ! command -v bc &> /dev/null
then
  echo we cant find the bc program -- does not seem to exist. Please install gnu coreutils.
  exit
fi

################################################################################
#
# General parameters
#
################################################################################

HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"

KEEP_TMP_IMAGES=0

USE_RANDOM_SEEDING=1

DIMENSION=3

ANATOMICAL_IMAGES=()

USE_FLOAT_PRECISION=0

# Intial affine supplied on command line
USER_INITIAL_AFFINE=""

################################################################################
#
# Programs and their parameters
#
################################################################################


ATROPOS=Atropos
ATROPOS_NUM_CLASSES=3
ATROPOS_CSF_CLASS_LABEL=1
ATROPOS_GM_CLASS_LABEL=2
ATROPOS_WM_CLASS_LABEL=3
ATROPOS_BRAIN_EXTRACTION_INITIALIZATION="kmeans[ ${ATROPOS_NUM_CLASSES} ]"
ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD="Gaussian"
ATROPOS_BRAIN_EXTRACTION_CONVERGENCE="[ 3,0.0 ]"

ANTS=antsRegistration
ANTS_MAX_ITERATIONS="100x100x70x20"
ANTS_TRANSFORMATION="SyN[ 0.1,3,0 ]"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[ 1000x500x250x100,1e-8,10 ]"
ANTS_METRIC="CC"
ANTS_METRIC_PARAMS="1,4"

# Default to search an arc fraction of 0.12 in 20 degree intervals
# In other words, roughly +/- 20 degrees from the initial angle
ANTS_AI_ROTATION_SEARCH_PARAMS="20,0.12"

# Default to search a translation of 40mm in the y and z direction
ANTS_AI_TRANSLATION_SEARCH_PARAMS="40,0x40x40"

WARP=antsApplyTransforms

N4=N4BiasFieldCorrection
N4_CONVERGENCE_1="[ 50x50x50x50,0.0000001 ]"
N4_CONVERGENCE_2="[ 50x50x50x50,0.0000001 ]"
N4_SHRINK_FACTOR_1=4
N4_SHRINK_FACTOR_2=2


function Usage {
    cat <<USAGE

`basename $0` performs template-based brain extraction.

Usage:

`basename $0` -d imageDimension
              -a anatomicalImage
              -e brainExtractionTemplate
              -m brainExtractionProbabilityMask
              <OPT_ARGS>
              -o outputPrefix

Example:

  bash $0 -d 3 -a t1.nii.gz -e brainWithSkullTemplate.nii.gz -m brainPrior.nii.gz -o output

Required arguments:

     -d:  Image dimension                       2 or 3 for 2- or 3-dimensional image (default = $DIMENSION)
     -a:  Anatomical image                      Structural image, typically T1.  If more than one
                                                anatomical image is specified, subsequently specified
                                                images are used during the segmentation process.  However,
                                                only the first image is used in the registration of priors.
                                                Our suggestion would be to specify the T1 as the first image.
     -e:  Brain extraction template             Anatomical template.
     -m:  Brain extraction probability mask     Brain probability mask, with intensity range 1 (definitely brain)
                                                to 0 (definitely background).
     -o:  Output prefix                         Output directory + file prefix.

Optional arguments:

     -c:  Tissue classification                 A k-means segmentation is run to find gray or white matter around
                                                the edge of the initial brain mask warped from the template.
                                                This produces a segmentation image with K classes, ordered by mean
                                                intensity in increasing order. With this option, you can control
                                                K and tell the script which classes represent CSF, gray and white matter.
                                                Format (\"KxcsfLabelxgmLabelxwmLabel\").
                                                Examples:
                                                         -c 3x1x2x3 for T1 with K=3, CSF=1, GM=2, WM=3 (default)
                                                         -c 3x3x2x1 for T2 with K=3, CSF=3, GM=2, WM=1
                                                         -c 3x1x3x2 for FLAIR with K=3, CSF=1 GM=3, WM=2
                                                         -c 4x4x2x3 uses K=4, CSF=4, GM=2, WM=3

     -f:  Brain extraction registration mask    Mask used for registration to limit the metric computation to
                                                a specific region.

     -k:  Keep temporary files                  Keep brain extraction/segmentation warps, etc (default = $KEEP_TMP_IMAGES).

     -q:  Use single floating point precision   Use antsRegistration with single (1) or double (0) floating point precision (default = $USE_FLOAT_PRECISION).

     -r:  Initial moving transform              An ITK affine transform (eg, from antsAI or ITK-SNAP) for the moving image.
                                                Without this option, this script calls antsAI to search for a good initial moving
                                                transform.

     -R:  Rotation search parameters            Rotation search parameters for antsAI in format step,arcFraction. The step is in
                                                degrees, the arc fraction goes from 0 (no search) to 1 (search -180 to 180
                                                degree rotations in increements of step). The search begins at -(180*arcFraction)
                                                in each dimension - users should choose parameters so that there is a search point
                                                near zero rotation. Default = $ANTS_AI_ROTATION_SEARCH_PARAMS.

     -s:  Image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz, mhd (default = $OUTPUT_SUFFIX)

     -T:  Translation search parameters         Translation search parameters for antsAI in format step,range. The step is in
                                                mm, -range to range will be tested in each dimension. The default does not search
                                                left-right translations because the brain is usually well-centered along this
                                                dimension in human images. Default = $ANTS_AI_TRANSLATION_SEARCH_PARAMS.

     -u:  Use random seeding                    Use random number generated from system clock (1) or a fixed seed (0). To produce identical
                                                results, multi-threading must also be disabled by setting the environment variable
                                                ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1. Default = $USE_RANDOM_SEEDING.

     -z:  Test / debug mode                     If > 0, runs a faster version of the script. Only for debugging, results will not be good.

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using antsBrainExtraction with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      extraction template     = ${EXTRACTION_TEMPLATE}
      extraction reg. mask    = ${EXTRACTION_REGISTRATION_MASK}
      extraction prior        = ${EXTRACTION_PRIOR}
      output prefix           = ${OUTPUT_PREFIX}
      output image suffix     = ${OUTPUT_SUFFIX}

    antsAI parameters (initial alignment):
      rotation search params  = ${ANTS_AI_ROTATION_SEARCH_PARAMS}
      translation search params = ${ANTS_AI_TRANSLATION_SEARCH_PARAMS}

    N4 parameters (pre brain extraction):
      convergence             = ${N4_CONVERGENCE_1}
      shrink factor           = ${N4_SHRINK_FACTOR_1}
      B-spline parameters     = ${N4_BSPLINE_PARAMS}

    Atropos parameters (extraction):
       convergence            = ${ATROPOS_BRAIN_EXTRACTION_CONVERGENCE}
       likelihood             = ${ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD}
       initialization         = ${ATROPOS_BRAIN_EXTRACTION_INITIALIZATION}
       mrf                    = ${ATROPOS_BRAIN_EXTRACTION_MRF}
       use clock random seed  = ${USE_RANDOM_SEEDING}

PARAMETERS
}


#    local  myresult='some value'
#    echo "$myresult"

# Echos a command to stdout, then runs it
# Will immediately exit on error unless you set debug flag here
DEBUG_MODE=0

function logCmd() {
  cmd="$@"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
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


if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:c:d:e:f:h:k:m:o:q:r:R:s:T:u:z:" OPT
    do
      case $OPT in
          d) #dimensions
       DIMENSION=$OPTARG
       if [[ ${DIMENSION} -gt 4 || ${DIMENSION} -lt 2 ]];
         then
           echo " Error:  ImageDimension must be 2, 3, or 4 "
           exit 1
         fi
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          a) #anatomical t1 image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
       ;;
          c) #k-means segmentation params
       # Check conventional ANTs vector designation (i.e.,  'x')
       kmeansParamsArr=(${OPTARG//x/ })
       if [[ ${#kmeansParamsArr[@]} -ne 4 ]];
         then
           # Check alternative form
           kmeansParamsArr=(${OPTARG//,/ })
           if [[ ${#kmeansParamsArr[@]} -ne 4 ]];
             then
               echo "ERROR:  unrecognized kmeans option (-c)."
               exit 1
             fi
         fi
       ATROPOS_NUM_CLASSES=${kmeansParamsArr[0]}
       ATROPOS_BRAIN_EXTRACTION_INITIALIZATION="kmeans[ ${ATROPOS_NUM_CLASSES} ]"
       ATROPOS_CSF_CLASS_LABEL=${kmeansParamsArr[1]}
       ATROPOS_GM_CLASS_LABEL=${kmeansParamsArr[2]}
       ATROPOS_WM_CLASS_LABEL=${kmeansParamsArr[3]}
       ;;
          k) #keep tmp images
       KEEP_TMP_IMAGES=$OPTARG
       ;;
          e) #brain extraction anatomical image
       EXTRACTION_TEMPLATE=$OPTARG
       ;;
          f) #brain extraction registration mask
       EXTRACTION_REGISTRATION_MASK=$OPTARG
       ;;
          m) #brain extraction prior probability mask
       EXTRACTION_PRIOR=$OPTARG
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          q)
       USE_FLOAT_PRECISION=$OPTARG
       ;;
          r)
       USER_INITIAL_AFFINE=$OPTARG
       ;;
          R)
       ANTS_AI_ROTATION_SEARCH_PARAMS=$OPTARG
       ;;
          s) #output suffix
       OUTPUT_SUFFIX=$OPTARG
       ;;
          T)
       ANTS_AI_TRANSLATION_SEARCH_PARAMS=$OPTARG
       ;;
          u) #use random seeding
       USE_RANDOM_SEEDING=$OPTARG
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

if [[ ${USE_RANDOM_SEEDING} -eq 0 ]]; then
  # Use random seed from Atropos unless one is already defined
  if [[ -z $ANTS_RANDOM_SEED ]] ; then
    export ANTS_RANDOM_SEED=19650218
  fi
fi

ATROPOS_BRAIN_EXTRACTION_MRF="[ 0.1,1x1x1 ]"
if [[ $DIMENSION -eq 2 ]];
  then
    ATROPOS_BRAIN_EXTRACTION_MRF="[ 0.1,1x1 ]"
  fi

if [[ -z "$ATROPOS_SEGMENTATION_MRF" ]];
  then
    ATROPOS_SEGMENTATION_MRF="[ 0.1,1x1x1 ]";
    if [[ $DIMENSION -eq 2 ]];
      then
        ATROPOS_SEGMENTATION_MRF="[ 0.1,1x1 ]"
      fi
  fi


N4_BSPLINE_PARAMS="[ 1x1x1, 3 ]"

if [[ $DIMENSION -eq 2 ]];
  then
    N4_BSPLINE_PARAMS="[ 1x1, 3 ]"
  fi


echo "
Will run Atropos segmentation with K=${ATROPOS_NUM_CLASSES}. Classes labeled in order of mean intensity. Assuming CSF=${ATROPOS_CSF_CLASS_LABEL}, GM=${ATROPOS_GM_CLASS_LABEL}, WM=${ATROPOS_WM_CLASS_LABEL}
"

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

if [[ $DEBUG_MODE -gt 0 ]];
  then

   echo "    WARNING - Running in test / debug mode. Results will be suboptimal "

   # Speed up by doing fewer its. Careful about changing this because
   # certain things are hard coded elsewhere, eg number of levels

   ANTS_MAX_ITERATIONS="40x40x20x0"
   ANTS_LINEAR_CONVERGENCE="[ 100x100x50x10,1e-8,10 ]"

   # Leave N4 / Atropos alone because they're pretty fast

  fi


echoParameters >&2

echo "---------------------  Running `basename $0` on $HOSTNAME  ---------------------"

time_start=`date +%s`

################################################################################
#
# Output image
#
################################################################################

BRAIN_EXTRACTION_MASK=${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}

################################################################################
#
# Brain extraction
#
################################################################################

N4_CORRECTED_IMAGES=()

BRAIN_EXTRACTION_OUTPUT=${OUTPUT_PREFIX}BrainExtraction
EXTRACTION_WARP_OUTPUT_PREFIX=${BRAIN_EXTRACTION_OUTPUT}Prior
EXTRACTION_WARP=${EXTRACTION_WARP_OUTPUT_PREFIX}1Warp.nii.gz
EXTRACTION_INVERSE_WARP=${EXTRACTION_WARP_OUTPUT_PREFIX}1InverseWarp.nii.gz
EXTRACTION_GENERIC_AFFINE=${EXTRACTION_WARP_OUTPUT_PREFIX}0GenericAffine.mat
EXTRACTION_MASK_PRIOR_WARPED=${EXTRACTION_WARP_OUTPUT_PREFIX}Warped.${OUTPUT_SUFFIX}
EXTRACTION_MASK=$BRAIN_EXTRACTION_MASK
EXTRACTION_SEGMENTATION=${BRAIN_EXTRACTION_OUTPUT}Segmentation.${OUTPUT_SUFFIX}
EXTRACTION_BRAIN=${BRAIN_EXTRACTION_OUTPUT}Brain.${OUTPUT_SUFFIX}
EXTRACTION_WM=${BRAIN_EXTRACTION_OUTPUT}WM.${OUTPUT_SUFFIX}
EXTRACTION_GM=${BRAIN_EXTRACTION_OUTPUT}GM.${OUTPUT_SUFFIX}
EXTRACTION_CSF=${BRAIN_EXTRACTION_OUTPUT}CSF.${OUTPUT_SUFFIX}
EXTRACTION_TMP=${BRAIN_EXTRACTION_OUTPUT}Tmp.${OUTPUT_SUFFIX}
EXTRACTION_INITIAL_AFFINE=${BRAIN_EXTRACTION_OUTPUT}InitialAffine.mat
EXTRACTION_INITIAL_AFFINE_FIXED=${BRAIN_EXTRACTION_OUTPUT}InitialAffineFixed.${OUTPUT_SUFFIX}
EXTRACTION_INITIAL_AFFINE_MOVING=${BRAIN_EXTRACTION_OUTPUT}InitialAffineMoving.${OUTPUT_SUFFIX}
EXTRACTION_LAPLACIAN=${BRAIN_EXTRACTION_OUTPUT}Laplacian.${OUTPUT_SUFFIX}
EXTRACTION_TEMPLATE_LAPLACIAN=${BRAIN_EXTRACTION_OUTPUT}TemplateLaplacian.${OUTPUT_SUFFIX}

TMP_FILES=( $EXTRACTION_MASK_PRIOR_WARPED $EXTRACTION_WARP $EXTRACTION_INVERSE_WARP $EXTRACTION_TMP $EXTRACTION_GM $EXTRACTION_CSF $EXTRACTION_SEGMENTATION $EXTRACTION_INITIAL_AFFINE $EXTRACTION_INITIAL_AFFINE_MOVING $EXTRACTION_INITIAL_AFFINE_FIXED $EXTRACTION_LAPLACIAN $EXTRACTION_TEMPLATE_LAPLACIAN $EXTRACTION_WM )

if [[ ! -f ${EXTRACTION_MASK} || ! -f ${EXTRACTION_WM} ]];
  then

    time_start_brain_extraction=`date +%s`

    ################################################################################
    #
    # N4 Correction (pre brain extraction)
    #
    ################################################################################

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Bias correction of anatomical images (pre brain extraction)"
    echo "   1) pre-process by truncating the image intensities"
    echo "   2) run N4"
    echo "--------------------------------------------------------------------------------------"
    echo

    time_start_n4_correction=`date +%s`

    for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
      do
        N4_TRUNCATED_IMAGE=${OUTPUT_PREFIX}N4Truncated${i}.${OUTPUT_SUFFIX}
        N4_CORRECTED_IMAGE=${OUTPUT_PREFIX}N4Corrected${i}.${OUTPUT_SUFFIX}

        TMP_FILES=( ${TMP_FILES[@]} $N4_TRUNCATED_IMAGE $N4_CORRECTED_IMAGE )
        N4_CORRECTED_IMAGES=( ${N4_CORRECTED_IMAGES[@]} ${N4_CORRECTED_IMAGE} )

        if [[ ! -f ${N4_CORRECTED_IMAGE} ]];
          then
            logCmd ImageMath ${DIMENSION} ${N4_TRUNCATED_IMAGE} TruncateImageIntensity ${ANATOMICAL_IMAGES[$i]} 0.01 0.999 256

            exe_n4_correction="${N4} -d ${DIMENSION} -i ${N4_TRUNCATED_IMAGE} -s ${N4_SHRINK_FACTOR_1} -c ${N4_CONVERGENCE_1} -b ${N4_BSPLINE_PARAMS} -o ${N4_CORRECTED_IMAGE} --verbose 1"
            logCmd $exe_n4_correction
          fi
      done

    time_end_n4_correction=`date +%s`
    time_elapsed_n4_correction=$((time_end_n4_correction - time_start_n4_correction))

    ## check if output was produced
    if [[ ! -f ${N4_CORRECTED_IMAGES[0]} ]];
      then
        echo "Expected output was not produced.  The N4 corrected image doesn't exist:"
        echo "   ${N4_CORRECTED_IMAGES[0]}"
        exit 1
      fi

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with N4 correction (pre brain extraction):  $(( time_elapsed_n4_correction / 3600 ))h $(( time_elapsed_n4_correction %3600 / 60 ))m $(( time_elapsed_n4_correction % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo

    if [[ ! -f ${EXTRACTION_INVERSE_WARP} ]];
      then
        if [[ ! -f ${N4_CORRECTED_IMAGES[0]} ]];
          then
            echo "The N4 corrected image doesn't exist:"
            echo "   ${N4_CORRECTED_IMAGES[0]}"
            exit 1
          fi

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Brain extraction using the following steps:"
        echo "   1) Register $EXTRACTION_TEMPLATE to ${N4_CORRECTED_IMAGES[0]}"
        echo "   2) Warp $EXTRACTION_PRIOR to ${ANATOMICAL_IMAGES[0]} using, from 1),"
        echo "      ${OUTPUT_PREFIX}BrainExtractionWarp/Affine"
        echo "   3) Refine segmentation results using Atropos"
        echo "--------------------------------------------------------------------------------------"
        echo

        ## Step 1 ##
          logCmd ImageMath ${DIMENSION} ${EXTRACTION_LAPLACIAN} Laplacian ${N4_CORRECTED_IMAGES[0]} 1.5 1
          logCmd ImageMath ${DIMENSION} ${EXTRACTION_TEMPLATE_LAPLACIAN} Laplacian ${EXTRACTION_TEMPLATE} 1.5 1

          if [[ ! -f "${USER_INITIAL_AFFINE}" ]]
            then

              # Smooth by 4 voxels
              antomical_spacing=($(PrintHeader ${N4_CORRECTED_IMAGES[0]} 1 | tr 'x' '\n'))
              template_spacing=($(PrintHeader ${EXTRACTION_TEMPLATE} 1 | tr 'x' '\n'))

              logCmd SmoothImage ${DIMENSION} ${EXTRACTION_TEMPLATE} 4 ${EXTRACTION_INITIAL_AFFINE_FIXED}
              logCmd SmoothImage ${DIMENSION} ${N4_CORRECTED_IMAGES[0]} 4 ${EXTRACTION_INITIAL_AFFINE_MOVING}

              # Downsample the template by a factor of 5, and resample both images to that resolution
              downsample_template_spacing=()

              for (( i = 0; i < ${DIMENSION}; i++ ))
                do
                  downsample_template_spacing[$i]=$(echo "${template_spacing[$i]} * 5" | bc -l)
                done

              logCmd ResampleImageBySpacing ${DIMENSION} ${EXTRACTION_INITIAL_AFFINE_FIXED} \
                ${EXTRACTION_INITIAL_AFFINE_FIXED} ${downsample_template_spacing[@]} 0

              logCmd ResampleImageBySpacing ${DIMENSION} ${EXTRACTION_INITIAL_AFFINE_MOVING} \
                ${EXTRACTION_INITIAL_AFFINE_MOVING} ${downsample_template_spacing[@]} 0

              exe_initial_align="antsAI -d ${DIMENSION} -v 1"
              exe_initial_align="${exe_initial_align} -m Mattes[ ${EXTRACTION_INITIAL_AFFINE_FIXED},${EXTRACTION_INITIAL_AFFINE_MOVING},32,Regular,0.2 ]"
              exe_initial_align="${exe_initial_align} -t Affine[ 0.1 ]"
              exe_initial_align="${exe_initial_align} -s [ ${ANTS_AI_ROTATION_SEARCH_PARAMS} ]"
              exe_initial_align="${exe_initial_align} -g [ ${ANTS_AI_TRANSLATION_SEARCH_PARAMS} ]"
              exe_initial_align="${exe_initial_align} -p 0"
              exe_initial_align="${exe_initial_align} -c 10"
              exe_initial_align="${exe_initial_align} -o ${EXTRACTION_INITIAL_AFFINE}"

              if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]];
                then
                  exe_initial_align="${exe_initial_align} -x ${EXTRACTION_REGISTRATION_MASK}"
              fi

              logCmd $exe_initial_align
            else
              logCmd antsApplyTransforms -d ${DIMENSION} -t ${USER_INITIAL_AFFINE} -o Linear[ ${EXTRACTION_INITIAL_AFFINE}, 0 ]
            fi

          basecall="${ANTS} -d ${DIMENSION} -u 1 -w [ 0.025,0.975 ] -o ${EXTRACTION_WARP_OUTPUT_PREFIX} -r ${EXTRACTION_INITIAL_AFFINE} -z 1 --float ${USE_FLOAT_PRECISION} --verbose 1"
          if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]];
            then
              basecall="${basecall} -x [ ${EXTRACTION_REGISTRATION_MASK} ]"
            fi
          stage1="-m MI[ ${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[ 0.1 ] -f 8x4x2x1 -s 4x2x1x0"
          stage2="-m MI[ ${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS} ] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[ 0.1 ] -f 8x4x2x1 -s 4x2x1x0"
          stage3="-m CC[ ${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},0.5,4 ] -m CC[ ${EXTRACTION_TEMPLATE_LAPLACIAN},${EXTRACTION_LAPLACIAN},0.5,4 ] -c [ 50x10x0,1e-9,15 ] -t ${ANTS_TRANSFORMATION} -f 4x2x1 -s 2x1x0"

          exe_brain_extraction_1="${basecall} ${stage1} ${stage2} ${stage3}"

          logCmd $exe_brain_extraction_1

        ## check to see if the output registration transforms exist
        if [[ ! -f ${EXTRACTION_GENERIC_AFFINE} ]];
          then
            echo "The registration component of the extraction step didn't complete properly."
            echo "The transform file ${EXTRACTION_GENERIC_AFFINE} does not exist."
            exit 1
          fi

        if [[ ! -f ${EXTRACTION_INVERSE_WARP} ]];
          then
            echo "The registration component of the extraction step didn't complete properly."
            echo "The transform file ${EXTRACTION_INVERSE_WARP} does not exist."
            exit 1
          fi

      fi

      if [[ ! -f ${EXTRACTION_SEGMENTATION} ]];
        then

        ## Step 2 ##

        exe_brain_extraction_2="${WARP} -d ${DIMENSION} -i ${EXTRACTION_PRIOR} -o ${EXTRACTION_MASK_PRIOR_WARPED} -r ${ANATOMICAL_IMAGES[0]} -n Gaussian -t [ ${EXTRACTION_GENERIC_AFFINE},1 ] -t ${EXTRACTION_INVERSE_WARP} --float ${USE_FLOAT_PRECISION} --verbose 1"
        logCmd $exe_brain_extraction_2

        ## superstep 1b ##
        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_MASK_PRIOR_WARPED} ${EXTRACTION_MASK_PRIOR_WARPED} 0.5 1 1 0
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK_PRIOR_WARPED} 2
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} GetLargestComponent ${EXTRACTION_MASK}

        ## superstep 6 ##
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE='';
        for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
          do
            ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${N4_CORRECTED_IMAGES[$i]}";
          done

        exe_brain_extraction_3="${ATROPOS} -d ${DIMENSION} -o ${EXTRACTION_SEGMENTATION} ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -x ${EXTRACTION_MASK} -i ${ATROPOS_BRAIN_EXTRACTION_INITIALIZATION} -c ${ATROPOS_BRAIN_EXTRACTION_CONVERGENCE} -m ${ATROPOS_BRAIN_EXTRACTION_MRF} -k ${ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD} -r ${USE_RANDOM_SEEDING} --verbose 1"
        logCmd $exe_brain_extraction_3
       fi

       # Pad image here to avoid errors from dilating into the edge of the image
        padVoxels=10

        logCmd ImageMath ${DIMENSION} ${EXTRACTION_SEGMENTATION} PadImage ${EXTRACTION_SEGMENTATION} $padVoxels
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK_PRIOR_WARPED} PadImage ${EXTRACTION_MASK_PRIOR_WARPED} $padVoxels

        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_WM} ${ATROPOS_WM_CLASS_LABEL} ${ATROPOS_WM_CLASS_LABEL} 1 0
        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_GM} ${ATROPOS_GM_CLASS_LABEL} ${ATROPOS_GM_CLASS_LABEL} 1 0
        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_CSF} ${ATROPOS_CSF_CLASS_LABEL} ${ATROPOS_CSF_CLASS_LABEL} 1 0

        logCmd ImageMath ${DIMENSION} ${EXTRACTION_WM} GetLargestComponent ${EXTRACTION_WM}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_GM} GetLargestComponent ${EXTRACTION_GM}

        logCmd ImageMath ${DIMENSION} ${EXTRACTION_TMP} FillHoles ${EXTRACTION_GM} 2
        logCmd MultiplyImages ${DIMENSION} ${EXTRACTION_GM} ${EXTRACTION_TMP} ${EXTRACTION_GM}

        logCmd MultiplyImages ${DIMENSION} ${EXTRACTION_WM} ${ATROPOS_WM_CLASS_LABEL} ${EXTRACTION_WM}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_TMP} ME ${EXTRACTION_CSF} 10

        logCmd ImageMath ${DIMENSION} ${EXTRACTION_GM} addtozero ${EXTRACTION_GM} ${EXTRACTION_TMP}
        logCmd MultiplyImages ${DIMENSION} ${EXTRACTION_GM} ${ATROPOS_GM_CLASS_LABEL} ${EXTRACTION_GM}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_SEGMENTATION} addtozero ${EXTRACTION_WM} ${EXTRACTION_GM}

        ## superstep 7 ##
        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_MASK} ${ATROPOS_WM_CLASS_LABEL} ${ATROPOS_WM_CLASS_LABEL} 1 0
        logCmd ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_TMP} ${ATROPOS_GM_CLASS_LABEL} ${ATROPOS_GM_CLASS_LABEL} 1 0
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} addtozero ${EXTRACTION_MASK} ${EXTRACTION_TMP}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} ME ${EXTRACTION_MASK} 2
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} GetLargestComponent ${EXTRACTION_MASK}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK} 4
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} FillHoles ${EXTRACTION_MASK} 2
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} addtozero ${EXTRACTION_MASK} ${EXTRACTION_MASK_PRIOR_WARPED}
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK} 5
        logCmd ImageMath ${DIMENSION} ${EXTRACTION_MASK} ME ${EXTRACTION_MASK} 5

        # De-pad
        for img in ${EXTRACTION_SEGMENTATION} ${EXTRACTION_MASK} ${EXTRACTION_WM} ${EXTRACTION_GM} ${EXTRACTION_CSF} ${EXTRACTION_MASK_PRIOR_WARPED}
          do
            logCmd ImageMath ${DIMENSION} ${img} PadImage ${img} -$padVoxels
          done


        logCmd MultiplyImages ${DIMENSION} ${N4_CORRECTED_IMAGES[0]} ${EXTRACTION_MASK} ${EXTRACTION_BRAIN}

        # Copy header information from original image into output
        logCmd CopyImageHeaderInformation ${ANATOMICAL_IMAGES[0]} ${EXTRACTION_BRAIN} ${EXTRACTION_BRAIN} 1 1 1 0
        logCmd CopyImageHeaderInformation ${ANATOMICAL_IMAGES[0]} ${EXTRACTION_MASK} ${EXTRACTION_MASK} 1 1 1 0


    if [[ ! -f ${EXTRACTION_MASK} ]];
      then
        echo "Expected output was not produced.  The brain mask doesn't exist:"
        echo "   $EXTRACTION_MASK"
        exit 1
      fi
    if [[ ! -f ${EXTRACTION_BRAIN} ]];
      then
        echo "Expected output was not produced.  The brain extracted image doesn't exist:"
        echo "   $EXTRACTION_BRAIN"
        exit 1
      fi

    time_end_brain_extraction=`date +%s`
    time_elapsed_brain_extraction=$((time_end_brain_extraction - time_start_brain_extraction))

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with brain extraction:  $(( time_elapsed_brain_extraction / 3600 ))h $(( time_elapsed_brain_extraction %3600 / 60 ))m $(( time_elapsed_brain_extraction % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo
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


################################################################################
#
# End of main routine
#
################################################################################

exit 0
