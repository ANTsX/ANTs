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
if [[ ! -s ${ANTSPATH}/N4BiasFieldCorrection ]]; then
  echo we cant find the N4 program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
if [[ ! -s ${ANTSPATH}/Atropos ]]; then
  echo we cant find the Atropos program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
if [[ ! -s ${ANTSPATH}/KellyKapowski ]]; then
  echo we cant find the DiReCT \(aka KellyKapowski\) program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

function Usage {
    cat <<USAGE

`basename $0` performs T1 anatomical brain processing where the following steps are currently applied:

  1. Brain extraction
  2. Brain 3-tissue segmentation
  3. Cortical thickness
  4. (Optional) registration to a template

Usage:

`basename $0` -d imageDimension
              -a anatomicalImage
              -e brainExtractionTemplate
              -m brainExtractionProbabilityMask
              -l brainSegmentationTemplate
              -p brainSegmentationPriors
              <OPTARGS>
              -o outputPrefix

Example:

  bash $0 -d 3 -i t1.nii.gz -e brainWithSkullTemplate.nii.gz -m brainPrior.nii.gz -l segmentationTemplate.nii.gz -p segmentationPriors%d.nii.gz -o output

Required arguments:

     -d:  Image dimension                       2 or 3 (for 2- or 3-dimensional image)
     -a:  Anatomical image                      Structural image, typically T1.  If more than one
                                                anatomical image is specified, subsequently specified
                                                images are used during the segmentation process.  However,
                                                only the first image is used in the registration of priors.
                                                Our suggestion would be to specify the T1 as the first image.
     -e:  Brain extraction template             Anatomical template created using e.g. LPBA40 data set with
                                                buildtemplateparallel.sh in ANTs.
     -m:  Brain extraction probability mask     Brain probability mask created using e.g. LPBA40 data set which
                                                have brain masks defined, and warped to anatomical template and
                                                averaged resulting in a probability image.
     -l:  Brain segmentation template           Anatomical template for brain segmentation.
     -p:  Brain segmentation priors             Label probability priors corresponding to the image specified
                                                with the -l option.  Specified using c-style formatting, e.g.
                                                -p labelsPriors%02d.nii.gz.
     -o:  Output prefix                         The following images are created:
                                                  * ${OUTPUT_PREFIX}N4Corrected.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}ExtractedBrain.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors1.${OUTPUT_SUFFIX}  CSF
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors2.${OUTPUT_SUFFIX}  GM
                                                  * ${OUTPUT_PREFIX}BrainSegmentationPosteriors3.${OUTPUT_SUFFIX}  WM
                                                  * ${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}

Optional arguments:

     -f:  Brain extraction registration mask    Mask used for registration to limit the metric computation to
                                                a specific region.
     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd
     -t:  template for t1 registration
     -k:  keep temporary files                  Keep brain extraction/segmentation warps, etc (default = false).
     -i:  max iterations for registration       ANTS registration max iterations (default = 100x100x70x20)
     -w:  Atropos prior segmentation weight     Atropos spatial prior probability weight for the segmentation (default = 0)
     -n:  number of segmentation iterations     N4 -> Atropos -> N4 iterations during segmentation (default = 15)

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using apb with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      extraction template     = ${EXTRACTION_TEMPLATE}
      extraction reg. mask    = ${EXTRACTION_REGISTRATION_MASK}
      extraction prior        = ${EXTRACTION_PRIOR}
      segmentation template   = ${SEGMENTATION_TEMPLATE}
      segmentation prior      = ${SEGMENTATION_PRIOR}
      output prefix           = ${OUTPUT_PREFIX}
      output image suffix     = ${OUTPUT_SUFFIX}
      registration template   = ${REGISTRATION_TEMPLATE}

    ANTs parameters:
      metric                  = ${ANTS_METRIC}[fixedImage,movingImage,${ANTS_METRIC_PARAMS}]
      regularization          = ${ANTS_REGULARIZATION}
      transformation          = ${ANTS_TRANSFORMATION}
      max iterations          = ${ANTS_MAX_ITERATIONS}

    N4 parameters (pre brain extraction):
      convergence             = ${N4_CONVERGENCE_1}
      shrink factor           = ${N4_SHRINK_FACTOR_1}
      B-spline parameters     = ${N4_BSPLINE_PARAMS}
    N4 parameters (segmentation):
      convergence             = ${N4_CONVERGENCE_2}
      shrink factor           = ${N4_SHRINK_FACTOR_2}
      B-spline parameters     = ${N4_BSPLINE_PARAMS}

    Atropos parameters (extraction):
       convergence            = ${ATROPOS_BRAIN_EXTRACTION_CONVERGENCE}
       likelihood             = ${ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD}
       initialization         = ${ATROPOS_BRAIN_EXTRACTION_INITIALIZATION}
       mrf                    = ${ATROPOS_BRAIN_EXTRACTION_MRF}
    Atropos parameters (segmentation):
       convergence            = ${ATROPOS_SEGMENTATION_CONVERGENCE}
       likelihood             = ${ATROPOS_SEGMENTATION_LIKELIHOOD}
       prior weight           = ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT}
       initialization         = ${ATROPOS_SEGMENTATION_INITIALIZATION}
       posterior formulation  = ${ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION}
       mrf                    = ${ATROPOS_SEGMENTATION_MRF}
       Max N4->Atropos iters. = ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS}
       N4->Atropos threshold  = ${ATROPOS_SEGMENTATION_DICE_THRESHOLD}

    DiReCT parameters:
      convergence             = ${DIRECT_CONVERGENCE}
      thickness prior         = ${DIRECT_THICKNESS_PRIOR}
      gradient step size      = ${DIRECT_GRAD_STEP_SIZE}
      smoothing sigma         = ${DIRECT_SMOOTHING_SIGMA}

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

KEEP_TMP_IMAGES='true'

DIMENSION=3

ANATOMICAL_IMAGES=()
REGISTRATION_TEMPLATE=""

EXTRACTION_TEMPLATE=""
EXTRACTION_REGISTRATION_MASK=""
EXTRACTION_PRIOR=""
SEGMENTATION_TEMPLATE=""
SEGMENTATION_PRIOR=""
WHITE_MATTER_LABEL=3
GRAY_MATTER_LABEL=2
CSF_MATTER_LABEL=1

################################################################################
#
# Programs and their parameters
#
################################################################################

ANTS=${ANTSPATH}antsRegistration
ANTS_MAX_ITERATIONS="100x100x70x20"
ANTS_TRANSFORMATION="SyN[0.1,3,0]"
ANTS_LINEAR_METRIC_PARAMS="1,32,Regular,0.25"
ANTS_LINEAR_CONVERGENCE="[1000x500x250x100,1e-8,10]"
ANTS_METRIC="CC"
ANTS_METRIC_PARAMS="1,4"

WARP=${ANTSPATH}antsApplyTransforms

N4=${ANTSPATH}N4BiasFieldCorrection
N4_CONVERGENCE_1="[50x50x50x50,0.0000001]"
N4_CONVERGENCE_2="[50x50x50x50,0.0000001]"
N4_SHRINK_FACTOR_1=4
N4_SHRINK_FACTOR_2=2
N4_BSPLINE_PARAMS="[200]"

ATROPOS=${ANTSPATH}Atropos
ATROPOS_BRAIN_EXTRACTION_INITIALIZATION="kmeans[3]"
ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD="Gaussian"
ATROPOS_BRAIN_EXTRACTION_CONVERGENCE="[3,0.0]"

ATROPOS_SEGMENTATION_INITIALIZATION="PriorProbabilityImages"
ATROPOS_SEGMENTATION_PRIOR_WEIGHT=0.0
ATROPOS_SEGMENTATION_LIKELIHOOD="Gaussian"
ATROPOS_SEGMENTATION_CONVERGENCE="[5,0.0]"
ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION="Socrates"
ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=15
ATROPOS_SEGMENTATION_DICE_THRESHOLD=0.99

DIRECT=${ANTSPATH}KellyKapowski
DIRECT_CONVERGENCE="[45,0.0,10]";
DIRECT_THICKNESS_PRIOR="10";
DIRECT_GRAD_STEP_SIZE="0.025";
DIRECT_SMOOTHING_SIGMA="1.5";
DIRECT_NUMBER_OF_DIFF_COMPOSITIONS="10";

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:d:e:f:h:i:k:l:m:n:p:o:s:t:v:w:" OPT
    do
      case $OPT in
          a) #anatomical t1 image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
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
       EXTRACTION_TEMPLATE=$OPTARG
       ;;
          f) #brain extraction registration mask
       EXTRACTION_REGISTRATION_MASK=$OPTARG
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          i) #max_iterations
       ANTS_MAX_ITERATIONS=$OPTARG
       ;;
          k) #keep tmp images
       KEEP_TMP_IMAGES=$OPTARG
       ;;
          l) #brain segmentation label anatomical image
       SEGMENTATION_TEMPLATE=$OPTARG
       ;;
          m) #brain extraction prior probability mask
       EXTRACTION_PRIOR=$OPTARG
       ;;
          n) #atropos segmentation iterations
       ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS=$OPTARG
       ;;
          p) #brain segmentation label prior image
       SEGMENTATION_PRIOR=$OPTARG
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          s) #output suffix
       OUTPUT_SUFFIX=$OPTARG
       ;;
          t) #template registration image
       REGISTRATION_TEMPLATE=$OPTARG
       ;;
          w) #atropos prior weight
       ATROPOS_SEGMENTATION_PRIOR_WEIGHT=$OPTARG
       ;;
          *) # getopts issues an error message
       echo "ERROR:  unrecognized option -$OPT $OPTARG"
       exit 1
       ;;
      esac
  done
fi

ATROPOS_BRAIN_EXTRACTION_MRF="[0.2,1x1x1]"
if [[ DIMENSION -eq 2 ]];
  then
    ATROPOS_BRAIN_EXTRACTION_MRF="[0.2,1x1]"
  fi
ATROPOS_SEGMENTATION_MRF="[0.1,1x1x1]";
if [[ DIMENSION -eq 2 ]];
  then
    ATROPOS_SEGMENTATION_MRF="[0.1,1x1]"
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

if [[ ! -f ${EXTRACTION_TEMPLATE} ]];
  then
    echo "The extraction template doesn't exist:"
    echo "   $EXTRACTION_TEMPLATE"
    exit 1
  fi
if [[ ! -f ${EXTRACTION_PRIOR} ]];
  then
    echo "The brain extraction prior doesn't exist:"
    echo "   $EXTRACTION_PRIOR"
    exit 1
  fi


## check segmentation inputs
if [[ ! -f ${SEGMENTATION_TEMPLATE} ]];
  then
    echo "The segmentation template doesn't exist:"
    echo "   $SEGMENTATION_TEMPLATE"
    exit 1
  fi

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

BRAIN_SEGMENTATION_POSTERIORS=${BRAIN_SEGMENTATION_OUTPUT}Posteriors%${FORMAT}d.${OUTPUT_SUFFIX}



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

if [[ ${NUMBER_OF_PRIOR_IMAGES} -ne 3 ]];
  then
    echo "Expected 3 prior images (${NUMBER_OF_PRIOR_IMAGES} are specified).  Check the command line specification."
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

# check to see if input extraction template and input segmentation are the same
# if so, we can initialize the segmentation registration with the affine transform
# derived from the extraction step.
INPUT_TEMPLATES_ARE_THE_SAME=0
if [[ `diff $EXTRACTION_TEMPLATE $SEGMENTATION_TEMPLATE >/dev/null` ]];
  then
    echo INPUT_TEMPLATES_ARE_THE_SAME=1
  fi

OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
if [[ ! -e $OUTPUT_PREFIX ]];
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

N4_CORRECTED_IMAGES=()
BRAIN_EXTRACTION_MASK=${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
BRAIN_SEGMENTATION=${OUTPUT_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_IMAGE=${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}

################################################################################
#
# Brain extraction
#
################################################################################

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

TMP_FILES=( $EXTRACTION_MASK_PRIOR_WARPED $EXTRACTION_WARP $EXTRACTION_INVERSE_WARP $EXTRACTION_TMP $EXTRACTION_GM $EXTRACTION_CSF $EXTRACTION_SEGMENTATION $EXTRACTION_INITIAL_AFFINE $EXTRACTION_INITIAL_AFFINE_MOVING $EXTRACTION_INITIAL_AFFINE_FIXED $EXTRACTION_LAPLACIAN $EXTRACTION_TEMPLATE_LAPLACIAN )

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
            logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${N4_TRUNCATED_IMAGE} TruncateImageIntensity ${ANATOMICAL_IMAGES[$i]} 0.025 0.975 256

            exe_n4_correction="${N4} -d ${DIMENSION} -i ${N4_TRUNCATED_IMAGE} -s ${N4_SHRINK_FACTOR_1} -c ${N4_CONVERGENCE_1} -b ${N4_BSPLINE_PARAMS} -o ${N4_CORRECTED_IMAGE}"
            logCmd $exe_n4_correction
          fi
      done

    time_end_n4_correction=`date +%s`
    time_elapsed_n4_correction=$((time_end_n4_correction - time_start_n4_correction))

    ## check if output was produced
    if [[ ! -f ${N4_CORRECTED_IMAGES[0]} ]];
      then
        echo "Expected output was not produce.  The N4 corrected image doesn't exist:"
        echo "   ${N4_CORRECTED_IMAGES[0]}"
        exit 1
      fi

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Done with N4 correction (pre brain extraction):  $(( time_elapsed_n4_correction / 3600 ))h $(( time_elapsed_n4_correction %3600 / 60 ))m $(( time_elapsed_n4_correction % 60 ))s"
    echo "--------------------------------------------------------------------------------------"
    echo

    if [[ ! -f ${EXTRACTION_MASK} ]];
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
        if [[ ! -f ${EXTRACTION_WARP} ]];
          then
            logCmd ${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${EXTRACTION_TEMPLATE} ${EXTRACTION_INITIAL_AFFINE_FIXED} 4 4 4 1
            logCmd ${ANTSPATH}/ResampleImageBySpacing ${DIMENSION} ${N4_CORRECTED_IMAGES[0]} ${EXTRACTION_INITIAL_AFFINE_MOVING} 4 4 4 1

            logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTION_LAPLACIAN} Laplacian ${N4_CORRECTED_IMAGES[0]} 1.5 1
            logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${EXTRACTION_TEMPLATE_LAPLACIAN} Laplacian ${EXTRACTION_TEMPLATE} 1.5 1

            exe_initial_align="${ANTSPATH}/antsAffineInitializer ${DIMENSION} ${EXTRACTION_INITIAL_AFFINE_FIXED} ${EXTRACTION_INITIAL_AFFINE_MOVING} ${EXTRACTION_INITIAL_AFFINE} 15 0.1 0 10"
            if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]];
              then
                exe_initial_align="${exe_initial_align} ${EXTRACTION_REGISTRATION_MASK}"
              fi
            logCmd $exe_initial_align

            basecall="${ANTS} -d ${DIMENSION} -u 1 -w [0.025,0.975] -o ${EXTRACTION_WARP_OUTPUT_PREFIX} -r ${EXTRACTION_INITIAL_AFFINE} -z 1"
            if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]];
              then
                basecall="${basecall} -x [${EXTRACTION_REGISTRATION_MASK}]"
              fi
            stage1="-m MI[${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0"
            stage2="-m MI[${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0"
            stage3="-m CC[${EXTRACTION_TEMPLATE},${N4_CORRECTED_IMAGES[0]},0.5,4] -m CC[${EXTRACTION_TEMPLATE_LAPLACIAN},${EXTRACTION_LAPLACIAN},0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0"

            exe_brain_extraction_1="${basecall} ${stage1} ${stage2} ${stage3}"
            logCmd $exe_brain_extraction_1
          fi

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

        ## Step 2 ##
        exe_brain_extraction_2="${WARP} -d ${DIMENSION} -i ${EXTRACTION_PRIOR} -o ${EXTRACTION_MASK_PRIOR_WARPED} -r ${ANATOMICAL_IMAGES[0]} -n Gaussian -t [${EXTRACTION_GENERIC_AFFINE},1] -t ${EXTRACTION_INVERSE_WARP}"
        logCmd $exe_brain_extraction_2

        ## superstep 1b ##
        logCmd ${ANTSPATH}ThresholdImage ${DIMENSION} ${EXTRACTION_MASK_PRIOR_WARPED} ${EXTRACTION_MASK_PRIOR_WARPED} 0.5 1 1 0
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK_PRIOR_WARPED} 2
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} GetLargestComponent ${EXTRACTION_MASK}

        ## superstep 6 ##
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE='';
        for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
          do
            ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${N4_CORRECTED_IMAGES[$i]}";
          done

        exe_brain_extraction_3="${ATROPOS} -d ${DIMENSION} -o ${EXTRACTION_SEGMENTATION} ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -x ${EXTRACTION_MASK} -i ${ATROPOS_BRAIN_EXTRACTION_INITIALIZATION} -c ${ATROPOS_BRAIN_EXTRACTION_CONVERGENCE} -m ${ATROPOS_BRAIN_EXTRACTION_MRF} -k ${ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD}"
        logCmd $exe_brain_extraction_3

        logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_WM} 3 3 1 0
        logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_GM} 2 2 1 0
        logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_CSF} 1 1 1 0

        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_WM} GetLargestComponent ${EXTRACTION_WM}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_GM} GetLargestComponent ${EXTRACTION_GM}

        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_TMP} FillHoles ${EXTRACTION_GM} 2
        logCmd ${ANTSPATH}MultiplyImages ${DIMENSION} ${EXTRACTION_GM} ${EXTRACTION_TMP} ${EXTRACTION_GM}

        logCmd ${ANTSPATH}MultiplyImages ${DIMENSION} ${EXTRACTION_WM} 3 ${EXTRACTION_WM}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_TMP} ME ${EXTRACTION_CSF} 10

        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_GM} addtozero ${EXTRACTION_GM} ${EXTRACTION_TMP}
        logCmd ${ANTSPATH}MultiplyImages ${DIMENSION} ${EXTRACTION_GM} 2 ${EXTRACTION_GM}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_SEGMENTATION} addtozero ${EXTRACTION_WM} ${EXTRACTION_GM}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_SEGMENTATION} addtozero ${EXTRACTION_SEGMENTATION} ${EXTRACTION_CSF}

        ## superstep 7 ##
        logCmd ${ANTSPATH}ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_MASK} 2 3
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} ME ${EXTRACTION_MASK} 2
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} GetLargestComponent ${EXTRACTION_MASK}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK} 4
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} FillHoles ${EXTRACTION_MASK} 2
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} addtozero ${EXTRACTION_MASK} ${EXTRACTION_MASK_PRIOR_WARPED}
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} MD ${EXTRACTION_MASK} 5
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${EXTRACTION_MASK} ME ${EXTRACTION_MASK} 5

      fi

    if [[ ! -f ${EXTRACTION_WM} ]];
      then
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE='';
        for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
          do
            ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${N4_CORRECTED_IMAGES[$i]}";
          done

        exe_brain_extraction_3="${ATROPOS} -d ${DIMENSION} -o ${EXTRACTION_SEGMENTATION} ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -x ${EXTRACTION_MASK} -i ${ATROPOS_BRAIN_EXTRACTION_INITIALIZATION} -c ${ATROPOS_BRAIN_EXTRACTION_CONVERGENCE} -m ${ATROPOS_BRAIN_EXTRACTION_MRF} -k ${ATROPOS_BRAIN_EXTRACTION_LIKELIHOOD}"
        logCmd $exe_brain_extraction_3

        logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${EXTRACTION_SEGMENTATION} ${EXTRACTION_WM} 3 3 1 0
      fi

    logCmd ${ANTSPATH}/MultiplyImages ${DIMENSION} ${EXTRACTION_MASK} ${N4_CORRECTED_IMAGES[0]} ${EXTRACTION_BRAIN}

    if [[ ! -f ${EXTRACTION_MASK} ]];
      then
        echo "Expected output was not produced.  The brain mask doesn't exist:"
        echo "   $EXTRACTION_MASK"
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

if [[ $KEEP_TMP_IMAGES = "false" || $KEEP_TMP_IMAGES = "0" ]];
  then
    for f in ${TMP_FILES[@]}
      do
        logCmd rm $f
      done
  fi

################################################################################
#
# Brain segmentation
#
################################################################################

SEGMENTATION_WARP=${SEGMENTATION_WARP_OUTPUT_PREFIX}1Warp.nii.gz
SEGMENTATION_INVERSE_WARP=${SEGMENTATION_WARP_OUTPUT_PREFIX}1InverseWarp.nii.gz
SEGMENTATION_GENERIC_AFFINE=${SEGMENTATION_WARP_OUTPUT_PREFIX}0GenericAffine.mat
SEGMENTATION_WHITE_MATTER_MASK=${EXTRACTION_WM}
SEGMENTATION_BRAIN=${EXTRACTION_BRAIN}
SEGMENTATION_MASK_DILATED=${BRAIN_SEGMENTATION_OUTPUT}MaskDilated.nii.gz
SEGMENTATION_BRAIN_WEIGHT_MASK=${BRAIN_SEGMENTATION_OUTPUT}WeightMask.nii.gz
SEGMENTATION_CONVERGENCE_FILE=${BRAIN_SEGMENTATION_OUTPUT}Convergence.txt
SEGMENTATION_PREVIOUS_ITERATION_LABELING=${BRAIN_SEGMENTATION_OUTPUT}PreviousIterationSegmentation.${OUTPUT_SUFFIX}

BRAIN_SEGMENTATION_POSTERIORS=${BRAIN_SEGMENTATION_OUTPUT}Posteriors%${FORMAT}d.${OUTPUT_SUFFIX}
POSTERIOR_IMAGE_FILENAMES=( ${BRAIN_SEGMENTATION_OUTPUT}Posteriors${CSF_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX} ${BRAIN_SEGMENTATION_OUTPUT}Posteriors${GRAY_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX} ${BRAIN_SEGMENTATION_OUTPUT}Posteriors${WHITE_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX} )

if [[ ! -f ${BRAIN_SEGMENTATION} ]];
  then

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Brain segmentation using the following steps:"
    echo "   1) Register $SEGMENTATION_TEMPLATE and $SEGMENTATION_PRIOR to ${N4_CORRECTED_IMAGES[0]}"
    echo "   2) Warp priors to ${N4_CORRECTED_IMAGES[0]}"
    echo "   3) N-tissue segmentation using Atropos and N4"
    echo "--------------------------------------------------------------------------------------"
    echo

    # Check inputs
    if [[ ! -f ${SEGMENTATION_TEMPLATE} ]];
      then
        echo "The segmentation template doesn't exist:"
        echo "   $SEGMENTATION_TEMPLATE"
        exit 1
      fi
    if [[ ! -f ${SEGMENTATION_BRAIN} ]];
      then
        echo "The extracted brain doesn't exist:"
        echo "   $SEGMENTATION_BRAIN"
        exit 1
      fi
    if [[ ${NUMBER_OF_PRIOR_IMAGES} -ne 3 ]];
      then
        echo "Expected 3 prior images (${NUMBER_OF_PRIOR_IMAGES} are specified).  Check the command line specification."
        exit 1
      fi

    time_start_brain_segmentation=`date +%s`

    ## Step 1 ##
    if [[ ! -f ${SEGMENTATION_WARP} ]];
      then
        logCmd ${ANTSPATH}ImageMath ${DIMENSION} ${SEGMENTATION_MASK_DILATED} MD ${EXTRACTION_MASK} 20

        basecall=''
        if [[ $INPUT_TEMPLATES_ARE_THE_SAME -eq 1 && -f ${EXTRACTION_AFFINE} ]];
          then
            basecall="${ANTS} -d ${DIMENSION} -u 1 -w [0.025,0.975] -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} -r [${EXTRACTION_AFFINE},1] -z 1"
          else
            basecall="${ANTS} -d ${DIMENSION} -u 1 -w [0.025,0.975] -o ${SEGMENTATION_WARP_OUTPUT_PREFIX} -r [${SEGMENTATION_BRAIN},${SEGMENTATION_TEMPLATE},1] -z 1"
          fi

        basecall="${basecall} -x [${SEGMENTATION_MASK_DILATED}]"
        stage1="-m MI[${SEGMENTATION_BRAIN},${SEGMENTATION_TEMPLATE},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0"
        stage2="-m MI[${SEGMENTATION_BRAIN},${SEGMENTATION_TEMPLATE},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0"
        stage3="-m CC[${SEGMENTATION_BRAIN},${SEGMENTATION_TEMPLATE},1,4] -c [${ANTS_MAX_ITERATIONS},1e-9,15] -t ${ANTS_TRANSFORMATION} -f 6x4x2x1 -s 3x2x1x0"

        exe_brain_segmentation_1="${basecall} ${stage1} ${stage2} ${stage3}"
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

        exe_brain_segmentation_2="${WARP} -d ${DIMENSION} -i ${PRIOR_IMAGE_FILENAMES[$i]} -o ${WARPED_PRIOR_IMAGE_FILENAMES[$i]} -r ${ANATOMICAL_IMAGES[0]} -n Gaussian  -t ${SEGMENTATION_WARP} -t ${SEGMENTATION_GENERIC_AFFINE}"
        logCmd $exe_brain_segmentation_2
      done

    ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE='';
    for (( j = 0; j < ${#ANATOMICAL_IMAGES[@]}; j++ ))
      do
        ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE="${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} -a ${ANATOMICAL_IMAGES[$j]}";
      done

    bash ${ANTSPATH}/antsAtroposN4.sh \
      -d ${DIMENSION} \
      ${ATROPOS_ANATOMICAL_IMAGES_COMMAND_LINE} \
      -x ${EXTRACTION_MASK} \
      -m ${ATROPOS_SEGMENTATION_NUMBER_OF_ITERATIONS} \
      -n 5 \
      -c 3 \
      -l 3 \
      -p ${SEGMENTATION_PRIOR_WARPED} \
      -k $KEEP_TMP_IMAGES \
      -o ${OUTPUT_PREFIX}Brain \
      -s ${OUTPUT_SUFFIX}

    ## Step 3 ###
    TMP_FILES=( $EXTRACTION_AFFINE $SEGMENTATION_WARP $SEGMENTATION_INVERSE_WARP $SEGMENTATION_GENERIC_AFFINE $SEGMENTATION_WHITE_MATTER_MASK $SEGMENTATION_BRAIN $SEGMENTATION_MASK_DILATED )
    TMP_FILES=( ${TMP_FILES[@]} ${WARPED_PRIOR_IMAGE_FILENAMES[@]} $SEGMENTATION_PREVIOUS_ITERATION_LABELING $EXTRACTION_GENERIC_AFFINE $SEGMENTATION_BRAIN_WEIGHT_MASK)

    if [[ $KEEP_TMP_IMAGES = "false" || $KEEP_TMP_IMAGES = "0" ]];
      then
        for f in ${TMP_FILES[@]}
          do
            logCmd rm $f
          done
      fi

     time_end_brain_segmentation=`date +%s`
     time_elapsed_brain_segmentation=$((time_end_brain_segmentation - time_start_brain_segmentation))

     echo
     echo "--------------------------------------------------------------------------------------"
     echo " Done with brain segmentation:  $(( time_elapsed_brain_segmentation / 3600 ))h $(( time_elapsed_brain_segmentation %3600 / 60 ))m $(( time_elapsed_brain_segmentation % 60 ))s"
     echo "--------------------------------------------------------------------------------------"
     echo

   fi

################################################################################
#
# Cortical thickness
#
################################################################################

CORTICAL_THICKNESS_GM=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${GRAY_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}
CORTICAL_THICKNESS_WM=${BRAIN_SEGMENTATION_OUTPUT}Posteriors${WHITE_MATTER_LABEL_FORMAT}.${OUTPUT_SUFFIX}

if [[ ! -f ${CORTICAL_THICKNESS_IMAGE} ]];
  then

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Cortical thickness using DiReCT (KellyKapowski)"
    echo "--------------------------------------------------------------------------------------"
    echo

    # Check inputs
    if [[ ! -f ${BRAIN_SEGMENTATION} ]];
      then
        echo "The brain segmentation image doesn't exist:"
        echo "   $BRAIN_SEGMENTATION"
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

    TMP_FILES=()

    exe_direct="${DIRECT} -d ${DIMENSION} -s [${BRAIN_SEGMENTATION},${GRAY_MATTER_LABEL},${WHITE_MATTER_LABEL}] -g ${CORTICAL_THICKNESS_GM} -w ${CORTICAL_THICKNESS_WM} -o ${CORTICAL_THICKNESS_IMAGE} -c ${DIRECT_CONVERGENCE} -t ${DIRECT_THICKNESS_PRIOR} -r ${DIRECT_GRAD_STEP_SIZE} -m ${DIRECT_SMOOTHING_SIGMA} -n ${DIRECT_NUMBER_OF_DIFF_COMPOSITIONS}"
    logCmd $exe_direct

    if [[ $KEEP_TMP_IMAGES = "false" || $KEEP_TMP_IMAGES = "0" ]];
      then
        for f in ${TMP_FILES[@]}
          do
            logCmd rm $f
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

################################################################################
#
# Registration to a template
#
################################################################################

if [[ -f ${REGISTRATION_TEMPLATE} ]];
  then

    REGISTRATION_TEMPLATE_OUTPUT_PREFIX=${OUTPUT_PREFIX}RegistrationToTemplate
    REGISTRATION_TEMPLATE_GENERIC_AFFINE=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}0GenericAffine.mat
    REGISTRATION_TEMPLATE_WARP=${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}1Warp.${OUTPUT_SUFFIX}

    if [[ ! -f ${REGISTRATION_TEMPLATE_GENERIC_AFFINE} || ! -f ${REGISTRATION_TEMPLATE_WARP} || ! -f "${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}Warped.${OUTPUT_SUFFIX}" ]];
      then

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " T1 registration to specified template"
        echo "--------------------------------------------------------------------------------------"
        echo

        TMP_FILES=()

        time_start_template_registration=`date +%s`

        basecall="${ANTS} -d ${DIMENSION} -u 1 -w [0.025,0.975] -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX} -r [${REGISTRATION_TEMPLATE},${SEGMENTATION_BRAIN_N4_IMAGES[0]},1] -z 1"
        stage1="-m MI[${REGISTRATION_TEMPLATE},${SEGMENTATION_BRAIN_N4_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Rigid[0.1] -f 4x2x1 -s 2x1x0"
        stage2="-m MI[${REGISTRATION_TEMPLATE},${SEGMENTATION_BRAIN_N4_IMAGES[0]},${ANTS_LINEAR_METRIC_PARAMS}] -c ${ANTS_LINEAR_CONVERGENCE} -t Affine[0.1] -f 4x2x1 -s 2x1x0"
        stage3="-m CC[${REGISTRATION_TEMPLATE},${SEGMENTATION_BRAIN_N4_IMAGES[0]},1,4] -c [${ANTS_MAX_ITERATIONS},1e-9,15] -t ${ANTS_TRANSFORMATION} -f 4x2x1 -s 2x1x0"

        exe_template_registration_1="${basecall} ${stage1} ${stage2} ${stage3}"

        if [[ ! -f ${REGISTRATION_TEMPLATE_WARP} ]];
          then
            logCmd $exe_template_registration_1
          fi

        exe_template_registration_2="${WARP} -d ${DIMENSION} -i ${SEGMENTATION_BRAIN_N4_IMAGES[0]} -o ${REGISTRATION_TEMPLATE_OUTPUT_PREFIX}Warped.${OUTPUT_SUFFIX} -r ${REGISTRATION_TEMPLATE} -n Gaussian -t ${REGISTRATION_TEMPLATE_WARP} -t ${REGISTRATION_TEMPLATE_GENERIC_AFFINE}"
        logCmd $exe_template_registration_2

        if [[ $KEEP_TMP_IMAGES = "false" || $KEEP_TMP_IMAGES = "0" ]];
          then
            for f in ${TMP_FILES[@]}
              do
                logCmd rm $f
              done
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

        time_end_template_registration=`date +%s`
        time_elapsed_template_registration=$((time_end_template_registration - time_start_template_registration))

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Done with registration to specified template:  $(( time_elapsed_template_registration / 3600 ))h $(( time_elapsed_template_registration %3600 / 60 ))m $(( time_elapsed_template_registration % 60 ))s"
        echo "--------------------------------------------------------------------------------------"
        echo

      fi
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
echo " Done with ANTs processing pipeline"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0

