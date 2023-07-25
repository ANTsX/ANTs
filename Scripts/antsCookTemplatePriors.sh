#!/bin/bash

VERSION="0.0"

# Check dependencies

PROGRAM_DEPENDENCIES=( 'ImageMath' 'SmoothImage' 'ThresholdImage' 'antsRegistration' 'antsApplyTransforms' 'N4BiasFieldCorrection' 'Atropos' 'KellyKapowski' )
SCRIPTS_DEPENDENCIES=( 'antsCorticalThickness.sh' 'antsBrainExtraction.sh' 'antsAtroposN4.sh' 'antsJointLabelFusion.sh' )

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

`basename $0` intializes tissue priors for a brain template.  The following steps
are performed:
  1. Run the template through antsCorticalThickness
  2. Create priors for the template
     a. The brain extraction prior is created by smoothing the brain extraction
        mask created during 1.
     b. If labeled atlases are not provided, we smooth the posteriors from 2a to create
        the new segmentation priors, otherwise we use antsJointLabelFusion to create a set of
        posteriors (https://github.com/ntustison/antsCookTemplatePriorsExample).

Usage:

`basename $0` -d imageDimension
              -e brainTemplate
              -m brainExtractionProbabilityMask
              -p brainSegmentationPriors
              <OPTARGS>
              -o outputPrefix
              \${templateImages[@]}

Example:

  bash $0 -d 3 -e brainWithSkullTemplate.nii.gz -m brainPrior.nii.gz -p segmentationPriors%d.nii.gz -o output \${templateImages[@]}

Required arguments:

     -d:  Image dimension                       2 or 3 (for 2- or 3-dimensional image)
     -e:  Brain reference template              Anatomical *intensity* template (possibly created using a population
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
     -o:  Output prefix                         The following subdirectory and images are created for the single
                                                subject template
                                                  * \${OUTPUT_PREFIX}SingleSubjectTemplate/
                                                  * \${OUTPUT_PREFIX}SingleSubjectTemplate/T_template*.nii.gz

     anatomical images                          Set of multimodal input data assumed to be specified ordered as
                                                follows:
                                                  \${time1_modality1} \${time1_modality2} ... \${time1_modalityN} \\
                                                  \${time2_modality1} \${time2_modality2} ...
                                                  .
                                                  .
                                                  .
                                                  \${timeN_modality1} ...

						A single modality is expected by default, in which case the input images
						are simply ordered by time:

						  \${time1_modality1} \${time2_modality1} ... \${timeN_modality1}

						If there are multiple modalities, use the -k option to specify how many.

Optional arguments:

     -s:  image file suffix                     Any of the standard ITK IO formats e.g. nrrd, nii.gz (default), mhd
     -c:  control type                          Control for parallel computation (default 0):
                                                  0 = run serially
                                                  1 = SGE qsub
                                                  2 = use PEXEC (localhost)
                                                  3 = Apple XGrid
                                                  4 = PBS qsub
                                                  5 = SLURM
     -a:                                        Atlases (assumed to be skull-stripped) used to cook template priors.  If atlases
                                                aren't used then we simply smooth the single-subject template posteriors after
                                                passing through antsCorticalThickness.sh. Example:

 						  -a atlas1.nii.gz -a atlas2.nii.gz ... -a atlasN.nii.gz

     -l:                                        Labels associated with each atlas, in the same order as they are specified
						                                    with the -a option. The number of labels in each image is assumed to be equal
                                                to the number of priors.
     -f:  extraction registration mask          Mask (defined in the template space) used during registration
                                                for brain extraction.
     -g:  denoise anatomical images             Denoise anatomical images (default = 0).
     -j:  number of cpu cores                   Number of cpu cores to use locally for pexec option (default 2; requires "-c 2")
     -k:  number of modalities                  Number of modalities used to construct the template (default 1):  For example,
                                                if one wanted to use multiple modalities consisting of T1, T2, and FA
                                                components ("-k 3").
     -u:  use floating-point precision          Use floating point precision in registrations (default = 0)
     -v:  Atropos segmentation weight           Atropos spatial prior *probability* weight for the segmentation for the
                                                template (default = 0.25)
     -w:  Atropos segmentation weight           Atropos spatial prior *probability* weight for the segmentations (default = 0.5)
     -q:  Use quick registration parameters     If 'yes' then we use antsRegistrationSyNQuick.sh as the basis for registration.
                                                Otherwise use antsRegistrationSyN.sh.  The options are as follows:
                                                '-q 0' = antsRegistrationSyN for everything (default)
                                                '-q 1' = Fast antsCorticalThickness to SST
                                                '-q 2' = Fast JLF cooking
                                                '-q 3' = Fast everything
     -z:  Test / debug mode                     If > 0, runs a faster version of the script. Only for testing. Implies -u 0
                                                in the antsCorticalThickness.sh script (i.e., no random seeding).
                                                Requires single thread computation for complete reproducibility.
USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using antsLongitudinalCorticalThickness with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      brain template          = ${BRAIN_TEMPLATE}
      extraction prior        = ${EXTRACTION_PRIOR}
      segmentation prior      = ${SEGMENTATION_PRIOR}
      output prefix           = ${OUTPUT_PREFIX}
      output image suffix     = ${OUTPUT_SUFFIX}
      registration template   = ${REGISTRATION_TEMPLATE}

    Other parameters:
      run quick               = ${RUN_QUICK}
      debug mode              = ${DEBUG_MODE}
      float precision         = ${USE_FLOAT_PRECISION}
      denoise                 = ${DENOISE}
      use random seeding      = ${USE_RANDOM_SEEDING}
      number of modalities    = ${NUMBER_OF_MODALITIES}
      number of cores         = ${CORES}
      control type            = ${DOQSUB}
      rigid alignment to SST  = ${RIGID_ALIGNMENT_TO_SST}

PARAMETERS
}

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


HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"

DIMENSION=3

NUMBER_OF_MODALITIES=1

TEMPLATE_IMAGES=()
RUN_QUICK=1
USE_RANDOM_SEEDING=1

BRAIN_TEMPLATE=""
EXTRACTION_PRIOR=""
EXTRACTION_REGISTRATION_MASK=""
SEGMENTATION_PRIOR=""
USE_SST_CORTICAL_THICKNESS_PRIOR=0
REGISTRATION_TEMPLATE=""
DO_REGISTRATION_TO_TEMPLATE=0
DENOISE=0

ATROPOS_SEGMENTATION_PRIOR_WEIGHT_SST=0.25
ATROPOS_SEGMENTATION_PRIOR_WEIGHT_TIMEPOINT=0.5

DOQSUB=0
CORES=2
RIGID_ALIGNMENT_TO_SST=0

MALF_ATLASES=()
MALF_LABELS=()
MALF_LABEL_STRINGS_FOR_PRIORS=()

################################################################################
#
# Programs and their parameters
#
################################################################################

USE_FLOAT_PRECISION=0

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:c:d:e:f:g:h:j:k:l:m:o:p:q:r:s:t:u:v:w:z:" OPT
    do
      case $OPT in
          a)
       MALF_ATLASES[${#MALF_ATLASES[@]}]=$OPTARG
       ;;
          b) # posterior formulation
       ATROPOS_SEGMENTATION_POSTERIOR_FORMULATION=$OPTARG
       ;;
          c)
       DOQSUB=$OPTARG
       if [[ $DOQSUB -gt 5 ]];
         then
           echo " DOQSUB must be an integer value (0=serial, 1=SGE qsub, 2=try pexec, 3=XGrid, 4=PBS qsub, 5=SLURM ) you passed  -c $DOQSUB "
           exit 1
         fi
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
          g) #denoise
       DENOISE=$OPTARG
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          j) #number of cpu cores to use (default = 2)
       CORES=$OPTARG
       ;;
          k) #number of modalities
       NUMBER_OF_MODALITIES=$OPTARG
       ;;
          l)
       MALF_LABELS[${#MALF_LABELS[@]}]=$OPTARG
       ;;
          m) #brain extraction prior probability mask
       EXTRACTION_PRIOR=$OPTARG
       ;;
          n) # use
       USE_SST_CORTICAL_THICKNESS_PRIOR=$OPTARG
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
          u) #use floating point precision
       USE_FLOAT_PRECISION=$OPTARG
       ;;
          v) #atropos prior weight for single subject template
       ATROPOS_SEGMENTATION_PRIOR_WEIGHT_SST=$OPTARG
       ;;
          w) #atropos prior weight for each individual time point
       ATROPOS_SEGMENTATION_PRIOR_WEIGHT_TIMEPOINT=$OPTARG
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

NUMBER_OF_PRIOR_IMAGES=${#WARPED_PRIOR_IMAGE_FILENAMES[*]}

# Shiftsize is calculated because a variable amount of arguments can be used on the command line.
# The shiftsize variable will give the correct number of arguments to skip. Issuing shift $shiftsize will
# result in skipping that number of arguments on the command line, so that only the input images remain.
shiftsize=$(($OPTIND - 1))
shift $shiftsize
# The invocation of $* will now read all remaining arguments into the variable IMAGESETVARIABLE
IMAGESETVARIABLE=$*
NINFILES=$(($nargs - $shiftsize))
IMAGESETARRAY=()

for IMG in $IMAGESETVARIABLE
  do
    TEMPLATE_IMAGES[${#TEMPLATE_IMAGES[@]}]=$IMG
  done

if [[ ${#TEMPLATE_IMAGES[@]} -eq 0 ]];
  then
    echo "Error:  no template images specified."
    exit 1
  fi

if [[ $NUMBER_OF_MODALITIES -gt 1 ]];
  then
    echo "--------------------------------------------------------------------------------------"
    echo " Cortical thickness using the following ${NUMBER_OF_MODALITIES}-tuples:  "
    echo "--------------------------------------------------------------------------------------"
    for (( i = 0; i < ${#TEMPLATE_IMAGES[@]}; i+=$NUMBER_OF_MODALITIES ))
      do
        IMAGEMETRICSET=""
        for (( j = 0; j < $TEMPLATE_IMAGES; j++ ))
          do
            k=0
            let k=$i+$j
            IMAGEMETRICSET="$IMAGEMETRICSET ${TEMPLATE_IMAGES[$k]}"
          done
        echo $IMAGEMETRICSET
      done
    echo "--------------------------------------------------------------------------------------"
fi

if [[ ${#MALF_ATLASES[@]} -ne ${#MALF_LABELS[@]} ]]
  then
    echo "Error:  The number of malf atlases and labels aren't equal."
  fi


# Set up various things related to RUN_QUICK

# Can't do everything fast and still get good results if there is large deformation.
# Initiate levels of fast:

# 0 - Fast SST (old ANTS) but everything else slower for quality
# 1 - + Fast antsct to SST
# 2 - + Fast MALF cooking
# 3 - + Fast everything

RUN_OLD_ANTS_SST_CREATION=1
RUN_ANTSCT_TO_SST_QUICK=0
RUN_FAST_MALF_COOKING=0
RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE=0

if [[ $RUN_QUICK -gt 0 ]];
  then
    RUN_ANTSCT_TO_SST_QUICK=1
  fi

if [[ $RUN_QUICK -gt 1 ]];
  then
    RUN_FAST_MALF_COOKING=1
  fi

if [[ $RUN_QUICK -gt 2 ]];
  then
    RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE=1
  fi

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#  3. See if $REGISTRATION_TEMPLATE is the same as $BRAIN_TEMPLATE
#
################################################################################

for (( i = 0; i < ${#TEMPLATE_IMAGES[@]}; i++ ))
  do
  if [[ ! -f ${TEMPLATE_IMAGES[$i]} ]];
    then
      echo "The specified image \"${TEMPLATE_IMAGES[$i]}\" does not exist."
      exit 1
    fi
  done





################################################################################
#
#  Run template through antsCorticalThickness.sh
#
################################################################################

TEMPLATE_EXTRACTION_MASK=${OUTPUT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
TEMPLATE_EXTRACTION_REGISTRATION_MASK=${OUTPUT_PREFIX}BrainExtractionRegistrationMask.${OUTPUT_SUFFIX}
TEMPLATE_PRIOR=${OUTPUT_PREFIX}Priors\%${FORMAT}d.${OUTPUT_SUFFIX}
TEMPLATE_EXTRACTION_PRIOR=${OUTPUT_PREFIX}BrainExtractionMaskPrior.${OUTPUT_SUFFIX}
TEMPLATE_CORTICAL_THICKNESS=${OUTPUT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}
TEMPLATE_SKULL_STRIPPED=${OUTPUT_PREFIX}BrainExtractionBrain.${OUTPUT_SUFFIX}

echo
echo "--------------------------------------------------------------------------------------"
echo " Creating template priors:  running template through antsCorticalThickness                 "
echo "--------------------------------------------------------------------------------------"
echo

time_start_priors=`date +%s`

TEMPLATE_IMAGES_LIST=''
for (( j=0; j < $NUMBER_OF_MODALITIES; j++ ))
  do
    TEMPLATE_IMAGES_LIST="${TEMPLATE_IMAGES_LIST} -a ${TEMPLATE_IMAGES[$j]}"
  done

REG_MASK=""
if [[ -f ${EXTRACTION_REGISTRATION_MASK} ]]; then
  REG_MASK="-f ${EXTRACTION_REGISTRATION_MASK}"
fi

if [[ ! -f ${TEMPLATE_CORTICAL_THICKNESS} ]];
  then
    logCmd antsCorticalThickness.sh \
      -d ${DIMENSION} \
      -q ${RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE} \
      $TEMPLATE_IMAGES_LIST \
      $REG_MASK \
      -e ${BRAIN_TEMPLATE} \
      -g ${DENOISE} \
      -m ${EXTRACTION_PRIOR} \
      -k 0 \
      -z ${DEBUG_MODE} \
      -p ${SEGMENTATION_PRIOR} \
      -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_SST} \
      -o ${OUTPUT_PREFIX}
  fi

TEMPLATE_POSTERIORS=( ${OUTPUT_PREFIX}BrainSegmentationPosteriors*.${OUTPUT_SUFFIX} )
TEMPLATE_POSTERIORS_EXIST=1
TEMPLATE_PRIORS_EXIST=1

TEMPLATE_PRIORS=()
for (( j = 0; j < ${#TEMPLATE_POSTERIORS[@]}; j++ ))
  do
    POSTERIOR=${TEMPLATE_POSTERIORS[$j]}

    if [[ ! -f ${POSTERIOR} ]];
      then
        TEMPLATE_POSTERIORS_EXIST=0
        TEMPLATE_PRIORS_EXIST=0
        break;
      fi

    TEMPLATE_PRIORS[$j]=${POSTERIOR/BrainSegmentationPosteriors/Priors}

    if [[ ! -f ${TEMPLATE_PRIORS[$j]} ]];
      then
        TEMPLATE_PRIORS_EXIST=0
      fi
  done

if [[ ${TEMPLATE_POSTERIORS_EXIST} -eq 0 ]];
  then
    echo "Error:  Posteriors for the template do not exist."
    exit 1
  fi

logCmd ImageMath ${DIMENSION} ${TEMPLATE_SKULL_STRIPPED} m ${TEMPLATE_IMAGES[0]} ${TEMPLATE_EXTRACTION_MASK}
logCmd SmoothImage ${DIMENSION} ${TEMPLATE_EXTRACTION_MASK} 1 ${TEMPLATE_EXTRACTION_PRIOR} 1
logCmd ImageMath ${DIMENSION} ${TEMPLATE_EXTRACTION_REGISTRATION_MASK} MD ${TEMPLATE_EXTRACTION_MASK} 40

if [[ ${TEMPLATE_PRIORS_EXIST} -eq 0 ]];
  then
    if [[ ${#MALF_ATLASES[@]} -eq 0 ]];
      then

        echo
        echo "   ---> Smoothing template posteriors as priors."
        echo

        for j in ${TEMPLATE_POSTERIORS[@]}
          do
            PRIOR=${j/BrainSegmentationPosteriors/Priors}
            logCmd SmoothImage ${DIMENSION} $j 1 $PRIOR 1
          done

      else

        echo
        echo "   ---> Cooking template priors using antsJointLabelFusion."
        echo

        TEMPLATE_MALF_LABELS_PREFIX=${OUTPUT_PREFIX}
        TEMPLATE_MALF_LABELS=${OUTPUT_PREFIX}Labels.nii.gz

        ATLAS_AND_LABELS_STRING=''
        for (( j=0; j < ${#MALF_ATLASES[@]}; j++ ))
          do
            ATLAS_AND_LABELS_STRING="${ATLAS_AND_LABELS_STRING} -g ${MALF_ATLASES[$j]} -l ${MALF_LABELS[$j]}"
          done

        if [[ ! -f ${TEMPLATE_MALF_LABELS} ]];
          then
            logCmd antsJointLabelFusion.sh \
              -d ${DIMENSION} \
              -q ${RUN_FAST_MALF_COOKING} \
              -x ${TEMPLATE_EXTRACTION_MASK} \
              -c ${DOQSUB} \
              -j ${CORES} \
              -t ${TEMPLATE_SKULL_STRIPPED} \
              -o ${TEMPLATE_MALF_LABELS_PREFIX} \
              ${ATLAS_AND_LABELS_STRING}
          fi

        TEMPLATE_PRIORS=()
        for (( j = 0; j < ${#TEMPLATE_POSTERIORS[@]}; j++ ))
          do
            POSTERIOR=${TEMPLATE_POSTERIORS[$j]}

            TEMPLATE_PRIORS[$j]=${POSTERIOR/BrainSegmentationPosteriors/Priors}

            let PRIOR_LABEL=$j+1
            logCmd ThresholdImage ${DIMENSION} ${TEMPLATE_MALF_LABELS} ${TEMPLATE_PRIORS[$j]} ${PRIOR_LABEL} ${PRIOR_LABEL} 1 0
            logCmd SmoothImage ${DIMENSION} ${TEMPLATE_PRIORS[$j]} 1 ${TEMPLATE_PRIORS[$j]} 1
          done

        TMP_CSF_POSTERIOR=${OUTPUT_PREFIX}BrainSegmentationCsfPosteriorTmp.${OUTPUT_SUFFIX}
        logCmd SmoothImage ${DIMENSION} ${TEMPLATE_POSTERIORS[0]} 1 ${TMP_CSF_POSTERIOR} 1
        logCmd ImageMath ${DIMENSION} ${TEMPLATE_PRIORS[0]} max ${TEMPLATE_PRIORS[0]} ${TMP_CSF_POSTERIOR}
        # Clip priors to remove precision errors
        logCmd ImageMath ${DIMENSION} ${TEMPLATE_PRIORS[0]} WindowImage ${TEMPLATE_PRIORS[0]} 0 1 0 1

        # Brian's finishing touches on "cooking"---subtract out CSF from all other priors
        for (( j = 1; j < ${#TEMPLATE_PRIORS[@]}; j++ ))
          do
            let PRIOR_LABEL=$j+1

            logCmd ImageMath ${DIMENSION} ${TEMPLATE_PRIORS[$j]} - ${TEMPLATE_PRIORS[$j]} ${TEMPLATE_PRIORS[0]}
            # Clip priors to range [0,1]
            logCmd ImageMath ${DIMENSION} ${TEMPLATE_PRIORS[$j]} WindowImage ${TEMPLATE_PRIORS[$j]} 0 1 0 1
          done

        logCmd rm -f $TMP_CSF_POSTERIOR
        logCmd rm -f ${TEMPLATE_MALF_LABELS_PREFIX}*log.txt
      fi
  fi

time_end_priors=`date +%s`
time_elapsed_priors=$((time_end_priors - time_start_priors))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with creating template priors:  $(( time_elapsed_priors / 3600 ))h $(( time_elapsed_priors %3600 / 60 ))m $(( time_elapsed_priors % 60 ))s"
echo "--------------------------------------------------------------------------------------"
echo
