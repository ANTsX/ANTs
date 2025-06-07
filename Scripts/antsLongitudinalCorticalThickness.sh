#!/bin/bash

VERSION="0.0"

# Check dependencies

PROGRAM_DEPENDENCIES=( 'antsRegistration' 'antsApplyTransforms' 'N4BiasFieldCorrection' 'Atropos' 'KellyKapowski' )
SCRIPTS_DEPENDENCIES=( 'antsBrainExtraction.sh' 'antsAtroposN4.sh' 'antsMultivariateTemplateConstruction2.sh' 'antsJointLabelFusion.sh' )

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

`basename $0` performs a longitudinal cortical thickness estimation.  The following steps
are performed:
  1. Create a single-subject template (SST) from all the anatomical images.
  2. Create priors for the SST
     a. Run the SST through the individual cortical thickness pipeline using the group
        template and associated priors.
     b. The brain extraction SST prior is created by smoothing the brain extraction
        mask created during 2a.
     c. If labeled atlases are not provided, we smooth the posteriors from 2a to create
        the SST segmentation priors, otherwise we use antsJointFusion to create a set of
        posteriors (https://github.com/ntustison/antsCookTemplatePriorsExample).
  3. Using the SST + priors, we run each time point image through the antsCorticalThickness
     pipeline.
  4. Compose (time point) -> SST and SST -> (group template) warps to deform images from the
     native anatomical space to the group template.

Usage:

`basename $0` -d imageDimension
              -e brainTemplate
              -m brainExtractionProbabilityMask
              -p brainSegmentationPriors
              <OPTARGS>
              -o outputPrefix
              \${anatomicalImages[@]}

Example:

  bash $0 -d 3 -e brainWithSkullTemplate.nii.gz -m brainPrior.nii.gz -p segmentationPriors%d.nii.gz -o output \${anatomicalImages[@]}

Required arguments:

     -d:  Image dimension                       2 or 3 (for 2- or 3-dimensional image)

     -e:  Brain segmentation template           Anatomical *intensity* template. This template is *not* skull-stripped.
                                                The following images must be in the same space as this template:
                                                    * Brain probability mask (-m)
                                                    * Segmentation priors (-p).
                                                If used, the following optional images must also be in the same space as
                                                this template:
                                                    * Registration metric mask (-f).

     -m:  Brain extraction probability mask     Brain *probability* mask in the segmentation template space. A binary mask
                                                is an acceptable probability image.

     -p:  Brain segmentation priors             Tissue *probability* priors corresponding to the image specified
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

     -o:  Output prefix                         The following subdirectory and images are created for the single
                                                subject template
                                                  * \${OUTPUT_PREFIX}SingleSubjectTemplate/
                                                  * \${OUTPUT_PREFIX}SingleSubjectTemplate/T_template*.nii.gz
                                                A subdirectory is created for each anatomical image.

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
                                                  5 = SLURM.

     -t:  template for t1 registration          Anatomical *intensity* template. This template *must* be skull-stripped.
                                                This template is used to produce a final, high-quality registration between
                                                the skull-stripped single-subject template and this template.

     -a:                                        Atlases (assumed to be skull-stripped) used to cook template priors. If atlases
                                                aren't used then we simply smooth the single-subject template posteriors after
                                                passing through antsCorticalThickness.sh. Example:

                                                  -a atlas1.nii.gz -a atlas2.nii.gz ... -a atlasN.nii.gz

     -l:                                        Labels associated with each atlas, in the same order as they are specified
                                                with the -a option. The number of labels in each image is assumed to be equal
                                                to the number of priors.

     -f:  extraction registration mask          Binary metric mask defined in the segmentation template space (-e). During the
                                                registration for brain extraction, the similarity metric is only computed within
                                                this mask.

     -g:  denoise anatomical images             Denoise anatomical images (default = 0).

     -j:  number of cpu cores                   Number of cpu cores to use locally for pexec option (default 2; requires "-c 2").

     -k:  number of modalities                  Number of modalities used to construct the template (default 1):  For example,
                                                if one wanted to use multiple modalities consisting of T1, T2, and FA
                                                components ("-k 3").

     -n:  use SST cortical thickness prior      If set to '1', the cortical thickness map from the single-subject template is used
                                                as a prior constraint for each of the individual calls to antsCorticalThickness.sh
                                                (default = 0).

     -u:  use single float precision            If 1, use single float precision in registrations (default = 0).

     -v:  Atropos segmentation weight (SST)     Atropos spatial prior *probability* weight for the segmentation of the
                                                single-subject template (default = 0.25).

     -w:  Atropos segmentation weight (Indiv.)  Atropos spatial prior *probability* weight for the segmentation of the individual
                                                time points (default = 0.5).

     -x:                                        Number of iterations within Atropos (default 5).

     -q:  Use quick registration parameters     Use antsRegistrationSyNQuick.sh for some or all of the registrations.
                                                The higher the number, the more registrations are performed quickly. The options are
                                                as follows:
                                                  '-q 0' = antsRegistrationSyN.sh for everything (default)
                                                  '-q 1' = Quick registration of time points to the SST
                                                  '-q 2' = Adds quick registrations for prior cooking
                                                  '-q 3' = Quick registrations throughout.

     -r:  rigid alignment to SST                If 1, register anatomical images to the single-subject template before processing
                                                with antsCorticalThickness. This potentially reduces bias caused by variable
                                                orientation and voxel spacing (default = 0).

     -b:  keep temporary files                  Keep brain extraction/segmentation warps, etc (default = 0).

     -y:  rigid template update component       Update the single-subject template with the full affine transform (default 0).
                                                If 1, the rigid component of the affine transform will be used to update the
                                                template. Using the rigid component is desireable to reduce bias, but variations
                                                in the origin or head position across time points can cause the template head to
                                                drift out of the field of view.

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
      average rigid component = ${AFFINE_UPDATE_FULL}

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

DIMENSION=3

NUMBER_OF_MODALITIES=1

ANATOMICAL_IMAGES=()
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
ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=5

AFFINE_UPDATE_FULL=0

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
KEEP_TMP_IMAGES=0

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:b:c:d:e:f:g:h:j:k:l:m:n:o:p:q:r:s:t:u:v:x:w:y:z:" OPT
    do
      case $OPT in
          a)
       MALF_ATLASES[${#MALF_ATLASES[@]}]=$OPTARG
       ;;
          b)
       KEEP_TMP_IMAGES=$OPTARG
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
          x) #atropos segmentation internal iterations
       ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS=$OPTARG
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          p) #brain segmentation label prior image
       SEGMENTATION_PRIOR=$OPTARG
       ;;
          r) #rigid alignment to SST
       RIGID_ALIGNMENT_TO_SST=$OPTARG
       ;;
          t) #template registration image
       REGISTRATION_TEMPLATE=$OPTARG
       DO_REGISTRATION_TO_TEMPLATE=1
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
          y) # 1 update with full affine, 0 for no rigid (default = 0)
       AFFINE_UPDATE_FULL=$OPTARG
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

# This is the numeric format used by template priors (if formatted numerically)
# and the numeric format used for warped priors / posteriors in antsCorticalThickness.sh.
#
# If template priors are named using BIDS convention, output priors and posteriors will use
# this format.
#
# If template priors are named numerically with a different format (eg, %d or $03d), then
# PRIOR_NUMERIC_FORMAT will be set to match template priors
#
PRIOR_NUMERIC_FORMAT=""

if [[ "${SEGMENTATION_PRIOR}" =~ %([0-9]*)d ]]; then
  # Numeric formatting - normally 1 through 6, though only the first four are mandatory
  PRIOR_NUMERIC_FORMAT="%${BASH_REMATCH[1]}d"
elif [[ "${SEGMENTATION_PRIOR}" =~ label-%s ]]; then
  # antsCorticalThickness.sh will use this format if the priors are named using BIDS convention
  PRIOR_NUMERIC_FORMAT="%02d"
else
  echo "Invalid prior specification, priors must use c-style printf formatting, see usage"
  exit 1
fi

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
    ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$IMG
  done

if [[ ${#ANATOMICAL_IMAGES[@]} -eq 0 ]];
  then
    echo "Error:  no anatomical images specified."
    exit 1
  fi

if [[ $NUMBER_OF_MODALITIES -gt 1 ]];
  then
    echo "--------------------------------------------------------------------------------------"
    echo " Cortical thickness using the following ${NUMBER_OF_MODALITIES}-tuples:  "
    echo "--------------------------------------------------------------------------------------"
    for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i+=$NUMBER_OF_MODALITIES ))
      do
        IMAGEMETRICSET=""
        for (( j = 0; j < $ANATOMICAL_IMAGES; j++ ))
          do
            k=0
            let k=$i+$j
            IMAGEMETRICSET="$IMAGEMETRICSET ${ANATOMICAL_IMAGES[$k]}"
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

for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do
  if [[ ! -f ${ANATOMICAL_IMAGES[$i]} ]];
    then
      echo "The specified image \"${ANATOMICAL_IMAGES[$i]}\" does not exist."
      exit 1
    fi
  done


if [[ ${#ANATOMICAL_IMAGES[@]} -eq ${NUMBER_OF_MODALITIES} ]];
  then
    #echo "Only one set of anatomical images.  Running antsCorticalThickness.sh on the single subject."

    #SUBJECT_ANATOMICAL_IMAGES=''
    #for (( j=0; j < $NUMBER_OF_MODALITIES; j++ ))
    # do
    #   SUBJECT_ANATOMICAL_IMAGES="${SUBJECT_ANATOMICAL_IMAGES} -a ${ANATOMICAL_IMAGES[$j]}"
    # done

    # Won't be quick unless -q 3 was specified
    # But if you are running a longitudinal script without longitudinal data, that may not be the only problem
    #if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 1 ]];
    #    then
    #      logCmd antsCorticalThickness.sh \
    #        -d ${DIMENSION} \
    #        -x ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
    #        -t ${REGISTRATION_TEMPLATE} \
    #        -q ${RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE} \
    #        -e ${BRAIN_TEMPLATE} \
    #        -f ${EXTRACTION_REGISTRATION_MASK} \
    #        -m ${EXTRACTION_PRIOR} \
    #        -k ${KEEP_TMP_IMAGES} \
    #        -g ${DENOISE} \
    #        -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_TIMEPOINT} \
    #        -z ${DEBUG_MODE} \
    #        -p ${SEGMENTATION_PRIOR} \
    #        -o ${OUTPUT_PREFIX} \
    #        ${SUBJECT_ANATOMICAL_IMAGES}
    #  fi

    #if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 0 ]];
    #  then
    #     logCmd antsCorticalThickness.sh \
    #       -d ${DIMENSION} \
    #       -x ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
    #       -q ${RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE} \
    #       -e ${BRAIN_TEMPLATE} \
    #       -f ${EXTRACTION_REGISTRATION_MASK} \
    #       -m ${EXTRACTION_PRIOR} \
    #       -k ${KEEP_TMP_IMAGES} \
    #       -g ${DENOISE} \
    #       -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_TIMEPOINT} \
    #       -z ${DEBUG_MODE} \
    #       -p ${SEGMENTATION_PRIOR} \
    #       -o ${OUTPUT_PREFIX} \
    #       ${SUBJECT_ANATOMICAL_IMAGES}
    # fi
    #exit 0

    # Use single timepoint as template
    echo "Only one set of anatomical images. Treating subject as its own template."
    OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE="${OUTPUT_PREFIX}SingleSubjectTemplate/"
    SINGLE_SUBJECT_TEMPLATE=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_template0.nii.gz
    logCmd mkdir -p ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}
    logCmd DenoiseImage -d ${DIMENSION} -i ${ANATOMICAL_IMAGES[0]} -o ${SINGLE_SUBJECT_TEMPLATE}  --verbose 1
  fi

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
# Single-subject template creation
#
################################################################################

echo
echo "--------------------------------------------------------------------------------------"
echo " Creating single-subject template                                                     "
echo "--------------------------------------------------------------------------------------"
echo

TEMPLATE_MODALITY_WEIGHT_VECTOR='1'
for(( i=1; i < ${NUMBER_OF_MODALITIES}; i++ ))
  do
    TEMPLATE_MODALITY_WEIGHT_VECTOR="${TEMPLATE_MODALITY_WEIGHT_VECTOR}x1"
  done

TEMPLATE_Z_IMAGES=''

OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE="${OUTPUT_PREFIX}SingleSubjectTemplate/"

logCmd mkdir -p ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}

# Pad initial template image to avoid problems with SST drifting out of FOV
for(( i=0; i < ${NUMBER_OF_MODALITIES}; i++ ))
  do
    TEMPLATE_INPUT_IMAGE="${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}initTemplateModality${i}.nii.gz"

    logCmd ImageMath 3 ${TEMPLATE_INPUT_IMAGE} PadImage ${ANATOMICAL_IMAGES[$i]} 5

    TEMPLATE_Z_IMAGES="${TEMPLATE_Z_IMAGES} -z ${TEMPLATE_INPUT_IMAGE}"
  done


SINGLE_SUBJECT_TEMPLATE=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_template0.nii.gz

time_start_sst_creation=`date +%s`

if [[ ! -f $SINGLE_SUBJECT_TEMPLATE ]];
  then

    if [[ $RUN_OLD_ANTS_SST_CREATION -gt 0 ]];
      then
        logCmd antsMultivariateTemplateConstruction.sh \
          -d ${DIMENSION} \
          -o ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_ \
          -b 0 \
          -g 0.25 \
          -i 4 \
          -c ${DOQSUB} \
          -j ${CORES} \
          -k ${NUMBER_OF_MODALITIES} \
          -w ${TEMPLATE_MODALITY_WEIGHT_VECTOR} \
          -m 100x70x30x3  \
          -n 1 \
          -r 1 \
          -s CC \
          -t GR \
          -y ${AFFINE_UPDATE_FULL} \
          ${TEMPLATE_Z_IMAGES} \
          ${ANATOMICAL_IMAGES[@]}
    else
       logCmd antsMultivariateTemplateConstruction2.sh \
         -d ${DIMENSION} \
         -o ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_ \
         -a 0 \
         -b 0 \
         -g 0.25 \
         -i 4 \
         -c ${DOQSUB} \
         -j ${CORES} \
         -e ${USE_FLOAT_PRECISION} \
         -k ${NUMBER_OF_MODALITIES} \
         -w ${TEMPLATE_MODALITY_WEIGHT_VECTOR} \
         -q 100x70x30x3  \
         -f 8x4x2x1 \
         -s 3x2x1x0 \
         -n 1 \
         -r 1 \
         -l 1 \
         -m CC[ 4 ] \
         -t SyN \
         -y ${AFFINE_UPDATE_FULL} \
         ${TEMPLATE_Z_IMAGES} \
         ${ANATOMICAL_IMAGES[@]}
    fi

  fi

if [[ ! -f ${SINGLE_SUBJECT_TEMPLATE} ]];
  then
    echo "Error:  The single subject template was not created.  Exiting."
    exit 1
  fi

# clean up

SINGLE_SUBJECT_ANTSCT_PREFIX=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/T_template

logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}job*.sh
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}job*.txt
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}rigid*
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}*Repaired*
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}*WarpedToTemplate*
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_template0warp.nii.gz
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_template0Affine.txt
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_templatewarplog.txt
logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}initTemplateModality*.nii.gz

# Also remove the warp files but we have to be careful to not remove the affine and
# warp files generated in subsequent steps (specifically from running the SST through
# the cortical thickness pipeline if somebody has to re-run the longitudinal pipeline

if [[ -f ${SINGLE_SUBJECT_ANTSCT_PREFIX}SubjectToTemplate1Warp.nii.gz ]];
  then

    logCmd mkdir -p ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/TmpFiles/
    logCmd mv -f ${SINGLE_SUBJECT_ANTSCT_PREFIX}SubjectToTemplate1*Warp.nii.gz \
                 ${SINGLE_SUBJECT_ANTSCT_PREFIX}SubjectToTemplate0GenericAffine.mat \
                 ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/TmpFiles/
    logCmd mv -f ${SINGLE_SUBJECT_ANTSCT_PREFIX}TemplateToSubject0*Warp.nii.gz \
                 ${SINGLE_SUBJECT_ANTSCT_PREFIX}TemplateToSubject1GenericAffine.mat \
                 ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/TmpFiles/

    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*Warp.nii.gz
    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*Affine.txt
    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*GenericAffine*

    logCmd mv -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/TmpFiles/* ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}
    logCmd rm -rf ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/TmpFiles/

  else

    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*Warp.nii.gz
    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*Affine.txt
    logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}T_*GenericAffine*

  fi

# Need to change the number of iterations to  -q \

time_end_sst_creation=`date +%s`
time_elapsed_sst_creation=$((time_end_sst_creation - time_start_sst_creation))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with single subject template:  $(( time_elapsed_sst_creation / 3600 ))h $(( time_elapsed_sst_creation %3600 / 60 ))m $(( time_elapsed_sst_creation % 60 ))s"
echo "--------------------------------------------------------------------------------------"
echo

################################################################################
#
#  Create template priors
#
################################################################################

SINGLE_SUBJECT_TEMPLATE_EXTRACTION_MASK=${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainExtractionMask.${OUTPUT_SUFFIX}
SINGLE_SUBJECT_TEMPLATE_EXTRACTION_REGISTRATION_MASK=${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainExtractionRegistrationMask.${OUTPUT_SUFFIX}
SINGLE_SUBJECT_TEMPLATE_PRIOR=${SINGLE_SUBJECT_ANTSCT_PREFIX}Priors${PRIOR_NUMERIC_FORMAT}.${OUTPUT_SUFFIX}
SINGLE_SUBJECT_TEMPLATE_EXTRACTION_PRIOR=${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainExtractionMaskPrior.${OUTPUT_SUFFIX}
SINGLE_SUBJECT_TEMPLATE_CORTICAL_THICKNESS=${SINGLE_SUBJECT_ANTSCT_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}
SINGLE_SUBJECT_TEMPLATE_SKULL_STRIPPED=${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainExtractionBrain.${OUTPUT_SUFFIX}

echo
echo "--------------------------------------------------------------------------------------"
echo " Creating template priors:  running SST through antsCorticalThickness                 "
echo "--------------------------------------------------------------------------------------"
echo

time_start_priors=`date +%s`

if [[ ! -f ${SINGLE_SUBJECT_TEMPLATE_CORTICAL_THICKNESS} ]];
  then
    if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 0 ]];
      then
        logCmd antsCorticalThickness.sh \
          -d ${DIMENSION} -x ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS} \
          -q ${RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE} \
          -a ${SINGLE_SUBJECT_TEMPLATE} \
          -e ${BRAIN_TEMPLATE} \
          -f ${EXTRACTION_REGISTRATION_MASK} \
          -g ${DENOISE} \
          -m ${EXTRACTION_PRIOR} \
          -k ${KEEP_TMP_IMAGES} \
          -z ${DEBUG_MODE} \
          -p ${SEGMENTATION_PRIOR} \
          -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_SST} \
          -o ${SINGLE_SUBJECT_ANTSCT_PREFIX}
      else
        logCmd antsCorticalThickness.sh \
          -d ${DIMENSION} -x ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS}  \
          -t ${REGISTRATION_TEMPLATE} \
          -q ${RUN_FAST_ANTSCT_TO_GROUP_TEMPLATE} \
          -a ${SINGLE_SUBJECT_TEMPLATE} \
          -e ${BRAIN_TEMPLATE} \
          -f ${EXTRACTION_REGISTRATION_MASK} \
          -g ${DENOISE} \
          -m ${EXTRACTION_PRIOR} \
          -k ${KEEP_TMP_IMAGES} \
          -z ${DEBUG_MODE} \
          -p ${SEGMENTATION_PRIOR} \
          -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_SST} \
          -o ${SINGLE_SUBJECT_ANTSCT_PREFIX}
      fi
  fi

SINGLE_SUBJECT_TEMPLATE_POSTERIORS=( ${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainSegmentationPosteriors*.${OUTPUT_SUFFIX} )

SINGLE_SUBJECT_TEMPLATE_POSTERIORS_EXIST=1
SINGLE_SUBJECT_TEMPLATE_PRIORS_EXIST=1

SINGLE_SUBJECT_TEMPLATE_PRIORS=()
for (( j = 0; j < ${#SINGLE_SUBJECT_TEMPLATE_POSTERIORS[@]}; j++ ))
  do
    POSTERIOR=${SINGLE_SUBJECT_TEMPLATE_POSTERIORS[$j]}

    if [[ ! -f ${POSTERIOR} ]];
      then
        SINGLE_SUBJECT_TEMPLATE_POSTERIORS_EXIST=0
        SINGLE_SUBJECT_TEMPLATE_PRIORS_EXIST=0
        break;
      fi

    SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]=${POSTERIOR/BrainSegmentationPosteriors/Priors}

    if [[ ! -f ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} ]];
      then
        SINGLE_SUBJECT_TEMPLATE_PRIORS_EXIST=0
      fi
  done

if [[ ${SINGLE_SUBJECT_TEMPLATE_POSTERIORS_EXIST} -eq 0 ]];
  then
    echo "Error:  Posteriors for the single-subject template do not exist."
    exit 1
  fi

logCmd ImageMath ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_SKULL_STRIPPED} m ${SINGLE_SUBJECT_TEMPLATE} ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_MASK}
logCmd SmoothImage ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_MASK} 1 ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_PRIOR} 1
logCmd ImageMath ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_REGISTRATION_MASK} MD ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_MASK} 40

if [[ ${SINGLE_SUBJECT_TEMPLATE_PRIORS_EXIST} -eq 0 ]];
  then
    if [[ ${#MALF_ATLASES[@]} -eq 0 ]];
      then

        echo
        echo "   ---> Smoothing single-subject template posteriors as priors."
        echo

        for j in ${SINGLE_SUBJECT_TEMPLATE_POSTERIORS[@]}
          do
            PRIOR=${j/BrainSegmentationPosteriors/Priors}
            logCmd SmoothImage ${DIMENSION} $j 1 $PRIOR 1
          done

      else

        echo
        echo "   ---> Cooking single-subject priors using antsJointLabelFusion."
        echo

        SINGLE_SUBJECT_TEMPLATE_MALF_LABELS_PREFIX=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_TEMPLATE}/T_template
        SINGLE_SUBJECT_TEMPLATE_MALF_LABELS=${SINGLE_SUBJECT_TEMPLATE_MALF_LABELS_PREFIX}Labels.nii.gz

        ATLAS_AND_LABELS_STRING=''
        for (( j=0; j < ${#MALF_ATLASES[@]}; j++ ))
          do
            ATLAS_AND_LABELS_STRING="${ATLAS_AND_LABELS_STRING} -g ${MALF_ATLASES[$j]} -l ${MALF_LABELS[$j]}"
          done

        if [[ ! -f ${SINGLE_SUBJECT_TEMPLATE_MALF_LABELS} ]];
          then
            logCmd antsJointLabelFusion.sh \
              -d ${DIMENSION} \
              -q ${RUN_FAST_MALF_COOKING} \
              -x ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_MASK} \
              -c ${DOQSUB} \
              -j ${CORES} \
              -t ${SINGLE_SUBJECT_TEMPLATE_SKULL_STRIPPED} \
              -o ${SINGLE_SUBJECT_TEMPLATE_MALF_LABELS_PREFIX} \
              ${ATLAS_AND_LABELS_STRING}
          fi

        SINGLE_SUBJECT_TEMPLATE_PRIORS=()
        for (( j = 0; j < ${#SINGLE_SUBJECT_TEMPLATE_POSTERIORS[@]}; j++ ))
          do
            POSTERIOR=${SINGLE_SUBJECT_TEMPLATE_POSTERIORS[$j]}

            SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]=${POSTERIOR/BrainSegmentationPosteriors/Priors}

            let PRIOR_LABEL=$j+1
            logCmd ThresholdImage ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_MALF_LABELS} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} ${PRIOR_LABEL} ${PRIOR_LABEL} 1 0
            logCmd SmoothImage ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} 1 ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} 1
          done

        TMP_CSF_POSTERIOR=${SINGLE_SUBJECT_ANTSCT_PREFIX}BrainSegmentationCsfPosteriorTmp.${OUTPUT_SUFFIX}
        logCmd SmoothImage ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_POSTERIORS[0]} 1 ${TMP_CSF_POSTERIOR} 1
        logCmd ImageMath ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[0]} max ${SINGLE_SUBJECT_TEMPLATE_PRIORS[0]} ${TMP_CSF_POSTERIOR}

        # Brian's finishing touches on "cooking"---subtract out CSF from all other priors
        for (( j = 1; j < ${#SINGLE_SUBJECT_TEMPLATE_PRIORS[@]}; j++ ))
          do
            let PRIOR_LABEL=$j+1

            logCmd ImageMath ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} - ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[0]}
            logCmd ImageMath ${DIMENSION} ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} WindowImage ${SINGLE_SUBJECT_TEMPLATE_PRIORS[$j]} 0 1 0 1
          done

        logCmd rm -f $TMP_CSF_POSTERIOR
        logCmd rm -f ${SINGLE_SUBJECT_TEMPLATE_MALF_LABELS_PREFIX}*log.txt
      fi
  fi

time_end_priors=`date +%s`
time_elapsed_priors=$((time_end_priors - time_start_priors))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with creating template priors:  $(( time_elapsed_priors / 3600 ))h $(( time_elapsed_priors %3600 / 60 ))m $(( time_elapsed_priors % 60 ))s"
echo "--------------------------------------------------------------------------------------"
echo

################################################################################
#
#  Run each individual subject through antsCorticalThickness
#
################################################################################

echo
echo "--------------------------------------------------------------------------------------"
echo " Run each individual through antsCorticalThickness                                    "
echo "--------------------------------------------------------------------------------------"
echo

time_start_antsct=`date +%s`


SUBJECT_COUNT=0
for (( i=0; i < ${#ANATOMICAL_IMAGES[@]}; i+=$NUMBER_OF_MODALITIES ))
  do

    BASENAME_ID=`basename ${ANATOMICAL_IMAGES[$i]}`
    BASENAME_ID=${BASENAME_ID/\.nii\.gz/}
    BASENAME_ID=${BASENAME_ID/\.nii/}

    OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS=${OUTPUT_DIR}/${BASENAME_ID}
    OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}_${SUBJECT_COUNT}

    echo $OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS

    if [[ ! -d $OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS ]];
      then
        echo "The output directory \"$OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS\" does not exist. Making it."
        mkdir -p $OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS
      fi

    let SUBJECT_COUNT=${SUBJECT_COUNT}+1

    ANATOMICAL_REFERENCE_IMAGE=${ANATOMICAL_IMAGES[$i]}

    SUBJECT_ANATOMICAL_IMAGES=''
    if [[ ${RIGID_ALIGNMENT_TO_SST} -ne 0 ]];
      then
        logCmd antsRegistrationSyN.sh \
          -d ${DIMENSION} \
          -o ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_ID}RigidToSST \
          -m ${ANATOMICAL_IMAGES[$i]} \
          -f ${SINGLE_SUBJECT_TEMPLATE} \
          -t r
        logCmd rm -f ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_ID}RigidToSSTInverseWarped.nii.gz

        ANATOMICAL_REFERENCE_IMAGE=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_ID}RigidToSSTWarped.nii.gz


        let k=$i+$NUMBER_OF_MODALITIES
        for (( j=$i; j < $k; j++ ))
          do
            BASENAME_LOCAL_ID=`basename ${ANATOMICAL_IMAGES[$j]}`
            BASENAME_LOCAL_ID=${BASENAME_LOCAL_ID/\.nii\.gz/}
            BASENAME_LOCAL_ID=${BASENAME_LOCAL_ID/\.nii/}

            logCmd antsApplyTransforms \
              -d ${DIMENSION} \
              -i ${ANATOMICAL_IMAGES[$j]} \
              -r ${SINGLE_SUBJECT_TEMPLATE} \
              -o ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_LOCAL_ID}RigidToSSTWarped.nii.gz \
              -n Linear \
              -t ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_ID}RigidToSST0GenericAffine.mat

            SUBJECT_ANATOMICAL_IMAGES="${SUBJECT_ANATOMICAL_IMAGES} -a ${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_LOCAL_ID}RigidToSSTWarped.nii.gz"
          done
      else
        let k=$i+$NUMBER_OF_MODALITIES
        for (( j=$i; j < $k; j++ ))
          do
            SUBJECT_ANATOMICAL_IMAGES="${SUBJECT_ANATOMICAL_IMAGES} -a ${ANATOMICAL_IMAGES[$j]}"
          done
      fi

    if [[ ${USE_SST_CORTICAL_THICKNESS_PRIOR} -ne 0 ]];
      then
        SUBJECT_ANATOMICAL_IMAGES="${SUBJECT_ANATOMICAL_IMAGES} -r ${SINGLE_SUBJECT_TEMPLATE_CORTICAL_THICKNESS}"
      fi

    OUTPUT_LOCAL_PREFIX=${OUTPUT_DIRECTORY_FOR_SINGLE_SUBJECT_CORTICAL_THICKNESS}/${BASENAME_ID}

    logCmd antsCorticalThickness.sh \
      -d ${DIMENSION} -x ${ATROPOS_SEGMENTATION_INTERNAL_ITERATIONS}  \
      -q ${RUN_ANTSCT_TO_SST_QUICK} \
      ${SUBJECT_ANATOMICAL_IMAGES} \
      -e ${SINGLE_SUBJECT_TEMPLATE} \
      -m ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_PRIOR} \
      -f ${SINGLE_SUBJECT_TEMPLATE_EXTRACTION_REGISTRATION_MASK} \
      -g ${DENOISE} \
      -k ${KEEP_TMP_IMAGES} \
      -z ${DEBUG_MODE} \
      -w ${ATROPOS_SEGMENTATION_PRIOR_WEIGHT_TIMEPOINT} \
      -p ${SINGLE_SUBJECT_TEMPLATE_PRIOR} \
      -t ${SINGLE_SUBJECT_TEMPLATE_SKULL_STRIPPED} \
      -o ${OUTPUT_LOCAL_PREFIX}

    if [[ $DO_REGISTRATION_TO_TEMPLATE -eq 1 ]];
      then

        if [[ ! -f ${OUTPUT_LOCAL_PREFIX}SubjectToGroupTemplateWarp.nii.gz ]];
          then
            logCmd antsApplyTransforms \
              -d ${DIMENSION} \
              -r ${REGISTRATION_TEMPLATE} \
              -o [ ${OUTPUT_LOCAL_PREFIX}SubjectToGroupTemplateWarp.nii.gz,1 ] \
              -t ${SINGLE_SUBJECT_ANTSCT_PREFIX}SubjectToTemplate1Warp.nii.gz \
              -t ${SINGLE_SUBJECT_ANTSCT_PREFIX}SubjectToTemplate0GenericAffine.mat \
              -t ${OUTPUT_LOCAL_PREFIX}SubjectToTemplate1Warp.nii.gz \
              -t ${OUTPUT_LOCAL_PREFIX}SubjectToTemplate0GenericAffine.mat
          fi

        if [[ ! -f ${OUTPUT_LOCAL_PREFIX}GroupTemplateToSubjectWarp.nii.gz ]];
          then
            logCmd antsApplyTransforms \
              -d ${DIMENSION} \
              -r ${ANATOMICAL_REFERENCE_IMAGE} \
              -o [ ${OUTPUT_LOCAL_PREFIX}GroupTemplateToSubjectWarp.nii.gz,1 ] \
              -t ${OUTPUT_LOCAL_PREFIX}TemplateToSubject1GenericAffine.mat \
              -t ${OUTPUT_LOCAL_PREFIX}TemplateToSubject0Warp.nii.gz \
              -t ${SINGLE_SUBJECT_ANTSCT_PREFIX}TemplateToSubject1GenericAffine.mat \
              -t ${SINGLE_SUBJECT_ANTSCT_PREFIX}TemplateToSubject0Warp.nii.gz
          fi

        if [[ -f ${CORTICAL_LABEL_IMAGE} ]];
          then

            SUBJECT_CORTICAL_LABELS=${OUTPUT_LOCAL_PREFIX}CorticalLabels.${OUTPUT_SUFFIX}
            SUBJECT_CORTICAL_THICKNESS=${OUTPUT_LOCAL_PREFIX}CorticalThickness.${OUTPUT_SUFFIX}
            SUBJECT_TMP=${OUTPUT_LOCAL_PREFIX}Tmp.${OUTPUT_SUFFIX}
            SUBJECT_STATS=${OUTPUT_LOCAL_PREFIX}LabelThickness.csv

            if [[ ! -f ${SUBJECT_CORTICAL_LABELS} ]];
              then
                logCmd antsApplyTransforms \
                  -d ${DIMENSION} \
                  -i ${CORTICAL_LABEL_IMAGE} \
                  -r ${ANATOMICAL_REFERENCE_IMAGE} \
                  -o ${SUBJECT_CORTICAL_LABELS} \
                  -n MultiLabel \
                  -t ${OUTPUT_LOCAL_PREFIX}GroupTemplateToSubjectWarp.nii.gz

                logCmd ThresholdImage ${DIMENSION} ${OUTPUT_LOCAL_PREFIX}BrainSegmentation.${OUTPUT_SUFFIX} ${SUBJECT_TMP} 2 2 1 0
                logCmd ImageMath ${DIMENSION} ${SUBJECT_CORTICAL_LABELS} m ${SUBJECT_TMP} ${SUBJECT_CORTICAL_LABELS}
                logCmd ImageMath ${DIMENSION} ${SUBJECT_STATS} LabelStats ${SUBJECT_CORTICAL_LABELS} ${SUBJECT_CORTICAL_THICKNESS}
              fi

            logCmd rm -f $SUBJECT_TMP
          fi
      fi

  done

time_end_antsct=`date +%s`
time_elapsed_antsct=$((time_end_antsct - time_start_antsct))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with individual cortical thickness:  $(( time_elapsed_antsct / 3600 ))h $(( time_elapsed_antsct %3600 / 60 ))m $(( time_elapsed_antsct % 60 ))s"
echo "--------------------------------------------------------------------------------------"
echo

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with ANTs longitudinal processing pipeline"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
