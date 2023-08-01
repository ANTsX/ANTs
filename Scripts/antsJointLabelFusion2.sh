#!/bin/bash

VERSION="0.0.0"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
JLF=antsJointFusion
ANTS=antsRegistration
WARP=antsApplyTransforms
SMOOTH=SmoothDisplacementField
PEXEC=ANTSpexec.sh
SGE=waitForSGEQJobs.pl
PBS=waitForPBSQJobs.pl
XGRID=waitForXGridJobs.pl
SLURM=waitForSlurmJobs.pl

fle_error=0
for FLE in $JLF $ANTS $WARP $SMOOTH $PEXEC $SGE $XGRID $PBS $SLURM
  do
  if ! command -v $FLE &> /dev/null
    then
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " FILE $FLE DOES NOT EXIST -- OR -- IS NOT EXECUTABLE !!! $0 will terminate."
      echo "--------------------------------------------------------------------------------------"
      echo " if the file is not executable, please change its permissions. "
      fle_error=1
    fi
  done

if [[ $fle_error = 1 ]];
  then
    echo "missing helper script"
    exit 1
  fi


# Assuming .nii.gz as default file type. This is the case for ANTS 1.7 and up

function Usage {
    cat <<USAGE

This adds some additional utilty to the original antsJointLabelFusion.sh script.
Specifically, it adds the parameters '-a', '-b', '-e', and '-s' for performing
multiple piecewise registration steps.  For example, suppose one wants to use
joint label fusion to label thoracic CT where the following atlases have been
supplied with the following labels:

  label 1:  left lung
  label 2:  right lung
  label 3:  spinal cord
  label 4:  heart

This script allows one to perform a piecewise registration where one can register
the two lungs separately and then "stitch" together the resulting transforms.  This
is followed up by a registration refinement step where other labeled structures
can be incorporated into driving the registration.  An additional benefit is that
one can use the segmentation labels of the atlases to create fixed masks for the
registration.  Example usage for the above would be

antsJointLabelFusion2.sh -d 3 ... \
                         -b 1 -b 2 -b 3 -b 4 \
                         -e 1 -e 2 \
                         -a 4 \
                         -s 10 \
                         ...

Usage:

`basename $0` -d ImageDimension -o OutputPrefix <other options> <images>

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3.

     -o:  OutputPrefix:   A prefix that is prepended to all output files.

     -t:  TargetImage:    Target image to be labeled.

     -g:  Atlas:          Atlas to be warped to target image.

     -l:  Labels:         Labels corresponding to atlas (cf -g).

Optional arguments:

     -m:  Majority vote:  Use majority vote instead of joint label fusion (default = 0).

     -k:  Keep files:     Keep warped atlas and label files (default = 0).

     -c:  Control for parallel computation (default 0) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM.

     -j:  Number of cpu cores to use (default 2; -- requires "-c 2").

     -r:  qsub options

     -q:  Use quick registration parameters:  Either 0 or 1 (default = 1).

     -a:  isotropic resampling in mm for specifying labels (default = 'NA' = "no resampling").

     -b:  specify labels to use for dilation and masking for the final registration

     -e:  subset of labels specified by the '-b' option to be used for the piecewise registration.
          If no labels are specified by '-e', it is assumed that all the '-b' labels are to be
          used for the initial piecewise registration.

     -s:  dilation radii for labels specified with option 'b'

     -p:  Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f:  Float precision: Use float precision (default = 1) -- 0 == double, 1 == float.

     -u:  Registration walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Registration memory limit (default = 8gb):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -w:  JLF walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          for the joint label fusion call.

     -z:  JLF Memory limit (default = 8gb):  Option for PBS/SLURM qsub specifying requested memory
          for the joint label fusion call.

     -y:  transform type (default = 's')
        t: translation
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        sr: rigid + deformable syn
        so: deformable syn only
        b: rigid + affine + deformable b-spline syn
        br: rigid + deformable b-spline syn
        bo: deformable b-spline syn only

     -x:  Target mask image (default = 'majorityvoting')
        otsu: use otsu thresholding to define foreground/background
        or: 'or' all the warped atlas images to defined foreground/background
        majorityvoting: perform a voxelwise label voting.  If >= 80% of the warped atlases agree at that
                        voxel, we keep that voted label at that voxel and *do not* perform JLF.  Note that
                        the 80% threshold is hard-coded but can be easily changed in the script.
        <filename>: a user-specified mask
        none: don't use a mask

Example:

`basename $0` -d 3 -t target.nii.gz -o malf \
              -p malfPosteriors%04d.nii.gz \
              -g atlas1.nii.gz -l labels1.nii.gz \
              -g atlas2.nii.gz -l labels2.nii.gz \
              -g atlas3.nii.gz -l labels3.nii.gz

--------------------------------------------------------------------------------------
JLF was created by:
--------------------------------------------------------------------------------------
Hongzhi Wang and Paul Yushkevich
Penn Image Computing And Science Laboratory
University of Pennsylvania

Please reference http://www.ncbi.nlm.nih.gov/pubmed/22732662 when employing this script
in your studies.

Wang H, Suh JW, Das SR, Pluta J, Craige C, Yushkevich PA.
Multi-Atlas Segmentation with Joint Label Fusion.
IEEE Trans Pattern Anal Mach Intell.

For lung-specific applications, cf https://www.ncbi.nlm.nih.gov/pubmed/26222827

Tustison NJ, Qing K, Wang C, Altes TA, Mugler JP 3rd.
Atlas-based estimation of lung and lobar anatomy in proton MRI.
Magn Reson Med. 2016 Jul;76(1):315-20.

--------------------------------------------------------------------------------------
script by Nick Tustison
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

function Help {
    cat <<HELP

This adds some additional utilty to the original antsJointLabelFusion.sh script.
Specifically, it adds the parameters '-a', '-b', '-e', and '-s' for performing
multiple piecewise registration steps.  For example, suppose one wants to use
joint label fusion to label thoracic CT where the following atlases have been
supplied with the following labels:

  label 1:  left lung
  label 2:  right lung
  label 3:  spinal cord
  label 4:  heart

This script allows one to perform a piecewise registration where one can register
the two lungs separately and then "stitch" together the resulting transforms.  This
is followed up by a registration refinement step where other labeled structures
can be incorporated into driving the registration.  An additional benefit is that
one can use the segmentation labels of the atlases to create fixed masks for the
registration.  Example usage for the above would be

antsJointLabelFusion2.sh -d 3 ... \
                         -b 1 -b 2 -b 3 -b 4 \
                         -e 1 -e 2 \
                         -a 4 \
                         -s 10 \
                         ...

`basename $0` will propagate labels from a set of pre-labeled atlases using the JLF
algorithm.

Usage:

`basename $0` -d ImageDimension -o OutputPrefix <other options> <images>

Example Case:

`basename $0` -d 3 -t target.nii.gz -o malf \
              -p malfPosteriors%04d.nii.gz \
              -g atlas1.nii.gz -l labels1.nii.gz \
              -g atlas2.nii.gz -l labels2.nii.gz \
              -g atlas3.nii.gz -l labels3.nii.gz


Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3.

     -o:  OutputPrefix:   A prefix that is prepended to all output files.

     -t:  TargetImage:    Target image to be labeled.

     -g:  Atlas:          Atlas to be warped to target image.

     -l:  Labels:         Labels corresponding to atlas (cf -g).

Optional arguments:

     -m:  Majority vote:  Use majority vote instead of joint label fusion (default = 0).

     -k:  Keep files:     Keep warped atlas and label files (default = 0).

     -c:  Control for parallel computation (default 0) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM.

     -j:  Number of cpu cores to use (default 2; -- requires "-c 2").

     -q:  Use quick registration parameters:  Either 0 or 1 (default = 1).

     -b:  specify labels to use for dilation and masking for the final registration

     -e:  subset of labels specified by the '-b' option to be used for the piecewise registration.
          If no labels are specified by '-e', it is assumed that all the '-b' labels are to be
          used for the initial piecewise registration.

     -s:  dilation radii (in mm) for labels specified with option 'b'

     -a:  isotropic resampling in mm for specifying labels (default = 4).

     -p:  Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f:  Float precision: Use float precision (default = 1) -- 0 == double, 1 == float.

     -u:  Registration walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Registration memory limit (default = 8gb):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -w:  JLF walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          for the joint label fusion call.

     -z:  JLF Memory limit (default = 8gb):  Option for PBS/SLURM qsub specifying requested memory
          for the joint label fusion call.

     -y:  Transform type (default = 's')
        t: translation
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        sr: rigid + deformable syn
        so: deformable syn only
        b: rigid + affine + deformable b-spline syn
        br: rigid + deformable b-spline syn
        bo: deformable b-spline syn only

     -x:  Target mask image (default = 'majorityvoting')
        otsu: use otsu thresholding to define foreground/background
        or: 'or' all the warped atlas images to defined foreground/background
        majorityvoting: perform a voxelwise label voting.  If >= 80% of the warped atlases agree at that
                        voxel, we keep that voted label at that voxel and *do not* perform JLF.  Note that
                        the 80% threshold is hard-coded but can be easily changed in the script.
        <filename>: a user-specified mask
        none: don't use a mask

Requirements:

This scripts relies on the following scripts in your $PATH directory. The script
will terminate prematurely if these files are not present or are not executable.
- pexec.sh
- waitForSGEQJobs.pl (only for use with Sun Grid Engine)
- waitForPBSQJobs.pl  (only for use with Portable Batch System)
- ANTSpexec.sh (only for use with localhost parallel execution)
- waitForXGridJobs.pl (only for use with Apple XGrid)
- waitForSlurmJobs.pl (only for use with SLURM)
- antsRegistrationSyN.sh
- antsRegistrationSyNQuick.sh ( quick parameters )

--------------------------------------------------------------------------------------
Get the latest ANTS version at:
--------------------------------------------------------------------------------------
https://github.com/stnava/ANTs/
--------------------------------------------------------------------------------------
Read the ANTS documentation at:
--------------------------------------------------------------------------------------
http://stnava.github.io/ANTs/

--------------------------------------------------------------------------------------
JLF was created by:
--------------------------------------------------------------------------------------
Hongzhi Wang and Paul Yushkevich
Penn Image Computing And Science Laboratory
University of Pennsylvania

Please reference http://www.ncbi.nlm.nih.gov/pubmed/22732662 when employing this script
in your studies.

Wang H, Suh JW, Das SR, Pluta J, Craige C, Yushkevich PA.
Multi-Atlas Segmentation with Joint Label Fusion.
IEEE Trans Pattern Anal Mach Intell.
--------------------------------------------------------------------------------------
script by Nick Tustison
--------------------------------------------------------------------------------------

HELP
    exit 1
}

function reportParameters {
    cat <<REPORTPARAMETERS

--------------------------------------------------------------------------------------
 Parameters
--------------------------------------------------------------------------------------
 Dimensionality:           $DIM
 Output prefix:            $OUTPUT_PREFIX
 Posteriors format:        $OUTPUT_POSTERIORS_FORMAT
 Target image:             $TARGET_IMAGE
 Atlas images:             ${ATLAS_IMAGES[@]}
 Atlas labels:             ${ATLAS_LABELS[@]}
 Transformation:           ${TRANSFORM_TYPE}

 Keep all images:          $KEEP_ALL_IMAGES
 Processing type:          $DOQSUB
 Number of cpu cores:      $CORES
--------------------------------------------------------------------------------------
REPORTPARAMETERS
}

cleanup()
{
  echo "\n*** Performing cleanup, please wait ***\n"

  runningANTSpids=$( ps --ppid $$ -o pid= )

  for thePID in $runningANTSpids
  do
      echo "killing:  ${thePID}"
      kill ${thePID}
  done

  return $?
}

control_c()
# run if user hits control-c
{
  echo -en "\n*** User pressed CTRL + C ***\n"
  cleanup
  exit $?
  echo -en "\n*** Script cancelled by user ***\n"
}

#initializing variables with global scope
time_start=`date +%s`
CURRENT_DIR=`pwd`/

DIM=3

OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"
OUTPUT_POSTERIORS_FORMAT=''

TARGET_IMAGE=''
ATLAS_IMAGES=()
ATLAS_LABELS=()
TRANSFORM='s'

PIECEWISE_LABELS=()
PIECEWISE_LABELS_SUBSET=()
PIECEWISE_DILATION_RADII=()
ISOTROPIC_RESAMPLING='NA'

KEEP_ALL_IMAGES=0
DOQSUB=0
CORES=1
PRECISION=0

XGRID_OPTS=""
SCRIPT_PREPEND=""
QSUB_OPTS=""
TARGET_MASK_IMAGE="majorityvoting"

REGISTRATION_WALLTIME="20:00:00"
REGISTRATION_MEMORY="8gb"
JLF_WALLTIME="20:00:00"
JLF_MEMORY="8gb"

##Getting system info from linux can be done with these variables.
# RAM=`cat /proc/meminfo | sed -n -e '/MemTotal/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# RAMfree=`cat /proc/meminfo | sed -n -e '/MemFree/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# cpu_free_ram=$((${RAMfree}/${cpu_count}))

if [[ ${OSTYPE:0:6} == 'darwin' ]];
  then
    cpu_count=`sysctl -n hw.physicalcpu`
  else
    cpu_count=`cat /proc/cpuinfo | grep processor | wc -l`
  fi

# Provide output for Help
if [[ "$1" == "-h" ]];
  then
    Help >&2
  fi
MAJORITYVOTE=0
RUNQUICK=1
TRANSFORM_TYPE="s"
# reading command line arguments
while getopts "a:b:c:d:e:f:g:h:j:k:l:m:o:p:q:r:s:t:u:v:w:x:y:z:" OPT
  do
  case $OPT in
      h) #help
   echo "$USAGE"
   exit 0
   ;;
      c) #use SGE cluster
   DOQSUB=$OPTARG
   if [[ $DOQSUB -gt 5 ]];
     then
       echo " DOQSUB must be an integer value (0=serial, 1=SGE qsub, 2=try pexec, 3=XGrid, 4=PBS qsub, 5=SLURM ) you passed  -c $DOQSUB "
       exit 1
     fi
   ;;
      d) #dimensions
   DIM=$OPTARG
   if [[ ${DIM} -ne 2 && $DIM -ne 3 ]];
     then
       echo " Dimensionality is only valid for 2 or 3.  You passed -d $DIM."
       exit 1
     fi
   ;;
      f)
      PRECISION=$OPTARG
   ;;
      g)
   ATLAS_IMAGES[${#ATLAS_IMAGES[@]}]=$OPTARG
   ;;
      j) #number of cpu cores to use
   CORES=$OPTARG
   ;;
      k)
   KEEP_ALL_IMAGES=$OPTARG
   ;;
      m) #majority voting option
   MAJORITYVOTE=$OPTARG
   ;;
      p)
   OUTPUT_POSTERIORS_FORMAT=$OPTARG
   ;;
      q)
   RUNQUICK=$OPTARG
   ;;
      a)
   ISOTROPIC_RESAMPLING=$OPTARG
   ;;
      b)
   PIECEWISE_LABELS[${#PIECEWISE_LABELS[@]}]=$OPTARG
   ;;
      e)
   PIECEWISE_LABELS_SUBSET[${#PIECEWISE_LABELS_SUBSET[@]}]=$OPTARG
   ;;
      s)
   PIECEWISE_DILATION_RADII[${#PIECEWISE_DILATION_RADII[@]}]=$OPTARG
   ;;
      j) #number of cpu cores to use
   CORES=$OPTARG
   ;;
      r)
   QSUB_OPTS=$OPTARG
   ;;
      l)
   ATLAS_LABELS[${#ATLAS_LABELS[@]}]=$OPTARG
   ;;
      o)
   OUTPUT_PREFIX=$OPTARG
   OUTPUT_DIR=`dirname ${OUTPUT_PREFIX}`
   ;;
      t)
   TARGET_IMAGE=$OPTARG
   ;;
      u)
   REGISTRATION_WALLTIME=$OPTARG
   ;;
      v)
   REGISTRATION_MEMORY=$OPTARG
   ;;
      w)
   JLF_WALLTIME=$OPTARG
   ;;
      z)
   JLF_MEMORY=$OPTARG
   ;;
      x)
   TARGET_MASK_IMAGE=$OPTARG
   ;;
      y)
   TRANSFORM_TYPE=$OPTARG
   ;;
      \?) # getopts issues an error message
      echo "$USAGE" >&2
      exit 1
      ;;
  esac
done

if [[ $DOQSUB -eq 1 || $DOQSUB -eq 4 ]];
  then
    qq=`which  qsub`
    if [[  ${#qq} -lt 1 ]];
      then
        echo "do you have qsub?  if not, then choose another c option ... if so, then check where the qsub alias points ..."
        exit 1
      fi
  fi
if [[ $DOQSUB -eq 5 ]];
  then
    qq=`which sbatch`
    if [[ ${#qq} -lt 1 ]];
      then
        echo "do you have sbatch?  if not, then choose another c option ... if so, then check where the sbatch alias points ..."
        exit
      fi
  fi

if [[ ! -f "$TARGET_IMAGE" ]];
  then
    echo "Target image '$TARGET_IMAGE' does not exist.  See usage: '$0 -h 1'"
    exit
  fi

if [[ ${#ATLAS_IMAGES[@]} -ne ${#ATLAS_LABELS[@]} ]];
  then
    echo "The number of atlas images does not equal the number of labels.  Ensure that a corresponding set of labels exist for each image."
    exit 1
  fi

if [[ ${#PIECEWISE_LABELS[@]} -gt 0 ]];
  then
    if [[ ${#PIECEWISE_DILATION_RADII[@]} -eq 1 ]];
      then
        for (( i = 0; i < ${#PIECEWISE_LABELS[@]}; i++ ))
          do
            PIECEWISE_DILATION_RADII[${i}]=${PIECEWISE_DILATION_RADII[0]}
          done
      elif [[ ${#PIECEWISE_LABELS[@]} -ne ${#PIECEWISE_DILATION_RADII[@]} ]];
      then
        echo "The number of dilation radii does not equal the number of piecewise labels."
        exit 1
      fi

    if [[ ${#PIECEWISE_LABELS_SUBSET[@]} -eq 0 ]];
      then
        for (( i = 0; i < ${#PIECEWISE_LABELS[@]}; i++ ))
          do
            PIECEWISE_LABELS_SUBSET[$i]=${PIECEWISE_LABELS[$i]}
          done
      else
        for (( i = 0; i < ${#PIECEWISE_LABELS_SUBSET[@]}; i++ ))
          do
            LABEL_IS_FOUND=0
            for (( j = 0; j < ${#PIECEWISE_LABELS[@]}; j++ ))
              do
                if [[ ${PIECEWISE_LABELS[$j]} -eq ${PIECEWISE_LABELS_SUBSET[$i]} ]];
                  then
                    LABEL_IS_FOUND=1
                    break
                  fi
              done
            if [[ $LABEL_IS_FOUND -eq 0 ]];
              then
                echo "Subset label ${PIECEWISE_LABELS_SUBSET[$i]} is not found in the labels specified by the -b option."
                exit 1
              fi
          done
      fi

  fi



ISOTROPIC_RESAMPLING_VECTOR="${ISOTROPIC_RESAMPLING}x${ISOTROPIC_RESAMPLING}"
if [[ $DIM -eq 3 ]];
  then
    ISOTROPIC_RESAMPLING_VECTOR="${ISOTROPIC_RESAMPLING}x${ISOTROPIC_RESAMPLING}x${ISOTROPIC_RESAMPLING}"
  fi

PRECISIONFLAG='f'
if [[ ${PRECISION} -eq 0 ]];
  then
    PRECISIONFLAG='d'
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


##########################################################################
#
# Perform JLF labeling by
#  1) registering all atlases to target image
#  2) call 'jointfusion'
#
##########################################################################

echo
echo "--------------------------------------------------------------------------------------"
echo " Start JLFization"
echo "--------------------------------------------------------------------------------------"
reportParameters

jobIDs=""

WARPED_ATLAS_IMAGES=()
WARPED_ATLAS_LABELS=()
AFFINE_FILES=()
WARP_FIELDS=()
INVERSE_WARP_FIELDS=()

for (( i = 0; i < ${#ATLAS_IMAGES[@]}; i++ ))
  do
    IMG_BASE=`basename ${ATLAS_IMAGES[$i]}`
    BASENAME=` echo ${IMG_BASE} | cut -d '.' -f 1 `

    qscript="${OUTPUT_DIR}/job_${BASENAME}_${i}.sh"

    WARPED_ATLAS_IMAGES[${#WARPED_ATLAS_IMAGES[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz"
    WARPED_ATLAS_LABELS[${#WARPED_ATLAS_LABELS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz"

    if [[ -f ${WARPED_ATLAS_LABELS[${i}]} ]];
      then
        echo "${WARPED_ATLAS_LABELS[${i}]} already exists."
        rm -f $qscript
        continue
      fi

    regcall=antsRegistrationSyN.sh
    if [[ $RUNQUICK -eq 1 ]];
      then
        regcall=antsRegistrationSyNQuick.sh
      fi

    if [[ ${#PIECEWISE_LABELS[@]} -eq 0 ]];
      then

        INVERSE_WARPED_ATLAS_IMAGES[${#INVERSE_WARPED_ATLAS_IMAGES[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_InverseWarped.nii.gz"
        WARP_FIELDS[${#WARP_FIELDS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_1Warp.nii.gz"
        INVERSE_WARP_FIELDS[${#INVERSE_WARP_FIELDS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_1InverseWarp.nii.gz"
        AFFINE_FILES[${#AFFINE_FILES[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_0GenericAffine.mat"

        registrationCall="$regcall \
                              -d ${DIM} \
                              -p ${PRECISIONFLAG} \
                              -j 1 \
                              -t ${TRANSFORM_TYPE} \
                              -f ${TARGET_IMAGE} \
                              -m ${ATLAS_IMAGES[$i]} \
                              -o ${OUTPUT_PREFIX}${BASENAME}_${i}_ > ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"

        labelXfrmCall="antsApplyTransforms \
                              -d ${DIM} \
                              --float 1 \
                              -i ${ATLAS_LABELS[$i]} \
                              -r ${TARGET_IMAGE} \
                              -o ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz \
                              -n GenericLabel \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_1Warp.nii.gz \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_0GenericAffine.mat >> ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"

        copyImageHeaderCall="CopyImageHeaderInformation \
                             ${TARGET_IMAGE} \
                             ${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz \
                             ${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz 1 1 1"

        copyLabelsHeaderCall="CopyImageHeaderInformation \
                              ${TARGET_IMAGE} \
                              ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz \
                              ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz 1 1 1"

        rm -f $qscript

        if [[ $DOQSUB -eq 5 ]];
          then
            # SLURM job scripts must start with a shebang
            echo '#!/bin/sh' > $qscript
          fi

        echo "$registrationCall" >> $qscript
        echo "$labelXfrmCall" >> $qscript
        echo "$copyImageHeaderCall" >> $qscript
        echo "$copyLabelsHeaderCall" >> $qscript

      else

        TMP_FILES=()

        logFile="${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"

        # Dilate regions and create masks

        SUBSET_CONFIDENCE_MASK="${OUTPUT_PREFIX}${BASENAME}_${i}_ConfidenceMask.nii.gz"
        CONFIDENCE_DILATED_MASK="${OUTPUT_PREFIX}${BASENAME}_${i}_ConfidenceDilatedMask.nii.gz"
        PIECEWISE_MASKS=()
        PIECEWISE_DILATED_MASKS=()
        PIECEWISE_OUTPUT_PREFIX="${OUTPUT_PREFIX}${BASENAME}_${i}_PiecewiseLabel"

        INITIAL_LABEL_IS_SPECIFIED=0

        for (( j = 0; j < ${#PIECEWISE_LABELS[@]}; j++ ))
          do

            PIECEWISE_MASKS[$j]="${PIECEWISE_OUTPUT_PREFIX}${PIECEWISE_LABELS[$j]}.nii.gz"
            PIECEWISE_DILATED_MASKS[$j]="${PIECEWISE_OUTPUT_PREFIX}Dilated${PIECEWISE_LABELS[$j]}.nii.gz"
            thresholdCall="ThresholdImage ${DIM} ${ATLAS_LABELS[$i]} ${PIECEWISE_MASKS[$j]} ${PIECEWISE_LABELS[$j]} ${PIECEWISE_LABELS[$j]} 1 0"

            echo "$thresholdCall" >> $qscript

            distanceCall="ImageMath ${DIM} ${PIECEWISE_DILATED_MASKS[$j]} MaurerDistance ${PIECEWISE_MASKS[$j]} 1"
            thresholdCall="ThresholdImage ${DIM} ${PIECEWISE_DILATED_MASKS[$j]} ${PIECEWISE_DILATED_MASKS[$j]} -1000000 ${PIECEWISE_DILATION_RADII[$j]} 1 0"

            echo "$distanceCall" >> $qscript
            echo "$thresholdCall" >> $qscript

            if [[ $j -eq 0 ]]
              then
                echo "cp ${PIECEWISE_DILATED_MASKS[$j]} ${CONFIDENCE_DILATED_MASK}" >> $qscript
              else
                addCall="ImageMath ${DIM} ${CONFIDENCE_DILATED_MASK} max ${CONFIDENCE_DILATED_MASK} ${PIECEWISE_DILATED_MASKS[$j]}"
                echo "$addCall" >> $qscript
              fi

            SUBSET_LABEL_IS_FOUND=0
            for (( k = 0; k < ${#PIECEWISE_LABELS_SUBSET[@]}; k++ ))
              do
                if [[ ${PIECEWISE_LABELS[$j]} -eq ${PIECEWISE_LABELS_SUBSET[$k]} ]];
                  then
                    SUBSET_LABEL_IS_FOUND=1
                    break
                  fi
              done

            if [[ $SUBSET_LABEL_IS_FOUND -eq 0 ]];
              then
                continue
              fi

            if [[ $INITIAL_LABEL_IS_SPECIFIED -eq 0 ]]
              then
                echo "cp ${PIECEWISE_MASKS[$j]} ${SUBSET_CONFIDENCE_MASK}" >> $qscript
                INITIAL_LABEL_IS_SPECIFIED=1
              else
                addCall="ImageMath ${DIM} ${SUBSET_CONFIDENCE_MASK} max ${SUBSET_CONFIDENCE_MASK} ${PIECEWISE_MASKS[$j]}"
                echo "$addCall" >> $qscript
              fi
          done

        if [[ ${ISOTROPIC_RESAMPLING} != 'NA' ]];
          then
            resampleCall="ResampleImage ${DIM} ${SUBSET_CONFIDENCE_MASK} ${SUBSET_CONFIDENCE_MASK} $ISOTROPIC_RESAMPLING_VECTOR 0 1"
            echo "$resampleCall" >> $qscript
          fi

        TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_MASKS[@]} ${PIECEWISE_DILATED_MASKS[@]} ${SUBSET_CONFIDENCE_MASK} ${CONFIDENCE_DILATED_MASK} )

        # Perform initial linear registration of the whole images
        # Note that we're registering the target image to the atlas image since the atlas
        # image has the masks (i.e., labels)

        INIT_OUTPUT_PREFIX=${OUTPUT_PREFIX}${BASENAME}_${i}_Init

        registrationCallInit="antsRegistration -d 3 -v 1 -o $INIT_OUTPUT_PREFIX -r [ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},0 ]"
        registrationCallInit="$registrationCallInit -t Rigid[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x40x10,1e-6,10 ] -f 6x4x3 -s 3x2x1"
        registrationCallInit="$registrationCallInit -t Similarity[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x40x10,1e-6,10 ] -f 6x4x3 -s 3x2x1"
        registrationCallInit="$registrationCallInit -t Affine[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x40x10,1e-6,10 ] -f 6x4x3 -s 3x2x1  >> $logFile"

        echo "$registrationCallInit" >> $qscript

        # Register piecewise labels separately

        PIECEWISE_TRANSFORM_WARPS=""

        for (( j = 0; j < ${#PIECEWISE_LABELS[@]}; j++ ))
          do

            SUBSET_LABEL_IS_FOUND=0
            for (( k = 0; k < ${#PIECEWISE_LABELS_SUBSET[@]}; k++ ))
              do
                if [[ ${PIECEWISE_LABELS[$j]} -eq ${PIECEWISE_LABELS_SUBSET[$k]} ]];
                  then
                    SUBSET_LABEL_IS_FOUND=1
                    break
                  fi
              done

            if [[ $SUBSET_LABEL_IS_FOUND -eq 0 ]];
              then
                continue
              fi

            PIECEWISE_OUTPUT_PREFIX=${OUTPUT_PREFIX}${BASENAME}_${i}_PiecewiseLabel${PIECEWISE_LABELS[$j]}_

            registrationCallPiecewise="antsRegistration -d 3 -v 1 -z 0 -o $PIECEWISE_OUTPUT_PREFIX -r ${INIT_OUTPUT_PREFIX}0GenericAffine.mat -x ${PIECEWISE_DILATED_MASKS[$j]}"
            registrationCallPiecewise="$registrationCallPiecewise -t Rigid[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x100x50,1e-6,10 ] -f 6x4x3 -s 3x2x1"
            registrationCallPiecewise="$registrationCallPiecewise -t Similarity[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x100x50,1e-6,10 ] -f 6x4x3 -s 3x2x1"
            registrationCallPiecewise="$registrationCallPiecewise -t Affine[ 0.15 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x100x50,1e-6,10 ] -f 6x4x3 -s 3x2x1"
            registrationCallPiecewise="$registrationCallPiecewise -t BSplineSyN[ 0.15,40,0 ] -m MI[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,32,Regular,0.25 ] -c [ 100x40x0,1e-6,10 ] -f 6x4x1 -s 3x2x0  >> $logFile"

            echo "$registrationCallPiecewise" >> $qscript

            composeCallPiecewise="$WARP \
                                 -d ${DIM} \
                                 -v 1 \
                                 -o [ ${PIECEWISE_OUTPUT_PREFIX}TotalWarp.nii.gz,1 ] \
                                 -r ${ATLAS_IMAGES[$i]} \
                                 -t ${PIECEWISE_OUTPUT_PREFIX}4Warp.nii.gz \
                                 -t ${PIECEWISE_OUTPUT_PREFIX}3Affine.mat \
                                 -t ${PIECEWISE_OUTPUT_PREFIX}2Similarity.mat \
                                 -t ${PIECEWISE_OUTPUT_PREFIX}1Rigid.mat >> $logFile"
            echo "$composeCallPiecewise" >> $qscript

            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}1Rigid.mat )
            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}2Similarity.mat )
            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}3Affine.mat )
            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}4Warp.nii.gz )
            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}4InverseWarp.nii.gz )
            TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}TotalWarp.nii.gz )

            componentCall="ConvertImage ${DIM} ${PIECEWISE_OUTPUT_PREFIX}TotalWarp.nii.gz ${PIECEWISE_OUTPUT_PREFIX}TotalWarp 10"
            echo "$componentCall" >> $qscript

            components=( 'xvec' 'yvec' 'zvec' )
            for (( d = 0; d < ${DIM}; d++ ))
              do
                maskingCall="ImageMath 3 ${PIECEWISE_OUTPUT_PREFIX}TotalWarp${components[$d]}.nii.gz m ${PIECEWISE_OUTPUT_PREFIX}TotalWarp${components[$d]}.nii.gz ${PIECEWISE_MASKS[$j]}"
                echo "$maskingCall" >> $qscript

                if [[ ${ISOTROPIC_RESAMPLING} != 'NA' ]];
                  then
                    resampleCall="ResampleImage 3 ${PIECEWISE_OUTPUT_PREFIX}TotalWarp${components[$d]}.nii.gz ${PIECEWISE_OUTPUT_PREFIX}TotalWarp${components[$d]}.nii.gz $ISOTROPIC_RESAMPLING_VECTOR 0 0"
                    echo "$resampleCall" >> $qscript
                  fi

                TMP_FILES=( ${TMP_FILES[@]} ${PIECEWISE_OUTPUT_PREFIX}TotalWarp${components[$d]}.nii.gz )
              done

            componentCall="ConvertImage ${DIM} ${PIECEWISE_OUTPUT_PREFIX}TotalWarp ${PIECEWISE_OUTPUT_PREFIX}TotalWarp.nii.gz 9"
            echo "$componentCall" >> $qscript

            PIECEWISE_TRANSFORM_WARPS="${PIECEWISE_TRANSFORM_WARPS} -t ${PIECEWISE_OUTPUT_PREFIX}TotalWarp.nii.gz"

          done

        # Compose all the piecewise warps

        composeCall="${WARP} -d ${DIM} -v 1 \
                         -o [ ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz,1 ] \
                         -r ${SUBSET_CONFIDENCE_MASK} \
                         ${PIECEWISE_TRANSFORM_WARPS} >> $logFile"
        echo "$composeCall" >> $qscript

        # Smooth the displacement field and estimate the inverse

        smoothCall="SmoothDisplacementField ${DIM} ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz 4x4x4 5 3 0 ${SUBSET_CONFIDENCE_MASK} >> $logFile"
        echo "$smoothCall" >> $qscript
        inverseSmoothCall="SmoothDisplacementField ${DIM} ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz 4x4x4 5 3 1 >> $logFile"
        echo "$inverseSmoothCall" >> $qscript

        # Resample to the full size

        xfrmForwardCall="${WARP} -d ${DIM} -v 1 \
                         -o [ ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz,1 ] \
                         -r ${ATLAS_IMAGES[$i]} \
                         -t ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz >> $logFile"
        echo "$xfrmForwardCall" >> $qscript

        xfrmInverseCall="${WARP} -d ${DIM} -v 1 \
                         -o [ ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz,1 ] \
                         -r ${ATLAS_IMAGES[$i]} \
                         -t ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz >> $logFile"
        echo "$xfrmInverseCall" >> $qscript

        # Now do a more refined registration

        registrationCallRefined="antsRegistration -d 3 -v 1 -z 0 -o ${OUTPUT_PREFIX}${BASENAME}_${i}_ -r ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz -r ${INIT_OUTPUT_PREFIX}0GenericAffine.mat -x ${CONFIDENCE_DILATED_MASK}"
        registrationCallRefined="$registrationCallRefined -t BSplineSyN[ 0.15,40,0 ] -m CC[ ${ATLAS_IMAGES[$i]},${TARGET_IMAGE},1,2 ] -c [ 100x40x0,1e-6,10 ] -f 6x4x1 -s 2x1x0 >> $logFile"

        echo "$registrationCallRefined" >> $qscript

        atlasXfrmCall="antsApplyTransforms \
                              -d ${DIM} \
                              --float 1 \
                              -i ${ATLAS_IMAGES[$i]} \
                              -r ${TARGET_IMAGE} \
                              -o ${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz \
                              -n Linear \
                              -t [ ${INIT_OUTPUT_PREFIX}0GenericAffine.mat,1 ] \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_2InverseWarp.nii.gz >> ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"
        echo "$atlasXfrmCall" >> $qscript

        labelXfrmCall="antsApplyTransforms \
                              -d ${DIM} \
                              --float 1 \
                              -i ${ATLAS_LABELS[$i]} \
                              -r ${TARGET_IMAGE} \
                              -o ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz \
                              -n GenericLabel \
                              -t [ ${INIT_OUTPUT_PREFIX}0GenericAffine.mat,1 ] \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz \
                              -t ${OUTPUT_PREFIX}${BASENAME}_${i}_2InverseWarp.nii.gz >> ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"
        echo "$labelXfrmCall" >> $qscript

        TMP_FILES=( ${TMP_FILES[@]} ${INIT_OUTPUT_PREFIX}0GenericAffine.mat )
        TMP_FILES=( ${TMP_FILES[@]} ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalWarp.nii.gz ${OUTPUT_PREFIX}${BASENAME}_${i}_TotalInverseWarp.nii.gz )
        TMP_FILES=( ${TMP_FILES[@]} ${OUTPUT_PREFIX}${BASENAME}_${i}_2Warp.nii.gz ${OUTPUT_PREFIX}${BASENAME}_${i}_2InverseWarp.nii.gz )

        copyImageHeaderCall="CopyImageHeaderInformation \
                             ${TARGET_IMAGE} \
                             ${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz \
                             ${OUTPUT_PREFIX}${BASENAME}_${i}_Warped.nii.gz 1 1 1"
        echo "$copyImageHeaderCall" >> $qscript

        copyLabelsHeaderCall="CopyImageHeaderInformation \
                              ${TARGET_IMAGE} \
                              ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz \
                              ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz 1 1 1"
        echo "$copyLabelsHeaderCall" >> $qscript

        if [[ $KEEP_ALL_IMAGES -eq 0 ]];
          then
            for f in ${TMP_FILES[@]}
              do
                echo "rm -f $f" >> $qscript
             done
          fi

      fi

    if [[ $DOQSUB -eq 1 ]];
      then
        id=`qsub -cwd -S /bin/bash -N antsJlfReg  $QSUB_OPTS $qscript | awk '{print $3}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 4 ]];
      then
        id=`qsub -N antsJlfReg  $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=${REGISTRATION_MEMORY} -l walltime=${REGISTRATION_WALLTIME} $qscript | awk '{print $1}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 3 ]];
      then
        id=`xgrid $XGRID_OPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
        jobIDs="$jobIDs $id"
    elif [[ $DOQSUB -eq 5 ]];
      then
        id=`sbatch --job-name=antsJlfReg${i}  $QSUB_OPTS --nodes=1 --cpus-per-task=1 --time=${REGISTRATION_WALLTIME} --mem=${REGISTRATION_MEMORY} $qscript | rev | cut -f1 -d\ | rev`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 0 ]];
      then
        echo $qscript
        bash $qscript
      fi
done

if [[ $DOQSUB -eq 2 ]];
  then
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF on max ${CORES} cpucores. "
    echo "--------------------------------------------------------------------------------------"
    chmod +x ${OUTPUT_DIR}/job_*.sh
    $PEXEC -j ${CORES} "sh" ${OUTPUT_DIR}/job_*.sh
  fi

jlfCall="antsJointFusion -d ${DIM} -t $TARGET_IMAGE --verbose 1 "

if [[ $DOQSUB -eq 0 ]];
  then
    # Run job locally
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF"
    echo "--------------------------------------------------------------------------------------"

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT} ]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "$maskCall" > $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    echo $qscript2
    bash $qscript2
  fi
if [[ $DOQSUB -eq 1 ]];
  then
    # Run jobs on SGE and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF on SGE cluster. "
    echo "--------------------------------------------------------------------------------------"

    waitForSGEQJobs.pl 1 600 $jobIDs

    # Returns 1 if there are errors
    if [[ ! $? -eq 0 ]];
      then
        echo "qsub submission failed - jobs went into error state"
        exit 1;
      fi

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT}]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "$maskCall" > $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    jobIDs=`qsub -cwd -S /bin/bash -N antsJlf  $QSUB_OPTS $qscript2 | awk '{print $3}'`
    waitForSGEQJobs.pl 1 600 $jobIDs
  fi
if [[ $DOQSUB -eq 4 ]];
  then
    # Run jobs on PBS and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF on PBS cluster. "
    echo "--------------------------------------------------------------------------------------"

    waitForPBSQJobs.pl 1 600 $jobIDs
    # Returns 1 if there are errors
    if [[ ! $? -eq 0 ]];
      then
        echo "qsub submission failed - jobs went into error state"
        exit 1;
      fi

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT}]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "$maskCall" > $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    jobIDs=`qsub -N antsJlf  $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=${JLF_MEMORY} -l walltime=${JLF_WALLTIME} $qscript2 | awk '{print $1}'`
    waitForPBSQJobs.pl 1 600 $jobIDs
  fi

if [[ $DOQSUB -eq 2 ]];
  then
    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT}]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "$maskCall" > $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    sh $qscript2
  fi
if [[ $DOQSUB -eq 3 ]];
  then
    # Run jobs on XGrid and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF on XGrid cluster. Submitted $count jobs "
    echo "--------------------------------------------------------------------------------------"

    waitForXGridJobs.pl -xgridflags "$XGRID_OPTS" -verbose -delay 30 $jobIDs
    # Returns 1 if there are errors
    if [[ ! $? -eq 0 ]];
      then
        echo "XGrid submission failed - jobs went into error state"
        exit 1;
      fi

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT}]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "$maskCall" > $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    sh $qscript2
  fi
if [[ $DOQSUB -eq 5 ]];
  then
    # Run jobs on SLURM and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting JLF on SLURM cluster. "
    echo "--------------------------------------------------------------------------------------"

    waitForSlurmJobs.pl 1 600 $jobIDs
    # Returns 1 if there are errors
    if [[ ! $? -eq 0 ]];
      then
        echo "SLURM submission failed - jobs went into error state"
        exit 1;
      fi

    # Remove the SLURM output files (which are likely to be empty)
    rm -f ${OUTPUT_DIR}/slurm-*.out

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
            EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 2 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    maskCall=''

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        jlfCall="ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      else

        for (( i = 0; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
          do
            jlfCall="${jlfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[$i]} -l ${EXISTING_WARPED_ATLAS_LABELS[$i]}"
          done

        if [[ -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
          then
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz ]"
          else
            jlfCall="${jlfCall} -o [ ${OUTPUT_PREFIX}Labels.nii.gz,${OUTPUT_PREFIX}Intensity.nii.gz,${OUTPUT_POSTERIORS_FORMAT} ]"
          fi

        if [[ ${TARGET_MASK_IMAGE} == 'otsu' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOtsu.nii.gz"
            maskCall="ThresholdImage ${DIM} ${TARGET_IMAGE} ${TARGET_MASK_IMAGE} Otsu 1;"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"

        elif [[ ${TARGET_MASK_IMAGE} == 'or' ]];
          then
            TARGET_MASK_IMAGE="${OUTPUT_PREFIX}TargetMaskImageOr.nii.gz"

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_LABELS[0]} ${EXISTING_WARPED_ATLAS_LABELS[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_LABELS[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_LABELS[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

        elif [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
          then
            MAJORITY_VOTING_IMAGE="${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
            maskCall="ImageMath ${DIM} ${MAJORITY_VOTING_IMAGE} MajorityVoting 0.8 ${EXISTING_WARPED_ATLAS_LABELS[@]};"
            jlfCall="${jlfCall} -x ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting_Mask.nii.gz"

        elif [[ -f ${TARGET_MASK_IMAGE} ]];
          then
            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
          fi
      fi

    qscript2="${OUTPUT_PREFIX}JLF.sh"

    echo "#!/bin/sh" > $qscript2
    echo "$maskCall" >> $qscript2
    echo "$jlfCall" >> $qscript2

    if [[ ${TARGET_MASK_IMAGE} == 'majorityvoting' ]];
      then
        combineCall="ImageMath ${DIM} ${OUTPUT_PREFIX}Labels.nii.gz max ${OUTPUT_PREFIX}Labels.nii.gz ${OUTPUT_PREFIX}TargetMaskImageMajorityVoting.nii.gz"
        echo "$combineCall" >> $qscript2
      fi

    jobIDs=`sbatch --job-name=antsJlf  $QSUB_OPTS --nodes=1 --cpus-per-task=1 --time=${JLF_WALLTIME} --mem=${JLF_MEMORY} $qscript2 | rev | cut -f1 -d\ | rev`
    waitForSlurmJobs.pl 1 600 $jobIDs
  fi

# clean up
rm -f ${OUTPUT_DIR}/job_*.sh
if [[ $KEEP_ALL_IMAGES -eq 0 ]];
  then
    rm -f ${WARPED_ATLAS_IMAGES[@]}
    rm -f ${INVERSE_WARPED_ATLAS_IMAGES[@]}
    rm -f ${WARPED_ATLAS_LABELS[@]}
    rm -f ${AFFINE_FILES[@]}
    rm -f ${WARP_FIELDS[@]}
    rm -f ${INVERSE_WARP_FIELDS[@]}
    rm -f $qscript
    rm -f $qscript2
    rm -f ${OUTPUT_DIR}/slurm-*.out
  fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))
echo
echo "--------------------------------------------------------------------------------------"
if [[ $MAJORITYVOTE -eq 1 ]];
  then
    echo " Done creating: ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz"
  else
    echo " Done creating: ${OUTPUT_PREFIX}Labels.nii.gz"
  fi
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
