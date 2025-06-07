#!/bin/bash

VERSION="0.0.0"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTS=antsRegistration
WARP=antsApplyTransforms
JLF=antsJointFusion
PEXEC=ANTSpexec.sh
SGE=waitForSGEQJobs.pl
PBS=waitForPBSQJobs.pl
XGRID=waitForXGridJobs.pl
SLURM=waitForSlurmJobs.pl

fle_error=0
for FLE in $JLF $ANTS $WARP $PEXEC $SGE $XGRID $PBS $SLURM
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


#assuming .nii.gz as default file type. This is the case for ANTS 1.7 and up

# Initialize variables with defaults, referred to in usage
DIM=3
KEEP_ALL_IMAGES=0
DOQSUB=0
CORES=2
PRECISION=1

XGRID_OPTS=""
SCRIPT_PREPEND=""
QSUB_OPTS=""
TARGET_MASK_IMAGE="majorityvoting"

REGISTRATION_WALLTIME="20:00:00"
REGISTRATION_MEMORY="8gb"
JLF_WALLTIME="20:00:00"
JLF_MEMORY="8gb"

MAJORITYVOTE=0
RUNQUICK=1
TRANSFORM_TYPE="s"


function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -o OutputPrefix <other options> <images>

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3.

     -o:  OutputPrefix:   A prefix that is prepended to all output files.

     -t:  TargetImage:    Target image to be labeled.

     -g:  Atlas:          Atlas to be warped to target image.

     -l:  Labels:         Labels corresponding to atlas (cf -g).

Optional arguments:

     -m:  Majority vote:  Use majority vote instead of joint label fusion (default = ${MAJORITYVOTE}).

     -k:  Keep files:     Keep warped atlas and label files (default = ${KEEP_ALL_IMAGES}).

     -c:  Control for parallel computation (default 0) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM.

     -j:  Number of cpu cores to use (default ${CORES}; -- requires "-c 2").

     -r:  qsub options

     -q:  Use quick registration parameters:  Either 0 or 1 (default = ${RUNQUICK}).

     -p:  Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f:  Float precision: Use float precision (default = ${PRECISION}) -- 0 == double, 1 == float.

     -u:  Registration walltime (default = ${REGISTRATION_WALLTIME}):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Registration memory limit (default = ${REGISTRATION_MEMORY}):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -w:  JLF walltime (default = ${JLF_WALLTIME}):  Option for PBS/SLURM qsub specifying requested time
          for the joint label fusion call.

     -z:  JLF Memory limit (default = ${JLF_MEMORY}):  Option for PBS/SLURM qsub specifying requested memory
          for the joint label fusion call.

     -y:  transform type (default = \'${TRANSFORM_TYPE}\')
        t: translation
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        sr: rigid + deformable syn
        so: deformable syn only
        b: rigid + affine + deformable b-spline syn
        br: rigid + deformable b-spline syn
        bo: deformable b-spline syn only

     -x:  Target mask image (default = \'${TARGET_MASK_IMAGE}\')
        otsu: use otsu thresholding to define foreground/background
        or: 'or' all the warped atlas images to defined foreground/background
        majorityvoting: perform a voxelwise label voting.  If >= 80% of the warped atlases agree at that
                        voxel, we keep that voted label at that voxel and *do not* perform JLF.  Note that
                        the 80% threshold is hard-coded but can be easily changed in the script.
        <filename>: a user-specified mask
        none: don't use a mask

Example:

`basename $0` -d 3 -t target.nii.gz -o malf \\
              -p malfPosteriors%04d.nii.gz \\
              -g atlas1.nii.gz -l labels1.nii.gz \\
              -g atlas2.nii.gz -l labels2.nii.gz \\
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
--------------------------------------------------------------------------------------
script by Nick Tustison
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

function Help {
    cat <<HELP

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

     -m:  Majority vote:  Use majority vote instead of joint label fusion (default = ${MAJORITYVOTE}).

     -k:  Keep files:     Keep warped atlas and label files (default = ${KEEP_ALL_IMAGES}).

     -c:  Control for parallel computation (default 0) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM.

     -j:  Number of cpu cores to use (default ${CORES}; -- requires "-c 2").

     -r:  qsub options

     -q:  Use quick registration parameters:  Either 0 or 1 (default = ${RUNQUICK}).

     -p:  Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f:  Float precision: Use float precision (default = ${PRECISION}) -- 0 == double, 1 == float.

     -u:  Registration walltime (default = ${REGISTRATION_WALLTIME}):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Registration memory limit (default = ${REGISTRATION_MEMORY}):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -w:  JLF walltime (default = ${JLF_WALLTIME}):  Option for PBS/SLURM qsub specifying requested time
          for the joint label fusion call.

     -z:  JLF Memory limit (default = ${JLF_MEMORY}):  Option for PBS/SLURM qsub specifying requested memory
          for the joint label fusion call.

     -y:  transform type (default = \'${TRANSFORM_TYPE}\')
        t: translation
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        sr: rigid + deformable syn
        so: deformable syn only
        b: rigid + affine + deformable b-spline syn
        br: rigid + deformable b-spline syn
        bo: deformable b-spline syn only

     -x:  Target mask image (default = \'${TARGET_MASK_IMAGE}\')
        otsu: use otsu thresholding to define foreground/background
        or: 'or' all the warped atlas images to defined foreground/background
        majorityvoting: perform a voxelwise label voting.  If >= 80% of the warped atlases agree at that
                        voxel, we keep that voted label at that voxel and *do not* perform JLF.  Note that
                        the 80% threshold is hard-coded but can be easily changed in the script.
        <filename>: a user-specified mask
        none: don't use a mask

Requirements:

This scripts relies on the following scripts in your $PATH. The script
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

OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"
OUTPUT_POSTERIORS_FORMAT=''

TARGET_IMAGE=''
ATLAS_IMAGES=()
ATLAS_LABELS=()
TRANSFORM='s'


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
if [[ $# -eq 0 ]];
  then
    Usage >&2
  fi


# Provide output for Help
if [[ "$1" == "-h" ]];
  then
    Help >&2
  fi

# reading command line arguments
while getopts "c:d:f:g:h:j:k:l:m:o:p:q:r:t:u:v:w:x:y:z:" OPT
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
    INVERSE_WARPED_ATLAS_IMAGES[${#INVERSE_WARPED_ATLAS_IMAGES[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_InverseWarped.nii.gz"
    WARPED_ATLAS_LABELS[${#WARPED_ATLAS_LABELS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz"
    WARP_FIELDS[${#WARP_FIELDS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_1Warp.nii.gz"
    INVERSE_WARP_FIELDS[${#INVERSE_WARP_FIELDS[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_1InverseWarp.nii.gz"
    AFFINE_FILES[${#AFFINE_FILES[@]}]="${OUTPUT_PREFIX}${BASENAME}_${i}_0GenericAffine.mat"

    if [[ -f "${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz" ]];
      then
        echo ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz already exists.
        rm -f $qscript
        continue
      fi

    regcall=antsRegistrationSyN.sh
    if [[ $RUNQUICK -eq 1 ]];
      then
        regcall=antsRegistrationSyNQuick.sh
      fi
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
                          -n NearestNeighbor \
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

    if [[ $DOQSUB -eq 1 ]];
      then
        id=`qsub -cwd -S /bin/bash -N antsJlfReg $QSUB_OPTS $qscript | awk '{print $3}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 4 ]];
      then
        id=`qsub -N antsJlfReg $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=${REGISTRATION_MEMORY} -l walltime=${REGISTRATION_WALLTIME} $qscript | awk '{print $1}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 3 ]];
      then
        id=`xgrid $XGRID_OPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
        jobIDs="$jobIDs $id"
    elif [[ $DOQSUB -eq 5 ]];
      then
        id=`sbatch --job-name=antsJlfReg${i} $QSUB_OPTS --nodes=1 --cpus-per-task=1 --time=${REGISTRATION_WALLTIME} --mem=${REGISTRATION_MEMORY} $qscript | rev | cut -f1 -d\ | rev`
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

            maskCall="ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${EXISTING_WARPED_ATLAS_IMAGES[0]} ${EXISTING_WARPED_ATLAS_IMAGES[1]};"
            for (( i = 2; i < ${#EXISTING_WARPED_ATLAS_IMAGES[@]}; i++ ))
              do
                maskCall="${maskCall} ImageMath ${DIM} ${TARGET_MASK_IMAGE} max ${TARGET_MASK_IMAGE} ${EXISTING_WARPED_ATLAS_IMAGES[$i]};"
              done
            maskCall="${maskCall} ThresholdImage ${DIM} ${TARGET_MASK_IMAGE} ${TARGET_MASK_IMAGE} 0 0 0 1"

            jlfCall="${jlfCall} -x ${TARGET_MASK_IMAGE}"
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

    jobIDs=`qsub -cwd -S /bin/bash -N antsJlf $QSUB_OPTS $qscript2 | awk '{print $3}'`
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

    jobIDs=`qsub -N antsJlf $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=${JLF_MEMORY} -l walltime=${JLF_WALLTIME} $qscript2 | awk '{print $1}'`
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

    jobIDs=`sbatch --job-name=antsJlf $QSUB_OPTS --nodes=1 --cpus-per-task=1 --time=${JLF_WALLTIME} --mem=${JLF_MEMORY} $qscript2 | rev | cut -f1 -d\ | rev`
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
