#!/bin/bash

VERSION="0.0.0"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

function setPath {
    cat <<SETPATH

--------------------------------------------------------------------------------------
Error locating ANTS
--------------------------------------------------------------------------------------
It seems that the ANTSPATH environment variable is not set. Please add the ANTSPATH
variable. This can be achieved by editing the .bash_profile in the home directory.
Add:

ANTSPATH=/home/yourname/bin/ants/

Or the correct location of the ANTS binaries.

Alternatively, edit this script ( `basename $0` ) to set up this parameter correctly.

SETPATH
    exit 1
}

# Uncomment the line below in case you have not set the ANTSPATH variable in your environment.
# export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS

#ANTSPATH=YOURANTSPATH
if [[ ${#ANTSPATH} -le 3 ]];
  then
    setPath >&2
  fi

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTS=${ANTSPATH}/antsRegistration
WARP=${ANTSPATH}/antsApplyTransforms
MALF=${ANTSPATH}/jointfusion
PEXEC=${ANTSPATH}ANTSpexec.sh
SGE=${ANTSPATH}waitForSGEQJobs.pl
PBS=${ANTSPATH}waitForPBSQJobs.pl
XGRID=${ANTSPATH}waitForXGridJobs.pl
SLURM=${ANTSPATH}/waitForSlurmJobs.pl

fle_error=0
for FLE in $MALF $ANTS $WARP $PEXEC $SGE $XGRID $PBS $SLURM
  do
  if [[ ! -x $FLE ]];
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

     -m:  Majority vote:  Use majority vote instead of joint label fusion (default = 0).

     -k:  Keep files:     Keep warped atlas and label files (default = 0).

     -c:  Control for parallel computation (default 0) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM.

     -j: Number of cpu cores to use (default 2; -- requires "-c 2").

     -q: Use quick registration parameters:  Either 0 or 1 (default = 1).

     -p: Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f: Float precision: Use float precision (default = 1) -- 0 == double, 1 == float.

     -x: Target mask image:  Used to check the quality of registrations, if available.

     -z: Dice threshold for target mask image and warped labels (default = 0.85).

Example:

`basename $0` -d 3 -t target.nii.gz -o malf \
              -p malfPosteriors%04d.nii.gz \
              -g atlas1.nii.gz -l labels1.nii.gz \
              -g atlas2.nii.gz -l labels2.nii.gz \
              -g atlas3.nii.gz -l labels3.nii.gz

--------------------------------------------------------------------------------------
MALF was created by:
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

`basename $0` will propagate labels from a set of pre-labeled atlases using the MALF
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

     -j: Number of cpu cores to use (default 2; -- requires "-c 2").

     -q: Use quick registration parameters:  Either 0 or 1 (default = 1).

     -p: Save posteriors:  Save posteriors in specified c-style format e.g. posterior%04d.nii.gz
                           Need to specify output directory.

     -f: Float precision: Use float precision (default = 1) -- 0 == double, 1 == float.

     -x: Target mask image:  Used to check the quality of registrations, if available.

     -z: Dice threshold for target mask image and warped labels (default = 0.85).

Requirements:

This scripts relies on the following scripts in your $ANTSPATH directory. The script
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
MALF was created by:
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
 ANTSPATH is $ANTSPATH

 Dimensionality:           $DIM
 Output prefix:            $OUTPUT_PREFIX
 Posteriors format:        $OUTPUT_POSTERIORS_FORMAT
 Target image:             $TARGET_IMAGE
 Atlas images:             ${ATLAS_IMAGES[@]}
 Atlas labels:             ${ATLAS_LABELS[@]}

 Keep all images:          $KEEP_ALL_IMAGES
 Processing type:          $DOQSUB
 Number of cpu cores:      $CORES
--------------------------------------------------------------------------------------
REPORTPARAMETERS
}

cleanup()
{
  echo "\n*** Performing cleanup, please wait ***\n"

# 1st attempt to kill all remaining processes
# put all related processes in array
runningANTSpids=( `ps -C antsRegistration -C ImageMath| awk '{ printf "%s\n", $1 ; }'` )

# kill these processes, skip the first since it is text and not a PID
for ((i = 1; i < ${#runningANTSpids[@]} ; i++))
  do
  echo "killing:  ${runningANTSpids[${i}]}"
  kill ${runningANTSpids[${i}]}
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


############
#
# Notice of eventual deprecation
#
echo "  ***********************************************************"
echo "  * Note that this script is slated to be deprecated.        "
echo "  * We recommend using the script antsJointLabelFusion.sh    "
echo "  * which has the same options.  Press any key to continue or"
echo "  * control-c to exit.                                       "
echo "  ***********************************************************"

sed -n q </dev/tty

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

KEEP_ALL_IMAGES=0
DOQSUB=0
CORES=1
PRECISION=0

XGRID_OPTS=""
SCRIPT_PREPEND=""
QSUB_OPTS=""
TARGET_MASK_IMAGE=""
DICE_THRESHOLD=0.85

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
# reading command line arguments
while getopts "c:d:f:g:h:j:k:l:m:o:p:t:q:x:z:" OPT
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
      x)
   TARGET_MASK_IMAGE=$OPTARG
   ;;
      z)
   DICE_THRESHOLD=$OPTARG
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

mkdir ${OUTPUT_DIR}

##########################################################################
#
# Perform MALF labeling by
#  1) registering all atlases to target image
#  2) call 'jointfusion'
#
##########################################################################

echo
echo "--------------------------------------------------------------------------------------"
echo " Start MALFization"
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

    regcall=${ANTSPATH}/antsRegistrationSyN.sh
    if [[ $RUNQUICK -eq 1 ]];
      then
        regcall=${ANTSPATH}/antsRegistrationSyNQuick.sh
      fi
    registrationCall="$regcall \
                          -d ${DIM} \
                          -p ${PRECISIONFLAG} \
                          -j 1 \
                          -f ${TARGET_IMAGE} \
                          -m ${ATLAS_IMAGES[$i]} \
                          -o ${OUTPUT_PREFIX}${BASENAME}_${i}_ > ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"

    labelXfrmCall="${ANTSPATH}/antsApplyTransforms \
                          -d ${DIM} \
                          --float 1 \
                          -i ${ATLAS_LABELS[$i]} \
                          -r ${TARGET_IMAGE} \
                          -o ${OUTPUT_PREFIX}${BASENAME}_${i}_WarpedLabels.nii.gz \
                          -n NearestNeighbor \
                          -t ${OUTPUT_PREFIX}${BASENAME}_${i}_1Warp.nii.gz \
                          -t ${OUTPUT_PREFIX}${BASENAME}_${i}_0GenericAffine.mat >> ${OUTPUT_PREFIX}${BASENAME}_${i}_log.txt"

    rm -f $qscript

    if [[ $DOQSUB -eq 5 ]];
      then
        # SLURM job scripts must start with a shebang
        echo '#!/bin/sh' > $qscript
      fi

    echo "$registrationCall" >> $qscript
    echo "$labelXfrmCall" >> $qscript

    if [[ $DOQSUB -eq 1 ]];
      then
        id=`qsub -cwd -S /bin/bash -N antsMalfReg -v ANTSPATH=$ANTSPATH $QSUB_OPTS $qscript | awk '{print $3}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 4 ]];
      then
        id=`qsub -N antsMalfReg -v ANTSPATH=$ANTSPATH $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=8gb -l walltime=20:00:00 $qscript | awk '{print $1}'`
        jobIDs="$jobIDs $id"
        sleep 0.5
    elif [[ $DOQSUB -eq 3 ]];
      then
        id=`xgrid $XGRID_OPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
        jobIDs="$jobIDs $id"
    elif [[ $DOQSUB -eq 5 ]];
      then
        id=`sbatch --job-name=antsMalfReg${i} --export=ANTSPATH=$ANTSPATH $QSUBOPTS --nodes=1 --cpus-per-task=1 --time=20:00:00 --mem=8192M $qscript | rev | cut -f1 -d\ | rev`
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
    echo " Starting MALF on max ${CORES} cpucores. "
    echo "--------------------------------------------------------------------------------------"
    chmod +x ${OUTPUT_DIR}/job_*.sh
    $PEXEC -j ${CORES} "sh" ${OUTPUT_DIR}/job_*.sh
  fi

malfCall="${ANTSPATH}/jointfusion ${DIM} 1 -m Joint[0.1,2] -tg $TARGET_IMAGE "
if [[ ! -z "${OUTPUT_POSTERIORS_FORMAT}" ]];
  then
    malfCall="${malfCall} -p ${OUTPUT_POSTERIORS_FORMAT}"
  fi

if [[ $DOQSUB -eq 0 ]];
  then
    # Run job locally
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting MALF"
    echo "--------------------------------------------------------------------------------------"

    EXISTING_WARPED_ATLAS_IMAGES=()
    EXISTING_WARPED_ATLAS_LABELS=()
    for (( i = 0; i < ${#WARPED_ATLAS_IMAGES[@]}; i++ ))
      do
        echo ${WARPED_ATLAS_IMAGES[$i]}
        if [[ -f ${WARPED_ATLAS_IMAGES[$i]} ]] && [[ -f ${WARPED_ATLAS_LABELS[$i]} ]];
          then
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then
                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
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

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "$malfCall" > $qscript2

    echo $qscript2
    bash $qscript2
  fi
if [[ $DOQSUB -eq 1 ]];
  then
    # Run jobs on SGE and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting MALF on SGE cluster. "
    echo "--------------------------------------------------------------------------------------"

    ${ANTSPATH}/waitForSGEQJobs.pl 1 600 $jobIDs

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
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then
                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 3 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "$malfCall" > $qscript2

    jobIDs=`qsub -cwd -S /bin/bash -N antsMalf -v ANTSPATH=$ANTSPATH $QSUB_OPTS $qscript2 | awk '{print $3}'`
    ${ANTSPATH}/waitForSGEQJobs.pl 1 600 $jobIDs
  fi
if [[ $DOQSUB -eq 4 ]];
  then
    # Run jobs on PBS and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting MALF on PBS cluster. "
    echo "--------------------------------------------------------------------------------------"

    ${ANTSPATH}/waitForPBSQJobs.pl 1 600 $jobIDs
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
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then
                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 3 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "$malfCall" > $qscript2

    jobIDs=`qsub -N antsMalf -v ANTSPATH=$ANTSPATH $QSUB_OPTS -q nopreempt -l nodes=1:ppn=1 -l mem=8gb -l walltime=30:00:00 $qscript2 | awk '{print $1}'`
    ${ANTSPATH}/waitForPBSQJobs.pl 1 600 $jobIDs
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
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then

                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 3 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "$malfCall" > $qscript2

    sh $qscript2
  fi
if [[ $DOQSUB -eq 3 ]];
  then
    # Run jobs on XGrid and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting MALF on XGrid cluster. Submitted $count jobs "
    echo "--------------------------------------------------------------------------------------"

    ${ANTSPATH}/waitForXGridJobs.pl -xgridflags "$XGRID_OPTS" -verbose -delay 30 $jobIDs
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
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then
                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 3 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "$malfCall" > $qscript2

    sh $qscript2
  fi
if [[ $DOQSUB -eq 5 ]];
  then
    # Run jobs on SLURM and wait to finish
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting MALF on SLURM cluster. "
    echo "--------------------------------------------------------------------------------------"

    ${ANTSPATH}/waitForSlurmJobs.pl 1 600 $jobIDs
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
            if [[ -f ${TARGET_MASK_IMAGE} ]];
              then
                TMP_WARPED_ATLAS_LABEL_MASK=${OUTPUT_PREFIX}WarpedAtlasLabelMask.nii.gz
                ${ANTSPATH}/ThresholdImage ${DIM} ${WARPED_ATLAS_LABELS[$i]} ${TMP_WARPED_ATLAS_LABEL_MASK} 0 0 0 1

                OVERLAP_MEASURES=( `${ANTSPATH}/LabelOverlapMeasures ${DIM} ${TMP_WARPED_ATLAS_LABEL_MASK} ${TARGET_MASK_IMAGE} 1` )
                TOKENS=( ${OVERLAP_MEASURES[1]//,/\ } )
                DICE_OVERLAP=${TOKENS[3]}

                if (( $(echo "${DICE_OVERLAP} >= ${DICE_THRESHOLD}" | bc -l) ));
                  then
                    EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                    EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
                  else
                    echo Not including ${WARPED_ATLAS_IMAGES[$i]} \(Dice = ${DICE_OVERLAP}\)
                  fi

                rm -f $TMP_WARPED_ATLAS_LABEL_MASK
              else
                EXISTING_WARPED_ATLAS_IMAGES[${#EXISTING_WARPED_ATLAS_IMAGES[@]}]=${WARPED_ATLAS_IMAGES[$i]}
                EXISTING_WARPED_ATLAS_LABELS[${#EXISTING_WARPED_ATLAS_LABELS[@]}]=${WARPED_ATLAS_LABELS[$i]}
              fi
          fi
      done

    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -lt 2 ]];
      then
        echo "Error:  At least 3 warped image/label pairs needs to exist for jointFusion."
        exit 1
      fi
    if [[ ${#EXISTING_WARPED_ATLAS_LABELS[@]} -ne ${#WARPED_ATLAS_LABELS[@]} ]];
      then
        echo "Warning:  One or more registrations failed."
      fi

    malfCall="${malfCall} -g ${EXISTING_WARPED_ATLAS_IMAGES[@]} -l ${EXISTING_WARPED_ATLAS_LABELS[@]} ${OUTPUT_PREFIX}MalfLabels.nii.gz"

    if [[ $MAJORITYVOTE -eq 1 ]];
      then
        malfCall="${ANTSPATH}/ImageMath ${DIM} ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz MajorityVoting ${EXISTING_WARPED_ATLAS_LABELS[@]} "
      fi

    qscript2="${OUTPUT_PREFIX}MALF.sh"
    echo "#!/bin/sh" > $qscript2
    echo "$malfCall" >> $qscript2

    jobIDs=`sbatch --job-name=antsMalf --export=ANTSPATH=$ANTSPATH $QSUBOPTS --nodes=1 --cpus-per-task=1 --time=30:00:00 --mem=8192M $qscript2 | rev | cut -f1 -d\ | rev`
    ${ANTSPATH}/waitForSlurmJobs.pl 1 600 $jobIDs
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
  fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))
echo
echo "--------------------------------------------------------------------------------------"
if [[ $MAJORITYVOTE -eq 1 ]];
  then
    echo " Done creating: ${OUTPUT_PREFIX}MajorityVotingLabels.nii.gz"
  else
    echo " Done creating: ${OUTPUT_PREFIX}MalfLabels.nii.gz"
  fi
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
