#!/bin/bash

VERSION="0.0.0 test"

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

ANTS=${ANTSPATH}/antsRegistration

if [[ ! -s ${ANTS} ]];
  then
    echo "antsRegistration program can't be found. Please (re)define \$ANTSPATH in your environment."
    exit
  fi

function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -f FixedImage -m MovingImage -o OutputPrefix

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)

     -f:  Fixed image or source image or reference image

     -m:  Moving image or target image

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = 1)

     -t:  transform type (default = 's')
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        b: rigid + affine + deformable b-spline syn

     -r:  radius for cross correlation metric used during SyN stage (default = 4)

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -p:  precision type (default = 'd')
        f: float
        d: double

     -j:  use histogram matching (default = 0)
        0: false
        1: true

     NB:  Multiple image pairs can be specified for registration during the SyN stage.
          Specify additional images using the '-m' and '-f' options.  Note that image
          pair correspondence is given by the order specified on the command line.
          Only the first fixed and moving image pair is used for the linear resgitration
          stages.

Example:

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -o output

--------------------------------------------------------------------------------------
ANTs was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

script by Nick Tustison

USAGE
    exit 1
}

function Help {
    cat <<HELP

Usage:

`basename $0` -d ImageDimension -f FixedImage -m MovingImage -o OutputPrefix

Example Case:

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -o output

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)

     -f:  Fixed image or source image or reference image

     -m:  Moving image or target image

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = 1)

     -t:  transform type (default = 's')
        r: rigid
        a: rigid + affine
        s: rigid + affine + deformable syn
        b: rigid + affine + deformable b-spline syn

     -r:  radius for cross correlation metric used during SyN stage (default = 4)

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -p:  precision type (default = 'd')
        f: float
        d: double

     -j:  use histogram matching (default = 0)
        0: false
        1: true

     NB:  Multiple image pairs can be specified for registration during the SyN stage.
          Specify additional images using the '-m' and '-f' options.  Note that image
          pair correspondence is given by the order specified on the command line.
          Only the first fixed and moving image pair is used for the linear resgitration
          stages.

--------------------------------------------------------------------------------------
Get the latest ANTs version at:
--------------------------------------------------------------------------------------
https://github.com/stnava/ANTs/

--------------------------------------------------------------------------------------
Read the ANTS documentation at:
--------------------------------------------------------------------------------------
http://stnava.github.io/ANTs/

--------------------------------------------------------------------------------------
ANTS was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

Relevent references for this script include:
   * http://www.ncbi.nlm.nih.gov/pubmed/20851191
   * http://www.frontiersin.org/Journal/10.3389/fninf.2013.00039/abstract
--------------------------------------------------------------------------------------
script by Nick Tustison
--------------------------------------------------------------------------------------

HELP
    exit 1
}

function reportMappingParameters {
    cat <<REPORTMAPPINGPARAMETERS

--------------------------------------------------------------------------------------
 Mapping parameters
--------------------------------------------------------------------------------------
 ANTSPATH is $ANTSPATH

 Dimensionality:           $DIM
 Output name prefix:       $OUTPUTNAME
 Fixed images:             ${FIXEDIMAGES[@]}
 Moving images:            ${MOVINGIMAGES[@]}
 Number of threads:        $NUMBEROFTHREADS
 Spline distance:          $SPLINEDISTANCE
 Transform type:           $TRANSFORMTYPE
 CC radius:                $CCRADIUS
 Precision:                $PRECISIONTYPE
 Use histogram matching    $USEHISTOGRAMMATCHING
======================================================================================
REPORTMAPPINGPARAMETERS
}

cleanup()
# example cleanup function
{

  cd ${currentdir}/

  echo "\n*** Performing cleanup, please wait ***\n"

# 1st attempt to kill all remaining processes
# put all related processes in array
runningANTSpids=( `ps -C antsRegistration | awk '{ printf "%s\n", $1 ; }'` )

# debug only
  #echo list 1: ${runningANTSpids[@]}

# kill these processes, skip the first since it is text and not a PID
for (( i = 1; i < ${#runningANTSpids[@]}; i++ ))
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


# Provide output for Help
if [[ "$1" == "-h" || $# -eq 0 ]];
  then
    Help >&2
  fi

#################
#
# default values
#
#################

DIM=3
FIXEDIMAGES=()
MOVINGIMAGES=()
OUTPUTNAME=output
NUMBEROFTHREADS=1
SPLINEDISTANCE=26
TRANSFORMTYPE='s'
PRECISIONTYPE='d'
USEHISTOGRAMMATCHING=0
CCRADIUS=4

# reading command line arguments
while getopts "d:f:h:j:m:n:o:p:r:s:t:" OPT
  do
  case $OPT in
      h) #help
   Help
   exit 0
   ;;
      d)  # dimensions
   DIM=$OPTARG
   ;;
      f)  # fixed image
   FIXEDIMAGES[${#FIXEDIMAGES[@]}]=$OPTARG
   ;;
      j)  # histogram matching
   USEHISTOGRAMMATCHING=$OPTARG
   ;;
      m)  # moving image
   MOVINGIMAGES[${#MOVINGIMAGES[@]}]=$OPTARG
   ;;
      n)  # number of threads
   NUMBEROFTHREADS=$OPTARG
   ;;
      o) #output name prefix
   OUTPUTNAME=$OPTARG
   ;;
      p)  # precision type
   PRECISIONTYPE=$OPTARG
   ;;
      r)  # cc radius
   CCRADIUS=$OPTARG
   ;;
      s)  # spline distance
   SPLINEDISTANCE=$OPTARG
   ;;
      t)  # transform type
   TRANSFORMTYPE=$OPTARG
   ;;
     \?) # getopts issues an error message
   echo "$USAGE" >&2
   exit 1
   ;;
  esac
done

###############################
#
# Check inputs
#
###############################
if [[ ${#FIXEDIMAGES[@]} -ne ${#MOVINGIMAGES[@]} ]];
  then
    echo "Number of fixed images is not equal to the number of moving images."
    exit 1
  fi

for(( i=0; i<${#FIXEDIMAGES[@]}; i++ ))
  do
    if [[ ! -f "${FIXEDIMAGES[$i]}" ]];
      then
        echo "Fixed image '${FIXEDIMAGES[$i]}' does not exist.  See usage: '$0 -h 1'"
        exit 1
      fi
    if [[ ! -f "${MOVINGIMAGES[$i]}" ]];
      then
        echo "Moving image '${MOVINGIMAGES[$i]}' does not exist.  See usage: '$0 -h 1'"
        exit 1
      fi
  done

###############################
#
# Set number of threads
#
###############################

ORIGINALNUMBEROFTHREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NUMBEROFTHREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

##############################
#
# Print out options
#
##############################

reportMappingParameters

##############################
#
# Infer the number of levels based on
# the size of the input fixed image.
#
##############################

ISLARGEIMAGE=0

SIZESTRING=$( ${ANTSPATH}/PrintHeader ${FIXEDIMAGES[0]} 2 )
SIZESTRING="${SIZESTRING%\\n}"
SIZE=( `echo $SIZESTRING | tr 'x' ' '` )

for (( i=0; i<${#SIZE[@]}; i++ ))
  do
    if [[ ${SIZE[$i]} -gt 256 ]];
      then
        ISLARGEIMAGE=1
        break
      fi
  done

##############################
#
# Construct mapping stages
#
##############################

RIGIDCONVERGENCE="[1000x500x250x100,1e-6,10]"
RIGIDSHRINKFACTORS="8x4x2x1"
RIGIDSMOOTHINGSIGMAS="3x2x1x0vox"

AFFINECONVERGENCE="[1000x500x250x100,1e-6,10]"
AFFINESHRINKFACTORS="8x4x2x1"
AFFINESMOOTHINGSIGMAS="3x2x1x0vox"

SYNCONVERGENCE="[100x70x50x20,1e-6,10]"
SYNSHRINKFACTORS="8x4x2x1"
SYNSMOOTHINGSIGMAS="3x2x1x0vox"

if [[ $ISLARGEIMAGE -eq 1 ]];
  then
    RIGIDCONVERGENCE="[1000x500x250x100,1e-6,10]"
    RIGIDSHRINKFACTORS="12x8x4x2"
    RIGIDSMOOTHINGSIGMAS="4x3x2x1vox"

    AFFINECONVERGENCE="[1000x500x250x100,1e-6,10]"
    AFFINESHRINKFACTORS="12x8x4x2"
    AFFINESMOOTHINGSIGMAS="4x3x2x1vox"

    SYNCONVERGENCE="[100x100x70x50x20,1e-6,10]"
    SYNSHRINKFACTORS="10x6x4x2x1"
    SYNSMOOTHINGSIGMAS="5x3x2x1x0vox"
  fi

RIGIDSTAGE="--initial-moving-transform [${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1] \
            --transform Rigid[0.1] \
            --metric MI[${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1,32,Regular,0.25] \
            --convergence $RIGIDCONVERGENCE \
            --shrink-factors $RIGIDSHRINKFACTORS \
            --smoothing-sigmas $RIGIDSMOOTHINGSIGMAS"

AFFINESTAGE="--transform Affine[0.1] \
             --metric MI[${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1,32,Regular,0.25] \
             --convergence $AFFINECONVERGENCE \
             --shrink-factors $AFFINESHRINKFACTORS \
             --smoothing-sigmas $AFFINESMOOTHINGSIGMAS"

SYNMETRICS=''
for(( i=0; i<${#FIXEDIMAGES[@]}; i++ ))
  do
    SYNMETRICS="$SYNMETRICS --metric CC[${FIXEDIMAGES[$i]},${MOVINGIMAGES[$i]},1,${CCRADIUS}]"
  done

SYNSTAGE="${SYNMETRICS} \
          --convergence $SYNCONVERGENCE \
          --shrink-factors $SYNSHRINKFACTORS \
          --smoothing-sigmas $SYNSMOOTHINGSIGMAS"

if [[ $TRANSFORMTYPE == 'b' ]];
  then
    SYNSTAGE="--transform BSplineSyN[0.1,${SPLINEDISTANCE},0,3] \
             $SYNSTAGE"
  fi
if [[ $TRANSFORMTYPE == 's' ]];
  then
    SYNSTAGE="--transform SyN[0.1,3,0] \
             $SYNSTAGE"
  fi

STAGES=''
case "$TRANSFORMTYPE" in
"r")
  STAGES="$RIGIDSTAGE"
  ;;
"a")
  STAGES="$RIGIDSTAGE $AFFINESTAGE"
  ;;
"b" | "s")
  STAGES="$RIGIDSTAGE $AFFINESTAGE $SYNSTAGE"
  ;;
*)
  echo "Transform type '$TRANSFORMTYPE' is not an option.  See usage: '$0 -h 1'"
  exit
  ;;
esac

PRECISION=''
case "$PRECISIONTYPE" in
"f")
  PRECISION="--float 1"
  ;;
"d")
  PRECISION="--float 0"
  ;;
*)
  echo "Precision type '$PRECISIONTYPE' is not an option.  See usage: '$0 -h 1'"
  exit
  ;;
esac

COMMAND="${ANTS} --verbose 1 \
                 --dimensionality $DIM $PRECISION \
                 --output [$OUTPUTNAME,${OUTPUTNAME}Warped.nii.gz] \
                 --interpolation Linear \
                 --winsorize-image-intensities [0.005,0.995] \
                 --use-histogram-matching ${USEHISTOGRAMMATCHING} \
                 $STAGES"

echo " antsRegistration call:"
echo "--------------------------------------------------------------------------------------"
echo ${COMMAND}
echo "--------------------------------------------------------------------------------------"

$COMMAND

###############################
#
# Restore original number of threads
#
###############################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ORIGINALNUMBEROFTHREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
