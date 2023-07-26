#!/bin/bash

VERSION="0.0.0 test"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

ANTS=antsRegistration

if ! command -v ${ANTS} &> /dev/null
  then
    echo "antsRegistration program can't be found. Please (re)define \$PATH in your environment."
    exit
  fi

function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -f FixedImage  -o OutputPrefix -m MovingImages*

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)

     -f:  ND fixed image

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = 1)

     -t:  timespacing

     -r:  number of repeats for each time point

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -x:  mask for the fixed image space

     -p:  precision type (default = 'd')
        f: float
        d: double

    -m:  Moving images

Example:

`basename $0` -d 3 -f fixedImage.nii.gz  -o output -m movingImage*.nii.gz

--------------------------------------------------------------------------------------
ANTs was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

script by Nick Tustison, BB Avants

USAGE
    exit 1
}

function Help {
    cat <<HELP

Usage:

`basename $0` -d ImageDimension -f FixedImage  -o OutputPrefix -m MovingImage*

Example Case:

`basename $0` -d 3 -f fixedImage.nii.gz -o output -m movingImage*.nii.gz

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)

     -f:  Fixed image or source image or reference image

     -m:  Moving images

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = 1)

     -t:  timespacing

     -j:  use histogram matching

     -r:  number of repeats for each time point

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -x:  mask for the fixed image space

     -p:  precision type (default = 'd')
        f: float
        d: double

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
 Dimensionality:           $DIM
 Output name prefix:       $OUTPUTNAME
 Fixed images:             ${FIXEDIMAGES[@]}
 Moving images:            ${MOVINGIMAGES[@]}
 Number of threads:        $NUMBEROFTHREADS
 Spline distance:          $SPLINEDISTANCE
 Transform type:           $TRANSFORMTYPE
 CC radius:                $nrepeats
 Precision:                $PRECISIONTYPE
======================================================================================
REPORTMAPPINGPARAMETERS
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
nrepeats=2
MASK=0
USEHISTOGRAMMATCHING=0
timespacing=1
# reading command line arguments
while getopts "d:f:h:m:j:n:o:p:r:s:t:x:" OPT
  do
  case $OPT in
      h) #help
   Help
   exit 0
   ;;
      d)  # dimensions
   DIM=$OPTARG
   ;;
      j)  # histogram matching
   USEHISTOGRAMMATCHING=$OPTARG
   ;;
      x)  # inclusive mask
   MASK=$OPTARG
   ;;
      f)  # fixed image
   FIXEDIMAGES[${#FIXEDIMAGES[@]}]=$OPTARG
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
      r)  # n repeats
   nrepeats=$OPTARG
   ;;
      s)  # spline distance
   SPLINEDISTANCE=$OPTARG
   ;;
      t)  # spacing in time
   timespacing=$OPTARG
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

for(( i=0; i<${#MOVINGIMAGES[@]}; i++ ))
  do
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


# what we do below
#
# stack the moving images together to produce ND+1 input
# antsMotionCorr to get affine to the fixed template
# stack fixed template to ND+1
# register with -r warp to the stacked template
# jacobian ...
let DIMP1=$DIM+1
echo $DIMP1
zero=${OUTPUTNAME}zero.nii.gz
zerob=${OUTPUTNAME}zerob.nii.gz
stack=${OUTPUTNAME}stack.nii.gz
nmov=${#MOVINGIMAGES[@]}
let nmov=$nmov-1
MultiplyImages $DIM ${MOVINGIMAGES[0]} 1 $zero
stackmovparam=" $zero $zero "
for mov in ${MOVINGIMAGES[@]} ; do
  for k in `seq 1 $nrepeats ` ; do
    stackmovparam=" $stackmovparam $mov "
  done
done
let nm1=${#MOVINGIMAGES[@]}-1
MultiplyImages $DIM ${MOVINGIMAGES[nm1]} 1 $zerob
stackmovparam=" $stackmovparam $zerob $zerob "
if [[ $DIM == 2 ]] ; then
  StackSlices $stack -1 -1 0 $stackmovparam | grep -v Slice
fi
if [[ $DIM == 3 ]] ; then
  ImageMath $DIMP1 $stack TimeSeriesAssemble $timespacing 0 ${stackmovparam}
fi
ImageMath $DIMP1 $stack SetTimeSpacing $stack $timespacing
nm=${OUTPUTNAME}
MultiplyImages $DIM ${FIXEDIMAGES} 1 $zero
stackparam=" $zero $zero  "
stacktemplate=${OUTPUTNAME}template.nii.gz
for mov in ${MOVINGIMAGES[@]} ; do
  for k in `seq 1 $nrepeats ` ; do
    stackparam=" $stackparam ${FIXEDIMAGES} "
  done
done
stackparam=" $stackparam $zero $zero  "
if [[ $DIM == 2 ]] ; then
  StackSlices $stacktemplate -1 -1 0 ${stackparam}
fi
if [[ $DIM == 3 ]] ; then
  ImageMath $DIMP1 $stacktemplate TimeSeriesAssemble  $timespacing 0  ${stackparam}
fi
ImageMath $DIMP1 $stacktemplate SetTimeSpacing $stacktemplate $timespacing
# echo $stackparam
# echo $stackmovparam
# echo $nrepeats $timespacing
if [[ $DIM == 2 ]] ; then rxt="1x1x0"; fi
if [[ $DIM == 3 ]] ; then rxt="1x1x1x0"; fi
antsMotionCorr  -d $DIM \
  -o [ ${nm}aff, ${nm}aff.nii.gz,${nm}_affavg.nii.gz ] \
  -m MI[ ${FIXEDIMAGES}, ${stack}, 1 , 20, Regular, 0.1 ] \
  -t rigid[ 0.1 ] -u 1 -e 1 -s 4x2x1x0 -f 6x4x2x1 \
  -i 100x100x0x0 \
  -m MI[ ${FIXEDIMAGES}, ${stack}, 1 , 20, Regular, 0.2 ] \
  -t Affine[ 0.1 ] -u 1 -e 1 -s 4x2x1x0 -f 6x4x2x1 \
  -i 100x100x100x15 \
  -n ${#MOVINGIMAGES[@]} -w 1 --verbose 1
ImageMath $DIMP1 ${nm}affWarp.nii.gz SetTimeSpacingWarp ${nm}affWarp.nii.gz $timespacing
ImageMath $DIMP1 ${nm}affInverseWarp.nii.gz SetTimeSpacingWarp ${nm}affInverseWarp.nii.gz $timespacing
# if below does not work - have to drop the -r
# and put the aff version into the metric
antsRegistration -d $DIMP1 -r ${nm}affWarp.nii.gz  \
 -c [ 100x70x50x10,1e-6,10 ] \
 -f 6x4x2x1               \
 -s 3x2x1x0vox            \
 -m CC[ $stacktemplate, ${nm}stack.nii.gz, 1, 2 ] \
 -t SyN[ 0.15,3,0 ] --restrict-deformation $rxt \
 -o [ ${nm},${nm}diffeoWarped.nii.gz,${nm}diffeoInvWarped.nii.gz ] -z 0
 CreateJacobianDeterminantImage $DIMP1 ${nm}1Warp.nii.gz ${nm}0logjacobian.nii.gz 1 1
###############################
#
# Restore original number of threads
#
###############################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ORIGINALNUMBEROFTHREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
