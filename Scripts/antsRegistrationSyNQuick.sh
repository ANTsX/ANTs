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

`basename $0` -d ImageDimension -f FixedImage -m MovingImage -o OutputPrefix

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)

     -f:  Fixed image(s) or source image(s) or reference image(s)

     -m:  Moving image(s) or target image(s)

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS if defined, otherwise 1)

     -i:  initial transform(s) --- order specified on the command line matters. If not specified, a
          default initialization is used, based on the transform type.

          For "syn only" and "b-spline syn only" transforms, the initial transform is the identity
          matrix. For other transforms, it is a translation that aligns the input centers of mass.

     -t:  transform type (default = 's')
        t: translation (1 stage)
        r: rigid (1 stage)
        a: rigid + affine (2 stages)
        s: rigid + affine + deformable syn (3 stages)
        sr: rigid + deformable syn (2 stages)
        so: deformable syn only (1 stage)
        b: rigid + affine + deformable b-spline syn (3 stages)
        br: rigid + deformable b-spline syn (2 stages)
        bo: deformable b-spline syn only (1 stage)

     -r:  histogram bins for mutual information in SyN stage (default = 32)

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -g:  gradient step size for SyN and B-spline SyN (default = 0.1)

     -x:  mask(s) for the fixed image space, or for the fixed and moving image space in the format
          "fixedMask,MovingMask". Use -x once to specify mask(s) to be used for all stages or use
          -x for each "stage" (cf -t option).  If no mask is to be used for a particular stage,
          the keyword 'NULL' should be used in place of file names.

     -p:  precision type (default = 'd')
        f: float
        d: double

     -j:  use histogram matching (default = 0)
        0: false
        1: true

     -y:  use 'repro' mode for exact reproducibility of output.  Uses GC metric for linear
          stages, CC metric for deformable stages, and a fixed random seed (default = 0).
        0: false
        1: true

     -z:  collapse output transforms (default = 1)
        0: false
        1: true

     -e:  Fix random seed to an int value

     NB:  Multiple image pairs can be specified for registration during the SyN stage.
          Specify additional images using the '-m' and '-f' options.  Note that image
          pair correspondence is given by the order specified on the command line.
          Only the first fixed and moving image pair is used for the linear registration
          stages.

Example:

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -o output

Example with masks:

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -x fixedMask.nii.gz -o output

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -x fixedMask.nii.gz,movingMask.nii.gz -o output

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -t sr -x NULL -x fixedMask.nii.gz  -t -o output


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

     -f:  Fixed image(s) or source image(s) or reference image(s)

     -m:  Moving image(s) or target image(s)

     -o:  OutputPrefix: A prefix that is prepended to all output files.

Optional arguments:

     -n:  Number of threads (default = ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS if defined, otherwise 1)

     -i:  initial transform(s) --- order specified on the command line matters. If not specified, a
          default initialization is used, based on the transform type.

          For "syn only" and "b-spline syn only" transforms, the initial transform is the identity
          matrix. For other transforms, it is a translation that aligns the input centers of mass.

     -t:  transform type (default = 's')
        t: translation (1 stage)
        r: rigid (1 stage)
        a: rigid + affine (2 stages)
        s: rigid + affine + deformable syn (3 stages)
        sr: rigid + deformable syn (2 stages)
        so: deformable syn only (1 stage)
        b: rigid + affine + deformable b-spline syn (3 stages)
        br: rigid + deformable b-spline syn (2 stages)
        bo: deformable b-spline syn only (1 stage)

     -r:  histogram bins for mutual information in SyN stage (default = 32)

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

     -g:  gradient step size for SyN and B-spline SyN (default = 0.1)

     -x:  mask(s) for the fixed image space, or for the fixed and moving image space in the format
          "fixedMask,MovingMask". Use -x once to specify mask(s) to be used for all stages or use
          -x for each "stage" (cf -t option).  If no mask is to be used for a particular stage,
          the keyword 'NULL' should be used in place of file names.

     -p:  precision type (default = 'd')
        f: float
        d: double

     -j:  use histogram matching (default = 0)
        0: false
        1: true

     -y:  use 'repro' mode for exact reproducibility of output.  Uses GC metric for linear
          stages, CC metric for deformable stages, and a fixed random seed (default = 0).
        0: false
        1: true

     -z:  collapse output transforms (default = 1)
        0: false
        1: true

     -e:  Fix random seed to an int value

     NB:  Multiple image pairs can be specified for registration during the SyN stage.
          Specify additional images using the '-m' and '-f' options.  Note that image
          pair correspondence is given by the order specified on the command line.
          Only the first fixed and moving image pair is used for the linear registration
          stages.

--------------------------------------------------------------------------------------
Get the latest ANTs version at:
--------------------------------------------------------------------------------------
https://github.com/ANTsX/ANTs/

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
 Mask images:              ${MASKIMAGES[@]}
 Initial transforms:       ${INITIALTRANSFORMS[@]}
 Number of threads:        $NUMBEROFTHREADS
 Spline distance:          $SPLINEDISTANCE
 Linear gradient step:     $LINEARGRADIENTSTEP
 SyN gradient step:        $SYNGRADIENTSTEP
 Transform type:           $TRANSFORMTYPE
 MI histogram bins:        $NUMBEROFBINS
 Precision:                $PRECISIONTYPE
 Use histogram matching:   $USEHISTOGRAMMATCHING
 Repro                     $REPRO
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
INITIALTRANSFORMS=()
OUTPUTNAME=output
NUMBEROFTHREADS=0
SPLINEDISTANCE=26
LINEARGRADIENTSTEP=0.1
SYNGRADIENTSTEP=0.2
TRANSFORMTYPE='s'
PRECISIONTYPE='d'
NUMBEROFBINS=32
MASKIMAGES=()
USEHISTOGRAMMATCHING=0
COLLAPSEOUTPUTTRANSFORMS=1
RANDOMSEED=0
REPRO=0

# reading command line arguments
while getopts "d:e:f:g:h:i:m:j:n:o:p:r:s:t:x:y:z:" OPT
  do
  case $OPT in
      h) #help
   Help
   exit 0
   ;;
      d)  # dimensions
   DIM=$OPTARG
   ;;
      e)  # seed
   RANDOMSEED=$OPTARG
   ;;
      x)  # inclusive mask
   MASKIMAGES[${#MASKIMAGES[@]}]=$OPTARG
   ;;
      f)  # fixed image
   FIXEDIMAGES[${#FIXEDIMAGES[@]}]=$OPTARG
   ;;
      g)  # SyN gradient step
   SYNGRADIENTSTEP=$OPTARG
   ;;
      j)  # histogram matching
   USEHISTOGRAMMATCHING=$OPTARG
   ;;
      m)  # moving image
   MOVINGIMAGES[${#MOVINGIMAGES[@]}]=$OPTARG
   ;;
      i)  # initial transform
   INITIALTRANSFORMS[${#INITIALTRANSFORMS[@]}]=$OPTARG
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
   NUMBEROFBINS=$OPTARG
   ;;
      s)  # spline distance
   SPLINEDISTANCE=$OPTARG
   ;;
      t)  # transform type
   TRANSFORMTYPE=$OPTARG
   ;;
      y)  # reproducibility
   REPRO=$OPTARG
   ;;
      z)  # collapse output transforms
   COLLAPSEOUTPUTTRANSFORMS=$OPTARG
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

##############################
#
# Mask stuff
#
##############################

NUMBEROFMASKSTAGES=${#MASKIMAGES[@]}

MASKCALL=""
if [[ ${#MASKIMAGES[@]} -gt 0 ]];
  then
    for (( i = 0; i < ${#MASKIMAGES[@]}; i++ ))
      do
        if [[ ${MASKIMAGES[$i]} =~ "," ]]; then
            MASKCALL="${MASKCALL} -x [ ${MASKIMAGES[$i]} ]"
        else
            MASKCALL="${MASKCALL} -x [ ${MASKIMAGES[$i]}, NULL ]"
        fi
      done
  fi


###############################
#
# Set number of threads
#
###############################

ORIGINALNUMBEROFTHREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}

# NUMBEROFTHREADS is > 0 if the option has been set to a positive value
if [[ $NUMBEROFTHREADS -lt 1 ]]
  then
    # Number of threads not set on the command line, try ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
    if [[ $ORIGINALNUMBEROFTHREADS -gt 0 ]]
      then
	NUMBEROFTHREADS=$ORIGINALNUMBEROFTHREADS
    else
	NUMBEROFTHREADS=1
    fi
  fi

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

SIZESTRING=$( PrintHeader ${FIXEDIMAGES[0]} 2 )
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

RIGIDCONVERGENCE="[ 1000x500x250x0,1e-6,10 ]"
RIGIDSHRINKFACTORS="8x4x2x1"
RIGIDSMOOTHINGSIGMAS="3x2x1x0vox"

AFFINECONVERGENCE="[ 1000x500x250x0,1e-6,10 ]"
AFFINESHRINKFACTORS="8x4x2x1"
AFFINESMOOTHINGSIGMAS="3x2x1x0vox"

SYNCONVERGENCE="[ 100x70x50x0,1e-6,10 ]"
SYNSHRINKFACTORS="8x4x2x1"
SYNSMOOTHINGSIGMAS="3x2x1x0vox"

if [[ $ISLARGEIMAGE -eq 1 ]];
  then
    RIGIDCONVERGENCE="[ 1000x500x250x0,1e-6,10 ]"
    RIGIDSHRINKFACTORS="12x8x4x2"
    RIGIDSMOOTHINGSIGMAS="4x3x2x1vox"

    AFFINECONVERGENCE="[ 1000x500x250x0,1e-6,10 ]"
    AFFINESHRINKFACTORS="12x8x4x2"
    AFFINESMOOTHINGSIGMAS="4x3x2x1vox"

    SYNCONVERGENCE="[ 100x100x70x50x0,1e-6,10 ]"
    SYNSHRINKFACTORS="10x6x4x2x1"
    SYNSMOOTHINGSIGMAS="5x3x2x1x0vox"
  fi

LINEARMETRIC="MI"
LINEARMETRICPARAMETER=32

# Precedence for random seeding
# 1. Command line option -e
# 2. Environment variable ANTS_RANDOM_SEED
# 3. Fixed seed = 1 if run in repro mode
# 4. ITK default (system time)

if [[ -n ${ANTS_RANDOM_SEED} ]] && [[ ${RANDOMSEED} -eq 0 ]];
  then
    RANDOMSEED=${ANTS_RANDOM_SEED}
  fi

if [[ $REPRO -eq 1 ]];
  then
    LINEARMETRIC="GC"
    LINEARMETRICPARAMETER=1
    if [[ ${RANDOMSEED} -eq 0 ]];
      then
        RANDOMSEED=1
      fi
  fi

INITIALSTAGE="--initial-moving-transform [ ${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1 ]"

if [[ ${TRANSFORMTYPE} == 'so' ]] || [[ ${TRANSFORMTYPE} == 'bo' ]];
  then
    # Could just set to an empty string but this is more explicit and also keeps transform
    # numbering consistent with other use cases
    #
    # Default to identity if we are doing "syn only" or "b-spline syn only"
    INITIALSTAGE="--initial-moving-transform Identity"
  fi

if [[ ${#INITIALTRANSFORMS[@]} -gt 0 ]];
  then
    INITIALSTAGE=""
    for(( i=0; i<${#INITIALTRANSFORMS[@]}; i++ ))
      do
        INITIALSTAGE="$INITIALSTAGE --initial-moving-transform ${INITIALTRANSFORMS[$i]}"
      done
  fi

tx=Rigid
if [[ $TRANSFORMTYPE == 't' ]] ; then
  tx=Translation
fi

RIGIDSTAGE="--transform ${tx}[ ${LINEARGRADIENTSTEP} ] \
            --metric ${LINEARMETRIC}[ ${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1,${LINEARMETRICPARAMETER},Regular,0.25 ] \
            --convergence $RIGIDCONVERGENCE \
            --shrink-factors $RIGIDSHRINKFACTORS \
            --smoothing-sigmas $RIGIDSMOOTHINGSIGMAS"

AFFINESTAGE="--transform Affine[ ${LINEARGRADIENTSTEP} ] \
             --metric ${LINEARMETRIC}[ ${FIXEDIMAGES[0]},${MOVINGIMAGES[0]},1,${LINEARMETRICPARAMETER},Regular,0.25 ] \
             --convergence $AFFINECONVERGENCE \
             --shrink-factors $AFFINESHRINKFACTORS \
             --smoothing-sigmas $AFFINESMOOTHINGSIGMAS"

SYNMETRICS=''
for(( i=0; i<${#FIXEDIMAGES[@]}; i++ ))
  do
    if [[ REPRO -eq 1 ]]
      then
      SYNMETRICS="$SYNMETRICS --metric CC[ ${FIXEDIMAGES[$i]},${MOVINGIMAGES[$i]},1,2]"
      else
      SYNMETRICS="$SYNMETRICS --metric MI[ ${FIXEDIMAGES[$i]},${MOVINGIMAGES[$i]},1,${NUMBEROFBINS}]"
      fi
  done

SYNSTAGE="${SYNMETRICS} \
          --convergence $SYNCONVERGENCE \
          --shrink-factors $SYNSHRINKFACTORS \
          --smoothing-sigmas $SYNSMOOTHINGSIGMAS"

if [[ $TRANSFORMTYPE == 'sr' ]] || [[ $TRANSFORMTYPE == 'br' ]];
  then
    SYNCONVERGENCE="[ 50x0,1e-6,10 ]"
    SYNSHRINKFACTORS="2x1"
    SYNSMOOTHINGSIGMAS="1x0vox"
          SYNSTAGE="${SYNMETRICS} \
          --convergence $SYNCONVERGENCE \
          --shrink-factors $SYNSHRINKFACTORS \
          --smoothing-sigmas $SYNSMOOTHINGSIGMAS"
  fi

if [[ $TRANSFORMTYPE == 'b' ]] || [[ $TRANSFORMTYPE == 'br' ]] || [[ $TRANSFORMTYPE == 'bo' ]];
  then
    SYNSTAGE="--transform BSplineSyN[ ${SYNGRADIENTSTEP},${SPLINEDISTANCE},0,3 ] \
             $SYNSTAGE"
  fi

if [[ $TRANSFORMTYPE == 's' ]] || [[ $TRANSFORMTYPE == 'sr' ]] || [[ $TRANSFORMTYPE == 'so' ]];
  then
    SYNSTAGE="--transform SyN[ ${SYNGRADIENTSTEP},3,0 ] \
             $SYNSTAGE"
  fi

NUMBEROFREGISTRATIONSTAGES=0
STAGES=''
case "$TRANSFORMTYPE" in
"r" | "t")
  STAGES="$INITIALSTAGE $RIGIDSTAGE"
  NUMBEROFREGISTRATIONSTAGES=1
  ;;
"a")
  STAGES="$INITIALSTAGE $RIGIDSTAGE $AFFINESTAGE"
  NUMBEROFREGISTRATIONSTAGES=2
  ;;
"b" | "s")
  STAGES="$INITIALSTAGE $RIGIDSTAGE $AFFINESTAGE $SYNSTAGE"
  NUMBEROFREGISTRATIONSTAGES=3
  ;;
"br" | "sr")
  STAGES="$INITIALSTAGE $RIGIDSTAGE  $SYNSTAGE"
  NUMBEROFREGISTRATIONSTAGES=2
  ;;
"bo" | "so")
  STAGES="$INITIALSTAGE $SYNSTAGE"
  NUMBEROFREGISTRATIONSTAGES=1
  ;;
*)
  echo "Transform type '$TRANSFORMTYPE' is not an option.  See usage: '$0 -h 1'"
  exit
  ;;
esac

if [[ $NUMBEROFMASKSTAGES -ne 0 && $NUMBEROFMASKSTAGES -ne 1 && $NUMBEROFMASKSTAGES -ne $NUMBEROFREGISTRATIONSTAGES ]];
  then
    echo "The specified number of mask images is not correct.  Please see help menu."
    exit
  fi

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

RANDOMOPT=""

if [[ ! $RANDOMSEED -eq 0 ]]; then
    RANDOMOPT=" --random-seed $RANDOMSEED "
fi

COMMAND="${ANTS} --verbose 1 $RANDOMOPT \
                 --dimensionality $DIM $PRECISION \
                 --collapse-output-transforms $COLLAPSEOUTPUTTRANSFORMS \
                 --output [ $OUTPUTNAME,${OUTPUTNAME}Warped.nii.gz,${OUTPUTNAME}InverseWarped.nii.gz ] \
                 --interpolation Linear \
                 --use-histogram-matching ${USEHISTOGRAMMATCHING} \
                 --winsorize-image-intensities [ 0.005,0.995 ] \
                 $MASKCALL \
                 $STAGES"

echo " antsRegistration call:"
echo "--------------------------------------------------------------------------------------"
echo ${COMMAND}
echo "--------------------------------------------------------------------------------------"

$COMMAND

echo " Registration finished. The antsRegistration call was:"
echo "--------------------------------------------------------------------------------------"
echo ${COMMAND}
echo "--------------------------------------------------------------------------------------"
echo "Moving image resampled into fixed space: ${OUTPUTNAME}Warped.nii.gz"
echo "Fixed image resampled into moving space: ${OUTPUTNAME}InverseWarped.nii.gz"
echo "--------------------------------------------------------------------------------------"

###############################
#
# Restore original number of threads
#
###############################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ORIGINALNUMBEROFTHREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
