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

PROGRAMS[0]=ExtractSliceFromImage
PROGRAMS[1]=PrintHeader
PROGRAMS[2]=TileImages
PROGRAMS[3]=DenoiseImage
PROGRAMS[4]=ImageMath
PROGRAMS[5]=N4BiasFieldCorrection

for (( i = 0; i < ${#PROGRAMS[@]}; i++ ))
  do
    if [[ ! -s "${ANTSPATH}/${PROGRAMS[$i]}" ]];
      then
        echo "The ${PROGRAMS[$i]} program can't be found. Please (re)define \$ANTSPATH in your environment."
        exit
      fi
  done

function Usage {
    cat <<USAGE

Usage:

`basename $0` outputImage operation whichDirection inputImage [optional parameters]

Operations:

     * ExtractAllSlices:

     * Convolve (ImageMath):

     * DenoiseImage:

     * N4BiasFieldCorrection:

     * RescaleImage (ImageMath):

     * TruncateImageIntensity (ImageMath):


Example:

`basename $0` output.nii.gz DenoiseImage 2 input.nii.gz

--------------------------------------------------------------------------------------
ANTs was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

script by Nick Tustison

--------------------------------------------------------------------------------------
Get the latest ANTs version at:
--------------------------------------------------------------------------------------
https://github.com/stnava/ANTs/

--------------------------------------------------------------------------------------
Read the ANTS documentation at:
--------------------------------------------------------------------------------------
http://stnava.github.io/ANTs/


USAGE
    exit 1
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
    Usage >&2
  fi

if [ $# -lt 4 ];
  then
    echo "Illegal number of parameters ($#)"
    Usage >&2
fi

OUTPUT_IMAGE=$1
OPERATION=$2
WHICH_DIRECTION=$3
INPUT_IMAGE=$4

PARAMETERS[0]=$5
PARAMETERS[1]=$6
PARAMETERS[2]=$7
PARAMETERS[3]=$8
PARAMETERS[4]=$9
PARAMETERS[5]=$10
PARAMETERS[6]=$11

SIZE_STRING=$( PrintHeader $INPUT_IMAGE 2 )
SIZE=( ${SIZE_STRING//x/ } )

if [[ ${#SIZE[@]} -ne 3 ]];
  then
    echo "Error:  The input image, $INPUT_IMAGE, is not 3-D."
    exit
  fi

NUMBER_OF_SLICES=${SIZE[$WHICH_DIRECTION]}
TMP_OUTPUT_DIR=$( mktemp -d )
TMP_OUTPUT_FILE=${TMP_OUTPUT_DIR}/tmp${RANDOM}.nii.gz

ALL_OUTPUT_SLICES=()
for (( i = 0; i < $NUMBER_OF_SLICES; i++ ))
  do
    echo "$OPERATION (direction $WHICH_DIRECTION, slice $i of $NUMBER_OF_SLICES)"

    OUTPUT_DIR=`dirname $OUTPUT_IMAGE`

    OUTPUT_SLICE=`basename $OUTPUT_IMAGE`
    OUTPUT_SLICE=${OUTPUT_SLICE/\.mha/}
    OUTPUT_SLICE=${OUTPUT_SLICE/\.nii\.gz/}
    OUTPUT_SLICE=${OUTPUT_SLICE/\.nii/}
    OUTPUT_SLICE=${OUTPUT_SLICE/\.nrrd/}

    OUTPUT_SLICE="${TMP_OUTPUT_DIR}/${OUTPUT_SLICE}Slice${i}.nii.gz"

    ${ANTSPATH}/ExtractSliceFromImage 3 $INPUT_IMAGE $OUTPUT_SLICE $WHICH_DIRECTION $i

    ALL_OUTPUT_SLICES[$i]=$OUTPUT_SLICE

    case "$OPERATION" in
    "ExtractAllSlices")
      mv ${OUTPUT_SLICE} ${OUTPUT_DIR}
      ;;
    "Convolve")
      ${ANTSPATH}/ImageMath 2 $OUTPUT_SLICE Convolve $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]}
      ;;
    "DenoiseImage")
      ${ANTSPATH}/DenoiseImage -d 2 -i $OUTPUT_SLICE -o $OUTPUT_SLICE -v 0
      ;;
    "N4BiasFieldCorrection")
      ${ANTSPATH}/N4BiasFieldCorrection -d 2 -i $OUTPUT_SLICE -o $OUTPUT_SLICE -v 0 -s 2
      ;;
    "RescaleImage")
      ${ANTSPATH}/ImageMath 2 $OUTPUT_SLICE RescaleImage $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]}
      ;;
    "TruncateImageIntensity")
      ${ANTSPATH}/ImageMath 2 $OUTPUT_SLICE TruncateImageIntensity $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]} ${PARAMETERS[2]} ${PARAMETERS[3]}
      ;;
    *)
      echo "The operation '$OPERATION' is not an option.  See usage: '$0 -h 1'"
      exit
      ;;
    esac

  done

PERMUTATION_ORDER[0]='2 0 1'
PERMUTATION_ORDER[1]='0 2 1'
PERMUTATION_ORDER[2]='0 1 2'

if [[ $OPERATION != "ExtractAllSlices" ]];
  then
    ${ANTSPATH}/TileImages 3 $TMP_OUTPUT_FILE 1x1x0 ${ALL_OUTPUT_SLICES[@]}
    ${ANTSPATH}/PermuteFlipImageOrientationAxes 3 $TMP_OUTPUT_FILE $TMP_OUTPUT_FILE ${PERMUTATION_ORDER[$WHICH_DIRECTION]} 0 0 0 0
    ${ANTSPATH}/CopyImageHeaderInformation $INPUT_IMAGE $TMP_OUTPUT_FILE $OUTPUT_IMAGE 1 1 1
  fi
