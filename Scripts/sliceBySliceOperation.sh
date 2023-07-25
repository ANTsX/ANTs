#!/bin/bash

VERSION="0.0.0 test"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

PROGRAMS[0]=ExtractSliceFromImage
PROGRAMS[1]=PrintHeader
PROGRAMS[2]=TileImages
PROGRAMS[3]=DenoiseImage
PROGRAMS[4]=ImageMath
PROGRAMS[5]=N4BiasFieldCorrection

for (( i = 0; i < ${#PROGRAMS[@]}; i++ ))
  do
    if ! command -v "${PROGRAMS[$i]}" &> /dev/null
      then
        echo "The ${PROGRAMS[$i]} program can't be found. Please (re)define \$PATH in your environment."
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

if [[ ${WHICH_DIRECTION} -gt 2 || ${WHICH_DIRECTION} -lt 0 ]];
  then
    echo "Error: Direction must be an integer in [0,2 ]"
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

    ExtractSliceFromImage 3 $INPUT_IMAGE $OUTPUT_SLICE $WHICH_DIRECTION $i

    ALL_OUTPUT_SLICES[$i]=$OUTPUT_SLICE

    case "$OPERATION" in
    "ExtractAllSlices")
      mv ${OUTPUT_SLICE} ${OUTPUT_DIR}
      ;;
    "Convolve")
      ImageMath 2 $OUTPUT_SLICE Convolve $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]}
      ;;
    "DenoiseImage")
      DenoiseImage -d 2 -i $OUTPUT_SLICE -o $OUTPUT_SLICE -v 0
      ;;
    "N4BiasFieldCorrection")
      N4BiasFieldCorrection -d 2 -i $OUTPUT_SLICE -o $OUTPUT_SLICE -v 0 -s 2
      ;;
    "RescaleImage")
      ImageMath 2 $OUTPUT_SLICE RescaleImage $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]}
      ;;
    "TruncateImageIntensity")
      ImageMath 2 $OUTPUT_SLICE TruncateImageIntensity $OUTPUT_SLICE ${PARAMETERS[0]} ${PARAMETERS[1]} ${PARAMETERS[2]} ${PARAMETERS[3]}
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
    TileImages 3 $TMP_OUTPUT_FILE 1x1x0 ${ALL_OUTPUT_SLICES[@]}
    PermuteFlipImageOrientationAxes 3 $TMP_OUTPUT_FILE $TMP_OUTPUT_FILE ${PERMUTATION_ORDER[$WHICH_DIRECTION]} 0 0 0 0
    CopyImageHeaderInformation $INPUT_IMAGE $TMP_OUTPUT_FILE $OUTPUT_IMAGE 1 1 1
  fi
