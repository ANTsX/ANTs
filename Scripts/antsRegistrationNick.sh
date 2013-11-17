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

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

Example:

`basename $0` -d 3 -f fixedImage.nii.gz -m movingImage.nii.gz -o output

--------------------------------------------------------------------------------------
ANTS was created by:
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

     -s:  spline distance for deformable B-spline SyN transform (default = 26)

--------------------------------------------------------------------------------------
Get the latest ANTS version at:
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

Please reference http://www.ncbi.nlm.nih.gov/pubmed/20851191 when employing this script
in your studies. A reproducible evaluation of ANTs similarity metric performance in
brain image registration:

* Avants BB, Tustison NJ, Song G, Cook PA, Klein A, Gee JC. Neuroimage, 2011.

Also see http://www.ncbi.nlm.nih.gov/pubmed/19818860 for more details.
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
 Fixed image:              $FIXEDIMAGE
 Moving image:             $MOVINGIMAGE
 Number of threads:        $NUMBEROFTHREADS
 Spline distance:          $SPLINEDISTANCE
 Transform type:           $TRANSFORMTYPE
--------------------------------------------------------------------------------------
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
FIXEDIMAGE=''
MOVINGIMAGE=''
OUTPUTNAME=output
NUMBEROFTHREADS=1
SPLINEDISTANCE=26
TRANSFORMTYPE='s'

# reading command line arguments
while getopts "d:f:h:m:n:o:s:t:" OPT
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
   FIXEDIMAGE=$OPTARG
   ;;
      m)  # moving image
   MOVINGIMAGE=$OPTARG
   ;;
      n)  # number of threads
   NUMBEROFTHREADS=$OPTARG
   ;;
      o) #output name prefix
   OUTPUTNAME=$OPTARG
   ;;
      s)  # spline distance
   SPLINEDISTANCE=$OPTARG
   ;;
      t)  # transform type
   TRANSFORMTYPE=$OPTARG
   ;;
  esac
done

###############################
#
# Check inputs
#
###############################

if [[ ! -f "$FIXEDIMAGE" ]];
  then
    echo "Fixed image '$FIXEDIMAGE' does not exist.  See usage: '$0 -h 1'"
    exit
  fi
if [[ ! -f "$MOVINGIMAGE" ]];
  then
    echo "Moving image '$MOVINGIMAGE' does not exist.  See usage: '$0 -h 1'"
    exit
  fi

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
# Construct mapping stages
#
##############################

RIGIDSTAGE="--initial-moving-transform [$FIXEDIMAGE,$MOVINGIMAGE,1] \
            --transform Rigid[0.1] \
            --metric MI[$FIXEDIMAGE,$MOVINGIMAGE,1,32,Regular,0.25] \
            --convergence 1000x500x250x100 \
            --shrink-factors 8x4x2x1 \
            --smoothing-sigmas 3x2x1x0"

AFFINESTAGE="--transform Affine[0.1] \
             --metric MI[$FIXEDIMAGE,$MOVINGIMAGE,1,32,Regular,0.25] \
             --convergence 1000x500x250x100 \
             --shrink-factors 8x4x2x1 \
             --smoothing-sigmas 3x2x1x0"

SYNSTAGE="--metric CC[$FIXEDIMAGE,$MOVINGIMAGE,1,4] \
          --convergence 100x70x50x20 \
          --shrink-factors 6x4x2x1 \
          --smoothing-sigmas 3x2x1x0"

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

${ANTS} --dimensionality $DIM \
								--output [$OUTPUTNAME,${OUTPUTNAME}Warped.nii.gz] \
								--interpolation Linear \
								--winsorize-image-intensities [0.005,0.995] \
								$STAGES

###############################
#
# Restore original number of threads
#
###############################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ORIGINALNUMBEROFTHREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
