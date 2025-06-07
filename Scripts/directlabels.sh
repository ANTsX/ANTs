#!/bin/bash

# Usage:  $0 imageDimension segmentationImage grayMatterProbabilityImage
# whiteMatterProbabilityImage outputImage

VERSION="0.0"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

if ! command -v ANTS &> /dev/null
then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$PATH in your environment.
  exit
fi

function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -s InputSegmentationImage -l CorticalLabelImage <OPTARGS> -o OutputImage

Example Case:

 echo " bash $0 -d 3 -s brainSegmentation.nii.gz -l corticalLabels.nii.gz -o thickness.nii.gz"

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional single image)
     -s:  InputSegmentationImage:  label image which has label = 2 for gray
            matter regions and label = 3 for white matter regions.
     -l:  CorticalLabelImage: subdivided cortical label image.
     -o:  OutputImage:  thickness output image file.

Optional arguments:

     -c:  Control for parallel computation (default=0) -- 0 == run serially,  1 == SGE qsub,  2 == use PEXEC (localhost)
     -j:  Number of cpu cores to use (default: 2; -- requires "-c 2")
     -g:  Gray matter probability image
     -w:  White matter probability image
     -t:  Thickness prior estimate (could be a .csv file, e.g. for label 3 with a thickness of 4.5 -> 3,4.5)
     -m:  Smoothing sigma
     -r:  Gradient step size (default=0.025) -- smaller in magnitude results in more cautious steps
     -i:  number of iterations (default=50)

--------------------------------------------------------------------------------------
ANTS was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using DiReCT with the following arguments:
      image dimension         = ${DIMENSION}
      segmentation image      = ${SEG_IMAGE}
      cortical label image    = ${LABEL_IMAGE}
        (labels -> ${LABELS[@]})
      output image            = ${OUTPUT_IMAGE}
      wm probability image    = ${WMPROB_IMAGE}
      gm probability image    = ${GMPROB_IMAGE}
      smoothing sigma         = ${SMOOTHING_SIGMA}
      thicknes prior estimate = ${THICKNESS_PRIOR_ESTIMATE}
      do qsub                 = ${DOQSUB}
      number of cores         = ${CORES}
      iteration limit         = ${ITERATION_LIMIT}
      gradient step           = ${GRADIENT_STEP}

PARAMETERS
}

getLabelsAndBoundingBoxes() {

  OUTPUT=(`LabelGeometryMeasures $DIMENSION $LABEL_IMAGE`);

  ## Get labels

  begin=10
  increment=`expr 5 + $DIMENSION + $DIMENSION + $DIMENSION + $DIMENSION`

  count=0;
  for (( i=begin; i<=${#OUTPUT[@]}; i+=$increment )); do
    LABELS[$count]=${OUTPUT[$i]}
    count=`expr $count + 1`
  done

  ## Get bounding boxes
  begin=`expr 15 + $DIMENSION + $DIMENSION`
  increment=`expr 5 + $DIMENSION + $DIMENSION + $DIMENSION + $DIMENSION`

  count=0;
  for (( i=begin; i<=${#OUTPUT[@]}; i+=$increment )); do
    for (( j=0; j<`expr $DIMENSION + $DIMENSION`; j++ )); do
      BOUNDING_BOXES[$count]=${BOUNDING_BOXES[$count]}${OUTPUT[`expr $i+$j`]}
    done
    count=`expr $count + 1`;
  done

  ## read thickness .csv file if it exists
  extension=`echo ${THICKNESS_PRIOR_ESTIMATE#*.}`

  if [[ $extension = 'csv' ]] || [[ $extension = 'txt' ]] ; then

    echo ${LABELS[@]}
    while read line
      do
      bar=(`echo $line | tr ',' ' '`)

      for (( i=0; i<${#LABELS[@]}; i++ )); do
        if [[ ${bar[0]} -eq ${LABELS[$i]} ]]; then
          LABEL_THICKNESSES[$i]=${bar[1]}
        fi
      done
    done < $THICKNESS_PRIOR_ESTIMATE
  else
    for (( count = 0; count<${#LABELS[@]}; count++ )); do
      LABEL_THICKNESSES[$count]=$THICKNESS_PRIOR_ESTIMATE
    done

  fi


}

writeSubimages() {

  if [ ! -d "$TMPDIR" ]; then
    mkdir $TMPDIR
    chmod ugo+rw $TMPDIR
  fi

  OUTFN=${POO%.*.*}

  for (( i=0; i<${#BOUNDING_BOXES[@]}; i++ )); do
    bbox=$(echo ${BOUNDING_BOXES[$i]}|sed 's/,/ /g')
    bbox=$(echo $bbox|sed 's/\[//')
    bbox=$(echo $bbox|sed 's/\]//')

    elements=( $bbox )

    minIndex=""
    maxIndex=""

    count=0
    for (( j=0; j<${#elements[@]}; j+=2 )); do
      minValue=${elements[$j]}

      jp1=`expr $j + 1`
      maxValue=${elements[$jp1]}

      minIndex=${minIndex}x${minValue}
      maxIndex=${maxIndex}x${maxValue}

      count=`expr $count + 1`
    done

   minIndex=${minIndex:1};
   maxIndex=${maxIndex:1};

   grayMatterMask=${TMPDIR}/grayMatterMask.nii.gz
   whiteMatterMask=${TMPDIR}/whiteMatterMask.nii.gz
   gmWmMask=${TMPDIR}/gmWmMask.nii.gz

   # For each label, we perform the following steps:
   #   1. Threshold out everything but the label region
   #   2. Threshold out everything but the white matter
   #   3. Combine the result from 1) and 2) to have an image with only
   #      the white matter and the current label (as the grey matter).
   #
   OUTPUT=(`ThresholdImage $DIMENSION $LABEL_IMAGE $grayMatterMask ${LABELS[$i]} ${LABELS[$i]} 1 0`)
   OUTPUT=(`ThresholdImage $DIMENSION $SEG_IMAGE $whiteMatterMask 3 3 3 0`)
   OUTPUT=(`ImageMath $DIMENSION $grayMatterMask m $grayMatterMask $SEG_IMAGE`)
   OUTPUT=(`ImageMath $DIMENSION $gmWmMask + $grayMatterMask $whiteMatterMask`)
   OUTPUT=(`ExtractRegionFromImage $DIMENSION $gmWmMask ${TMPDIR}seg_${LABELS[$i]}.nii.gz $minIndex $maxIndex`)
   OUTPUT=(`ImageMath $DIMENSION ${TMPDIR}seg_${LABELS[$i]}.nii.gz PadImage ${TMPDIR}seg_${LABELS[$i]}.nii.gz $PADDING`)
   if [ -f $WMPROB_IMAGE ]; then
     OUTPUT=(`ImageMath $DIMENSION $whiteMatterMask m $whiteMatterMask $WMPROB_IMAGE`)
     OUTPUT=(`ExtractRegionFromImage $DIMENSION $whiteMatterMask ${TMPDIR}wm_${LABELS[$i]}.nii.gz $minIndex $maxIndex`);
     OUTPUT=(`ImageMath $DIMENSION ${TMPDIR}wm_${LABELS[$i]}.nii.gz PadImage ${TMPDIR}wm_${LABELS[$i]}.nii.gz $PADDING`);
   fi
   if [ -f $GMPROB_IMAGE ]; then
     OUTPUT=(`ImageMath $DIMENSION $grayMatterMask m $grayMatterMask $GMPROB_IMAGE`)
     OUTPUT=(`ExtractRegionFromImage $DIMENSION $grayMatterMask ${TMPDIR}gm_${LABELS[$i]}.nii.gz $minIndex $maxIndex`);
     OUTPUT=(`ImageMath $DIMENSION ${TMPDIR}gm_${LABELS[$i]}.nii.gz PadImage ${TMPDIR}gm_${LABELS[$i]}.nii.gz $PADDING`);
   fi

   rm -rf $grayMatterMask
   rm -rf $whiteMatterMask
   rm -rf $gmWmMask

  done
}

function jobfnamepadding {

    files=`ls ${TMPDIR}job*.sh`
    BASENAME1=`echo $files[1] | cut -d 'b' -f 1`

    for file in ${files}
      do

      if [ "${#file}" -eq "9" ]
	  then
	  BASENAME2=`echo $file | cut -d 'b' -f 2 `
	  mv "$file" "${BASENAME1}b_000${BASENAME2}"

      elif [ "${#file}" -eq "10" ]
	  then
	  BASENAME2=`echo $file | cut -d 'b' -f 2 `
	  mv "$file" "${BASENAME1}b_00${BASENAME2}"

      elif [ "${#file}" -eq "11" ]
	  then
	  BASENAME2=`echo $file | cut -d 'b' -f 2 `
	  mv "$file" "${BASENAME1}b_0${BASENAME2}"
      fi
   done
}

function pasteImages {

  output=( `CreateImage $DIMENSION $SEG_IMAGE $OUTPUT_IMAGE 0` );

  for (( i=0; i<${#BOUNDING_BOXES[@]}; i++ )); do
    bbox=$(echo ${BOUNDING_BOXES[$i]}|sed 's/,/ /g')
    bbox=$(echo $bbox|sed 's/\[//')
    bbox=$(echo $bbox|sed 's/\]//')

    elements=( $bbox )

    minIndex=""

    count=0
    for (( j=0; j<${#elements[@]}; j+=2 )); do
      minValue=${elements[$j]}
      minValue=`expr $minValue - $PADDING`;

      minIndex=${minIndex}x${minValue}

      count=`expr $count + 1`
    done

    minIndex=${minIndex:1};

    output=(`PasteImageIntoImage $DIMENSION $OUTPUT_IMAGE $TMPDIR/direct_${LABELS[$i]}.nii.gz $OUTPUT_IMAGE $minIndex 0 2 -1`);
  done

}

cleanup() {
  rm -rf ${TMPDIR}
}

################################################################################
#
# Main routine
#
################################################################################

time_start=`date +%s`
CURRENTDIR=`pwd`/
TMPDIR=${CURRENTDIR}/tmp$RANDOM/

DIRECT=KellyKapowski
DIMENSION=3
SEG_IMAGE=""
GMPROB_IMAGE=""
WMPROB_IMAGE=""
OUTPUT_IMAGE=""
LABEL_IMAGE=""
LABELS=()
LABEL_THICKNESSES=()
BOUNDING_BOXES=()
PADDING=5;

GRADIENT_STEP=0.025
SMOOTHING_SIGMA=1.5
THICKNESS_PRIOR_ESTIMATE=8.0
ITERATION_LIMIT=50

DOQSUB=0
CORES=2

# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS

PEXEC=ANTSpexec.sh
SGE=waitForSGEQJobs.pl

for FLE in $PEXEC $SGE
  do
  if [ ! -x $FLE ] ;
      then
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " FILE $FLE DOES NOT EXIST -- OR -- IS NOT EXECUTABLE !!! $0 will terminate."
      echo "--------------------------------------------------------------------------------------"
      echo " if the file is not executable, please change its permissions. "
      exit 1
  fi
done

if [[ $# -eq 0 ]] ; then
  Usage >&2
else
  while getopts "c:d:g:i:j:w:m:o:s:r:t:h:l:" OPT
    do
    case $OPT in
        h) #help
     Usage >&2
     exit 0
     ;;
        c) #use SGE cluster
     DOQSUB=$OPTARG
     if [[ ${#DOQSUB} -gt 2 ]] ; then
         echo " DOQSUB must be an integer value (0=serial, 1=SGE qsub, 2=try pexec ) you passed  -c $DOQSUB "
         exit 1
     fi
     ;;
        d) #dimensions
     DIMENSION=$OPTARG
     if [[ ${DIMENSION} -gt 3 || ${DIMENSION} -lt 2 ]] ; then
         echo " Error:  ImageDimension must be 2 or 3 "
         exit 1
     fi
     ;;
        r) #gradient stepsize (default = 0.025)
     GRADIENT_STEP=$OPTARG
     ;;
        i) #iteration limit (default = 100)
     ITERATION_LIMIT=$OPTARG
     ;;
        j) #number of cpu cores to use (default = 2)
     CORES=$OPTARG
     ;;
        o) #output name prefix
     OUTPUT_IMAGE=$OPTARG
     ;;
        s) #segmentation image
     SEG_IMAGE=$OPTARG
     ;;
        g) #gray matter probability image
     GMPROB_IMAGE=$OPTARG
     ;;
        w) #white matter probability image
     WMPROB_IMAGE=$OPTARG
     ;;
        l) #label image
     LABEL_IMAGE=$OPTARG
     ;;
        m) #smoothing sigma
     SMOOTHING_SIGMA=$OPTARG
     ;;
        t) #thickness prior estimate
     THICKNESS_PRIOR_ESTIMATE=$OPTARG
     ;;
        *) # getopts issues an error message
     Usage >&2
     exit 1
     ;;
    esac
  done
fi

# Get label information and subimages

getLabelsAndBoundingBoxes

echoParameters >&2

writeSubimages

# Job IDs of jobs submitted to queue in loop below
jobIDs=""

# Reinitialize count to 0
count=0

# Submit registration of each input to volume template to SGE or run locally.
for LABEL in ${LABELS[@]}
  do
  # prepare DiReCT command
  segLabelImage=${TMPDIR}seg_${LABEL}.nii.gz
  gmLabelImage=${TMPDIR}gm_${LABEL}.nii.gz
  wmLabelImage=${TMPDIR}wm_${LABEL}.nii.gz

  exe="${DIRECT} -d ${DIMENSION} -c [ $ITERATION_LIMIT,0.000001,10 ] -t ${LABEL_THICKNESSES[$count]} -r $GRADIENT_STEP -m $SMOOTHING_SIGMA -s $segLabelImage -o ${TMPDIR}direct_${LABEL}.nii.gz"
  if [[ -f "$gmLabelImage" ]]; then
    exe=${exe}" -g $gmLabelImage"
  fi
  if [[ -f "$wmLabelImage" ]]; then
    exe=${exe}" -w $wmLabelImage"
  fi
  pexe=" $exe >> ${TMPDIR}job_label_${LABEL}_metriclog.txt "

  # 6 submit to SGE or else run locally
  if [ $DOQSUB -eq 1 ]; then
    id=`qsub -cwd -N DiReCT_${LABEL} -S /bin/bash $QSUBOPTS $exe | awk '{print $3}'`
    jobIDs="$jobIDs $id"
    sleep 0.5
  elif [ $DOQSUB -eq 2 ] ; then
    echo $pexe
    echo $pexe >> ${TMPDIR}job_label_${LABEL}.sh
  elif  [ $DOQSUB -eq 0 ] ; then
    echo "  Performing DiReCT for label ${LABEL}."
    output=(`$exe`)
  fi

  # SGE wait for script to finish
  if [ $DOQSUB -eq 1 ];
    then
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Starting DiReCT on SGE cluster. Label: $LABEL of ${#LABELS[@]} total labels."
    echo "--------------------------------------------------------------------------------------"

    # now wait for the stuff to finish - this will take a while so poll queue every 5 mins
    $SGE 1 300 $jobIDs

    if [ ! $? -eq 0 ]; then
      echo "qsub submission failed - jobs went into error state"
      exit 1;
    fi
  fi
  count=`expr $count + 1`
done

# Run jobs on localhost and wait to finish
if [ $DOQSUB -eq 2 ];
  then
  echo
  echo "--------------------------------------------------------------------------------------"
  echo " Starting DiReCT on max ${CORES} cpucores."
  echo " Progress can be viewed in ${TMPDIR}job_label*_metriclog.txt"
  echo "--------------------------------------------------------------------------------------"
  jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
  chmod +x ${TMPDIR}job*.sh
  $PEXEC -j ${CORES} sh ${TMPDIR}job*.sh
fi

pasteImages

cleanup

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with DiReCT"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
