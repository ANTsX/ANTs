#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

shopt -s extglob

VERSION="0.0.0"

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTS=antsRegistration
WARP=antsApplyTransforms
N4=N4BiasFieldCorrection
PEXEC=ANTSpexec.sh
SGE=waitForSGEQJobs.pl
PBS=waitForPBSQJobs.pl
XGRID=waitForXGridJobs.pl
SLURM=waitForSlurmJobs.pl

fle_error=0
for FLE in $ANTS $WARP $N4 $PEXEC $SGE $XGRID $PBS $SLURM
  do
    if ! command -v $FLE &> /dev/null
      then
        echo
        echo "-----------------------------------------------------------------------------"
        echo " FILE $FLE DOES NOT EXIST -- OR -- IS NOT EXECUTABLE !!! $0 will terminate."
        echo "-----------------------------------------------------------------------------"
        echo " if the file is not executable, please change its permissions. "
        fle_error=1
      fi
  done

if [[ $fle_error = 1 ]];
  then
    echo "missing helper script"
    exit 1
  fi

function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -o OutputPrefix <other options> <images>

Compulsory arguments (minimal command line requires SGE/PBS cluster, otherwise use -c and
-j options):

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
          ImageDimension: 4 (for template generation of time-series data)

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.  Optionally, one can specify a .csv or .txt file where each
          line is the location of the input image.  One can also specify more than
          one file for each image for multi-modal template construction (e.g. t1 and t2).
          For the multi-modal case, the templates will be consecutively numbered (e.g.
          {OutputPrefix}template0.nii.gz, {OutputPrefix}template1.nii.gz, ...).

NB: All images to be added to the template should be in the same directory, and this
    script should be invoked from that directory.

Optional arguments:

     -a   image statistic used to summarize images (default 1)
          0 = mean
          1 = mean of normalized intensities
          2 = median

          Normalization here means dividing each image by its mean intensity.

     -A   sharpening applied to template at each iteration (default 1)
          0 = none
          1 = Laplacian
          2 = Unsharp mask

     -b:  Backup images and results from all iterations (default = 0):  Boolean to save
          the transform files, bias corrected inputs, templates, transforms, and warped images
          for each iteration. By default, only the templates and the shape update warp field
          are saved.

     -c:  Control for parallel computation (default 0):
          0 = run serially
          1 = SGE qsub
          2 = use PEXEC (localhost)
          3 = Apple XGrid
          4 = PBS qsub
          5 = SLURM

     -e   use single precision ( default 1 )

     -g:  Gradient step size (default 0.25): smaller in magnitude results in more
          cautious steps.  Use smaller steps to refine template details.
          0.25 is an upper (aggressive) limit for this parameter.

     -i:  Iteration limit (default 4): iterations of the template construction
          (Iteration limit)*NumImages registrations.

     -j:  Number of cpu cores to use locally for pexec option (default 2; requires "-c 2")

     -k:  Number of modalities used to construct the template (default 1):  For example,
          if one wanted to create a multimodal template consisting of T1,T2,and FA
          components ("-k 3").

     -w:  Modality weights used in the similarity metric (default = 1): specified as
          e.g. 1x0.5x0.75.

     -q:  Max iterations for each pairwise registration (default = 100x100x70x20): specified in
          the form JxK...xF where
            J = max iterations at first (coarsest) resolution
            K = max iterations at next resolution
            F = max iterations at the final resolution
          Finer resolutions take much more time per iteration than coarser resolutions.
          The resolution for each level is controlled by the shrink factors "-f", so
          "-q JxKxL -f 4x2x1" does J iterations at factor 4, K at factor 2, L at factor 1.

     -f:  Shrink factors in pairwise registration (default = 6x4x2x1): in the same form as "-q"
          max iterations. Must have the same number of components as the iterations "-q" and smoothing
          "-s". The shrink factors are integer factors for downsampling the the virtual space
          (usually the template image) during registration.

     -s:  Smoothing kernels in pairwise registration (default = 3x2x1x0): also in the same form as
          "-q" and "-f", with the same number of components. Standard deviation of a Gaussian smoothing
          kernel applied to the images before downsampling at each level. Needs to have the same number
          of components as the number of iterations and shrink factors. The kernel may be specified in
          mm units or voxels with "AxBxCmm" or "AxBxCvox". Missing units implies vox.

     -n:  N4BiasFieldCorrection of moving image: 0 == off, 1 == on (default 1).

     -o:  OutputPrefix; A prefix that is prepended to all output files (default = "antsBTP").

     -p:  Commands to prepend to job scripts (e.g., change into appropriate directory, set
          paths, etc)

     -r:  Do rigid-body registration of inputs to the initial template, before doing the main
          pairwise registration. 0 == off 1 == on (default 0). If you are trying to refine or update
          an existing template, you would use '-r 0'.
          Rigid initialization is useful when you do not have an initial template, or you want to use
          a single image as a reference for rigid alignment only. For example,
            "-z tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz -y 0 -r 1"
          will rigidly align the inputs to the MNI template, and then use their average to begin the
          template building process.

     -l:  Use linear image registration stages during the pairwise (template/subject)
          deformable registration.  Otherwise, registration is limited to SyN or
          B-spline SyN (see '-t' option).  This is '1' by default.

     -m:  Type of similarity metric used for pairwise registration (default = CC). Options are case
          sensitive.
            CC = cross-correlation
            MI = mutual information
            MSQ = mean square difference
            DEMONS = demon's metric
          A similarity metric per modality can be specified.  If the CC metric is chosen,
          one can also specify the radius in brackets, e.g. '-m "CC[4]"'. This option controls
          the metric for the transformation specified with "-t", so "-m CC -t SyN" means the SyN
          stage is run with the CC metric; preceding linear stages use MI.

     -t:  Type of transformation model used for registration (default = SyN):  Options are case
          sensitive.
            SyN = Greedy SyN
            BSplineSyN = Greedy B-spline SyN
            TimeVaryingVelocityField = Time-varying velocity field
            TimeVaryingBSplineVelocityField = Time-varying B-spline velocity field

          Transformation parameters may be specified with brackets, eg '-t "SyN[0.1,3,0]"'.

          The transformations above are used after Rigid + Affine linear registration, unless linear
          registration is disabled with "-l". To use linear registration only, options are:
            Affine = Affine (runs Rigid + Affine)
            Rigid = Rigid (runs Rigid only).

     -u:  Walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Memory limit (default = 8gb):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -x:  XGrid arguments (e.g., -x "-p password -h controlhost")

     -y:  Update the template with the full affine transform (default 1). If 0, the rigid
          component of the affine transform will not be used to update the template. If your
          template drifts in translation or orientation try "-y 0".

     -z:  Use this this volume as the target of all inputs. When not used, the script will create an unbiased
          starting point by averaging all inputs, then aligning the center of mass of all inputs to that of
          the initial average. If you do not use -z, it is recommended to use "-r 1". Use the full path.
          For multiple modalities, specify -z modality1.nii.gz -z modality2.nii.gz ...
          in the same modality order as the input images.

Example:

`basename $0` -d 3 -i 3 -k 1 -f 4x2x1 -s 2x1x0vox -q 30x20x4 -t SyN -m CC -c 0 -r 1 -o MY sub*avg.nii.gz

Multimodal example:

`basename $0` -d 3 -i 3 -k 2 -f 4x2x1 -s 2x1x0vox -q 30x20x4 -t SyN -z t1.nii.gz -z t2.nii.gz \\
 -m CC -c 0 -o MY templateInput.csv

where templateInput.csv contains

subjectA_t1.nii.gz,subjectA_t2.nii.gz
subjectB_t1.nii.gz,subjectB_t2.nii.gz
...

Example of first building an affine template, then refining it with SyN:

`basename $0` -d 3 -i 2 -A 2 -k 1 -f 6x4x2x1 -s 4x2x1x0vox -q 200x100x50x0 -t Affine -m MI -c 0 \\
  -r 1 -o MY_Affine_ sub*avg.nii.gz

`basename $0` -d 3 -i 5 -k 1 -f 6x4x2x1 -s 4x2x1x0vox -q 50x100x70x20 -t SyN -m CC -c 0 \\
  -r 0 -o MY_Deformable_ -z MY_Affine_template0.nii.gz sub*avg.nii.gz


Output

{OutputPrefix}template{m}.nii.gz
  final template for each modality m.

{OutputPrefix}template{m}{inputFile}{n}WarpedToTemplate.nii.gz
{OutputPrefix}template{m}{inputFile}{n}0GenericAffine.mat
{OutputPrefix}template{m}{inputFile}{n}1Warp.nii.gz
{OutputPrefix}template{m}{inputFile}{n}1InverseWarp.nii.gz
  each of n input images warped to the penultimate template m, with transforms. If the template has converged,
  these should be well aligned to {OutputPrefix}template{m}.nii.gz.

intermediateTemplates/
                     initial_{OutputPrefix}template{m}.nii.gz :
                       initial template
                     initialRigid_{OutputPrefix}template{m}.nii.gz :
                       initial rigid template if requested with "-r 1"
                     {transform}_iteration{i}_{OutputPrefix}template{m}.nii.gz
                       Template computed with {transform} (-t) for each iteration (-i) and modality.
                     {transform}_iteration{i}_shapeUpdateWarp.nii.gz
                       Shape update warp applied to the template at iteration i. As the template converges,
                       the magnitude of the update warp will converge to a minimal value.

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

The script has been updated and improved since this publication.

--------------------------------------------------------------------------------------
Script by Nick Tustison
--------------------------------------------------------------------------------------
Apple XGrid support by Craig Stark
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

function reportMappingParameters {
    cat <<REPORTMAPPINGPARAMETERS

--------------------------------------------------------------------------------------
 Mapping parameters
--------------------------------------------------------------------------------------
 Dimensionality:           $DIM
 Do N4 bias correction:    $N4CORRECT
 Back up each iteration:   $BACKUPEACHITERATION
 Similarity metric:        ${METRICTYPE[@]}
 Gradient step:            $GRADIENTSTEP
 Transformation:           $TRANSFORMATIONTYPE
 Max iterations:           $MAXITERATIONS
 Smoothing factors:        $SMOOTHINGFACTORS
 Shrink factors:           $SHRINKFACTORS
 Output prefix:            $OUTPUTNAME
 Template:                 $TEMPLATENAME
 Template update steps:    $ITERATIONLIMIT
 Template population:      $IMAGESETVARIABLE
 Number of modalities:     $NUMBEROFMODALITIES
 Modality weights:         $MODALITYWEIGHTSTRING
 Image statistic:          $STATSMETHOD
 Sharpening method:        $SHARPENMETHOD
 Shape update full affine: $AFFINE_UPDATE_FULL
--------------------------------------------------------------------------------------
REPORTMAPPINGPARAMETERS
}

function summarizeimageset() {

  local dim=$1
  shift
  local output=$1
  shift
  local summarizemethod=$1
  shift
  local sharpenmethod=$1
  shift
  local images=( "${@}" )

  if [[ ${#images[@]} -ne ${IMAGESPERMODALITY} ]]
    then
      echo "ERROR summarizeimageset - imagelist length is ${#images[@]}, expected ${IMAGESPERMODALITY}"
      exit 1
    fi

  rm -f "$output"

  case $summarizemethod in
    0) #mean
      AverageImages $dim $output 0 ${images[@]}
      ;;
    1) #mean of normalized images
      AverageImages $dim $output 2 ${images[@]}
      ;;
    2) #median
      local image
      for image in "${images[@]}";
        do
          echo $image >> ${output}_list.txt
        done
      ImageSetStatistics $dim ${output}_list.txt ${output} 0
      rm ${output}_list.txt
      ;;
  esac

  if [[ ! -f "$output" ]];
    then
      echo "summarizeimageset: ERROR - output file $output could not be created"
      exit 1
    fi

  case $sharpenmethod in
    0)
      echo "Sharpening method none"
      ;;
    1)
      echo "Laplacian sharpening"
      ImageMath $dim $output Sharpen $output 0
      ;;
    2)
      echo "Unsharp mask sharpening"
      ImageMath $dim $output UnsharpMask $output 0.5 1 0 0
      ;;
  esac

  local sharpenExit=$?

  if [[ $? -ne 0 ]]
    then
      echo "summarizeimageset: ERROR - template sharpening failed with status $?"
      exit 1
    fi

}

function shapeupdatetotemplate() {

   echo "shapeupdatetotemplate()"

    # local declaration of values
    dim=$1
    template=$2
    templatename=$3
    outputname=$4
    gradientstep=-$5
    whichtemplate=$6
    statsmethod=$7
    sharpenmethod=$8

# debug only
# echo $dim
# echo ${template}
# echo ${templatename}
# echo ${outputname}
# echo ${outputname}*WarpedToTemplate.nii*
# echo ${gradientstep}

# We find the average warp to the template and apply its inverse to the template image
# This keeps the template shape stable over multiple iterations of template building

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate---voxel-wise averaging of the warped images to the current template"
    date
    #echo "   AverageImages $dim ${template} 1 ${templatename}${whichtemplate}*WarpedToTemplate.nii.gz    "
    #echo "    ImageSetStatistics $dim ${whichtemplate}WarpedToTemplateList.txt ${template} 0"
    echo "--------------------------------------------------------------------------------------"

    imagelist=(`ls ${outputname}*-modality${whichtemplate}-*-WarpedToTemplate.nii.gz`)
    if [[ ${#imagelist[@]} -ne ${IMAGESPERMODALITY} ]]
      then
        echo "ERROR shapeupdatetotemplate - imagelist length is ${#imagelist[@]}, expected ${IMAGESPERMODALITY}"
        exit 1
      fi

    summarizeimageset $dim $template $statsmethod $sharpenmethod ${imagelist[@]}

    WARPLIST=( `ls ${outputname}input*-[0-9]Warp.nii.gz 2> /dev/null` ) || true
    NWARPS=${#WARPLIST[@]}
    echo "number of warps = $NWARPS"
    if [[ $NWARPS -ne 0 ]];
      then
        echo "${WARPLIST[@]}"
      fi

    if [[ $whichtemplate -eq 0 ]];
      then

        if [[ $NWARPS -ne 0 ]]; then
          echo "$NWARPS does not equal 0"
          echo
          echo "--------------------------------------------------------------------------------------"
          echo " shapeupdatetotemplate---voxel-wise averaging of the inverse warp fields (from subject to template)"
          echo "   AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 ${WARPLIST[@]}"
          date
          echo "--------------------------------------------------------------------------------------"
          AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 ${WARPLIST[@]}

          echo
          echo "--------------------------------------------------------------------------------------"
          echo " shapeupdatetotemplate---scale the averaged inverse warp field by the gradient step"
          echo "   MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz"
          date
          echo "--------------------------------------------------------------------------------------"
          MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz
        fi

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---average the affine transforms (template <-> subject)"
        echo "                      ---transform the inverse field by the resulting average affine transform"
        echo "   ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0GenericAffine.mat ${outputname}*-input*GenericAffine.mat"
        echo "   ${WARP} -d ${dim} -e vector -i ${templatename}0warp.nii.gz -o ${templatename}0warp.nii.gz -t [ ${templatename}0GenericAffine.mat,1 ] -r ${template} --verbose 1"
        echo "--------------------------------------------------------------------------------------"

        ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0GenericAffine.mat ${outputname}input*GenericAffine.mat

        if [[ $NWARPS -ne 0 ]];
          then
            ${WARP} -d ${dim} -e vector -i ${templatename}0warp.nii.gz -o ${templatename}0warp.nii.gz -t [ ${templatename}0GenericAffine.mat,1 ] -r ${template} --verbose 1
            MeasureMinMaxMean ${dim} ${templatename}0warp.nii.gz ${templatename}warplog.txt 1
          fi
      fi


    if [[ -f "${templatename}0warp.nii.gz" ]];
      then
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---warp each template by the resulting transforms"
        echo "   ${WARP} -d ${dim} --float $USEFLOAT --verbose 1 -i ${template} -o ${template} -t [ ${templatename}0GenericAffine.mat,1 ] -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -r ${template}"
        echo "--------------------------------------------------------------------------------------"
        ${WARP} -d ${dim} --float $USEFLOAT --verbose 1 -i ${template} -o ${template} -t [ ${templatename}0GenericAffine.mat,1 ] -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -r ${template}
      else
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---warp each template by the resulting transform"
        echo "   ${WARP} -d ${dim} --float $USEFLOAT --verbose 1 -i ${template} -o ${template} -t [ ${templatename}0GenericAffine.mat,1 ] -r ${template}"
        echo "--------------------------------------------------------------------------------------"
        ${WARP} -d ${dim} --float $USEFLOAT --verbose 1 -i ${template} -o ${template} -t [ ${templatename}0GenericAffine.mat,1 ] -r ${template}
      fi

}

function jobfnamepadding {

    outdir=`dirname ${TEMPLATES[0]}`
    if [[ ${#outdir} -eq 0 ]]
        then
        outdir=`pwd`
    fi

    files=`ls ${outdir}/job*.sh`
    BASENAME1=`echo $files[1] | cut -d 'b' -f 1`

    for file in ${files}
      do

      if [[ "${#file}" -eq "9" ]];
       then
         BASENAME2=`echo $file | cut -d 'b' -f 2 `
         mv "$file" "${BASENAME1}b_000${BASENAME2}"

      elif [[ "${#file}" -eq "10" ]];
        then
          BASENAME2=`echo $file | cut -d 'b' -f 2 `
          mv "$file" "${BASENAME1}b_00${BASENAME2}"

      elif [[ "${#file}" -eq "11" ]];
        then
          BASENAME2=`echo $file | cut -d 'b' -f 2 `
          mv "$file" "${BASENAME1}b_0${BASENAME2}"
      fi
    done
}

function setCurrentImageSet() {

WHICHMODALITY=$1

CURRENTIMAGESET=()
COUNT=0

for (( g = $WHICHMODALITY; g < ${#IMAGESETARRAY[@]}; g+=$NUMBEROFMODALITIES ))
  do
    CURRENTIMAGESET[$COUNT]=${IMAGESETARRAY[$g]}
    COUNT=$(( COUNT + 1 ))
  done
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
currentdir=`pwd`
nargs=$#

STATSMETHOD=1
SHARPENMETHOD=1
USEFLOAT=1
BACKUPEACHITERATION=0
MAXITERATIONS=100x100x70x20
SMOOTHINGFACTORS=3x2x1x0
SHRINKFACTORS=6x4x2x1
METRICTYPE=()
TRANSFORMATIONTYPE="SyN"
NUMBEROFMODALITIES=1
MODALITYWEIGHTSTRING=""
N4CORRECT=1
DOLINEAR=1
NOWARP=0
DOQSUB=0
GRADIENTSTEP=0.25
ITERATIONLIMIT=4
CORES=2
TDIM=0
RIGID=0
range=0
REGTEMPLATES=()
TEMPLATES=()
CURRENTIMAGESET=()
XGRIDOPTS=""
SCRIPTPREPEND=""
WALLTIME="20:00:00"
MEMORY="8gb"
# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS
OUTPUTNAME=antsBTP
TEMPLATENAME=${OUTPUTNAME}template
AFFINE_UPDATE_FULL=1

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
if [[ $# -eq 0 || "$1" == "-h" ]];
  then
    Usage >&2
  fi

# reading command line arguments
while getopts "A:a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:s:r:t:u:v:w:x:y:z:" OPT
  do
  case $OPT in
      h) #help
      Usage >&2
      exit 0
   ;;
      A) # Sharpening method
      SHARPENMETHOD=$OPTARG
   ;;
      a) # summarizing statistic
      STATSMETHOD=$OPTARG
   ;;
      b) #backup each iteration
   BACKUPEACHITERATION=$OPTARG
   ;;
      e) #float boolean
   USEFLOAT=$OPTARG
   ;;
      c) #use SGE cluster
   DOQSUB=$OPTARG
   if [[ $DOQSUB -gt 5 ]];
     then
       echo " DOQSUB must be an integer value (0=serial, 1=SGE qsub, 2=try pexec, 3=XGrid, 4=PBS qsub, 5=SLURM) you passed  -c $DOQSUB "
       exit 1
     fi
   ;;
      d) #dimensions
   DIM=$OPTARG
   if [[ ${DIM} -eq 4 ]];
     then
       DIM=3
       TDIM=4
     fi
   ;;
      g) #gradient stepsize (default = 0.25)
   GRADIENTSTEP=$OPTARG
   ;;
      i) #iteration limit (default = 3)
   ITERATIONLIMIT=$OPTARG
   ;;
      j) #number of cpu cores to use (default = 2)
   CORES=$OPTARG
   ;;
      k) #number of modalities used to construct the template (default = 1)
   NUMBEROFMODALITIES=$OPTARG
   ;;
      l) #do linear (rigid + affine) for deformable registration
   DOLINEAR=$OPTARG
   ;;
      w) #modality weights (default = 1)
   MODALITYWEIGHTSTRING=$OPTARG
   ;;
      q) #max iterations other than default
   MAXITERATIONS=$OPTARG
   ;;
      f) #shrink factors
   SHRINKFACTORS=$OPTARG
   ;;
      s) #smoothing factors
   SMOOTHINGFACTORS=$OPTARG
   ;;
      n) #apply bias field correction
   N4CORRECT=$OPTARG
   ;;
      o) #output name prefix
   OUTPUTNAME=$OPTARG
   TEMPLATENAME=${OUTPUTNAME}template
   ;;
      p) #Script prepend
   SCRIPTPREPEND=$OPTARG
   ;;
      m) #similarity model
   METRICTYPE[${#METRICTYPE[@]}]=$OPTARG
   ;;
      r) #start with rigid-body registration
   RIGID=$OPTARG
   ;;
      t) #transformation model
   TRANSFORMATIONTYPE=$OPTARG
   ;;
      u)
   WALLTIME=$OPTARG
   ;;
      v)
   MEMORY=$OPTARG
   ;;
      x) #initialization template
   XGRIDOPTS=$OPTARG
   ;;
      y) # update with full affine, 0 for no rigid (default = 1)
   AFFINE_UPDATE_FULL=$OPTARG
   ;;
      z) #initialization template
   REGTEMPLATES[${#REGTEMPLATES[@]}]=$OPTARG
   ;;
      \?) # getopts issues an error message
      echo "$USAGE" >&2
      exit 1
      ;;
  esac
done

# Provide different output for Usage and Help
if [[ ${TDIM} -eq 4 && $nargs -lt 5 ]];
  then
    Usage >&2
elif [[ ${TDIM} -eq 4 && $nargs -eq 5 ]];
  then
    echo ""
    # This option is required to run 4D template creation on SGE with a minimal command line
elif [[ $nargs -lt 6 ]]
  then
    Usage >&2
fi

if [[ -z ${DIM} ]]
  then
    echo "Image dimension (-d) is required"
    exit 1
  fi

if [[ ${OUTPUTNAME} == */ ]];
  then
    OUTPUT_DIR=${OUTPUTNAME%/}
  else
    OUTPUT_DIR=$(dirname $OUTPUTNAME)
  fi

if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi

# Intermediate template output. Keep the template for each iteration and also the average warp if defined.
# Useful for debugging and monitoring convergence
intermediateTemplateDir=${OUTPUT_DIR}/intermediateTemplates
mkdir -p $intermediateTemplateDir

if [[ $DOQSUB -eq 1 || $DOQSUB -eq 4 ]];
  then
    qq=`which  qsub`
    if [[ ${#qq} -lt 1 ]];
      then
        echo "do you have qsub?  if not, then choose another c option ... if so, then check where the qsub alias points ..."
        exit
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

for (( i = 0; i < $NUMBEROFMODALITIES; i++ ))
  do
    TEMPLATES[$i]=${TEMPLATENAME}${i}.nii.gz
  done

if [[ ${#METRICTYPE[@]} -eq 0 ]];
  then
    METRICTYPE[0]=CC
  fi

if [[ ${#METRICTYPE[@]} -eq 1 ]];
  then
    for (( i = 1; i < $NUMBEROFMODALITIES; i++ ))
      do
        METRICTYPE[${#METRICTYPE[@]}]=${METRICTYPE[0]}
      done
  fi

if [[ ${#METRICTYPE[@]} -ne $NUMBEROFMODALITIES ]];
  then
    echo "The number of similarity metrics does not match the number of specified modalities (see -s option)"
    exit
  fi

if [[ ! -n "$MODALITYWEIGHTSTRING" ]];
  then
    for (( i = 0; i < $NUMBEROFMODALITIES; i++ ))
      do
        MODALITYWEIGHTS[$i]=1
      done
  else
    MODALITYWEIGHTS=(`echo $MODALITYWEIGHTSTRING | tr 'x' "\n"`)
    if [[ ${#MODALITYWEIGHTS[@]} -ne $NUMBEROFMODALITIES ]];
      then
        echo "The number of weights (specified e.g. -w 1x1x1) does not match the number of specified modalities (see -k option)";
        exit 1
      fi
  fi


TRANSFORMATION=''
# Allow SyN[step] or SyN[ step ]
# Users should use the latter to avoid glob problems, but but don't want force that
#
if [[ $TRANSFORMATIONTYPE == "BSplineSyN"* ]];
  then
    if [[ $TRANSFORMATIONTYPE == "BSplineSyN["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="BSplineSyN"
      else
        TRANSFORMATION="BSplineSyN[ 0.2,26,0,3 ]"
    fi
elif [[ $TRANSFORMATIONTYPE == "SyN"* ]];
  then
    if [[ $TRANSFORMATIONTYPE == "SyN["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="SyN"
      else
        TRANSFORMATION="SyN[ 0.2,3,0 ]"
    fi
elif [[ $TRANSFORMATIONTYPE == "TimeVaryingVelocityField"* ]];
  then
    if [[ $TRANSFORMATIONTYPE == "TimeVaryingVelocityField["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="TimeVaryingVelocityField"
      else
        TRANSFORMATION="TimeVaryingVelocityField[ 0.5,4,3,0,0,0 ]"
    fi
elif [[ $TRANSFORMATIONTYPE == "TimeVaryingBSplineVelocityField"* ]];
  then
    if [[ $TRANSFORMATIONTYPE == "TimeVaryingBSplineVelocityField["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="TimeVaryingBSplineVelocityField"
      else
        TRANSFORMATION="TimeVaryingVelocityField[ 0.5,12x12x12x2,4,3 ]"
    fi
elif [[ $TRANSFORMATIONTYPE == "Affine"* ]];
  then
    echo "Linear transforms only!!!!"
    NOWARP=1
    if [[ $TRANSFORMATIONTYPE == "Affine["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="Affine"
      else
        TRANSFORMATION="Affine[ 0.1 ]"
    fi
elif [[ $TRANSFORMATIONTYPE == "Rigid"* ]];
  then
    echo "Rigid transforms only!!!!"
    NOWARP=1
    if [[ $TRANSFORMATIONTYPE == "Rigid["*"]" ]]
      then
        TRANSFORMATION=${TRANSFORMATIONTYPE}
        TRANSFORMATIONTYPE="Rigid"
      else
        TRANSFORMATION="Rigid[ 0.1 ]"
    fi
else
  echo "Invalid transformation. See `basename $0` -h for help menu."
  exit 1
fi

if [[ $STATSMETHOD -lt 0 ]] || [[ $STATSMETHOD -gt 2 ]];
  then
  echo "Invalid stats type: using normalized mean (1)"
  STATSMETHOD=1
fi

if [[ $SHARPENMETHOD -lt 0 ]] || [[ $SHARPENMETHOD -gt 2 ]];
  then
  echo "Invalid sharpening method: using Laplacian (1)"
  SHARPENMETHOD=1
fi

AVERAGE_AFFINE_PROGRAM="AverageAffineTransform"

if [[ $AFFINE_UPDATE_FULL -eq 0 ]];
  then
    AVERAGE_AFFINE_PROGRAM="AverageAffineTransformNoRigid"
  fi

# Creating the file list of images to make a template from.
# Shiftsize is calculated because a variable amount of arguments can be used on the command line.
# The shiftsize variable will give the correct number of arguments to skip. Issuing shift $shiftsize will
# result in skipping that number of arguments on the command line, so that only the input images remain.
shiftsize=$(($OPTIND - 1))
shift $shiftsize
# The invocation of $* will now read all remaining arguments into the variable IMAGESETVARIABLE
IMAGESETVARIABLE=$*
NINFILES=$(($nargs - $shiftsize))
IMAGESETARRAY=()



if [[ ${NINFILES} -eq 0 ]];
    then
    echo "Please provide at least 2 filenames for the template inputs."
    echo "Use `basename $0` -h for help"
    exit 1
elif [[ ${NINFILES} -eq 1 ]];
    then
    extension=`echo ${IMAGESETVARIABLE##*.}`
    if [[ $extension = 'csv' || $extension = 'txt' ]];
        then
        IMAGESFILE=$IMAGESETVARIABLE
        IMAGECOUNT=0
        while read line
            do
            line=$(echo "$line" | tr -d '\r') # remove carriage return from python / windows line-endings
            files=(`echo $line | tr "," "\n"`)
            if [[ ${#files[@]} -ne $NUMBEROFMODALITIES ]];
                then
                echo "The number of files in the csv file does not match the specified number of modalities."
                echo "See the -k option."
                exit 1
            fi
            for (( i = 0; i < ${#files[@]}; i++ ));
                do
                IMAGESETARRAY[$IMAGECOUNT]=${files[$i]}
                IMAGECOUNT=$(( IMAGECOUNT + 1 ))
            done
         done < $IMAGESFILE
    else
        range=`ImageMath $TDIM abs nvols ${IMAGESETVARIABLE} | tail -1 | cut -d "," -f 4 | cut -d " " -f 2 | cut -d "]" -f 1 `
        if [[ ${range} -eq 1 && ${TDIM} -ne 4 ]];
          then
            echo "Please provide at least 2 filenames for the template."
            echo "Use `basename $0` -h for help"
            exit 1
        elif [[ ${range} -gt 1 && ${TDIM} -ne 4 ]]
          then
            echo "This is a multivolume file. Use -d 4"
            echo "Use `basename $0` -h for help"
            exit 1
        elif [[ ${range} -gt 1 && ${TDIM} -eq 4 ]];
          then
            echo
            echo "--------------------------------------------------------------------------------------"
            echo " Creating template of 4D input. "
            echo "--------------------------------------------------------------------------------------"

             #splitting volume
             #setting up working dirs
             tmpdir=${currentdir}/tmp_${RANDOM}_${RANDOM}_${RANDOM}_$$
             (umask 077 && mkdir ${tmpdir}) || {
                 echo "Could not create temporary directory! Exiting." 1>&2
                 exit 1
                 }

             mkdir ${tmpdir}/selection

             #split the 4D file into 3D elements
             cp ${IMAGESETVARIABLE} ${tmpdir}/
             cd ${tmpdir}/
             # ImageMath $TDIM vol0.nii.gz TimeSeriesSubset ${IMAGESETVARIABLE} ${range}
             # rm -f ${IMAGESETVARIABLE}

             # selecting 16 volumes randomly from the timeseries for averaging, placing them in tmp/selection folder.
             # the script will automatically divide timeseries into $total_volumes/16 bins from wich to take the random volumes;
             # if there are more than 32 volumes in the time-series (in case they are smaller

             nfmribins=16
            if [[ ${range} -gt 31 ]];
              then
                BINSIZE=$((${range} / ${nfmribins}))
                j=1 # initialize counter j
                for ((i = 0; i < ${nfmribins}; i++))
                    do
                    FLOOR=$((${i} * ${BINSIZE}))
                    BINrange=$((${j} * ${BINSIZE}))
                    # Retrieve random number between two limits.
                    number=0   #initialize
                    while [[ "$number" -le $FLOOR ]];
                        do
                        number=$RANDOM
                        if [[ $i -lt 15 ]];
                            then
                            number=$((number %= $BINrange))  # Scales $number down within $range.
                        elif [[ $i -eq 15 ]];
                            then
                            number=$((number %= $range))  # Scales $number down within $range.
                        fi
                    done
                    #debug only
                    echo
                    echo "Random number between $FLOOR and $BINrange ---  $number"
                    #  echo "Random number between $FLOOR and $range ---  $number"

                    if [[ ${number} -lt 10 ]];
                        then
                        ImageMath $TDIM selection/vol000${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #   cp vol000${number}.nii.gz selection/
                    elif [[ ${number} -ge 10 && ${number} -lt 100 ]];
                        then
                        ImageMath $TDIM selection/vol00${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #   cp vol00${number}.nii.gz selection/
                    elif [[ ${number} -ge 100 && ${number} -lt 1000 ]];
                        then
                        ImageMath $TDIM selection/vol0${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #   cp vol0${number}.nii.gz selection/
                    fi
                    j=$(( j + 1 ))
                done
            fi
        elif [[ ${range} -gt ${nfmribins} && ${range} -lt 32 ]];
            then
            for ((i = 0; i < ${nfmribins} ; i++))
                do
                number=$RANDOM
                number=$((number %= $range))
                if [[ ${number} -lt 10 ]];
                    then
                    ImageMath $TDIM selection/vol0.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                    #   cp vol000${number}.nii.gz selection/
                elif [[ ${number} -ge 10 && ${number} -lt 100 ]];
                    then
                    ImageMath $TDIM selection/vol0.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                    #   cp vol00${number}.nii.gz selection/
                fi
            done
        elif [[ ${range} -le ${nfmribins} ]];
            then
            ImageMath selection/$TDIM vol0.nii.gz TimeSeriesSubset ${IMAGESETVARIABLE} ${range}
            # cp *.nii.gz selection/
        fi
        # set filelist variable
        rm -f ${IMAGESETVARIABLE}
        cd selection/
        IMAGESETVARIABLE=`ls *.nii.gz`

        IMAGESETARRAY=()
        for IMG in $IMAGESETVARIABLE
          do
            IMAGESETARRAY[${#IMAGESETARRAY[@]}]=$IMG
          done
    fi
else
    IMAGESETARRAY=()
    for IMG in $IMAGESETVARIABLE
      do
        IMAGESETARRAY[${#IMAGESETARRAY[@]}]=$IMG
      done
fi

if [[ $NUMBEROFMODALITIES -gt 1 ]];
    then
    echo "--------------------------------------------------------------------------------------"
    echo " Multivariate template construction using the following ${NUMBEROFMODALITIES}-tuples:  "
    echo "--------------------------------------------------------------------------------------"
    for (( i = 0; i < ${#IMAGESETARRAY[@]}; i+=$NUMBEROFMODALITIES ))
      do
        IMAGEMETRICSET=""
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
          do
            k=0
            k=$((i+j))
            IMAGEMETRICSET="$IMAGEMETRICSET ${IMAGESETARRAY[$k]}"
          done
        echo $IMAGEMETRICSET
      done
    echo "--------------------------------------------------------------------------------------"
fi

# Useful to check the right number of images exist for various ops
IMAGESPERMODALITY=$(( ${#IMAGESETARRAY[@]} / ${NUMBEROFMODALITIES} ))

# check for initial template images
for (( i = 0; i < $NUMBEROFMODALITIES; i++ ))
  do
    setCurrentImageSet $i

    if [[ ${#REGTEMPLATES[@]} -gt 0 && -n "${REGTEMPLATES[$i]}" ]];
      then
        if [[ ! -r "${REGTEMPLATES[$i]}" ]];
          then
            echo "Initial template ${REGTEMPLATES[$i]} cannot be read"
            exit 1
          fi
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Initial template $i found.  This will be used for guiding the registration. use : ${REGTEMPLATES[$i]} and ${TEMPLATES[$i]} "
        echo "--------------------------------------------------------------------------------------"
     # now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.
        cp ${REGTEMPLATES[$i]} ${TEMPLATES[$i]}
      else
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Creating template ${TEMPLATES[$i]} from a population average image from the inputs."
        echo "   ${CURRENTIMAGESET[@]}"
        echo "--------------------------------------------------------------------------------------"
        # Normalized mean, no sharpening
        # This forces a call to AverageImages, which resizes images to match the largest input
        summarizeimageset $DIM ${TEMPLATES[$i]} 1 0 ${CURRENTIMAGESET[@]}
        # Quickly align COM of input images to average, and then recompute average
        IMAGECOMSET=()
        for (( j = 0; j < ${#CURRENTIMAGESET[@]}; j+=1 ))
          do
            IMGbase=`basename ${CURRENTIMAGESET[$j]}`
            BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
            COM="${OUTPUT_DIR}/initialCOM${i}_${j}_${IMGbase}"
            COMTRANSFORM="${OUTPUT_DIR}/initialCOM${i}_${j}_${BASENAME}.mat"
            antsAI -d ${DIM} --convergence 0 --verbose 1 -m Mattes[${TEMPLATES[$i]},${CURRENTIMAGESET[$j]},32,None] -o ${COMTRANSFORM} -t AlignCentersOfMass
            antsApplyTransforms -d ${DIM} -r ${TEMPLATES[$i]} -i ${CURRENTIMAGESET[$j]} -t ${COMTRANSFORM} -o ${COM} --verbose
            rm -f $COMTRANSFORM
            IMAGECOMSET[${#IMAGECOMSET[@]}]=$COM
          done
        # Now safe to let user control stat method
        summarizeimageset $DIM ${TEMPLATES[$i]} $STATSMETHOD 0 ${IMAGECOMSET[@]}
        # Clean up
        rm -f ${IMAGECOMSET[@]}
      fi

    if [[ ! -s ${TEMPLATES[$i]} ]];
      then
        echo "Your initial template : $TEMPLATES[$i] was not created.  This indicates trouble!  You may want to check correctness of your input parameters. exiting."
        exit 1
      fi

    # Back up template
    intermediateTemplateBase=`basename ${TEMPLATES[$i]}`
    cp ${TEMPLATES[$i]} ${intermediateTemplateDir}/initial_${intermediateTemplateBase}

done


# remove old job bash scripts
outdir=`dirname ${TEMPLATES[0]}`
if [[ ${#outdir} -eq 0 ]];
    then
    outdir=`pwd`
fi
rm -f ${outdir}/job*.sh

##########################################################################
#
# perform rigid body registration if requested
#
##########################################################################
if [[ "$RIGID" -eq 1 ]];
  then
    count=0
    jobIDs=""

    for (( i = 0; i < ${#IMAGESETARRAY[@]}; i+=$NUMBEROFMODALITIES ))
      do

        basecall="${ANTS} -d ${DIM} --float $USEFLOAT --verbose 1 -u 0 -w [ 0.01,0.99 ] -z 1 -r [ ${TEMPLATES[0]},${IMAGESETARRAY[$i]},1 ]"

        IMAGEMETRICSET=""
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
          do
            k=0
            k=$((i+j))
            IMAGEMETRICSET="$IMAGEMETRICSET -m MI[ ${TEMPLATES[$j]},${IMAGESETARRAY[$k]},${MODALITYWEIGHTS[$j]},32,Regular,0.25 ]"
          done

        stage1="-t Rigid[ 0.1 ] ${IMAGEMETRICSET} -c [ 1000x500x250x0,1e-6,10 ] -f 6x4x2x1 -s 3x2x1x0 -o ${outdir}/rigid${i}_"
        #stage1="-t Rigid[ 0.1 ] ${IMAGEMETRICSET} -c [ 10x10x10x10,1e-8,10 ] -f 8x4x2x1 -s 4x2x1x0 -o ${outdir}/rigid${i}_"
        exe="${basecall} ${stage1}"

        qscript="${outdir}/job_${count}_qsub.sh"
        rm -f $qscript

        if [[ $DOQSUB -eq 5 ]];
            then
            # SLURM job scripts must start with a shebang
            echo '#!/bin/sh' > $qscript
            fi

        echo "$SCRIPTPREPEND" >> $qscript

        IMGbase=`basename ${IMAGESETARRAY[$i]}`
        BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
        RIGID="${outdir}/rigid${i}_0_${IMGbase}"

        echo "$exe" >> $qscript

        exe2='';
        pexe2='';
        pexe=" $exe > ${outdir}/job_${count}_metriclog.txt "
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
          do
            k=0
            k=$((i+j))
            IMGbase=`basename ${IMAGESETARRAY[$k]}`
            BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
            RIGID="${outdir}/rigid${i}_${j}_${IMGbase}"
            IMGbaseBASE=`basename ${IMAGESETARRAY[$i]}`
            BASENAMEBASE=` echo ${IMGbaseBASE} | cut -d '.' -f 1 `
            exe2="$exe2 ${WARP} -d $DIM --float $USEFLOAT --verbose 1 -i ${IMAGESETARRAY[$k]} -o $RIGID -t ${outdir}/rigid${i}_0GenericAffine.mat -r ${TEMPLATES[$j]}\n"
            pexe2="$exe2 ${WARP} -d $DIM --float $USEFLOAT --verbose 1 -i ${IMAGESETARRAY[$k]} -o $RIGID -t ${outdir}/rigid${i}_0GenericAffine.mat -r ${TEMPLATES[$j]} >> ${outdir}/job_${count}_metriclog.txt\n"
          done

        echo -e "$exe2" >> $qscript;

        if [[ $DOQSUB -eq 1 ]];
          then
            id=`qsub -cwd -S /bin/bash -N antsBuildTemplate_rigid -v  $QSUBOPTS $qscript | awk '{print $3}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
          then
            id=`qsub -N antsrigid -v  $QSUBOPTS -q nopreempt -l nodes=1:ppn=1 -l mem=${MEMORY} -l walltime=${WALLTIME} $qscript | awk '{print $1}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 2 ]];
          then
            # Send pexe and exe2 to same job file so that they execute in series
            echo $pexe >> ${outdir}/job${count}_r.sh
            echo -e $pexe2 >> ${outdir}/job${count}_r.sh
        elif [[ $DOQSUB -eq 3 ]];
          then
            id=`xgrid $XGRIDOPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
            #echo "xgrid $XGRIDOPTS -job submit /bin/bash $qscript"
            jobIDs="$jobIDs $id"
        elif [[ $DOQSUB -eq 5 ]];
            then
            id=`sbatch --job-name=antsrigid --export=${QSUBOPTS} --nodes=1 --cpus-per-task=1 --time=${WALLTIME} --mem=${MEMORY} $qscript | rev | cut -f1 -d\ | rev`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
          then
            echo $qscript
             # execute jobs in series
            bash $qscript
          fi
        count=$(( count + 1 ))
    done
    if [[ $DOQSUB -eq 1 ]];
      then
        # Run jobs on SGE and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on SGE cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
        waitForSGEQJobs.pl 1 60 $jobIDs
        # Returns 1 if there are errors
        if [[ ! $? -eq 0 ]];
          then
            echo "qsub submission failed - jobs went into error state"
            exit 1;
          fi
      fi
    if [[ $DOQSUB -eq 4 ]];
      then
        # Run jobs on PBS and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on PBS cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
               # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
        waitForPBSQJobs.pl 1 60 $jobIDs
        # Returns 1 if there are errors
        if [[ ! $? -eq 0 ]];
          then
            echo "qsub submission failed - jobs went into error state"
            exit 1;
          fi
      fi
    # Run jobs on localhost and wait to finish
    if [[ $DOQSUB -eq 2 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on max ${CORES} cpucores. "
        echo " Progress can be viewed in ${outdir}/job*_metriclog.txt"
        echo "--------------------------------------------------------------------------------------"
        jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
        chmod +x ${outdir}/job*_r.sh
        $PEXEC -j ${CORES} "sh" ${outdir}/job*_r.sh
      fi
    if [[ $DOQSUB -eq 3 ]];
      then
        # Run jobs on XGrid and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on XGrid cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
        waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 30 $jobIDs
        # Returns 1 if there are errors
        if [[ ! $? -eq 0 ]];
          then
            echo "XGrid submission failed - jobs went into error state"
            exit 1;
          fi
      fi
    if [[ $DOQSUB -eq 5 ]];
      then
        # Run jobs on SLURM and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on SLURM cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
               # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
        waitForSlurmJobs.pl 1 60 $jobIDs
        # Returns 1 if there are errors
        if [[ ! $? -eq 0 ]];
          then
            echo "SLURM submission failed - jobs went into error state"
            exit 1;
          fi
      fi

    for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
      do
        IMAGERIGIDSET=()
        for (( i = $j; i < ${#IMAGESETARRAY[@]}; i+=$NUMBEROFMODALITIES ))
          do
            k=0
            k=$((i-j))
            IMGbase=`basename ${IMAGESETARRAY[$i]}`
            BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
            RIGID="${outdir}/rigid${k}_${j}_${IMGbase}"

            IMAGERIGIDSET[${#IMAGERIGIDSET[@]}]=$RIGID
          done
        echo
        echo  "Building rigid template from ${IMAGERIGIDSET[@]}"

        # No sharpening at rigid stage
        summarizeimageset $DIM ${TEMPLATES[$j]} $STATSMETHOD 0 ${IMAGERIGIDSET[@]}
        intermediateTemplateBase=`basename ${TEMPLATES[$j]}`
        cp ${TEMPLATES[$j]} ${intermediateTemplateDir}/initialRigid_${intermediateTemplateBase}

      done

    if [[ BACKUPEACHITERATION -eq 1 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from rigid iteration"
        echo "--------------------------------------------------------------------------------------"

        mkdir ${outdir}/rigid
        mv ${outdir}/rigid*.nii.gz ${outdir}/*GenericAffine.mat ${outdir}/rigid/
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
          then
            mv ${outdir}/antsBuildTemplate_rigid* ${outdir}/rigid/
            # Remove qsub scripts
            rm -f ${outdir}/job_${count}_qsub.sh
        elif [[ $DOQSUB -eq 4 ]];
          then
            mv ${outdir}/antsrigid* ${outdir}/job* ${outdir}/rigid/
        elif [[ $DOQSUB -eq 2 ]];
          then
            mv ${outdir}/job*.txt ${outdir}/rigid/
        elif [[ $DOQSUB -eq 3 ]];
          then
            rm -f ${outdir}/job_*_qsub.sh
        elif [[ $DOQSUB -eq 5 ]];
          then
            mv ${outdir}/slurm-*.out ${outdir}/rigid/
            mv ${outdir}/job*.txt ${outdir}/rigid/

            # Remove submission scripts
            rm -f ${outdir}/job_${count}_qsub.sh
        fi
      else
        rm -f  ${outdir}/rigid*.* ${outdir}/job*.txt ${outdir}/slurm-*.out
    fi
fi # endif RIGID

##########################################################################
#
# begin main level
#
##########################################################################

ITERATLEVEL=( $(echo $MAXITERATIONS | tr 'x' '\n') )
NUMLEVELS=${#ITERATLEVEL[@]}

SHRINKLEVEL=( $(echo $SHRINKFACTORS | tr 'x' '\n') )
NUMSHRINK=${#SHRINKLEVEL[@]}

SMOOTHLEVEL=( $(echo $SMOOTHINGFACTORS | tr 'x' '\n') )
NUMSMOOTH=${#SMOOTHLEVEL[@]}

if [[ $NUMLEVELS -ne $NUMSHRINK ]]
  then
    echo "Number of shrink factors in [ $SHRINKFACTORS ] does not match number of iteration levels in [ $MAXITERATIONS ]"
    exit 1
  fi

if [[ $NUMLEVELS -ne $NUMSMOOTH ]]
  then
    echo "Number of smoothing levels in [ $SMOOTHINGFACTORS ] does not match number of iteration levels in [ $MAXITERATIONS ]"
    exit 1
  fi

#
# debugging only
#echo $ITERATLEVEL
#echo $NUMLEVELS
#echo ${ITERATIONLIMIT}
#
echo
echo "--------------------------------------------------------------------------------------"
echo " Start to build templates: ${TEMPLATES[@]}"
echo "--------------------------------------------------------------------------------------"
#


reportMappingParameters

i=0
while [[ $i -lt ${ITERATIONLIMIT} ]];
  do
    itdisplay=$((i+1))
    rm -f ${OUTPUTNAME}*WarpedToTemplate.nii.gz
    rm -f ${OUTPUTNAME}*Warp.nii*
    rm -f ${OUTPUTNAME}*warp.nii*
    rm -f ${OUTPUTNAME}*GenericAffine.mat
    rm -f ${outdir}/job*.sh
    # Used to save time by only running coarse registration for the first couple of iterations
    # This may also help convergence, but because there's no way to turn it off, it makes it harder
    # to refine templates with multiple calls to this script.
    # If you uncomment this, replace MAXITERATIONS with ITERATIONS in the call to ants below
    #
    # if [[ $i -gt $((NUMLEVELS - 1)) ]];
    #    then
    #    ITERATIONS=$MAXITERATIONS
    # else
    #
    #    ITERATIONS=${ITERATLEVEL[0]}
    #
    #    for (( n = 1 ; n < ${NUMLEVELS}; n++ ))
    #      do
    #      ITERATIONS=${ITERATIONS}x$((${ITERATLEVEL[n]} * $((n <= i)) ))
    #    done
    # fi
    # Job IDs of jobs submitted to queue in loop below
    jobIDs=""
    # Reinitialize count to 0
    count=0
    # Submit registration of each input to volume template to SGE or run locally.

    for (( j = 0; j < ${#IMAGESETARRAY[@]}; j+=$NUMBEROFMODALITIES ))
      do
        basecall="${ANTS} -d ${DIM} --float $USEFLOAT --verbose 1  -u 0 -w [ 0.01,0.99 ] -z 1"

        IMAGEMETRICLINEARSET=''
        IMAGEMETRICSET=''
        exe=''
        warpexe=''
        pexe=''
        warppexe=''

        for (( k = 0; k < $NUMBEROFMODALITIES; k++ ))
          do
            l=0
            l=$((j+k))

            if [[ "${METRICTYPE[$k]}" == "DEMONS" ]];
              then
                # Mapping Parameters
                METRIC="Demons[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},4 ]"
            elif [[ "${METRICTYPE[$k]}" == CC* ]];
              then
                METRIC="CC[ "
                RADIUS=4
                if [[ "${METRICTYPE[$k]}" == CC[* ]]
                  then
                    RADIUS=${METRICTYPE[$k]%]*}
                    RADIUS=${RADIUS##*[}
                  fi
                METRICPARAMS="${MODALITYWEIGHTS[$k]},${RADIUS} ]"
            elif [[ "${METRICTYPE[$k]}" == "MI" ]];
              then
                # Mapping Parameters
                METRIC="MI[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},32 ]"
            elif [[ "${METRICTYPE[$k]}" == "MSQ" ]];
              then
                # Mapping Parameters
                METRIC="MeanSquares[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},0 ]"
            else
              echo "Invalid similarity metric. Use CC, MI, MSQ, DEMONS or type bash `basename $0` -h."
              exit 1
            fi
            TEMPLATEbase=`basename ${TEMPLATES[$k]}`
            indir=`dirname ${IMAGESETARRAY[$j]}`
            if [[ ${#indir} -eq 0 ]];
              then
                indir=`pwd`
              fi
            IMGbase=`basename ${IMAGESETARRAY[$l]}`
            OUTFN=${OUTPUTNAME}input$(printf "%04d" $l)-modality${k}-${IMGbase/%?(.nii.gz|.nii)}
            OUTFN=`basename ${OUTFN}`
            DEFORMED="${outdir}/${OUTFN}-WarpedToTemplate.nii.gz"

            IMGbase=`basename ${IMAGESETARRAY[$j]}`
            OUTWARPFN=${OUTPUTNAME}input$(printf "%04d" $j)-${IMGbase/%?(.nii.gz|.nii)}-
            OUTWARPFN=`basename ${OUTWARPFN}`

            if [[ $NOWARP -eq 0 ]];
              then
                OUTPUTTRANSFORMS="-t ${outdir}/${OUTWARPFN}1Warp.nii.gz -t ${outdir}/${OUTWARPFN}0GenericAffine.mat"
              else
                OUTPUTTRANSFORMS="-t ${outdir}/${OUTWARPFN}0GenericAffine.mat"
              fi

            if [[ $N4CORRECT -eq 1 ]];
              then
                REPAIRED="${outdir}/${OUTFN}Repaired.nii.gz"
                if [[ ! -s ${REPAIRED} ]]; then
                  exe=" $exe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${IMAGESETARRAY[$l]} -o ${REPAIRED} -r 0 -s 2 --verbose 1\n"
                  pexe=" $pexe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${IMAGESETARRAY[$l]} -o ${REPAIRED} -r 0 -s 2 --verbose 1  >> ${outdir}/job_${count}_metriclog.txt >> ${outdir}/job_${count}_metriclog.txt\n"
                fi
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${TEMPLATES[$k]},${REPAIRED},${METRICPARAMS}"
                IMAGEMETRICLINEARSET="$IMAGEMETRICLINEARSET -m MI[ ${TEMPLATES[$k]},${REPAIRED},${MODALITYWEIGHTS[$k]},32,Regular,0.25 ]"

                warpexe=" $warpexe ${WARP} -d ${DIM} --float $USEFLOAT --verbose 1 -i ${REPAIRED} -o ${DEFORMED} -r ${TEMPLATES[$k]} ${OUTPUTTRANSFORMS}\n"
                warppexe=" $warppexe ${WARP} -d ${DIM} --float $USEFLOAT --verbose 1 -i ${REPAIRED} -o ${DEFORMED} -r ${TEMPLATES[$k]} ${OUTPUTTRANSFORMS} >> ${outdir}/job_${count}_metriclog.txt\n"
              else
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${TEMPLATES[$k]},${IMAGESETARRAY[$l]},${METRICPARAMS}"
                IMAGEMETRICLINEARSET="$IMAGEMETRICLINEARSET -m MI[ ${TEMPLATES[$k]},${IMAGESETARRAY[$l]},${MODALITYWEIGHTS[$k]},32,Regular,0.25 ]"

                warpexe=" $warpexe ${WARP} -d ${DIM} --float $USEFLOAT --verbose 1 -i ${IMAGESETARRAY[$l]} -o ${DEFORMED} -r ${TEMPLATES[$k]} ${OUTPUTTRANSFORMS}\n"
                warppexe=" $warppexe ${WARP} -d ${DIM} --float $USEFLOAT --verbose 1 -i ${IMAGESETARRAY[$l]} -o ${DEFORMED} -r ${TEMPLATES[$k]} ${OUTPUTTRANSFORMS} >> ${outdir}/job_${count}_metriclog.txt\n"
              fi

        done


      #  Already defined above
      #  IMGbase=`basename ${IMAGESETARRAY[$j]}`
      #  OUTWARPFN=${OUTPUTNAME}${IMGbase/%?(.nii.gz|.nii)}
      #  OUTWARPFN=`basename ${OUTWARPFN}${j}`

        stage0="-r [ ${TEMPLATES[0]},${IMAGESETARRAY[$j]},1 ]"
        stage1="-t Rigid[ 0.1 ] ${IMAGEMETRICLINEARSET} -c [ 1000x500x250x0,1e-6,10 ] -f 6x4x2x1 -s 4x2x1x0"
        stage2="-t Affine[ 0.1 ] ${IMAGEMETRICLINEARSET} -c [ 1000x500x250x0,1e-6,10 ] -f 6x4x2x1 -s 4x2x1x0"
        #stage1="-t Rigid[ 0.1 ] ${IMAGEMETRICLINEARSET} -c [ 10x10x10x10,1e-8,10 ] -f 8x4x2x1 -s 4x2x1x0"
        #stage2="-t Affine[ 0.1 ] ${IMAGEMETRICLINEARSET} -c [ 10x10x10x10,1e-8,10 ] -f 8x4x2x1 -s 4x2x1x0"
        stage3="-t ${TRANSFORMATION} ${IMAGEMETRICSET} -c [ ${MAXITERATIONS},1e-9,10 ] -f ${SHRINKFACTORS} -s ${SMOOTHINGFACTORS} -o ${outdir}/${OUTWARPFN}"

        stageId="-t Rigid[ 0.1 ] ${IMAGEMETRICLINEARSET} -c [ 0,1e-8,10 ] -f 1 -s 0"
        exebase=$exe
        pexebase=$pexe

        if [[ $DOLINEAR -eq 0 ]];
          then
            exe="$exe ${basecall} ${stageId} ${stage3}\n"
            pexe="$pexe ${basecall} ${stageId} ${stage3} >> ${outdir}/job_${count}_metriclog.txt\n"
          elif [[ $NOWARP -eq 1 ]];
            then
    	      if [[ ${TRANSFORMATION} == "Affine"* ]];
	        then
          	  # If affine, do standard rigid, then affine with levels, etc from command line
		  exe="$exebase ${basecall} ${stage0} ${stage1} ${stage3}\n";
		  pexe="$pexebase ${basecall} ${stage0} ${stage1} ${stage3} >> ${outdir}/job_${count}_metriclog.txt\n"
	        else
		  # Rigid only - just use command line params
		  exe="$exebase ${basecall} ${stage0} ${stage3}\n";
		  pexe="$pexebase ${basecall} ${stage0} ${stage3} >> ${outdir}/job_${count}_metriclog.txt\n"
                fi
          else
            exe="$exe ${basecall} ${stage0} ${stage1} ${stage2} ${stage3}\n"
            pexe="$pexe ${basecall} ${stage0} ${stage1} ${stage2} ${stage3} >> ${outdir}/job_${count}_metriclog.txt\n"
          fi

        exe="$exe $warpexe"
        pexe="$pexe $warppexe"

        qscript="${outdir}/job_${count}_${i}.sh"

        echo -e $exe >> ${outdir}/job_${count}_${i}_metriclog.txt
        # 6 submit to SGE (DOQSUB=1), PBS (DOQSUB=4), PEXEC (DOQSUB=2), XGrid (DOQSUB=3), SLURM (DOQSUB=5) or else run locally (DOQSUB=0)
        if [[ $DOQSUB -eq 1 ]];
          then
            echo "$SCRIPTPREPEND" > $qscript
            echo -e "$exe" >> $qscript
            id=`qsub -cwd -N antsBuildTemplate_deformable_${i} -S /bin/bash -v  $QSUBOPTS $qscript | awk '{print $3}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
          then
            echo -e "$SCRIPTPREPEND" > $qscript
            echo -e "$exe" >> $qscript
            id=`qsub -N antsdef${i} -v  -q nopreempt -l nodes=1:ppn=1 -l mem=${MEMORY} -l walltime=${WALLTIME} $QSUBOPTS $qscript | awk '{print $1}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 2 ]];
          then
            echo -e $pexe >> ${outdir}/job${count}_r.sh
        elif [[ $DOQSUB -eq 3 ]];
          then
            echo -e "$SCRIPTPREPEND" > $qscript
            echo -e "$exe" >> $qscript
            id=`xgrid $XGRIDOPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
            jobIDs="$jobIDs $id"
        elif [[ $DOQSUB -eq 5 ]];
          then
            echo '#!/bin/sh' > $qscript
            echo -e "$SCRIPTPREPEND" >> $qscript
            echo -e "$exe" >> $qscript
            id=`sbatch --job-name=antsdef${i} --export=${QSUBOPTS} --nodes=1 --cpus-per-task=1 --time=${WALLTIME} --mem=${MEMORY} $QSUBOPTS $qscript | rev | cut -f1 -d\ | rev`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
          then
            echo -e $exe > $qscript
            bash $qscript
        fi

        # counter updated, but not directly used in this loop
        count=`expr $count + 1`;
    # echo " submitting job number $count " # for debugging only
    done
    # SGE wait for script to finish
    if [[ $DOQSUB -eq 1 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on SGE cluster. Iteration: $itdisplay of $ITERATIONLIMIT"
        echo "--------------------------------------------------------------------------------------"
        # now wait for the stuff to finish - this will take a while so poll queue every 10 mins
        waitForSGEQJobs.pl 1 600 $jobIDs
        if [[ ! $? -eq 0 ]];
          then
            echo "qsub submission failed - jobs went into error state"
            exit 1;
          fi
    elif [[ $DOQSUB -eq 4 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on PBS cluster. Iteration: $itdisplay of $ITERATIONLIMIT"
        echo "--------------------------------------------------------------------------------------"
        # now wait for the stuff to finish - this will take a while so poll queue every 10 mins
        waitForPBSQJobs.pl 1 600 $jobIDs
        if [[ ! $? -eq 0 ]];
          then
            echo "qsub submission failed - jobs went into error state"
            exit 1;
          fi
      fi
    # Run jobs on localhost and wait to finish
    if [[ $DOQSUB -eq 2 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on max ${CORES} cpucores. Iteration: $itdisplay of $ITERATIONLIMIT"
        echo " Progress can be viewed in job*_${i}_metriclog.txt"
        echo "--------------------------------------------------------------------------------------"
        jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
        chmod +x ${outdir}/job*.sh
        $PEXEC -j ${CORES} sh ${outdir}/job*.sh
      fi

    if [[ $DOQSUB -eq 3 ]];
      then
        # Run jobs on XGrid and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on XGrid cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. This is slow, so poll less often
        waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 300 $jobIDs
        # Returns 1 if there are errors
        if [[ ! $? -eq 0 ]];
          then
            echo "XGrid submission failed - jobs went into error state"
            exit 1;
          fi
      fi

    if [[ $DOQSUB -eq 5 ]];
      then
        # Run jobs on SLURM and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on SLURM cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the stuff to finish - this will take a while so poll queue every 10 mins
        waitForSlurmJobs.pl 1 600 $jobIDs
        if [[ ! $? -eq 0 ]];
          then
            echo "SLURM submission failed - jobs went into error state"
            exit 1;
          fi
      fi

    NUM_AFFINEFILES=0
    NUM_AFFINEFILES=$(ls ${OUTPUTNAME}*GenericAffine.mat | wc -l) || true

    if [[ $NOWARP -eq 1 ]];
      then
        if [[ ${NUM_AFFINEFILES} -eq 0 ]];
          then
            echo "The registrations did not terminate properly.  There are no warp files or affine files."
            exit 1
          fi
      if [[ ${NUM_AFFINEFILES} -ne $IMAGESPERMODALITY ]];
        then
          echo "The registrations did not terminate properly.  The number of affine transform files"
          echo "does not match the number of input images."
          exit 1
        fi
      else
        NUM_WARPFILES=0
        NUM_WARPFILES=$(ls ${OUTPUTNAME}*Warp.nii.gz | grep -v InverseWarp | wc -l) || true
        if [[ ${NUM_WARPFILES} -ne $IMAGESPERMODALITY ]];
          then
            echo "The registrations did not terminate properly.  The number of warp files"
            echo "does not match the number of input images."
            exit 1
          fi
      fi

    for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
      do
        shapeupdatetotemplate ${DIM} ${TEMPLATES[$j]} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP} ${j} ${STATSMETHOD} ${SHARPENMETHOD}
        # Back up templates and average warps if defined
        intermediateTemplateBase=`basename ${TEMPLATES[$j]}`
        cp ${TEMPLATES[$j]} ${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_${intermediateTemplateBase}
      done

    if [[ -f "${TEMPLATENAME}0warp.nii.gz" ]]
      then
        cp ${TEMPLATENAME}0warp.nii.gz ${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_shapeUpdateWarp.nii.gz
      fi

    if [[ $BACKUPEACHITERATION -eq 1 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from iteration $itdisplay"
        echo "--------------------------------------------------------------------------------------"
        mkdir ${outdir}/ANTs_iteration_${i}
        cp -f ${TEMPLATENAME}*warplog.txt ${OUTPUTNAME}*.nii.gz ${OUTPUTNAME}*.mat ${outdir}/ANTs_iteration_${i}/
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
            then
            mv ${outdir}/antsBuildTemplate_deformable_* ${outdir}/ANTs_iteration_${i}
        elif [[ $DOQSUB -eq 4 ]];
            then
            mv ${outdir}/antsdef* ${outdir}/ANTs_iteration_${i}
        elif [[ $DOQSUB -eq 2 ]];
            then
            mv ${outdir}/job*.txt ${outdir}/ANTs_iteration_${i}
        elif [[ $DOQSUB -eq 3 ]];
            then
            rm -f ${outdir}/job_*.sh
        elif [[ $DOQSUB -eq 5 ]];
            then
            mv ${outdir}/slurm-*.out ${outdir}/ANTs_iteration_${i}
            mv ${outdir}/job*.txt ${outdir}/ANTs_iteration_${i}
        fi
      else
        rm -f ${outdir}/job*.txt ${outdir}/slurm-*.out
    fi
    echo "Iteration $itdisplay completed"
    i=$(( i + 1 ))
done

# end main loop

rm -f job*.sh
#cleanup of 4D files
if [[ "${range}" -gt 1 && "${TDIM}" -eq 4 ]];
  then
    mv ${tmpdir}/selection/${TEMPLATES[@]} ${currentdir}/
    cd ${currentdir}
    rm -rf ${tmpdir}/
  fi
time_end=`date +%s`
time_elapsed=$((time_end - time_start))
echo
echo "--------------------------------------------------------------------------------------"
echo " Done creating: ${TEMPLATES[@]}"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo " Intermediate templates and warps stored in: ${intermediateTemplateDir}/"
echo "--------------------------------------------------------------------------------------"
echo " Note: The template has been updated after the final iteration."
echo " Pairwise warps from the final iteration are left in ${outdir}."
echo " They register input to the penultimate template(s), stored in ${intermediateTemplateDir}/"
echo "--------------------------------------------------------------------------------------"

exit 0
