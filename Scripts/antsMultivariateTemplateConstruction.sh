#!/bin/bash

set -euo pipefail

# Use extended globbing to allow for more flexible pattern matching
shopt -s extglob

VERSION="0.0.0"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTS=ANTS
WARP=WarpImageMultiTransform
N4=N4BiasFieldCorrection
PEXEC=ANTSpexec.sh
SGE=waitForSGEQJobs.pl
PBS=waitForPBSQJobs.pl
XGRID=waitForXGridJobs.pl
SLURM=waitForSlurmJobs.pl
PRINTHEADER=PrintHeader

fle_error=0
for FLE in "$ANTS" "$WARP" "$N4" "$PEXEC" "$SGE" "$XGRID" "$PBS" "$SLURM" "$PRINTHEADER"
  do
    if ! command -v "$FLE" &> /dev/null
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

function Usage {
    cat <<USAGE

Usage:

`basename $0` -d ImageDimension -o OUTPREFIX <other options> <images>

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
          ImageDimension: 4 (for template generation of time-series data).

          With -d 4, input should be a single 4D NIFTI image. Up to 16 3D volumes from
          the time series are selected and written to the output directory as
          time_index_XXXX.nii.gz. The template is then built from these 3D images.

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.  Optionally, one can specify a .csv or .txt file where each
          line is the location of the input image.  One can also specify more than
          one file for each image for multi-modal template construction (e.g. t1 and t2).
          For the multi-modal case, the templates will be consecutively numbered (e.g.
          \${OUTPUTPREFIX}template0.nii.gz, \${OUTPUTPREFIX}template1.nii.gz, ...).

NB: All images to be added to the template should be in the same directory, and this script
should be invoked from that directory.

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

     -c:  Control for parallel computation (default 1) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost, see -j), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM

          Use -T to control how many threads are used by each parallel job. For slurm,
          --cpus-per-task is set automatically. For SGE and PBS, you must specify the parallel
          environment options with -p in order to reserve the correct number of CPUs.

     -g:  Gradient step size (default 0.25) -- smaller in magnitude results in
          more cautious steps. Use smaller steps to refine template details.
          0.25 is an upper (aggressive) limit for this parameter.

     -i:  Iteration limit (default 4) -- iterations of the template construction (Iteration limit)*NumImages registrations.

     -j:  Number of parallel processes to use for pexec execution (-c 2). The default is 2. Each process may use
          multiple threads, see the -T option. For example, if you have 8 cores, you could use -j 4 -T 2 to run
          4 parallel processes, each using 2 threads, or -j 8 -T 1 to run 8 parallel processes, each using 1 thread.

          If not using pexec, this option is ignored.

     -T:  Number of ITK threads to use in registration and other ANTs programs. For cluster jobs, this value will be
          set inside the job script.

          A reasonable range for this value is 1 to 8. The optimal setting depends on the system and data but in general,
          prefer more parallel processes and fewer threads.

          For slurm, this script will set the '--cpus-per-task' option in the sbatch command line.
          For SGE and PBS, you must specify the parallel environment options with -p in order to reserve the
          correct number of CPUs.

          For parallel execution, the default is to use the environment variable ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
          if it is defined, or 1 otherwise.

          For serial execution (-c 0), the default is to use the environment variable ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
          if it is defined, otherwise the min(8, number of system CPU cores).

          If the number of threads is set to 0, or if ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=0, then the script will use
          min(8, number of system CPU cores).

     -k:  Number of modalities used to construct the template (default 1)

     -w:  Modality weights used in the similarity metric (default = 1) --- specified as e.g. 1x0.5x0.75

     -m:  Max-iterations in each registration

     -n:  N4BiasFieldCorrection of moving image (default 1) -- 0 == off, 1 == on

     -p:  A command or directive to prepend to job scripts. Repeat -p to add multiple lines; entries
          are written to job scripts in the order given. Job options for SGE/PBS/SLURM can be included here,
          but command-line args in this script will take precedence for wall time (see -u), memory (see -v),
          and multithreading (see -T).

          For slurm, a shebang for '/bin/sh' is always the first line in the job script. For SGE and PBS,
          the shell interpreter is set to '/bin/bash'.

     -r:  Do rigid-body registration of inputs to the initial template, before doing the main
          pairwise registration. 0 == off 1 == on (default 0). If you are trying to refine or update
          an existing template, you would use '-r 0'.
          Rigid initialization is useful when you do not have an initial template, or you want to use
          a single image as a reference for rigid alignment only. For example,
            "-z tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz -y 0 -r 1"
          will rigidly align the inputs to the MNI template, and then use their average to begin the
          template building process.

     -s:  Type of similarity metric used for nonlinear registration (affine is always MI). Default = CC.
          Options are case sensitive.
             CC  : Cross-correlation
             MI  : Mutual information
             MSQ : Mean squared differences
             PR  : CC after subtraction of local mean from the image (deprecated)

     -t:  Type of transformation model used for nonlinear registration. Options are case sensitive.
             GR             : Greedy SyN (default for scalar data)
             GR_Constrained : Greedy SyN with regularization on the total deformation (default for time series)
             EL             : Elastic
             EX             : Exponential
             DD             : Greedy exponential, diffemorphic-demons-style optimization
             SY             : LDDMM-style SyN with symmetric time-dependent gradient estimation
             LDDMM          : Like SY, but with asymmetric time-dependent gradient estimation
             S2             : Like SY, but with no time-dependent gradient estimation

     -u:  Walltime (default = 20:00:00):  Option for PBS/SLURM qsub specifying requested time
          per pairwise registration.

     -v:  Memory limit (default = 8G):  Option for PBS/SLURM qsub specifying requested memory
          per pairwise registration.

     -x:  XGrid arguments (e.g., -x "-p password -h controlhost")

     -y:  Update the template with the full affine transform (default 1). If 0, the rigid
          component of the affine transform will not be used to update the template. If your
          template drifts in translation or orientation try -y 0.

     -z:  Use this this volume as the target of all inputs. When not used, the script will create an unbiased
          starting point by averaging all inputs, then aligning the center of mass of all inputs to that of
          the initial average. If you do not use -z, it is recommended to use "-r 1". Use the full path.
          For multiple modalities, specify -z modality1.nii.gz -z modality2.nii.gz ...
          in the same modality order as the input images.

     -b:  Boolean for saving full iteration output to directories (default = 0). If 1, images and warps
          are saved for each pairwise registration at each iteration. Otherwise, only templates and the shape
          update warps are saved.


Example:

`basename $0` -d 3 -m 30x50x20 -t GR -s CC -c 1 -o MY -z InitialTemplate.nii.gz  *RF*T1x.nii.gz

- In this example 30x50x20 iterations per registration are used for template creation (that is the default)
- Greedy-SyN and CC are the metrics to guide the mapping.
- Output is prepended with MY and the initial template is InitialTemplate.nii.gz (optional).
- The -c option is set to 1, which will result in using the Sun Grid Engine (SGE) to distribute the computation.
- if you do not have SGE, read the help for multi-core computation on the local machine, or Apple X-grid options.

Output:

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
script adapted by N.M. van Strien, http://www.mri-tutorial.com | NTNU MR-Center
multivariate template adaption by Nick Tustison
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
 Dimensionality:                    $DIM
 N4BiasFieldCorrection:             $N4CORRECT
 Similarity Metric:                 $METRICTYPE
 Transformation:                    $TRANSFORMATIONTYPE
 Regularization:                    $REGULARIZATION
 MaxIterations:                     $MAXITERATIONS
 Number Of MultiResolution Levels:  $NUMLEVELS
 OutputName prefix:                 $OUTPUTNAME
 Template:                          $TEMPLATENAME
 Template Update Steps:             $ITERATIONLIMIT
 Template population:               $IMAGESETVARIABLE
 Number of Modalities:              $NUMBEROFMODALITIES
 Modality weights:                  $MODALITYWEIGHTSTRING
 Image statistic:                   $STATSMETHOD
 Sharpening method:                 $SHARPENMETHOD
 Shape update full affine:          $AFFINE_UPDATE_FULL
--------------------------------------------------------------------------------------
REPORTMAPPINGPARAMETERS
}

function write_script_prepend() {
  if (( ${#SCRIPTPREPEND[@]} > 0 )); then
    printf '%s\n' "${SCRIPTPREPEND[@]}"
  fi
}

function get_image_stem() {
  local filename=${1##*/}

  # Treat .gz as a compression suffix, then remove the image-format suffix.
  # This preserves internal periods for both compressed and uncompressed formats.
  if [[ $filename == *.gz ]]; then
    filename=${filename%.gz}
  fi
  printf '%s\n' "${filename%.*}"
}

# Set RANDOM_INDEX to an unbiased random integer in [0, upper_bound).
# Combining two values of RANDOM gives a uniform 30-bit source value; rejection
# sampling avoids the modulo bias that occurs when upper_bound does not divide
# the source range evenly.
function random_index_below() {
  local upper_bound=$1
  local random_value
  local rejection_limit
  local random_source_size=$((1 << 30))

  if (( upper_bound < 1 || upper_bound > random_source_size )); then
    printf 'random_index_below: upper bound must be between 1 and %d\n' "$random_source_size" >&2
    return 1
  fi

  rejection_limit=$((random_source_size - random_source_size % upper_bound))
  while true; do
    random_value=$((RANDOM << 15 | RANDOM))
    if (( random_value < rejection_limit )); then
      RANDOM_INDEX=$((random_value % upper_bound))
      return 0
    fi
  done
}

function cleanup_4d_tempdir() {
  if [[ -n ${tmpdir:-} && -d $tmpdir ]]; then
    rm -rf "$tmpdir"
    tmpdir=""
  fi
}
trap cleanup_4d_tempdir EXIT

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
      AverageImages "$dim" "$output" 0 "${images[@]}"
      ;;
    1) #mean of normalized images
      AverageImages "$dim" "$output" 2 "${images[@]}"
      ;;
    2) #median
      local image
      for image in "${images[@]}";
        do
          printf '%s\n' "$image" >> "${output}_list.txt"
        done
      ImageSetStatistics "$dim" "${output}_list.txt" "$output" 0
      rm "${output}_list.txt"
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
      ImageMath "$dim" "$output" Sharpen "$output" 0
      ;;
    2)
      echo "Unsharp mask sharpening"
      ImageMath "$dim" "$output" UnsharpMask "$output" 0.5 1 0 0
      ;;
  esac

}


function shapeupdatetotemplate() {

   echo "shapeupdatetotemplate()"

    # local declaration of values
    dim=$1
    template=$2
    templatename=$3
    outputname=$4
    gradientstep=-$5
    summarizemethod=$6
    sharpenmethod=$7
    whichtemplate=$8

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
    echo "--------------------------------------------------------------------------------------"

    imagelist=( "${outputname}"template-modality"${whichtemplate}"-*WarpedToTemplate.nii.gz )
    if [[ ${#imagelist[@]} -eq 0 ]]; then
      imagelist=()
    elif [[ ! -e ${imagelist[0]} ]]; then
      imagelist=()
    fi
    if [[ ${#imagelist[@]} -ne ${IMAGESPERMODALITY} ]]
      then
        echo "ERROR shapeupdatedtotemplate - imagelist length is ${#imagelist[@]}, expected ${IMAGESPERMODALITY}"
        exit 1
      fi

    summarizeimageset "${dim}" "${template}" "${summarizemethod}" "${sharpenmethod}" "${imagelist[@]}"

    if [[ $whichtemplate -eq 0 ]] ;
      then
        WARPFILES=()
        for warpFile in "${outputname}"*Warp.nii.gz; do
          if [[ -e $warpFile && $warpFile != *InverseWarp* ]]; then
            WARPFILES+=( "$warpFile" )
          fi
        done
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---voxel-wise averaging of the inverse warp fields (from subject to template)"
        echo "   AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 ${WARPFILES[*]}"
        echo "--------------------------------------------------------------------------------------"

        AverageImages "$dim" "${templatename}${whichtemplate}warp.nii.gz" 0 "${WARPFILES[@]}"

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---scale the averaged inverse warp field by the gradient step"
        echo "   MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz"
        echo "--------------------------------------------------------------------------------------"

        MultiplyImages "$dim" "${templatename}${whichtemplate}warp.nii.gz" "${gradientstep}" "${templatename}${whichtemplate}warp.nii.gz"

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---average the affine transforms (template <-> subject)"
        echo "                      ---transform the inverse field by the resulting average affine transform"
        echo "   ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0Affine.txt ${outputname}*Affine.txt"
        echo "   WarpImageMultiTransform ${dim} ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -i  ${templatename}0Affine.txt -R ${template}"
        echo "--------------------------------------------------------------------------------------"

        "${AVERAGE_AFFINE_PROGRAM}" "${dim}" "${templatename}0Affine.txt" "${outputname}"*Affine.txt
        WarpImageMultiTransform "${dim}" "${templatename}0warp.nii.gz" "${templatename}0warp.nii.gz" -i "${templatename}0Affine.txt" -R "${template}"

        MeasureMinMaxMean "${dim}" "${templatename}0warp.nii.gz" "${templatename}warplog.txt" 1
      fi

    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate---warp each template by the resulting transforms"
    echo "   WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}0Affine.txt ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -R ${template}"
    echo "--------------------------------------------------------------------------------------"
    WarpImageMultiTransform "${dim}" "${template}" "${template}" -i "${templatename}0Affine.txt" "${templatename}0warp.nii.gz" "${templatename}0warp.nii.gz" "${templatename}0warp.nii.gz" "${templatename}0warp.nii.gz" -R "${template}"
}

function jobfnamepadding {

    outdir=`dirname "${TEMPLATES[0]}"`
    if [[ ${#outdir} -eq 0 ]]
        then
        outdir=`pwd`
    fi

    files=( "${outdir}"/job*.sh )
    if [[ ${#files[@]} -eq 0 ]]; then
      return 0
    elif [[ ! -e ${files[0]} ]]; then
      return 0
    fi
    BASENAME1=`printf '%s\n' "${files[0]}" | cut -d 'b' -f 1`

    for file in "${files[@]}"
      do

      if [[ "${#file}" -eq "9" ]];
       then
         BASENAME2=`printf '%s\n' "$file" | cut -d 'b' -f 2 `
         mv "$file" "${BASENAME1}b_000${BASENAME2}"

      elif [[ "${#file}" -eq "10" ]];
        then
          BASENAME2=`printf '%s\n' "$file" | cut -d 'b' -f 2 `
          mv "$file" "${BASENAME1}b_00${BASENAME2}"

      elif [[ "${#file}" -eq "11" ]];
        then
          BASENAME2=`printf '%s\n' "$file" | cut -d 'b' -f 2 `
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

control_c()
# run if user hits control-c
{
  trap - SIGINT
  printf '\n*** User pressed CTRL + C ***\n'
  cleanup_4d_tempdir
  printf '\n*** Script cancelled by user ***\n'
  exit 130
}

#initializing variables with global scope
time_start=`date +%s`
currentdir=`pwd`
nargs=$#

MAXITERATIONS=30x90x20
LABELIMAGE=0 # initialize optional parameter
METRICTYPE=()
TRANSFORMATIONTYPE="GR" # initialize optional parameter
DIM=""
NUMBEROFMODALITIES=1
MODALITYWEIGHTSTRING=""
N4CORRECT=1 # initialize optional parameter
DOQSUB=1 # By default, antsMultivariateTemplateConstruction tries to do things in parallel
GRADIENTSTEP=0.25 # Gradient step size, smaller in magnitude means more smaller (more cautious) steps
ITERATIONLIMIT=4
# Number of threads to use, we set this to -1 as a default, after arg parsing we will set
# this to a sensible value if it is not specified.
NUMBER_OF_THREADS=-1
PEXEC_PARALLEL_PROCESSES=2
TDIM=0
RIGID=0
RIGIDTYPE="" # set to an empty string to use affine initialization
LINEARTRANSFORMPARAMS="--number-of-affine-iterations 10000x10000x1000 --MI-option 32x16000"
range=0
REGTEMPLATES=()
TEMPLATES=()
CURRENTIMAGESET=()
XGRIDOPTS=""
SCRIPTPREPEND=()
WALLTIME="20:00:00"
MEMORY="8G"
# System specific queue options for SGE/PBS/SLURM, eg "-q name" to submit to a specific queue.
# It can be set to an empty string if you do not need any special cluster options.
# Some queue options are additive (eg, -l for resources) so check that your options are compatible
# with what's hard-coded below
QSUBOPTS="" # EDIT THIS
OUTPUTNAME=antsBTP
TEMPLATENAME=${OUTPUTNAME}template

BACKUP_EACH_ITERATION=0

AFFINE_UPDATE_FULL=1

# Methods for averaging warped images and sharpening next template
STATSMETHOD=1
SHARPENMETHOD=1

##Getting system info from linux can be done with these variables.
# RAM=`cat /proc/meminfo | sed -n -e '/MemTotal/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# RAMfree=`cat /proc/meminfo | sed -n -e '/MemFree/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# cpu_free_ram=$((${RAMfree}/${cpu_count}))

if [[ ${OSTYPE:0:6} == 'darwin' ]];
 then
 cpu_count=$(sysctl -n hw.physicalcpu 2>/dev/null) || cpu_count=1
else
 cpu_count=$(grep -c processor /proc/cpuinfo) || cpu_count=1
fi

# Provide output for Help
if [[ $# -eq 0 || ${1:-} == "-h" ]];
    then
    Usage >&2

fi

# reading command line arguments
while getopts "A:T:a:b:c:d:g:h:i:j:k:m:n:o:p:s:r:t:u:v:w:x:y:z:" OPT
  do

  case $OPT in
      A|T|a|b|c|d|e|i|j|k|l|n|r|y)
      if [[ ! $OPTARG =~ ^[0-9]+$ ]];
        then
          echo "Option -$OPT requires a non-negative integer, but received '$OPTARG'." >&2
          exit 1
        fi
      # Force base 10 so values with leading zeroes are safe in arithmetic expressions.
      OPTARG=$((10#$OPTARG))
      ;;
      g)
      if [[ ! $OPTARG =~ ^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]];
        then
          echo "Option -$OPT requires a numeric value, but received '$OPTARG'." >&2
          exit 1
        fi
      ;;
  esac

  case $OPT in
      h) #help
   Usage >&2
   exit 0
   ;;
      A) # Sharpening method
      SHARPENMETHOD=$OPTARG
   ;;
      T) # number of threads to use for each process
   NUMBER_OF_THREADS=$OPTARG
   ;;
      a) # summarizing statistic
      STATSMETHOD=$OPTARG
   ;;
      b) #backup each iteration (default = 0)
   BACKUP_EACH_ITERATION=$OPTARG
   ;;
      c) #use SGE cluster
   DOQSUB=$OPTARG
   if [[ ${#DOQSUB} -gt 2 ]]; then
       echo " DOQSUB must be an integer value (0=serial, 1=SGE qsub, 2=try pexec, 3=XGrid, 4=PBS qsub, 5=SLURM) you passed  -c $DOQSUB "
       exit 1
   fi
   ;;
      d) #dimensions
   DIM=$OPTARG
   if [[ ${DIM} -eq 4 ]]; then
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
      j) #number of parallel processes to use for pexec execution (default = 2)
   PEXEC_PARALLEL_PROCESSES=$OPTARG
   if [[ $PEXEC_PARALLEL_PROCESSES -lt 2 ]];
     then
       echo " Number of parallel processes for pexec must be > 1. Use -c 0 to run serially."
       exit 1
     fi
   ;;
      k) #number of modalities used to construct the template (default = 1)
   NUMBEROFMODALITIES=$OPTARG
   ;;
      w) #modality weights (default = 1)
   MODALITYWEIGHTSTRING=$OPTARG
   ;;
      m) #max iterations other than default
   MAXITERATIONS=$OPTARG
   ;;
      n) #apply bias field correction
   N4CORRECT=$OPTARG
   ;;
      o) #output name prefix
   OUTPUTNAME=$OPTARG
   TEMPLATENAME=${OUTPUTNAME}template
   ;;
      p) #Script prepend
   SCRIPTPREPEND+=( "$OPTARG" )
   ;;
      s) #similarity model
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
      x) # XGrid options (deprecated)
   XGRIDOPTS=$OPTARG
   ;;
      y) # update with full affine, 0 for no rigid (default = 1)
   AFFINE_UPDATE_FULL=$OPTARG
   ;;
      z) #initialization template
   REGTEMPLATES[${#REGTEMPLATES[@]}]=$OPTARG
   ;;
      \?) # getopts issues an error message
      Usage >&2
      exit 1
  esac
done

# Set up multi-threading
if [[ $NUMBER_OF_THREADS -lt 0 ]];
  then
    # Number of threads not set on the command line, try ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
    if [[ -n ${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS:-} ]];
      then
	    NUMBER_OF_THREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}
        echo "Number of threads per process not set on the command line, using environment value."
        echo "Environment has ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}"
      else
        if [[ $DOQSUB -eq 0 ]];
          then
            NUMBER_OF_THREADS=$(( cpu_count < 8 ? cpu_count : 8 ))
            echo "Number of threads not set on the command line or via ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS."
            echo "Using min(cpu_count, 8) = ${NUMBER_OF_THREADS} threads per process for serial execution"
          else
            NUMBER_OF_THREADS=1
            echo "Number of threads not set on the command line or via ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS."
            echo "Defaulting to 1 thread per parallel job"
          fi
      fi
  else
    echo "Number of threads per process set to $NUMBER_OF_THREADS on the command line"
  fi

if [[ $NUMBER_OF_THREADS -eq 0 ]];
  then
    NUMBER_OF_THREADS=$(( cpu_count < 8 ? cpu_count : 8 ))
    echo "Number of threads set to 0. Using min(cpu_count, 8) = ${NUMBER_OF_THREADS} threads per process"
  fi

# This sets an appropriate number of threads for serial / pexec jobs and also things that run locally
# like AverageImages, ImageMath, ImageSetStatistics, MeasureMinMaxMean, etc
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NUMBER_OF_THREADS

if [[ $TDIM -eq 4 && $TRANSFORMATIONTYPE == "GR" ]]; then
  # Use stronger regularization for 4D mapping, where deformations should be small and local.
  echo "Using GR_Constrained for 4D template construction"
  TRANSFORMATIONTYPE="GR_Constrained"
fi

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
    OUTPUT_DIR=$(dirname "$OUTPUTNAME")
  fi

if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p "$OUTPUT_DIR"
  fi

# Intermediate template output. Keep the template for each iteration and also the average warp if defined.
# Useful for debugging and monitoring convergence
intermediateTemplateDir=${OUTPUT_DIR}/intermediateTemplates
mkdir -p "$intermediateTemplateDir"

if [[ $DOQSUB -eq 1 || $DOQSUB -eq 4 ]];
  then
    qq=`which  qsub`
    if [[  ${#qq} -lt 1 ]];
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
  IFS=x read -r -a MODALITYWEIGHTS <<< "$MODALITYWEIGHTSTRING"
  if [[ ${#MODALITYWEIGHTS[@]} -ne $NUMBEROFMODALITIES ]];
    then
      echo "The number of weights (specified e.g. -w 1x1x1) does not match the number of specified modalities (see -k option)";
      exit
    fi
  fi


# Creating the file list of images to make a template from.
# Shiftsize is calculated because a variable amount of arguments can be used on the command line.
# The shiftsize variable will give the correct number of arguments to skip. Issuing shift $shiftsize will
# result in skipping that number of arguments on the command line, so that only the input images remain.
shiftsize=$(($OPTIND - 1))
shift "$shiftsize"
NINFILES=$#
IMAGESETARRAY=( "$@" )
IMAGESETVARIABLE=${IMAGESETARRAY[*]}

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

if [[ ${NINFILES} -eq 0 ]];
  then
    echo "Please provide at least 2 filenames for the template."
    echo "Use `basename $0` -h for help"
    exit 1
  elif [[ ${NINFILES} -eq 1 ]];
    then
    IMAGESETVARIABLE=${IMAGESETARRAY[0]}
    extension=${IMAGESETVARIABLE##*.}
    if [[ $extension = 'csv' || $extension = 'txt' ]];
        then
        IMAGESFILE=$IMAGESETVARIABLE
        IMAGESETARRAY=()
        IMAGECOUNT=0
        while IFS= read -r line || [[ -n $line ]]
            do
            line=${line%$'\r'} # remove carriage return from python / windows line-endings
            IFS=',' read -r -a files <<< "$line"
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
         done < "$IMAGESFILE"
    else
        if [[ ${TDIM} -eq 0 ]];
            then
            echo "Please provide at least 2 filenames for the template."
            echo "If building a template from a 4D image, use -d 4"
            echo "Use `basename $0` -h for help"
            exit 1
        fi
        if ! size_string=$("$PRINTHEADER" "$IMAGESETVARIABLE" 2); then
            printf 'Could not read image dimensions from %s\n' "$IMAGESETVARIABLE" >&2
            exit 1
        fi
        IFS='x' read -r -a image_size <<< "$size_string"

        if (( ${#image_size[@]} != TDIM )); then
            printf 'Expected a %dD time-series image, but %s has %d dimensions\n' \
                "$TDIM" "$IMAGESETVARIABLE" "${#image_size[@]}" >&2
            exit 1
        fi

        range=${image_size[$((TDIM - 1))]}
        if [[ ! $range =~ ^[1-9][0-9]*$ ]]; then
            printf 'Invalid time dimension reported for %s: %s\n' "$IMAGESETVARIABLE" "$range" >&2
            exit 1
        fi

        if [[ ${range} -eq 1 ]];
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
            echo " Creating template from 4D input. "
            echo "--------------------------------------------------------------------------------------"

             # Stage a private copy of the input. Selected volumes are written
             # to OUTPUT_DIR and remain there after template construction.
             tmpdir=${currentdir}/tmp_${RANDOM}_${RANDOM}_${RANDOM}_$$
             (umask 077 && mkdir "${tmpdir}") || {
                 echo "Could not create temporary directory! Exiting." 1>&2
                 exit 1
                 }

             staged_input="${tmpdir}/$(basename "$IMAGESETVARIABLE")"
             if ! cp "$IMAGESETVARIABLE" "$staged_input"; then
                 cleanup_4d_tempdir
                 echo "Could not stage 4D input image. Exiting." 1>&2
                 exit 1
             fi

             nfmribins=16
             selected_indices=()

             if (( range <= nfmribins )); then
                 # Use every time point when the series contains at most 16.
                 for ((i = 0; i < range; i++)); do
                     selected_indices+=( "$i" )
                 done
             elif (( range < 32 )); then
                 # Select 16 unique time points uniformly without replacement.
                 candidate_indices=()
                 for ((i = 0; i < range; i++)); do
                     candidate_indices+=( "$i" )
                 done
                 for ((i = 0; i < nfmribins; i++)); do
                     random_index_below $((range - i))
                     swap_index=$((i + RANDOM_INDEX))
                     number=${candidate_indices[$swap_index]}
                     candidate_indices[$swap_index]=${candidate_indices[$i]}
                     candidate_indices[$i]=$number
                     selected_indices+=( "$number" )
                 done
             else
                 # Divide the complete time series into 16 non-overlapping bins
                 # and select every member of each bin with equal probability.
                 for ((i = 0; i < nfmribins; i++)); do
                     bin_start=$((i * range / nfmribins))
                     bin_end=$(((i + 1) * range / nfmribins - 1))
                     random_index_below $((bin_end - bin_start + 1))
                     selected_indices+=( "$((bin_start + RANDOM_INDEX))" )
                 done
             fi

             IMAGESETARRAY=()
             selected_volume_dir=$(cd "$OUTPUT_DIR" && pwd -P)
             for number in "${selected_indices[@]}"; do
                 printf -v selected_volume '%s/time_index_%04d.nii.gz' "$selected_volume_dir" "$number"
                 printf 'Selecting time point %d as %s\n' "$number" "$selected_volume"
                 if ! ImageMath "$TDIM" "$selected_volume" ExtractSlice "$staged_input" "$number"; then
                     cleanup_4d_tempdir
                     echo "Could not extract time point $number. Exiting." 1>&2
                     exit 1
                 fi
                 IMAGESETARRAY+=( "$selected_volume" )
             done

             cleanup_4d_tempdir
             unset candidate_indices selected_indices selected_volume_dir staged_input selected_volume
             # Continue in the invocation directory with absolute input paths.
             cd "$currentdir"
        fi
    fi
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
        echo "$IMAGEMETRICSET"
    done
    echo "--------------------------------------------------------------------------------------"
fi

# Useful to check the right number of images exist for various ops
IMAGESPERMODALITY=$(( ${#IMAGESETARRAY[@]} / ${NUMBEROFMODALITIES} ))

# check for initial template images
for (( i = 0; i < $NUMBEROFMODALITIES; i++ ))
    do
    setCurrentImageSet $i

    if [[ ${#REGTEMPLATES[@]} -gt 0 && -n "${REGTEMPLATES[$i]:-}" ]];
      then
        if [[ ! -r "${REGTEMPLATES[$i]}" ]];
          then
            echo "Initial template {REGTEMPLATES[$i]} cannot be read"
            exit 1
          fi
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Initial template $i found.  This will be used for guiding the registration. use : ${REGTEMPLATES[$i]} and ${TEMPLATES[$i]} "
        echo "--------------------------------------------------------------------------------------"
     # now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.
        cp "${REGTEMPLATES[$i]}" "${TEMPLATES[$i]}"
    else
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Creating template ${TEMPLATES[$i]} from a population average image from the inputs."
        echo "   ${CURRENTIMAGESET[*]}"
        echo "--------------------------------------------------------------------------------------"
        # Normalized mean, no sharpening
        # This forces a call to AverageImages, which resizes images to match the largest input
        summarizeimageset "$DIM" "${TEMPLATES[$i]}" 1 0 "${CURRENTIMAGESET[@]}"
        # Quickly align COM of input images to average, and then recompute average
        IMAGECOMSET=()
        for (( j = 0; j < ${#CURRENTIMAGESET[@]}; j+=1 ))
          do
            IMGbase=`basename "${CURRENTIMAGESET[$j]}"`
            BASENAME=$(get_image_stem "$IMGbase")
            COM="${OUTPUT_DIR}/initialCOM${i}_${j}_${IMGbase}"
            COMTRANSFORM="${OUTPUT_DIR}/initialCOM${i}_${j}_${BASENAME}.mat"
            antsAI -d "${DIM}" --convergence 0 --verbose 1 -m "Mattes[${TEMPLATES[$i]},${CURRENTIMAGESET[$j]},32,None]" -o "${COMTRANSFORM}" -t AlignCentersOfMass
            antsApplyTransforms -d "${DIM}" -r "${TEMPLATES[$i]}" -i "${CURRENTIMAGESET[$j]}" -t "${COMTRANSFORM}" -o "${COM}" --verbose
            rm -f "$COMTRANSFORM"
            IMAGECOMSET[${#IMAGECOMSET[@]}]=$COM
          done
          # Now safe to let user control stat method
          summarizeimageset "$DIM" "${TEMPLATES[$i]}" "${STATSMETHOD}" 0 "${IMAGECOMSET[@]}"
          # Clean up
          rm -f "${IMAGECOMSET[@]}"
    fi

    if [[ ! -s ${TEMPLATES[$i]} ]];
        then
        echo "Your template : ${TEMPLATES[$i]} was not created.  This indicates trouble!  You may want to check correctness of your input parameters. exiting."
        exit
    fi

    # Back up template
    intermediateTemplateBase=`basename "${TEMPLATES[$i]}"`
    cp "${TEMPLATES[$i]}" "${intermediateTemplateDir}/initial_${intermediateTemplateBase}"

done

# remove old job bash scripts
outdir=`dirname "${TEMPLATES[0]}"`
if [[ ${#outdir} -eq 0 ]];
    then
    outdir=`pwd`
fi
rm -f "${outdir}"/job*.sh

##########################################################################
#
# perform rigid body registration if requested
#
##########################################################################
if [[ "$RIGID" -eq 1 ]];
    then
    count=0
    jobIDs=()
    for (( i = 0; i < ${#IMAGESETARRAY[@]}; i+=$NUMBEROFMODALITIES ))
        do
        IMAGEMETRICSET=""
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
            do
            k=0
            k=$((i+j))
            printf -v template_j_q '%q' "${TEMPLATES[$j]}"
            printf -v input_k_q '%q' "${IMAGESETARRAY[$k]}"
            IMAGEMETRICSET="$IMAGEMETRICSET -m MI[ ${template_j_q},${input_k_q},${MODALITYWEIGHTS[$j]},32 ]"
        done

        qscript="${outdir}/job_${count}_qsub.sh"
        rm -f "$qscript"

        if [[ $DOQSUB -eq 5 ]];
            then
            # SLURM job scripts must start with a shebang
            printf '%s\n' '#!/bin/sh' > "$qscript"
            fi

        write_script_prepend >> "$qscript"
        printf 'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%s\n' "$NUMBER_OF_THREADS" >> "$qscript"

        IMGbase=`basename "${IMAGESETARRAY[$i]}"`
        BASENAME=$(get_image_stem "$IMGbase")
        RIGID="${outdir}/rigid${i}_0_${IMGbase}"
        printf -v rigid_q '%q' "$RIGID"

        exe="$ANTS $DIM $IMAGEMETRICSET -o $rigid_q -i 0 $LINEARTRANSFORMPARAMS $RIGIDTYPE"

        printf '%s\n' "$exe" >> "$qscript"

        exe2='';
        pexe2='';
        printf -v metriclog_q '%q' "${outdir}/job_${count}_metriclog.txt"
        pexe=" $exe > ${metriclog_q} "
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
            do
            k=0
            k=$((i+j))
            IMGbase=`basename "${IMAGESETARRAY[$k]}"`
            BASENAME=$(get_image_stem "$IMGbase")
            RIGID="${outdir}/rigid${i}_${j}_${IMGbase}"
            IMGbaseBASE=`basename "${IMAGESETARRAY[$i]}"`
            BASENAMEBASE=$(get_image_stem "$IMGbaseBASE")
            printf -v input_k_q '%q' "${IMAGESETARRAY[$k]}"
            printf -v rigid_q '%q' "$RIGID"
            printf -v affine_q '%q' "${outdir}/rigid${i}_0_${BASENAMEBASE}Affine.txt"
            printf -v template_j_q '%q' "${TEMPLATES[$j]}"
            exe2+=" ${WARP} $DIM ${input_k_q} ${rigid_q} ${affine_q} -R ${template_j_q}"$'\n'
            pexe2="$exe2 ${WARP} $DIM ${input_k_q} ${rigid_q} ${affine_q} -R ${template_j_q} >> ${metriclog_q}"$'\n'
        done

        printf '%s' "$exe2" >> "$qscript";

        if [[ $DOQSUB -eq 1 ]];
            then
            id=`qsub -cwd -S /bin/bash -N antsBuildTemplate_rigid $QSUBOPTS "$qscript" | awk '{print $3}'`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
            then
            id=`qsub -N antsrigid -l mem=${MEMORY} -l walltime=${WALLTIME} $QSUBOPTS "$qscript" | awk '{print $1}'`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 2 ]];
            then
            # Send pexe and exe2 to same job file so that they execute in series
            printf '%s\n' "$pexe" >> "${outdir}/job${count}_r.sh"
            printf '%s' "$pexe2" >> "${outdir}/job${count}_r.sh"
        elif [[ $DOQSUB -eq 3 ]];
            then
            id=`xgrid $XGRIDOPTS -job submit /bin/bash "$qscript" | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
            #echo "xgrid $XGRIDOPTS -job submit /bin/bash $qscript"
            jobIDs+=("$id")
        elif [[ $DOQSUB -eq 5 ]];
            then
            id=`sbatch --job-name=antsrigid --nodes=1 --ntasks=1 --cpus-per-task=$NUMBER_OF_THREADS --time=${WALLTIME} --mem=${MEMORY} $QSUBOPTS "$qscript" | rev | cut -f1 -d\ | rev`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
            then
             # execute jobs in series
             bash "$qscript"
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
        # Returns 1 if there are errors
        if ! waitForSGEQJobs.pl 1 60 "${jobIDs[@]}";
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
        # Returns 1 if there are errors
        if ! waitForPBSQJobs.pl 1 60 "${jobIDs[@]}";
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
        echo " Starting ANTS rigid registration with ${PEXEC_PARALLEL_PROCESSES} parallel processes. "
        echo " Progress can be viewed in ${outdir}/job*_metriclog.txt"
        echo "--------------------------------------------------------------------------------------"
        jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
        chmod +x "${outdir}"/job*_r.sh
        "$PEXEC" -j "${PEXEC_PARALLEL_PROCESSES}" "sh" "${outdir}"/job*_r.sh
    fi
    if [[ $DOQSUB -eq 3 ]];
        then
        # Run jobs on XGrid and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS rigid registration on XGrid cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
        # Returns 1 if there are errors
        if ! waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 30 "${jobIDs[@]}";
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
        # Returns 1 if there are errors
        if ! waitForSlurmJobs.pl 1 60 "${jobIDs[@]}";
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
            IMGbase=`basename "${IMAGESETARRAY[$i]}"`
            BASENAME=$(get_image_stem "$IMGbase")
            RIGID="${outdir}/rigid${k}_${j}_${IMGbase}"

            IMAGERIGIDSET[${#IMAGERIGIDSET[@]}]=$RIGID
        done
        echo
        echo  "AverageImages $DIM ${TEMPLATES[$j]} 2 ${IMAGERIGIDSET[*]}"

      # Don't sharpen after rigid alignment
      summarizeimageset "$DIM" "${TEMPLATES[$j]}" "${STATSMETHOD}" 0 "${IMAGERIGIDSET[@]}"
      intermediateTemplateBase=`basename "${TEMPLATES[$j]}"`
      cp "${TEMPLATES[$j]}" "${intermediateTemplateDir}/initialRigid_${intermediateTemplateBase}"

    done

    # cleanup and save output in seperate folder
    if [[ BACKUP_EACH_ITERATION -eq 1 ]];
      then

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from rigid iteration"
        echo "--------------------------------------------------------------------------------------"

        mkdir "${outdir}/rigid"
        mv "${outdir}"/rigid*.nii.gz "${outdir}"/*Affine.txt "${outdir}"/*GenericAffine.mat "${outdir}/rigid/"
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
          then
            mv "${outdir}"/antsBuildTemplate_rigid* "${outdir}/rigid/"
            # Remove qsub scripts
            rm -f "${outdir}/job_${count}_qsub.sh"
        elif [[ $DOQSUB -eq 4 ]];
          then
            mv "${outdir}"/antsrigid* "${outdir}/rigid/"
            # Remove qsub scripts
            rm -f "${outdir}/job_${count}_qsub.sh"
        elif [[ $DOQSUB -eq 2 ]];
          then
            mv "${outdir}"/job*.txt "${outdir}/rigid/"
        elif [[ $DOQSUB -eq 3 ]];
          then
            rm -f "${outdir}"/job_*_qsub.sh
        elif [[ $DOQSUB -eq 5 ]];
          then
            mv "${outdir}"/slurm-*.out "${outdir}/rigid/"
            mv "${outdir}"/job*.txt "${outdir}/rigid/"

            # Remove qsub scripts
            rm -f "${outdir}/job_${count}_qsub.sh"
        fi
      else
        rm -f "${outdir}"/rigid*.* "${outdir}"/job*.txt "${outdir}"/slurm-*.out
      fi
fi # endif RIGID

##########################################################################
#
# begin main level
#
##########################################################################

IFS=x read -r -a ITERATLEVEL <<< "$MAXITERATIONS"
NUMLEVELS=${#ITERATLEVEL[@]}
#
# debugging only
#echo $ITERATLEVEL
#echo $NUMLEVELS
#echo ${ITERATIONLIMIT}
#
echo
echo "--------------------------------------------------------------------------------------"
echo " Start to build templates: ${TEMPLATES[*]}"
echo "--------------------------------------------------------------------------------------"

TRANSFORMATION=''
REGULARIZATION=''
if [[ "${TRANSFORMATIONTYPE}" == "EL" ]];
    then
    # Mapping Parameters
    TRANSFORMATION="Elast[ 1 ]"
    REGULARIZATION="Gauss[ 3,0.5 ]"
    # Gauss[3,x ] is usually the best option.    x is usually 0 for SyN --- if you want to reduce flexibility/increase mapping smoothness, the set x > 0.
    # We did a large scale evaluation of SyN gradient parameters in normal brains and found 0.25 => 0.5 to perform best when
    # combined with default Gauss[3,0 ] regularization.    You would increase the gradient step in some cases, though, to make
    # the registration converge faster --- though oscillations occur if the step is too high and other instability might happen too.
elif [[ "${TRANSFORMATIONTYPE}" == "S2" ]];
    then
    # Mapping Parameters for the LDDMM style SyN --- the params are SyN[ GradientStepLength,NTimeDiscretizationPoints,IntegrationTimeStep]
    # increasing IntegrationTimeStep increases accuracy in the diffeomorphism integration and takes more computation time.
    # NTimeDiscretizationPoints is set to 2 here
    TRANSFORMATION="SyN[ 1,2,0.05 ]"
    REGULARIZATION="Gauss[ 3,0. ]"
elif [[ "${TRANSFORMATIONTYPE}" == "SY" ]];
    then
    # Mapping Parameters for the LDDMM style SyN --- the params are SyN[ GradientStepLength,NTimeDiscretizationPoints,IntegrationTimeStep]
    # increasing IntegrationTimeStep increases accuracy in the diffeomorphism integration and takes more computation time.
    # NTimeDiscretizationPoints is the number of spatial indices in the time dimension (the 4th dim when doing 3D registration)
    # increasing NTimeDiscretizationPoints increases flexibility and takes more computation time.
    # the --geodesic option enables either 1 asymmetric gradient estimation or 2 symmetric gradient estimation (the default here )
    TRANSFORMATION="SyN[ 1,2,0.05 ] --geodesic 2"
    REGULARIZATION="Gauss[ 3,0. ]"
elif [[ "${TRANSFORMATIONTYPE}" == "LDDMM" ]];
   then
   # Mapping Parameters for the LDDMM style SyN --- the params are SyN[ GradientStepLength,NTimeDiscretizationPoints,IntegrationTimeStep]
   # increasing IntegrationTimeStep increases accuracy in the diffeomorphism integration and takes more computation time.
   # NTimeDiscretizationPoints is the number of spatial indices in the time dimension (the 4th dim when doing 3D registration)
   # increasing NTimeDiscretizationPoints increases flexibility and takes more computation time.
   # the --geodesic option enables either 1 asymmetric gradient estimation or 2 symmetric gradient estimation (the default here )
   TRANSFORMATION="SyN[1,2,0.05 ] --geodesic 1"
   REGULARIZATION="Gauss[ 3,0. ]"
elif [[ "${TRANSFORMATIONTYPE}" == "GR" ]];
    then
    # Mapping Parameters for the greedy gradient descent (fast) version of SyN -- only needs GradientStepLength
    TRANSFORMATION="SyN[ 0.25 ]"
    REGULARIZATION="Gauss[ 3,0 ]"
elif [[ "${TRANSFORMATIONTYPE}" == "GR_Constrained" ]];
    then
    # Mapping Parameters for the greedy gradient descent (fast) version of SyN -- only needs GradientStepLength
    TRANSFORMATION="SyN[ 0.25 ]"
    REGULARIZATION="Gauss[ 3,0.5 ]"

elif [[ "${TRANSFORMATIONTYPE}" == "EX" ]];
    then
    # Mapping Parameters
    TRANSFORMATION="Exp[ 0.5,10 ]"
    REGULARIZATION="Gauss[ 3,0.5 ]"
elif [[ "${TRANSFORMATIONTYPE}" == "DD" ]];
    then
    # Mapping Parameters for diffemorphic demons style optimization Exp[GradientStepLength,NTimePointsInIntegration]
    #  NTimePointsInIntegration controls the number of compositions in the transformation update , see the DD paper
    TRANSFORMATION="GreedyExp[ 0.5,10 ]"
    REGULARIZATION="Gauss[ 3,0.5 ]"
else
    echo "Invalid transformation metric. Use EL, SY, S2, GR , DD or EX or type bash `basename $0` -h."
    exit 1
fi

reportMappingParameters
#

i=0
while [[ $i -lt ${ITERATIONLIMIT} ]];
    do
    itdisplay=$((i+1))
    rm -f "${OUTPUTNAME}"*Warp*.nii*
    rm -f "${outdir}"/job*.sh
    # Used to save time by only running coarse registration for the first couple of iterations
    # This may also help convergence, but because there's no way to turn it off, it makes it harder
    # to refine templates with multiple calls to this script.
    # If you uncomment this, replace MAXITERATIONS with ITERATIONS in the call to ants below
    #
    # # For the first couple of iterations, use high-level registration only
    # # eg if MAXITERATIONS=30x90x20, then for iteration 0, do 30x0x0
    # # for iteration 1 do 30x90x0, then do 30x90x20 on subsequent iterations
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
    jobIDs=()
    # Reinitialize count to 0
    count=0
    # Submit registration of each input to volume template to SGE or run locally.

    for (( j = 0; j < ${#IMAGESETARRAY[@]}; j+=$NUMBEROFMODALITIES ))
        do
        IMAGEMETRICSET=''
        exe=''
        warpexe=''
        pexe=''
        warppexe=''

        for (( k = 0; k < $NUMBEROFMODALITIES; k++ ))
          do
            l=0
            l=$((j+k))

            if [[ "${METRICTYPE[$k]}" == "PR" ]];
                then
                # Mapping Parameters
                METRIC="PR[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},4 ]"
            elif [[ "${METRICTYPE[$k]}" == "CC" ]];
                then
                # Mapping Parameters
                METRIC="CC[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},5 ]"
            elif [[ "${METRICTYPE[$k]}" == "MI" ]];
                then
                # Mapping Parameters
                METRIC="MI[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},32 ]"
            elif [[ "${METRICTYPE[$k]}" == "MSQ" ]];
                then
                # Mapping Parameters
                METRIC="MSQ[ "
                METRICPARAMS="${MODALITYWEIGHTS[$k]},0 ]"
            else
                echo "Invalid similarity metric. Use CC, MI, MSQ or PR or type bash `basename $0` -h."
                exit 1
            fi
            TEMPLATEbase=`basename "${TEMPLATES[$k]}"`
            indir=`dirname "${IMAGESETARRAY[$j]}"`
            if [[ ${#indir} -eq 0 ]];
                then
                indir=`pwd`
            fi
            IMGbase=`basename "${IMAGESETARRAY[$l]}"`
            IMGstem=$(get_image_stem "$IMGbase")
            OUTFN=${OUTPUTNAME}template-modality${k}-${IMGstem}
            OUTFN=`basename "${OUTFN}"`
            OUTFN="${OUTFN}${l}"
            DEFORMED="${outdir}/${OUTFN}${l}WarpedToTemplate.nii.gz"

            IMGbase=`basename "${IMAGESETARRAY[$j]}"`
            IMGstem=$(get_image_stem "$IMGbase")
            OUTWARPFN=${OUTPUTNAME}${IMGstem}
            OUTWARPFN=`basename "${OUTWARPFN}"`
            OUTWARPFN="${OUTWARPFN}${j}"

            printf -v template_k_q '%q' "${TEMPLATES[$k]}"
            printf -v input_l_q '%q' "${IMAGESETARRAY[$l]}"
            printf -v deformed_q '%q' "$DEFORMED"
            printf -v warp_q '%q' "${outdir}/${OUTWARPFN}Warp.nii.gz"
            printf -v affine_q '%q' "${outdir}/${OUTWARPFN}Affine.txt"
            printf -v metriclog_q '%q' "${outdir}/job_${count}_metriclog.txt"

            if [[ $N4CORRECT -eq 1 ]];
              then
                REPAIRED="${outdir}/${OUTFN}Repaired.nii.gz"
                printf -v repaired_q '%q' "$REPAIRED"
                if [[ ! -s ${REPAIRED} ]]; then
                  exe=" $exe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${input_l_q} -o ${repaired_q} -r 0 -s 2"$'\n'
                  pexe=" $pexe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${input_l_q} -o ${repaired_q} -r 0 -s 2  >> ${metriclog_q}"$'\n'
                fi
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${template_k_q},${repaired_q},${METRICPARAMS}"
                warpexe=" $warpexe ${WARP} ${DIM} ${repaired_q} ${deformed_q} -R ${template_k_q} ${warp_q} ${affine_q}"$'\n'
                warppexe=" $warppexe ${WARP} ${DIM} ${repaired_q} ${deformed_q} -R ${template_k_q} ${warp_q} ${affine_q} >> ${metriclog_q}"$'\n'
              else
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${template_k_q},${input_l_q},${METRICPARAMS}";
                warpexe=" $warpexe ${WARP} ${DIM} ${input_l_q} ${deformed_q} -R ${template_k_q} ${warp_q} ${affine_q}"$'\n'
                warppexe=" $warppexe ${WARP} ${DIM} ${input_l_q} ${deformed_q} -R ${template_k_q} ${warp_q} ${affine_q} >> ${metriclog_q}"$'\n'
              fi

        done

        IMGbase=`basename "${IMAGESETARRAY[$j]}"`
        IMGstem=$(get_image_stem "$IMGbase")
        OUTWARPFN=${OUTPUTNAME}${IMGstem}
        OUTWARPFN=`basename "${OUTWARPFN}${j}"`
        printf -v output_prefix_q '%q' "${outdir}/${OUTWARPFN}"

        LINEARTRANSFORMPARAMS="--number-of-affine-iterations 10000x10000x1000 --MI-option 32x16000"

        exe="$exe $ANTS ${DIM} $IMAGEMETRICSET -i ${MAXITERATIONS} -t ${TRANSFORMATION} -r $REGULARIZATION -o ${output_prefix_q} $LINEARTRANSFORMPARAMS"$'\n'
        exe="$exe $warpexe"

        pexe="$pexe $ANTS ${DIM} $IMAGEMETRICSET -i ${MAXITERATIONS} -t ${TRANSFORMATION} -r $REGULARIZATION -o ${output_prefix_q} $LINEARTRANSFORMPARAMS >> ${metriclog_q}"$'\n'
        pexe="$pexe $warppexe"

        qscript="${outdir}/job_${count}_${i}.sh"

        printf '%s' "$exe" >> "${outdir}/job_${count}_${i}_metriclog.txt"
        # 6 submit to SGE (DOQSUB=1), PBS (DOQSUB=4), PEXEC (DOQSUB=2), XGrid (DOQSUB=3) or else run locally (DOQSUB=0)
        if [[ $DOQSUB -eq 1 ]];
            then
            write_script_prepend > "$qscript"
            printf 'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%s\n' "$NUMBER_OF_THREADS" >> "$qscript"
            printf '%s' "$exe" >> "$qscript"
            id=`qsub -cwd -N antsBuildTemplate_deformable_${i} -S /bin/bash $QSUBOPTS "$qscript" | awk '{print $3}'`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
            then
            write_script_prepend > "$qscript"
            printf 'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%s\n' "$NUMBER_OF_THREADS" >> "$qscript"
            printf '%s' "$exe" >> "$qscript"
            id=`qsub -N antsdef${i} -l walltime=${WALLTIME} -l mem=${MEMORY} $QSUBOPTS "$qscript" | awk '{print $1}'`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 2 ]];
            then
            printf '%s' "$pexe" >> "${outdir}/job${count}_r.sh"
        elif [[ $DOQSUB -eq 3 ]];
            then
            write_script_prepend > "$qscript"
            printf 'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%s\n' "$NUMBER_OF_THREADS" >> "$qscript"
            printf '%s' "$exe" >> "$qscript"
            id=`xgrid $XGRIDOPTS -job submit /bin/bash "$qscript" | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
            jobIDs+=("$id")
        elif [[ $DOQSUB -eq 5 ]];
            then
            printf '%s\n' '#!/bin/sh' > "$qscript"
            write_script_prepend >> "$qscript"
            printf 'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%s\n' "$NUMBER_OF_THREADS" >> "$qscript"
            printf '%s' "$exe" >> "$qscript"
            id=`sbatch --job-name=antsdef${i} --nodes=1 --ntasks=1 --cpus-per-task=$NUMBER_OF_THREADS --time=${WALLTIME} --mem=${MEMORY} $QSUBOPTS "$qscript" | rev | cut -f1 -d\ | rev`
            jobIDs+=("$id")
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
            then
            printf '%s' "$exe" > "$qscript"
            bash "$qscript"
        fi

        # counter updated, but not directly used in this loop
        count=`expr $count + 1`;
    #  echo " submitting job number $count " # for debugging only
    done
    # SGE wait for script to finish
    if [[ $DOQSUB -eq 1 ]];
        then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on SGE cluster. Iteration: $itdisplay of $ITERATIONLIMIT"
        echo "--------------------------------------------------------------------------------------"
        # now wait for the stuff to finish - this will take a while so poll queue every 10 mins
        if ! waitForSGEQJobs.pl 1 600 "${jobIDs[@]}";
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
        if ! waitForPBSQJobs.pl 1 600 "${jobIDs[@]}";
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
        echo " Starting ANTS registration using ${PEXEC_PARALLEL_PROCESSES} parallel processes. Iteration: $itdisplay of $ITERATIONLIMIT"
        echo " Progress can be viewed in job*_${i}_metriclog.txt"
        echo "--------------------------------------------------------------------------------------"
        jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
        chmod +x "${outdir}"/job*.sh
        "$PEXEC" -j "${PEXEC_PARALLEL_PROCESSES}" sh "${outdir}"/job*.sh
    fi
    if [[ $DOQSUB -eq 3 ]];
        then
        # Run jobs on XGrid and wait to finish
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Starting ANTS registration on XGrid cluster. Submitted $count jobs "
        echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. This is slow, so poll less often
        # Returns 1 if there are errors
        if ! waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 300 "${jobIDs[@]}";
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
        if ! waitForSlurmJobs.pl 1 600 "${jobIDs[@]}";
            then
            echo "SLURM submission failed - jobs went into error state"
            exit 1;
        fi
    fi

    for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
      do
        shapeupdatetotemplate "${DIM}" "${TEMPLATES[$j]}" "${TEMPLATENAME}" "${OUTPUTNAME}" "${GRADIENTSTEP}" "${STATSMETHOD}" "${SHARPENMETHOD}" "${j}"
        intermediateTemplateBase=`basename "${TEMPLATES[$j]}"`
        cp "${TEMPLATES[$j]}" "${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_${intermediateTemplateBase}"
      done

    if [[ -f "${TEMPLATENAME}0warp.nii.gz" ]]
      then
        cp "${TEMPLATENAME}0warp.nii.gz" "${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_shapeUpdateWarp.nii.gz"
      fi

    if [[ BACKUP_EACH_ITERATION -eq 1 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from iteration $itdisplay"
        echo "--------------------------------------------------------------------------------------"
        mkdir "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
        cp "${TEMPLATENAME}${j}warplog.txt" "${outdir}"/*.cfg "${outdir}"/*Affine.txt "${OUTPUTNAME}"*.nii.gz "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}/"
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
            then
            mv "${outdir}"/antsBuildTemplate_deformable_* "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
        elif [[ $DOQSUB -eq 4 ]];
            then
            mv "${outdir}"/antsdef* "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
        elif [[ $DOQSUB -eq 2 ]];
            then
            mv "${outdir}"/job*.txt "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
        elif [[ $DOQSUB -eq 3 ]];
            then
            rm -f "${outdir}"/job_*.sh
        elif [[ $DOQSUB -eq 5 ]];
            then
            mv "${outdir}"/slurm-*.out "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
            mv "${outdir}"/job*.txt "${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}"
        fi
      else
        rm -f "${outdir}"/job*.txt "${outdir}"/slurm-*.out
      fi
    i=$(( i + 1 ))
done
# end main loop
rm -f job*.sh
time_end=`date +%s`
time_elapsed=$((time_end - time_start))
echo
echo "--------------------------------------------------------------------------------------"
echo " Done creating: ${TEMPLATES[*]}"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
