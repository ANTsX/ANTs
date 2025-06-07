#!/bin/bash

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

fle_error=0
for FLE in $ANTS $WARP $N4 $PEXEC $SGE $XGRID $PBS $SLURM
  do
    if ! command -v $FLE &> /dev/null
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
   ImageDimension: 4 (for template generation of time-series data)

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.  Optionally, one can specify a .csv or .txt file where each
          line is the location of the input image.  One can also specify more than
          one file for each image for multi-modal template construction (e.g. t1 and t2).
          For the multi-modal case, the templates will be consecutively numbered (e.g.
          ${OUTPUTPREFIX}template0.nii.gz, ${OUTPUTPREFIX}template1.nii.gz, ...).

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
          2 == use PEXEC (localhost), 3 == Apple XGrid, 4 == PBS qsub, 5 == SLURM

     -g:  Gradient step size (default 0.25) -- smaller in magnitude results in
          more cautious steps. Use smaller steps to refine template details.
          0.25 is an upper (aggressive) limit for this parameter.

     -i:  Iteration limit (default 4) -- iterations of the template construction (Iteration limit)*NumImages registrations.

     -j:  Number of cpu cores to use (default 2; -- requires "-c 2")

     -k:  Number of modalities used to construct the template (default 1)

     -w:  Modality weights used in the similarity metric (default = 1) --- specified as e.g. 1x0.5x0.75

     -m:  Max-iterations in each registration

     -n:  N4BiasFieldCorrection of moving image (default 1) -- 0 == off, 1 == on

     -p:  Commands to prepend to job scripts (e.g., change into appropriate directory, set paths, etc)

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

    imagelist=(`ls ${outputname}template-modality${whichtemplate}-*WarpedToTemplate.nii.gz`)
    if [[ ${#imagelist[@]} -ne ${IMAGESPERMODALITY} ]]
      then
        echo "ERROR shapeupdatedtotemplate - imagelist length is ${#imagelist[@]}, expected ${IMAGESPERMODALITY}"
        exit 1
      fi

    summarizeimageset ${dim} ${template} ${summarizemethod} ${sharpenmethod} ${imagelist[@]}

    if [[ $whichtemplate -eq 0 ]] ;
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---voxel-wise averaging of the inverse warp fields (from subject to template)"
        echo "   AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`"
        echo "--------------------------------------------------------------------------------------"

        AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---scale the averaged inverse warp field by the gradient step"
        echo "   MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz"
        echo "--------------------------------------------------------------------------------------"

        MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---average the affine transforms (template <-> subject)"
        echo "                      ---transform the inverse field by the resulting average affine transform"
        echo "   ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0Affine.txt ${outputname}*Affine.txt"
        echo "   WarpImageMultiTransform ${dim} ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -i  ${templatename}0Affine.txt -R ${template}"
        echo "--------------------------------------------------------------------------------------"

        ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0Affine.txt ${outputname}*Affine.txt
        WarpImageMultiTransform ${dim} ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -i  ${templatename}0Affine.txt -R ${template}

        MeasureMinMaxMean ${dim} ${templatename}0warp.nii.gz ${templatename}warplog.txt 1
      fi

    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate---warp each template by the resulting transforms"
    echo "   WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}0Affine.txt ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -R ${template}"
    echo "--------------------------------------------------------------------------------------"
    WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}0Affine.txt ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz ${templatename}0warp.nii.gz -R ${template}
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
    (( COUNT++ ))
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

MAXITERATIONS=30x90x20
LABELIMAGE=0 # initialize optional parameter
METRICTYPE=()
TRANSFORMATIONTYPE="GR" # initialize optional parameter
if [[ $dim == 4 ]]; then
  # we use a more constrained regularization for 4D mapping b/c we expect deformations to be relatively small and local
  TRANSFORMATIONTYPE="GR_Constrained"
fi
NUMBEROFMODALITIES=1
MODALITYWEIGHTSTRING=""
N4CORRECT=1 # initialize optional parameter
DOQSUB=1 # By default, antsMultivariateTemplateConstruction tries to do things in parallel
GRADIENTSTEP=0.25 # Gradient step size, smaller in magnitude means more smaller (more cautious) steps
ITERATIONLIMIT=4
CORES=2
TDIM=0
RIGID=0
RIGIDTYPE="" # set to an empty string to use affine initialization
range=0
REGTEMPLATES=()
TEMPLATES=()
CURRENTIMAGESET=()
XGRIDOPTS=""
SCRIPTPREPEND=""
# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS
OUTPUTNAME=antsBTP

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
 cpu_count=`sysctl -n hw.physicalcpu`
else
 cpu_count=`cat /proc/cpuinfo | grep processor | wc -l`
fi

# Provide output for Help
if [[ "$1" == "-h" ]];
    then
    Usage >&2

fi

# reading command line arguments
while getopts "A:a:b:c:d:g:h:i:j:k:m:n:o:p:s:r:t:w:x:y:z:" OPT
  do
  case $OPT in
      h) #help
   echo "$USAGE"
   exit 0
   ;;
      A) # Sharpening method
      SHARPENMETHOD=$OPTARG
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
      j) #number of cpu cores to use (default = 2)
   CORES=$OPTARG
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
   SCRIPTPREPEND=$OPTARG
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
      x) #initialization template
   XGRIDOPTS=$XGRIDOPTS
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
  MODALITYWEIGHTS=(`echo $MODALITYWEIGHTSTRING | tr 'x' "\n"`)
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
shift $shiftsize
# The invocation of $* will now read all remaining arguments into the variable IMAGESETVARIABLE
IMAGESETVARIABLE=$*
NINFILES=$(($nargs - $shiftsize))
IMAGESETARRAY=()

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

# FSL not needed anymore, all dependent on ImageMath
# #test if FSL is available in case of 4D, exit if not
# if [[  ${TDIM} -eq 4 && ${#FSLDIR} -le 0 ]];
#     then
#     setFSLPath >&2
# fi

if [[ ${NINFILES} -eq 0 ]];
  then
    echo "Please provide at least 2 filenames for the template."
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
            files=(`echo $line | tr ',' ' '`)
            if [[ ${#files[@]} -ne $NUMBEROFMODALITIES ]];
                then
                echo "The number of files in the csv file does not match the specified number of modalities."
                echo "See the -k option."
                exit 1
            fi
            for (( i = 0; i < ${#files[@]}; i++ ));
                do
                IMAGESETARRAY[$IMAGECOUNT]=${files[$i]}
                ((IMAGECOUNT++))
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
                            let "number %= $BINrange"  # Scales $number down within $range.
                        elif [[ $i -eq 15 ]];
                            then
                            let "number %= $range"  # Scales $number down within $range.
                        fi
                    done
                    #debug only
                    echo
                    echo "Random number between $FLOOR and $BINrange ---  $number"
                    #   echo "Random number between $FLOOR and $range ---  $number"

                    if [[ ${number} -lt 10 ]];
                        then
                        ImageMath $TDIM selection/vol000${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #     cp vol000${number}.nii.gz selection/
                    elif [[ ${number} -ge 10 && ${number} -lt 100 ]];
                        then
                        ImageMath $TDIM selection/vol00${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #     cp vol00${number}.nii.gz selection/
                    elif [[ ${number} -ge 100 && ${number} -lt 1000 ]];
                        then
                        ImageMath $TDIM selection/vol0${number}.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                        #     cp vol0${number}.nii.gz selection/
                    fi
                    let j++
                done
            fi
        elif [[ ${range} -gt ${nfmribins} && ${range} -lt 32 ]];
            then
            for ((i = 0; i < ${nfmribins} ; i++))
                do
                number=$RANDOM
                let "number %= $range"
                if [[ ${number} -lt 10 ]];
                    then
                    ImageMath $TDIM selection/vol0.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                    #     cp vol000${number}.nii.gz selection/
                elif [[ ${number} -ge 10 && ${number} -lt 100 ]];
                    then
                    ImageMath $TDIM selection/vol0.nii.gz ExtractSlice ${IMAGESETVARIABLE} ${number}
                    #     cp vol00${number}.nii.gz selection/
                fi
            done
        elif [[ ${range} -le ${nfmribins} ]];
            then
            ImageMath selection/$TDIM vol0.nii.gz TimeSeriesSubset ${IMAGESETVARIABLE} ${range}
            #  cp *.nii.gz selection/
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
            let k=$i+$j
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

    if [[ -n "${REGTEMPLATES[$i]}" ]];
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
          summarizeimageset $DIM ${TEMPLATES[$i]} ${STATSMETHOD} 0 ${IMAGECOMSET[@]}
          # Clean up
          rm -f ${IMAGECOMSET[@]}
    fi

    if [[ ! -s ${TEMPLATES[$i]} ]];
        then
        echo "Your template : $TEMPLATES[$i] was not created.  This indicates trouble!  You may want to check correctness of your input parameters. exiting."
        exit
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
        IMAGEMETRICSET=""
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
            do
            k=0
            let k=$i+$j
            IMAGEMETRICSET="$IMAGEMETRICSET -m MI[ ${TEMPLATES[$j]},${IMAGESETARRAY[$k]},${MODALITYWEIGHTS[$j]},32 ]"
        done

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

        exe="$ANTS $DIM $IMAGEMETRICSET -o $RIGID -i 0 $LINEARTRANSFORMPARAMS $RIGIDTYPE"

        echo "$exe" >> $qscript

        exe2='';
        pexe2='';
        pexe=" $exe > ${outdir}/job_${count}_metriclog.txt "
        for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
            do
            k=0
            let k=$i+$j
            IMGbase=`basename ${IMAGESETARRAY[$k]}`
            BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
            RIGID="${outdir}/rigid${i}_${j}_${IMGbase}"
            IMGbaseBASE=`basename ${IMAGESETARRAY[$i]}`
            BASENAMEBASE=` echo ${IMGbaseBASE} | cut -d '.' -f 1 `
            exe2="$exe2 ${WARP} $DIM ${IMAGESETARRAY[$k]} $RIGID ${outdir}/rigid${i}_0_${BASENAMEBASE}Affine.txt -R ${TEMPLATES[$j]}\n"
            pexe2="$exe2 ${WARP} $DIM ${IMAGESETARRAY[$k]} $RIGID ${outdir}/rigid${i}_0_${BASENAMEBASE}Affine.txt -R ${TEMPLATES[$j]} >> ${outdir}/job_${count}_metriclog.txt\n"
        done

        echo -e "$exe2" >> $qscript;

        if [[ $DOQSUB -eq 1 ]];
            then
            id=`qsub -cwd -S /bin/bash -N antsBuildTemplate_rigid $QSUBOPTS $qscript | awk '{print $3}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
            then
            id=`qsub -N antsrigid $QSUBOPTS -q nopreempt -l nodes=1:ppn=1 -l walltime=20:00:00 -l mem=8gb $qscript | awk '{print $1}'`
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
            id=`sbatch --job-name=antsrigid $QSUBOPTS --nodes=1 --cpus-per-task=1 --time=20:00:00 --mem=8192M $qscript | rev | cut -f1 -d\ | rev`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
            then
             # execute jobs in series
             bash $qscript
        fi
        ((count++))
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
            let k=$i-$j
            IMGbase=`basename ${IMAGESETARRAY[$i]}`
            BASENAME=` echo ${IMGbase} | cut -d '.' -f 1 `
            RIGID="${outdir}/rigid${k}_${j}_${IMGbase}"

            IMAGERIGIDSET[${#IMAGERIGIDSET[@]}]=$RIGID
        done
        echo
        echo  "AverageImages $DIM ${TEMPLATES[$j]} 2 ${IMAGERIGIDSET[@]}"

      # Don't sharpen after rigid alignment
      summarizeimageset $DIM ${TEMPLATES[$j]} ${STATSMETHOD} 0 ${IMAGERIGIDSET[@]}
      intermediateTemplateBase=`basename ${TEMPLATES[$j]}`
      cp ${TEMPLATES[$j]} ${intermediateTemplateDir}/initialRigid_${intermediateTemplateBase}

    done

    # cleanup and save output in seperate folder
    if [[ BACKUP_EACH_ITERATION -eq 1 ]];
      then

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from rigid iteration"
        echo "--------------------------------------------------------------------------------------"

        mkdir ${outdir}/rigid
        mv ${outdir}/rigid*.nii.gz ${outdir}/*Affine.txt ${outdir}/*GenericAffine.mat ${outdir}/rigid/
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
          then
            mv ${outdir}/antsBuildTemplate_rigid* ${outdir}/rigid/
            # Remove qsub scripts
            rm -f ${outdir}/job_${count}_qsub.sh
        elif [[ $DOQSUB -eq 4 ]];
          then
            mv ${outdir}/antsrigid* ${outdir}/rigid/
            # Remove qsub scripts
            rm -f ${outdir}/job_${count}_qsub.sh
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

            # Remove qsub scripts
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

ITERATLEVEL=(` echo $MAXITERATIONS | tr 'x' ' ' `)
NUMLEVELS=${#ITERATLEVEL[@]}
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
reportMappingParameters
#

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

i=0
while [[ $i -lt ${ITERATIONLIMIT} ]];
    do
    itdisplay=$((i+1))
    rm -f ${OUTPUTNAME}*Warp*.nii*
    rm -f ${outdir}/job*.sh
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
    jobIDs=""
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
            let l=$j+$k

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
            TEMPLATEbase=`basename ${TEMPLATES[$k]}`
            indir=`dirname ${IMAGESETARRAY[$j]}`
            if [[ ${#indir} -eq 0 ]];
                then
                indir=`pwd`
            fi
            IMGbase=`basename ${IMAGESETARRAY[$l]}`
            POO=${OUTPUTNAME}template-modality${k}-${IMGbase}
            OUTFN=${POO%.*.*}
            OUTFN=`basename ${OUTFN}`
            OUTFN="${OUTFN}${l}"
            DEFORMED="${outdir}/${OUTFN}${l}WarpedToTemplate.nii.gz"

            IMGbase=`basename ${IMAGESETARRAY[$j]}`
            POO=${OUTPUTNAME}${IMGbase}
            OUTWARPFN=${POO%.*.*}
            OUTWARPFN=`basename ${OUTWARPFN}`
            OUTWARPFN="${OUTWARPFN}${j}"

            if [[ $N4CORRECT -eq 1 ]];
              then
                REPAIRED="${outdir}/${OUTFN}Repaired.nii.gz"
                if [[ ! -s ${REPAIRED} ]]; then
                  exe=" $exe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${IMAGESETARRAY[$l]} -o ${REPAIRED} -r 0 -s 2\n"
                  pexe=" $pexe $N4 -d ${DIM} -b [ 200 ] -c [ 50x50x40x30,0.00000001 ] -i ${IMAGESETARRAY[$l]} -o ${REPAIRED} -r 0 -s 2  >> ${outdir}/job_${count}_metriclog.txt\n"
                fi
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${TEMPLATES[$k]},${REPAIRED},${METRICPARAMS}"
                warpexe=" $warpexe ${WARP} ${DIM} ${REPAIRED} ${DEFORMED} -R ${TEMPLATES[$k]} ${outdir}/${OUTWARPFN}Warp.nii.gz ${outdir}/${OUTWARPFN}Affine.txt\n"
                warppexe=" $warppexe ${WARP} ${DIM} ${REPAIRED} ${DEFORMED} -R ${TEMPLATES[$k]} ${outdir}/${OUTWARPFN}Warp.nii.gz ${outdir}/${OUTWARPFN}Affine.txt >> ${outdir}/job_${count}_metriclog.txt\n"
              else
                IMAGEMETRICSET="$IMAGEMETRICSET -m ${METRIC}${TEMPLATES[$k]},${IMAGESETARRAY[$l]},${METRICPARAMS}";
                warpexe=" $warpexe ${WARP} ${DIM} ${IMAGESETARRAY[$l]} ${DEFORMED} -R ${TEMPLATES[$k]} ${outdir}/${OUTWARPFN}Warp.nii.gz ${outdir}/${OUTWARPFN}Affine.txt\n"
                warppexe=" $warppexe ${WARP} ${DIM} ${IMAGESETARRAY[$l]} ${DEFORMED} -R ${TEMPLATES[$k]} ${outdir}/${OUTWARPFN}Warp.nii.gz ${outdir}/${OUTWARPFN}Affine.txt >> ${outdir}/job_${count}_metriclog.txt\n"
              fi

        done

        IMGbase=`basename ${IMAGESETARRAY[$j]}`
        POO=${OUTPUTNAME}${IMGbase}
        OUTWARPFN=${POO%.*.*}
        OUTWARPFN=`basename ${OUTWARPFN}${j}`

        LINEARTRANSFORMPARAMS="--number-of-affine-iterations 10000x10000x1000 --MI-option 32x16000"

        exe="$exe $ANTS ${DIM} $IMAGEMETRICSET -i ${MAXITERATIONS} -t ${TRANSFORMATION} -r $REGULARIZATION -o ${outdir}/${OUTWARPFN} $LINEARTRANSFORMPARAMS\n"
        exe="$exe $warpexe"

        pexe="$pexe $ANTS ${DIM} $IMAGEMETRICSET -i ${MAXITERATIONS} -t ${TRANSFORMATION} -r $REGULARIZATION -o ${outdir}/${OUTWARPFN} $LINEARTRANSFORMPARAMS >> ${outdir}/job_${count}_metriclog.txt\n"
        pexe="$pexe $warppexe"

        qscript="${outdir}/job_${count}_${i}.sh"

        echo -e $exe >> ${outdir}/job_${count}_${i}_metriclog.txt
        # 6 submit to SGE (DOQSUB=1), PBS (DOQSUB=4), PEXEC (DOQSUB=2), XGrid (DOQSUB=3) or else run locally (DOQSUB=0)
        if [[ $DOQSUB -eq 1 ]];
            then
            echo -e "$SCRIPTPREPEND" > $qscript
            echo -e "$exe" >> $qscript
            id=`qsub -cwd -N antsBuildTemplate_deformable_${i} -S /bin/bash $QSUBOPTS $qscript | awk '{print $3}'`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 4 ]];
            then
            echo -e "$SCRIPTPREPEND" > $qscript
            echo -e "$exe" >> $qscript
            id=`qsub -N antsdef${i} -q nopreempt -l nodes=1:ppn=1 -l walltime=20:00:00 -l mem=8gb $QSUBOPTS $qscript | awk '{print $1}'`
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
            id=`sbatch --job-name=antsdef${i} --nodes=1 --cpus-per-task=1 --time=20:00:00 --mem=8192M $QSUBOPTS $qscript | rev | cut -f1 -d\ | rev`
            jobIDs="$jobIDs $id"
            sleep 0.5
        elif [[ $DOQSUB -eq 0 ]];
            then
            echo -e $exe > $qscript
            bash $qscript
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

    for (( j = 0; j < $NUMBEROFMODALITIES; j++ ))
      do
        shapeupdatetotemplate ${DIM} ${TEMPLATES[$j]} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP} ${STATSMETHOD} ${SHARPENMETHOD} ${j}
        intermediateTemplateBase=`basename ${TEMPLATES[$j]}`
        cp ${TEMPLATES[$j]} ${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_${intermediateTemplateBase}
      done

    if [[ -f "${TEMPLATENAME}0warp.nii.gz" ]]
      then
        cp ${TEMPLATENAME}0warp.nii.gz ${intermediateTemplateDir}/${TRANSFORMATIONTYPE}_iteration${i}_shapeUpdateWarp.nii.gz
      fi

    if [[ BACKUP_EACH_ITERATION -eq 1 ]];
      then
        echo
        echo "--------------------------------------------------------------------------------------"
        echo " Backing up results from iteration $itdisplay"
        echo "--------------------------------------------------------------------------------------"
        mkdir ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
        cp ${TEMPLATENAME}${j}warplog.txt ${outdir}/*.cfg ${outdir}/*Affine.txt ${OUTPUTNAME}*.nii.gz ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}/
        # backup logs
        if [[ $DOQSUB -eq 1 ]];
            then
            mv ${outdir}/antsBuildTemplate_deformable_* ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
        elif [[ $DOQSUB -eq 4 ]];
            then
            mv ${outdir}/antsdef* ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
        elif [[ $DOQSUB -eq 2 ]];
            then
            mv ${outdir}/job*.txt ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
        elif [[ $DOQSUB -eq 3 ]];
            then
            rm -f ${outdir}/job_*.sh
        elif [[ $DOQSUB -eq 5 ]];
            then
            mv ${outdir}/slurm-*.out ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
            mv ${outdir}/job*.txt ${outdir}/${TRANSFORMATIONTYPE}_iteration_${i}
        fi
      else
        rm -f ${outdir}/job*.txt ${outdir}/slurm-*.out
      fi
    ((i++))
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
echo "--------------------------------------------------------------------------------------"

exit 0
