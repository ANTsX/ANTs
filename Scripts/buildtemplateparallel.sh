#!/bin/bash


VERSION="0.0.13"

# trap keyboard interrupt (control-c)
trap control_c SIGINT

# Uncomment the line below in case you have not set the ANTSPATH variable in your environment.
if [ ${#ANTSPATH} -le 3 ] ; then
  echo we guess at your ants path
  export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS
fi
if [ ! -s ${ANTSPATH}/ANTS ] ; then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

#assuming .nii.gz as default file type. This is the case for ANTS 1.7 and up

function Usage {
    cat <<USAGE

Please reference http://www.ncbi.nlm.nih.gov/pubmed/20851191 when employing this script in your studies.
A reproducible evaluation of ANTs similarity metric performance in brain image registration.
Avants BB, Tustison NJ, Song G, Cook PA, Klein A, Gee JC. Neuroimage, 2011.

Also see http://www.ncbi.nlm.nih.gov/pubmed/19818860 for more details.  The script has been updated and improved since this publication.

Usage:

`basename $0` -d ImageDimension -o OUTPREFIX <other options> <images>

Example Case:

 echo " bash $0 -d 3 -m 30x50x20 -t GR  -s CC -c 1 -o MY -z InitialTemplate.nii.gz  *RF*T1x.nii.gz "
 echo " in this case you use 30x50x20 iterations per registration "
 echo " 4 iterations over template creation (that is the default) "
 echo " with Greedy-SyN and CC metrics to guide the mapping. "
 echo " Output is prepended with MY and the initial template is InitialTemplate.nii.gz (optional). "
 echo " the -c option is set to 1 which will try to use SGE to distribute the computation. "
 echo " if you do not have SGE, use -c 0 or -c 2 --- read the help."

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options):

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
	  ImageDimension: 4 (for 3 dimensional registration of time-series; requires FSL)

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.

NB: All images to be added to the template should be in the same directory, and this script
should be invoked from that directory.

Optional arguments:

     -c:  Control for parallel computation (default 1) -- 0 == run serially,  1 == SGE qsub,
          2 == use PEXEC (localhost),  3 == Apple XGrid

     -g:  Gradient step size (default 0.25) -- smaller in magnitude results in more cautious steps

     -i:  Iteration limit (default 4) -- iterations of the template construction (Iteration limit)*NumImages registrations.

     -j:  Number of cpu cores to use (default: 2; -- requires "-c 2")

     -m:  Max-iterations in each registration

     -n:  N4BiasFieldCorrection of moving image (default 1) -- 0 == off, 1 == on

     -p:  Commands to prepend to job scripts (e.g., change into appropriate directory, set paths, etc)

     -r:  Do rigid-body registration of inputs before creating template (default 0) -- 0 == off 1 == on. Only useful when
          you do not have an initial template

     -s:  Type of similarity metric used for registration.

     -t:  Type of transformation model used for registration.

     -x:  XGrid arguments (e.g., -x "-p password -h controlhost")

     -z:  Use this this volume as the target of all inputs. When not used, the script
          will create an unbiased starting point by averaging all inputs. Use the full path!

--------------------------------------------------------------------------------------
ANTS was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

--------------------------------------------------------------------------------------
script adapted by N.M. van Strien, http://www.mri-tutorial.com | NTNU MR-Center
--------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------
Apple XGrid support by Craig Stark
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

function Help {
    cat <<HELP
Please reference http://www.ncbi.nlm.nih.gov/pubmed/20851191 when employing this script in your studies.
A reproducible evaluation of ANTs similarity metric performance in brain image registration.
Avants BB, Tustison NJ, Song G, Cook PA, Klein A, Gee JC. Neuroimage, 2011.

Also see http://www.ncbi.nlm.nih.gov/pubmed/19818860 for more details.  The script has been updated and improved since this publication.

./buildtemplateparallel.sh  will make a template out of the input files using an elastic
or diffeomorphic transformation. This script builds a template iteratively from the input
images and uses Sun Grid Engine (SGE) or multiple cpu cores on the localhost (min 2) to
parallelize the registration of each subject to the template.

Usage:

`basename $0` -d ImageDimension -o OUTPREFIX <other options> <images>

Compulsory arguments (minimal command line requires SGE cluster, otherwise use -c & -j options)::

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
	  ImageDimension: 4 (for 3 dimensional registration of time-series; requires FSL)

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.

NB: All files to be added to the template should be in the same directory.

Optional arguments:

     -c:  Control for parallel computation   --- if set to zero, run serially, if set to 2 , use PEXEC (localhost), if set to 1 , use SGE qsub.

     -g:  Gradient step size; smaller in magnitude results in more cautious steps (default 0.25)

     -i:  Iteration limit (default = 4) for template construction. requires 4*NumImages registrations.

     -j:  Number of cpu cores to use (default: 2; --- set -c option to 2 to use this .

	  The optimal number of cpu cores to use for template generation depends on the availability of cores, the amount of
	  free working memory (RAM) and the resolution of the data. High resolution datasets typically require more RAM during
	  processing. Running out of RAM during a calculation will slow down all processing on your computer.

     -m:  Max-iterations

          Max-Iterations in form: JxKxL where
	     J = max iterations at coarsest resolution (here, reduce by power of 2^2)
	     K = middle resolution iterations (here,reduce by power of 2)
	     L = fine resolution iterations (here, full resolution) !!this level takes much
                 more time per iteration!!

	  Adding an extra value before JxKxL (i.e. resulting in IxJxKxL) would add another
	  iteration level.

     -n:  N4BiasFieldCorrection of moving image ( 0 = off; 1 = on (default) )

     -p:  Commands to prepend to job scripts (e.g., change into appropriate directory, set paths, etc)

     -r:  Do rigid-body registration of inputs before creating template (default 0) -- 0 == off 1 == on. Only useful when
          you do not have an initial template

          In case a template is specified (-z option), all inputs are registered to that template. If
          no template is specified, the inputs will be registered to the averaged input.

     -s:  Type of similarity metric used for registration.

	     For intramodal image registration, use:
	     CC = cross-correlation
	     MI = mutual information
	     PR = probability mapping (default)
      MSQ = mean square difference (Demons-like)
      SSD = sum of squared differences

	     For intermodal image registration, use:
	     MI = mutual information
	     PR = probability mapping (default)

     -t:  Type of transformation model used for registration.

	     For elastic image registration, use:
	     EL = elastic transformation model (less deformation possible)

	     For diffeomorphic image registration, use:
	     SY = SyN with time (default) with arbitrary number of time points in time discretization
	     S2 = SyN with time optimized specifically for 2 time points in the time discretization
	     GR = Greedy SyN
	     EX = Exponential
             DD = Diffeomorphic Demons style exponential mapping

     -x:  XGrid arguments (e.g., -x "-p password -h controlhost")

     -z:  Use this this volume as the target of all inputs. When not used, the script
          will create an unbiased starting point by averaging all inputs. Use the full path!

Requirements:

This scripts relies on the following scripts in your $ANTSPATH directory. The script
will terminate prematurely if these files are not present or are not executable.
- antsIntroduction.sh
- pexec.sh
- waitForSGEQJobs.pl (only for use with Sun Grid Engine)
- ANTSpexec.sh (only for use with localhost parallel execution)
- waitForXGridJobs.pl (only for use with Apple XGrid)

For 4D template building FSL is also needed

--------------------------------------------------------------------------------------
Get the latest ANTS version at:
--------------------------------------------------------------------------------------
http://sourceforge.net/projects/advants/

--------------------------------------------------------------------------------------
Read the ANTS documentation at:
--------------------------------------------------------------------------------------
http://picsl.upenn.edu/ANTS/

--------------------------------------------------------------------------------------
ANTS was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

--------------------------------------------------------------------------------------
script adapted by N.M. van Strien, http://www.mri-tutorial.com | NTNU MR-Center
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
Apple XGrid support by Craig Stark
--------------------------------------------------------------------------------------

HELP
    exit 1
}

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

Alternatively, edit this script ( $0 ) to set up this parameter correctly.

SETPATH
    exit 1
}

setFSLPath() {
    cat <<SETFSLPATH
--------------------------------------------------------------------------------------
Error locating FSL
--------------------------------------------------------------------------------------
The FSLDIR environment variable is not set. Please add the FSLDIR variable.

see the FSL website for more information about installation:

http://www.fmrib.ox.ac.uk/fsl/fsl/downloading.html

SETFSLPATH
#    exit 1
}


function reportMappingParameters {
    cat <<REPORTMAPPINGPARAMETERS

--------------------------------------------------------------------------------------
 Mapping parameters
--------------------------------------------------------------------------------------
 ANTSPATH is $ANTSPATH

 Dimensionality:			$DIM
 N4BiasFieldCorrection:			$N4CORRECT
 Similarity Metric:			$METRICTYPE
 Transformation:			$TRANSFORMATIONTYPE
 Regularization:			$REGULARIZATION
 MaxIterations:				$MAXITERATIONS
 Number Of MultiResolution Levels:	$NUMLEVELS
 OutputName prefix:			$OUTPUTNAME
 Template:  				$TEMPLATE
 Template Update Steps:			$ITERATIONLIMIT
 Template population:	   		$IMAGESETVARIABLE
--------------------------------------------------------------------------------------
REPORTMAPPINGPARAMETERS
}

function shapeupdatetotemplate {

    # local declaration of values
    dim=${DIM}
    template=${TEMPLATE}
    templatename=${TEMPLATENAME}
    outputname=${OUTPUTNAME}
    gradientstep=-${GRADIENTSTEP}

# debug only
# echo $dim
# echo ${template}
# echo ${templatename}
# echo ${outputname}
# echo ${outputname}*formed.nii*
# echo ${gradientstep}

# We find the average warp to the template and apply its inverse to the template image
# This keeps the template shape stable over multiple iterations of template building

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 1"
    echo "--------------------------------------------------------------------------------------"
    ${ANTSPATH}AverageImages $dim ${template} 1 ${outputname}*formed.nii.gz

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 2"
    echo "--------------------------------------------------------------------------------------"

	${ANTSPATH}AverageImages $dim ${templatename}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 3"
    echo "--------------------------------------------------------------------------------------"
	${ANTSPATH}MultiplyImages $dim ${templatename}warp.nii.gz ${gradientstep} ${templatename}warp.nii.gz


    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 4"
    echo "--------------------------------------------------------------------------------------"
    rm -f ${templatename}Affine.txt

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 5"
    echo "--------------------------------------------------------------------------------------"

    # Averaging and inversion code --- both are 1st order estimates.
#    if [ ${dim} -eq 2   ] ; then
#      ANTSAverage2DAffine ${templatename}Affine.txt ${outputname}*Affine.txt
#    elif [ ${dim} -eq 3  ] ; then
#      ANTSAverage3DAffine ${templatename}Affine.txt ${outputname}*Affine.txt
#    fi
    ${ANTSPATH}AverageAffineTransform ${dim} ${templatename}Affine.txt ${outputname}*Affine.txt

    ${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warp.nii.gz ${templatename}warp.nii.gz -i  ${templatename}Affine.txt -R ${template}
    ${ANTSPATH}WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}Affine.txt ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz -R ${template}

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 6"
    ${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warp.nii.gz ${templatename}warplog.txt 1

}

function ANTSAverage2DAffine {

    OUTNM=${templatename}Affine.txt
    FLIST=${outputname}*Affine.txt
    NFILES=0
    PARAM1=0
    PARAM2=0
    PARAM3=0
    PARAM4=0
    PARAM5=0
    PARAM6=0
    PARAM7=0
    PARAM8=0
    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM1=` awk -v a=$PARAM1 -v b=$x 'BEGIN{print (a + b)}' ` ;  let NFILES=$NFILES+1  ; done
    PARAM1=` awk -v a=$PARAM1 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM2=` awk -v a=$PARAM2 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM2=` awk -v a=$PARAM2 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM3=` awk -v a=$PARAM3 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM3=` awk -v a=$PARAM3 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 5  `
    for x in $LL ; do PARAM4=` awk -v a=$PARAM4 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM4=` awk -v a=$PARAM4 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 6  `
    for x in $LL ; do PARAM5=` awk -v a=$PARAM5 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM5=0 # ` awk -v a=$PARAM5 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 7  `
    for x in $LL ; do PARAM6=` awk -v a=$PARAM6 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM6=0 # ` awk -v a=$PARAM6 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 2  `
    for x in $LL ; do PARAM7=` awk -v a=$PARAM7 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM7=` awk -v a=$PARAM7 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM8=` awk -v a=$PARAM8 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM8=` awk -v a=$PARAM8 -v b=$NFILES 'BEGIN{print (a / b)}' `

    echo "# Insight Transform File V1.0 " > $OUTNM
    echo "# Transform 0 " >> $OUTNM
    echo "Transform: MatrixOffsetTransformBase_double_2_2  " >> $OUTNM
    echo "Parameters:  $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6  " >> $OUTNM
    echo "FixedParameters: $PARAM7 $PARAM8 " >> $OUTNM


}

function ANTSAverage3DAffine {

    OUTNM=${templatename}Affine.txt
    FLIST=${outputname}*Affine.txt
    NFILES=0
    PARAM1=0
    PARAM2=0
    PARAM3=0
    PARAM4=0
    PARAM5=0
    PARAM6=0
    PARAM7=0
    PARAM8=0
    PARAM9=0
    PARAM10=0
    PARAM11=0
    PARAM12=0
    PARAM13=0
    PARAM14=0
    PARAM15=0
    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM1=` awk -v a=$PARAM1 -v b=$x 'BEGIN{print (a + b)}' ` ;  let NFILES=$NFILES+1  ; done
    PARAM1=` awk -v a=$PARAM1 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM2=` awk -v a=$PARAM2 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM2=` awk -v a=$PARAM2 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM3=` awk -v a=$PARAM3 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM3=` awk -v a=$PARAM3 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 5  `
    for x in $LL ; do PARAM4=` awk -v a=$PARAM4 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM4=` awk -v a=$PARAM4 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 6  `
    for x in $LL ; do PARAM5=` awk -v a=$PARAM5 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM5=` awk -v a=$PARAM5 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 7  `
    for x in $LL ; do PARAM6=` awk -v a=$PARAM6 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM6=` awk -v a=$PARAM6 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 8  `
    for x in $LL ; do PARAM7=` awk -v a=$PARAM7 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM7=` awk -v a=$PARAM7 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 9  `
    for x in $LL ; do PARAM8=` awk -v a=$PARAM8 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM8=` awk -v a=$PARAM8 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 10  `
    for x in $LL ; do PARAM9=` awk -v a=$PARAM9 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM9=` awk -v a=$PARAM9 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 11  `
    for x in $LL ; do PARAM10=` awk -v a=$PARAM10 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM10=0 # ` awk -v a=$PARAM10 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 12  `
    for x in $LL ; do PARAM11=` awk -v a=$PARAM11 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM11=0 # ` awk -v a=$PARAM11 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 13  `
    for x in $LL ; do PARAM12=` awk -v a=$PARAM12 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM12=0 # ` awk -v a=$PARAM12 -v b=$NFILES 'BEGIN{print (a / b)}' `

# origin params below

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM13=` awk -v a=$PARAM13 -v b=$x 'BEGIN{print (a + b)}' ` ;  done
    PARAM13=` awk -v a=$PARAM13 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM14=` awk -v a=$PARAM14 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM14=` awk -v a=$PARAM14 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM15=` awk -v a=$PARAM15 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM15=` awk -v a=$PARAM15 -v b=$NFILES 'BEGIN{print (a / b)}' `

    echo "# Insight Transform File V1.0 " > $OUTNM
    echo "# Transform 0 " >> $OUTNM
    echo "Transform: MatrixOffsetTransformBase_double_3_3  " >> $OUTNM
    echo "Parameters:  $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6  $PARAM7 $PARAM8 $PARAM9 $PARAM10 $PARAM11 $PARAM12  " >> $OUTNM
    echo "FixedParameters: $PARAM13 $PARAM14 $PARAM15 " >> $OUTNM

}

function jobfnamepadding {

    files=`ls job*.sh`
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

cleanup()
# example cleanup function
{

  cd ${currentdir}/

  echo -en "\n*** Performing cleanup, please wait ***\n"

# 1st attempt to kill all remaining processes
# put all related processes in array
  runningANTSpids=( `ps -C ANTS -C N4BiasFieldCorrection -C fslsplit | awk '{ printf "%s\n", $1 ; }'` )

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

#initializing variables with global scope
time_start=`date +%s`
currentdir=`pwd`
nargs=$#

MAXITERATIONS=30x90x20
LABELIMAGE=0 # initialize optional parameter
METRICTYPE=CC # initialize optional parameter
TRANSFORMATIONTYPE="GR" # initialize optional parameter
N4CORRECT=1 # initialize optional parameter
DOQSUB=1 # By default, buildtemplateparallel tries to do things in parallel
GRADIENTSTEP=0.25 # Gradient step size, smaller in magnitude means more smaller (more cautious) steps
ITERATIONLIMIT=4
CORES=2
TDIM=0
RIGID=0
RIGIDTYPE=" --do-rigid" # set to an empty string to use affine initialization
range=0
REGTEMPLATE=target
XGRIDOPTS=""
SCRIPTPREPEND=""
# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS



if [ $OSTYPE == 'darwin' ]
	then
	cpu_count=`sysctl -n hw.physicalcpu`
else
	cpu_count=`cat /proc/cpuinfo | grep processor | wc -l`
fi

##Getting system info from linux can be done with these variables.
# RAM=`cat /proc/meminfo | sed -n -e '/MemTotal/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# RAMfree=`cat /proc/meminfo | sed -n -e '/MemFree/p' | awk '{ printf "%s %s\n", $2, $3 ; }' | cut -d " " -f 1`
# cpu_free_ram=$((${RAMfree}/${cpu_count}))

# Provide output for Help
if [ "$1" == "-h" ]
    then
    Help >&2

fi
OUTPUTNAME=antsBTP
# reading command line arguments
while getopts "c:d:g:i:j:h:m:n:o:p:s:r:t:z:" OPT
  do
  case $OPT in
      h) #help
	  echo "$USAGE"
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
	  DIM=$OPTARG
	  if [[ ${DIM} -eq 4 ]] ; then
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
      m) #max iterations other than default
	  MAXITERATIONS=$OPTARG
	  ;;
      n) #apply bias field correction
	  N4CORRECT=$OPTARG
	  ;;
      o) #output name prefix
	  OUTPUTNAME=$OPTARG
	  TEMPLATENAME=${OUTPUTNAME}template
	  TEMPLATE=${TEMPLATENAME}.nii.gz
	  ;;
      p) #Script prepend
	  SCRIPTPREPEND=$OPTARG
	  ;;
      s) #similarity model
	  METRICTYPE=$OPTARG
	  ;;
      r) #start with rigid-body registration
	  RIGID=$OPTARG
	  ;;
      t) #transformation model
	  TRANSFORMATIONTYPE=$OPTARG
	  ;;
      z) #initialization template
	  REGTEMPLATE=$OPTARG
	  ;;
      \?) # getopts issues an error message
      echo "$USAGE" >&2
      exit 1
      ;;
  esac
done

# Provide different output for Usage and Help
if [ ${TDIM} -eq 4 ] && [ $nargs -lt 5 ]
    then
    Usage >&2
elif [ ${TDIM} -eq 4 ] && [ $nargs -eq 5 ]
    then
    echo ""
    # This option is required to run 4D template creation on SGE with a minimal command line
elif [ $nargs -lt 6 ]
    then
    Usage >&2
fi

if [ $DOQSUB -eq 1 ] ; then
  qq=`which  qsub`
  if [  ${#qq} -lt 1 ] ; then
    echo do you have qsub?  if not, then choose another c option ... if so, then check where the qsub alias points ...
    exit
  fi
fi

#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
    then
    setPath >&2
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
#test if FSL is available in case of 4D, exit if not
if [  ${TDIM} -eq 4 ] && [  ${#FSLDIR} -le 0 ]
    then
    setFSLPath >&2
fi

if [ ${NINFILES} -eq 0 ]
    then
    echo "Please provide at least 2 filenames for the template."
    echo "Use `basename $0` -h for help"
    exit 1
elif [[ ${NINFILES} -eq 1 ]] # && [[ -s ${FSLDIR}/bin/fslnvols ]]
    then

#    range=`fslnvols ${IMAGESETVARIABLE}`
    range=` ${ANTSPATH}ImageMath $TDIM a nvols ${IMAGESETVARIABLE}  `

    if [ ${range} -eq 1 ] && [ ${TDIM} -ne 4 ]
	then
	echo "Please provide at least 2 filenames for the template."
	echo "Use `basename $0` -h for help"
	exit 1
    elif [ ${range} -gt 1 ] && [ ${TDIM} -ne 4 ]
	then
	echo "This is a multivolume file. Use -d 4"
	echo "Use `basename $0` -h for help"
	exit 1
    elif [ ${range} -gt 1 ] && [ ${TDIM} -eq 4 ]
	then
	echo "--------------------------------------------------------------------------------------"
	echo " Creating template of 4D input. "
	echo "--------------------------------------------------------------------------------------"
		#splitting volume
		#setting up working dirs
	mkdir tmp
	mkdir tmp/selection

		#split the 4D file into 3D elements
	cp ${IMAGESETVARIABLE} tmp/
	cd tmp/
        ${ANTSPATH}ImageMath $TDIM selection/vol.nii.gz TimeSeriesSubset ${IMAGESETVARIABLE} 16

#	fslsplit ${IMAGESETVARIABLE}
#	rm -f ${IMAGESETVARIABLE}

		# selecting 16 volumes randomly from the timeseries for averaging, placing them in tmp/selection. folder
#	for ((i = 0; i < 16 ; i++))
#	  do
#	  number=$RANDOM
#	  let "number %= $range"
#
#	  if [ ${number} -lt 10 ]
#	      then
#	      cp vol000${number}.nii.gz selection/
#	  elif [ ${number} -ge 10 ] && [ ${number} -lt 100 ]
#	      then
#	      cp vol00${number}.nii.gz selection/
#	  elif [ ${number} -ge 100 ] && [ ${number} -lt 1000 ]
#	      then
#	      cp vol0${number}.nii.gz selection/
#	  fi
#	done

		# set filelist variable
	cd selection/
	IMAGESETVARIABLE=`ls *.nii.gz`

    fi
fi

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTSSCRIPTNAME=${ANTSPATH}antsIntroduction.sh
PEXEC=${ANTSPATH}ANTSpexec.sh
SGE=${ANTSPATH}waitForSGEQJobs.pl
XGRID=${ANTSPATH}waitForXGridJobs.pl

for FLE in $ANTSSCRIPTNAME $PEXEC $SGE $XGRID
  do
  if [ ! -x $FLE  ] ;
      then
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " FILE $FLE DOES NOT EXIST -- OR -- IS NOT EXECUTABLE !!! $0 will terminate."
      echo "--------------------------------------------------------------------------------------"
      echo " if the file is not executable, please change its permissions. "
      exit 1
  fi
done

# exit
# check for an initial template image and perform rigid body registration if requested
if [ ! -s $REGTEMPLATE ]
    then
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " No initial template exists. Creating a population average image from the inputs."
    echo "--------------------------------------------------------------------------------------"
    ${ANTSPATH}AverageImages $DIM $TEMPLATE 1 $IMAGESETVARIABLE
else
    echo
    echo "--------------------------------------------------------------------------------------"
    echo " Initial template found.  This will be used for guiding the registration. use : $REGTEMPLATE and $TEMPLATE "
    echo "--------------------------------------------------------------------------------------"
	# now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.

    cp ${REGTEMPLATE} ${TEMPLATE}

fi


if [ ! -s $TEMPLATE ] ; then
  echo Your template : $TEMPLATE was not created.  This indicates trouble!  You may want to check correctness of your input parameters. exiting.
  exit
fi
# remove old job bash scripts
rm -f job*.sh

if [ "$RIGID" -eq 1 ] ;
    then
    count=0
    jobIDs=""

    RIGID_IMAGESET=""

    for IMG in $IMAGESETVARIABLE
      do

      RIGID_IMAGESET="$RIGID_IMAGESET rigid_${IMG}"

      BASENAME=` echo ${IMG} | cut -d '.' -f 1 `

      exe=" ${ANTSPATH}ANTS $DIM -m MI[${TEMPLATE},${IMG},1,32] -o rigid_${IMG} -i 0 --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000 $RIGIDTYPE"
      exe2="${ANTSPATH}WarpImageMultiTransform $DIM ${IMG} rigid_${IMG} rigid_${BASENAME}Affine.txt -R ${TEMPLATE}"
      pexe=" $exe >> job_${count}_metriclog.txt "

      qscript="job_${count}_qsub.sh"

	  echo "$SCRIPTPREPEND" > $qscript

      echo "$exe" >> $qscript

      echo "$exe2" >> $qscript

      if [ $DOQSUB -eq 1 ] ; then
	  id=`qsub -cwd -S /bin/bash -N antsBuildTemplate_rigid -v ANTSPATH=$ANTSPATH $QSUBOPTS $qscript | awk '{print $3}'`
	  jobIDs="$jobIDs $id"
	  sleep 0.5
      elif  [ $DOQSUB -eq 2 ] ; then
		# Send pexe and exe2 to same job file so that they execute in series
		echo $pexe >> job${count}_r.sh
		echo $exe2 >> job${count}_r.sh
      elif  [ $DOQSUB -eq 3 ] ; then
	id=`xgrid $XGRIDOPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
	#echo "xgrid $XGRIDOPTS -job submit /bin/bash $qscript"
		jobIDs="$jobIDs $id"
	  elif  [ $DOQSUB -eq 0 ] ; then
	  # execute jobs in series
	  $exe
	  $exe2
      fi

      ((count++))
    done


    if [ $DOQSUB -eq 1 ];
	then
	# Run jobs on SGE and wait to finish
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Starting ANTS rigid registration on cluster. Submitted $count jobs "
	echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
	${ANTSPATH}waitForSGEQJobs.pl 1 60 $jobIDs

	# Returns 1 if there are errors
	if [ ! $? -eq 0 ]; then
	    echo "qsub submission failed - jobs went into error state"
	    exit 1;
	fi
    fi

    # Run jobs on localhost and wait to finish
    if [ $DOQSUB -eq 2 ];
	then
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Starting ANTS rigid registration on max ${CORES} cpucores. "
	echo " Progress can be viewed in job*_metriclog.txt"
	echo "--------------------------------------------------------------------------------------"
	jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
	chmod +x job*.sh
	$PEXEC -j ${CORES} "sh" job*.sh
    fi

    if [ $DOQSUB -eq 3 ];
	then
	# Run jobs on XGrid and wait to finish
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Starting ANTS rigid registration on XGrid cluster. Submitted $count jobs "
	echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. Rigid registration is quick, so poll queue every 60 seconds
	${ANTSPATH}waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 30 $jobIDs
	# Returns 1 if there are errors
	if [ ! $? -eq 0 ]; then
	    echo "XGrid submission failed - jobs went into error state"
	    exit 1;
	fi
    fi

    # Update template
    ${ANTSPATH}AverageImages $DIM $TEMPLATE 1 $RIGID_IMAGESET

    # cleanup and save output in seperate folder

    mkdir rigid
    mv *.cfg rigid_*.nii.gz rigid_*.nii *Affine.txt rigid/

    # backup logs
    if [ $DOQSUB -eq 1 ];
	then
	mv antsBuildTemplate_rigid* rigid/

        # Remove qsub scripts
	rm -f job_${count}_qsub.sh

    elif [ $DOQSUB -eq 2 ];
		then
		mv job*.txt rigid/
	elif [ $DOQSUB -eq 3 ];
		then
		rm -f job_*_qsub.sh
    fi


fi # endif RIGID



# Begin Main Loop
ITERATLEVEL=(` echo $MAXITERATIONS | tr 'x' ' ' `)
NUMLEVELS=${#ITERATLEVEL[@]}

# debugging only
#echo $ITERATLEVEL
#echo $NUMLEVELS
#echo ${ITERATIONLIMIT}

echo "--------------------------------------------------------------------------------------"
echo " Start to build template: ${TEMPLATE}"
echo "--------------------------------------------------------------------------------------"
reportMappingParameters

i=0
while [  $i -lt ${ITERATIONLIMIT} ]
  do

  itdisplay=$((i+1))

  rm -f  ${OUTPUTNAME}*Warp*.nii*
  rm -f job*.sh

# Used to save time by only running coarse registration for the first couple of iterations
# But with decent initialization, this is probably not worthwhile.
# If you uncomment this, replace MAXITERATIONS with ITERATIONS in the call to ants below
#
# # For the first couple of iterations, use high-level registration only
# # eg if MAXITERATIONS=30x90x20, then for iteration 0, do 30x0x0
# # for iteration 1 do 30x90x0, then do 30x90x20 on subsequent iterations
# if [ $i -gt $((NUMLEVELS - 1)) ]
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
  for IMG in $IMAGESETVARIABLE
    do
    # 1 determine working dir
    dir=`pwd`

    # 2 determine new filename
    POO=${OUTPUTNAME}${IMG}


    # 3 Make variable OUTFILENAME and remove anything behind . ; for example .nii.gz.gz
    OUTFN=${POO%.*.*}

    # 4 Test if outputfilename has only a single extention and remove that
    if [ ${#OUTFN} -eq ${#POO} ]
	then
	OUTFN=${OUTPUTNAME}${IMG%.*}
    fi

    # 5 prepare registration command
    exe="${ANTSSCRIPTNAME} -d ${DIM} -r ${dir}/${TEMPLATE} -i ${dir}/${IMG} -o ${dir}/${OUTFN} -m ${MAXITERATIONS} -n ${N4CORRECT} -s ${METRICTYPE} -t ${TRANSFORMATIONTYPE} "
    pexe=" $exe >> job_${count}_${i}_metriclog.txt "

    # 6 submit to SGE or else run locally
    if [ $DOQSUB -eq 1 ]; then
	id=`qsub -cwd -N antsBuildTemplate_deformable_${i} -S /bin/bash -v ANTSPATH=$ANTSPATH $QSUBOPTS $exe | awk '{print $3}'`
	jobIDs="$jobIDs $id"
	sleep 0.5
    elif [ $DOQSUB -eq 2 ] ; then
	echo $pexe
	echo $pexe >> job${count}_${i}.sh
    elif [ $DOQSUB -eq 3 ] ; then
      qscript="job_${count}_${i}.sh"
      exe="${ANTSSCRIPTNAME} -d ${DIM} -r ./${TEMPLATE} -i ./${IMG} -o ./${OUTFN} -m ${MAXITERATIONS} -n ${N3CORRECT} -s ${METRICTYPE} -t ${TRANSFORMATIONTYPE} "
	  echo "$SCRIPTPREPEND" > $qscript
	  echo "$exe" >> $qscript
      id=`xgrid $XGRIDOPTS -job submit /bin/bash $qscript | awk '{sub(/;/,"");print $3}' | tr '\n' ' ' | sed 's:  *: :g'`
	  jobIDs="$jobIDs $id"
    elif  [ $DOQSUB -eq 0 ] ; then
	bash $exe
    fi

    # counter updated, but not directly used in this loop
    count=`expr $count + 1`;

#		echo " submitting job number $count " # for debugging only
  done

  # SGE wait for script to finish
  if [ $DOQSUB -eq 1 ];
      then
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " Starting ANTS registration on SGE cluster. Iteration: $itdisplay of $ITERATIONLIMIT"
      echo "--------------------------------------------------------------------------------------"

      # now wait for the stuff to finish - this will take a while so poll queue every 10 mins
      ${ANTSPATH}waitForSGEQJobs.pl 1 600 $jobIDs

      if [ ! $? -eq 0 ]; then
	  echo "qsub submission failed - jobs went into error state"
	  exit 1;
      fi

  fi


  # Run jobs on localhost and wait to finish
  if [ $DOQSUB -eq 2 ];
      then
      echo
      echo "--------------------------------------------------------------------------------------"
      echo " Starting ANTS registration on max ${CORES} cpucores. Iteration: $itdisplay of $ITERATIONLIMIT"
      echo " Progress can be viewed in job*_${i}_metriclog.txt"
      echo "--------------------------------------------------------------------------------------"
      jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
      chmod +x job*.sh
      $PEXEC -j ${CORES} sh job*.sh
  fi

  if [ $DOQSUB -eq 3 ];
	then
	# Run jobs on XGrid and wait to finish
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Starting ANTS rigid registration on XGrid cluster. Submitted $count jobs "
	echo "--------------------------------------------------------------------------------------"
        # now wait for the jobs to finish. This is slow, so poll less often
	${ANTSPATH}waitForXGridJobs.pl -xgridflags "$XGRIDOPTS" -verbose -delay 300 $jobIDs
	# Returns 1 if there are errors
	if [ ! $? -eq 0 ]; then
	    echo "XGrid submission failed - jobs went into error state"
	    exit 1;
	fi
  fi

  shapeupdatetotemplate ${DIM} ${TEMPLATE} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP}

  echo
  echo "--------------------------------------------------------------------------------------"
  echo " Backing up results from iteration $itdisplay"
  echo "--------------------------------------------------------------------------------------"

  mkdir ${TRANSFORMATIONTYPE}_iteration_${i}
  cp ${TEMPLATENAME}warp*log.txt *.cfg *${OUTPUTNAME}*.nii.gz ${TRANSFORMATIONTYPE}_iteration_${i}/

  # backup logs
  if [ $DOQSUB -eq 1 ];
      then
      mv antsBuildTemplate_deformable_* ${TRANSFORMATIONTYPE}_iteration_${i}

  elif [ $DOQSUB -eq 2 ];
      then
      mv job*.txt ${TRANSFORMATIONTYPE}_iteration_${i}
  elif [ $DOQSUB -eq 3 ];
      then
      rm -f job_*.sh
  fi



  ((i++))

done

# end main loop

rm -f job*.sh

#cleanup of 4D files
if [ "${range}" -gt 1 ] && [ "${TDIM}" -eq 4 ]
    then
    mv ${currentdir}/tmp/selection/${TEMPLATE} ${currentdir}/
    cd ${currentdir}
    rm -rf ${currentdir}/tmp/
fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done creating: ${TEMPLATE}"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
