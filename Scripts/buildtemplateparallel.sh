# !/bin/sh

VERSION="0.0.10"

# Uncomment the line below in case you have not set the ANTSPATH variable in your environment.
# ANTSPATH=

#assuming .nii.gz as default file type. This is the case for ANTS 1.7 and up

function Usage {
    cat <<USAGE

Usage:

$0 -d ImageDimension -o OUTPREFIX <other options> <images>

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
	  ImageDimension: 4 (for 3 dimensional registration of time-series; requires FSL)

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.

NB: All files to be added to the template should be in the same directory.

Optional arguments:

     -c:  Use SGE cluster (cannot be used in combination with -j; requires SGE)

     -g:  Gradient step size; smaller in magnitude results in more cautious steps (default 0.25)

     -i:  Iteration limit (default = 4) for template construction. requires 4*NumImages registrations.

     -j:  Number of cpu cores to use (default 2)

     -m:  Max-iterations

     -n:  N3BiasFieldCorrection of moving image ( 0 = off; 1 = on (default) )

     -r:  Do rigid-body registration of inputs before creating template ( 0 = off (default); 1 = on )

     -s:  Type of similarity metric used for registration.

     -t:  Type of transformation model used for registration.

     -z:  Use this this volume as the target of all inputs. When not used, the script
          will create an unbiased starting point by averaging all inputs.

--------------------------------------------------------------------------------------
ANTS was created by:
--------------------------------------------------------------------------------------
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

--------------------------------------------------------------------------------------
script adapted by N.M. van Strien, http://www.mri-tutorial.com | NTNU MR-Center
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

function Help {
    cat <<Help

./buildtemplateparallel.sh  will make a template out of the input files using an elastic
or diffeomorphic transformation. This script builds a template iteratively from the input
images and uses Sun Grid Engine (SGE) or multiple cpu cores on the localhost (min 2) to
parallelize the registration of each subject to the template.

Usage:

sh buildtemplateparallel.sh -d ImageDimension -o OUTPREFIX <other options> <images>

Compulsory arguments:

     -d:  ImageDimension: 2 or 3 (for 2 or 3 dimensional registration of single volume)
	  ImageDimension: 4 (for 3 dimensional registration of time-series; requires FSL)

     -o:  OUTPREFIX; A prefix that is prepended to all output files.

<images>  List of images in the current directory, eg *_t1.nii.gz. Should be at the end
          of the command.

NB: All files to be added to the template should be in the same directory.

Optional arguments:

     -c:  Use SGE cluster (cannot be used in combination with -j; requires SGE)

     -g:  Gradient step size; smaller in magnitude results in more cautious steps (default 0.25)

     -i:  Iteration limit (default = 4) for template construction. requires 4*NumImages registrations.

     -j:  Number of cpu cores to use (default 2)

     -m:  Max-iterations

          Max-Iterations in form: JxKxL where
	     J = max iterations at coarsest resolution (here, reduce by power of 2^2)
	     K = middle resolution iterations (here,reduce by power of 2)
	     L = fine resolution iterations (here, full resolution) !!this level takes much
                 more time per iteration!!

	  Adding an extra value before JxKxL (i.e. resulting in IxJxKxL) would add another
	  iteration level.

     -n:  N3BiasFieldCorrection of moving image ( 0 = off; 1 = on (default) )

     -r:  Do rigid-body registration of inputs before creating template ( 0 = off (default); 1 = on )

          In case a template is specified (-z option), all inputs are registered to that template. If
          no template is specified, the inputs will be registered to the averaged input.

     -s:  Type of similarity metric used for registration.

	     For intramodal image registration, use:
	     CC = cross-correlation
	     MI = mutual information
	     PR = probability mapping (default)
	     MSQ = mean square difference

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

     -z:  Use this this volume as the target of all inputs. When not used, the script
          will create an unbiased starting point by averaging all inputs.

Requirements:

This scripts relies on the following scripts in your $ANTSPATH directory. The script
will terminate prematurely if these files are not present.
- antsIntroduction.sh
- pexec.sh
- waitForSGEQJobs.pl (only for use with Sun Grid Engine)

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

Help
    exit 1
}

function setPath {
    cat <<setPath

--------------------------------------------------------------------------------------
Error locating ANTS
--------------------------------------------------------------------------------------
It seems that the ANTSPATH environment variable is not set. Please add the ANTSPATH
variable. This can be achieved by editing the .bash_profile in the home directory.
Add:

ANTSPATH=/home/yourname/bin/ants/

Or the correct location of the ANTS binaries.

Alternatively, edit this script (ants.sh) to set up this parameter correctly.

setPath
    exit 1
}

setFSLPath() {
    cat <<setFSLPath
--------------------------------------------------------------------------------------
Error locating FSL
--------------------------------------------------------------------------------------
The FSLDIR environment variable is not set. Please add the FSLDIR variable.

see the FSL website for more information about installation:

http://www.fmrib.ox.ac.uk/fsl/fsl/downloading.html

setFSLPath
    exit 1
}

jobfnamepadding() {

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

function reportMappingParameters {
    cat <<reportMappingParameters

--------------------------------------------------------------------------------------
 Mapping parameters
--------------------------------------------------------------------------------------
 ANTSPATH is $ANTSPATH

 Dimensionality:			$DIM
 N3BiasFieldCorrection:			$N3CORRECT
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
reportMappingParameters
}

function shapeupdatetotemplate {

#local declaration of values
dim=${DIM}
template=${TEMPLATE}
templatename=${TEMPLATENAME}
outputname=${OUTPUTNAME}
gradientstep=-${GRADIENTSTEP}

#debug only
# echo $dim
# echo ${template}
# echo ${templatename}
# echo ${outputname}
# echo ${outputname}*formed.nii*
# echo ${gradientstep}

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 1"
echo "--------------------------------------------------------------------------------------"
${ANTSPATH}AverageImages $dim ${template} 1 ${outputname}*formed.nii.gz

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 2"
echo "--------------------------------------------------------------------------------------"
if [ $dim -eq 2  ]
then
	${ANTSPATH}AverageImages $dim ${templatename}warpxvec.nii.gz 0 ${outputname}*Warpxvec.nii.gz
	${ANTSPATH}AverageImages $dim ${templatename}warpyvec.nii.gz 0 ${outputname}*Warpyvec.nii.gz

elif [ $dim -eq 3  ]
then
	${ANTSPATH}AverageImages $dim ${templatename}warpxvec.nii.gz 0 ${outputname}*Warpxvec.nii.gz
	${ANTSPATH}AverageImages $dim ${templatename}warpyvec.nii.gz 0 ${outputname}*Warpyvec.nii.gz
	${ANTSPATH}AverageImages $dim ${templatename}warpzvec.nii.gz 0 ${outputname}*Warpzvec.nii.gz
fi

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 3"
echo "--------------------------------------------------------------------------------------"
if [ $dim -eq 2  ]
then
	${ANTSPATH}MultiplyImages $dim ${templatename}warpxvec.nii.gz ${gradientstep} ${templatename}warpxvec.nii.gz
	${ANTSPATH}MultiplyImages $dim ${templatename}warpyvec.nii.gz ${gradientstep} ${templatename}warpyvec.nii.gz

elif [ $dim -eq 3  ]
then
	${ANTSPATH}MultiplyImages $dim ${templatename}warpxvec.nii.gz ${gradientstep} ${templatename}warpxvec.nii.gz
	${ANTSPATH}MultiplyImages $dim ${templatename}warpyvec.nii.gz ${gradientstep} ${templatename}warpyvec.nii.gz
	${ANTSPATH}MultiplyImages $dim ${templatename}warpzvec.nii.gz ${gradientstep} ${templatename}warpzvec.nii.gz
fi

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 4"
echo "--------------------------------------------------------------------------------------"
rm -f ${templatename}Affine.txt

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 5"
echo "--------------------------------------------------------------------------------------"
# Averaging and inversion code
if [ ${dim} -eq 2   ]
then
	ANTSAverage2DAffine ${templatename}Affine.txt ${outputname}*Affine.txt

	${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warpxvec.nii.gz ${templatename}warpxvec.nii.gz -i  ${templatename}Affine.txt -R ${template}
	${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warpyvec.nii.gz ${templatename}warpyvec.nii.gz -i  ${templatename}Affine.txt -R ${template}

	${ANTSPATH}WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}Affine.txt ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz -R ${template}

elif [ ${dim} -eq 3  ]
then
	ANTSAverage3DAffine ${templatename}Affine.txt ${outputname}*Affine.txt

	${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warpxvec.nii.gz ${templatename}warpxvec.nii.gz -i  ${templatename}Affine.txt -R ${template}
	${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warpyvec.nii.gz ${templatename}warpyvec.nii.gz -i  ${templatename}Affine.txt -R ${template}
	${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warpzvec.nii.gz ${templatename}warpzvec.nii.gz -i  ${templatename}Affine.txt -R ${template}

	${ANTSPATH}WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}Affine.txt ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz -R ${template}
fi

echo
echo "--------------------------------------------------------------------------------------"
echo " shapeupdatetotemplate 6"
echo "--------------------------------------------------------------------------------------"
if [ ${dim} -eq 2  ]
then
	${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warpxvec.nii.gz ${templatename}warpxlog.txt 1
	${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warpyvec.nii.gz ${templatename}warpylog.txt 1
elif [ ${dim} -eq 3  ]
then
	${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warpxvec.nii.gz ${templatename}warpxlog.txt 1
	${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warpyvec.nii.gz ${templatename}warpylog.txt 1
	${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warpzvec.nii.gz ${templatename}warpzlog.txt 1
fi

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

#initializing variables with global scope
time_start=`date +%s`
currentdir=`pwd`
nargs=$#

MAXITERATIONS=30x90x20
LABELIMAGE=0 # initialize optional parameter
METRICTYPE="PR" # initialize optional parameter
TRANSFORMATIONTYPE="GR" # initialize optional parameter
N3CORRECT=1 # initialize optional parameter
DOQSUB=0 # run locally by default
GRADIENTSTEP="0.25" # Gradient step size, smaller in magnitude means more smaller (more cautious) steps
ITERATIONLIMIT=4
CORES=2
TDIM=0
RIGID=0
RIGIDTYPE=" --do-rigid" # you can edit this to your liking; alternatively use --rigid-affine for registration with more (12?) DOF
range=0
REGTEMPLATE=target

# Provide different output for Usage and Help
if [ "$1" == "-h" ]
then
Help >&2
elif [ $nargs -lt 6 ]
then
Usage >&2
fi

# reading command line arguments
while getopts "c:d:i:j:h:m:n:o:s:r:t:z:" OPT
do
    case $OPT in
    h) #help
        echo "$USAGE"
        exit 0
        ;;
    c) #use SGE cluster
        DOQSUB=$OPTARG
	if [[ ${#DOQSUB} -gt 1 ]] ; then
	    echo " DOQSUB must be an integer value (0 or 1 ) you passed  -c $DOQSUB "
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
# will develop when I have time
#	testproc=`cat /proc/cpuinfo | grep processor | cut -d ':' -f 2 `
# 	if [[ ${#testproc} -gt ${CORES} ]] ;
# 	then
# 		echo " You have spefied more threads than that there are cpu cores in your system. "
# 		echo " This may slow down template generation"
# 		echo " Maximum number of threads is: ${#testproc}"
# 		echo " Press CTRL + C to exit and respecify "
# 		sleep 5
#
# 	elif [[ ${#testproc} -lt ${CORES} ]] ;
# 	then
# 		echo " You have spefied less threads than that there are cpu cores in your system. "
# 		echo " This may slow down template generation"
# 		echo " Maximum number of threads is: ${#testproc}"
# 	fi
        ;;
    m) #max iterations other than default
        MAXITERATIONS=$OPTARG
        ;;
    n) #apply bias field correction
        N3CORRECT=$OPTARG
        ;;
    o) #output name prefix
        OUTPUTNAME=$OPTARG
	TEMPLATENAME=${OUTPUTNAME}template
	TEMPLATE=${TEMPLATENAME}.nii.gz
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

#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
then
setPath >&2
fi

# Creating the file list of images to make a template from.
# Shiftsize is calculated because a variable amount of arguments can be used on the command line.
# The shiftsize variable will give the correct number of arguments to skip. Issueing shift $shiftsize will
# results in skipping that number of arguments on the command line, so that only the input images remain.
shiftsize=`expr $OPTIND - 1`
shift $shiftsize
# The invocation of $* will now read all remaining arguments into the variable IMAGESETVARIABLE
IMAGESETVARIABLE=$*
NINFILES=`expr $nargs - $shiftsize`

#test if FSL is available in case of 4D, exit if not
if [  ${TDIM} -eq 0 ] && [  ${#FSLDIR} -le 0 ]
then
setFSLPath >&2
fi

if [ ${NINFILES} -eq 0 ]
then
echo "Please provide at least 2 filenames for the template."
echo "Use $0 -h for help"
exit 1
elif [ ${NINFILES} -eq 1 ]
then
range=`fslnvols ${IMAGESETVARIABLE}`
	if [ ${range} -eq 1 ] && [ ${TDIM} -ne 4 ]
	then
		echo "Please provide at least 2 filenames for the template."
		echo "Use $0 -h for help"
		exit 1
	elif [ ${range} -gt 1 ] && [ ${TDIM} -ne 4 ]
	then
		echo "This is a multivolume file. Use -d 4"
		echo "Use $0 -h for help"
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
		fslsplit ${IMAGESETVARIABLE}
		rm ${IMAGESETVARIABLE}

		# selecting 16 volumes randomly from the timeseries for averaging, placing them in tmp/selection. folder
		for ((i = 0; i < 16 ; i++))
		do
			number=$RANDOM
			let "number %= $range"

			if [ ${number} -lt 10 ]
			then
			cp vol000${number}.nii.gz selection/
			elif [ ${number} -ge 10 ] && [ ${number} -lt 100 ]
			then
			cp vol00${number}.nii.gz selection/
			elif [ ${number} -ge 100 ] && [ ${number} -lt 1000 ]
			then
			cp vol0${number}.nii.gz selection/
			fi
		done

		# set filelist variable
		cd selection/
		IMAGESETVARIABLE=`ls *.nii.gz`

	fi
fi

# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS

# Test availability of helper scripts.
# No need to test this more than once. Can reside outside of the main loop.
ANTSSCRIPTNAME=${ANTSPATH}antsIntroduction.sh
PEXEC=${ANTSPATH}ANTSpexec.sh

for FLE in $ANTSSCRIPTNAME $PEXEC
do
	if [ ! -f $FLE  ] ;
	then
		echo
		echo "--------------------------------------------------------------------------------------"
		echo " FILE $FLE DOES NOT EXIST !!! $0 will terminate."
		echo "--------------------------------------------------------------------------------------"
		exit 1
	fi
done

# This script is not available in the current ANTS distribution.
# For SGE to work, please include this in the next ANTS update.

SGE=${ANTSPATH}waitForSGEQJobs.pl

if [ $DOQSUB -gt 0 ];
then
	for FLE in $SGE
	do
		if [ ! -f $FLE  ] ;
		then
			echo
			echo "--------------------------------------------------------------------------------------"
			echo " FILE $FLE DOES NOT EXIST !!! $0 will terminate."
			echo "--------------------------------------------------------------------------------------"
			exit 1
		fi
	done
fi

# check for an initial template image and perform rigid body registration if requested
if [ ! -s $REGTEMPLATE ] && [ $RIGID == 0 ];
then
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " No initial template exists. Creating a population average image from the inputs."
	echo "--------------------------------------------------------------------------------------"
	${ANTSPATH}AverageImages $DIM populationmean.nii.gz 1 $IMAGESETVARIABLE
	cp populationmean.nii.gz $TEMPLATE

elif [ ! -s $REGTEMPLATE ] && [ "$RIGID" == 1 ];
then
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " No initial template exists. Creating a population average image from the inputs."
	echo "--------------------------------------------------------------------------------------"
	${ANTSPATH}AverageImages $DIM populationmean.nii.gz 1 $IMAGESETVARIABLE
	cp populationmean.nii.gz $TEMPLATE

	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Rigid body registration to the population average image from the inputs."
	echo "--------------------------------------------------------------------------------------"

elif [ -s $REGTEMPLATE ] && [ "$RIGID" == 0 ] ;
then
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Initial template found.  This will be used for guiding the registration."
	echo "--------------------------------------------------------------------------------------"
	# now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.
	cp ${REGTEMPLATE} ${TEMPLATE}

elif [ -s $REGTEMPLATE ] && [ "$RIGID" == 1 ] ;
then
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Initial template found.  This will be used for guiding the registration."
	echo "--------------------------------------------------------------------------------------"
	# now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.
	cp ${REGTEMPLATE} ${TEMPLATE}
	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Rigid body registration to ${REGTEMPLATE} of all inputs."
	echo "--------------------------------------------------------------------------------------"
	# now move the initial registration template to OUTPUTNAME, otherwise this input gets overwritten.

fi

if [ "$RIGID" == 1 ] ;
then
	count=0
	for IMG in $IMAGESETVARIABLE
	do

		BASENAME=` echo ${IMG} | cut -d '.' -f 1 `

		exe=" ${ANTSPATH}ANTS $DIM -m  MI[${TEMPLATE},${IMG},1,32] -o rigid_${IMG} -i 0 --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000 $RIGIDTYPE"
		exe2="${ANTSPATH}WarpImageMultiTransform $DIM ${IMG} rigid_${IMG} rigid_${BASENAME}Affine.txt -R ${TEMPLATE}"
		pexe=" $exe >> job_${count}_metriclog.txt "

# this SGE call will go wrong, since the exe2 is not executed. Is there a way to sequentially execute exe and exe 2 without repeating the SGE loop? I assume that exe and exe2 need to be run on the same SGE node.
		if [ $DOQSUB -gt 0 ]; then
			id=`qsub -S /bin/bash -v ANTSPATH=$ANTSPATH $QSUBOPTS $exe | awk '{print $3}'`
			jobIDs="$jobIDs $id"
			sleep 0.3
		else
			# here comes the pexec call
			# sh $exe
			echo $pexe >> job${count}_r.sh
			echo $exe2 >> job${count}_r.sh
		fi

		((count++))
	done

	# Run jobs on SGE and wait to finish
	if [ $DOQSUB -gt 0 ];
	then
		echo " submitted $count jobs "
		# now wait for the stuff to finish; this script is absent
		${ANTSPATH}waitForSGEQJobs.pl 1 120 $jobIDs

		if [ ! $? -eq 0 ]; then
		echo "qsub submission failed - jobs went into error state"
		exit 1;
		fi

	fi

	# Run jobs on localhost and wait to finish
	if [ $DOQSUB -eq 0 ];
	then

		echo
		echo "--------------------------------------------------------------------------------------"
		echo " Starting ANTS rigid registration on max ${CORES} cpucores. "
		echo " Progress can be viewed in job*_metriclog.txt"
		echo "--------------------------------------------------------------------------------------"
		jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
		chmod +x job*.sh
		$PEXEC -j ${CORES} sh job*.sh
	fi

	# cleanup and save output in seperate folder
	rm job*.sh
	IMAGESETVARIABLE=`ls rigid_*.nii.gz`
	mkdir rigid
	cp rigid_*.nii.gz rigid_*.cfg *Affine.txt rigid/

fi



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
	rm -f  ${OUTPUTNAME}*warp*nii*

	# iteration 1
	if [  $i -eq 0 ]
	then
		for (( n = 0 ; n < ${NUMLEVELS}; n++ ))
		do
			val=0;
			if [ $n  -eq  0  ]
			then
				val=${ITERATLEVEL}
				ITERATIONS=$val;
			fi

			if [ $n  -gt  0  ]
			then
				ITERATIONS=${ITERATIONS}x${val}
			fi
		done
	fi

	# iteration 2
	if [  $i -eq 1 ]
	then
		for (( n = 0; n<${NUMLEVELS}; n++ ));
		do
			val=0;
			if [ $n  -eq  0  ]
			then
				val=${ITERATLEVEL}
				ITERATIONS=$val
			fi

			if [ $n  -eq  1  ]
			then
				val=${ITERATLEVEL[1]}
				ITERATIONS=${ITERATIONS}x${val}
			fi

			if [ $n  -gt  1  ]
			then
				ITERATIONS=${ITERATIONS}x${val}
			fi
		done
	fi

	# iteration 3 and further
	if [  $i -gt 1 ]
	then
		ITERATIONS=$MAXITERATIONS
	fi

	# Job IDs of jobs submitted to queue in loop below
	jobIDs=""

	# Reinitialize count to 0
	count=0

	# Submit registration of each input to volume template to SGE or run locally.
	for IMG in $IMAGESETVARIABLE
	do
		#1 determine working dir
		dir=`pwd`

		#2 determine new filename
		POO=${OUTPUTNAME}${IMG}


		#3 Make variable OUTFILENAME and remove anything behind . ; for example .nii.gz.gz
		OUTFN=${POO%.*.*}

		#4 Test if outputfilename has only a single extention and remove that
		if [ ${#OUTFN} -eq ${#POO} ]
		then
			OUTFN=${OUTPUTNAME}${IMG%.*}
		fi

		#5 prepare registration command
		exe="${ANTSSCRIPTNAME} -d ${DIM} -r ${dir}/${TEMPLATE} -i ${dir}/${IMG} -o ${dir}/${OUTFN} -m ${MAXITERATIONS} -n ${N3CORRECT} -s ${METRICTYPE} -t ${TRANSFORMATIONTYPE} "
                pexe=" $exe >> job_${count}_${i}_metriclog.txt "

		#6 submit to SGE or else run locally
		if [ $DOQSUB -gt 0 ]; then
			id=`qsub -S /bin/bash -v ANTSPATH=$ANTSPATH $QSUBOPTS $exe | awk '{print $3}'`
			jobIDs="$jobIDs $id"
			sleep 0.3
		else
			# here comes the pexec call
			# sh $exe
			echo $pexe >> job${count}_${i}.sh
		fi

		# counter updated, but not directly used in this loop
		count=`expr $count + 1`;
#		echo " submitting job number $count " # for debugging only
	done
	# SGE wait for script to finish
	if [ $DOQSUB -gt 0 ];
	then
		echo " submitted $count jobs "
		# now wait for the stuff to finish; this script is absent
		${ANTSPATH}waitForSGEQJobs.pl 1 120 $jobIDs

		if [ ! $? -eq 0 ]; then
		echo "qsub submission failed - jobs went into error state"
		exit 1;
		fi

	fi

#	echo "finished $i " >> ${TEMPLATENAME}metriclog.txt


	# Run jobs on localhost and wait to finish
	if [ $DOQSUB -eq 0 ];
	then
		itdisplay=`expr $i + 1`;
		echo
		echo "--------------------------------------------------------------------------------------"
		echo " Starting ANTS registration on max ${CORES} cpucores. Iteration: $itdisplay of $ITERATIONLIMIT"
		echo " Progress can be viewed in job*_${i}_metriclog.txt"
		echo "--------------------------------------------------------------------------------------"
		jobfnamepadding #adds leading zeros to the jobnames, so they are carried out chronologically
		chmod +x job*.sh
		$PEXEC -j ${CORES} sh job*_${i}.sh
	fi

	shapeupdatetotemplate ${DIM} ${TEMPLATE} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP}

	if [ $DIM -eq 2  ]
	then
		${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpxvec.nii.gz ${TEMPLATENAME}warpxlog.txt 1
		${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpyvec.nii.gz ${TEMPLATENAME}warpylog.txt 1

	elif [ $DIM -eq 3  ]
	then
		${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpxvec.nii.gz ${TEMPLATENAME}warpxlog.txt 1
		${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpyvec.nii.gz ${TEMPLATENAME}warpylog.txt 1
		${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpzvec.nii.gz ${TEMPLATENAME}warpzlog.txt 1
	fi

	echo
	echo "--------------------------------------------------------------------------------------"
	echo " Backing up results from iteration $itdisplay"
	echo "--------------------------------------------------------------------------------------"

	mkdir ${TRANSFORMATIONTYPE}_iteration_${i}
	cp j* *.cfg *${OUTPUTNAME}*.nii.gz ${TRANSFORMATIONTYPE}_iteration_${i}

	((i++))

	done

# end main loop

rm job*.sh

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
