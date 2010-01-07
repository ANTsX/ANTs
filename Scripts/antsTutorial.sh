# !/bin/sh

VERSION="0.0.2"

function Usage {
    cat <<USAGE

ants.sh will perform an elastic or diffeomorphic transformation of the moving image
into the fixed image space.

this script encodes some reasonable defaults for the parameters and is useful for testing and comparing methods

Usage:

sh ants.sh -d ImageDimension -r fixed.ext -i moving.ext

Compulsory arguments:

-d:  ImageDimension: 2 or 3 (for 2 or 3 Dimensional registration)

-r:  reference image

-i:  input image

Optional arguments

-o:  OUTPREFIX; A prefix that is prepended to all output files.

-m:  max-iterations

	Max-Iterations in form: JxKxL where
	J = max iterations at coarsest resolution (here, reduce by power of 2^2)
	K = middle resolution iterations (here,reduce by power of 2)
	L = fine resolution iterations (here, full resolution) !!this level takes much
            more time per iteration!!

	Adding an extra value before JxKxL (i.e. resulting in IxJxKxL) would add another
	iteration level.

-n3  N3BiasFieldCorrection of moving image ( 0 = off; 1 = on (default) )

-l:  Labels-In-Fixed-Image-Space-To-Deform-To-Moving-Image (default is off)

-tm: type of similarity metric used for registration.

	For intramodal image registration, use:
	CC = cross-correlation
	MI = mutual information
	PR = probability mapping (default)
	MSQ = mean square difference

	For intermodal image registration, use:
	MI = mutual information
	PR = probability mapping (default)

	You can change the values for each of these methods in this script.

-tt: type of transformation model used for registration.

	For elastic image registration, use:
	EL = elastic transformation model (less deformation possible)

	For diffeomorphic image registration, use:
	SY = SyN with time (default) with arbitrary number of time points in time discretization
	S2 = SyN with time optimized specifically for 2 time points in the time discretization
	GR = Greedy SyN
	EX = Exponential
        DD = Diffeomorphic Demons style exponential mapping

	You can change the values for each of these methods in this script.

-q:  Perform a Quality Check (QC) of the result ( 0 = off (default); 1 = on ).
--------------------------------------------------------------------------------------

Default number of iterations is $MAXITERATIONS. Most often this is sufficient.

The other parameters in this script are updates to the defaults that were used in the
Neuroimage 2009 paper; for more information read:

- Klein A, et al. 2009. Evaluation of 14 nonlinear deformation algorithms applied to
human brain MRI registration. NeuroImage, 46(3):786-802

http://www.ncbi.nlm.nih.gov/pubmed/19195496

--------------------------------------------------------------------------------------
I/O Data Formats In ANTS
--------------------------------------------------------------------------------------
see: http://picsl.upenn.edu/ANTS/ioants.php

--------------------------------------------------------------------------------------
Get the latest ANTS version at:
http://sourceforge.net/projects/advants/

Read the ANTS documentation at:
http://picsl.upenn.edu/ANTS/
--------------------------------------------------------------------------------------
ANTS was created by:
Brian B. Avants, Nick Tustison and Gang Song
Penn Image Computing And Science Laboratory
University of Pennsylvania

script adapted by N.M. van Strien, www.mri-tutorial.com
--------------------------------------------------------------------------------------

USAGE
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

function reportMappingParameters {
    cat <<reportMappingParameters
--------------------------------------------------------------------------------------
Mapping parameters
--------------------------------------------------------------------------------------
ANTSPATH is $ANTSPATH

Dimensionality:				$DIM
Fixed image:				$FIXED
Moving image:				$MOVING
Outputname prefix:			$OUTPUTNAME
Use label image:			$LABELIMAGE
N3BiasFieldCorrection:			$N3CORRECT
DoANTS Quality Check: 			$DoANTSQC
Similarity Metric:			$METRICTYPE
Transformation:				$TRANSFORMATIONTYPE
Regularization:				$REGULARIZATION
MaxIterations:				$MAXITERATIONS
Number Of MultiResolution Levels:	$NUMLEVELS
OutputName:				$OUTPUTNAME
--------------------------------------------------------------------------------------
reportMappingParameters
}

function dataCheck {
    cat <<dataCheck
--------------------------------------------------------------------------------------
Data integrity check failed
--------------------------------------------------------------------------------------
There seems to be a problem with the header definitions of the input files. This
script has tried to repair this issue automatically by invoking:

${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}repaired.nii CompareHeadersAndImages $FIXED $MOVING

The repaired image is:

${OUTPUTNAME}repaired.nii

and should be located in:

$currentdir

Check if the corrected file is OK and use this as input to the ants.sh script. Otherwise
call ImageMath's CompareHeadersAndImages on your Fixed and Moving image, to troubleshoot
further.
--------------------------------------------------------------------------------------

dataCheck
    exit 1
}

currentdir=`pwd`
nargs=$#

MAXITERATIONS=30x90x20
LABELIMAGE=0 # initialize optional parameter
DoANTSQC=0 # initialize optional parameter
METRICTYPE="PR" # initialize optional parameter
TRANSFORMATIONTYPE="GR" # initialize optional parameter
N3CORRECT=1 # initialize optional parameter

if [ $nargs -lt 3 ]
then
Usage >&2
fi

# reading command line arguments
while getopts "d:r:i:h:o:m:n:l:q:s:t:" OPT
do
    case $OPT in
    h) #help
        echo "$USAGE"
        exit 0
        ;;
    d) #dimensions
        DIM=$OPTARG
        ;;
    i) #input or moving image
        MOVING=$OPTARG
        OUTPUTNAME=` echo $MOVING | cut -d '.' -f 1 `
        ;;
    l) #use label image
        LABELIMAGE=$OPTARG
        ;;
    m) #max iterations other than default
        MAXITERATIONS=$OPTARG
        ;;
    n) #apply bias field correction
        N3CORRECT=$OPTARG
        ;;
    o) #output name prefix
        OUTPUTNAME=$OPTARG
        ;;
    q) #perform quality check after running ANTS
        DoANTSQC=$OPTARG
        ;;
    r) #ref image
        FIXED=$OPTARG
        ;;
    s) #similarity model
        METRICTYPE=$OPTARG
        ;;
    t) #transformation model
        TRANSFORMATIONTYPE=$OPTARG
        ;;
    \?) # getopts issues an error message
        echo "$USAGE" >&2
        exit 1
        ;;
    esac
done

export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS
#ANTSPATH=YOURANTSPATH
if [  ${#ANTSPATH} -le 0 ]
then
setPath >&2
fi

# Uncomment the line below in case you have not set the ANTSPATH variable in your environment.
# ANTSPATH=

# test input parameter before starting
if  [ ${#DIM} -gt 1 ]
then
echo "Problem with specified ImageDimension => User Specified Value = $DIM"
exit
fi

if  [ ${#FIXED} -lt 1 -o  ! -f $FIXED ]
then
echo "Problem with specified Fixed Image => User Specified Value = $FIXED"
exit
fi

if  [ ${#MOVING} -lt 1 -o  ! -f $MOVING ]
then
echo "Problem with specified Moving Image => User Specified Value = $MOVING"
exit
fi

# check the image headers
echo "input files: checking"
compareheaders=`${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}repaired.nii CompareHeadersAndImages $FIXED $MOVING  | grep FailureState | cut -d ' ' -f 4 `

if [ $compareheaders -ne 0 ]
then
dataCheck
else
echo "input files: check passed"
fi

ITERATLEVEL=(`echo $MAXITERATIONS | tr 'x' ' '`)
NUMLEVELS=${#ITERATLEVEL[@]}

# Transformation model
if [ "${TRANSFORMATIONTYPE}" == "EL" ]
then
# Mapping Parameters
TRANSFORMATION=Elast[1]
REGULARIZATION=Gauss[0,3]

elif [ "${TRANSFORMATIONTYPE}" == "S2"  ]
then
# Mapping Parameters for the LDDMM style SyN --- the params are SyN[GradientStepLength,NTimeDiscretizationPoints,IntegrationTimeStep]
# increasing IntegrationTimeStep increases accuracy in the diffeomorphism integration and takes more computation time.
# NTimeDiscretizationPoints is set to 2 here
TRANSFORMATION=SyN[1,2,0.05]
REGULARIZATION=Gauss[3,0.]

elif [ "${TRANSFORMATIONTYPE}" == "SY"  ]
then
# Mapping Parameters for the LDDMM style SyN --- the params are SyN[GradientStepLength,NTimeDiscretizationPoints,IntegrationTimeStep]
# increasing IntegrationTimeStep increases accuracy in the diffeomorphism integration and takes more computation time.
# NTimeDiscretizationPoints is the number of spatial indices in the time dimension (the 4th dim when doing 3D registration)
# increasing NTimeDiscretizationPoints increases flexibility and takes more computation time.
# the --geodesic option enables either 1 asymmetric gradient estimation or 2 symmetric gradient estimation (the default here )
TRANSFORMATION=" SyN[1,2,0.05] --geodesic 2 "
REGULARIZATION=Gauss[3,0.25]

elif [ "${TRANSFORMATIONTYPE}" == "GR" ]
then
# Mapping Parameters for the greedy gradient descent (fast) version of SyN -- only needs GradientStepLength
TRANSFORMATION=SyN[0.25]
REGULARIZATION=Gauss[3,0]

elif [ "${TRANSFORMATIONTYPE}" == "EX" ]
then
# Mapping Parameters
TRANSFORMATION=Exp[0.5,10]
REGULARIZATION=Gauss[3,0.5]

elif [ "${TRANSFORMATIONTYPE}" == "DD" ]
then
# Mapping Parameters for diffemorphic demons style optimization Exp[GradientStepLength,NTimePointsInIntegration]
#  NTimePointsInIntegration controls the number of compositions in the transformation update , see the DD paper
TRANSFORMATION=GreedyExp[0.5,10]
REGULARIZATION=Gauss[3,0.5]
else
echo "Invalid transformation metric. Use EL, SY, S2, GR , DD or EX or type sh ants.sh -h."
exit 1
fi

# Similarity Measures
if [ "${METRICTYPE}" == "PR" ]
then
# Mapping Parameters
METRIC=PR[
METRICPARAMS=1,4]

elif [ "${METRICTYPE}" == "CC"  ]
then
# Mapping Parameters
METRIC=CC[
METRICPARAMS=1,5]

elif [ "${METRICTYPE}" == "MI" ]
then
# Mapping Parameters
METRIC=MI[
METRICPARAMS=1,32]

elif [ "${METRICTYPE}" == "MSQ" ]
then
# Mapping Parameters
METRIC=MSQ[
METRICPARAMS=1,0]

else
echo "Invalid similarity metric. Use CC, MI, MSQ or PR or type sh ants.sh -h."
exit 1
fi

# main script
reportMappingParameters

if  [ ${N3CORRECT} -eq 0 ]
then
# Apply ANTS mapping command without N3BiasFieldCorrection
exe="${ANTSPATH}ANTS $DIM -m  ${METRIC}${FIXED},${MOVING},${METRICPARAMS} -t $TRANSFORMATION -r $REGULARIZATION -o ${OUTPUTNAME} -i $MAXITERATIONS --use-Histogram-Matching "
echo
echo "--------------------------------------------------------------------------------------"
echo "ANTS command:"
echo "$exe "
echo "--------------------------------------------------------------------------------------"
$exe

# Apply forward transformation to MOVING
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying forward transformation to ${MOVING}"
echo "--------------------------------------------------------------------------------------"
${ANTSPATH}WarpImageMultiTransform $DIM ${MOVING} ${OUTPUTNAME}deformed.nii ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt -R ${FIXED}

# Apply inverse transformation to FIXED
if [ "${TRANSFORMATIONTYPE}" == "EL" ]
then
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying inverse transformation to ${FIXED} is not possible for elastic registration."
echo "--------------------------------------------------------------------------------------"
else
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying inverse transformation to ${FIXED}"
echo "--------------------------------------------------------------------------------------"
FIXEDBASE=` echo $MOVING | cut -d '.' -f 1 `
${ANTSPATH}WarpImageMultiTransform $DIM ${FIXED} ${FIXEDBASE}_InverseWarp.nii -R ${MOVING} -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii
fi

elif [ ${N3CORRECT} -eq 1 ]
then
# Apply N3BiasFieldCorrection
exe="${ANTSPATH}N3BiasFieldCorrection $DIM $MOVING ${OUTPUTNAME}.nii 4"
echo
echo "--------------------------------------------------------------------------------------"
echo "N3BiasFieldCorrection command:"
echo "$exe "
echo "--------------------------------------------------------------------------------------"
$exe

# Apply ANTS mapping command on N3 corrected image
exe="${ANTSPATH}ANTS $DIM -m  ${METRIC}${FIXED},${OUTPUTNAME}.nii,${METRICPARAMS} -t $TRANSFORMATION -r $REGULARIZATION -o ${OUTPUTNAME} -i $MAXITERATIONS --use-Histogram-Matching "
echo
echo "--------------------------------------------------------------------------------------"
echo "ANTS command:"
echo "$exe "
echo "--------------------------------------------------------------------------------------"
$exe

# Apply forward transformation to MOVING
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying forward transformation to ${MOVING}"
echo "--------------------------------------------------------------------------------------"
${ANTSPATH}WarpImageMultiTransform $DIM ${OUTPUTNAME}.nii ${OUTPUTNAME}deformed.nii ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt -R ${FIXED}

# Apply inverse transformation to FIXED
if [ "${TRANSFORMATIONTYPE}" == "EL" ]
then
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying inverse transformation to ${FIXED} is not possible for elastic registration."
echo "--------------------------------------------------------------------------------------"
else
echo
echo "--------------------------------------------------------------------------------------"
echo "Applying inverse transformation to ${FIXED}"
echo "--------------------------------------------------------------------------------------"
FIXEDBASE=` echo $MOVING | cut -d '.' -f 1 `
${ANTSPATH}WarpImageMultiTransform $DIM ${FIXED} ${FIXEDBASE}_InverseWarp.nii -R ${OUTPUTNAME}.nii -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii
fi

fi

# Apply transformation on labeled image?
if  [ ${#LABELIMAGE} -gt 3 ]
then
${ANTSPATH}WarpImageMultiTransform $DIM  $LABELIMAGE ${OUTPUTNAME}labeled.nii -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii  -R ${MOVING}   --use-NN
fi

if [ $DoANTSQC -eq 1 ]
then
#  measure image similarity - 0 = intensity difference, 1 = correlation , 2 = mutual information
for SIM in 0 1 2 ; do
${ANTSPATH}MeasureImageSimilarity $DIM $SIM $FIXED ${OUTPUTNAME}deformed.nii
done

# the lines below  measure dice overlap and min-dist-sum from the approximately segmented brain (3 tissues and background)
${ANTSPATH}ThresholdImage $DIM $FIXED ${OUTPUTNAME}fixthresh.nii.gz Otsu 4
${ANTSPATH}ThresholdImage $DIM $MOVING ${OUTPUTNAME}movthresh.nii.gz Otsu 4

${ANTSPATH}WarpImageMultiTransform $DIM ${OUTPUTNAME}movthresh.nii.gz ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt -R ${FIXED} --use-NN

${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}dicestats.txt DiceAndMinDistSum  ${OUTPUTNAME}fixthresh.nii.gz   ${OUTPUTNAME}movthresh.nii.gz ${OUTPUTNAME}mindistsum.nii.gz

# below not used unless you want to compare segmentation volumes to those estimated by the jacobian
# labelstats for jacobian wrt segmenation below, to compose
# ${ANTSPATH}ComposeMultiTransform $DIM   ${OUTPUTNAME}CompWarp.nii  -R $FIXED ${OUTPUTNAME}Warp.nii ${OUTPUTNAME}Affine.txt
# ${ANTSPATH}CreateJacobianDeterminantImage $DIM ${OUTPUTNAME}CompWarp.nii ${OUTPUTNAME}jacobian.nii  0
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}movlabstat.txt LabelStats ${OUTPUTNAME}movthresh.nii.gz ${OUTPUTNAME}movthresh.nii.gz
# ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}jaclabstat.txt LabelStats ${OUTPUTNAME}defthresh.nii.gz ${OUTPUTNAME}jacobian.nii
# we compare the output of these last two lines:
#  the Volume of the movlabstat computation vs. the mass of the jaclabstat
fi

exit