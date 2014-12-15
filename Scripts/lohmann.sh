#!/bin/bash

NUMPARAMS=$#
MAXITERATIONS=30x90x20

if [ $NUMPARAMS -lt 2  ]
then
echo " USAGE ::  "
echo "  sh   $0 ImageDimension  image.ext   "
echo " be sure to set ANTSPATH environment variable "
exit
fi

ANTSPATH=/mnt/aibs1/avants/bin/ants/
if [  ${#ANTSPATH} -le 0 ]
then
echo " Please set ANTSPATH=LocationOfYourAntsBinaries "
echo " Either set this in your profile or directly, here in the script. "
echo " For example : "
echo " ANTSPATH=/home/yourname/bin/ants/ "
exit
else
echo " ANTSPATH is $ANTSPATH "
fi

#initialization, here, is unbiased
DIM=$1

if  [ ${#DIM} -gt 1 ]
then
echo " Problem with specified ImageDimension => User Specified Value = $DIM "
exit
fi

FIXED=$2

if  [ ${#FIXED} -lt 1 -o  ! -f $FIXED ]
then
echo " Problem with specified Fixed Image => User Specified Value = $FIXED "
exit
fi

OUTPUTNAME=` echo $FIXED | cut -d '.' -f 1 `

if [[ ! -s ${OUTPUTNAME}repaired.nii.gz ]] ; then
${ANTSPATH}/N3BiasFieldCorrection $DIM $FIXED ${OUTPUTNAME}repaired.nii.gz 4
fi
${ANTSPATH}/ThresholdImage $DIM  ${OUTPUTNAME}repaired.nii.gz ${OUTPUTNAME}mask.nii.gz 0.1 99999

WM=${OUTPUTNAME}seg.nii.gz
if [[ ! -s $WM ]] ; then
  ${ANTSPATH}/Apocrita $DIM -x ${OUTPUTNAME}mask.nii.gz  -m [0.5,1,0,0] -o [${OUTPUTNAME}seg.nii.gz] -i Otsu[${OUTPUTNAME}repaired.nii.gz,3] -h 0
  ${ANTSPATH}/ThresholdImage $DIM $WM $WM 3 3
  ${ANTSPATH}/ImageMath $DIM $WM GetLargestComponent $WM
fi
OUT=${OUTPUTNAME}
if [[ ! -s  ${OUT}max.nii.gz ]] ; then
# the closing operation - arachnoid surface
${ANTSPATH}/ImageMath 3 ${OUT}max.nii.gz MD $WM  15
${ANTSPATH}/ImageMath 3 ${OUT}max.nii.gz ME ${OUT}max.nii.gz  15

# unfolded or lissencephalic surface
# we would prefer to do this in a topology preserving way
${ANTSPATH}/ImageMath 3 ${OUT}min.nii.gz ME $WM 5
${ANTSPATH}/ImageMath 3 ${OUT}min.nii.gz MD ${OUT}min.nii.gz 5

# distance transform from arachnoid surface to inside of brain
# perhaps should be a geodesic distance transform
${ANTSPATH}/ImageMath 3 ${OUT}max.nii.gz Neg ${OUT}max.nii
${ANTSPATH}/ImageMath 3 ${OUT}dist.nii.gz D ${OUT}max.nii
fi
# multiply the distance transform by the "unfolded" surface
# this gives the "feature" surface/image that could be input
# to an image registration similarity metric
${ANTSPATH}/ImageMath 3 ${OUT}lohmann.nii.gz m ${OUT}dist.nii.gz ${OUT}min.nii.gz

# surface stuff for visualization
${ANTSPATH}/ImageMath 3 ${OUT}surf.nii.gz ME ${OUT}min.nii.gz 2
${ANTSPATH}/ImageMath 3 ${OUT}surf.nii.gz - ${OUT}min.nii.gz ${OUT}surf.nii.gz
${ANTSPATH}/ImageMath 3 ${OUT}maxvals.nii.gz m ${OUT}lohmann.nii.gz ${OUT}surf.nii.gz

# find sulci using curvature of the WM/lohmann surface
	${ANTSPATH}/SurfaceCurvature ${OUT}lohmann.nii.gz ${OUT}lohmannk.nii.gz 1.5

${ANTSPATH}/ThresholdImage 3 ${OUT}lohmann.nii.gz ${OUT}lohmann.nii.gz 0.001 9999
${ANTSPATH}/ImageMath 3 ${OUT}lohmannk.nii.gz Normalize ${OUT}lohmannk.nii.gz
${ANTSPATH}/ImageMath 3 ${OUT}lohmannk.nii.gz m ${OUT}lohmannk.nii.gz ${OUT}lohmann.nii.gz


exit
