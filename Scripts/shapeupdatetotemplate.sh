#!/bin/sh

NUMPARAMS=$#

if [ $NUMPARAMS -lt 4  ]
then
echo " USAGE ::  "
echo "   sh $0  ImageDimension OutputRoot GradientStep <images>"
echo " example : "
echo " sh $0  2 TEST 0.5 YFace*nii "
echo " above, dimension is 2,  TEST is the OutputRoot ,  1 is the grad-step and the images you originally mapped are in the list YFace<var>.nii "
echo " we assume you used ants.sh naming conventions to deform your images and the deformations are in the same directory as the template and the call to the script . "
echo " "
echo " ImageDimension  -  Dimension of your image, eg 3 for 3D."
echo " OutputRoot      -  Root file name of the output. Should be a file root only, no path information. "
echo " GradientStep -   the size of the shape update gradient = 0.25 is typical "
echo " <images>        -  List of images in the current directory, eg *_t1.nii.gz. "
echo
echo " This script performs a shape update to the template assuming that  a registration to the template exists. "
echo
echo " We assume all files to be added to the template are in the current directory. You can modify the "
echo " script if you want to relax this assumption, but be sure that the qsubbed jobs get the correct "
echo " absolute path to the images."
echo
echo " Things within the script that you may need to change for your needs are highlighted by EDIT THIS "
echo
echo " The template will be written to [OutputRoot]template.nii. If the template file exists, it is used as the starting point for "
echo " the new template creation. Otherwise, we create an unbiased starting point by averaging the input dataset. "
exit
fi

#initialization, here, is unbiased
DIM=$1

# Root of the output name, will produce ${OUTPUTNAME}template.nii
# If this already exists, it will be used to initialize the template building
OUTPUTNAME=$2

# ANTSPATH - you will need to edit this if it is not set before the script is called
# note trailing slash - this is needed
export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS

# System specific queue options, eg "-q name" to submit to a specific queue
# It can be set to an empty string if you do not need any special cluster options
QSUBOPTS="" # EDIT THIS

  TEMPLATENAME=${OUTPUTNAME}template
  TEMPLATE=${TEMPLATENAME}.nii
  # Gradient step size, smaller in magnitude means more smaller (more cautious) steps
  GRADIENTSTEP=-${3}

shift 3

# Optionally disable qsub for debugging purposes - runs jobs in series
DOQSUB=0


IMAGESETVARIABLE=$*

if [ ! -s $TEMPLATE ] ; then
echo " No initial template exists. Cannot update the template."
exit
fi

echo  " ANTSPATH  $ANTSPATH "
echo " OutputName :  $OUTPUTNAME "
echo " template  $TEMPLATE "
echo " Template Update Steps $ITERATIONLIMIT "
echo " Template population :   $IMAGESETVARIABLE "
echo " grad step -- $GRADIENTSTEP "
echo " if the files and parameters are all ok then uncomment the exit call below this line  "
echo " "
#exit

   rm -f  ${OUTPUTNAME}*InverseWarp*vec.nii

deformedimages=(` ls ${OUTPUTNAME}*formed.nii*     `)
NUM=${#deformedimages[@]}
if [ $NUM -le 1 ] ; then
echo " you do not have images of the type :  ${OUTPUTNAME}*formed.nii* "
echo " if they dont exist, you need to run ants.sh to map your input images $IMAGESETVARIABLE to the template : $TEMPLATE "
exit
fi
deformx=(` ls ${OUTPUTNAME}*Warpxvec.nii     `)
if [ ${#deformx[@]} -le 1 ] ; then
echo " you do not have images of the type :  ${OUTPUTNAME}*Warpxvec.nii  "
echo " if they dont exist, you need to run ants.sh to map your input images $IMAGESETVARIABLE to the template : $TEMPLATE "
exit
fi
deformy=(` ls ${OUTPUTNAME}*Warpyvec.nii     `)
if [ ${#deformy[@]} -le 1 ] ; then
echo " you do not have images of the type :  ${OUTPUTNAME}*Warpyvec.nii  "
echo " if they dont exist, you need to run ants.sh to map your input images $IMAGESETVARIABLE to the template : $TEMPLATE "
exit
fi
if [ $DIM -eq 3 ] ; then
deformz=(` ls ${OUTPUTNAME}*Warpzvec.nii     `)
if [ ${#deformz[@]} -le 1 ] ; then
echo " you do not have images of the type :  ${OUTPUTNAME}*Warpzvec.nii  "
echo " if they dont exist, you need to run ants.sh to map your input images $IMAGESETVARIABLE to the template : $TEMPLATE "
exit
fi
fi


${ANTSPATH}AverageImages $DIM ${TEMPLATE} 1 ${OUTPUTNAME}*formed.nii*
#sh sygnccavg.sh 0.1  $TEMPLATE   # uncomment this for sygn template.

# below, a cheap approach to integrating the negative velocity field
# in the absence of other code and saving the velocity fields.
# additionally, this works for all types of registration algorithms

     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpxvec.nii 0 ${OUTPUTNAME}*Warpxvec.nii
     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpyvec.nii 0 ${OUTPUTNAME}*Warpyvec.nii
if [ $DIM -gt 2  ]
then
     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpzvec.nii 0 ${OUTPUTNAME}*Warpzvec.nii
fi
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpxvec.nii $GRADIENTSTEP ${TEMPLATENAME}warpxvec.nii
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpyvec.nii $GRADIENTSTEP  ${TEMPLATENAME}warpyvec.nii
if [ $DIM -gt 2  ]
then
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpzvec.nii $GRADIENTSTEP  ${TEMPLATENAME}warpzvec.nii
fi
    ${ANTSPATH}WarpImageMultiTransform $DIM  ${TEMPLATE}   ${TEMPLATE} ${TEMPLATENAME}warp.nii ${TEMPLATENAME}warp.nii ${TEMPLATENAME}warp.nii  ${TEMPLATENAME}warp.nii  -R ${TEMPLATE}
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpxvec.nii  ${TEMPLATENAME}warpxlog.txt  1
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpyvec.nii  ${TEMPLATENAME}warpylog.txt  1
if [ $DIM -gt 2  ]
then
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpzvec.nii  ${TEMPLATENAME}warpzlog.txt  1
fi
