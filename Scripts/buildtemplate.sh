#!/bin/sh

NUMPARAMS=$#

if [ $NUMPARAMS -lt 4  ]
then
echo " USAGE ::  "
echo "  sh   buildtemplate.sh  ImageDimension WildcardingSide1 WildcardingSide2 Naming  OptionalIterationLimit-Default=4  OptionalDoQSub "
echo " "
echo " the files used for template construction will  be  of the type     WildcardingSide1*WildcardingSide2  "
echo " "
echo " We assume all files to be added to the template are in the same directory -- you can modify the script if you prefer something different "
echo "  "
echo " things you should change for your needs are highlighted by   EDIT THIS "
echo " if the template file does not exist yet, we create an unbiased starting point by averaging the input dataset "
exit
fi

#initialization, here, is unbiased
DIM=$1

OUTPUTNAME=''
if [ $NUMPARAMS -gt 2  ]
then
OUTPUTNAME=$4
fi

ANTSPATH="/mnt/aibs1/avants/bin/ants/"  # EDIT THIS -- location of ANTS programs
# Mapping Parameters
  TRANSFORMATION=SyN[0.25]   # EDIT THIS  -- transformation model
  MAXITERATIONS=50x10x0    # EDIT THIS -- max iterations per pairwise registration
  ITERATLEVEL=(`echo $MAXITERATIONS | tr 'x' ' '`)
  NUMLEVELS=${#ITERATLEVEL[@]}
  echo $NUMLEVELS
  REGULARIZATION=Gauss[3,0]   # EDIT THIS  --- regularization of the mapping
#SEE ALSO!!! -- METRIC variable Below!!!
  TEMPLATENAME=${OUTPUTNAME}template
  TEMPLATE=${TEMPLATENAME}.nii
  METRIC=PR[${TEMPLATE}    # EDIT THIS  -- similarity metric
  METRICPARAMS=1,2]    # EDIT THIS

IMAGESETVARIABLE=${2}*${3}

ITERATIONLIMIT=4
if [ $NUMPARAMS -gt 4  ]
then
ITERATIONLIMIT=$5
fi

DOQSUB=0
if [ $NUMPARAMS -gt 5 ]
then
DOQSUB=$6
fi

LONG=0
count=0
for x in `ls -tr  $IMAGESETVARIABLE  `
do
echo " Image $count is:   $x  "
BASE=${x%.*}
NAMING=${OUTPUTNAME}${BASE}
echo " Base Name is $BASE "
echo " Full Naming is  $NAMING "
echo " "
if [ $count -eq 0 ] && [ $LONG -eq 1 ]
then
echo " Using biased initialization from first image! "
 cp  $x  ${TEMPLATE}
fi

count=`expr $count + 1`
done

if [  $LONG -eq 0  ]
then
if [ ! -s $TEMPLATE ] ; then
echo " Averaging :   $IMAGESETVARIABLE "
echo " "
 ${ANTSPATH}AverageImages $DIM $TEMPLATE  $IMAGESETVARIABLE
fi
fi

echo  " ANTSPATH  $ANTSPATH "
echo " Mapping Parameters  :: "
echo  " Transformation is:  $TRANSFORMATION "
echo " MaxIterations :   $MAXITERATIONS "
echo " Number Of MultiResolution Levels   $NUMLEVELS "
echo " Regularization :  $REGULARIZATION "
echo " Metric :  ${METRIC},File,${METRICPARAMS} "
echo " OutputName :  $OUTPUTNAME "
echo " The template will be called  $TEMPLATE "
echo " Template Update Steps $ITERATIONLIMIT "
echo " Averaging :   $IMAGESETVARIABLE "
echo " "
echo " if the files and parameters are all ok then uncomment the exit call below this line  "
echo " "
# exit



#${ANTSPATH}AverageImages $DIM ${TEMPLATE} 0 $IMAGESETVARIABLE
#ConvertToJpg template.nii template0.jpg

\rm metriclog.txt ${TEMPLATENAME}warpxlog.txt  ${TEMPLATENAME}warpylog.txt   ${TEMPLATENAME}warpzlog.txt

# begin loop
i=0

while [  $i -lt ${ITERATIONLIMIT} ]
  do
  rm -f  ${OUTPUTNAME}*registered*nii
  rm -f ${OUTPUTNAME}*Warp*nii ${OUTPUTNAME}*warp*nii
  echo " update $i "
count=100

if [ $DOQSUB -eq 0  ]
then
for x in `ls -tr  $IMAGESETVARIABLE  `
do					# Increment the counter
   echo " count is $count  and file is $x at it $i "
BASE=${x%.*}
    NAMING=${OUTPUTNAME}${BASE}
echo " BASE $NAMING "

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

  if [  $i -eq 1 ]
  then
    for (( n = 0 ; n < ${NUMLEVELS}; n++ ))
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

  if [  $i -gt 1 ]
  then
    ITERATIONS=$MAXITERATIONS
  fi
    echo ITERATIONS $ITERATIONS
  LOCALMETRIC=${METRIC},$x,${METRICPARAMS}

     ${ANTSPATH}ANTS $DIM -m  $LOCALMETRIC  -t $TRANSFORMATION  -r $REGULARIZATION -o ${NAMING}   -i $ITERATIONS
#--MI-option 16x8000 #-a InitAffine.txt --continue-affine 0

      ${ANTSPATH}WarpImageMultiTransform $DIM  $x    ${NAMING}registered.nii ${NAMING}Warp.nii ${NAMING}Affine.txt  -R ${TEMPLATE}

#    ${ANTSPATH}ComposeMultiTransform $DIM   ${NAMING}Warp.nii  -R ${TEMPLATE} ${NAMING}Warp.nii ${NAMING}Affine.txt   -R ${TEMPLATE}

if [ $DIM -lt 3  ]
then
      ${ANTSPATH}ConvertToJpg  ${NAMING}registered.nii ${NAMING}registered.jpg
      ${ANTSPATH}MeasureImageSimilarity $DIM 2 ${TEMPLATE} ${NAMING}registered.nii ${TEMPLATENAME}metriclog.txt
fi

    count=`expr $count + 1`

done
fi


if [ $DOQSUB -eq 1  ]
then
echo " perl registerimages.pl \" $IMAGESETVARIABLE \"   ${TEMPLATE}  $DIM "
MAPTHESE=$IMAGESETVARIABLE
echo " NOT IMPLEMENTED "
#perl registerimages.pl  \"$IMAGESETVARIABLE\" ${TEMPLATE} $DIM
exit
#perl warpimages.pl  ${IMAGESETVARIABLE}  ${TEMPLATE}   $DIM
fi

echo " finished $i " >> ${TEMPLATENAME}metriclog.txt

    ${ANTSPATH}AverageImages $DIM ${TEMPLATE} 1 ${OUTPUTNAME}*registered.nii

LESSLIMIT=$((${ITERATIONLIMIT} - 1))
if [  $i  -lt  ${LESSLIMIT} ]
 then
   rm -f  ${OUTPUTNAME}*InverseWarp*vec.nii

# below, a cheap approach to integrating the negative velocity field
# in the absence of other code and saving the velocity fields.
# additionally, this works for all types of registration algorithms

     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpxvec.nii 0 ${OUTPUTNAME}*Warpxvec.nii
     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpyvec.nii 0 ${OUTPUTNAME}*Warpyvec.nii
if [ $DIM -gt 2  ]
then
     ${ANTSPATH}AverageImages $DIM ${TEMPLATENAME}warpzvec.nii 0 ${OUTPUTNAME}*Warpzvec.nii
fi
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpxvec.nii -0.15 ${TEMPLATENAME}warpxvec.nii
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpyvec.nii -0.15  ${TEMPLATENAME}warpyvec.nii
if [ $DIM -gt 2  ]
then
     ${ANTSPATH}MultiplyImages  $DIM ${TEMPLATENAME}warpzvec.nii -0.15  ${TEMPLATENAME}warpzvec.nii
fi
    ${ANTSPATH}WarpImageMultiTransform $DIM  ${TEMPLATE}   ${TEMPLATE} ${TEMPLATENAME}warp.nii ${TEMPLATENAME}warp.nii ${TEMPLATENAME}warp.nii  ${TEMPLATENAME}warp.nii  -R ${TEMPLATE}
if [ $DIM -lt 3  ]
then
    ${ANTSPATH}ConvertToJpg ${TEMPLATE} ${TEMPLATENAME}$i.jpg
fi
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpxvec.nii  ${TEMPLATENAME}warpxlog.txt  1
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpyvec.nii  ${TEMPLATENAME}warpylog.txt  1
if [ $DIM -gt 2  ]
then
    ${ANTSPATH}MeasureMinMaxMean $DIM ${TEMPLATENAME}warpzvec.nii  ${TEMPLATENAME}warpylog.txt  1
fi
fi

i=$((i + 1))

done


# end loop


cat ${TEMPLATENAME}metriclog.txt
cat ${TEMPLATENAME}warp*log.txt
