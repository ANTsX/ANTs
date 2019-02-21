#!/bin/bash
usage=" $0 -d 3 -f symmetric_template.nii.gz -m moving.nii.gz -o output_prefix "
:<<supercalifragilisticexpialidocious

here is a first pass at registration-based asymmetry based on mapping an image A to a symmetric-template.

note - care must be taken to look at prefix_L_aff.nii.gz to make sure that it's correct - we dont yet have an automated way to verify this - though we know the theory.  if prefix_L_aff.nii.gz is not well aligned, results will be invalid.  

supercalifragilisticexpialidocious

A=A ; B=B ; prefix=J ; dim=3 ; a=1
if [[ $# -eq 0 ]] ; then echo $usage ; exit 0 ; fi
while getopts ":d:f:m:o:h:a:" opt; do
  case $opt in
    d)
      echo "-d $OPTARG" >&2
      dim=$OPTARG
      ;;
    f)
      echo "-f $OPTARG" >&2
      A=$OPTARG
      ;;
    h)
      echo $usage
      exit 0;
      ;;
    m)
      echo "-m $OPTARG" >&2
      B=$OPTARG
      ;;
    o)
      echo "-o $OPTARG " >&2
      prefix=$OPTARG
      ;;
    a)
      echo "-a $OPTARG " >&2
      a=$OPTARG
      ;;
    \?)
      echo "Usage: $usage " >&2
      exit 0
      ;;
  esac
done
echo inputs: $A $B $prefix $dim
if [[ ${#dim} -lt 1 ]] ; then echo must provide input dimension $dim ; echo $usage ; exit 0 ; fi
if [[ ${#prefix} -lt 3 ]] ; then echo must provide output prefix $prefix ; echo $usage ; exit 0 ; fi
if [[ ! -s $A ]] || [[  ! -s $B ]]  ; then echo inputs: $A $B $prefix ; echo $usage ; exit 1 ; fi
#####################################################
reg=antsRegistration
uval=0
affits=999x999x1550x200
rig=" -t rigid[ 0.2 ]  -c [ $affits ,1.e-7,20 ]  -s 3x2x1x0 -f 8x4x2x1 -u $uval -l 0 "
aff=" -t affine[ 0.2 ]  -c [ $affits ,1.e-7,20 ]  -s 3x2x1x0 -f 8x4x2x1 -u $uval -l 0 "
metparams=" 1 , 32, regular , 0.5 "
synits=220x220x100x50  #BA 
# synits=0x0x0x0  #BA 
dtx="syn[ 0.25, 3, 0. ]"
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
#####################################################
ImageMath $dim ${prefix}_reflection.mat ReflectionMatrix $A $a
antsApplyTransforms -d $dim -t ${prefix}_reflection.mat -i $B -o ${prefix}_reflect.nii.gz -r $B
imgs=" $A, $B "
antsAffineInitializer ${dim} $A $B ${prefix}_init.mat 5 0.25 0 3
$reg -d $dim -r ${prefix}_init.mat\
    -m mattes[ $imgs, $metparams ] $rig \
    -m mattes[ $imgs, $metparams ] $aff \
    -u $uval -b 0 -z 1 \
    -o [ ${prefix}_L,${prefix}_L_aff.nii.gz]
$reg -d $dim -r ${prefix}_L0GenericAffine.mat \
    -m mattes[ $imgs , 1 , 32 ] \
    -t $dtx \
    -c [ ${synits},1.e-8,10 ]  \
    -s 3x2x1x0 \
    -f 8x4x2x1 \
    -u $uval -b 0 -z 1 \
    -o [ ${prefix}_L,${prefix}_L.nii.gz]

$reg -d $dim   -r  ${prefix}_reflection.mat -r ${prefix}_L0GenericAffine.mat \
    -m mattes[ $imgs , 1 , 32 ] \
    -t $dtx \
    -c [ ${synits},1.e-8,10 ]  \
    -s 3x2x1x0 \
    -f 8x4x2x1 \
    -u $uval -b 0 -z 1 \
    -o [ ${prefix}_R,${prefix}_R.nii.gz]

##########################################################
ANTSJacobian $dim ${prefix}_R1Warp.nii.gz  ${prefix}_R 1
ANTSJacobian $dim ${prefix}_L1Warp.nii.gz  ${prefix}_L 1
ImageMath $dim ${prefix}_asym.nii.gz  - ${prefix}_Llogjacobian.nii.gz ${prefix}_Rlogjacobian.nii.gz
