#!/bin/bash
usage=" $0 -d 3 -f fixed.nii.gz -m moving.nii.gz  -o output_prefix "
:<<supercalifragilisticexpialidocious

here is a first pass at an unbiased registration between an image pair, A ( fixed ) and B.

the reference frame for the registration uses the header of A.

however, the image content ---- call it AB ----- is actually half-way between A and B.

so we have ( in geometry )

A . . . . AB . . . . B

or very close to it.

we then compute :

AB =>  A  which gives   Aff_A
AB =>  B  which gives   Aff_B

and , finally      A( Aff_A )   <=>   B( Aff_ B )

where    <=>    is SyN.

so briefly, this implements:

A => midAffA => midWarpAtoMID * MID * midWarpBtoMID <= midAffB <= B

i don't worry about the header bias too much or the fact that my first transform depends on mapping B to A --- the interpolation is still symmetric ( i think ) with this approach.

would be nice to use the CompositeTransformUtil to convert the output of this to just a fwd/inv tx.

if this turns out to be biased, i suppose we need header tricks.

supercalifragilisticexpialidocious

A=A ; B=B ; prefix=J ; dim=3
if [[ $# -eq 0 ]] ; then echo $usage ; exit 0 ; fi
while getopts ":d:f:m:o:t:h:" opt; do
  case $opt in
    d)
      echo "-d $OPTARG" >&2
      dim=$OPTARG
      ;;
    f)
      echo "-f $OPTARG" >&2
      A=$OPTARG
      ;;
    m)
      echo "-m $OPTARG" >&2
      B=$OPTARG
      ;;
    o)
      echo "-o $OPTARG " >&2
      prefix=$OPTARG
      ;;
    h)
      echo "Usage: $usage " >&2
      exit 0
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
reg=antsRegistration
uval=0
aff=" -t affine[ 0.25 ]  -c [ 1009x200x20,1.e-8,20 ]  -s 4x2x0 -f 4x2x1 "
synits=20x20x10
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
nmA=${prefix}_A_norm
nmB=${prefix}_B_norm
nm=${prefix}_B_to_A
nminv=${prefix}_A_to_B
initAmat=${prefix}_initA0GenericAffine.mat
initBmat=${prefix}_initB0GenericAffine.mat
# register in both directions, then average the result
initA=${prefix}_initA
$reg -d $dim -r [ $A, $B, 1 ]   \
                        -m mattes[  $A, $B , 1 , 32, regular , 0.25 ] $aff -z 1 \
                       -o [ ${initA}]
initB=${prefix}_initB
$reg -d $dim -r [ $B, $A, 1 ] \
                        -m mattes[  $B, $A , 1 , 32, regular , 0.25 ] $aff -z 1 \
                       -o [ ${initB}]
# get the identity map
ComposeMultiTransform $dim ${initA}_id.mat -R ${initA}0GenericAffine.mat  ${initA}0GenericAffine.mat -i ${initA}0GenericAffine.mat
# invert the 2nd affine registration map
ComposeMultiTransform $dim ${initB}_inv.mat -R ${initA}0GenericAffine.mat -i ${initB}0GenericAffine.mat
# get the average affine map
AverageAffineTransform $dim ${prefix}_avg.mat  ${initB}_inv.mat ${initA}0GenericAffine.mat
# get the midpoint affine map
AverageAffineTransform $dim ${prefix}_mid.mat   ${initA}_id.mat  ${prefix}_avg.mat
#.........#
# this applies, to B, a map from B to midpoint(B,A)
antsApplyTransforms -d $dim -i $B -o ${prefix}_mid.nii.gz -t  ${prefix}_mid.mat  -r  $A
# compute the map from A to midpoint(B,A) --- "fair" interpolation
$reg -d $dim -r  [ ${prefix}_mid.nii.gz, $A, 1 ] \
                        -m mattes[  ${prefix}_mid.nii.gz, $A, 1 , 32, random , 0.25 ] $aff \
                       -o [ ${nmA},${nmA}_aff.nii.gz]
# compute the map from B to midpoint(B,A) --- "fair" interpolation
$reg -d $dim  -r [ ${nmA}_aff.nii.gz, $B, 1 ] \
                        -m mattes[  ${nmA}_aff.nii.gz, $B, 1 , 32, random , 0.25 ] $aff \
                       -o [ ${nmB},${nmB}_aff.nii.gz]
midaaff=${nmA}0GenericAffine.mat
midbaff=${nmB}0GenericAffine.mat
# now we can do a symmetric deformable mapping - not sure
echo now do deformable expecting $initB and $initA to exist
$reg -d $dim  --initial-fixed-transform $midaaff  --initial-moving-transform $midbaff \
                         -m mattes[  $A, $B , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.0 ] \
                         -c [ ${synits},1.e-8,10 ]  \
                        -s 2x1x0 \
                        -f 4x2x1 \
                       -o [ ${nm},${nm}_diff_symm.nii.gz]
# this is the composite mapping - via the mid-point - to the A image
antsApplyTransforms -d $dim -i $B -o ${nm}_diffX.nii.gz    -t [ $midaaff, 1 ] -t ${nm}1Warp.nii.gz        -t  $midbaff -r $A
# one can gain the inverse in a similar way
antsApplyTransforms -d $dim -i $A -o ${nminv}_diffX.nii.gz -t [ $midbaff, 1 ] -t ${nm}1InverseWarp.nii.gz -t  $midaaff -r $B
