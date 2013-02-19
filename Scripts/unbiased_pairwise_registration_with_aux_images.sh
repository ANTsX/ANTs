
#!/bin/bash
usage=" $0 -d 3 -f fixed.nii.gz -m moving.nii.gz -t grouptemplate.nii.gz -b templatebrainmask.nii.gz -g f_auximage_tensor -n m_auximage_tensor -o output_prefix -h f_auximage_scalar -k m_auximage_scalar "
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

i don't worry about the header bias too much or the fact that my first transform depends on mapping B to A --- the interpolation is still symmetric ( i think ) with this approach.

would be nice to use the CompositeTransformUtil to convert the output of this to just a fwd/inv tx.

if this turns out to be biased, i suppose we need header tricks.

supercalifragilisticexpialidocious

A=A ; B=B ; prefix=J ; dim=3
if [[ $# -eq 0 ]] ; then echo $usage ; exit 0 ; fi
while getopts ":d:f:m:b:o:t:g:h:n:k:h:" opt; do
  case $opt in
    d)
      echo "-d $OPTARG" >&2
      dim=$OPTARG
      ;;
    f)
      echo "-f $OPTARG" >&2
      A=$OPTARG
      ;;
    g)
      echo "-g $OPTARG" >&2
      G=$OPTARG
      ;;
    h)
      echo "-h $OPTARG" >&2
      H=$OPTARG
      ;;
    k)
      echo "-k $OPTARG" >&2
      K=$OPTARG
      ;;
    m)
      echo "-m $OPTARG" >&2
      B=$OPTARG
      ;;
    n)
      echo "-n $OPTARG" >&2
      N=$OPTARG
      ;;
    o)
      echo "-o $OPTARG " >&2
      prefix=$OPTARG
      ;;
    t)
      echo "-t $OPTARG " >&2
      template=$OPTARG
      ;;
    b)
      echo "-b $OPTARG " >&2
      templatebm=$OPTARG
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
affits=999x550x20
aff=" -t rigid[ 0.2 ]  -c [ $affits ,1.e-7,20]  -s 3x2x0 -f 4x2x1 -u $uval -l 0 "
synits=100x50  #BA 
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
nmA=${prefix}_A_norm
nmB=${prefix}_B_norm
nm=${prefix}_B_to_A
nminv=${prefix}_A_to_B

################
#  T1 Mapping  #
################
initAmat=${prefix}_A_norm0GenericAffine.mat
initBmat=${prefix}_B_norm0GenericAffine.mat
FWD=" -t [ $initAmat , 1 ] -t ${nm}1Warp.nii.gz -t  $initB -r $A "
INV=" -t [ $initBmat , 1 ] -t ${nm}1InverseWarp.nii.gz -t  $initA -r $B  "
echo  FWD   $FWD
echo  INV   $INV
initA=${prefix}_initA
initB=${prefix}_initB
#####################
if [[ ! -s ${nm}1Warp.nii.gz ]]  ; then
$reg -d $dim -r [ $A, $B, 1 ]   \
                        -m mattes[  $A, $B , 1 , 32, random , 0.25 ] $aff -z 1 \
                       -o [${initA}]
$reg -d $dim -r [ $B, $A, 1 ] \
                        -m mattes[  $B, $A , 1 , 32, random , 0.25 ] $aff -z 1 \
                       -o [${initB}]
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
$reg -d $dim  \
                        -m mattes[  ${prefix}_mid.nii.gz, $A, 1 , 32, random , 0.25 ] $aff \
                       -o [${nmA},${nmA}_aff.nii.gz]
# compute the map from B to midpoint(B,A) --- "fair" interpolation
$reg -d $dim  \
                        -m mattes[  ${nmA}_aff.nii.gz, $B, 1 , 32, random , 0.25 ] $aff \
                       -o [${nmB},${nmB}_aff.nii.gz]
# now we can do a symmetric deformable mapping
antsApplyTransforms -d $dim -i $B -o ${nm}_aff.nii.gz -t [ $initAmat, 1 ] -t  $initBmat -r $A
N3BiasFieldCorrection $dim $A ${nm}_n3_a.nii.gz 4 
N3BiasFieldCorrection $dim ${nm}_n3_a.nii.gz ${nm}_n3_a.nii.gz 2 
N3BiasFieldCorrection $dim $B ${nm}_n3_b.nii.gz 4 
N3BiasFieldCorrection $dim ${nm}_n3_b.nii.gz ${nm}_n3_b.nii.gz 2 
  echo now do deformable expecting $initB and $initA to exist
  if [[ -s $initAmat  ]] && [[ -s $initBmat ]] ; then
    $reg -d $dim  --initial-fixed-transform $initAmat  --initial-moving-transform $initBmat \
                         -m mattes[  ${nm}_n3_a.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [${synits},1.e-8,10]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [${nm},${nm}_diff_symm.nii.gz]
#    $reg -d $dim  --initial-fixed-transform $initBmat  --initial-moving-transform $initAmat \
#                         -m cc[  $B, $A , 1 , 4 ] \
#                         -t syn[ 0.25, 3, 0.5 ] \
#                         -c [${synits},1.e-8,10]  \
#                        -s 1x0 \
#                        -f 2x1 \
#                       -u $uval \
#                       -o [${nminv},${nminv}_diff_symm.nii.gz]
  else 
    echo  $initBmat and $initAmat DO NOT exist 
    exit
  fi
fi
antsApplyTransforms -d $dim -i $B -o ${nm}_diff.nii.gz -t [ $initAmat, 1 ] -t ${nm}1Warp.nii.gz -t  $initBmat -r $A
ANTSJacobian $dim ${nm}1Warp.nii.gz ${nm} 1 no 0
# antsApplyTransforms -d $dim -i ${nm}logjacobian.nii.gz -o  ${nm}logjacobian.nii.gz -t [ $initAmat, 1 ] -r $A
# CompositeTransformUtil --assemble ${prefix}_fwd.mat $FWD
# CompositeTransformUtil --assemble ${prefix}_inv.mat $INV
# ImageMath $dim ${prefix}invid.nii.gz InvId ${nminv}1Warp.nii.gz ${nm}1Warp.nii.gz


initafffn=${nm}_init_aff.mat 
if [[ -s $template ]] && [[ ! -s  ${nm}_gt_0GenericAffine.mat ]] ; then 
  imgsmall=${nm}_diffsmall.nii.gz
  temsmall=${nm}_temsmall.nii.gz
  ResampleImageBySpacing $dim ${nm}_diff_symm.nii.gz $imgsmall 4 4 4 
  ResampleImageBySpacing $dim $template $temsmall 4 4 4 
  antsAffineInitializer 3 $temsmall  $imgsmall  $initafffn 15 0.1 0 10 
  gf=$template
  gm=${nm}_diff_symm.nii.gz 
  imgs=" $gf, $gm "
  $reg -d $dim  -r $initafffn  \
                        -m mattes[  $imgs , 1 , 32, regular, 0.25 ] \
                         -t affine[ 0.1 ] \
                         -c [ $affits ,1.e-7,20]  \
                        -s 4x2x1vox  \
                        -f 4x2x1 -l 1 \
                        -m cc[  $imgs , 1 , 4 ] \
                         -t syn[ .2, 3, 0.0 ] \
                         -c [100x50x20,1.e-8,20]  \
                        -s 2x1x0vox  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [${nm}_gt_,${nm}_gt.nii.gz]

# map brain mask to subject space T1
trans=" -t [ $initAmat, 1 ] -t [ ${nm}_gt_0GenericAffine.mat, 1] -t  ${nm}_gt_1InverseWarp.nii.gz "
antsApplyTransforms -d $dim -i $templatebm -o  ${nm}_bm_A.nii.gz -n NearestNeighbor -r $A $trans
trans=" -t [ $initBmat, 1 ] -t ${nm}1InverseWarp.nii.gz -t [ ${nm}_gt_0GenericAffine.mat, 1] -t  ${nm}_gt_1InverseWarp.nii.gz "
antsApplyTransforms -d $dim -i $templatebm -o  ${nm}_bm_B.nii.gz -n NearestNeighbor -r $B $trans
MultiplyImages $dim  ${nm}_bm_A.nii.gz $A  ${nm}_A_brain.nii.gz
MultiplyImages $dim  ${nm}_bm_B.nii.gz $B  ${nm}_B_brain.nii.gz


  ######### now redo bias correction & syn #########
  N3BiasFieldCorrection $dim ${nm}_A_brain.nii.gz ${nm}_n3_a.nii.gz 4 
  N3BiasFieldCorrection $dim ${nm}_n3_a.nii.gz ${nm}_n3_a.nii.gz 2  
  N3BiasFieldCorrection $dim ${nm}_n3_a.nii.gz ${nm}_n3_a.nii.gz 2 
  N3BiasFieldCorrection $dim ${nm}_B_brain.nii.gz ${nm}_n3_b.nii.gz 4 
  N3BiasFieldCorrection $dim ${nm}_n3_b.nii.gz ${nm}_n3_b.nii.gz 2 
  N3BiasFieldCorrection $dim ${nm}_n3_b.nii.gz ${nm}_n3_b.nii.gz 2 
  echo now do deformable expecting $initB and $initA to exist
  if [[ -s $initAmat  ]] && [[ -s $initBmat ]] ; then
    $reg -d $dim  --initial-fixed-transform $initAmat  --initial-moving-transform $initBmat \
                         -m mattes[  ${nm}_n3_a.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [${synits},1.e-8,10]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [${nm},${nm}_diff_symm.nii.gz]
  else 
    echo  $initBmat and $initAmat DO NOT exist 
    exit
  fi
  antsApplyTransforms -d $dim -i $B -o ${nm}_diff.nii.gz -t [ $initAmat, 1 ] -t ${nm}1Warp.nii.gz -t  $initBmat -r $A
  ANTSJacobian $dim ${nm}1Warp.nii.gz ${nm} 1 no 0
fi 
if [[ -s $G ]] && [[ -s $N ]] && [[ ! -s ${nm}_fadiff.nii.gz  ]] ; then 
  echo deal with auxiliary images ... here DTI
  ffa=${nm}_ffa.nii.gz
  mfa=${nm}_mfa.nii.gz
  ImageMath 3 $ffa TensorFA $G
  ImageMath 3 $mfa TensorFA $N
  synits=15x25x3
  imgs=" ${nm}_A_brain.nii.gz, $ffa "
  $reg -d $dim  \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [1000x1000x5,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [1000x25,1.e-7,20]  \
                        -s 2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [$synits,1.e-7,20]  \
                        -s 4x2x1mm -x [${nm}_bm_A.nii.gz] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [${nm}_ffa,${nm}_ffa_distcorr.nii.gz]
  imgs=" ${nm}_B_brain.nii.gz, $mfa "
  $reg -d $dim  \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [1000x1000x5,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [1000x25,1.e-7,20]  \
                        -s 2x1mm  -x [${nm}_bm_B.nii.gz]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [$synits,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_B.nii.gz] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -o [${nm}_mfa,${nm}_mfa_distcorr.nii.gz]
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat                       -t  $initAmat -t ${nm}_ffa1Warp.nii.gz -t ${nm}_ffa0GenericAffine.mat "
antsApplyTransforms -d $dim -i $ffa -o  ${nm}_ffanorm.nii.gz  -r $template $trans
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}1Warp.nii.gz  -t  $initBmat -t ${nm}_mfa1Warp.nii.gz -t ${nm}_mfa0GenericAffine.mat "
antsApplyTransforms -d $dim -i $mfa -o  ${nm}_mfanorm.nii.gz  -r $template $trans
ImageMath $dim ${nm}_fadiff.nii.gz - ${nm}_ffanorm.nii.gz  ${nm}_mfanorm.nii.gz
fi 
echo done with aux images --- now get final jacobians 



if [[ -s $H ]] && [[ -s $K ]] && [[ ! -s ${nm}_cbfdiff.nii.gz  ]] ; then 
  echo deal with auxiliary images ... here DTI
  fcbf=${nm}_fcbf.nii.gz
  mcbf=${nm}_mcbf.nii.gz
  ImageMath 3 $fcbf Normalize $H
  ImageMath 3 $mcbf Normalize $K
  synits=15x25x3
  imgs=" ${nm}_A_brain.nii.gz, $fcbf "
  $reg -d $dim  \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [1000x1000x5,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [1000x25,1.e-7,20]  \
                        -s 2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [$synits,1.e-7,20]  \
                        -s 4x2x1mm -x [${nm}_bm_A.nii.gz] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [${nm}_fcbf,${nm}_fcbf_distcorr.nii.gz]
  imgs=" ${nm}_B_brain.nii.gz, $mcbf "
  $reg -d $dim  \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [1000x1000x5,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_A.nii.gz]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [1000x25,1.e-7,20]  \
                        -s 2x1mm  -x [${nm}_bm_B.nii.gz]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [$synits,1.e-7,20]  \
                        -s 4x2x1mm  -x [${nm}_bm_B.nii.gz] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -o [${nm}_mcbf,${nm}_mcbf_distcorr.nii.gz]
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat                       -t  $initAmat -t ${nm}_fcbf1Warp.nii.gz -t ${nm}_fcbf0GenericAffine.mat "
antsApplyTransforms -d $dim -i $fcbf -o  ${nm}_fcbfnorm.nii.gz  -r $template $trans
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}1Warp.nii.gz  -t  $initBmat -t ${nm}_mcbf1Warp.nii.gz -t ${nm}_mcbf0GenericAffine.mat "
antsApplyTransforms -d $dim -i $mcbf -o  ${nm}_mcbfnorm.nii.gz  -r $template $trans
ImageMath $dim ${nm}_cbfdiff.nii.gz - ${nm}_fcbfnorm.nii.gz  ${nm}_mcbfnorm.nii.gz
fi 
echo done with 2nd aux images --- now get final jacobians 


# get final jacobian values 
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat                       -t  $initAmat"
antsApplyTransforms -d $dim -i ${nm}_A_brain.nii.gz -o  [${nm}_A_fullWarp.nii.gz, 1 ]  -r $template $trans
ANTSJacobian $dim ${nm}_A_fullWarp.nii.gz ${nm}_A_full 1 no 0
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}1Warp.nii.gz  -t  $initBmat"
antsApplyTransforms -d $dim -i ${nm}_B_brain.nii.gz -o  [${nm}_B_fullWarp.nii.gz, 1 ]  -r $template $trans
ANTSJacobian $dim ${nm}_B_fullWarp.nii.gz ${nm}_B_full 1 no 0
