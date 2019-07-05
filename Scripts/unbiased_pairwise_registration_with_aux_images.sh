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
echo THIS SCRIPT MAY OR MAY NOT WORK - IT HAS NOT BEEN TESTED RECENTLY 
echo USE AT OWN RISK BY COMMENTING OUT THE EXIT CALL BELOW
exit 0
if [[ ${#dim} -lt 1 ]] ; then echo must provide input dimension $dim ; echo $usage ; exit 0 ; fi
if [[ ${#prefix} -lt 3 ]] ; then echo must provide output prefix $prefix ; echo $usage ; exit 0 ; fi
if [[ ! -s $A ]] || [[  ! -s $B ]]  ; then echo inputs: $A $B $prefix ; echo $usage ; exit 1 ; fi

reg=antsRegistration
uval=0
affits=999x550x20
aff=" -t affine[ 0.2 ]  -c [ $affits ,1.e-7,20 ]  -s 3x2x0 -f 4x2x1 -u $uval -l 0 "
metparams=" 1 , 32, random , 0.25 "
synits=100x50  #BA 
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
nmA=${prefix}_A_norm
nmB=${prefix}_B_norm
nm=${prefix}_Long
nminv=${prefix}_ILong
inA=${prefix}_A.nii.gz
inB=${prefix}_B.nii.gz 
if [[ ! -s $inA ]] && [[ ! -s $inB ]] ; then 
  cp $A $inA
  cp $B $inB 
  A=$inA
  B=$inB
  # map the input images to a similar intensity space
  ImageMath $dim $A Normalize $A 
  ImageMath $dim $B Normalize $B
  if [[ -s $template ]] ; then 
      ImageMath $dim $A HistogramMatch $A $template
      ImageMath $dim $B HistogramMatch $B $template
  fi 
  N3BiasFieldCorrection $dim $A $A 4
  N3BiasFieldCorrection $dim $B $B 4
fi 
A=$inA
B=$inB
################
#  T1 Mapping  #
################
initAmat=${prefix}_A_norm0GenericAffine.mat
initBmat=${prefix}_B_norm0GenericAffine.mat
initA=${prefix}_initA
initB=${prefix}_initB
#####################
if [[ ! -s ${nm}1Warp.nii.gz ]]  ; then
$reg -d $dim -r [ $A, $B, 1 ]   \
                        -m mattes[  $A, $B , $metparams ] $aff -z 1 \
                       -o [ ${initA}]
$reg -d $dim -r [ $B, $A, 1 ] \
                        -m mattes[  $B, $A , $metparams ] $aff -z 1 \
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
ImageMath $dim  ${prefix}_mid.nii.gz PadImage $A 10 
antsApplyTransforms -d $dim -i $B -o ${prefix}_mid.nii.gz -t  ${prefix}_mid.mat  -r  $A
# compute the map from A to midpoint(B,A) --- "fair" interpolation
$reg -d $dim  \
                        -m mattes[  ${prefix}_mid.nii.gz, $A, $metparams ] $aff \
                       -o [ ${nmA},${nmA}_aff.nii.gz]
# compute the map from B to midpoint(B,A) --- "fair" interpolation
$reg -d $dim  \
                        -m mattes[  ${nmA}_aff.nii.gz, $B, $metparams ] $aff \
                       -o [ ${nmB},${nmB}_aff.nii.gz]

# now we can do a symmetric deformable mapping
N3BiasFieldCorrection $dim $A ${nm}_n3_a.nii.gz 4 
N3BiasFieldCorrection $dim ${nm}_n3_a.nii.gz ${nm}_n3_a.nii.gz 2 
N3BiasFieldCorrection $dim $B ${nm}_n3_b.nii.gz 4 
N3BiasFieldCorrection $dim ${nm}_n3_b.nii.gz ${nm}_n3_b.nii.gz 2 
  echo now do deformable expecting $initB and $initA to exist
  if [[ -s $initAmat ]] && [[ -s $initBmat ]] ; then
    $reg -d $dim  --initial-fixed-transform $initAmat  --initial-moving-transform $initBmat \
                         -m mattes[  ${nm}_n3_a.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [ ${synits},1.e-8,10 ]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [ ${nm},${nm}_diff_symm.nii.gz]
#    $reg -d $dim  --initial-fixed-transform $initBmat  --initial-moving-transform $initAmat \
#                         -m mattes[  ${nm}_n3_a.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
#                         -t syn[ 0.25, 3, 0.5 ] \
#                         -c [ ${synits},1.e-8,10 ]  \
#                        -s 1x0 \
#                        -f 2x1 \
#                       -u $uval -b 0 -z 1 \
#                       -o [ ${nminv},${nminv}_diff_symm.nii.gz]
  else
    echo  $initBmat and $initAmat DO NOT exist
    exit
  fi
fi
AverageImages  $dim ${nm}_avg.nii.gz 0 ${nmA}_aff.nii.gz ${nm}_diff_symm.nii.gz
MultiplyImages $dim ${nm}1InverseWarp.nii.gz 0.5 ${nm}tempWarp.nii.gz 
antsApplyTransforms -d $dim -i ${nm}_avg.nii.gz -o ${nm}_avg.nii.gz -t ${nm}tempWarp.nii.gz -r ${nm}_avg.nii.gz
rm ${nm}tempWarp.nii.gz 

# recompute the mappings 
  if [[ -s $initAmat ]] && [[ -s $initBmat ]] ; then
    $reg -d $dim  --initial-moving-transform $initBmat \
                         -m mattes[  ${nm}_avg.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [ ${synits},1.e-8,10 ]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [ ${nm}_B,${nm}_B_symm.nii.gz]
    $reg -d $dim  --initial-moving-transform $initAmat \
                         -m mattes[  ${nm}_avg.nii.gz, ${nm}_n3_a.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [ ${synits},1.e-8,10 ]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [ ${nm}_A,${nm}_A_symm.nii.gz]
  else 
    echo  $initBmat and $initAmat DO NOT exist 
    exit
  fi

initafffn=${nm}_init_aff.mat
if [[ -s $template ]] && [[ ! -s  ${nm}_gt_0GenericAffine.mat ]] ; then
  imgsmall=${nm}_diffsmall.nii.gz
  temsmall=${nm}_temsmall.nii.gz
  ResampleImageBySpacing $dim ${nm}_avg.nii.gz $imgsmall 4 4 4 
  ResampleImageBySpacing $dim $template $temsmall 4 4 4 
  antsAffineInitializer $dim $temsmall  $imgsmall  $initafffn 15 0.1 0 10 
  gf=$template
  gm=${nm}_avg.nii.gz
  imgs=" $gf, $gm "
  $reg -d $dim  -r $initafffn  \
                        -m mattes[  $imgs , 1 , 32, regular, 0.25 ] \
                         -t affine[ 0.1 ] \
                         -c [ $affits ,1.e-7,20 ]  \
                        -s 4x2x1vox  \
                        -f 4x2x1 -l 1 \
                        -m cc[  $imgs , 1 , 4 ] \
                         -t syn[ .2, 3, 0.0 ] \
                         -c [ 100x50x20,1.e-8,20 ]  \
                        -s 2x1x0vox  \
                        -f 4x2x1 -l 1 -u 0 -z 1 \
                       -o [ ${nm}_gt_,${nm}_gt.nii.gz]

  # map brain mask to subject space T1
  trans=" -t [ $initAmat, 1 ] -t ${nm}_A1InverseWarp.nii.gz -t [ ${nm}_gt_0GenericAffine.mat, 1 ] -t  ${nm}_gt_1InverseWarp.nii.gz "
  antsApplyTransforms -d $dim -i $templatebm -o  ${nm}_bm_A.nii.gz -n NearestNeighbor -r $A $trans
  trans=" -t [ $initBmat, 1 ] -t ${nm}_B1InverseWarp.nii.gz -t [ ${nm}_gt_0GenericAffine.mat, 1 ] -t  ${nm}_gt_1InverseWarp.nii.gz "
  antsApplyTransforms -d $dim -i $templatebm -o  ${nm}_bm_B.nii.gz -n NearestNeighbor -r $B $trans
  MultiplyImages $dim  ${nm}_bm_A.nii.gz $A  ${nm}_A_brain.nii.gz
  MultiplyImages $dim  ${nm}_bm_B.nii.gz $B  ${nm}_B_brain.nii.gz
fi 
echo done with brain extraction 

if [[ -s $G ]] && [[ -s $N ]] && [[ ! -s ${nm}_fadiff.nii.gz ]] ; then 
# map brain mask to subject space T1
trans=" -t [ $initAmat, 1 ] -t [ ${nm}_gt_0GenericAffine.mat, 1 ] -t  ${nm}_gt_1InverseWarp.nii.gz "
antsApplyTransforms -d $dim -i $templatebm -o  ${nm}_bm_A.nii.gz -n NearestNeighbor -r $A $trans
trans=" -t [ $initBmat, 1 ] -t ${nm}1InverseWarp.nii.gz -t [ ${nm}_gt_0GenericAffine.mat, 1 ] -t  ${nm}_gt_1InverseWarp.nii.gz "
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
  if [[ -s $initAmat ]] && [[ -s $initBmat ]] ; then
    $reg -d $dim  --initial-fixed-transform $initAmat  --initial-moving-transform $initBmat \
                         -m mattes[  ${nm}_n3_a.nii.gz, ${nm}_n3_b.nii.gz , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0.5 ] \
                         -c [ ${synits},1.e-8,10 ]  \
                        -s 1x0 \
                        -f 2x1 \
                       -u $uval -b 0 -z 1 \
                       -o [ ${nm},${nm}_diff_symm.nii.gz]
  else
    echo  $initBmat and $initAmat DO NOT exist
    exit
  fi
  antsApplyTransforms -d $dim -i $B -o ${nm}_diff.nii.gz -t [ $initAmat, 1 ] -t ${nm}1Warp.nii.gz -t  $initBmat -r $A
  ANTSJacobian $dim ${nm}1Warp.nii.gz ${nm} 1 no 0
fi
if [[ -s $G ]] && [[ -s $N ]] && [[ ! -s ${nm}_fadiff.nii.gz ]] ; then
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
                         -c [ 1000x1000x5,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [ 1000x25,1.e-7,20 ]  \
                        -s 2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [ $synits,1.e-7,20 ]  \
                        -s 4x2x1mm -x [ ${nm}_bm_A.nii.gz ] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [ ${nm}_ffa,${nm}_ffa_distcorr.nii.gz]
  imgs=" ${nm}_B_brain.nii.gz, $mfa "
  $reg -d $dim  \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [ 1000x1000x5,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [ 1000x25,1.e-7,20 ]  \
                        -s 2x1mm  -x [ ${nm}_bm_B.nii.gz ]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [ $synits,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_B.nii.gz ] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -o [ ${nm}_mfa,${nm}_mfa_distcorr.nii.gz]
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}_A1Warp.nii.gz -t  $initAmat -t ${nm}_ffa1Warp.nii.gz -t ${nm}_ffa0GenericAffine.mat "
antsApplyTransforms -d $dim -i $ffa -o  ${nm}_ffanorm.nii.gz  -r $template $trans
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}_B1Warp.nii.gz  -t  $initBmat -t ${nm}_mfa1Warp.nii.gz -t ${nm}_mfa0GenericAffine.mat "
antsApplyTransforms -d $dim -i $mfa -o  ${nm}_mfanorm.nii.gz  -r $template $trans
ImageMath $dim ${nm}_fadiff.nii.gz - ${nm}_ffanorm.nii.gz  ${nm}_mfanorm.nii.gz
fi
echo done with aux images --- now get final jacobians

 

if [[ -s $H ]] && [[ -s $K ]] && [[ ! -s ${nm}_cbfdiff.nii.gz ]] ; then 
  echo deal with auxiliary images ... here a scalar image pair 
# if [[ -s $H ]] && [[ -s $K ]] && [[ ! -s ${nm}_cbfdiff.nii.gz ]] ; then
#  echo deal with auxiliary images ... here DTI
  fcbf=${nm}_fcbf.nii.gz
  mcbf=${nm}_mcbf.nii.gz
  ImageMath $dim $fcbf Normalize $H
  ImageMath $dim $mcbf Normalize $K
  synits=15x25x3
  imgs=" ${nm}_A_brain.nii.gz, $fcbf "
  $reg -d $dim  \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [ 1000x1000x5,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [ 1000x25,1.e-7,20 ]  \
                        -s 2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [ $synits,1.e-7,20 ]  \
                        -s 4x2x1mm -x [ ${nm}_bm_A.nii.gz ] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [ ${nm}_fcbf,${nm}_fcbf_distcorr.nii.gz]
  imgs=" ${nm}_B_brain.nii.gz, $mcbf "
  $reg -d $dim  \
                         -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t rigid[ 0.1 ] \
                         -c [ 1000x1000x5,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_A.nii.gz ]  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -m mattes[  $imgs , 1 , 32, regular, 0.2 ] \
                         -t affine[ 0.1 ] \
                         -c [ 1000x25,1.e-7,20 ]  \
                        -s 2x1mm  -x [ ${nm}_bm_B.nii.gz ]  \
                        -f 2x1 -l 1 -u 1 -z 1 \
                         -m mattes[  $imgs , 1 , 32 ] \
                         -m cc[  $imgs , 1 , 2 ] \
                         -t SyN[ 0.2, 3, 0.5 ] \
                         -c [ $synits,1.e-7,20 ]  \
                        -s 4x2x1mm  -x [ ${nm}_bm_B.nii.gz ] \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -o [ ${nm}_mcbf,${nm}_mcbf_distcorr.nii.gz]
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}_A1Warp.nii.gz  -t  $initAmat -t ${nm}_fcbf1Warp.nii.gz -t ${nm}_fcbf0GenericAffine.mat "
exe1="antsApplyTransforms -d $dim -i $fcbf -o  ${nm}_fcbfnorm.nii.gz  -r $template $trans"
$exe1
#####
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}_B1Warp.nii.gz  -t  $initBmat -t ${nm}_mcbf1Warp.nii.gz -t ${nm}_mcbf0GenericAffine.mat "
exe2="antsApplyTransforms -d $dim -i $mcbf -o  ${nm}_mcbfnorm.nii.gz  -r $template $trans"
$exe2
echo $exe1 > ${nm}_cbf_map.txt 
echo $exe2 >> ${nm}_cbf_map.txt 
#####
ImageMath $dim ${nm}_cbfdiff.nii.gz - ${nm}_fcbfnorm.nii.gz  ${nm}_mcbfnorm.nii.gz
fi
echo done with 2nd aux images --- now get final jacobians

# get final jacobian values 
# trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat  -t ${nm}_A1Warp.nii.gz  -t  $initAmat"
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat                       -t  $initAmat"
antsApplyTransforms -d $dim -i ${nm}_A_brain.nii.gz -o  [ ${nm}_A_fullWarp.nii.gz, 1 ]  -r $template $trans
ANTSJacobian $dim ${nm}_A_fullWarp.nii.gz ${nm}_A_full 1 no 0
#
trans=" -t ${nm}_gt_1Warp.nii.gz -t ${nm}_gt_0GenericAffine.mat -t ${nm}_B1Warp.nii.gz  -t  $initBmat"
antsApplyTransforms -d $dim -i ${nm}_B_brain.nii.gz -o  [ ${nm}_B_fullWarp.nii.gz, 1 ]  -r $template $trans
ANTSJacobian $dim ${nm}_B_fullWarp.nii.gz ${nm}_B_full 1 no 0
