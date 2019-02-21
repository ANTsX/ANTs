#!/bin/bash
#
dim=3 # image dimensionality
AP="" # /home/yourself/code/ANTS/bin/bin/  # path to ANTs binaries
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2    # fixed and moving image file names
mysetting=$3
if [[ ! -s $f ]] ; then echo no fixed $f ; exit; fi
if [[ ! -s $m ]] ; then echo no moving $m ;exit; fi
if [[ ${#mysetting} -eq 0 ]] ; then
echo usage is
echo $0 fixed.nii.gz moving.nii.gz mysetting
echo  where mysetting is either forproduction or fastfortesting
exit
fi
nm1=` basename $f | cut -d '.' -f 1 `
nm2=` basename $m | cut -d '.' -f 1 `
reg=${AP}antsRegistration           # path to antsRegistration
if [[ $mysetting == "fastfortesting" ]] ; then
  its=10000x0x0
  percentage=0.1
  syn="100x0x0,0,5"
else
  its=10000x111110x11110
  percentage=0.3
  syn="100x100x50,-0.01,5"
  mysetting=forproduction
fi
echo affine $m $f outname is $nm am using setting $mysetting
nm=${D}${nm1}_fixed_${nm2}_moving_setting_is_${mysetting}   # construct output prefix
$reg -d $dim -r [ $f, $m ,1 ]  \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t translation[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 6x4x2 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t rigid[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t affine[ 0.1 ] \
                         -c [ $its,1.e-8,20 ]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 0.5 , 32 ] \
                        -m cc[  $f, $m , 0.5 , 4 ] \
                         -t SyN[ .20, 3, 0 ] \
                         -c [ $syn ]  \
                        -s 1x0.5x0vox  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                       -o [ ${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz]

${AP}antsApplyTransforms -d $dim -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -o ${nm}_warped.nii.gz
