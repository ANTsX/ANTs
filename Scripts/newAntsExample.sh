#!/bin/bash
#
dim=3  # image dimensionality
AP="" # /home/yourself/code/ANTS/bin/bin/  # path to ANTs binaries
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2    # fixed and moving image file names
if [[ ! -s $f ]] ; then echo no $f ; exit; fi
if [[ ! -s $m ]] ; then echo no $m ;exit; fi
nm1=` basename $f | cut -d '.' -f 1 `
nm2=` basename $m | cut -d '.' -f 1 `
nm=${D}${nm1}_fixed_${nm2}_moving   # construct output prefix
reg=${AP}antsRegistration           # path to antsRegistration
echo affine $m $f outname is $nm
$reg -d $dim  \
                        -m mattes[  $f, $m , 1 , 32, random , 0.05 ] \
                         -t affine[ 2 ] \
                         -c [1800x1000x1500x20,1.e-8,20]  \
                        -s 4x2x1x0  \
                        -f 8x4x2x1 -l 1 \
                        -m mattes[  $f, $m , 1 , 32 ] \
                         -t syn[ 0.25, 3, 0 ] \
                         -c [30x20x0,1.e-8,20]  \
                        -s 2x1x0  \
                        -f 4x2x1 -l 1 \
                       -u 1 \
                       -o [${nm},${nm}_diff.nii.gz]

${AP}antsApplyTransforms -d 3 -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0Affine.mat
