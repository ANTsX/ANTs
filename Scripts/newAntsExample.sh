#!/bin/bash
#
dim=2 # image dimensionality
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
                        -m mattes[  $f, $m , 1 , 32 ] \
                         -t exponential[ .5, 3, 0.25 ] \
                         -c [10x10,1.e-8,20]  \
                        -s 2x2  \
                        -f 2x2 -l 1 \
                       -u 1 \
                       -o [${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz]

${AP}antsApplyTransforms -d $dim -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0Affine.mat -o ${nm}_warped.nii.gz
