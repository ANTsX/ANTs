#!/bin/bash
#
dim=2 # image dimensionality
AP="" # /home/yourself/code/ANTS/bin/bin/  # path to ANTs binaries
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2    # fixed and moving image file names
if [[ ! -s $f ]] ; then echo no fixed $f ; exit; fi
if [[ ! -s $m ]] ; then echo no moving $m ;exit; fi
nm1=` basename $f | cut -d '.' -f 1 `
nm2=` basename $m | cut -d '.' -f 1 `
nm=${D}${nm1}_fixed_${nm2}_moving   # construct output prefix
reg=${AP}antsRegistration           # path to antsRegistration
echo affine $m $f outname is $nm
$reg -d $dim -r [ $f, $m ,1]  \
                        -m mattes[  $f, $m , 1 , 32, regular, 0.1 ] \
                         -t affine[ 1 ] \
                         -c [10000x1000x1000,1.e-8,20]  \
                        -s 4x2x1  \
                        -f 6x4x2 -l 1 \
                        -m cc[  $f, $m , 1 , 2 ] \
                         -t syn[ .3, 3, 0.0 ] \
                         -c [100x50x50,1.e-8,20]  \
                        -s 2x1x0  \
                        -f 4x2x1 -l 1 -u 1 -b 0 -z 1 \
                       -o [${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz]
echo CORRELATION `ImageMath $dim x NeighborhoodCorrelation $f ${nm}_diff.nii.gz`
${AP}antsApplyTransforms -d $dim -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0MatrixOffset.mat -o ${nm}_warped.nii.gz
