#!/bin/bash
NUMPARAMS=$#
if [ $NUMPARAMS -lt 3 ]
then
echo " USAGE ::  "
echo " $0   inWarp WhichTypeOfMagnitude outputmagnitudeimage.nii "
echo " WhichTypeOfMagnitude ==  Laplacian Grad  Euclidean "
echo " the numbers that come out of this depend on some itk filters being dimensionally correct and consistent. "
echo " this may not be the case -- however, relative values should be ok , e.g comparing Grad values to Grad values etc "
exit
fi

inwarp=$1
whichmag=$2
magimage=$3
dimarr=(`PrintHeader ${inwarp}xvec.nii.gz  | grep Spacing | grep , | wc `)
DIM=${dimarr[0]}
let DIM=${DIM}+1
# multiply some image of the right size by 0 to make the initial magnitude image
ImageMath  $DIM  $magimage m  ${inwarp}xvec.nii.gz 0
ImageMath  $DIM  ${magimage}temp.nii.gz m  ${inwarp}xvec.nii.gz 0
ok=0
for  x in  x y z
do
if [ -s  ${inwarp}${x}vec.nii.gz ]
then
# compute the magnitude of the x-component deformation
if [ $whichmag == Laplacian ]
then
echo " Laplacian $x "
ok=1
ImageMath  $DIM  ${magimage}temp.nii.gz Laplacian  ${inwarp}${x}vec.nii.gz 1.0 0
elif [ $whichmag == Grad ]
then
echo " Grad $x "
ok=1
ImageMath  $DIM  ${magimage}temp.nii.gz Grad ${inwarp}${x}vec.nii.gz 1.0  0
#ImageMath  $DIM  temp.nii.gz m temp.nii.gz temp.nii.gz
elif [ $whichmag == Euclidean ]
then
echo " Euclidean $x "
ok=1
ImageMath  $DIM  ${magimage}temp.nii.gz m  ${inwarp}${x}vec.nii.gz ${inwarp}${x}vec.nii.gz
fi
# add the xvec magnitude to the magimage
ImageMath  $DIM  $magimage  +  $magimage ${magimage}temp.nii.gz
fi # ifz
done #loop
if [ $ok -eq 1 ]
then
# take the square root
ImageMath  $DIM  $magimage ^ $magimage 0.5
MeasureMinMaxMean $DIM $magimage
#thus,
#magimage=sqrt( x*x + y*y + z*z )
else
echo " User Requested deformation measure ::   $whichmag "
echo " that is not available -- user needs to choose another type of deformation measure "
fi

rm -f  ${magimage}temp.nii.gz
