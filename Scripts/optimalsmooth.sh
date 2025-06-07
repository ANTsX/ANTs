#!/bin/bash

echo  " sh $0 outputprefix "
if [ $# -lt 1 ] ; then
exit
fi

SUB=$1
TH=${SUB}thicknorm.nii.gz
JA=${SUB}logjacobian.nii
MK=${SUB}mask.nii.gz
OUTT=${SUB}smooththick.nii.gz
OUTJ=${SUB}smoothjac.nii.gz
REPS=5
  SurfaceBasedSmoothing $TH 1.5 $MK $OUTT $REPS
  SurfaceBasedSmoothing $JA 1.5 $MK $OUTJ $REPS
