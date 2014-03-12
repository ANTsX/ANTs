#!/bin/bash

echo  " sh $0 outputprefix "
if [ $# -lt 1 ] ; then
exit
fi

export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"}
SUB=$1
TH=${SUB}thicknorm.nii.gz
JA=${SUB}logjacobian.nii
MK=${SUB}mask.nii.gz
OUTT=${SUB}smooththick.nii.gz
OUTJ=${SUB}smoothjac.nii.gz
REPS=5
  ${ANTSPATH}SurfaceBasedSmoothing $TH 1.5 $MK $OUTT $REPS
  ${ANTSPATH}SurfaceBasedSmoothing $JA 1.5 $MK $OUTJ $REPS
