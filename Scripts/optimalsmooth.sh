

# usage :  SurfaceBasedSmoothing ImageToSmooth  sigma SurfaceImage  outname  {numrepeatsofsmoothing}

export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"}
SUB=$1
TH=${SUB}thicknorm.nii.gz
JA=${SUB}logjacobian.nii
MK=${SUB}mask.nii.gz
OUTT=s${SUB}thick.nii.gz
OUTJ=s${SUB}jac.nii.gz
REPS=5
  ${ANTSPATH}SurfaceBasedSmoothing $TH 1.5 $MK $OUTT $REPS
  ${ANTSPATH}SurfaceBasedSmoothing $JA 1.5 $MK $OUTJ $REPS
