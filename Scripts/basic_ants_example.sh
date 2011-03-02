#!/bin/bash
#
#
if [ $# -lt 5 ] ; then
echo this script will allow you to run \( Reg \) a pairwise segmentation,  \( RegSegPSE or RegSegMSQ \) or  run segmentation, then use segmentation output in registration or \( BTP \) build a template.
echo you call the script like this
echo $0 ImageDimensionality Opt OutputPrefix  r16slice.nii.gz r64slice.nii.gz
echo where Opt is Reg , RegSeg , BTP
echo
echo should work for 2D or 3D images.
echo
echo these are BASIC examples --- not intended to illustrate optimal usage but provide general, reasonable settings.
exit
fi
#
if [ ${#ANTSPATH} -le 3 ] ; then
  echo we guess at your ants path
  export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS
fi
if [ ! -s ${ANTSPATH}/ANTS ] ; then
  echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
#
# in general, you should call ANTS -h for help
#
DIM=$1
OPT=$2
OUTPUTNAME=$3 # choose the output name
II=$4 # change to your own images
JJ=$5
if [[ ! -s $II ]] ; then
  echo cannot find $II --- please change the filename by editing this script.
  exit
fi
if [[ ! -s $JJ ]] ; then
  echo cannot find $JJ --- please change the filename by editing this script.
  exit
fi

ITS=" -i 100x100x30 " # 3 optimization levels

# different transformation models you can choose
TSYNWITHTIME=" -t SyN[0.25,5,0.01] -r DMFFD[10x10,0,3] " # spatiotemporal (full) diffeomorphism
GGREEDYSYN=" -t SyN[0.15]  -r Gauss[3,0] " # fast symmetric normalization gaussian regularization
BGREEDYSYN=" -t SyN[0.15]  -r DMFFD[4x4,0,3] " # fast symmetric normalization dmffd regularization
TELAST=" -t Elast[1] -r Gauss[0.5,3] "             # elastic
TEXP=" -t Exp[0.001,100] -r DMFFD[3,0] "            # exponential

# different metric choices for the user
INTMSQ=" -m MSQ[${II},${JJ},1,0] "
INTMI=" -m MI[${II},${JJ},1,16] "
INTCC=" -m CC[${II},${JJ},1,4] "
#
#
# these are the forward and backward warps.
#
INVW=" -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz "
FWDW=" ${OUTPUTNAME}Warp.nii.gz ${OUTPUTNAME}Affine.txt "
#
INT=$INTCC  # choose a  metric ( here , cross-correlation )
TRAN=$BGREEDYSYN  # choose a transformation
# TRAN=$TEXP  # choose a transformation
if [[ $OPT == "Reg" ]] ; then

# run the registration
${ANTSPATH}/ANTS $DIM -o $OUTPUTNAME $ITS $TRAN $INT
# this is how you apply the output transformation
${ANTSPATH}/WarpImageMultiTransform $DIM ${II} ${OUTPUTNAME}IItoJJ.nii.gz -R ${JJ} $INVW
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJ} ${OUTPUTNAME}JJtoII.nii.gz -R ${II} $FWDW

elif [[ $OPT == "RegSegPSE" ]] ; then

# try something different -- segment the images and use them in ANTS
# to run multivariate registration.
IISEG=${OUTPUTNAME}A_seg.nii.gz
JJSEG=${OUTPUTNAME}B_seg.nii.gz
#
# get masks for atropos
#
${ANTSPATH}ThresholdImage $DIM $II $IISEG 1 1.e9
${ANTSPATH}ThresholdImage $DIM $JJ $JJSEG 1 1.e9
AtroposParams=" -d $DIM -m [0.1,1x1] -c [5,0] -i kmeans[3] "
${ANTSPATH}/Atropos $AtroposParams -a $II -o $IISEG  -x $IISEG
${ANTSPATH}/Atropos $AtroposParams -a $JJ -o $JJSEG  -x $JJSEG
# compute some segmentations and use them in labelguided mapping
LABELGUIDED=" -m PSE[${II},${JJ},${IISEG},${JJSEG},0.75,0.1,25,0,10] "
#
${ANTSPATH}/ANTS $DIM -o ${OUTPUTNAME} $ITS $TRAN $INT $LABELGUIDED
${ANTSPATH}/WarpImageMultiTransform $DIM ${II} ${OUTPUTNAME}IItoJJ.nii.gz -R ${JJ} $INVW
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJ} ${OUTPUTNAME}JJtoII.nii.gz -R ${II} $FWDW

# now warp the labels in both directions
${ANTSPATH}/WarpImageMultiTransform $DIM ${IISEG} ${OUTPUTNAME}IIsegtoJJseg.nii.gz -R ${JJ} $INVW --use-NN
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJSEG} ${OUTPUTNAME}JJsegtoIIseg.nii.gz -R ${II} $FWDW --use-NN


elif [[ $OPT == "RegSegMSQ" ]] ; then

# try something different -- segment the images and use them in ANTS
# to run multivariate registration.
IISEG=${OUTPUTNAME}A_seg.nii.gz
JJSEG=${OUTPUTNAME}B_seg.nii.gz
#
# get masks for atropos
#
${ANTSPATH}ThresholdImage $DIM $II $IISEG 1 1.e9
${ANTSPATH}ThresholdImage $DIM $JJ $JJSEG 1 1.e9
AtroposParams=" -d $DIM -m [0.1,1x1] -c [5,0] -i kmeans[3] "
${ANTSPATH}/Atropos $AtroposParams -a $II -o $IISEG  -x $IISEG
${ANTSPATH}/Atropos $AtroposParams -a $JJ -o $JJSEG  -x $JJSEG
# compute some segmentations and use them in labelguided mapping
LABELGUIDED=" -m MSQ[${IISEG},${JJSEG},0.75] "
#
${ANTSPATH}/ANTS $DIM -o ${OUTPUTNAME} $ITS $TRAN $INT $LABELGUIDED
${ANTSPATH}/WarpImageMultiTransform $DIM ${II} ${OUTPUTNAME}IItoJJ.nii.gz -R ${JJ} $INVW
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJ} ${OUTPUTNAME}JJtoII.nii.gz -R ${II} $FWDW

# now warp the labels in both directions
${ANTSPATH}/WarpImageMultiTransform $DIM ${IISEG} ${OUTPUTNAME}IIsegtoJJseg.nii.gz -R ${JJ} $INVW --use-NN
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJSEG} ${OUTPUTNAME}JJsegtoIIseg.nii.gz -R ${II} $FWDW --use-NN


elif [[ $OPT == "BTP" ]] ; then
#
# now build a template from the two input images .
# in a real application, you would modify the call to btp
# such that it's more appropriate for your data.
# call buildtemplateparallel.sh -h to get usage.
#
${ANTSPATH}buildtemplateparallel.sh -d 2  -o ${OUTPUTNAME}BTP -c 0  $II $JJ
rm -r GR* *cfg

TEM=${OUTPUTNAME}BTPtemplate.nii.gz
NM1=` echo $JJ | cut -d '.' -f 1 `
INVW=" -i ${OUTPUTNAME}BTP${NM1}Affine.txt ${OUTPUTNAME}BTP${NM1}InverseWarp.nii.gz "
FWDW=" ${OUTPUTNAME}BTP${NM1}Warp.nii.gz ${OUTPUTNAME}BTP${NM1}Affine.txt "
${ANTSPATH}/WarpImageMultiTransform $DIM ${JJ} ${OUTPUTNAME}JJtoTemplate.nii.gz -R $TEM $FWDW
${ANTSPATH}/WarpImageMultiTransform $DIM $TEM ${OUTPUTNAME}TemplatetoJJ.nii.gz -R ${JJ} $INVW

else
 echo unrecognized option $OPT
fi
