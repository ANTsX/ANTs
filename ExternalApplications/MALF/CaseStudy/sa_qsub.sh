#!/bin/bash -x
#$ -S /bin/bash

### change these paths to your local paths
ROOT=/home/hwang3/MA2012 # path of the folder holding the 2012 multi atlas label challenge data
BIN=/home/hwang3/MA2012/PICSL_MALF # path to the PICSL_MALF software

ID=$1

$BIN/sa $ROOT/malf/${ID}_JointLabel_R2S3Sig2Lamda0.1.nii.gz $ROOT/malf/BL/JointLabel_BL $ROOT/malf/${ID}_JointLabel_BC.nii.gz -f Testing/$ID.nii.gz malf/${ID}_JointLabel_R2S3Sig2Lamda0.1_posterior%04d.nii.gz 
##-p $ROOT/malf/${ID}_JointLabel_R2S3Sig2Lamda0.1_BC_posterior%04d.nii.gz


