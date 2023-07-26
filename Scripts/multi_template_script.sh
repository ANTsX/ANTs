#! /bin/bash
ANTSPATH=$1
if [ ${#ANTSPATH} -lt 3 ] ; then
 echo usage
 echo $0 PATH_TO_ANTS_BINARIES image_to_label.nii.gz
 exit
fi
if ! command -v ants.sh &> /dev/null
then
  echo you need the file ants.sh - exiting
  exit
fi

Image_To_Be_Labeled=$2
if [[ ! -s $Image_To_Be_Labeled ]] || [[ ${#Image_To_Be_Labeled} -lt 3 ]] ; then
  echo you need to pass an image you want to label as the 2nd argument
  echo " you tried to pass $Image_To_Be_Labeled  with name length ${#Image_To_Be_Labeled} --- sure that's right? "
  exit
fi

LONGITERATIONS=30x50x30
FASTITERATIONS=10x0x0
ct=0
LIST_OF_IMAGES=" image1.nii.gz image2.nii.gz "
LIST_OF_LABELS=( image1_labels.nii.gz image2_labels.nii.gz )
LIST_OF_OUTPUT=( image1_output image2_output  )
for  labeled_img in $LIST_OF_IMAGES ; do
  sh ants.sh 3  $labeled_img  $Image_To_Be_Labeled ${LIST_OF_OUTPUT[${ct}]} $FASTITERATIONS ${LIST_OF_LABELS[${ct}]}
  let ct=$ct+1
done
ls *labeled.nii.gz > labeled_list.txt
ImageSetStatistics 3 labeled_list.txt  new_subject_labels.nii.gz  1
