#!/bin/bash

# trap keyboard interrupt (control-c)
trap control_c SIGINT

function setPath {
    cat <<SETPATH

--------------------------------------------------------------------------------------
Error locating ANTS
--------------------------------------------------------------------------------------
It seems that the ANTSPATH environment variable is not set. Please add the ANTSPATH
variable. This can be achieved by editing the .bash_profile in the home directory.
Add:

ANTSPATH=/home/yourname/bin/ants/

Or the correct location of the ANTS binaries.

Alternatively, edit this script ( `basename $0` ) to set up this parameter correctly.

SETPATH
    exit 1
}

# Uncomment the line below in case you have not set the ANTSPATH variable in your environment.
# export ANTSPATH=${ANTSPATH:="$HOME/bin/ants/"} # EDIT THIS

#ANTSPATH=YOURANTSPATH
if [[ ${#ANTSPATH} -le 3 ]];
  then
    setPath >&2
  fi

WARP=${ANTSPATH}/antsApplyTransforms
AVERAGE_AFFINE_PROGRAM=${ANTSPATH}/AverageAffineTransform # NoRigid

if [[ ! -s ${WARP} ]];
  then
    echo "antsApplyTransforms program can't be found. Please (re)define \$ANTSPATH in your environment."
    exit
  fi

if [[ ! -s ${AVERAGE_AFFINE_PROGRAM} ]];
  then
    echo "AverageAffineTransform* program can't be found. Please (re)define \$ANTSPATH in your environment."
    exit
  fi


DIM=3
TEMPLATE=template0.nii.gz
OUTPUTNAME=T_
GRADIENTSTEP=0.25
whichtemplate=0
statsmethod=1
useaff=0
USAGE="$0 -d 3 -t template0.nii.gz -o T_ -g 0.25 -s 1 -w 0 -y 0"
while getopts "d:t:o:g:w:s:y:h:" OPT
  do
  case $OPT in
      h) #help
   echo $USAGE
   exit 0
   ;;
      d)  # dimensions
   DIM=$OPTARG
   ;;
      t)  # name of image
   TEMPLATE=$OPTARG
   ;;
      o)  # output prefix
   OUTPUTNAME=$OPTARG
   ;;
      g)  # moving image
   GRADIENTSTEP=$OPTARG
   ;;
      w)  # for multivar templates
   whichtemplate=$OPTARG
   ;;
    y)  # affine only in update?
    useaff=$OPTARG
    ;;
      s)  # median, mean, etc
   statsmethod=$OPTARG
   ;;
     \?) # getopts issues an error message
   echo "$USAGE" >&2
   exit 1
   ;;
  esac
done

if [[ $useaff -eq 1 ]] ; then
  AVERAGE_AFFINE_PROGRAM=${ANTSPATH}/AverageAffineTransformNoRigid
fi


function summarizeimageset() {

  local dim=$1
  shift
  local output=$1
  shift
  local method=$1
  shift
  local images=( "${@}" "" )

  case $method in
    0) #mean
      AverageImages $dim $output 0 ${images[*]}
      ;;
    1) #mean of normalized images
      AverageImages $dim $output 1 ${images[*]}
      ;;
    2) #median
      for i in "${images[@]}";
        do
          echo $i >> ${output}_list.txt
        done
      ImageSetStatistics $dim ${output}_list.txt ${output} 0
      rm ${output}_list.txt
      ;;
    3) #median+sharpen
      for i in "${images[@]}";
        do
          echo $i >> ${output}_list.txt
        done
      ImageSetStatistics $dim ${output}_list.txt ${output} 0
      ImageMath $dim ${output} Sharpen ${output}
      rm ${output}_list.txt
      ;;
      4) #mean+sharpen
        AverageImages $dim $output 1 ${images[*]}
        ImageMath $dim ${output} Sharpen ${output}
        ;;
  esac

  }


function shapeupdatetotemplate() {

   echo "shapeupdatetotemplate()"

    # local declaration of values
    dim=$1
    template=$2
    templatename=$3
    outputname=$4
    gradientstep=-$5
    whichtemplate=$6
    statsmethod=$7

# debug only
# echo $dim
# echo ${template}
# echo ${templatename}
# echo ${outputname}
# echo ${outputname}*WarpedToTemplate.nii*
# echo ${gradientstep}

# We find the average warp to the template and apply its inverse to the template image
# This keeps the template shape stable over multiple iterations of template building

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate---voxel-wise averaging of the warped images to the current template"
    date
    #echo "   ${ANTSPATH}/AverageImages $dim ${template} 1 ${templatename}${whichtemplate}*WarpedToTemplate.nii.gz    "
    #echo "    ${ANTSPATH}/ImageSetStatistics $dim ${whichtemplate}WarpedToTemplateList.txt ${template} 0"
    echo "--------------------------------------------------------------------------------------"
    imagelist=(`ls ${outputname}template${whichtemplate}*WarpedToTemplate.nii.gz`)
    if [[ ${#imagelist[@]} -eq 0  ]] ; then
      echo ERROR shapeupdatedtotemplate - imagelist length is 0
      exit 1
    fi

    summarizeimageset $dim $template $statsmethod ${imagelist[@]}

    WARPLIST=( `ls ${outputname}*[0-9]Warp.nii.gz 2> /dev/null` )
    NWARPS=${#WARPLIST[*]}
    echo "number of warps = $NWARPS"
    echo "$WARPLIST"

    if [[ $whichtemplate -eq 0 ]];
      then

        if [[ $NWARPS -ne 0 ]]; then
          echo "$NWARPS does not equal 0"
          echo
          echo "--------------------------------------------------------------------------------------"
          echo " shapeupdatetotemplate---voxel-wise averaging of the inverse warp fields (from subject to template)"
          echo "   ${ANTSPATH}/AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`"
          date
          echo "--------------------------------------------------------------------------------------"
          ${ANTSPATH}/AverageImages $dim ${templatename}${whichtemplate}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`

          echo
          echo "--------------------------------------------------------------------------------------"
          echo " shapeupdatetotemplate---scale the averaged inverse warp field by the gradient step"
          echo "   ${ANTSPATH}/MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz"
          date
          echo "--------------------------------------------------------------------------------------"
          ${ANTSPATH}/MultiplyImages $dim ${templatename}${whichtemplate}warp.nii.gz ${gradientstep} ${templatename}${whichtemplate}warp.nii.gz
        fi

        echo
        echo "--------------------------------------------------------------------------------------"
        echo " shapeupdatetotemplate---average the affine transforms (template <-> subject)"
        echo "                      ---transform the inverse field by the resulting average affine transform"
        echo "   ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0GenericAffine.mat ${outputname}*GenericAffine.mat"
        echo "   ${WARP} -d ${dim} -e vector -i ${templatename}0warp.nii.gz -o ${templatename}0warp.nii.gz -t [${templatename}0GenericAffine.mat,1] -r ${template}"
        echo "--------------------------------------------------------------------------------------"

        ${AVERAGE_AFFINE_PROGRAM} ${dim} ${templatename}0GenericAffine.mat ${outputname}*GenericAffine.mat

        if [[ $NWARPS -ne 0 ]];
          then
            ${WARP} -d ${dim} -e vector -i ${templatename}0warp.nii.gz -o ${templatename}0warp.nii.gz -t [${templatename}0GenericAffine.mat,1] -r ${template}
            ${ANTSPATH}/MeasureMinMaxMean ${dim} ${templatename}0warp.nii.gz ${templatename}warplog.txt 1
          fi
      fi

    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate---warp each template by the resulting transforms"
    echo "   ${WARP} -d ${dim} --float $USEFLOAT -i ${template} -o ${template} -t [${templatename}0GenericAffine.mat,1] -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -r ${template}"
    echo "--------------------------------------------------------------------------------------"

    if [ -f "${templatename}0warp.nii.gz" ];
      then
        ${WARP} -d ${dim} --float $USEFLOAT -i ${template} -o ${template} -t [${templatename}0GenericAffine.mat,1] -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -t ${templatename}0warp.nii.gz -r ${template}
      else
        ${WARP} -d ${dim} --float $USEFLOAT -i ${template} -o ${template} -t [${templatename}0GenericAffine.mat,1] -r ${template}
      fi

}

TEMPLATENAME=$OUTPUTNAME
shapeupdatetotemplate ${DIM} ${TEMPLATE} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP} ${whichtemplate} $statsmethod
