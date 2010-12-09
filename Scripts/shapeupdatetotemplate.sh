#!/bin/bash

function shapeupdatetotemplate {

    # local declaration of values
    dim=${DIM}
    template=${TEMPLATE}
    templatename=${TEMPLATENAME}
    outputname=${OUTPUTNAME}
    gradientstep=-${GRADIENTSTEP}

# debug only
# echo $dim
# echo ${template}
# echo ${templatename}
# echo ${outputname}
# echo ${outputname}*formed.nii*
# echo ${gradientstep}

# We find the average warp to the template and apply its inverse to the template image
# This keeps the template shape stable over multiple iterations of template building

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 1"
    echo "--------------------------------------------------------------------------------------"
    ${ANTSPATH}AverageImages $dim ${template} 1 ${outputname}*formed.nii.gz

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 2"
    echo "--------------------------------------------------------------------------------------"

	${ANTSPATH}AverageImages $dim ${templatename}warp.nii.gz 0 `ls ${outputname}*Warp.nii.gz | grep -v "InverseWarp"`

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 3"
    echo "--------------------------------------------------------------------------------------"
	${ANTSPATH}MultiplyImages $dim ${templatename}warp.nii.gz ${gradientstep} ${templatename}warp.nii.gz

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 4"
    echo "--------------------------------------------------------------------------------------"
    rm -f ${templatename}Affine.txt

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 5"
    echo "--------------------------------------------------------------------------------------"

    # Averaging and inversion code --- both are 1st order estimates.
    if [ ${dim} -eq 2   ] ; then
      ANTSAverage2DAffine ${templatename}Affine.txt ${outputname}*Affine.txt
    elif [ ${dim} -eq 3  ] ; then
      ANTSAverage3DAffine ${templatename}Affine.txt ${outputname}*Affine.txt
    fi
    ${ANTSPATH}WarpImageMultiTransform ${dim} ${templatename}warp.nii.gz ${templatename}warp.nii.gz -i  ${templatename}Affine.txt -R ${template}
    ${ANTSPATH}WarpImageMultiTransform ${dim} ${template} ${template} -i ${templatename}Affine.txt ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz ${templatename}warp.nii.gz -R ${template}

    echo
    echo "--------------------------------------------------------------------------------------"
    echo " shapeupdatetotemplate 6"
    ${ANTSPATH}MeasureMinMaxMean ${dim} ${templatename}warp.nii.gz ${templatename}warplog.txt 1

}

function ANTSAverage2DAffine {

    OUTNM=${templatename}Affine.txt
    FLIST=${outputname}*Affine.txt
    NFILES=0
    PARAM1=0
    PARAM2=0
    PARAM3=0
    PARAM4=0
    PARAM5=0
    PARAM6=0
    PARAM7=0
    PARAM8=0
    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM1=` awk -v a=$PARAM1 -v b=$x 'BEGIN{print (a + b)}' ` ;  let NFILES=$NFILES+1  ; done
    PARAM1=` awk -v a=$PARAM1 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM2=` awk -v a=$PARAM2 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM2=` awk -v a=$PARAM2 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM3=` awk -v a=$PARAM3 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM3=` awk -v a=$PARAM3 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 5  `
    for x in $LL ; do PARAM4=` awk -v a=$PARAM4 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM4=` awk -v a=$PARAM4 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 6  `
    for x in $LL ; do PARAM5=` awk -v a=$PARAM5 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM5=0 # ` awk -v a=$PARAM5 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 7  `
    for x in $LL ; do PARAM6=` awk -v a=$PARAM6 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM6=0 # ` awk -v a=$PARAM6 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 2  `
    for x in $LL ; do PARAM7=` awk -v a=$PARAM7 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM7=` awk -v a=$PARAM7 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM8=` awk -v a=$PARAM8 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM8=` awk -v a=$PARAM8 -v b=$NFILES 'BEGIN{print (a / b)}' `

    echo "# Insight Transform File V1.0 " > $OUTNM
    echo "# Transform 0 " >> $OUTNM
    echo "Transform: MatrixOffsetTransformBase_double_2_2  " >> $OUTNM
    echo "Parameters:  $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6  " >> $OUTNM
    echo "FixedParameters: $PARAM7 $PARAM8 " >> $OUTNM


}

function ANTSAverage3DAffine {

    OUTNM=${templatename}Affine.txt
    FLIST=${outputname}*Affine.txt
    NFILES=0
    PARAM1=0
    PARAM2=0
    PARAM3=0
    PARAM4=0
    PARAM5=0
    PARAM6=0
    PARAM7=0
    PARAM8=0
    PARAM9=0
    PARAM10=0
    PARAM11=0
    PARAM12=0
    PARAM13=0
    PARAM14=0
    PARAM15=0
    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM1=` awk -v a=$PARAM1 -v b=$x 'BEGIN{print (a + b)}' ` ;  let NFILES=$NFILES+1  ; done
    PARAM1=` awk -v a=$PARAM1 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM2=` awk -v a=$PARAM2 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM2=` awk -v a=$PARAM2 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM3=` awk -v a=$PARAM3 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM3=` awk -v a=$PARAM3 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 5  `
    for x in $LL ; do PARAM4=` awk -v a=$PARAM4 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM4=` awk -v a=$PARAM4 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 6  `
    for x in $LL ; do PARAM5=` awk -v a=$PARAM5 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM5=` awk -v a=$PARAM5 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 7  `
    for x in $LL ; do PARAM6=` awk -v a=$PARAM6 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM6=` awk -v a=$PARAM6 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 8  `
    for x in $LL ; do PARAM7=` awk -v a=$PARAM7 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM7=` awk -v a=$PARAM7 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 9  `
    for x in $LL ; do PARAM8=` awk -v a=$PARAM8 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM8=` awk -v a=$PARAM8 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 10  `
    for x in $LL ; do PARAM9=` awk -v a=$PARAM9 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM9=` awk -v a=$PARAM9 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 11  `
    for x in $LL ; do PARAM10=` awk -v a=$PARAM10 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM10=0 # ` awk -v a=$PARAM10 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 12  `
    for x in $LL ; do PARAM11=` awk -v a=$PARAM11 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM11=0 # ` awk -v a=$PARAM11 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 13  `
    for x in $LL ; do PARAM12=` awk -v a=$PARAM12 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM12=0 # ` awk -v a=$PARAM12 -v b=$NFILES 'BEGIN{print (a / b)}' `

# origin params below

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 2  `
    for x in $LL ; do  PARAM13=` awk -v a=$PARAM13 -v b=$x 'BEGIN{print (a + b)}' ` ;  done
    PARAM13=` awk -v a=$PARAM13 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 3  `
    for x in $LL ; do PARAM14=` awk -v a=$PARAM14 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM14=` awk -v a=$PARAM14 -v b=$NFILES 'BEGIN{print (a / b)}' `

    LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 4  `
    for x in $LL ; do PARAM15=` awk -v a=$PARAM15 -v b=$x 'BEGIN{print (a + b)}' `  ; done
    PARAM15=` awk -v a=$PARAM15 -v b=$NFILES 'BEGIN{print (a / b)}' `

    echo "# Insight Transform File V1.0 " > $OUTNM
    echo "# Transform 0 " >> $OUTNM
    echo "Transform: MatrixOffsetTransformBase_double_3_3  " >> $OUTNM
    echo "Parameters:  $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6  $PARAM7 $PARAM8 $PARAM9 $PARAM10 $PARAM11 $PARAM12  " >> $OUTNM
    echo "FixedParameters: $PARAM13 $PARAM14 $PARAM15 " >> $OUTNM

}


DIM=3
TEMPLATE=t2template.nii.gz
TEMPLATENAME=t2template
OUTPUTNAME=t2
GRADIENTSTEP=0.1
  shapeupdatetotemplate ${DIM} ${TEMPLATE} ${TEMPLATENAME} ${OUTPUTNAME} ${GRADIENTSTEP}
