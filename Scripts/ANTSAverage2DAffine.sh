#!/bin/bash
NUMPARAM=$#
if [ $NUMPARAM -lt 2 ] ; then
echo " Usage "
echo " $0  OUTNameAffine.txt   *Affine.txt "
echo " assumes close to idetity affine transforms "
exit
fi
OUTNM=$1
shift 1
FLIST=$*
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
PARAM5=` awk -v a=$PARAM5 -v b=$NFILES 'BEGIN{print (a / b)}' `

LL=` head -n 4 $FLIST | grep Paramet | cut -d ' ' -f 7  `
for x in $LL ; do PARAM6=` awk -v a=$PARAM6 -v b=$x 'BEGIN{print (a + b)}' `  ; done
PARAM6=` awk -v a=$PARAM6 -v b=$NFILES 'BEGIN{print (a / b)}' `

LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 2  `
for x in $LL ; do PARAM7=` awk -v a=$PARAM7 -v b=$x 'BEGIN{print (a + b)}' `  ; done
PARAM7=` awk -v a=$PARAM7 -v b=$NFILES 'BEGIN{print (a / b)}' `

LL=` cat $FLIST | grep FixedParamet | cut -d ' ' -f 3  `
for x in $LL ; do PARAM8=` awk -v a=$PARAM8 -v b=$x 'BEGIN{print (a + b)}' `  ; done
PARAM8=` awk -v a=$PARAM8 -v b=$NFILES 'BEGIN{print (a / b)}' `

echo "#Insight Transform File V1.0 " > $OUTNM
echo "# Transform 0 " >> $OUTNM
echo "Transform: MatrixOffsetTransformBase_double_2_2  " >> $OUTNM
echo "Parameters:  $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6  " >> $OUTNM
echo "FixedParameters: $PARAM7 $PARAM8 " >> $OUTNM


exit
