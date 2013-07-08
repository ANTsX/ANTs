#!/bin/bash
#$ -S /bin/bash

dir=$1
cmd=$2
para=$3
naming=$4
fixname=$5
movingname=$6
echo $6
opta=$7
echo $7
optb=$8
echo $8
optc=$9
echo $9
optd=${10}
echo ${10}
opte=${11}
echo ${11}
optf=${12}
echo ${12}
optg=${13}
echo ${13}
opth=${14}
echo ${14}
opti=${15}
echo ${15}
optj=${16}
echo ${16}
optk=${17}
echo ${17}
optl=${18}
echo ${18}
optm=${19}
echo ${19}
optn=${20}
echo ${20}
opto=${21}
echo ${21}
optp=${22}
echo ${22}
optq=${23}
echo ${23}
optr=${24}
echo ${24}
opts=${25}
echo ${25}
optt=${26}
echo ${26}

echo " dir "
echo $1

cd $1
echo $cmd $para $naming $fixname $movingname $opta $optb $optc $optd $opte  $optf $optg $opth $opti $optj $optk $optl $optm $optn $opto $optp $optq $optr $opts $optt
$cmd $para $naming $fixname $movingname $opta $optb $optc $optd $opte   $optf $optg $opth $opti $optj $optk $optl $optm  $optn $opto $optp $optq $optr $opts $optt ${27} ${28} ${29} ${30} ${31} ${32} ${33} ${34} ${35} ${36} ${37} ${38} ${39}  ${40} ${41} ${42} ${43} ${44} ${45} ${46} ${47} ${48} ${49}
cd -
