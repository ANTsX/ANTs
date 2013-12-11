#!/bin/bash

count=0
for x in  phantomA phantomB phantomC phantomD phantomE phantomF phantomG phantomH
do
    count=`expr $count + 1`					# Increment the counter
    echo " $x "
  echo " LaplacianThickness ${x}wm.hdr ${x}gm.hdr 15 TEST{$count}thick.hdr  1 0.5   "
    ThresholdImage 2 ${x}wmgm.jpg ${x}wm.nii 245 290
    ThresholdImage 2 ${x}wmgm.jpg ${x}gm.nii 111  133
    LaplacianThickness ${x}wm.nii ${x}gm.nii 15 TEST{$count}thick.hdr  5 3
   WarpImageMultiTransform 2  TEST{$count}thick.hdr  TEST{$count}thickreg.hdr   TEST{$count}Warp.nii   TEST{$count}Affine.txt
#    SmoothImage  2 TEST{$count}thickreg.hdr  1  TEST{$count}thickreg.hdr
done

StudentsTestOnImages 2  THICK2.hdr 4 4


 ~/Code/bin/ants/GLM 2 mask.nii designmat.txt contrastvec.txt temp.hdr 1000   TEST{1}thickreg.hdr TEST{2}thickreg.hdr TEST{3}thickreg.hdr TEST{4}thickreg.hdr TEST{5}thickreg.hdr TEST{6}thickreg.hdr TEST{7}thickreg.hdr TEST{8}thickreg.hdr

MultiplyImages 2  THICK2.hdr phantomtemplate.jpg THICK2.hdr
