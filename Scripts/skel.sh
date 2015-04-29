#!/bin/bash
if [[ $# -lt 3 ]] ; then 
  echo read script before using 
  echo $0 segmentationimage outputname thresholdlevel
  echo e.g. $0 seg.nii.gz gm  2 
  exit 
fi
TEM=$1
nm=$2
ThresholdImage 3 $TEM ${nm}_topo.nii.gz $3 $3
ImageMath 3 ${nm}_topo.nii.gz MD ${nm}_topo.nii.gz 1
ImageMath 3 ${nm}_topo.nii.gz ME ${nm}_topo.nii.gz 1
ImageMath 3 ${nm}_topo.nii.gz Neg ${nm}_topo.nii.gz
ImageMath 3 ${nm}_topod.nii.gz D ${nm}_topo.nii.gz
ImageMath 3 ${nm}_topodl.nii.gz Laplacian ${nm}_topod.nii.gz 1.0 1
ThresholdImage 3 ${nm}_topodl.nii.gz speed.nii.gz 0.0 0.25
ImageMath 3 speed.nii.gz Neg speed.nii.gz
ImageMath 3 ${nm}_topo.nii.gz GetLargestComponent ${nm}_topo.nii.gz
ImageMath 3 ${nm}_topo_skel.nii.gz PropagateLabelsThroughMask speed.nii.gz ${nm}_topo.nii.gz 20000 1
ImageMath 3 ${nm}_topo_skel.nii.gz Neg ${nm}_topo_skel.nii.gz
ThresholdImage 3 ${nm}_topo_skel.nii.gz ${nm}_topo_skel.nii.gz 1 10000
MultiplyImages 3 ${nm}_topo_skel.nii.gz $3 ${nm}_topo_skel.nii.gz 
