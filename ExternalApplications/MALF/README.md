PICSL Multi-Atlas Segmentation Tool
===================================

This package contains source code for joint label fusion [1] and corrective learning [2], which were
applied in MICCAI 2012 Grand Challenge on Multi-Atlas Labeling and finished in the first place.

Joint label fusion is for combining candidate segmentations produced by registering and warping multiple atlases for a target image. Corrective learning can be applied to further reduce systematic errors produced by joint label fusion (see [2] for detail). In general, corrective learning can be applied to correct systematic errors produced by other segmentation methods as well. If you use this software to produce results for a publication, please cite the
following paper(s) accordingly:

Suggested Citation Language:
============================

Multi-atlas label fusion is a segmentation strategy that has been applied to a number of medical image segmentation problems. The method makes use of a set of expert-labeled atlases, where each atlas consists of a sample image and a set of labels for the anatomic structures in that image.  When a new target image is presented for segmentation, each atlas image is registered to the target image.  The deformation fields obtained by registration are then used to propagate the atlas labels to the target image.  Depending on dissimilarities in anatomy and appearance between the atlas and target images, each atlas produces a different segmentation of the target image.  Multi-atlas label fusion combines these results to produce a consensus segmentation of the target image.

Joint label fusion is an extension of multi-atlas label fusion that aims to reduce segmentation bias produced by redundancies in the atlas set. To capture the redundancies, the method accounts for both similarity between each atlas and the target as well as similarity between atlases for label fusion. The expected label error produced by one atlas is large when the image intensity difference between the warped atlas and target image is large. The expectation that any two atlases both produce a label error is large only when the two atlases are similar and both have large intensity differences from the target image. Given the estimated pairwise atlas correlations, joint label fusion computes spatially-varying voting weights for each atlas in a closed form. The final consensus segmentation is produced by weighted voting.

[1] Wang, Hongzhi; Suh, Jung W.; Das, Sandhitsu R.; Pluta, John B.; Craige, Caryne; Yushkevich, Paul A.; , "Multi-Atlas Segmentation with Joint Label Fusion," Pattern Analysis and Machine Intelligence, IEEE Transactions on , vol.35, no.3, pp.611-623, March 2013
doi: 10.1109/TPAMI.2012.143
keywords: {Accuracy;Biomedical imaging;Educational institutions;Image segmentation;Indexes;Joints;Radiology;Multi-atlas label fusion segmentation;dependence;hippocampal segmentation;}
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6226425&isnumber=6461861


[2] H. Wang, S. R. Das, J. W. Suh, M. Altinay, J. Pluta, C. Craige, B. B. Avants, and P. A. Yushkevich,
"A Learning-Based Wrapper Method to Correct Systematic Errors in Automatic Image Segmentation:
Consistently Improved Performance in Hippocampus, Cortex and Brain," Neuroimage, vol. 55, iss. 3,
pp. 968-985, 2011.



INSTRUCTIONS FOR COMPILING THE CODE

A cmake file is provided to facilitate compiling the code. The user needs to set up configurations for
compiling. To do so, go to the folder of the software package. Then type

ccmake .

to setup the enviroment. After this is done, type

make

to complie and build the executable files. Our program uses ITK's I/O functions to handle image input
and output. Hence, it requries prebuilt ITK. The following three executable files will be built,

jointfusion : joint label fusion
bl          : learning classifiers for correcting systematic errors
sa          : apply the learned classifiers to correct systematic segmentation errors on testing images

Type each command to see details on how to use them.


09/04/2012
Hongzhi Wang
wanghongzhi78@gmail.com


===============================

