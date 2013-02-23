## This work is described in detail in http://dx.doi.org/10.1109/TPAMI.2012.143

## ========= Reference text that can be used to describe these methods
#Multi-atlas label fusion is a segmentation strategy that has been applied to a number of medical image segmentation problems. The method makes use of a set of expert-labeled atlases, where each atlas consists of a sample image and a set of labels for the anatomic structures in that image.  When a new target image is presented for segmentation, each atlas image is registered to the target image.  The deformation fields obtained by registration are then used to propagate the atlas labels to the target image.  Depending on dissimilarities in anatomy and appearance between the atlas and target images, each atlas produces a different segmentation of the target image.  Multi-atlas label fusion combines these results to produce a consensus segmentation of the target image.  
#
#Joint label fusion is an extension of multi-atlas label fusion that aims to reduce segmentation bias produced by redundancies in the atlas set. To capture the redundancies, the method accounts for both similarity between each atlas and the target as well as similarity between atlases for label fusion. The expected label error produced by one atlas is large when the image intensity difference between the warped atlas and target image is large. The expectation that any two atlases both produce a label error is large only when the two atlases are similar and both have large intensity differences from the target image. Given the estimated pairwise atlas correlations, joint label fusion computes spatially-varying voting weights for each atlas in a closed form. The final consensus segmentation is produced by weighted voting.    

project(PICSL_MALF)

add_executable(jointfusion ./JointFusion/LabelFusion.cxx)
add_executable(sa ./CorrectiveLearning/segAdapter.cxx)
add_executable(bl ./CorrectiveLearning/BiasLearn.cxx)

target_link_libraries(jointfusion ${ITK_LIBRARIES} )
target_link_libraries(sa ${ITK_LIBRARIES} )
target_link_libraries(bl ${ITK_LIBRARIES} )


install(TARGETS jointfusion sa bl
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
