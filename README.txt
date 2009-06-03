Release 1.0

Homepage: http://www.picsl.upenn.edu/ANTS/


# directory guide:
Examples    -- the executable programs and test data in Examples/Data
Scripts -- user-friendly scripts for template building and running studies
Utilities --- basic utilities
ImageRegistration	-- base code for ImageRegistration
Temporary  -- where temporary code lives
Tensor  -- base code for diffusion tensor operations


Use cmake (cmake.org) to set up compilation.
To build ANTS, do the following:

1.  get ANTS

cvs -d:pserver:anonymous@advants.cvs.sourceforge.net:/cvsroot/advants co ANTS


2.  get itk :

 cvs -d :pserver:anoncvs@www.itk.org:/cvsroot/Insight co Insight

3.  compile itk and ANTS -- link ANTS to itk build directory

  ccmake ANTS/Examples/

4.   call  ctest in the compile directory and verify that the tests pass

5.  in general, to perform a mapping :

# include the mask, if desired.  mask in inclusive.

ANTS 3 -m PR[tp22_s1.nii,template.nii.gz,1,4] -i 50x20x10 -o tp22map   -t SyN[0.25]  -x mask.nii.gz  -r Gauss[3,0]


# The ANTS executable reflects the variational optimization problem
# which balances regularization of the transformation model's parameters
# and the quality of matchins as driven by a similarity (or data) term
#
# explanation :   -m  PR  -- the similarity metric =>  PR[fixed.nii,moving,nii,weight,metric-radius]
#             :   -i 50x20 --  the number of iterations and number of resolution levels
#             :   -o  tp22map --  the output naming convention (can add an extension)
#             :   -t  SyN/Elast/Exp/Syn[time] --- transformation model
#             :   -r  Gauss/Bspline  -- the regularization models
#                                       Gauss[gradient-regularize,deformation-regularize]
#             :   -x  mask -- an inclusive mask -- dictates what information to use in registration
#                          -- defined in the fixed domain but works on both domains
#             :    -m  other metrics :  PSE MSQ MI  etc -- some are label-image (or point-set) metrics
#                                       and some are intensity metrics
#
# Call ANTS with no params to get detailed help
#
# warp the tp22 to template image

WarpImageMultiTransform 3 tp22_s1.nii tp22totemplate.nii  -R template.nii.gz -i  tp22mapAffine.txt tp22mapInverseWarp.nii

#
#  use CreateJacobianDeterminantImage to get log-Jacobian (volumetric change) images
#  and programs StudentsTestOnImages or  GLM to peform a statistical study
#  one might also use SurfaceCurvature to do a curvature study
#  or LaplacianThickness to do a thickness study
#
