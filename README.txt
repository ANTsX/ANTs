Release 1.8

Homepage: http://www.picsl.upenn.edu/ANTS/

Introduction -- ANTS is a tool for computational neuroanatomy based on
medical images.  ANTS reads any image type that can be read by ITK
(www.itk.org), that is, jpg, tiff, hdr, nii, nii.gz, mha/d and more
image types as well.  For the most part, ANTS will output float images
which you can convert to other types with the ANTS
ConvertImagePixelType tool.  ImageMath has a bunch of basic utilities
such as multiplication, inversion and many more advanced tools such as
computation of the Lipschitz norm of a deformation field.  ANTS
programs may be called from the command line on almost any platform
.... you can compile the code yourself or use the precompiled binaries
for Windows (Vista), OSX (Darwin) or linux (32 bit or 64 bit).
Download the binaries that are correct for you.  If you are a naive
user (not from a medical imaging background) then you might still find
the tools here useful.  Many of the operations available, for
instance, in PhotoShop are available in ANTS and many more are
available as well.  For instance, ANTS tools may be useful in face
mapping / morphing and also in generating animations from two
different images, for instance, interpolating between frames in a
movie.  But, mainly, ANTS is useful for brain mapping, segmentation,
measuring cortical thickness and in generating automated or
semi-automated labeling of three-dimensional imagery (e.g. labeling
hippocampus or cortical regions or lobes of the lung).  Many
prior-based segmentation possibilities are available in the Atropos
tool, including three tissue segmentation, structure-specific
segmentation and brain extracton.

The ants.pdf file has more details and examples.

New Stuff 1.8 :

Sped up CC metric -- comparable to PR but faster.

Fixed the MI metric -- fast and functional.

WarpTimeSeries --- for deforming 4D or vector images.

New Stuff 1.7 :

Now using SymmetricSecondOrderPixelType -- fixes some DT bugs with Nifti I/O etc.

This means our nii tensors have SYMMATRIX as intent (as with standard)

Update in parameters for convergence and inversion -- aids performance.

TensorToVector coloring --- also preliminary integration of vector field

Speed up number one for the CC metric (number two coming later).

Atropos !  new tool for segmentation.

Nick's N4 bias correction tool (publication on the way).

Updates to buildtemplateparallel that allow parallel use on multicore machines.

Various utilities and a few improvements in usage.

New Stuff 1.6 :

Check DT tensor ordering in DTI Read/Write

HistogramMatching in ImageMath

ConvertImagePixelType utility

Affine averaging in buildtemplateparallel

Updated time-dependent diffeomorphic mapping (option --geodesic  1 / 2 ).

Updated with a greedy exponential mapping diffeomorphic approach akin to DiffeomorphicDemons.

Bug fix in checking ANTS convergence.

Other miscellaneous including minor Apocrita changes and allowing spaces in command line interface.

# directory guide:
Documentation -- pdf / tex describing ANTS
Examples    -- the executable programs and test data in Examples/Data
Scripts -- user-friendly scripts for template building and running studies
Utilities --- basic utilities
ImageRegistration	-- base code for ImageRegistration
Temporary  -- where temporary code lives
Tensor  -- base code for diffusion tensor operations


Use cmake (cmake.org) to set up compilation.
To build ANTS, do the following:

1.  get ANTS

svn checkout https://advants.svn.sourceforge.net/svnroot/advants ANTS


2.  get itk :

 cvs -d :pserver:anoncvs@www.itk.org:/cvsroot/Insight co Insight

3.  compile itk and ANTS -- link ANTS to itk build directory

  ccmake ANTS/Examples/

4.   call  ctest in the compile directory and verify that the tests pass

5.  in general, to perform a mapping :

# include the mask, if desired.  mask is inclusive.

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

# warp the template image to tp22 -- note reversal of order from above

WarpImageMultiTransform 3 template.nii.gz templatetotp22.nii  -R tp22_s1.nii   tp22mapWarp.nii tp22mapAffine.txt

# or call ants.sh for a standard approach.

#
#  use CreateJacobianDeterminantImage to get log-Jacobian (volumetric change) images
#  and programs StudentsTestOnImages or  GLM to peform a statistical study
#  one might also use SurfaceCurvature to do a curvature study
#  or LaplacianThickness to do a thickness study
#
