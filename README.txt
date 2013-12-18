See the *new* ants [website](http://stnava.github.io/ANTs/ "ANTs")

[Installation/Compilation](http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/ "Build")

[Old Homepage](http://www.picsl.upenn.edu/ANTS/)

Experimental transition to github for development: Wed Jan 23 10:46:13 EST 2013

Release 1.9.x --- final svn release = 1781, now moved to git ....

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

New Stuff 1.9.x :  Requires git-itk.  must be compiled with USE_REVIEW on.

official release of atropos segmentation tool.

checked compilation on windows os.

we now save vector nii.gz files and do not store component images.

improves ImageMath tensor functions.

many operations may be performed on vector images.

some enhancements include:

you can warp a vector field represented as  a nii.gz   via WarpImageMultiTransform

MeasureMinMaxMean and MultiplyImages support vector images.

the only functionality that we "lost" is the ability to use a bspline interpolator with WIMT.   the BSplineInterpolateImageFunction doesnt support vector valued pixels.

Atropos updated and validated.

New Stuff 1.9.2 :

New atropos interface + ROIStatistics in ImageMath

New Stuff 1.9.1 :

Atropos refactored , vtk dependencies allowed , additional tools for surface-based mapping (not much tested), augmented warping for vtk files

Must compile ITK with USE_REVIEW_STATISTICS ON if you want Atropos functionality
Should compile ITK with USE_REVIEW ON
Should compile ITK with USE_OPTIMIZED_REGISTRATION ON
Should have ITK v 3.20 or greater.

New Stuff 1.9 :

Atropos revisions + various utilities.

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

Nick's N4 bias correction tool.

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
ImageRegistration    -- base code for ImageRegistration
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

References:

ANTs registration

  Principal references
    http://www.ncbi.nlm.nih.gov/pubmed/20851191
    Avants BB, Tustison NJ, Song G, Cook PA, Klein A, Gee JC.
    A reproducible evaluation of ANTs similarity metric performance in brain
    image registration.
    Neuroimage. 2011 Feb 1;54(3):2033-44. doi: 10.1016/j.neuroimage.2010.09.025.

    http://link.springer.com/chapter/10.1007/978-3-642-31340-0_28
    Brian B. Avants, Nicholas J. Tustison, Gang Song, Baohua Wu, Michael Stauffer,
    Matthew McCormick, Hans J. Johnson, James C. Gee.
    A Unified Image Registration Framework for ITK.
    Proceedings of the Fifth Workshop on Biomedical Image Registration 2012:266-275.

  Symmetric Normalization (SyN)
    http://www.ncbi.nlm.nih.gov/pubmed/17659998
    Avants BB, Epstein CL, Grossman M, Gee JC.
    Symmetric diffeomorphic image registration with cross-correlation: evaluating
    automated labeling of elderly and neurodegenerative brain.
    Med Image Anal. 2008 Feb;12(1):26-41.

  B-spline-based
    http://www.ncbi.nlm.nih.gov/pubmed/19171516
    Tustison NJ, Avants BB, Gee JC.
    Directly manipulated free-form deformation image registration.
    IEEE Trans Image Process. 2009 Mar;18(3):624-35. doi: 10.1109/TIP.2008.2010072.

    http://link.springer.com/chapter/10.1007%2F978-3-642-31340-0_4
    Nicholas J. Tustison, Brian B. Avants:
    Diffeomorphic Directly Manipulated Free-Form Deformation Image Registration
    via Vector Field Flows.
    Proceedings of the Fifth Workshop on Biomedical Image Registration 2012:31-39.

  Point set registration
     Point-set expectation (PSE)
     http://www.ncbi.nlm.nih.gov/pubmed/19437413
     Pluta J, Avants BB, Glynn S, Awate S, Gee JC, Detre JA.
     Appearance and incomplete label matching for diffeomorphic template based
     hippocampus segmentation.
     Hippocampus. 2009 Jun;19(6):565-71. doi: 10.1002/hipo.20619.

     Havrda-Charvat-Tsallis (JTB)
     http://www.ncbi.nlm.nih.gov/pubmed/20937578
     Tustison NJ, Awate SP, Song G, Cook TS, Gee JC.
     Point set registration using Havrda-Charvat-Tsallis entropy measures.
     IEEE Trans Med Imaging. 2011 Feb;30(2):451-60. doi: 10.1109/TMI.2010.2086065.

Template construction
  http://www.ncbi.nlm.nih.gov/pubmed/15501083
  Avants B, Gee JC.
  Geodesic estimation for large deformation anatomical shape averaging and interpolation.
  Neuroimage. 2004;23 Suppl 1:S139-50.

  http://www.ncbi.nlm.nih.gov/pubmed/19818860
  Avants BB, Yushkevich P, Pluta J, Minkoff D, Korczykowski M, Detre J, Gee JC.
  The optimal template effect in hippocampus studies of diseased populations.
  Neuroimage. 2010 Feb 1;49(3):2457-66. doi: 10.1016/j.neuroimage.2009.09.062.

  http://www.ncbi.nlm.nih.gov/pubmed/18995188
  Avants B, Duda JT, Kim J, Zhang H, Pluta J, Gee JC, Whyte J.
  Multivariate analysis of structural and diffusion imaging in traumatic
  brain injury.
  Acad Radiol. 2008 Nov;15(11):1360-75. doi: 10.1016/j.acra.2008.07.007.

Atropos (n-tissue multivariate segmentation)
  http://www.ncbi.nlm.nih.gov/pubmed/21373993
  Avants BB, Tustison NJ, Wu J, Cook PA, Gee JC.
  An open source multivariate framework for n-tissue segmentation with
  evaluation on public data.
  Neuroinformatics. 2011 Dec;9(4):381-400. doi: 10.1007/s12021-011-9109-y.

N4 bias correction
  http://www.ncbi.nlm.nih.gov/pubmed/20378467
  Tustison NJ, Avants BB, Cook PA, Zheng Y, Egan A, Yushkevich PA, Gee JC.
  N4ITK: improved N3 bias correction.
  IEEE Trans Med Imaging. 2010 Jun;29(6):1310-20. doi: 10.1109/TMI.2010.2046908.

DiReCT aka KellyKapowski/KellySlater
  http://www.ncbi.nlm.nih.gov/pubmed/19150502
  Das SR, Avants BB, Grossman M, Gee JC.
  Registration based cortical thickness measurement.
  Neuroimage. 2009 Apr 15;45(3):867-79. doi: 10.1016/j.neuroimage.2008.12.016.

SCCAN
  http://www.ncbi.nlm.nih.gov/pubmed/20083207
  Avants BB, Cook PA, Ungar L, Gee JC, Grossman M.
  Dementia induces correlated reductions in white matter integrity and cortical
  thickness: a multivariate neuroimaging study with sparse canonical correlation
  analysis.
  Neuroimage. 2010 Apr 15;50(3):1004-16. doi: 10.1016/j.neuroimage.2010.01.041.

SurfaceCurvature
  http://www.ncbi.nlm.nih.gov/pubmed/15344450
  Avants B, Gee J.
  The shape operator for differential analysis of images.
  Inf Process Med Imaging. 2003 Jul;18:101-13.

Topological well-composedness
  http://www.ncbi.nlm.nih.gov/pubmed/21118779
  Tustison NJ, Avants BB, Siqueira M, Gee JC.
  Topological well-composedness and glamorous glue: a digital gluing algorithm
  for topologically constrained front propagation.
  IEEE Trans Image Process. 2011 Jun;20(6):1756-61. doi: 10.1109/TIP.2010.2095021.

ANTs-related Studies
  http://www.ncbi.nlm.nih.gov/pubmed/15948659
  Avants BB, Schoenemann PT, Gee JC.
  Lagrangian frame diffeomorphic image registration: Morphometric comparison of
  human and chimpanzee cortex.
  Med Image Anal. 2006 Jun;10(3):397-412. Epub 2005 Jun 3.

  http://www.ncbi.nlm.nih.gov/pubmed/19195496
  Klein A, Andersson J, Ardekani BA, Ashburner J, Avants B, Chiang MC,
  Christensen GE, Collins DL, Gee J, Hellier P, Song JH, Jenkinson M, Lepage C,
  Rueckert D, Thompson P, Vercauteren T, Woods RP, Mann JJ, Parsey RV.
  Evaluation of 14 nonlinear deformation algorithms applied to human brain MRI
  registration.
  Neuroimage. 2009 Jul 1;46(3):786-802. doi: 10.1016/j.neuroimage.2008.12.037.

  http://www.ncbi.nlm.nih.gov/pubmed/20123029
  Klein A, Ghosh SS, Avants B, Yeo BT, Fischl B, Ardekani B, Gee JC, Mann JJ,
  Parsey RV.
  Evaluation of volume-based and surface-based brain image registration methods.
  Neuroimage. 2010 May 15;51(1):214-20. doi: 10.1016/j.neuroimage.2010.01.091.

  http://www.ncbi.nlm.nih.gov/pubmed/23151955
  Tustison NJ, Avants BB, Cook PA, Kim J, Whyte J, Gee JC, Stone JR.
  Logical circularity in voxel-based analysis: Normalization strategy may induce
  statistical bias.
  Hum Brain Mapp. 2012 Nov 14. doi: 10.1002/hbm.22211.

  http://www.ncbi.nlm.nih.gov/pubmed/17999940
  Kim J, Avants B, Patel S, Whyte J, Coslett BH, Pluta J, Detre JA, Gee JC.
  Structural consequences of diffuse traumatic brain injury: a large deformation
  tensor-based morphometry study.
  Neuroimage. 2008 Feb 1;39(3):1014-26. Epub 2007 Oct 13.


# gource visualization
gource ./ -s 0.05 --stop-at-end --output-ppm-stream ants.ppm
ffmpeg -y  -f image2pipe -vcodec ppm -i ants.ppm -vcodec mpeg4 -preset slow -crf 2 -b:v 4M ./ants_gource.mp4
