/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>
#include <vnl/vnl_inverse.h>
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "antsAllocImage.h"
#include "antsSCCANObject.h"
#include "itkAlternatingValueDifferenceImageFilter.h"
#include "itkAlternatingValueSimpleSubtractionImageFilter.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkArray.h"
#include "itkAverageOverDimensionImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkBlackTopHatImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCompositeValleyFunction.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConvolutionImageFilter.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkCyclicShiftImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkFastMarchingExtensionImageFilterBase.h"
#include "itkFastMarchingExtensionImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHistogram.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageGaussianModelEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkListSample.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkMRIBiasFieldCorrectionFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkMultivariateLegendrePolynomial.h"
#include "itkNeighborhood.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodFirstOrderStatisticsImageFilter.h"
#include "itkNormalVariateGenerator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkSampleToHistogramFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkSize.h"
#include "itkSliceTimingCorrectionImageFilter.h"
#include "itkSphereSpatialFunction.h"
#include "itkSplitAlternatingTimeSeriesImageFilter.h"
#include "itkSTAPLEImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSumProjectionImageFilter.h"
#include "itkTDistribution.h"
#include "itkTileImageFilter.h"
#include "itkTimeProbe.h"
#include "itkTranslationTransform.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkWhiteTopHatImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkTransformFactory.h"
#include "itkSurfaceImageCurvature.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"

#include <fstream>
#include <iostream>
#include <map> // Here I'm using a map but you could choose even other containers
#include <sstream>
#include <string>

#include "ReadWriteData.h"
#include "TensorFunctions.h"
#include "antsMatrixUtilities.h"
#include "antsFastMarchingImageFilter.h"
#include "itkFastMarchingImageFilterBase.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"

namespace ants
{

// External functions in separate files for more module compilation
// These functions are defined in independant compilation units in
// ImageMathHelper2D.cpp, ImageMathHelper3D.cpp, and ImageMathHelper4D.cpp
extern int
ImageMathHelper2D(int argc, char ** argv);
extern int
ImageMathHelper3D(int argc, char ** argv);
extern int
ImageMathHelper4D(int argc, char ** argv);

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ImageMath(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ImageMath");

  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);
  // antscout->set_stream( out_stream );

  if (argc < 5)
  {
    std::cout << "\nUsage: " << argv[0]
              << " ImageDimension <OutputImage.ext> [operations and inputs] <Image1.ext> <Image2.ext>" << std::endl;

    std::cout << "\nUsage Information " << std::endl;
    std::cout << " ImageDimension: 2 or 3 (for 2 or 3 dimensional operations)." << std::endl;
    std::cout << " ImageDimension: 4 (for operations on 4D file, e.g. time-series data)." << std::endl;
    std::cout << " Operator: See list of valid operators below." << std::endl;
    std::cout << " The last two arguments can be an image or float value " << std::endl;
    std::cout << " NB: Some options output text files" << std::endl;

    std::cout << "\nMathematical Operations:" << std::endl;
    std::cout << "  m            : Multiply ---  use vm for vector multiply " << std::endl;
    std::cout << "  +             : Add ---  use v+ for vector add " << std::endl;
    std::cout << "  -             : Subtract ---  use v- for vector subtract " << std::endl;
    std::cout << "  /             : Divide" << std::endl;
    std::cout << "  ^            : Power" << std::endl;
    std::cout << "  max            : voxelwise max" << std::endl;
    std::cout << "  exp            : Take exponent exp(imagevalue*value)" << std::endl;
    std::cout << "  addtozero        : add image-b to image-a only over points where image-a has zero values"
              << std::endl;
    std::cout << "  overadd        : replace image-a pixel with image-b pixel if image-b pixel is non-zero"
              << std::endl;
    std::cout << "  abs            : absolute value " << std::endl;
    std::cout << "  total            : Sums up values in an image or in image1*image2 (img2 is the probability mask)"
              << std::endl;
    std::cout << "  mean            :  Average of values in an image or in image1*image2 (img2 is the probability mask)"
              << std::endl;
    std::cout << "  vtotal            : Sums up volumetrically weighted values in an image or in image1*image2 (img2 "
                 "is the probability mask)"
              << std::endl;
    std::cout << "  Decision        : Computes result=1./(1.+exp(-1.0*( pix1-0.25)/pix2))" << std::endl;
    std::cout << "  Neg            : Produce image negative" << std::endl;

    std::cout << "\nSpatial Filtering:" << std::endl;
    std::cout << "  Project Image1.ext axis-a which-projection   : Project an image along axis a, "
                 "which-projection=0(sum, 1=max, 2=min)"
              << std::endl;
    std::cout << "  G Image1.ext s    : Smooth with Gaussian of sigma = s" << std::endl;
    std::cout << "  MD Image1.ext s    : Morphological Dilation with radius s" << std::endl;
    std::cout << "  ME Image1.ext s    : Morphological Erosion with radius s" << std::endl;
    std::cout << "  MO Image1.ext s    : Morphological Opening with radius s" << std::endl;
    std::cout << "  MC Image1.ext s    : Morphological Closing with radius s" << std::endl;
    std::cout << "  GD Image1.ext s    : Grayscale Dilation with radius s" << std::endl;
    std::cout << "  GE Image1.ext s    : Grayscale Erosion with radius s" << std::endl;
    std::cout << "  GO Image1.ext s    : Grayscale Opening with radius s" << std::endl;
    std::cout << "  GC Image1.ext s    : Grayscale Closing with radius s" << std::endl;
    std::cout << "  Extract contours: extract contours from a label image" << std::endl;
    std::cout << "    Usage: ExtractContours inputImage [doFullyConnected=1]" << std::endl;
    std::cout
      << "  BlobDetector Image1.ext NumberOfBlobsToExtract  Optional-Input-Image2 Blob-2-out.nii.gz N-Blobs-To-Match  "
         ":  blob detection by searching for local extrema of the Laplacian of the Gassian (LoG) "
      << std::endl;
    std::cout << "    Example matching 6 best blobs from 2 images: " << std::endl;
    std::cout << "    ImageMath 2 blob.nii.gz BlobDetector image1.nii.gz 1000  image2.nii.gz blob2.nii.gz 6 "
              << std::endl;
    std::cout << "  MatchBlobs Image1.ext Image1LM.ext Image2.ext" << std::endl;
    std::cout << std::endl;

    std::cout << "\nTransform Image: " << std::endl;
    std::cout << "Translate InImage.ext x [ y z ] " << std::endl;

    std::cout << "\nTime Series Operations:" << std::endl;
    std::cout << " CompCorrAuto : Outputs a csv file containing global signal vector and N comp-corr eigenvectors "
                 "determined from PCA of the high-variance voxels.  Also outputs a comp-corr + global signal corrected "
                 "4D image as well as a 3D image measuring the time series variance.  Requires a label image with "
                 "label 1 identifying voxels in the brain."
              << std::endl;
    std::cout << "   ImageMath 4 ${out}compcorr.nii.gz ThreeTissueConfounds ${out}.nii.gz  ${out}seg.nii.gz 1 3  "
              << " : Outputs average global, CSF and WM signals.  Requires a label image with 3 labels , csf, gm , wm ."
              << std::endl;
    std::cout << "    Usage        : ThreeTissueConfounds 4D_TimeSeries.nii.gz LabeLimage.nii.gz  csf-label wm-label "
              << std::endl;
    std::cout
      << " TimeSeriesSubset : Outputs n 3D image sub-volumes extracted uniformly from the input time-series 4D image."
      << std::endl;
    std::cout << "    Usage        : TimeSeriesSubset 4D_TimeSeries.nii.gz n " << std::endl;
    std::cout << " TimeSeriesDisassemble : Outputs n 3D image volumes for each time-point in time-series 4D image."
              << std::endl;
    std::cout << "    Usage        : TimeSeriesDisassemble 4D_TimeSeries.nii.gz " << std::endl << std::endl;
    std::cout << " TimeSeriesAssemble : Outputs a 4D time-series image from a list of 3D volumes." << std::endl;
    std::cout << "    Usage        : TimeSeriesAssemble time_spacing time_origin *images.nii.gz " << std::endl;
    std::cout << " TimeSeriesToMatrix : Converts a 4D image + mask to matrix (stored as csv file) where rows are time "
                 "and columns are space ."
              << std::endl;
    std::cout << "    Usage        : TimeSeriesToMatrix 4D_TimeSeries.nii.gz mask " << std::endl;
    std::cout << " TimeSeriesSimpleSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes."
              << std::endl;
    std::cout << "    Usage        : TimeSeriesSimpleSubtraction image.nii.gz " << std::endl;
    std::cout << " TimeSeriesSurroundSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes."
              << std::endl;
    std::cout << "    Usage        : TimeSeriesSurroundSubtraction image.nii.gz " << std::endl;
    std::cout << " TimeSeriesSincSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes." << std::endl;
    std::cout << "    Usage        : TimeSeriesSincSubtraction image.nii.gz " << std::endl;
    std::cout << " SplitAlternatingTimeSeries : Outputs 2 3D time series" << std::endl;
    std::cout << "    Usage        : SplitAlternatingTimeSeries image.nii.gz " << std::endl;

    std::cout
      << " ComputeTimeSeriesLeverage : Outputs a csv file that identifies the raw leverage and normalized leverage for "
         "each time point in the 4D image.  leverage, here, is the difference of the time-point image from the average "
         "of the n images.  the normalized leverage is =  average( sum_k abs(Leverage(t)-Leverage(k)) )/Leverage(t). "
      << std::endl;
    std::cout << "    Usage        : ComputeTimeSeriesLeverage 4D_TimeSeries.nii.gz k_neighbors " << std::endl;

    std::cout << " SliceTimingCorrection : Outputs a slice-timing corrected 4D time series" << std::endl;
    std::cout << "    Usage        : SliceTimingCorrection image.nii.gz sliceTiming [sinc / bspline] [sincRadius=4 / "
                 "bsplineOrder=3]"
              << std::endl;

    std::cout << " PASL : computes the PASL model of CBF  " << std::endl
              << "f =  \frac{      lambda DeltaM        } " << std::endl
              << " {     2 alpha M_0 TI_1 exp( - TI_2 / T_{1a} )  } " << std::endl;
    std::cout << "    Usage        : PASL 3D/4D_TimeSeries.nii.gz BoolFirstImageIsControl M0Image parameter_list.txt "
              << std::endl;

    std::cout << " pCASL : computes the pCASL model of CBF  " << std::endl
              << " f =  \frac{      lambda DeltaM R_{1a}        }  " << std::endl
              << "  {     2 alpha M_0 [ exp( - w R_{1a} ) - exp( -w ( \tau + w ) R_{1a}) ]     } " << std::endl;
    std::cout << "    Usage        : pCASL 3D/4D_TimeSeries.nii.gz parameter_list.txt " << std::endl;
    std::cout << " PASLQuantifyCBF : Outputs a 3D CBF image in ml/100g/min from a magnetization ratio image"
              << std::endl;
    std::cout << "    Usage        : PASLQuantifyCBF mag_ratio.nii.gz [TI1=700] [TI2=1900] [T1blood=1664] [Lambda=0.9] "
                 "[Alpha=0.95] [SliceDelay-45] "
              << std::endl;

    std::cout << "\nTensor Operations:" << std::endl;
    std::cout << "  4DTensorTo3DTensor    : Outputs a 3D_DT_Image with the same information. " << std::endl;
    std::cout << "    Usage        : 4DTensorTo3DTensor 4D_DTImage.ext" << std::endl;
    std::cout << "  ComponentTo3DTensor    : Outputs a 3D_DT_Image with the same information as component images. "
              << std::endl;
    std::cout << "    Usage        : ComponentTo3DTensor component_image_prefix[xx,xy,xz,yy,yz,zz] extension"
              << std::endl;
    std::cout << "  ExtractComponentFrom3DTensor    : Outputs a component images. " << std::endl;
    std::cout << "    Usage        : ExtractComponentFrom3DTensor dtImage.ext which={xx,xy,xz,yy,yz,zz}" << std::endl;
    std::cout << "  ExtractVectorComponent: Produces the WhichVec component of the vector " << std::endl;
    std::cout << "    Usage        : ExtractVectorComponent VecImage WhichVec" << std::endl;
    std::cout << "  TensorColor        : Produces RGB values identifying principal directions " << std::endl;
    std::cout << "    Usage        : TensorColor DTImage.ext" << std::endl;
    std::cout << "  TensorFA        : " << std::endl;
    std::cout << "    Usage        : TensorFA DTImage.ext" << std::endl;
    std::cout << "  TensorFADenominator    : " << std::endl;
    std::cout << "    Usage        : TensorFADenominator DTImage.ext" << std::endl;
    std::cout << "  TensorFANumerator    : " << std::endl;
    std::cout << "    Usage        : TensorFANumerator DTImage.ext" << std::endl;
    std::cout << "  TensorIOTest    : Will write the DT image back out ... tests I/O processes for consistency. "
              << std::endl;
    std::cout << "    Usage        : TensorIOTest DTImage.ext" << std::endl;
    std::cout << "  TensorMeanDiffusion      : Mean of the eigenvalues" << std::endl;
    std::cout << "    Usage        : TensorMeanDiffusion DTImage.ext" << std::endl;
    std::cout << "  TensorRadialDiffusion    : Mean of the two smallest eigenvalues" << std::endl;
    std::cout << "    Usage        : TensorRadialDiffusion DTImage.ext" << std::endl;
    std::cout << "  TensorAxialDiffusion     : Largest eigenvalue, equivalent to TensorEigenvalue DTImage.ext 2"
              << std::endl;
    std::cout << "    Usage        : TensorAxialDiffusion DTImage.ext" << std::endl;
    std::cout << "  TensorEigenvalue         : Gets a single eigenvalue 0-2, where 0 = smallest, 2 = largest"
              << std::endl;
    std::cout << "    Usage        : TensorEigenvalue DTImage.ext WhichInd" << std::endl;
    std::cout << "  TensorToVector    : Produces vector field identifying one of the principal directions, 2 = largest "
                 "eigenvalue"
              << std::endl;
    std::cout << "    Usage        : TensorToVector DTImage.ext WhichVec" << std::endl;
    std::cout << "  TensorToVectorComponent: 0 => 2 produces component of the principal vector field (largest "
                 "eigenvalue). 3 = 8 => gets values from the tensor "
              << std::endl;
    std::cout << "    Usage        : TensorToVectorComponent DTImage.ext WhichVec" << std::endl;
    std::cout << "  TensorMask     : Mask a tensor image, sets background tensors to zero or to isotropic tensors with "
                 "specified mean diffusivity "
              << std::endl;
    std::cout << "    Usage        : TensorMask DTImage.ext mask.ext [ backgroundMD = 0 ] " << std::endl;
    std::cout << "  FuseNImagesIntoNDVectorField     : Create ND field from N input scalar images" << std::endl;
    std::cout << "    Usage        : FuseNImagesIntoNDVectorField imagex imagey imagez" << std::endl;

    std::cout << "\nLabel Fusion:" << std::endl;
    std::cout << "  MajorityVoting : Select label with most votes from candidates" << std::endl;
    std::cout << "    Usage: MajorityVoting LabelImage1.nii.gz .. LabelImageN.nii.gz" << std::endl;
    std::cout << "  CorrelationVoting : Select label with local correlation weights" << std::endl;
    std::cout << "    Usage: CorrelationVoting Template.ext IntenistyImages* LabelImages* {Optional-Radius=5}"
              << std::endl;
    std::cout << "  STAPLE : Select label using STAPLE method" << std::endl;
    std::cout << "    Usage: STAPLE confidence-weighting LabelImages*" << std::endl;
    std::cout << "    Note:  Gives probabilistic output (float)" << std::endl;
    std::cout << "  MostLikely : Select label from from maximum probabilistic segmentations" << std::endl;
    std::cout << "    Usage: MostLikely probabilityThreshold ProbabilityImages*" << std::endl;
    std::cout << "  AverageLabels : Select label using STAPLE method" << std::endl;
    std::cout << "    Usage: AverageLabels LabelImages*" << std::endl;
    std::cout << "    Note:  Gives probabilistic output (float)" << std::endl;

    std::cout << "\nImage Metrics & Info:" << std::endl;
    std::cout << "  PearsonCorrelation: r-value from intesities of two images" << std::endl;
    std::cout << "    Usage: PearsonCorrelation image1.ext image2.ext {Optional-mask.ext}" << std::endl;
    std::cout << "  NeighborhoodCorrelation: local correlations" << std::endl;
    std::cout << "    Usage: NeighborhoodCorrelation image1.ext image2.ext {Optional-radius=5} {Optional-image-mask}"
              << std::endl;
    std::cout << "  NormalizedCorrelation: r-value from intesities of two images" << std::endl;
    std::cout << "    Usage: NormalizedCorrelation image1.ext image2.ext {Optional-image-mask}" << std::endl;
    std::cout << "  Demons: " << std::endl;
    std::cout << "    Usage: Demons image1.ext image2.ext" << std::endl;
    std::cout << "  Mattes: mutual information" << std::endl;
    std::cout << "    Usage: Mattes image1.ext image2.ext {Optional-number-bins=32} {Optional-image-mask}" << std::endl;

    std::cout << "\nUnclassified Operators:" << std::endl;

    std::cout << "  ReflectionMatrix : Create a reflection matrix about an axis" << std::endl;
    std::cout << " out.mat ReflectionMatrix image_in axis " << std::endl << std::endl;

    std::cout << "  MakeAffineTransform : Create an itk affine transform matrix " << std::endl;

    std::cout << "  ClosestSimplifiedHeaderMatrix : does what it says ... image-in, image-out" << std::endl;

    std::cout << "  Byte            : Convert to Byte image in [0,255]" << std::endl;

    std::cout
      << "\n  CompareHeadersAndImages: Tries to find and fix header errors. Outputs a repaired image with new header. "
      << std::endl;
    std::cout << "                Never use this if you trust your header information. " << std::endl;
    std::cout << "      Usage        : CompareHeadersAndImages Image1 Image2" << std::endl;

    std::cout << "\n  ConvertImageSetToMatrix: Each row/column contains image content extracted from mask applied to "
                 "images in *img.nii "
              << std::endl;
    std::cout << "      Usage        : ConvertImageSetToMatrix rowcoloption Mask.nii *images.nii" << std::endl;
    std::cout << " ConvertImageSetToMatrix output can be an image type or csv file type." << std::endl;

    std::cout << "\n  RandomlySampleImageSetToCSV: N random samples are selected from each image in a list "
              << std::endl;
    std::cout << "      Usage        : RandomlySampleImageSetToCSV N_samples *images.nii" << std::endl;
    std::cout << " RandomlySampleImageSetToCSV outputs a csv file type." << std::endl;

    std::cout << "\n  FrobeniusNormOfMatrixDifference: take the difference between two itk-transform matrices and then "
                 "compute the frobenius norm"
              << std::endl;
    std::cout << "      Usage        : FrobeniusNormOfMatrixDifference mat1 mat2 " << std::endl;
    std::cout << "\n  ConvertImageSetToEigenvectors: Each row/column contains image content extracted from mask "
                 "applied to images in *img.nii "
              << std::endl;
    std::cout << "      Usage        : ConvertImageSetToEigenvectors N_Evecs Mask.nii *images.nii" << std::endl;
    std::cout << " ConvertImageSetToEigenvectors output will be a csv file for each label value > 0 in the mask."
              << std::endl;

    std::cout << "\n  ConvertImageToFile    : Writes voxel values to a file  " << std::endl;
    std::cout << "      Usage        : ConvertImageToFile imagevalues.nii {Optional-ImageMask.nii}" << std::endl;

    std::cout
      << "\n  ConvertLandmarkFile    : Converts landmark file between formats. See ANTS.pdf for description of formats."
      << std::endl;
    std::cout << "      Usage        : ConvertLandmarkFile InFile.txt" << std::endl;
    std::cout << "      Example 1        : ImageMath 3  outfile.vtk  ConvertLandmarkFile  infile.txt" << std::endl;

    std::cout << "\n  ConvertToGaussian    : " << std::endl;
    std::cout << "      Usage        : ConvertToGaussian  TValueImage  sigma-float" << std::endl;

    std::cout << "\n  ConvertVectorToImage    : The vector contains image content extracted from a mask. Here the "
                 "vector is returned to its spatial origins as image content "
              << std::endl;
    std::cout << "      Usage        : ConvertVectorToImage Mask.nii vector.nii" << std::endl;

    std::cout << "\n  CorrelationUpdate    : In voxels, compute update that makes Image2 more like Image1."
              << std::endl;
    std::cout << "      Usage        : CorrelationUpdate Image1.ext Image2.ext RegionRadius" << std::endl;

    std::cout << "\n  CountVoxelDifference    : The where function from IDL " << std::endl;
    std::cout << "      Usage        : CountVoxelDifference Image1 Image2 Mask" << std::endl;

    std::cout << "\n  CorruptImage        : " << std::endl;
    std::cout << "      Usage        : CorruptImage Image NoiseLevel Smoothing" << std::endl;

    std::cout << "\n  D             : Danielson Distance Transform" << std::endl;

    std::cout << "\n  MaurerDistance : Maurer distance transform (much faster than Danielson)" << std::endl;
    std::cout << "      Usage        : MaurerDistance inputImage {foreground=1}" << std::endl;

    std::cout << "\n  DiceAndMinDistSum    : Outputs DiceAndMinDistSum and Dice Overlap to text log file + optional "
                 "distance image"
              << std::endl;
    std::cout << "      Usage        : DiceAndMinDistSum LabelImage1.ext LabelImage2.ext OptionalDistImage"
              << std::endl;

    std::cout << "\n  EnumerateLabelInterfaces: " << std::endl;
    std::cout << "      Usage        : EnumerateLabelInterfaces ImageIn ColoredImageOutname NeighborFractionToIgnore"
              << std::endl;

    std::cout << "\n  ClusterThresholdVariate        :  for sparse estimation " << std::endl;
    std::cout << "      Usage        : ClusterThresholdVariate image mask  MinClusterSize" << std::endl;

    std::cout << "\n  ExtractSlice        : Extracts slice number from last dimension of volume (2,3,4) dimensions "
              << std::endl;
    std::cout << "      Usage        : ExtractSlice volume.nii.gz slicetoextract" << std::endl;

    std::cout << "\n  FastMarchingSegmentation: final output is the propagated label image. Optional stopping value: "
                 "higher values allow more distant propagation "
              << std::endl;
    std::cout << "      Usage        : FastMarchingSegmentation speed/binaryimagemask.ext initiallabelimage.ext "
                 "Optional-Stopping-Value"
              << std::endl;

    std::cout << "\n  FillHoles        : Parameter = ratio of edge at object to edge at background;  --  " << std::endl;
    std::cout << "                Parameter = 1 is a definite hole bounded by object only, 0.99 is close" << std::endl;
    std::cout << "                Default of parameter > 1 will fill all holes" << std::endl;
    std::cout << "      Usage        : FillHoles Image.ext parameter" << std::endl;

    std::cout << "\n  InPaint        : very simple inpainting --- assumes zero values should be inpainted  "
              << std::endl;
    std::cout << "      Usage        : InPaint #iterations" << std::endl;

    std::cout << "\n  PeronaMalik       : anisotropic diffusion w/varying conductance param (0.25 in example below)"
              << std::endl;
    std::cout << "      Usage        : PeronaMalik image #iterations conductance " << std::endl;

    std::cout << "\n  Convolve       : convolve input image with kernel image" << std::endl;
    std::cout << "      Usage        : Convolve inputImage kernelImage {normalize=1} " << std::endl;

    std::cout << "  Finite            : replace non-finite values with finite-value (default = 0)" << std::endl;
    std::cout << "      Usage        : Finite Image.exdt {replace-value=0}" << std::endl;

    std::cout << "\n  LabelSurfaceArea        : " << std::endl;
    std::cout << "      Usage        : LabelSurfaceArea ImageIn {MaxRad-Default=1}" << std::endl;

    std::cout << "\n  FlattenImage        : Replaces values greater than %ofMax*Max to the value %ofMax*Max "
              << std::endl;
    std::cout << "      Usage        : FlattenImage Image %ofMax" << std::endl;

    std::cout << "\n  GetLargestComponent    : Get the largest object in an image" << std::endl;
    std::cout << "      Usage        : GetLargestComponent InputImage {MinObjectSize}" << std::endl;

    std::cout << "\n  Grad            : Gradient magnitude with sigma s (if normalize, then output in range [0, 1])"
              << std::endl;
    std::cout << "      Usage        : Grad Image.ext s normalize?" << std::endl;

    std::cout << "\n  HistogramMatch    : " << std::endl;
    std::cout << "      Usage        : HistogramMatch SourceImage ReferenceImage {NumberBins-Default=255} "
                 "{NumberPoints-Default=64} {useThresholdAtMeanIntensity=false}"
              << std::endl;

    std::cout << "\n  RescaleImage    : " << std::endl;
    std::cout << "      Usage        : RescaleImage InputImage min max" << std::endl;
    std::cout << "\n  WindowImage    : " << std::endl;
    std::cout << "      Usage        : WindowImage InputImage windowMinimum windowMaximum outputMinimum outputMaximum"
              << std::endl;
    std::cout << "\n  NeighborhoodStats    : " << std::endl;
    std::cout
      << "      Usage        : NeighborhoodStats inputImage whichStat radius"
         "             whichStat:  1 = min, 2 = max, 3 = variance, 4 = sigma, 5 = skewness, 6 = kurtosis, 7 = entropy"
      << std::endl;

    std::cout << "\n  InvId            : computes the inverse-consistency of two deformations and write the inverse "
                 "consistency error image "
              << std::endl;
    std::cout << "      Usage        : InvId VectorFieldName VectorFieldName" << std::endl;

    std::cout << "\n  ReplicateDisplacement            : replicate a ND displacement to a ND+1 image" << std::endl;
    std::cout << "      Usage        : ReplicateDisplacement VectorFieldName TimeDims TimeSpacing TimeOrigin"
              << std::endl;

    std::cout << "\n  ReplicateImage            : replicate a ND image to a ND+1 image" << std::endl;
    std::cout << "      Usage        : ReplicateImage ImageName TimeDims TimeSpacing TimeOrigin" << std::endl;

    std::cout << "\n  ShiftImageSlicesInTime            : shift image slices by one " << std::endl;
    std::cout
      << "      Usage        : ShiftImageSlicesInTime ImageName shift-amount-default-1 shift-dim-default-last-dim"
      << std::endl;

    std::cout << "\n  LabelStats        : Compute volumes / masses of objects in a label image. Writes to text file"
              << std::endl;
    std::cout << "      Usage        : LabelStats labelimage.ext valueimage.nii" << std::endl;

    std::cout << "\n  Laplacian        : Laplacian computed with sigma s (if normalize, then output in range [0, 1])"
              << std::endl;
    std::cout << "      Usage        : Laplacian Image.ext s normalize?" << std::endl;

    std::cout << "\n  Canny        : Canny edge detector" << std::endl;
    std::cout << "      Usage        : Canny Image.ext sigma lowerThresh upperThresh" << std::endl;

    std::cout << "\n  Lipschitz        : Computes the Lipschitz norm of a vector field " << std::endl;
    std::cout << "      Usage        : Lipschitz VectorFieldName" << std::endl;

    std::cout << "\n  MakeImage        : " << std::endl;
    std::cout << "      Usage        : MakeImage SizeX  SizeY {SizeZ};" << std::endl;

    std::cout
      << "\n  MTR        : Computes the magnetization transfer ratio ( (M0-M1)/M0 ) and truncates values to [0,1]"
      << std::endl;
    std::cout << "      Usage        : MTR M0Image M1Image [MaskImage];" << std::endl;

    std::cout << "\n  Normalize        : Normalize to [0,1]. Option instead divides by average value.  If opt is a "
                 "mask image, then we normalize by mean intensity in the mask ROI."
              << std::endl;
    std::cout << "      Usage        : Normalize Image.ext opt" << std::endl;

    std::cout << "\n  PadImage       : If Pad-Number is negative, de-Padding occurs" << std::endl;
    std::cout << "      Usage        : PadImage ImageIn PaddingSize [PaddingVoxelValue=0]" << std::endl;

    std::cout << "\n  SigmoidImage   : " << std::endl;
    std::cout << "      Usage        : SigmoidImage ImageIn [alpha=1.0] [beta=0.0]" << std::endl;

    std::cout << "\n  Sharpen        : Apply a Laplacian sharpening filter" << std::endl;
    std::cout << "      Usage        : Sharpen ImageIn [useImageSpacing=(1)/0]" << std::endl;

    std::cout << "\n  UnsharpMask     Apply an Unsharp Mask filter" << std::endl;
    std::cout
      << "      Usage        : UnsharpMask ImageIn [amount=0.5] [radius=1] [threshold=0] [radius in spacing unit (0)/1]"
      << std::endl;

    std::cout << "\n  CoordinateComponentImages   : " << std::endl;
    std::cout << "      Usage        : CoordinateComponentImages domainImage" << std::endl;

    std::cout << "\n  CenterImage2inImage1        : " << std::endl;
    std::cout << "      Usage       : ReferenceImageSpace ImageToCenter " << std::endl;

    std::cout << "\n  PH            : Print Header" << std::endl;

    std::cout << "\n  PoissonDiffusion        : Solves Poisson's equation in a designated region using non-zero sources"
              << std::endl;
    std::cout << "      Usage        : PoissonDiffusion inputImage labelImage [sigma=1.0] [regionLabel=1] "
                 "[numberOfIterations=500] [convergenceThreshold=1e-10]"
              << std::endl;

    std::cout << "\n  PropagateLabelsThroughMask: Final output is the propagated label image. Optional stopping value: "
                 "higher values allow more distant propagation"
              << std::endl;
    std::cout << "      Usage        : PropagateLabelsThroughMask speed/binaryimagemask.nii.gz "
                 "initiallabelimage.nii.gz Optional-Stopping-Value  0/1/2"
              << std::endl;
    std::cout << "      0/1/2  =>  0, no topology constraint, 1 - strict topology constraint, 2 - no handles "
              << std::endl;

    std::cout << "\n  PValueImage        : " << std::endl;
    std::cout << "      Usage        : PValueImage TValueImage dof" << std::endl;

    std::cout << "\n  RemoveLabelInterfaces: " << std::endl;
    std::cout << "      Usage        : RemoveLabelInterfaces ImageIn" << std::endl;

    std::cout << "\n  ReplaceVoxelValue: replace voxels in the range [a,b] in the input image with c" << std::endl;
    std::cout << "      Usage        : ReplaceVoxelValue inputImage a b c" << std::endl;
    std::cout << "\n  ROIStatistics        : computes anatomical locations, cluster size and mass of a stat image "
                 "which should be in the same physical space (but not nec same resolution) as the label image."
              << std::endl;
    std::cout << "      Usage        : ROIStatistics LabelNames.txt labelimage.ext valueimage.nii" << std::endl;

    std::cout << "\n  SetOrGetPixel    : " << std::endl;
    std::cout << "      Usage        : SetOrGetPixel ImageIn Get/Set-Value IndexX IndexY {IndexZ}" << std::endl;
    std::cout
      << "      Example 1        : ImageMath 2 outimage.nii SetOrGetPixel Image Get 24 34; Gets the value at 24, 34"
      << std::endl;
    std::cout << "      Example 2        : ImageMath 2 outimage.nii SetOrGetPixel Image 1.e9 24 34; This sets 1.e9 as "
                 "the value at 23 34"
              << std::endl;
    std::cout << "                You can also pass a boolean at the end to force the physical space to be used"
              << std::endl;

    std::cout << "\n  SetTimeSpacing            : sets spacing for last dimension" << std::endl;
    std::cout << "      Usage        : SetTimeSpacing Image.ext tspacing" << std::endl;

    std::cout << "\n  SetTimeSpacingWarp            : sets spacing for last dimension" << std::endl;
    std::cout << "      Usage        : SetTimeSpacingWarp Warp.ext tspacing" << std::endl;

    std::cout << "\n  stack            : Will put 2 images in the same volume" << std::endl;
    std::cout << "      Usage        : Stack Image1.ext Image2.ext" << std::endl;

    std::cout << "\n  ThresholdAtMean    : See the code" << std::endl;
    std::cout << "      Usage        : ThresholdAtMean Image %ofMean" << std::endl;

    std::cout << "\n  TileImages    : " << std::endl;
    std::cout << "      Usage        : TileImages NumColumns ImageList*" << std::endl;

    std::cout << "\n  TriPlanarView    : " << std::endl;
    std::cout << "      Usage        : TriPlanarView  ImageIn.nii.gz PercentageToClampLowIntensity "
                 "PercentageToClampHiIntensity x-slice y-slice z-slice"
              << std::endl;

    std::cout << "\n  TruncateImageIntensity: " << std::endl;
    std::cout << "      Usage        : TruncateImageIntensity InputImage.ext {lowerQuantile=0.05} {upperQuantile=0.95} "
                 "{numberOfBins=65} {binary-maskImage}"
              << std::endl;

    std::cout << "\n  Where            : The where function from IDL" << std::endl;
    std::cout << "      Usage        : Where Image ValueToLookFor maskImage-option tolerance" << std::endl;
    std::cout << "\n  KinematicTensor    : Evaluation kinematic tensor from a displacement field" << std::endl;
    std::cout << "      Usage        : KinematicTensor displacementField whichTensor ["
              << "'d'=DeformationFieldGradient, "
              << "'l'=Lagrangian, "
              << "'e'=Eulerian, "
              << "'rc'=RightCauchyGreen, "
              << "'lc'=LeftCauchyGreen, "
              << "'rs'=RightStretch, "
              << "'ls'=LeftStretch]" << std::endl;

    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  int returnvalue = EXIT_SUCCESS;

  std::string operation = std::string(argv[3]);

  unsigned int imageDimension = std::stoi(argv[1]);

  switch (imageDimension)
  {
    case 2:
      returnvalue = ImageMathHelper2D(argc, argv);
      break;
    case 3:
      returnvalue = ImageMathHelper3D(argc, argv);
      break;
    case 4:
      returnvalue = ImageMathHelper4D(argc, argv);
      break;
    default:
      std::cout << " Dimension " << imageDimension << " is not supported " << std::endl;
      return EXIT_FAILURE;
  }

  if (returnvalue == EXIT_FAILURE)
  {
    std::cout << " Operation " << operation << " not found or not supported for dimension " << imageDimension
              << std::endl;
  }

  return returnvalue;
}

} // namespace ants
