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
#include "antsAllocImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkArray.h"
#include "itkGradientImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkCompositeValleyFunction.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHistogram.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImage.h"
#include "itkImageClassifierBase.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageGaussianModelEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkListSample.h"
#include "itkMRFImageFilter.h"
#include "itkMRIBiasFieldCorrectionFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMultivariateLegendrePolynomial.h"
#include "itkMultiStartOptimizerv4.h"
#include "itkNeighborhood.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkNormalVariateGenerator.h"
#include "itkOptimizerParameterScalesEstimator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSampleToHistogramFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSimilarity3DTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkSize.h"
#include "itkSphereSpatialFunction.h"
#include "itkSTAPLEImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTDistribution.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkTransformFactory.h"
#include "itkSurfaceImageCurvature.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkCompositeTransform.h"

#include <fstream>
#include <iostream>
#include <map> // Here I'm using a map but you could choose even other containers
#include <sstream>
#include <string>

#include "ReadWriteData.h"
#include "TensorFunctions.h"
#include "antsMatrixUtilities.h"

namespace ants
{

template <typename TComputeType, unsigned int ImageDimension>
class SimilarityTransformTraits
{
  // Don't worry about the fact that the default option is the
  // affine Transform, that one will not actually be instantiated.
public:
  using TransformType = itk::AffineTransform<TComputeType, ImageDimension>;
};


template <>
class SimilarityTransformTraits<double, 2>
{
public:
  using TransformType = itk::Similarity2DTransform<double>;
};

template <>
class SimilarityTransformTraits<float, 2>
{
public:
  using TransformType = itk::Similarity2DTransform<float>;
};

template <>
class SimilarityTransformTraits<double, 3>
{
public:
  using TransformType = itk::Similarity3DTransform<double>;
};

template <>
class SimilarityTransformTraits<float, 3>
{
public:
  using TransformType = itk::Similarity3DTransform<float>;
};

using RealType = double;

// Specializations try to rotate around tertiary and secondary axis

static void l_rotate_around_tertiatry_and_secondary_axis(vnl_vector_fixed<RealType, 3> &       evec_tert,
                                                         const vnl_vector_fixed<RealType, 3> & evec1_primary,
                                                         vnl_vector_fixed<RealType, 3> &       evec1_2ndary)
{
  evec_tert = vnl_cross_3d(evec1_primary, evec1_2ndary);
}

static void l_rotate_around_tertiatry_and_secondary_axis(vnl_vector_fixed<RealType, 2> &       evec_tert,
                                                         const vnl_vector_fixed<RealType, 2> & evec1_primary,
                                                         vnl_vector_fixed<RealType, 2> &       evec1_2ndary)
{
  evec_tert = evec1_2ndary;
  evec1_2ndary = evec1_primary;
}

static void l_rotate_around_tertiatry_and_secondary_axis(vnl_vector_fixed<RealType, 4> &,
                                                         const vnl_vector_fixed<RealType, 4> &,
                                                         vnl_vector_fixed<RealType, 4> &)
{
  return; // Do nothing in the case of 4D
}

template <unsigned int ImageDimension>
int
antsAffineInitializerImp(int argc, char * argv[])
{

  using PixelType = float;

  /** Define All Parameters Here */
  double       pi = itk::Math::pi;      // probably a vnl alternative
  RealType     searchfactor = 10;       // in degrees, passed by user
  unsigned int mibins = 32;             // for mattes MI metric
  RealType     degtorad = 0.0174532925; // to convert degrees to radians
  // NOT USED: unsigned int localoptimizeriterations = 20;    // for local search via conjgrad
  // piover4 is (+/-) for cross-section of the sphere to multi-start search in increments
  // of searchfactor ( converted from degrees to radians ).
  // the search is centered +/- from the principal axis alignment of the images.
  RealType piover4 = pi / 4; // works in preliminary practical examples in 3D, in 2D use pi.
  bool     useprincaxis = false;
  using maskimagetype = typename itk::ImageMaskSpatialObject<ImageDimension>::ImageType;
  std::string  whichMetric = std::string("MI");
  unsigned int localSearchIterations = 20;
  using TransformWriterType = itk::TransformFileWriter;
  using VectorType = itk::Vector<float, ImageDimension>;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using ImageCalculatorType = typename itk::ImageMomentsCalculator<ImageType>;
  using AffineType = itk::AffineTransform<RealType, ImageDimension>;
  using MatrixType = typename ImageCalculatorType::MatrixType;
  if (argc < 2)
  {
    return 0;
  }
  int         argct = 2;
  std::string fn1 = std::string(argv[argct]);
  argct++;
  std::string fn2 = std::string(argv[argct]);
  argct++;
  std::string outname = std::string(argv[argct]);
  argct++;
  if (argc > argct)
  {
    searchfactor = atof(argv[argct]);
    argct++;
  }
  if (argc > argct)
  {
    RealType temp = atof(argv[argct]);
    argct++;
    if (temp > 1)
    {
      temp = 1;
    }
    if (temp < 0.01)
    {
      temp = 0.01;
    }
    piover4 = pi * temp;
  }
  if (argc > argct)
  {
    useprincaxis = std::stoi(argv[argct]);
    argct++;
  }
  if (argc > argct)
  {
    // NOT USED: localoptimizeriterations = std::stoi( argv[argct] );
    argct++;
  }
  typename ImageType::Pointer     image1 = nullptr;
  typename ImageType::Pointer     image2 = nullptr;
  typename maskimagetype::Pointer mask = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str());
  ReadImage<ImageType>(image2, fn2.c_str());
  std::string maskfn = "";
  if (argc > argct)
  {
    maskfn = std::string(argv[argct]);
    argct++;
    ReadImage<maskimagetype>(mask, maskfn.c_str());
  }
  searchfactor *= degtorad; // convert degrees to radians
  VectorType ccg1;
  VectorType cpm1;
  MatrixType cpa1;
  VectorType ccg2;
  VectorType cpm2;
  MatrixType cpa2;

  typename ImageCalculatorType::Pointer calculator1 = ImageCalculatorType::New();
  typename ImageCalculatorType::Pointer calculator2 = ImageCalculatorType::New();
  calculator1->SetImage(image1);
  calculator2->SetImage(image2);
  typename ImageCalculatorType::VectorType fixed_center;
  fixed_center.Fill(0);
  typename ImageCalculatorType::VectorType moving_center;
  moving_center.Fill(0);
  try
  {
    calculator1->Compute();
    fixed_center = calculator1->GetCenterOfGravity();
    ccg1 = calculator1->GetCenterOfGravity();
    cpm1 = calculator1->GetPrincipalMoments();
    cpa1 = calculator1->GetPrincipalAxes();
    try
    {
      calculator2->Compute();
      moving_center = calculator2->GetCenterOfGravity();
      ccg2 = calculator2->GetCenterOfGravity();
      cpm2 = calculator2->GetPrincipalMoments();
      cpa2 = calculator2->GetPrincipalAxes();
    }
    catch (...)
    {
      std::cerr << " zero image2 error ";
      fixed_center.Fill(0);
    }
  }
  catch (...)
  {
    std::cerr << " zero image1 error ";
  }
  RealType bestscale = calculator2->GetTotalMass() / calculator1->GetTotalMass();
  RealType powlev = 1.0 / static_cast<RealType>(ImageDimension);
  bestscale = std::pow(bestscale, powlev);
  bestscale = 1;
  unsigned int eigind1 = 1;
  unsigned int eigind2 = 1;
  if (ImageDimension == 3)
  {
    eigind1 = 2;
  }
  vnl_vector_fixed<RealType, ImageDimension> evec1_primary{ cpa1.GetVnlMatrix().get_row(eigind1).as_vector() };
  vnl_vector_fixed<RealType, ImageDimension> evec2_primary{ cpa2.GetVnlMatrix().get_row(eigind1).as_vector() };
  vnl_vector_fixed<RealType, ImageDimension> evec1_2ndary{ cpa1.GetVnlMatrix().get_row(eigind2).as_vector() };
  vnl_vector_fixed<RealType, ImageDimension> evec2_2ndary{ cpa2.GetVnlMatrix().get_row(eigind2).as_vector() };
  /** Solve Wahba's problem --- http://en.wikipedia.org/wiki/Wahba%27s_problem */
  vnl_matrix_fixed<RealType, ImageDimension, ImageDimension> B = outer_product(evec2_primary, evec1_primary);
  if (ImageDimension == 3)
  {
    B = outer_product(evec2_2ndary, evec1_2ndary) + outer_product(evec2_primary, evec1_primary);
  }
  vnl_svd_fixed<RealType, ImageDimension, ImageDimension>    wahba(B);
  vnl_matrix_fixed<RealType, ImageDimension, ImageDimension> A_solution =
    vnl_inverse(wahba.V() * wahba.U().transpose());
  const RealType det = vnl_determinant(A_solution);
  if (det < 0)
  {
    std::cerr << " bad det " << det << " v " << vnl_determinant(wahba.V()) << " u " << vnl_determinant(wahba.U())
              << std::endl;
    vnl_matrix_fixed<RealType, ImageDimension, ImageDimension> id(A_solution);
    id.set_identity();
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      if (A_solution(i, i) < 0)
      {
        id(i, i) = -1.0;
      }
    }
    A_solution = A_solution * id.transpose();
    std::cerr << " bad det " << det << " v " << vnl_determinant(wahba.V()) << " u " << vnl_determinant(wahba.U())
              << " new " << vnl_determinant(A_solution) << std::endl;
  }

  typename AffineType::MatrixType AA_solution = typename AffineType::MatrixType(A_solution);

  typename AffineType::Pointer       affine1 = AffineType::New(); // translation to center
  typename AffineType::OffsetType    trans = affine1->GetOffset();
  itk::Point<double, ImageDimension> trans2;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    trans[i] = moving_center[i] - fixed_center[i];
    trans2[i] = fixed_center[i] * (1);
  }
  affine1->SetIdentity();
  affine1->SetOffset(trans);
  if (useprincaxis)
  {
    affine1->SetMatrix(AA_solution);
  }
  affine1->SetCenter(trans2);
  {
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetInput(affine1);
    transformWriter->SetFileName(outname.c_str());
#if ITK_VERSION_MAJOR >= 5
    transformWriter->SetUseCompression(true);
#endif
    transformWriter->Update();
  }
  if (ImageDimension > 3)
  {
    return EXIT_SUCCESS;
  }
  vnl_vector_fixed<RealType, ImageDimension> evec_tert;
  l_rotate_around_tertiatry_and_secondary_axis(evec_tert, evec1_primary, evec1_2ndary);
  itk::Vector<RealType, ImageDimension> axis2;
  itk::Vector<RealType, ImageDimension> axis1;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    axis1[d] = evec_tert[d];
    axis2[d] = evec1_2ndary[d];
  }
  typename AffineType::Pointer affinesearch = AffineType::New();
  using OptimizerType = itk::MultiStartOptimizerv4;
  typename OptimizerType::MetricValuesListType metricvalues;
  typename OptimizerType::Pointer              mstartOptimizer = OptimizerType::New();
  using GCMetricType = itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType>;
  using MetricType = itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType>;
  typename MetricType::ParametersType newparams(affine1->GetParameters());
  typename GCMetricType::Pointer      gcmetric = GCMetricType::New();
  gcmetric->SetFixedImage(image1);
  gcmetric->SetVirtualDomainFromImage(image1);
  gcmetric->SetMovingImage(image2);
  gcmetric->SetMovingTransform(affinesearch);
  gcmetric->SetParameters(newparams);
  typename MetricType::Pointer mimetric = MetricType::New();
  mimetric->SetNumberOfHistogramBins(mibins);
  mimetric->SetFixedImage(image1);
  mimetric->SetMovingImage(image2);
  mimetric->SetMovingTransform(affinesearch);
  mimetric->SetParameters(newparams);
  if (mask.IsNotNull())
  {
    typename itk::ImageMaskSpatialObject<ImageDimension>::Pointer so =
      itk::ImageMaskSpatialObject<ImageDimension>::New();
    so->SetImage(const_cast<maskimagetype *>(mask.GetPointer()));
    mimetric->SetFixedImageMask(so);
    gcmetric->SetFixedImageMask(so);
  }
  using LocalOptimizerType = itk::ConjugateGradientLineSearchOptimizerv4;
  typename LocalOptimizerType::Pointer localoptimizer = LocalOptimizerType::New();
  RealType                             localoptimizerlearningrate = 0.1;
  localoptimizer->SetLearningRate(localoptimizerlearningrate);
  localoptimizer->SetMaximumStepSizeInPhysicalUnits(localoptimizerlearningrate);
  localoptimizer->SetNumberOfIterations(localSearchIterations);
  localoptimizer->SetLowerLimit(0);
  localoptimizer->SetUpperLimit(2);
  localoptimizer->SetEpsilon(0.1);
  localoptimizer->SetMaximumLineSearchIterations(10);
  localoptimizer->SetDoEstimateLearningRateOnce(true);
  localoptimizer->SetMinimumConvergenceValue(1.e-6);
  localoptimizer->SetConvergenceWindowSize(5);
  if (true)
  {
    using PointSetType = typename MetricType::FixedSampledPointSetType;
    using PointType = typename PointSetType::PointType;
    typename PointSetType::Pointer               pset(PointSetType::New());
    unsigned int                                 ind = 0;
    unsigned int                                 ct = 0;
    itk::ImageRegionIteratorWithIndex<ImageType> It(image1, image1->GetLargestPossibleRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      // take every N^th point
      if (ct % 10 == 0)
      {
        PointType pt;
        image1->TransformIndexToPhysicalPoint(It.GetIndex(), pt);
        pset->SetPoint(ind, pt);
        ind++;
      }
      ct++;
    }
    mimetric->SetFixedSampledPointSet(pset);
    mimetric->SetUseSampledPointSet(true);
    gcmetric->SetFixedSampledPointSet(pset);
    gcmetric->SetUseSampledPointSet(true);
  }
  if (whichMetric.compare("MI") == 0)
  {
    mimetric->Initialize();
    using RegistrationParameterScalesFromPhysicalShiftType =
      itk::RegistrationParameterScalesFromPhysicalShift<MetricType>;
    typename RegistrationParameterScalesFromPhysicalShiftType::Pointer shiftScaleEstimator =
      RegistrationParameterScalesFromPhysicalShiftType::New();
    shiftScaleEstimator->SetMetric(mimetric);
    shiftScaleEstimator->SetTransformForward(true);
    typename RegistrationParameterScalesFromPhysicalShiftType::ScalesType movingScales(
      affinesearch->GetNumberOfParameters());
    shiftScaleEstimator->EstimateScales(movingScales);
    mstartOptimizer->SetScales(movingScales);
    mstartOptimizer->SetMetric(mimetric);
    localoptimizer->SetMetric(mimetric);
    localoptimizer->SetScales(movingScales);
  }
  if (whichMetric.compare("MI") != 0)
  {
    gcmetric->Initialize();
    using RegistrationParameterScalesFromPhysicalShiftType =
      itk::RegistrationParameterScalesFromPhysicalShift<GCMetricType>;
    typename RegistrationParameterScalesFromPhysicalShiftType::Pointer shiftScaleEstimator =
      RegistrationParameterScalesFromPhysicalShiftType::New();
    shiftScaleEstimator->SetMetric(gcmetric);
    shiftScaleEstimator->SetTransformForward(true);
    typename RegistrationParameterScalesFromPhysicalShiftType::ScalesType movingScales(
      affinesearch->GetNumberOfParameters());
    shiftScaleEstimator->EstimateScales(movingScales);
    mstartOptimizer->SetScales(movingScales);
    mstartOptimizer->SetMetric(gcmetric);
    localoptimizer->SetMetric(gcmetric);
    localoptimizer->SetScales(movingScales);
  }
  typename OptimizerType::ParametersListType parametersList = mstartOptimizer->GetParametersList();
  for (double ang1 = (piover4 * (-1)); ang1 <= (piover4 + searchfactor); ang1 = ang1 + searchfactor)
  {
    if (useprincaxis)
    {
      affinesearch->SetMatrix(AA_solution);
    }
    if (ImageDimension == 3)
    {
      for (double ang2 = (piover4 * (-1)); ang2 <= (piover4 + searchfactor); ang2 = ang2 + searchfactor)
      {
        affinesearch->SetIdentity();
        affinesearch->SetCenter(trans2);
        affinesearch->SetOffset(trans);
        if (useprincaxis)
        {
          affinesearch->SetMatrix(AA_solution);
        }
        affinesearch->Rotate3D(axis1, ang1, 1);
        affinesearch->Rotate3D(axis2, ang2, 1);
        affinesearch->Scale(bestscale);
        parametersList.push_back(affinesearch->GetParameters());
      }
    }
    if (ImageDimension == 2)
    {
      affinesearch->SetIdentity();
      affinesearch->SetCenter(trans2);
      affinesearch->SetOffset(trans);
      if (useprincaxis)
      {
        affinesearch->SetMatrix(AA_solution);
      }
      affinesearch->Rotate2D(ang1, 1);
      affinesearch->Scale(bestscale);
      parametersList.push_back(affinesearch->GetParameters());
    }
  }
  mstartOptimizer->SetParametersList(parametersList);
  if (localSearchIterations > 0)
  {
    mstartOptimizer->SetLocalOptimizer(localoptimizer);
  }
  mstartOptimizer->StartOptimization();
  typename AffineType::Pointer bestaffine = AffineType::New();
  bestaffine->SetCenter(trans2);
  bestaffine->SetParameters(mstartOptimizer->GetBestParameters());


  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput(bestaffine);
  transformWriter->SetFileName(outname.c_str());
#if ITK_VERSION_MAJOR >= 5
  transformWriter->SetUseCompression(true);
#endif
  transformWriter->Update();
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsAffineInitializer(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsAffineInitializer");

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

  if (argc < 3)
  {
    std::cerr
      << "\nUsage: " << argv[0]
      << " ImageDimension <Image1.ext> <Image2.ext> TransformOutput.mat Optional-SearchFactor Optional-Radian-Fraction "
         "Optional-bool-UsePrincipalAxes Optional-uint-UseLocalSearch Optional-Image1Mask "
      << std::endl;
    std::cerr << " Optional-SearchFactor is in degrees --- e.g. 10 = search in 10 degree increments ." << std::endl;
    std::cerr << " Radian-Fraction should be between 0 and 1 --- will search this arc +/- around principal axis."
              << std::endl;
    std::cerr << " Optional-bool-UsePrincipalAxes determines whether the rotation is searched around an initial "
                 "principal axis alignment.  Default = false. "
              << std::endl;
    std::cerr << " Optional-uint-UseLocalSearch determines if a local optimization is run at each search point for the "
                 "set number of iterations. Default = 20."
              << std::endl;
    return 0;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return antsAffineInitializerImp<2>(argc, argv);
    }
    case 3:
    {
      return antsAffineInitializerImp<3>(argc, argv);
    }
    case 4:
      return antsAffineInitializerImp<4>(argc, argv);
  }
  return antsAffineInitializerImp<2>(argc, argv);
}
} // namespace ants
