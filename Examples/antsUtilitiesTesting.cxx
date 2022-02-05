#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include <cstdio>

#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkSimilarity2DTransform.h"

#include <string>
#include <vector>

namespace ants
{
int
antsUtilitiesTesting(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsUtilitiesTesting");

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
    std::cerr << "Usage: " << argv[0]
              << " testImage whichMetric[GCorMIorMattes] numberOfIterations lowerScalexUpperScalexNumberOfScaleSamples "
                 "numberOfRotationSamples doPrintMatchingFileName transformName setOfTrainingImages"
              << std::endl;
    std::cerr << "Notes:  if transformName='none', no transform is printed" << std::endl;
    std::cerr << "Example call:  " << std::endl;
    std::cerr << "   " << argv[0] << "test.nii.gz CC 10 0.25x3x11 15 testTransform.txt 0 training*.nii.gz" << std::endl;

    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  constexpr unsigned int ImageDimension = 2;

  using PixelType = double;

  using ImageType = itk::Image<PixelType, ImageDimension>;

  // read in test image

  ImageType::Pointer testImage = ImageType::New();
  ReadImage<ImageType>(testImage, argv[1]);

  // read in training image file names

  std::vector<std::string> trainingImageFileNames;
  for (int n = 8; n < argc; n++)
  {
    trainingImageFileNames.emplace_back(argv[n]);
  }

  // Get number of iterations

  unsigned int numberOfIterations = std::stoi(argv[3]);

  // Get scale parameters

  std::vector<double> scaleParameters = ConvertVector<double>(std::string(argv[4]));
  if (scaleParameters.size() != 3)
  {
    std::cerr << "The scale parameters were improperly specified.  See usage." << std::endl;
    return EXIT_FAILURE;
  }
  double scaleLowerBoundLog = std::log(scaleParameters[0]);
  double scaleUpperBoundLog = std::log(scaleParameters[1]);
  auto   scaleNumberOfSamples = static_cast<unsigned int>(scaleParameters[2]);
  double scaleDelta = (scaleUpperBoundLog - scaleLowerBoundLog) / static_cast<double>(scaleNumberOfSamples - 1);

  // Get rotation sampling resolution

  unsigned int rotationNumberOfSamples = static_cast<unsigned int>(std::stoi(argv[5]));
  double       rotationDelta = (2.0 * itk::Math::pi - 0.0) / static_cast<double>(rotationNumberOfSamples - 1);

  // Set up metric
  using ImageMetricType = itk::ImageToImageMetricv4<ImageType, ImageType, ImageType>;
  using PointSetType = ImageMetricType::FixedSampledPointSetType;

  ImageMetricType::Pointer imageMetric = nullptr;

  if (strcmp(argv[2], "Mattes") == 0)
  {
    using MattesMetricType = itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType>;
    MattesMetricType::Pointer mattesMetric = MattesMetricType::New();
    mattesMetric->SetNumberOfHistogramBins(20);

    imageMetric = mattesMetric;
  }
  else if (strcmp(argv[2], "GC") == 0)
  {
    using GCMetricType = itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType>;
    GCMetricType::Pointer gcMetric = GCMetricType::New();

    imageMetric = gcMetric;
  }
  else if (strcmp(argv[2], "MI") == 0)
  {
    using MIMetricType = itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType>;
    MIMetricType::Pointer miMetric = MIMetricType::New();
    miMetric->SetNumberOfHistogramBins(20);

    imageMetric = miMetric;
  }
  else
  {
    std::cerr << "Unrecognized metric option." << std::endl;
    return EXIT_FAILURE;
  }

  imageMetric->SetFixedImage(testImage);
  imageMetric->SetVirtualDomainFromImage(testImage);
  imageMetric->SetUseMovingImageGradientFilter(false);
  imageMetric->SetUseFixedImageGradientFilter(false);

  // identity transform for fixed image

  using IdentityTransformType = itk::IdentityTransform<double, ImageDimension>;
  IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
  identityTransform->SetIdentity();

  imageMetric->SetFixedTransform(identityTransform);

  // Do a random sampling

  unsigned int index = 0;
  unsigned int count = 0;

  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  itk::ImageRegionIteratorWithIndex<ImageType> It(testImage, testImage->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    // take every N^th point
    if (count % 20 == 0)
    {
      PointSetType::PointType point;
      testImage->TransformIndexToPhysicalPoint(It.GetIndex(), point);
      pointSet->SetPoint(index++, point);
    }
    count++;
  }
  imageMetric->SetFixedSampledPointSet(pointSet);
  imageMetric->SetUseSampledPointSet(true);

  // Now go through the rotations + scalings to find the optimal pose.

  using AffineTransformType = itk::AffineTransform<double, ImageDimension>;
  using SimilarityTransformType = itk::Similarity2DTransform<double>;

  double                           optimalMetricValue = itk::NumericTraits<double>::max();
  SimilarityTransformType::Pointer optimalTransform = nullptr;
  unsigned int                     optimalMetricIndex = 0;

  for (unsigned int n = 0; n < trainingImageFileNames.size(); n++)
  {
    double optimalMetricValuePerImage = itk::NumericTraits<double>::max();

    ImageType::Pointer trainingImage = ImageType::New();
    ReadImage<ImageType>(trainingImage, trainingImageFileNames[n].c_str());

    imageMetric->SetMovingImage(trainingImage);

    // Initialize centered transform (based on the center of the image)

    AffineTransformType::Pointer initialTransform = AffineTransformType::New();

    using TransformInitializerType = itk::CenteredTransformInitializer<AffineTransformType, ImageType, ImageType>;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform(initialTransform);
    initializer->SetFixedImage(testImage);
    initializer->SetMovingImage(trainingImage);
    initializer->GeometryOn();

    initializer->InitializeTransform();

    for (double angle = 0.0; angle < 2.0 * itk::Math::pi; angle += rotationDelta)
    {
      AffineTransformType::MatrixType rotationMatrix;
      rotationMatrix(0, 0) = rotationMatrix(1, 1) = std::cos(angle);
      rotationMatrix(1, 0) = std::sin(angle);
      rotationMatrix(0, 1) = -rotationMatrix(1, 0);

      for (double scaleLog = scaleLowerBoundLog; scaleLog <= scaleUpperBoundLog; scaleLog += scaleDelta)
      {
        double scale = std::exp(scaleLog);

        SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();
        similarityTransform->SetCenter(initialTransform->GetCenter());
        similarityTransform->SetTranslation(initialTransform->GetTranslation());
        similarityTransform->SetMatrix(initialTransform->GetMatrix() * rotationMatrix);
        similarityTransform->SetScale(scale);

        imageMetric->SetMovingTransform(similarityTransform);
        imageMetric->Initialize();

        if (numberOfIterations > 0)
        {
          using ScalesEstimatorType = itk::RegistrationParameterScalesFromPhysicalShift<ImageMetricType>;
          ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
          scalesEstimator->SetMetric(imageMetric);
          scalesEstimator->SetTransformForward(true);

          using ConjugateGradientDescentOptimizerType = itk::ConjugateGradientLineSearchOptimizerv4Template<double>;
          ConjugateGradientDescentOptimizerType::Pointer optimizer = ConjugateGradientDescentOptimizerType::New();
          optimizer->SetLowerLimit(0);
          optimizer->SetUpperLimit(2);
          optimizer->SetEpsilon(0.2);
          optimizer->SetLearningRate(0.15);
          optimizer->SetMaximumStepSizeInPhysicalUnits(0.15);
          optimizer->SetNumberOfIterations(numberOfIterations);
          optimizer->SetScalesEstimator(scalesEstimator);
          optimizer->SetMinimumConvergenceValue(1e-10);
          optimizer->SetConvergenceWindowSize(10);
          optimizer->SetDoEstimateLearningRateAtEachIteration(false);
          optimizer->SetDoEstimateLearningRateOnce(true);
          optimizer->SetMetric(imageMetric);

          //         typedef itk::GradientDescentOptimizerv4Template<double> GradientDescentOptimizerType;
          //         GradientDescentOptimizerType::Pointer optimizer2 = GradientDescentOptimizerType::New();
          //         optimizer2->SetLearningRate( 0.15 );
          //         optimizer2->SetMaximumStepSizeInPhysicalUnits( 0.15 );
          //         optimizer2->SetNumberOfIterations( numberOfIterations );
          //         optimizer2->SetScalesEstimator( scalesEstimator );
          //         optimizer2->SetMinimumConvergenceValue( 1e-6 );
          //         optimizer2->SetConvergenceWindowSize( 10 );
          //         optimizer2->SetDoEstimateLearningRateAtEachIteration( true );
          //         optimizer2->SetDoEstimateLearningRateOnce( false );
          //         optimizer2->SetMetric( imageMetric );

          try
          {
            optimizer->StartOptimization();
          }
          catch (...)
          {
            continue;
          }
        }

        //         typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResamplerType;
        //         ResamplerType::Pointer resampleFilter = ResamplerType::New();
        //         resampleFilter->SetInput( trainingImage );
        //         resampleFilter->SetOutputParametersFromImage( testImage );
        //         resampleFilter->SetTransform( similarityTransform );
        //         resampleFilter->SetDefaultPixelValue( 255 );
        //         resampleFilter->Update();
        //         ANTs::WriteImage<ImageType>( resampleFilter->GetOutput(), "training.nii.gz" );

        double metricValue = imageMetric->GetValue();

        if (metricValue < optimalMetricValue)
        {
          optimalMetricValue = metricValue;
          optimalTransform = similarityTransform;
          optimalMetricIndex = n;
        }

        if (metricValue < optimalMetricValuePerImage)
        {
          optimalMetricValuePerImage = metricValue;
        }
      }
    }
    //     std::cout << trainingImageFileNames[n] << " -> " << optimalMetricValuePerImage << std::endl;
  }

  if (strcmp(argv[7], "none") != 0)
  {
    using TransformWriterType = itk::TransformFileWriter;
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetInput(optimalTransform);
    transformWriter->SetFileName(argv[7]);
#if ITK_VERSION_MAJOR >= 5
    transformWriter->SetUseCompression(true);
#endif
    transformWriter->Update();
  }

  if (static_cast<bool>(std::stoi(argv[6])))
  {
    std::cout << trainingImageFileNames[optimalMetricIndex] << std::endl;
  }
  std::cout << optimalMetricValue << std::endl;

  return EXIT_SUCCESS;
}

} // namespace ants
