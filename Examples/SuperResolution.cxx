#include "antsUtilities.h"
#include <algorithm>

#include "itkAddImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkPointSet.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int
SuperResolution(unsigned int argc, char * argv[])
{
  using RealType = float;
  using ImageType = itk::Image<RealType, ImageDimension>;

  using ScalarType = itk::Vector<RealType, 1>;
  using ScalarImageType = itk::Image<ScalarType, ImageDimension>;
  using PointSetType = itk::PointSet<ScalarType, ImageDimension>;
  using BSplineFilterType = itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, ScalarImageType>;

  typename ImageType::Pointer domainImage = nullptr;
  ReadImage<ImageType>(domainImage, argv[3]);

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  typename PointSetType::Pointer bsplinePoints = PointSetType::New();
  bsplinePoints->Initialize();

  typename BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();
  weights->Initialize();

  unsigned int                          splineOrder = 3;
  typename BSplineFilterType::ArrayType numberOfLevels;
  typename BSplineFilterType::ArrayType ncps;

  bool     useGradientWeighting = true;
  RealType gradientSigma = atof(argv[4]);
  if (gradientSigma < itk::NumericTraits<RealType>::ZeroValue())
  {
    useGradientWeighting = false;
  }

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>(std::string(argv[5]));
  if (meshSize.size() == 1)
  {
    ncps.Fill(meshSize[0] + splineOrder);
  }
  else if (meshSize.size() == ImageDimension)
  {
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      ncps[d] = meshSize[d] + splineOrder;
    }
  }
  else
  {
    std::cerr << "Invalid ncps format." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<unsigned int> nlevels = ConvertVector<unsigned int>(std::string(argv[6]));
  if (nlevels.size() == 1)
  {
    numberOfLevels.Fill(nlevels[0]);
  }
  else if (nlevels.size() == ImageDimension)
  {
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      numberOfLevels[d] = nlevels[d];
    }
  }
  else
  {
    std::cerr << "Invalid nlevels format." << std::endl;
    return EXIT_FAILURE;
  }

  // Find the average for the offset

  typename ImageType::IndexType domainBeginIndex = domainImage->GetRequestedRegion().GetIndex();
  typename ImageType::IndexType domainEndIndex;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    domainEndIndex[d] = domainBeginIndex[d] + domainImage->GetRequestedRegion().GetSize()[d] - 1;
  }

  RealType     averageIntensity = 0.0;
  unsigned int N = 0;
  for (unsigned int n = 7; n < argc; n++)
  {
    typename ImageType::Pointer inputImage = nullptr;
    ReadImage<ImageType>(inputImage, argv[n]);

    itk::ImageRegionConstIteratorWithIndex<ImageType> It(inputImage, inputImage->GetRequestedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      typename ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint(It.GetIndex(), imagePoint);

      itk::ContinuousIndex<RealType, ImageDimension> cidx;
      bool isInside = domainImage->TransformPhysicalPointToContinuousIndex(imagePoint, cidx);

      if (isInside)
      {
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          if (cidx[d] <= domainBeginIndex[d] || cidx[d] >= domainEndIndex[d])
          {
            isInside = false;
            break;
          }
        }
      }

      if (!isInside)
      {
        continue;
      }

      averageIntensity += It.Get();

      N++;
    }
  }
  averageIntensity /= static_cast<RealType>(N);


  typename ScalarImageType::DirectionType identity;
  identity.SetIdentity();

  const typename ScalarImageType::RegionType & bufferedRegion = domainImage->GetBufferedRegion();
  const itk::SizeValueType                     numberOfPixels = bufferedRegion.GetNumberOfPixels();
  const bool                                   filterHandlesMemory = false;

  using ImporterType = itk::ImportImageFilter<RealType, ImageDimension>;
  typename ImporterType::Pointer importer = ImporterType::New();
  importer->SetImportPointer(domainImage->GetBufferPointer(), numberOfPixels, filterHandlesMemory);
  importer->SetRegion(domainImage->GetBufferedRegion());
  importer->SetOrigin(domainImage->GetOrigin());
  importer->SetSpacing(domainImage->GetSpacing());
  importer->SetDirection(identity);
  importer->Update();

  const typename ImporterType::OutputImageType * parametricDomainImage = importer->GetOutput();

  N = 0;
  for (unsigned int n = 7; n < argc; n++)
  {
    typename ImageType::Pointer inputImage = nullptr;
    ReadImage<ImageType>(inputImage, argv[n]);

    typename ImageType::Pointer gradientImage = nullptr;
    if (useGradientWeighting)
    {
      using GradientFilterType = itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType>;
      typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
      gradientFilter->SetSigma(gradientSigma);
      gradientFilter->SetInput(inputImage);

      using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
      typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
      rescaler->SetOutputMinimum(0.0);
      rescaler->SetOutputMaximum(1.0);
      rescaler->SetInput(gradientFilter->GetOutput());

      gradientImage = rescaler->GetOutput();
      gradientImage->Update();
      gradientImage->DisconnectPipeline();
    }

    itk::ImageRegionConstIteratorWithIndex<ImageType> It(inputImage, inputImage->GetRequestedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      typename ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint(It.GetIndex(), imagePoint);

      itk::ContinuousIndex<RealType, ImageDimension> cidx;
      bool isInside = domainImage->TransformPhysicalPointToContinuousIndex(imagePoint, cidx);

      if (isInside)
      {
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          if (cidx[d] <= domainBeginIndex[d] || cidx[d] >= domainEndIndex[d])
          {
            isInside = false;
            break;
          }
        }
      }

      if (!isInside)
      {
        continue;
      }

      RealType weight = 1.0;
      if (gradientImage.IsNotNull())
      {
        weight = gradientImage->GetPixel(It.GetIndex());
      }
      if (itk::Math::FloatAlmostEqual(weight, itk::NumericTraits<RealType>::ZeroValue()))
      {
        continue;
      }

      //       domainImage->SetPixel( index, 1 );

      parametricDomainImage->TransformContinuousIndexToPhysicalPoint(cidx, imagePoint);

      ScalarType scalar;
      scalar[0] = It.Get() - averageIntensity;

      bsplinePoints->SetPointData(N, scalar);
      bsplinePoints->SetPoint(N, imagePoint);
      weights->InsertElement(N, weight);

      N++;
    }
  }

  itk::TimeProbe timer;
  timer.Start();

  typename ScalarImageType::PointType parametricOrigin = domainImage->GetOrigin();
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    parametricOrigin[d] += (domainImage->GetSpacing()[d] * domainImage->GetLargestPossibleRegion().GetIndex()[d]);
  }

  bspliner->SetOrigin(parametricOrigin);
  bspliner->SetSpacing(domainImage->GetSpacing());
  bspliner->SetSize(domainImage->GetRequestedRegion().GetSize());
  bspliner->SetDirection(domainImage->GetDirection());
  bspliner->SetGenerateOutputImage(true);
  bspliner->SetNumberOfLevels(numberOfLevels);
  bspliner->SetSplineOrder(splineOrder);
  bspliner->SetNumberOfControlPoints(ncps);
  bspliner->SetInput(bsplinePoints);
  bspliner->SetPointWeights(weights);
  bspliner->Update();

  timer.Stop();

  std::cout << "Elapsed Time:  " << timer.GetMean() << std::endl;

  using SelectorType = itk::VectorIndexSelectionCastImageFilter<ScalarImageType, ImageType>;
  typename SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput(bspliner->GetOutput());
  selector->SetIndex(0);
  selector->Update();

  using AdderType = itk::AddImageFilter<ImageType>;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput(selector->GetOutput());
  adder->SetConstant2(averageIntensity);
  adder->Update();

  ANTs::WriteImage<ImageType>(adder->GetOutput(), argv[2]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
SuperResolution(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "SuperResolution");

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

  if (argc < 8)
  {
    std::cerr
      << argv[0]
      << " imageDimension outputImage domainImage gradientSigma meshSize numberOfLevels inputImage1 ... inputImageN"
      << std::endl;
    std::cerr << std::endl;
    std::cerr
      << "    N.B. The \'gradientSigma\' parameter is used in calculating the gradient magnitude of the input images "
      << std::endl;
    std::cerr
      << "       for weighting the voxel points during fitting.  If a negative \'gradient\' sigma is specified then no "
      << std::endl;
    std::cerr << "       weighting is used." << std::endl;

    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  const int ImageDimension = static_cast<int>(std::stoi(argv[1]));

  switch (ImageDimension)
  {
    case 2:
      return SuperResolution<2>(argc, argv);
      break;
    case 3:
      return SuperResolution<3>(argc, argv);
      break;
    case 4:
      return SuperResolution<4>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
} // namespace ants
