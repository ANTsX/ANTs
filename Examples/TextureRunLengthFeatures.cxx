#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include <cstdio>

#include "itkBoundingBox.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"

#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkDenseFrequencyContainer2.h"

namespace ants
{
template <unsigned int ImageDimension>
int
TextureRunLengthFeatures(int argc, char * argv[])
{
  using PixelType = float;
  using RealType = float;

  using ImageType = itk::Image<PixelType, ImageDimension>;
  using RealImageType = itk::Image<RealType, ImageDimension>;

  typename ImageType::Pointer inputImage = ImageType::New();
  ReadImage<ImageType>(inputImage, argv[2]);

  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;

  using RunLengthFilterType =
    itk::Statistics::ScalarImageToRunLengthFeaturesFilter<RealImageType, HistogramFrequencyContainerType>;
  typename RunLengthFilterType::Pointer runLengthFilter = RunLengthFilterType::New();
  runLengthFilter->SetInput(inputImage);

  typename ImageType::Pointer mask = nullptr;
  PixelType                   label = itk::NumericTraits<PixelType>::OneValue();
  if (argc > 4)
  {
    ReadImage<ImageType>(mask, argv[4]);
    runLengthFilter->SetMaskImage(mask);

    if (argc > 5)
    {
      label = static_cast<PixelType>(std::stoi(argv[5]));
    }
    runLengthFilter->SetInsidePixelValue(label);
  }


  unsigned int numberOfBins = 256;
  if (argc > 3)
  {
    numberOfBins = static_cast<PixelType>(std::stoi(argv[3]));
  }
  runLengthFilter->SetNumberOfBinsPerAxis(numberOfBins);


  itk::ImageRegionIteratorWithIndex<ImageType> ItI(inputImage, inputImage->GetLargestPossibleRegion());

  PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  using BoundingBoxType = itk::BoundingBox<unsigned long, ImageDimension, RealType>;
  typename BoundingBoxType::Pointer                bbox = BoundingBoxType::New();
  typename BoundingBoxType::PointsContainerPointer points = BoundingBoxType::PointsContainer::New();
  itk::Point<RealType, ImageDimension>             point;

  unsigned int idx = 0;

  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (!mask || (itk::Math::FloatAlmostEqual(mask->GetPixel(ItI.GetIndex()), label)))
    {
      if (ItI.Get() < minValue)
      {
        minValue = ItI.Get();
      }
      else if (ItI.Get() > maxValue)
      {
        maxValue = ItI.Get();
      }
      inputImage->TransformIndexToPhysicalPoint(ItI.GetIndex(), point);
      points->InsertElement(idx++, point);
    }
  }
  bbox->SetPoints(points);
  bbox->ComputeBoundingBox();
  typename BoundingBoxType::PointType pointMin = bbox->GetMinimum();
  typename BoundingBoxType::PointType pointMax = bbox->GetMaximum();

  runLengthFilter->SetPixelValueMinMax(minValue, maxValue);
  runLengthFilter->SetDistanceValueMinMax(0, pointMin.EuclideanDistanceTo(pointMax));
  runLengthFilter->SetNumberOfBinsPerAxis(numberOfBins);
  runLengthFilter->FastCalculationsOff();

  runLengthFilter->Update();

  typename RunLengthFilterType::FeatureValueVectorPointer means = runLengthFilter->GetFeatureMeans();
  const typename RunLengthFilterType::FeatureNameVector * names = runLengthFilter->GetRequestedFeatures();

  typename RunLengthFilterType::FeatureValueVector::ConstIterator mIt = means->Begin();
  typename RunLengthFilterType::FeatureNameVector::ConstIterator  nIt = names->Begin();

  while (mIt != means->End())
  {
    //    std::cout << nIt.Value() << ": " << mIt.Value() << std::endl;
    std::cout << mIt.Value() << " ";
    ++mIt;
    ++nIt;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
TextureRunLengthFeatures(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "TextureRunLengthFeatures");

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

  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
              << "[numberOfBinsPerAxis=256] [maskImage] [maskLabel=1]" << std::endl;
    std::cout << "Features: ShortRunEmphasis,LongRunEmphasis,GreyLevelNonuniformity,";
    std::cout << "RunLengthNonuniformity,LowGreyLevelRunEmphasis,HighGreyLevelRunEmphasis,";
    std::cout << "ShortRunLowGreyLevelEmphasis,ShortRunHighGreyLevelEmphasis,";
    std::cout << "LongRunLowGreyLevelEmphasis,LongRunHighGreyLevelEmphasis" << std::endl;

    exit(1);
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      TextureRunLengthFeatures<2>(argc, argv);
      break;
    case 3:
      TextureRunLengthFeatures<3>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
} // namespace ants
