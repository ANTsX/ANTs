#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include <cstdio>

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkDenseFrequencyContainer2.h"

namespace ants
{
template <unsigned int ImageDimension>
int
TextureCooccurrenceFeatures(int argc, char * argv[])
{

  using PixelType = float;
  using RealType = float;


  using ImageType = itk::Image<PixelType, ImageDimension>;
  using RealImageType = itk::Image<RealType, ImageDimension>;

  typename ImageType::Pointer inputImage = ImageType::New();
  ReadImage<ImageType>(inputImage, argv[2]);

  // we have to rescale the input image since the cooccurrence filter
  // adds 1 to the upper bound of the joint histogram and small ranges
  // will be greatly affected by this.

  using RescalerType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput(inputImage);
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(10000);
  rescaler->Update();

  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;

  using TextureFilterType =
    itk::Statistics::ScalarImageToTextureFeaturesFilter<RealImageType, HistogramFrequencyContainerType>;
  typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  textureFilter->SetInput(rescaler->GetOutput());

  typename ImageType::Pointer mask = nullptr;
  PixelType                   label = itk::NumericTraits<PixelType>::OneValue();
  if (argc > 4)
  {
    ReadImage<ImageType>(mask, argv[4]);
    textureFilter->SetMaskImage(mask);

    if (argc > 5)
    {
      label = static_cast<PixelType>(std::stoi(argv[5]));
    }
    textureFilter->SetInsidePixelValue(label);
  }

  unsigned int numberOfBins = 256;
  if (argc > 3)
  {
    numberOfBins = static_cast<PixelType>(std::stoi(argv[3]));
  }
  textureFilter->SetNumberOfBinsPerAxis(numberOfBins);

  itk::ImageRegionIteratorWithIndex<ImageType> ItI(rescaler->GetOutput(),
                                                   rescaler->GetOutput()->GetLargestPossibleRegion());

  PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (!mask || itk::Math::FloatAlmostEqual(mask->GetPixel(ItI.GetIndex()), label))
    {
      if (ItI.Get() < minValue)
      {
        minValue = ItI.Get();
      }
      else if (ItI.Get() > maxValue)
      {
        maxValue = ItI.Get();
      }
    }
  }

  textureFilter->SetPixelValueMinMax(minValue, maxValue);
  textureFilter->SetNumberOfBinsPerAxis(numberOfBins);
  textureFilter->FastCalculationsOff();

  /**
   * Second order measurements
   * These include:
   *   1. energy (f1) *cth, *amfm
   *   2. entropy (f2) *cth, *amfm
   *   3. correlation (f3) *amfm
   *   4. inverse difference moment (f4) *cth, *amfm
   *   5. inertia (f5) *cth, *amfm
   *   6. cluster shade (f6) *cth
   *   7. cluster prominence (f7) *cth
   *   8. haralick's correlation (f8)
   */

  textureFilter->Update();

  typename TextureFilterType::FeatureValueVectorPointer means = textureFilter->GetFeatureMeans();
  const typename TextureFilterType::FeatureNameVector * names = textureFilter->GetRequestedFeatures();

  typename TextureFilterType::FeatureValueVector::ConstIterator mIt = means->Begin();
  typename TextureFilterType::FeatureNameVector::ConstIterator  nIt = names->Begin();


  while (mIt != means->End())
  {
    //    std::cout << nIt.Value() << ": " << mIt.Value() << std::endl;
    std::cout << mIt.Value() << " ";
    ++mIt;
    ++nIt;
  }
  std::cout << std::endl;

  //  RealType entropy = textureFilter->GetFeatureMeansOutput()[1];
  //  RealType inverseDifferenceMoment = textureFilter->GetFeatureMeansOutput()[2];
  //  RealType inertia = textureFilter->GetFeatureMeansOutput()[3];
  //  RealType clusterShade = textureFilter->GetFeatureMeansOutput()[4];
  //  RealType clusterProminence = textureFilter->GetFeatureMeansOutput()[5];
  //
  //  std::cout << energy << " "
  //            << entropy << " "
  //            << inverseDifferenceMoment << " "
  //            << inertia << " "
  //            << clusterShade << " "
  //            << clusterProminence << std::endl;

  /*
    std::cout << "energy             : " << energy << std::endl;
    std::cout << "entropy            : " << entropy << std::endl;
    std::cout << "inverse diff moment: " << inverseDifferenceMoment << std::endl;
    std::cout << "inertia            : " << inertia << std::endl;
    std::cout << "cluster Shade      : " << clusterShade << std::endl;
    std::cout << "cluster prominence : " << clusterProminence << std::endl;
  */

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
TextureCooccurrenceFeatures(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "TextureCooccurrenceFeatures");

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
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage "
              << "[numberOfBinsPerAxis=256] [maskImage] [maskLabel=1]" << std::endl;
    std::cerr << "Features: Energy,Entropy,InverseDifferenceMoment,Inertia,ClusterShade,ClusterProminence" << std::endl;
    exit(1);
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      TextureCooccurrenceFeatures<2>(argc, argv);
      break;
    case 3:
      TextureCooccurrenceFeatures<3>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
} // namespace ants
