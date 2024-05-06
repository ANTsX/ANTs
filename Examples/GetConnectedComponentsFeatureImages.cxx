#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include <string>
#include <vector>

namespace ants
{
template <unsigned int ImageDimension>
int
GetConnectedComponentsFeatureImages(int itkNotUsed(argc), char * argv[])
{
  using PixelType = int;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using RealImageType = itk::Image<float, ImageDimension>;

  typename ImageType::Pointer inputImage = nullptr;
  ReadImage<ImageType>(inputImage, argv[2]);

  // Output images:
  // [0] = volume (in physical coordinates)
  // [1] = volume / surface area
  // [2] = eccentricity
  // [3] = elongation

  std::vector<typename RealImageType::Pointer> outputImages;
  for (unsigned int n = 0; n < 4; n++)
  {
    typename RealImageType::Pointer output = RealImageType::New();
    output->CopyInformation(inputImage);
    output->SetRegions(inputImage->GetRequestedRegion());
    output->AllocateInitialized();

    outputImages.push_back(output);
  }

  typename ImageType::SpacingType spacing = inputImage->GetSpacing();

  using RelabelerType = itk::RelabelComponentImageFilter<ImageType, ImageType>;
  typename RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput(inputImage);
  relabeler->Update();

  for (unsigned int i = 1; i <= relabeler->GetNumberOfObjects(); i++)
  {
    using ThresholderType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput(relabeler->GetOutput());
    thresholder->SetLowerThreshold(i);
    thresholder->SetUpperThreshold(i);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);
    thresholder->Update();

    using ConnectedComponentType = itk::ConnectedComponentImageFilter<ImageType, ImageType>;
    typename ConnectedComponentType::Pointer filter = ConnectedComponentType::New();
    filter->SetInput(thresholder->GetOutput());
    filter->Update();

    typename RelabelerType::Pointer relabeler2 = RelabelerType::New();
    relabeler2->SetInput(filter->GetOutput());
    relabeler2->Update();

    using GeometryFilterType = itk::LabelImageToShapeLabelMapFilter<ImageType>;
    typename GeometryFilterType::Pointer geometry = GeometryFilterType::New();
    geometry->SetInput(relabeler2->GetOutput());
    geometry->ComputeOrientedBoundingBoxOff();
    geometry->ComputePerimeterOn();

    geometry->Update();

    itk::ImageRegionIteratorWithIndex<ImageType> It(relabeler->GetOutput(),
                                                    relabeler->GetOutput()->GetRequestedRegion());
    itk::ImageRegionIterator<ImageType> It2(relabeler2->GetOutput(), relabeler2->GetOutput()->GetRequestedRegion());

    for (It.GoToBegin(), It2.GoToBegin(); !It.IsAtEnd(); ++It, ++It2)
    {
      int label = It2.Get();
      if (label != 0)
      {
        typename ImageType::IndexType index = It.GetIndex();

        // Output images:
        // [0] = volume (in physical coordinates)
        // [1] = volume / surface area
        // [2] = eccentricity
        // [3] = elongation

        try
        {
          auto labelObject = geometry->GetOutput()->GetLabelObject(label);
          float volume = labelObject->GetPhysicalSize();

          outputImages[0]->SetPixel(index, volume);
          outputImages[1]->SetPixel(index, volume / static_cast<float>(labelObject->GetPerimeter()));

          auto principalMoments = labelObject->GetPrincipalMoments();

          float lambda1 = principalMoments[0];
          float lambdaN = principalMoments[ImageDimension - 1];
          float eccentricity = 0.0;

          if (!itk::Math::FloatAlmostEqual(lambda1, 0.0f))
          {
            eccentricity = std::sqrt(1.0 - (lambda1 * lambda1) / (lambdaN * lambdaN));
          }

          outputImages[2]->SetPixel(index, eccentricity);
          outputImages[3]->SetPixel(index, labelObject->GetElongation());
        }
        catch (itk::ExceptionObject & e)
        {
          std::cerr << "Could not find label " << label << std::endl;
        }
      }
    }
  }

  {
    std::string filename = std::string(argv[3]) + std::string("PHYSICAL_VOLUME.nii.gz");
    ANTs::WriteImage<RealImageType>(outputImages[0], filename.c_str());
  }

  {
    std::string filename = std::string(argv[3]) + std::string("VOLUME_TO_SURFACE_AREA_RATIO.nii.gz");
    ANTs::WriteImage<RealImageType>(outputImages[1], filename.c_str());
  }

  {
    std::string filename = std::string(argv[3]) + std::string("ECCENTRICITY.nii.gz");
    ANTs::WriteImage<RealImageType>(outputImages[2], filename.c_str());
  }

  {
    std::string filename = std::string(argv[3]) + std::string("ELONGATION.nii.gz");
    ANTs::WriteImage<RealImageType>(outputImages[3], filename.c_str());
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
GetConnectedComponentsFeatureImages(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "GetConnectedComponentsFeatureImages");

  const int argc = args.size();
  char **   argv = new char *[args.size() + 1];
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

  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " imageDimension "
              << "inputSegmentationImage outputImagePrefix" << std::endl;
    return EXIT_FAILURE;
  }

  int returnStatus = EXIT_FAILURE;
  switch (std::stoi(argv[1]))
  {
    case 2:
      returnStatus = GetConnectedComponentsFeatureImages<2>(argc, argv);
      break;
    case 3:
      returnStatus = GetConnectedComponentsFeatureImages<3>(argc, argv);
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
  }
  return returnStatus;
}
} // namespace ants
