
#include "antsUtilities.h"
#include "ReadWriteData.h"
#include <algorithm>

#include <cstdio>

#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <string>
#include <vector>

namespace ants
{

template <unsigned int ImageDimension>
int
ExtractRegionFromImage(int argc, char * argv[])
{
  using PixelType = float;

  using ImageType = itk::Image<PixelType, ImageDimension>;

  typename ImageType::Pointer inputImage = ImageType::New();
  ReadImage<ImageType>(inputImage, argv[2]);

  typename ImageType::RegionType            region;
  typename ImageType::RegionType::SizeType  regionSize;
  typename ImageType::RegionType::IndexType regionIndex;
  if (argc == 6)
  {
    std::vector<int> minIndex;
    std::vector<int> maxIndex;
    minIndex = ConvertVector<int>(std::string(argv[4]));
    maxIndex = ConvertVector<int>(std::string(argv[5]));
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      regionIndex[i] = minIndex[i];
      regionSize[i] = maxIndex[i] - minIndex[i] + 1;
    }
    region.SetSize(regionSize);
    region.SetIndex(regionIndex);
  }
  else if (argc == 7)
  {
    typename ImageType::Pointer labimg;
    ReadImage<ImageType>(labimg, argv[5]);
    using ShortImageType = itk::Image<unsigned short, ImageDimension>;
    using CasterType = itk::CastImageFilter<ImageType, ShortImageType>;
    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput(labimg);
    caster->Update();

    using StatsFilterType = itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetLabelInput(caster->GetOutput());
    stats->SetInput(caster->GetOutput());
    stats->Update();

    region = stats->GetRegion(std::stoi(argv[4]));
  }
  else
  {
    typename ImageType::Pointer domainImage = nullptr;
    ReadImage<ImageType>(domainImage, argv[4]);

    if (domainImage.IsNotNull())
    {
      typename ImageType::IndexType maxIndex;
      typename ImageType::IndexType minIndex;

      minIndex.Fill(itk::NumericTraits<int>::max());
      maxIndex.Fill(itk::NumericTraits<int>::NonpositiveMin());

      itk::ImageRegionConstIteratorWithIndex<ImageType> It(inputImage, inputImage->GetLargestPossibleRegion());

      for (It.GoToBegin(); !It.IsAtEnd(); ++It)
      {
        typename ImageType::IndexType index = It.GetIndex();

        typename ImageType::PointType point;
        inputImage->TransformIndexToPhysicalPoint(index, point);

        typename ImageType::IndexType trashIndex;
        bool                          isInside = domainImage->TransformPhysicalPointToIndex(point, trashIndex);

        if (isInside)
        {
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            if (index[d] < minIndex[d])
            {
              minIndex[d] = index[d];
            }
            if (index[d] > maxIndex[d])
            {
              maxIndex[d] = index[d];
            }
          }
        }
      }

      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        regionIndex[i] = minIndex[i];
        regionSize[i] = maxIndex[i] - regionIndex[i] + 1;
      }
      region.SetSize(regionSize);
      region.SetIndex(regionIndex);
    }
    else
    {
      using ShortImageType = itk::Image<unsigned short, ImageDimension>;
      using CasterType = itk::CastImageFilter<ImageType, ShortImageType>;
      typename CasterType::Pointer caster = CasterType::New();
      caster->SetInput(inputImage);
      caster->Update();

      using StatsFilterType = itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>;
      typename StatsFilterType::Pointer stats = StatsFilterType::New();
      stats->SetLabelInput(caster->GetOutput());
      stats->SetInput(caster->GetOutput());
      stats->Update();

      region = stats->GetRegion(std::stoi(argv[4]));
    }
  }

  using CropperType = itk::ExtractImageFilter<ImageType, ImageType>;
  typename CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput(inputImage);
  cropper->SetExtractionRegion(region);
  cropper->SetDirectionCollapseToSubmatrix();
  cropper->Update();

  ANTs::WriteImage<ImageType>(cropper->GetOutput(), argv[3]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ExtractRegionFromImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ExtractRegionFromImage");

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

  if (argc < 5 || argc > 7)
  {
    std::cerr << "Usage 1: " << argv[0] << " ImageDimension "
              << "inputImage outputImage minIndex maxIndex " << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " ImageDimension "
              << "inputImage outputImage label " << std::endl;
    std::cerr << "Usage 3: " << argv[0] << " ImageDimension "
              << "inputImage outputImage domainImage " << std::endl;
    std::cerr << "Usage 4: " << argv[0] << " ImageDimension "
              << "inputImage outputImage label labelImage 1 " << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return ExtractRegionFromImage<2>(argc, argv);
    }
    break;
    case 3:
    {
      return ExtractRegionFromImage<3>(argc, argv);
    }
    break;
    case 4:
    {
      return ExtractRegionFromImage<4>(argc, argv);
    }
    break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
