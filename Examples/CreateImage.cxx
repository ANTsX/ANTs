
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include <cstdio>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include <string>
#include <vector>

namespace ants
{

template <int ImageDimension>
int
CreateZeroImage(int argc, char * argv[])
{
  using PixelType = float;
  using ImageType = itk::Image<PixelType, ImageDimension>;

  using GeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
  typename GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();
  generator->SetSeed();

  std::string                     which = std::string(argv[3]);
  typename std::string::size_type pos = which.find(".");

  typename std::string::size_type pos3 = std::string::npos;
  if (argc > 6)
  {
    std::string pixelValues = std::string(argv[6]);
    pos3 = pixelValues.find("x");
  }
  typename ImageType::Pointer image;
  if (pos3 != std::string::npos)
  {
    std::vector<float> og = ConvertVector<float>(std::string(argv[3]));
    std::vector<float> sp = ConvertVector<float>(std::string(argv[4]));
    std::vector<int>   sz = ConvertVector<int>(std::string(argv[5]));
    std::vector<float> values = ConvertVector<float>(std::string(argv[6]));

    unsigned long numberOfPixels = 1;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      numberOfPixels *= sz[d];
    }
    if (values.size() > numberOfPixels)
    {
      std::cout << "Number of specified pixel values is greater than "
                << "the size of the image." << std::endl;
      return EXIT_FAILURE;
    }

    if (og.size() != ImageDimension)
    {
      std::cout << "Invalid origin size." << std::endl;
      return EXIT_FAILURE;
    }
    if (sp.size() != ImageDimension)
    {
      std::cout << "Invalid spacing size." << std::endl;
      return EXIT_FAILURE;
    }
    if (sz.size() != ImageDimension)
    {
      std::cout << "Invalid Size size." << std::endl;
      return EXIT_FAILURE;
    }

    typename ImageType::PointType   origin;
    typename ImageType::SpacingType spacing;
    typename ImageType::SizeType    size;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      origin[d] = og[d];
      spacing[d] = sp[d];
      size[d] = sz[d];
    }
    typename ImageType::RegionType region;
    region.SetSize(size);
    typename ImageType::DirectionType direction;
    direction.SetIdentity();
    image = AllocImage<ImageType>(region, spacing, origin, direction, 0.0);

    unsigned long                       count = 0;
    itk::ImageRegionIterator<ImageType> It(image, image->GetRequestedRegion());

    It.GoToBegin();
    while (!It.IsAtEnd() && count < values.size())
    {
      It.Set(values[count]);

      ++It;
      ++count;
    }

    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(image);
    writer->Update();
  }
  else if (pos != std::string::npos)
  {
    using ReaderType = itk::ImageFileReader<ImageType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[2]);
    reader->Update();
    // ORIENTATION ALERT  -- the original code here
    // set the region, spacing, and origin without setting directions.
    image = AllocImage<ImageType>(reader->GetOutput(), atof(argv[4]));

    if (argc > 5)
    {
      switch (std::stoi(argv[5]))
      {
        case 1:
        default:
        {
          itk::ImageRegionIterator<ImageType> It(image, image->GetLargestPossibleRegion());
          for (It.GoToBegin(); !It.IsAtEnd(); ++It)
          {
            It.Set(static_cast<PixelType>(generator->GetIntegerVariate(static_cast<int>(It.Get()))));
          }
          break;
        }
          //        case 2:
          //         {
          //                                    itk::ImageRegionIteratorWithIndex<ImageType> ItI( image,
          //                                            image->GetLargestPossibleRegion() );
          //                                    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
          //                                            {
          //                                            ItI.Set( constant - ItI.GetIndex()[d] );
          //                                            }
          //          break;
          //          }
          //        default:
          //          std::cout << "Incorrect choice" << std::endl;
          //          return EXIT_FAILURE;
          //          break;
      }
    }

    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[3]);
    writer->SetInput(image);
    writer->Update();
  }
  else
  {
    std::vector<float> og = ConvertVector<float>(std::string(argv[3]));
    std::vector<float> sp = ConvertVector<float>(std::string(argv[4]));
    std::vector<int>   sz = ConvertVector<int>(std::string(argv[5]));

    if (og.size() != ImageDimension)
    {
      std::cout << "Invalid origin size." << std::endl;
      return EXIT_FAILURE;
    }
    if (sp.size() != ImageDimension)
    {
      std::cout << "Invalid spacing size." << std::endl;
      return EXIT_FAILURE;
    }
    if (sz.size() != ImageDimension)
    {
      std::cout << "Invalid Size size." << std::endl;
      return EXIT_FAILURE;
    }

    typename ImageType::PointType     origin;
    typename ImageType::SpacingType   spacing;
    typename ImageType::SizeType      size;
    typename ImageType::DirectionType direction;
    direction.SetIdentity();
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      origin[d] = og[d];
      spacing[d] = sp[d];
      size[d] = sz[d];
    }
    typename ImageType::RegionType region;
    region.SetSize(size);
    image = AllocImage<ImageType>(region, spacing, origin, direction, atof(argv[6]));
    if (argc > 7)
    {
      switch (std::stoi(argv[7]))
      {
        case 1:
        default:
        {
          itk::ImageRegionIterator<ImageType> It(image, image->GetLargestPossibleRegion());
          for (It.GoToBegin(); !It.IsAtEnd(); ++It)
          {
            It.Set(static_cast<PixelType>(generator->GetIntegerVariate(static_cast<int>(It.Get()))));
          }
          break;
        }
          //        case 2:
          //         {
          //                                    itk::ImageRegionIteratorWithIndex<ImageType> ItI( image,
          //                                            image->GetLargestPossibleRegion() );
          //                                    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
          //                                            {
          //                                            ItI.Set( constant - ItI.GetIndex()[d] );
          //                                            }
          //          break;
          //          }
          //        default:
          //          {
          //          std::cout << "Incorrect choice" << std::endl;
          //          return EXIT_FAILURE;
          //          break;
          //          }
      }
    }

    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(image);
    writer->Update();
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
CreateImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "CreateImage");

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
    std::cout << "Usage 1: " << argv[0] << " imageDimension referenceImage outputImage constant [random?]" << std::endl;
    std::cout << "Usage 2: " << argv[0] << " imageDimension outputImage origin spacing size constant [random?]"
              << std::endl;
    std::cout << "Usage 3: " << argv[0] << " imageDimension outputImage origin spacing size pixelValues" << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 1:
    {
      return CreateZeroImage<1>(argc, argv);
    }
    break;
    case 2:
    {
      return CreateZeroImage<2>(argc, argv);
    }
    break;
    case 3:
    {
      return CreateZeroImage<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
