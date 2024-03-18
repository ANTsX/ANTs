#include "antsUtilities.h"

#include "ReadWriteData.h"

#include "itkCastImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include <string>

namespace ants
{
template <typename TPixel, unsigned int ImageDimension>
int
ConvertImage(int argc, char * argv[])
{
  using OutputPixelType = TPixel;

  if (argc > 4 && std::stoi(argv[4]) == 9)
  {
    using VectorType = itk::Vector<OutputPixelType, ImageDimension>;
    using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
    using ComponentImageType = itk::Image<OutputPixelType, ImageDimension>;

    typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();


    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      std::string filename = std::string(argv[2]);
      if (d == 0)
      {
        filename += std::string("xvec.nii.gz");
      }
      else if (d == 1)
      {
        filename += std::string("yvec.nii.gz");
      }
      else if (d == 2)
      {
        filename += std::string("zvec.nii.gz");
      }

      typename ComponentImageType::Pointer inputImage;
      ReadImage<ComponentImageType>(inputImage, filename.c_str());

      if (d == 0)
      {
        displacementField->CopyInformation(inputImage);
        displacementField->SetRegions(inputImage->GetRequestedRegion());
        displacementField->AllocateInitialized();
      }

      itk::ImageRegionConstIterator<ComponentImageType> It(inputImage, inputImage->GetLargestPossibleRegion());
      itk::ImageRegionIterator<DisplacementFieldType>   ItD(displacementField,
                                                          displacementField->GetLargestPossibleRegion());
      for (It.GoToBegin(), ItD.GoToBegin(); !It.IsAtEnd(); ++ItD, ++It)
      {
        VectorType V = ItD.Get();
        V[d] = It.Get();
        ItD.Set(V);
      }
    }
    ANTs::WriteImage<DisplacementFieldType>(displacementField, argv[3]);
  }
  else if (argc > 4 && std::stoi(argv[4]) == 10)
  {
    using VectorType = itk::Vector<OutputPixelType, ImageDimension>;
    using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
    using ComponentImageType = itk::Image<OutputPixelType, ImageDimension>;

    typename DisplacementFieldType::Pointer inputImage;
    ReadImage<DisplacementFieldType>(inputImage, argv[2]);

    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      using SelectorType = itk::VectorIndexSelectionCastImageFilter<DisplacementFieldType, ComponentImageType>;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput(inputImage);
      selector->SetIndex(d);
      selector->Update();

      std::string filename = std::string(argv[3]);
      if (d == 0)
      {
        filename += std::string("xvec.nii.gz");
      }
      else if (d == 1)
      {
        filename += std::string("yvec.nii.gz");
      }
      else if (d == 2)
      {
        filename += std::string("zvec.nii.gz");
      }
      ANTs::WriteImage<ComponentImageType>(selector->GetOutput(), filename.c_str());
    }
  }
  else if (argc > 4 && std::stoi(argv[4]) == 11)
  {
    using VectorType = itk::Vector<OutputPixelType, ImageDimension>;
    using VelocityFieldType = itk::Image<VectorType, ImageDimension + 1>;
    using ComponentImageType = itk::Image<OutputPixelType, ImageDimension + 1>;

    using ReaderType = itk::ImageFileReader<VelocityFieldType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[2]);
    reader->Update();

    typename VelocityFieldType::Pointer inputImage;
    ReadImage<VelocityFieldType>(inputImage, argv[2]);

    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      using SelectorType = itk::VectorIndexSelectionCastImageFilter<VelocityFieldType, ComponentImageType>;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput(reader->GetOutput());
      selector->SetIndex(d);
      selector->Update();

      std::string filename = std::string(argv[3]);
      if (d == 0)
      {
        filename += std::string("xvec.nii.gz");
      }
      else if (d == 1)
      {
        filename += std::string("yvec.nii.gz");
      }
      else if (d == 2)
      {
        filename += std::string("zvec.nii.gz");
      }
      ANTs::WriteImage<ComponentImageType>(selector->GetOutput(), filename.c_str());
    }
  }
  else if (argc > 4 && std::stoi(argv[4]) == 12)
  {
    using VectorType = itk::Vector<OutputPixelType, ImageDimension>;
    using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;

    typename DisplacementFieldType::Pointer displacementField;
    ReadImage<DisplacementFieldType>(displacementField, argv[2]);

    ANTs::WriteImage<DisplacementFieldType>(displacementField, argv[3]);
  }
  else if (argc == 4 || (argc > 4 && std::stoi(argv[4]) < 9))
  {
    using RealType = typename itk::NumericTraits<OutputPixelType>::RealType;
    using InputImageType = itk::Image<RealType, ImageDimension>;
    using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

    typename InputImageType::Pointer inputImage;
    ReadImage<InputImageType>(inputImage, argv[2]);

    std::vector<std::string> rescaleFileTypes;
    rescaleFileTypes.emplace_back(".png");
    rescaleFileTypes.emplace_back(".jpeg");
    rescaleFileTypes.emplace_back(".jpg");
    rescaleFileTypes.emplace_back(".tiff");
    rescaleFileTypes.emplace_back(".tif");
    rescaleFileTypes.emplace_back(".bmp");

    bool isRescaleType = false;
    for (auto & rescaleFileType : rescaleFileTypes)
    {
      if (strstr(argv[3], rescaleFileType.c_str()) != nullptr)
      {
        isRescaleType = true;
        break;
      }
    }

    if (isRescaleType)
    {
      using FilterType = itk::RescaleIntensityImageFilter<InputImageType, OutputImageType>;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput(inputImage);
      filter->SetOutputMinimum(itk::NumericTraits<OutputPixelType>::min());
      filter->SetOutputMaximum(itk::NumericTraits<OutputPixelType>::max());
      filter->Update();

      ANTs::WriteImage<OutputImageType>(filter->GetOutput(), argv[3]);
    }
    else
    {
      using CasterType = itk::CastImageFilter<InputImageType, OutputImageType>;
      typename CasterType::Pointer caster = CasterType::New();
      caster->SetInput(inputImage);
      caster->Update();

      ANTs::WriteImage<OutputImageType>(caster->GetOutput(), argv[3]);
    }
  }

  return EXIT_SUCCESS;
}

int
ConvertImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ConvertImage");

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

  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
              << "inputImage outputImage <pixelType>" << std::endl;
    std::cerr << "pixelType:  0 -> float (default)" << std::endl
              << "            1 -> unsigned char" << std::endl
              << "            2 -> unsigned short" << std::endl
              << "            3 -> unsigned int" << std::endl
              << "            4 -> unsigned long" << std::endl
              << "            5 -> char" << std::endl
              << "            6 -> short" << std::endl
              << "            7 -> int" << std::endl
              << "            8 -> long" << std::endl
              << "            9 -> component images to a float vector image" << std::endl
              << "           10 -> vector image to component images" << std::endl
              << "           11 -> time-varying velocity field image to component images (ImageDimension is the "
                 "dimensionality of the displacement vector)"
              << std::endl
              << "           12 -> float vector image" << std::endl;

    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      if (argc > 4 && std::stoi(argv[4]) == 1)
      {
        ConvertImage<unsigned char, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 2)
      {
        ConvertImage<unsigned short, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 3)
      {
        ConvertImage<unsigned int, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 4)
      {
        ConvertImage<unsigned long, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 5)
      {
        ConvertImage<char, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 6)
      {
        ConvertImage<short, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 7)
      {
        ConvertImage<int, 2>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 8)
      {
        ConvertImage<long, 2>(argc, argv);
      }
      else
      {
        ConvertImage<float, 2>(argc, argv);
      }
      break;
    case 3:
      if (argc > 4 && std::stoi(argv[4]) == 1)
      {
        ConvertImage<unsigned char, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 2)
      {
        ConvertImage<unsigned short, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 3)
      {
        ConvertImage<unsigned int, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 4)
      {
        ConvertImage<unsigned long, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 5)
      {
        ConvertImage<char, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 6)
      {
        ConvertImage<short, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 7)
      {
        ConvertImage<int, 3>(argc, argv);
      }
      else if (argc > 4 && std::stoi(argv[4]) == 8)
      {
        ConvertImage<long, 3>(argc, argv);
      }
      else
      {
        ConvertImage<float, 3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
