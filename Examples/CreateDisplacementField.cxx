
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace ants
{
/* See usage output below */

template <unsigned int ImageDimension, unsigned int NumberOfComponents>
void
CreateDisplacementField(int argc, char * argv[])
{
  using ValueType = float;
  using ComponentImageType = itk::Image<ValueType, ImageDimension>;
  using VectorPixelType = itk::Vector<ValueType, NumberOfComponents>;
  using VectorImageType = itk::Image<VectorPixelType, ImageDimension>;
  using FileReaderType = itk::ImageFileReader<ComponentImageType>;

  typename FileReaderType::Pointer     reader = FileReaderType::New();
  typename ComponentImageType::Pointer componentImage;

  // EnforceZeroBoundaryFlag
  bool enforceZeroBoundaryFlag = static_cast<bool>(std::stoi(argv[2]));

  // Get image info
  reader->SetFileName(argv[3]);
  reader->Update();
  componentImage = reader->GetOutput();
  using RegionType = typename ComponentImageType::RegionType;
  RegionType regionOfFirstComponent = componentImage->GetLargestPossibleRegion();

  // Create output vector image
  typename VectorImageType::Pointer vectorImage = AllocImage<VectorImageType>(componentImage.GetPointer());
  vectorImage->SetRegions(componentImage->GetLargestPossibleRegion());

  itk::ImageRegionIteratorWithIndex<ComponentImageType> comIt(componentImage, regionOfFirstComponent);
  using VecItType = typename itk::ImageRegionIteratorWithIndex<VectorImageType>;
  using VecItIndexValueType = typename VecItType::IndexValueType;

  VecItType vecIt(vectorImage, regionOfFirstComponent);
  for (itk::SizeValueType n = 0; n < NumberOfComponents; n++)
  {
    reader->SetFileName(argv[3 + n]);
    reader->Update();
    componentImage = reader->GetOutput();
    if (componentImage->GetLargestPossibleRegion() != regionOfFirstComponent)
    {
      std::cout << "LargestPossibleRegion of component " << n << " does not match 1st component image." << std::endl;
      throw std::exception();
    }
    // Walk the images
    comIt.GoToBegin();
    vecIt.GoToBegin();
    for (; !comIt.IsAtEnd(); ++comIt, ++vecIt)
    {
      VectorPixelType vec = vecIt.Get();
      vec[n] = comIt.Get();
      vecIt.Set(vec);
    }
  }

  // Set zero boundary vectors
  if (enforceZeroBoundaryFlag)
  {
    VectorPixelType zeros(ImageDimension);
    zeros.Fill(itk::NumericTraits<ValueType>::ZeroValue());
    vecIt.GoToBegin();
    while (!vecIt.IsAtEnd())
    {
      bool isBoundary = false;
      for (itk::SizeValueType dim = 0; dim < ImageDimension; dim++)
      {
        if (vecIt.GetIndex()[dim] == 0 ||
            vecIt.GetIndex()[dim] == static_cast<VecItIndexValueType>(regionOfFirstComponent.GetSize()[dim] - 1))
        {
          isBoundary = true;
          break;
        }
      }
      if (isBoundary)
      {
        vecIt.Set(zeros);
      }
      ++vecIt;
    }
  }

  // Write out the output
  typename itk::ImageFileWriter<VectorImageType>::Pointer writer = itk::ImageFileWriter<VectorImageType>::New();
  writer->SetFileName(argv[argc - 1]);
  writer->SetInput(vectorImage);
  writer->Update();

  // debug
#if 0
  itk::ImageRegionIteratorWithIndex<VectorImageType> it( vectorImage, regionOfFirstComponent );
  it.GoToBegin();
  itk::OffsetValueType col = regionOfFirstComponent.GetSize()[0] / 2;
  itk::OffsetValueType row = regionOfFirstComponent.GetSize()[1] / 2;
  typename VectorImageType::IndexType index;
  index[0] = col;
  index[1] = 0;
  std::cout << "Debug output of field: " << std::endl;
  std::cout << "Middle column: " << std::endl;
  while( it.GetIndex()[1] < row * 2 )
    {
    it.SetIndex( index );
    std::cout << it.Get() << " ";
    index[1]++;
    }

  std::cout << std::endl;

  index[0] = 0;
  index[1] = row;
  std::cout << "Middle row: " << std::endl;
  while( it.GetIndex()[0] < col * 2 )
    {
    it.SetIndex( index );
    std::cout << it.Get() << " ";
    index[0]++;
    }

  std::cout << std::endl;
#endif // debug
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
CreateDisplacementField(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "CreateDisplacementField");

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
    std::cout << "Create an itkImage of itkVector pixels (NOT an itkVectorImage), using each scalar input component "
                 "image for each vector component. An itkImage of itkVectors is the standard type for displacement "
                 "fields in ITK. All component images (up to 8) are assumed to have the same size, offset, origin, and "
                 "spacing. The 'EnforceZeroBoundaryFlag' option will create zero-valued vectors along the borders when "
                 "enabled (pass 1), and is recommended for better displacement field behavior."
              << std::endl;
    std::cout << "Usage: " << argv[0]
              << " ImageDimension EnforceZeroBoundaryFlag{0/1} ComponentImage1 [ ComponentImage2 [...ComponentImageN] "
                 "] OutputImage "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }
  itk::SizeValueType imageDimension = std::stoi(argv[1]);

  itk::SizeValueType numberOfComponents = argc - 4;

  switch (imageDimension)
  {
    case 2:
    {
      switch (numberOfComponents)
      {
        case 1:
        {
          CreateDisplacementField<2, 1>(argc, argv);
        }
        break;
        case 2:
        {
          CreateDisplacementField<2, 2>(argc, argv);
        }
        break;
        case 3:
        {
          CreateDisplacementField<2, 3>(argc, argv);
        }
        break;
        case 4:
        {
          CreateDisplacementField<2, 4>(argc, argv);
        }
        break;
        case 5:
        {
          CreateDisplacementField<2, 5>(argc, argv);
        }
        break;
        case 6:
        {
          CreateDisplacementField<2, 6>(argc, argv);
        }
        break;
        case 7:
        {
          CreateDisplacementField<2, 7>(argc, argv);
        }
        break;
        case 8:
        {
          CreateDisplacementField<2, 8>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported number of components: " << numberOfComponents << std::endl;
          return EXIT_FAILURE;
      }
    }
    break;
    case 3:
    {
      switch (numberOfComponents)
      {
        case 1:
        {
          CreateDisplacementField<3, 1>(argc, argv);
        }
        break;
        case 2:
        {
          CreateDisplacementField<3, 2>(argc, argv);
        }
        break;
        case 3:
        {
          CreateDisplacementField<3, 3>(argc, argv);
        }
        break;
        case 4:
        {
          CreateDisplacementField<3, 4>(argc, argv);
        }
        break;
        case 5:
        {
          CreateDisplacementField<3, 5>(argc, argv);
        }
        break;
        case 6:
        {
          CreateDisplacementField<3, 6>(argc, argv);
        }
        break;
        case 7:
        {
          CreateDisplacementField<3, 7>(argc, argv);
        }
        break;
        case 8:
        {
          CreateDisplacementField<3, 8>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported number of components: " << numberOfComponents << std::endl;
          return EXIT_FAILURE;
      }
    }
    break;
    default:
      std::cout << "Unsupported number of dimensions: " << imageDimension << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
