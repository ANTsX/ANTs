#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMatrixOffsetTransformBase.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkGridImageSource.h"

namespace ants
{

template <unsigned int ImageDimension>
int
CreateWarpedGridImage(int argc, char * argv[])
{
  using RealType = float;
  using RealImageType = itk::Image<RealType, ImageDimension>;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using VectorImageType = itk::Image<VectorType, ImageDimension>;

  /**
   * Read in vector field
   */
  using ReaderType = itk::ImageFileReader<VectorImageType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  reader->Update();

  using GridSourceType = itk::GridImageSource<RealImageType>;
  typename GridSourceType::Pointer gridder = GridSourceType::New();
  gridder->SetSpacing(reader->GetOutput()->GetSpacing());
  gridder->SetOrigin(reader->GetOutput()->GetOrigin());
  gridder->SetSize(reader->GetOutput()->GetLargestPossibleRegion().GetSize());

  using ArrayType = typename GridSourceType::ArrayType;
  using BoolArrayType = typename GridSourceType::BoolArrayType;

  ArrayType     gridSpacing;
  ArrayType     gridSigma;
  BoolArrayType which;

  which.Fill(false);
  for (itk::SizeValueType i = 0; i < 2; i++)
  {
    which[i] = true;
  }

  if (argc > 4)
  {
    std::vector<unsigned int> directions = ConvertVector<unsigned int>(std::string(argv[4]));
    if (directions.size() != ImageDimension)
    {
      std::cout << "Incorrect direction size." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      for (itk::SizeValueType i = 0; i < ImageDimension; i++)
      {
        which[i] = static_cast<bool>(directions[i]);
      }
    }
  }
  for (itk::SizeValueType i = 0; i < ImageDimension; i++)
  {
    gridSpacing[i] =
      static_cast<typename ArrayType::ValueType>(reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i]) *
      static_cast<typename ArrayType::ValueType>(reader->GetOutput()->GetSpacing()[i]) /
      static_cast<typename ArrayType::ValueType>(25.0);
    gridSigma[i] = gridSpacing[i] / static_cast<typename ArrayType::ValueType>(10.0);
  }
  if (argc > 5)
  {
    std::vector<RealType> spacing = ConvertVector<RealType>(std::string(argv[5]));
    if (spacing.size() != ImageDimension)
    {
      std::cout << "Incorrect spacing size." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      for (itk::SizeValueType i = 0; i < ImageDimension; i++)
      {
        gridSpacing[i] = spacing[i];
        gridSigma[i] = gridSpacing[i] / 10.0;
      }
    }
  }
  if (argc > 6)
  {
    std::vector<typename ArrayType::ValueType> sigma =
      ConvertVector<typename ArrayType::ValueType>(std::string(argv[6]));
    if (sigma.size() != ImageDimension)
    {
      std::cout << "Incorrect sigma size." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      for (itk::SizeValueType i = 0; i < ImageDimension; i++)
      {
        gridSigma[i] = sigma[i] / static_cast<typename ArrayType::ValueType>(10.0);
      }
    }
  }

  gridder->SetGridSpacing(gridSpacing);
  gridder->SetSigma(gridSigma);
  gridder->SetWhichDimensions(which);
  gridder->Update();
  typename RealImageType::Pointer grid = gridder->GetOutput();
  grid->SetDirection(reader->GetOutput()->GetDirection());
  grid->SetOrigin(reader->GetOutput()->GetOrigin());
  grid->SetSpacing(reader->GetOutput()->GetSpacing());

  using TransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  using WarperType = itk::WarpImageMultiTransformFilter<RealImageType, RealImageType, VectorImageType, TransformType>;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput(grid);
  warper->SetEdgePaddingValue(0);
  warper->SetSmoothScale(1);
  warper->PushBackDisplacementFieldTransform(reader->GetOutput());
  warper->SetOutputParametersFromImage(reader->GetOutput());
  warper->Update();

  std::string file = std::string(argv[3]);
  using ImageWriterType = itk::ImageFileWriter<RealImageType>;
  typename ImageWriterType::Pointer gridWriter = ImageWriterType::New();
  gridWriter->SetFileName(file.c_str());
  gridWriter->SetInput(warper->GetOutput());
  gridWriter->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
CreateWarpedGridImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "CreateWarpedGridImage");

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
    std::cout << "Usage: " << argv[0] << " ImageDimension deformationField "
              << "outputImage [directions, e.g. 1x0x0] [gridSpacing, e.g. 10x10x10] [gridSigma, e.g. 1x1x1]"
              << std::endl;
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
      return CreateWarpedGridImage<2>(argc, argv);
    }
    break;
    case 3:
    {
      return CreateWarpedGridImage<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
