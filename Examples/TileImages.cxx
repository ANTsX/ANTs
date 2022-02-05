#include "antsUtilities.h"
#include <algorithm>

#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"
#include "ReadWriteData.h"

#include <string>
#include <vector>

namespace ants
{

template <unsigned int ImageDimension>
int
TileImages(unsigned int argc, char * argv[])
{
  using PixelType = float;
  using ImageType = itk::Image<PixelType, ImageDimension>;

  using FilterType = itk::TileImageFilter<ImageType, ImageType>;
  typename FilterType::Pointer         filter = FilterType::New();
  typename FilterType::LayoutArrayType array;

  std::vector<unsigned int> layout = ConvertVector<unsigned int>(std::string(argv[3]));
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    array[d] = layout[d];
  }
  filter->SetLayout(array);
  for (unsigned int n = 4; n < argc; n++)
  {
    typename ImageType::Pointer inputImage;
    ReadImage<ImageType>(inputImage, argv[n]);

    filter->SetInput(n - 4, inputImage);
  }
  filter->Update();

  ANTs::WriteImage<ImageType>(filter->GetOutput(), argv[2]);

  return EXIT_SUCCESS;
}

int
CreateMosaic(unsigned int argc, char * argv[])
{
  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0] << " imageDimension outputImage layout inputImage1" << std::endl;
    return EXIT_FAILURE;
  }

  constexpr unsigned int ImageDimension = 3;

  using PixelType = float;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using SliceType = itk::Image<PixelType, ImageDimension - 1>;

  ImageType::Pointer inputImage;
  ReadImage<ImageType>(inputImage, argv[4]);

  std::vector<int> layout = ConvertVector<int>(std::string(argv[3]));
  if (layout.size() != 3)
  {
    std::cerr << "Layout for CreateMosaic is DxRxC where" << std::endl;
    std::cerr << "  D is direction, i.e. 0, 1, or 2.  If not any of those numbers, we pick the coarsest spacing."
              << std::endl;
    std::cerr << "  R is number of rows." << std::endl;
    std::cerr << "  C is number of cols." << std::endl;
    std::cerr << "  If R < 0 and C > 0 (or vice versa), the negative value is selected based on D" << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::SpacingType spacing = inputImage->GetSpacing();
  ImageType::SizeType    size = inputImage->GetRequestedRegion().GetSize();

  if (layout[0] < 0 || layout[0] > 2)
  {
    float        maxSpacing = spacing[0];
    unsigned int maxIndex = 0;
    for (unsigned int d = 1; d < ImageDimension; d++)
    {
      if (spacing[d] > static_cast<double>(maxSpacing))
      {
        maxSpacing = spacing[d];
        maxIndex = d;
      }
    }
    layout[0] = maxIndex;
  }

  unsigned long numberOfSlices = size[layout[0]];

  int numberOfRows = std::min(static_cast<int>(layout[1]), static_cast<int>(numberOfSlices));
  int numberOfColumns = std::min(static_cast<int>(layout[2]), static_cast<int>(numberOfSlices));

  if (numberOfRows <= 0 && numberOfColumns > 0)
  {
    numberOfRows = std::ceil(static_cast<float>(numberOfSlices) / static_cast<float>(numberOfColumns));
  }
  else if (numberOfColumns <= 0 && numberOfRows > 0)
  {
    numberOfColumns = std::ceil(static_cast<float>(numberOfSlices) / static_cast<float>(numberOfRows));
  }
  else if (numberOfColumns <= 0 && numberOfRows <= 0)
  {
    numberOfRows = static_cast<int>(std::sqrt(static_cast<float>(numberOfSlices)));
    numberOfColumns = std::ceil(static_cast<float>(numberOfSlices) / static_cast<float>(numberOfRows));
  }

  std::cout << "Slices[" << layout[0] << "]: " << numberOfSlices << std::endl;
  std::cout << "Rows:  " << numberOfRows << std::endl;
  std::cout << "Columns:  " << numberOfColumns << std::endl;

  using FilterType = itk::TileImageFilter<SliceType, SliceType>;
  FilterType::LayoutArrayType array;

  array[0] = numberOfColumns;
  array[1] = numberOfRows;

  ImageType::RegionType region;
  size[layout[0]] = 0;

  FilterType::Pointer filter = FilterType::New();
  filter->SetLayout(array);

  for (unsigned int n = 0; n < numberOfSlices; n++)
  {
    ImageType::IndexType index;
    index.Fill(0);
    index[layout[0]] = static_cast<int>(n);
    region.SetIndex(index);
    region.SetSize(size);

    using ExtracterType = itk::ExtractImageFilter<ImageType, SliceType>;
    ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput(inputImage);
    extracter->SetExtractionRegion(region);
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    filter->SetInput(n, extracter->GetOutput());
  }
  filter->Update();

  ANTs::WriteImage<SliceType>(filter->GetOutput(), argv[2]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
TileImages(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "TileImages");

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
    std::cout << argv[0] << " imageDimension outputImage layout inputImage1 ... inputImageN" << std::endl;
    std::cout << "  The layout has the same dimension as the output image. If all entries of " << std::endl;
    std::cout << "  the layout are positive, the tiled output will contain the exact number  " << std::endl;
    std::cout << "  of tiles. If the layout contains a 0 in the last dimension, the filter " << std::endl;
    std::cout << "  will compute a size that will accommodate all of the images. " << std::endl;
    std::cout << "  The input images must have a dimension less than or equal to the output " << std::endl;
    std::cout << "  image. The output image could have a larger dimension than the input. " << std::endl;
    std::cout << "  For example, This filter can be used to create a 3-d volume from a series " << std::endl;
    std::cout << "  of 2-d inputs by specifying a layout of 1x1x0. " << std::endl << std::endl;

    std::cout << "  In addition to the above functionality, there is another usage option" << std::endl;
    std::cout << "  for creating a 2-d tiled mosaic from a 3-D image.  The command line options" << std::endl;
    std::cout << "  are the same except only 1 input is expected and the layout for this option" << std::endl;
    std::cout << "  is DxRxC where:" << std::endl;
    std::cout << "      D is direction, i.e. 0, 1, or 2.  If not any of those numbers, we pick the coarsest spacing."
              << std::endl;
    std::cout << "      R is number of rows." << std::endl;
    std::cout << "      C is number of cols." << std::endl;
    std::cout << "      If R < 0 and C > 0 (or vice versa), the negative value is selected based on D" << std::endl;

    // Should add the following options:
    //    * add rgb overlay (with alpha value?)
    //    * number of slices to skip
    //    * beginning and ending slice
    //    * add or subtract border around each slice/tile
    //    * if adding, set pad constant value


    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  const int ImageDimension = static_cast<int>(std::stoi(argv[1]));

  if (ImageDimension == 3 && argc == 5)
  {
    CreateMosaic(argc, argv);
  }
  else
  {
    switch (ImageDimension)
    {
      case 2:
      {
        return TileImages<2>(argc, argv);
      }
      break;
      case 3:
      {
        return TileImages<3>(argc, argv);
      }
      break;
      case 4:
      {
        return TileImages<4>(argc, argv);
      }
      break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
} // namespace ants
