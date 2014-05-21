
#include "antsUtilities.h"
#include <algorithm>

#include <stdio.h>

#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <string>
#include <vector>

namespace ants
{

template <unsigned int ImageDimension>
int ExtractRegionFromImageByMask(int argc, char *argv[])
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>       ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  reader->Update();

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size;
  typename ImageType::RegionType::IndexType index;

  if( 0 )
    {
    std::vector<int> minIndex;
    std::vector<int> maxIndex;
    minIndex = ConvertVector<int>(std::string(argv[4]) );
    maxIndex = ConvertVector<int>(std::string(argv[5]) );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      index[i] = minIndex[i];
      size[i] = maxIndex[i] - minIndex[i] + 1;
      }
    region.SetSize(size);
    region.SetIndex(index);
    }
  else
    {
    typedef itk::Image<unsigned short, ImageDimension> ShortImageType;
//    typedef itk::CastImageFilter<ImageType, ShortImageType> CasterType;
//    typename CasterType::Pointer caster = CasterType::New();
//    caster->SetInput(reader->GetOutput());
//    caster->Update();

    typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
    typename ShortImageReaderType::Pointer shortReader = ShortImageReaderType::New();
    shortReader->SetFileName(argv[4]);
    shortReader->Update();

    // typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>
    typedef itk::LabelStatisticsImageFilter<ImageType, ShortImageType>
      StatsFilterType;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
//    stats->SetLabelInput(caster->GetOutput());
    stats->SetLabelInput(shortReader->GetOutput() );
//    stats->SetInput(caster->GetOutput());
    stats->SetInput(reader->GetOutput() );
    stats->Update();

    const unsigned int label = (argc >= 6) ? atoi(argv[5]) : 1;
    region = stats->GetRegion(label);

    std::cout << "bounding box of label=" << label
             << " : " << region << std::endl;

    const unsigned int padWidth = (argc >= 7) ? atoi(argv[6]) : 0;

    region.PadByRadius(padWidth);

    std::cout << "padding radius = " << padWidth
             << " : " << region << std::endl;

    region.Crop(reader->GetOutput()->GetBufferedRegion() );

    std::cout << "crop with original image region " << reader->GetOutput()->GetBufferedRegion()
             << " : " << region << std::endl;
    }

  std::cout << "final cropped region: " << region << std::endl;

  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput(reader->GetOutput() );
  cropper->SetExtractionRegion(region);
  cropper->SetDirectionCollapseToSubmatrix();
  cropper->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(cropper->GetOutput() );
  writer->SetFileName(argv[3]);
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ExtractRegionFromImageByMask( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ExtractRegionFromImageByMask" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );

  if( argc < 6 || argc > 7 )
    {
    std::cout << "Extract a sub-region from image using the bounding"
      " box from a label image, with optional padding radius."
             << std::endl << "Usage : " << argv[0] << " ImageDimension "
             << "inputImage outputImage labelMaskImage [label=1] [padRadius=0]"
             << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      ExtractRegionFromImageByMask<2>(argc, argv);
      }
      break;
    case 3:
      {
      ExtractRegionFromImageByMask<3>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
