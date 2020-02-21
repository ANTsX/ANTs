
#include "antsUtilities.h"
#include <algorithm>
#include <algorithm>
#include "stdio.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "ReadWriteData.h"
#include "TensorFunctions.h"

namespace ants
{
template <unsigned int ImageDimension>
int AverageTensorImages(unsigned int argc, char *argv[])
{
  // typedef itk::Vector<float,6> TensorType;
  typedef itk::SymmetricSecondRankTensor<float, 3> TensorType;

  typedef itk::Image<TensorType, ImageDimension>       ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

  char * outputName = argv[2];
  int    mathtype = std::stoi(argv[3]);
  float  numberofimages = static_cast<float>( argc ) - 4.0f;

  std::cout << "Averaging " << numberofimages << " images " << std::endl;

  typename ImageType::Pointer averageimage = nullptr;
  typename ImageType::Pointer image2 = nullptr;

  typename ImageType::SizeType size;
  size.Fill(0);
  unsigned int bigimage = 0;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    std::cout << " fn " << fn << std::endl;
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        bigimage = j;
        std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  std::cout << " largest image " << size << std::endl;

  bool logeuc = true;
  if( mathtype == 1 )
    {
    logeuc = false;
    }

  TensorType nullTensor;
  nullTensor[0] = nullTensor[1] = nullTensor[2] = nullTensor[3]
          = nullTensor[4] = nullTensor[5] = 0;

  ReadTensorImage<ImageType>(averageimage, argv[bigimage], logeuc);
  averageimage->FillBuffer(nullTensor);
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::string fn = std::string(argv[j]);
    ReadTensorImage<ImageType>(image2, fn.c_str(), logeuc);

    IteratorType vfIter( image2,  image2->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      TensorType val =  vfIter.Get() / numberofimages;
      averageimage->SetPixel(vfIter.GetIndex(), val + averageimage->GetPixel(vfIter.GetIndex() ) );
      }
    }

  WriteTensorImage<ImageType>(averageimage, outputName, logeuc);

  return EXIT_SUCCESS;
}

// Main Program

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int AverageTensorImages( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "AverageTensorImages" );
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
  argv[argc] = nullptr;
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

  try
    {
    if( argc - 4 < 1 )
      {
      std::cerr << "Basic useage ex: " << std::endl;
      std::cerr << argv[0] << " ImageDimension  average.nii mathtype list-of-files-via-wildcard " << std::endl;
      std::cerr << " e.g. \n   AverageTensorImages 3  average.nii  1  *registered.nii " << std::endl;
      std::cerr << " mathtype=[0=log-euclidean, 1=euclidean] " << std::endl;
      if( argc >= 2 &&
          ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
        return EXIT_SUCCESS;
        }
      return EXIT_FAILURE;
      }

    int dim = std::stoi(argv[1]);
    // char * outputName = argv[2];
    // int mathtype = std::stoi(argv[3]);
    // int numberofimages = argc - 4;

    // Get the image dimension
    switch( dim )
      {
      case 2:
        {
        return AverageTensorImages<2>(argc, argv);
        }
        break;
      case 3:
        {
        return AverageTensorImages<3>(argc, argv);
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }

    return EXIT_SUCCESS;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
}
} // namespace ants
