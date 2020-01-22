/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension, unsigned int NVectorComponents>
int MultiplyImages(int argc, char *argv[])
{
  typedef itk::Vector<float, NVectorComponents>        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;

  if( argc < 3 )
    {
    std::cout << "missing 1st filename" << std::endl;
    throw;
    }
  if( argc < 4 )
    {
    std::cout << "missing 2nd filename" << std::endl;
    throw;
    }
  if( argc < 5 )
    {
    std::cout << "missing output filename" << std::endl;
    throw;
    }

  std::string fn1 = std::string(argv[2]);
  std::string fn2 = std::string(argv[3]);

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer varimage = nullptr;

  typename readertype::Pointer reader2 = readertype::New();
  typename readertype::Pointer reader1 = readertype::New();
  reader2->SetFileName(fn2.c_str() );

  bool isfloat = false;
  try
    {
    reader2->UpdateLargestPossibleRegion();
    }
  catch( ... )
    {
//    std::cout << " Rather than opening " << fn2
//             <<
//      " as an image file, this program has decided, in its great wisdom, to consider it to be a floating point numerical value, and has acted accordingly -- i.e. read this as a number. "
//             << std::endl;
    isfloat = true;
    }

  float floatval = 1.0;
  if( isfloat )
    {
    floatval = atof(argv[3]);
    }
  else
    {
    image2 = reader2->GetOutput();
    }

  reader1->SetFileName(fn1.c_str() );
  try
    {
    reader1->UpdateLargestPossibleRegion();
    image1 = reader1->GetOutput();
    }
  catch( ... )
    {
    std::cout << " read 1 error ";
    }

  varimage = AllocImage<ImageType>(image1);

  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    typename ImageType::IndexType ind = vfIter2.GetIndex();
    PixelType pix1 = image1->GetPixel(ind);
    if( isfloat )
      {
      vfIter2.Set(pix1 * floatval);
      }
    else
      {
      vfIter2.Set(pix1 * image2->GetPixel(ind) );
      }
    }
  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[4]);
  writer->SetInput( varimage );
  writer->Write();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int MultiplyImages( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "MultiplyImages" );

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

  if( argc < 4 )
    {
    std::cout << "Usage:  " << std::endl;
    std::cout << argv[0] << " ImageDimension img1.nii img2.nii product.nii {smoothing}" << std::endl;
    std::cout <<     " 2nd image file may also be floating point numerical value, and program will act accordingly -- i.e. read as a number. " << std::endl;
    std::cout << " Program handles vector and tensor images as well " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  int                       dim = std::stoi( argv[1] );
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(argv[2], itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(argv[2]);
  imageIO->ReadImageInformation();
  unsigned int ncomponents = imageIO->GetNumberOfComponents();

  int returnValue = EXIT_FAILURE;

  // Get the image dimension
  switch( std::stoi(argv[1]) )
    {
    case 1:
      {
      switch( ncomponents )
        {
        case 3:
          {
          returnValue = MultiplyImages<1, 3>(argc, argv);
          }
          break;
        case 2:
          {
          returnValue = MultiplyImages<1, 2>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MultiplyImages<1, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 2:
      {
      switch( ncomponents )
        {
        case 3:
          {
          returnValue = MultiplyImages<2, 3>(argc, argv);
          }
          break;
        case 2:
          {
          returnValue = MultiplyImages<2, 2>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MultiplyImages<2, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 3:
      {
      switch( ncomponents )
        {
        case 7:
          {
          returnValue = MultiplyImages<3, 7>(argc, argv);
          }
          break;
        case 6:
          {
          returnValue = MultiplyImages<3, 6>(argc, argv);
          }
          break;
        case 3:
          {
          returnValue = MultiplyImages<3, 3>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MultiplyImages<3, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 4:
      {
      switch( ncomponents )
        {
        case 7:
          {
          returnValue = MultiplyImages<4, 7>(argc, argv);
          }
          break;
        case 6:
          {
          returnValue = MultiplyImages<4, 6>(argc, argv);
          }
          break;
        case 4:
          {
          returnValue = MultiplyImages<4, 4>(argc, argv);
          }
          break;
        case 3:
          {
          returnValue = MultiplyImages<4, 3>(argc, argv);
          }
          break;
        case 2:
          {
          returnValue = MultiplyImages<4, 2>(argc, argv);
          }
          break;
        default:
          {
          returnValue = MultiplyImages<4, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    default:
      std::cout << " not supported " << dim  << std::endl;
      return EXIT_FAILURE;
    }

  return returnValue;
}
} // namespace ants
