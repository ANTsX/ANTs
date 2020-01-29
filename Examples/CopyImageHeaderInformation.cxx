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
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cstdio>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "ReadWriteData.h"
#include "TensorFunctions.h"

namespace ants
{
template <unsigned int ImageDimension>
int CopyImageHeaderInformation(int argc, char *argv[])
{
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;

  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  //  std::cout << " Spacing " << reader->GetOutput()->GetSpacing() << std::endl;
  // std::cout << " Origin " << reader->GetOutput()->GetOrigin() << std::endl;
  // std::cout << " Direction " << std::endl << reader->GetOutput()->GetDirection() << std::endl;
  // std::cout << " Size " << std::endl << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

  bool istensor = false;
  if( argc > 7 )
    {
    if( std::stoi(argv[7]) )
      {
      istensor = true;
      }
    }
  if( istensor )
    {
    typedef itk::Vector<float, 6>                  TensorType;
    typedef itk::Image<TensorType, ImageDimension> TensorFieldType;
    typename TensorFieldType::Pointer timage;
    ReadTensorImage<TensorFieldType>(timage, argv[2], false);
    //      std::cout<< " tim dir " << timage->GetDirection() << std::endl;
    if( argc > 6 )
      {
      if( std::stoi(argv[6]) )
        {
        timage->SetSpacing(  reader->GetOutput()->GetSpacing()  );
        }
      }
    if( argc > 5 )
      {
      if( std::stoi(argv[5]) )
        {
        timage->SetOrigin(  reader->GetOutput()->GetOrigin()  );
        }
      }
    if( argc > 4 )
      {
      if( std::stoi(argv[4]) )
        {
        timage->SetDirection(  reader->GetOutput()->GetDirection()  );
        }
      }

    //      std::cout<< " tim dir " << timage->GetDirection() << std::endl;
    WriteTensorImage<TensorFieldType>( timage, argv[3], false);

    return EXIT_SUCCESS;
    }

  typename readertype::Pointer reader2 = readertype::New();
  reader2->SetFileName(argv[2]);
  reader2->Update();

  // MakeNewImage(typename TImage::Pointer image1, typename TImage::PixelType initval)
  typename ImageType::Pointer newimage = MakeNewImage<ImageType>(reader2->GetOutput(), -1);

  if( argc > 6 )
    {
    if( std::stoi(argv[6]) )
      {
      newimage->SetSpacing(  reader->GetOutput()->GetSpacing()  );
      }
    }
  if( argc > 5 )
    {
    if( std::stoi(argv[5]) )
      {
      newimage->SetOrigin(  reader->GetOutput()->GetOrigin()  );
      }
    }
  if( argc > 4 )
    {
    if( std::stoi(argv[4]) )
      {
      newimage->SetDirection(  reader->GetOutput()->GetDirection()  );
      }
    }

  WriteImage<ImageType>(newimage, argv[3]);

  return EXIT_FAILURE;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int CopyImageHeaderInformation( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CopyImageHeaderInformation" );

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

  if( argc < 4  )
    {
    std::cout << "Usage:  " << argv[0]
             <<
      " refimage.ext imagetocopyrefimageinfoto.ext imageout.ext   boolcopydirection  boolcopyorigin boolcopyspacing  {bool-Image2-IsTensor}"
             << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim = imageIO->GetNumberOfDimensions();

  switch( dim  )
    {
    case 2:
      {
      CopyImageHeaderInformation<2>(argc, argv);
      }
      break;
    case 3:
      {
      CopyImageHeaderInformation<3>(argc, argv);
      }
      break;
    case 4:
      {
      CopyImageHeaderInformation<4>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension : " << dim << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
