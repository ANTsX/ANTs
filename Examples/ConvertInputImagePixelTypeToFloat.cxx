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
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <climits>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension, typename TPIXELTYPE>
int ConvertTypeToFloat( int argc, char *argv[] )
{
  typedef  TPIXELTYPE                                inPixelType;
  typedef  float                                     floatPixelType;
  typedef  float                                     outPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::Image<floatPixelType, ImageDimension> IntermediateType;
  typedef itk::Image<outPixelType, ImageDimension>   OutImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;
  typedef itk::ImageFileWriter<OutImageType>         writertype;

  typename readertype::Pointer reader = readertype::New();
  if( argc < 2 )
    {
    std::cout << "Missing input filename" << std::endl;
    throw;
    }
  reader->SetFileName(argv[1]);
  reader->Update();

  typedef itk::CastImageFilter<ImageType, IntermediateType> castertype;
  typename   castertype::Pointer caster = castertype::New();
  caster->SetInput(reader->GetOutput() );
  caster->Update();
  typedef itk::CastImageFilter<IntermediateType, OutImageType> castertype2;
  typename castertype2::Pointer caster2 = castertype2::New();
  caster2->SetInput(caster->GetOutput() );
  caster2->Update();

  typename   OutImageType::Pointer outim = caster2->GetOutput();
  typename   writertype::Pointer writer = writertype::New();
  if( argc < 3 )
    {
    std::cout << "Missing output filename" << std::endl;
    throw;
    }
  writer->SetFileName(argv[2]);
  writer->SetInput(outim);
  writer->Update();
  writer->Write();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ConvertInputImagePixelTypeToFloat( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ConvertInputImagePixelTypeToFloat" );
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

  if( argc < 3 )
    {
    std::cout << "Usage:   " << argv[0] << " infile.nii out.ext TYPE-OPTION " << std::endl;
    std::cout << " ext is the extension you want, e.g. tif.  " << std::endl;
    std::cout << " TYPE-OPTION  :  TYPE " << std::endl;
    std::cout << "  0  :  char   " << std::endl;
    std::cout << "  1  :  unsigned char   " << std::endl;
    std::cout << "  2  :  short   " << std::endl;
    std::cout << "  3  :  unsigned short   " << std::endl;
    std::cout << "  4  :  int   " << std::endl;
    std::cout << "  5  :  unsigned int   " << std::endl;
    std::cout
      << " Note that some pixel types are not supported by some image formats. e.g.  int is not supported by jpg. "
      << std::endl;
    std::cout
      << " You can easily extend this for other pixel types with a few lines of code and adding usage info. "
      << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }
  unsigned int typeoption = 0;
  if( argc > 3 )
    {
    typeoption = std::stoi(argv[3]);
    }
  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();

  if( typeoption == 0 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, char>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, char>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, char>(argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  else if( typeoption == 1 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, unsigned char>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, unsigned char>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, unsigned char>(argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  else if( typeoption == 2 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, short>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, short>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, short>(argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  else if( typeoption == 3 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, unsigned short>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, unsigned short>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, unsigned short>(argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  else if( typeoption == 4 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, int>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, int>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, int>(argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  else if( typeoption == 5 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertTypeToFloat<2, unsigned int>(argc, argv );
        }
        break;
      case 3:
        {
        ConvertTypeToFloat<3, unsigned int>(argc, argv );
        }
        break;
      case 4:
        {
        ConvertTypeToFloat<4, unsigned int>(argc, argv );
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
