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
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension>
int ConvertType(int argc, char *argv[])
{
  typedef  unsigned char                             outPixelType;
  typedef  float                                     floatPixelType;
  typedef  float                                     inPixelType;
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
  std::cout << " Updated reader " << std::endl;

  typedef itk::CastImageFilter<ImageType, IntermediateType> castertype;
  typename   castertype::Pointer caster = castertype::New();
  caster->SetInput(reader->GetOutput() );
  caster->Update();

  // Rescale the image intensities so that they fall between 0 and 255
  typedef itk::RescaleIntensityImageFilter<IntermediateType, IntermediateType> FilterType;
  typename   FilterType::Pointer fixedrescalefilter = FilterType::New();
  fixedrescalefilter->SetInput(caster->GetOutput() );
  constexpr double desiredMinimum = 0.0;
  constexpr double desiredMaximum = 255.0;
  fixedrescalefilter->SetOutputMinimum( desiredMinimum );
  fixedrescalefilter->SetOutputMaximum( desiredMaximum );
  fixedrescalefilter->UpdateLargestPossibleRegion();

  typedef itk::CastImageFilter<IntermediateType, OutImageType> castertype2;
  typename castertype2::Pointer caster2 = castertype2::New();
  caster2->SetInput(fixedrescalefilter->GetOutput() );
  caster2->Update();

  typename   OutImageType::Pointer outim = caster2->GetOutput();
  typename   OutImageType::SpacingType spc = outim->GetSpacing();
  outim->SetSpacing(spc);
  std::cout << " Dire in " << reader->GetOutput()->GetDirection() << std::endl;
  std::cout << " Dire out " << outim->GetDirection() << std::endl;
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
int ConvertToJpg( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  try
    {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert( args.begin(), "ConvertToJpg" );
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
      std::cout << "Usage:   ConvertToJpg infile.nii out.jpg " << std::endl;
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
        fn.c_str(), itk::ImageIOFactory::FileModeType::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();

    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2>(argc, argv);
        }
        break;
      case 3:
        {
        ConvertType<3>(argc, argv);
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  catch( const itk::ExceptionObject & e )
    {
    e.Print( std::cerr );
    return EXIT_FAILURE;
    }
  catch( const std::exception & e )
    {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
