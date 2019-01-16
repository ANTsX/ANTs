/*
 * SetSpacing.cxx
 *
 *  Created on: Mar 5, 2009
 *      Author: songgang
 */

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
int SetSpacing(int argc, char *argv[])
{
  typedef  float                                     outPixelType;
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::Image<outPixelType, ImageDimension>   OutImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;
  typedef itk::ImageFileWriter<OutImageType>         writertype;

  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  typename OutImageType::Pointer outim = reader->GetOutput();
  typename OutImageType::SpacingType spacing = outim->GetSpacing();

  std::cout << " Old Spacing " <<  outim->GetSpacing();
  if( argc > 3 )
    {
    spacing[0] = atof(argv[3]);
    }
  if( argc > 4 )
    {
    spacing[1] = atof(argv[4]);
    }
  if( argc > 5 )
    {
    spacing[2] = atof(argv[5]);
    }
  std::cout << "  New Spacing " << spacing << std::endl;

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typename ImageType::Pointer varimage = AllocImage<ImageType>(outim);
  varimage->SetSpacing(spacing);

  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    vfIter2.Set(outim->GetPixel(vfIter2.GetIndex() ) );
    }

  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(  varimage );
  writer->Update();
  writer->Write();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SetSpacing( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SetSpacing" );

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
    std::cout << "Usage:   " << argv[0] << "  Dimension infile.hdr outfile.nii  SpacingX SpacingY {SpacingZ} "
             << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension
  switch( std::stoi(argv[1]) )
    {
    case 2:
      {
      return SetSpacing<2>(argc - 1, argv + 1);
      }
      break;
    case 3:
      {
      return SetSpacing<3>(argc - 1, argv + 1);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // namespace ants
