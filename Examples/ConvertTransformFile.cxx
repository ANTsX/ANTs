/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ConvertTransformFile.cxx,v $
  Language:  C++
  Date:      $Date: 2012/08/02 20:09:47 $
  Version:   $Revision: 1.4 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include <iostream>
#include <sys/stat.h>

#include "itksys/SystemTools.hxx"

#include <fstream>
#include <stdio.h>

#include "itkantsReadWriteTransform.h"

/* Utility to read in a transform file (presumed to be in binary format) and output
 * it in legacy text format for human reading.
 * The option '--matrix' will instead output only the transform matrix to a text file,
 * one row per line with space-delimited values. This option works only for
 * transform of MatrixOffsetTranformBase or derived. */

namespace ants
{
using namespace std;

template <unsigned int ImageDimension>
int ConvertTransformFile(std::string inputFilename, std::string outFilename, bool outputMatrixOnly)
{
  typedef itk::Transform<double, ImageDimension, ImageDimension> TransformType;
  typename TransformType::Pointer transform;
  transform = itk::ants::ReadTransform<ImageDimension>( inputFilename );
  if( transform.IsNull() )
    {
    antscout << "Error while reading transform file. " << std::endl;
    return EXIT_FAILURE;
    }

  if( outputMatrixOnly )
    {
    typedef itk::MatrixOffsetTransformBase<typename TransformType::ScalarType, ImageDimension,
                                           ImageDimension> OffsetType;
    typename OffsetType::Pointer matrixOffsetTransform = dynamic_cast<OffsetType *>(transform.GetPointer() );
    if( matrixOffsetTransform.IsNull() )
      {
      antscout << "The transfrom read from file is not derived from MatrixOffsetTransformBase. Cannot output matrix."
               << std::endl;
      return EXIT_FAILURE;
      }
    typename OffsetType::MatrixType matrix = matrixOffsetTransform->GetMatrix();

    std::ofstream outputStream;
    outputStream.open(outFilename.c_str(), std::ios::out);
    if( outputStream.fail() )
      {
      outputStream.close();
      antscout << "Failed opening the output file " << outFilename << std::endl;
      return EXIT_FAILURE;
      }
    outputStream << matrix << std::endl;
    outputStream.close();
    }
  else
    {
    // Write it out as a text file using the legacy txt transform format
    int result = itk::ants::WriteTransform<ImageDimension>( transform, outFilename );
    if( result == EXIT_FAILURE )
      {
      antscout << "Failed writing transform to text format." << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

bool FileExists(string strFilename)
{
  struct stat stFileInfo;
  bool        blnReturn;
  int         intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(), &stFileInfo);
  if( intStat == 0 )
    {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
    }
  else
    {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
    }

  return blnReturn;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ConvertTransformFile( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ConvertTransformFile" );

  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  antscout->set_stream( out_stream );

  if( argc < 4  || ( strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ) )
    {
    antscout << " Utility to read in a transform file (presumed to be in binary format) " << std::endl
             << " and output it in legacy text format for human reading. " << std::endl
             << " The option '--matrix' will instead output only the transform matrix " << std::endl
             << " to a text file, one row per line with space-delimited values. " << std::endl
             << " This option works only for transforms of type MatrixOffsetTranformBase or derived." << std::endl
             << std::endl;
    antscout << "Usage:  " << argv[0]
             << " dimensions inputTransfromFile.ext outputTransformFile['.txt'|'.tfm'] [--matrix | -m] " << std::endl;
    if( argc < 4 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }

  int  dimensionPos = 1;
  int  inputFilenamePos = 2;
  int  outFilenamePos = 3;
  bool outputMatrixOnly = false;

  // User option
  if( argc > 4 )
    {
    if( strcmp(argv[4], "--matrix") == 0 || strcmp(argv[4], "-m") == 0 )
      {
      // User has requested outputting matrix information only.
      outputMatrixOnly = true;
      }
    else
      {
      antscout << "Unrecognized option: " << argv[4] << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Get the image dimension
  unsigned int dimension = atoi( argv[dimensionPos] );

  // Check the filename
  std::string inputFilename = std::string( argv[inputFilenamePos] );
  if( !FileExists(inputFilename) )
    {
    antscout << " file " << inputFilename << " does not exist . " << std::endl;
    return EXIT_FAILURE;
    }

  // Check the output filename
  std::string outFilename = std::string( argv[outFilenamePos] );
  if( itksys::SystemTools::GetFilenameLastExtension(outFilename) != ".txt" &&
      itksys::SystemTools::GetFilenameLastExtension(outFilename) != ".tfm" )
    {
    antscout << "Output filename '" << outFilename << "' must end in '.txt' or '.tfm' " << std::endl;
    return EXIT_FAILURE;
    }

  switch( dimension )
    {
    case 1:
      {
      ConvertTransformFile<1>(inputFilename, outFilename, outputMatrixOnly);
      }
      break;
    case 2:
      {
      ConvertTransformFile<2>(inputFilename, outFilename, outputMatrixOnly);
      }
      break;
    case 3:
      {
      ConvertTransformFile<3>(inputFilename, outFilename, outputMatrixOnly);
      }
      break;
    case 4:
      {
      ConvertTransformFile<4>(inputFilename, outFilename, outputMatrixOnly);
      }
      break;
    default:
      antscout << "Unsupported dimension " <<  dimension << "." << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // namespace ants
