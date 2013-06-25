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

#include "itkTransformFactory.h"
#include "itkAffineTransform.h"
#include "itkTranslationTransform.h"
#include "itkIdentityTransform.h"

#include "itkantsReadWriteTransform.h"

/* Utility to read in a transform file (presumed to be in binary format) and output
 * it in one of several different formats, defaulting to legacy text format for human reading.
 * Options are available to instead output only a transform matrix to a text file,
 * one row per dimension with space-delimited values. This option works only for
 * transforms of MatrixOffsetTranformBase or derived, Translation and Identity transforms. */

namespace ants
{
using namespace std;

/*
 *
 */
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

/*
 *
 */
template <class TTransform>
bool GetMatrix( const typename TTransform::Pointer & transform, typename TTransform::MatrixType & matrix,
                bool outputRAS )
{
  const unsigned int ImageDimension = TTransform::InputSpaceDimension;

  typedef typename TTransform::ScalarType ScalarType;

  matrix.Fill( itk::NumericTraits<typename TTransform::ScalarType>::Zero );
  bool done = false;

  // Matrix-offset derived
    {
    typedef itk::MatrixOffsetTransformBase<ScalarType, ImageDimension, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      matrix = castTransform->GetMatrix();
      done = true;
      }
    }

  // Translation
  if( !done )
    {
    typedef itk::TranslationTransform<ScalarType, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        matrix(i, i) = itk::NumericTraits<typename TTransform::ScalarType>::One;
        }
      done = true;
      }
    }

  // Identity
  if( !done )
    {
    typedef itk::IdentityTransform<ScalarType, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        matrix(i, i) = itk::NumericTraits<typename TTransform::ScalarType>::One;
        }
      done = true;
      }
    }

  if( !done )
    {
    // Unsupported transform type
    return false;
    }

  if( outputRAS )
    {
    // Convert to RAS coordinate system. ITK uses LPS.
    // x and y dimensions are flipped.
    // This code is from c3d app.
    vnl_vector<ScalarType> v_lps_to_ras(ImageDimension, 1.0);
    v_lps_to_ras[0] = -1.0;
    if( ImageDimension > 1 )
      {
      v_lps_to_ras[1] = -1.0;
      }
    vnl_diag_matrix<ScalarType> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<ScalarType>      mold = matrix.GetVnlMatrix();
    matrix.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    }

  return true;
}

/*
 *
 */
template <class TTransform, class TMatrix>
bool GetHomogeneousMatrix( const typename TTransform::Pointer & transform, TMatrix & hMatrix, bool outputRAS )
{
  const unsigned int ImageDimension = TTransform::InputSpaceDimension;

  typedef typename TTransform::ScalarType ScalarType;

  hMatrix.Fill( itk::NumericTraits<ScalarType>::Zero );

  bool done = false;

  // Get the NxN matrix
  typename TTransform::MatrixType matrix;
  if( !GetMatrix<TTransform>( transform, matrix, outputRAS ) )
    {
    return false;
    }
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      hMatrix(i, j) = matrix(i, j);
      }
    }

  // Set the lower-right corner to 1
  unsigned int corner = ImageDimension;
  hMatrix(corner, corner) = itk::NumericTraits<typename TTransform::ScalarType>::One;

  //
  // Get the offset
  //

  // Identity
    {
    typedef itk::IdentityTransform<ScalarType, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      // Nothing more to do here.
      return true;
      }
    }

  typename TTransform::OutputVectorType offset;
  offset.Fill( itk::NumericTraits<typename TTransform::ScalarType>::Zero );

  // Matrix-offset derived
    {
    typedef itk::MatrixOffsetTransformBase<ScalarType, ImageDimension, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      offset = castTransform->GetOffset();
      done = true;
      }
    }

  // Translation
  if( !done )
    {
    typedef itk::TranslationTransform<ScalarType, ImageDimension> CastTransformType;
    typename CastTransformType::Pointer castTransform = dynamic_cast<CastTransformType *>(transform.GetPointer() );

    if( castTransform.IsNotNull() )
      {
      offset = castTransform->GetOffset();
      done = true;
      }
    }

  if( !done )
    {
    // Unsupported transform type
    return false;
    }

  if( outputRAS )
    {
    // Convert to RAS coordinate system. ITK uses LPS.
    // x and y dimensions are flipped.
    // This code is from c3d app.
    offset[0] *= -1.0;
    if( ImageDimension > 1 )
      {
      offset[1] *= -1.0;
      }
    }
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    hMatrix(i, ImageDimension) = offset[i];
    }

  return true;
}

/*
 *
 */
template <unsigned int ImageDimension>
int ConvertTransformFile(int argc, char* argv[])
{
  int  inputFilenamePos = 2;
  int  outFilenamePos = 3;
  bool outputMatrix = false;
  bool outputHomogeneousMatrix = false;
  bool outputAffine = false;
  bool outputRAS = false;

  // User options
  if( argc > 4 )
    {
    for( int n = 4; n < argc; n++ )
      {
      if( strcmp(argv[n], "--matrix") == 0 || strcmp(argv[n], "-m") == 0 )
        {
        // User has requested outputting matrix information only.
        outputMatrix = true;
        }
      else if( strcmp(argv[n], "--homogeneousMatrix") == 0 || strcmp(argv[n], "--hm") == 0 )
        {
        // User has requested outputting homogeneous matrix information only.
        outputHomogeneousMatrix = true;
        }
      else if( strcmp(argv[n], "--convertToAffineType") == 0 )
        {
        outputAffine = true;
        }
      else if( strcmp(argv[n], "--RAS") == 0 || strcmp(argv[n], "--ras") == 0 )
        {
        outputRAS = true;
        }
      else
        {
        antscout << "Unrecognized option: " << argv[n] << std::endl;
        return EXIT_FAILURE;
        }
      }
    if( outputRAS && !outputMatrix && !outputHomogeneousMatrix )
      {
      antscout << " '--RAS' option must be used with either of 'matrix' or 'homongeneousMatrix' options." << std::endl;
      return EXIT_FAILURE;
      }
    if( (outputMatrix &&
         outputHomogeneousMatrix) || (outputMatrix && outputAffine) || (outputHomogeneousMatrix && outputAffine) )
      {
      antscout << "Only one primary output option allowed at once." << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Check the filename
  std::string inputFilename = std::string( argv[inputFilenamePos] );
  if( !FileExists(inputFilename) )
    {
    antscout << " file " << inputFilename << " does not exist . " << std::endl;
    return EXIT_FAILURE;
    }

  // Get the output filename
  std::string outFilename = std::string( argv[outFilenamePos] );

  // Read the transform
  typedef itk::Transform<double, ImageDimension, ImageDimension> TransformType;
  typename TransformType::Pointer transform;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> baseTransformType;
  itk::TransformFactory<baseTransformType>::RegisterTransform();
  transform = itk::ants::ReadTransform<double, ImageDimension>( inputFilename );
  if( transform.IsNull() )
    {
    antscout << "Error while reading transform file. Did you specify the correct dimension?" << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Outputs
  //
  if( outputMatrix || outputHomogeneousMatrix )
    {
    std::ofstream outputStream;
    outputStream.open(outFilename.c_str(), std::ios::out);
    if( outputStream.fail() )
      {
      outputStream.close();
      antscout << "Failed opening the output file " << outFilename << std::endl;
      return EXIT_FAILURE;
      }

    if( outputMatrix )
      {
      typedef itk::Matrix<typename TransformType::ScalarType, ImageDimension, ImageDimension> MatrixType;
      MatrixType matrix;
      if( GetMatrix<TransformType>( transform, matrix, outputRAS ) )
        {
        outputStream << matrix;
        }
      else
        {
        antscout << "Error. Transform type is unsupported for getting matrix: " << transform->GetNameOfClass()
                 << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      // Homogeneous matrix
      typedef itk::Matrix<typename TransformType::ScalarType, ImageDimension + 1, ImageDimension + 1> MatrixType;
      MatrixType hMatrix;

      if( GetHomogeneousMatrix<TransformType, MatrixType>( transform, hMatrix, outputRAS ) )
        {
        outputStream << hMatrix;
        }
      else
        {
        antscout << "Error. Transform type is unsupported for getting matrix: " << transform->GetNameOfClass()
                 << std::endl;
        return EXIT_FAILURE;
        }
      }
    outputStream.close();
    return EXIT_SUCCESS;
    }

  if( outputAffine )
    {
    // Convert to Affine and output as binary.
    // This is done by taking the matrix and offset from the transform
    // and assigning them to a new affine transform.
    typedef itk::MatrixOffsetTransformBase<typename TransformType::ScalarType, ImageDimension,
                                           ImageDimension> CastTransformType;
    typename CastTransformType::Pointer matrixOffsetTransform =
      dynamic_cast<CastTransformType *>(transform.GetPointer() );
    if( matrixOffsetTransform.IsNull() )
      {
      antscout
        << "The transfrom read from file is not derived from MatrixOffsetTransformBase. Cannot convert to Affine."
        << std::endl;
      return EXIT_FAILURE;
      }
    if( itksys::SystemTools::GetFilenameLastExtension(outFilename) != ".mat" )
      {
      antscout << "Output filename '" << outFilename << "' must end in '.mat' for binary output." << std::endl;
      return EXIT_FAILURE;
      }
    typedef itk::AffineTransform<typename TransformType::ScalarType, ImageDimension> AffineTransformType;
    typename AffineTransformType::Pointer newAffineTransform = AffineTransformType::New();
    newAffineTransform->SetMatrix( matrixOffsetTransform->GetMatrix() );
    newAffineTransform->SetOffset( matrixOffsetTransform->GetOffset() );
    transform = dynamic_cast<TransformType *>(newAffineTransform.GetPointer() );
    if( transform.IsNull() )
      {
      antscout << "Unexpected error casting from affine transform to transform type." << std::endl;
      return EXIT_FAILURE;
      }
    int result = itk::ants::WriteTransform<double, ImageDimension>( transform, outFilename );
    if( result == EXIT_FAILURE )
      {
      antscout << "Failed writing converted transform to binary format." << std::endl;
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }

  // Default behavior.
  // Write it out as a text file using the legacy txt transform format
  if( itksys::SystemTools::GetFilenameLastExtension(outFilename) != ".txt" &&
      itksys::SystemTools::GetFilenameLastExtension(outFilename) != ".tfm" )
    {
    antscout << "Output filename '" << outFilename << "' must end in '.txt' or '.tfm' for text-format output."
             << std::endl;
    return EXIT_FAILURE;
    }
  int result = itk::ants::WriteTransform<double, ImageDimension>( transform, outFilename );
  if( result == EXIT_FAILURE )
    {
    antscout << "Failed writing transform to text format." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
/*
 *
 */
int ConvertTransformFile( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ConvertTransformFile" );

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
    antscout << "USAGE:  " << std::endl
             << " " << argv[0] << " dimensions inputTransfromFile.ext outputTransformFile.ext [OPTIONS]" << std::endl
             << std::endl;
    antscout << "COMMAND: " << std::endl
             << " Utility to read in a transform file (presumed to be in binary format) " << std::endl
             << " and output it in various formats. Default output is legacy human-readable" << std::endl
             << " text format.  Without any options, the output filename extension must be " << std::endl
             << " .txt or .tfm to signify a text-formatted transform file. " << std::endl
             << std::endl
             << " OPTIONS: " << std::endl
             << std::endl
             << " --matrix, -m " << std::endl
             << "   Output only the transform matrix (from transform::GetMatrix() )" << std::endl
             << "   to a text file, one row per line with space-delimited values. " << std::endl
             << "   Only works for transforms of type identity, translation or " << std::endl
             << "   MatrixOffsetTranformBase and its derived types." << std::endl
             << "   The output filename must end in '.mat'." << std::endl
             << std::endl
             << " --homogeneousMatrix, --hm" << std::endl
             << "   Output an N+1 square homogeneous matrix from the transform matrix and offset." << std::endl
             << "   Only works for transforms of type identity, translation or " << std::endl
             << "   MatrixOffsetTranformBase and its derived types." << std::endl
             << "   The output filename must end in '.mat'." << std::endl
             << std::endl
             << " --RAS, --ras" << std::endl
             << "   Combined with the 'matrix' or 'homogeneousMatrix' options, this will convert" << std::endl
             << "   the output into the RAS coordinate system (Right, Anterior, Superior)." << std::endl
             << "   Otherwise, the output is in the LPS coordinate system (Left, Posterior," << std::endl
             << "   Superior), which is used by ITK. RAS is used, for example, by Slicer. " << std::endl
             << std::endl
             << " --convertToAffineType" << std::endl
             << "   Convert the input transform type to AffineTransform using the transform's " << std::endl
             << "   matrix and offset, and output again as as a binary transform file." << std::endl
             << "   This is useful for using transforms in programs" << std::endl
             << "   that do not register all available Transform factory types." << std::endl
             << std::endl;
    if( argc < 4 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  if( argc > 6 )
    {
    antscout << "Too many arguments." << std::endl;
    return EXIT_FAILURE;
    }

  // Get the image dimension
  unsigned int dimension = atoi( argv[1] );

  switch( dimension )
    {
    case 1:
      {
      return ConvertTransformFile<1>(argc, argv);
      }
      break;
    case 2:
      {
      return ConvertTransformFile<2>(argc, argv);
      }
      break;
    case 3:
      {
      return ConvertTransformFile<3>(argc, argv);
      }
      break;
    case 4:
      {
      return ConvertTransformFile<4>(argc, argv);
      }
      break;
    default:
      antscout << "Unsupported dimension " <<  dimension << "." << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // namespace ants
