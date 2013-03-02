/*=========================================================================1

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ReorientTensorImage.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/17 18:55:26 $
  Version:   $Revision: 1.2 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "antsUtilities.h"
#include <algorithm>

#include "ReadWriteImage.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int RebaseTensorImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "RebaseTensorImage" );

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

  if( argc < 5 )
    {
    antscout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <PHYSICAL/LOCAL/reference.nii.gz> "
             << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  int dim = atoi(argv[1]);

  char * moving_image_filename = argv[2];
  char * output_image_filename = argv[3];

  if( dim != 3 )
    {
    antscout << "RebaseTensorImage only supports 3D image volumes" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::DiffusionTensor3D<double> PixelType;
  typedef itk::Image<PixelType, 3>       TensorImageType;
  typedef itk::Image<float, 3>           ImageType;

  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  TensorImageType::Pointer img_mov =TensorImageType::New();

  // No reason to use log-euclidean space
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  TensorImageType::DirectionType::InternalMatrixType direction = img_mov->GetDirection().GetVnlMatrix();
  direction.set_identity();

  antscout << "Transforming space of " << moving_image_filename;

  // antscout << i << " = " << argv[i-1] << std::endl;
  char * convert = argv[4];

  if( strcmp( convert, "PHYSICAL" ) == 0 )
    {
    antscout << " -> physical space";
    direction = img_mov->GetDirection().GetVnlMatrix();
    }
  else if( strcmp( convert, "LOCAL" ) == 0 )
    {
    antscout << " -> local space";
    direction = img_mov->GetDirection().GetTranspose();
    }
  else
    {
    antscout << " -> " << convert << " space";
    ImageType::Pointer target;
    ReadImage<ImageType>( target, convert );
    direction =  img_mov->GetDirection().GetTranspose() * target->GetDirection().GetVnlMatrix();
    }

  // direction = direction.transpose(); // to accomodate for how
  // eigenvectors are stored

  antscout << std::endl;
  antscout << "Final rebasing matrix: " << std::endl << direction << std::endl;

  if( !direction.is_identity(0.00001) )
    {
    itk::ImageRegionIteratorWithIndex<TensorImageType> it( img_mov, img_mov->GetLargestPossibleRegion() );
    while( !it.IsAtEnd() )
      {
      /*
      PixelType dt = it.Value();
      PixelType::EigenValuesArrayType evalues;
      PixelType::EigenVectorsMatrixType evectors;
      dt.ComputeEigenAnalysis( evalues, evectors );

      evectors = evectors * direction;

      PixelType::EigenVectorsMatrixType emat;
      emat.Fill( 0.0 );
      for (unsigned int i=0; i<3; i++)
        {
        emat(i,i) = evalues[i];
        }

      PixelType::EigenVectorsMatrixType::InternalMatrixType matrixDT = evectors.GetTranspose() * emat.GetVnlMatrix() * evectors.GetVnlMatrix();
      */

      PixelType::EigenVectorsMatrixType::InternalMatrixType dt;
      dt(0, 0) = it.Value()[0];
      dt(0, 1) = dt(1, 0) = it.Value()[1];
      dt(0, 2) = dt(2, 0) = it.Value()[2];
      dt(1, 1) = it.Value()[3];
      dt(1, 2) = dt(2, 1) = it.Value()[4];
      dt(2, 2) = it.Value()[5];

      if( ( it.Value()[0] + it.Value()[3] + it.Value()[5] ) > 0.00001 )
        {
        PixelType::EigenVectorsMatrixType::InternalMatrixType matrixDT = direction * dt * direction.transpose();

        PixelType outDT;
        outDT[0] = matrixDT(0, 0);
        outDT[1] = matrixDT(1, 0);
        outDT[2] = matrixDT(2, 0);
        outDT[3] = matrixDT(1, 1);
        outDT[4] = matrixDT(2, 1);
        outDT[5] = matrixDT(2, 2);
        it.Set( outDT );
        }

      ++it;
      }
    }
  else
    {
    antscout << "Identity transform detected.. image unmodified" << std::endl;
    }

  // No reason to use log-euclidean space here
  WriteTensorImage<TensorImageType>(img_mov, output_image_filename, false);

  return EXIT_SUCCESS;
}
} // namespace ants
