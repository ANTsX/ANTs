/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ANTS.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:05 $
  Version:   $Revision: 1.19 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkPICSLAdvancedNormalizationToolKit.h"
#include "itkANTSImageTransformation.h"
#include "itkANTSImageRegistrationOptimizer.h"

#include <string>

template <unsigned int ImageDimension>
int ANTSex(int argc, char *argv[])
{
  typedef itk::PICSLAdvancedNormalizationToolKit<ImageDimension> RegistrationType;
  typename RegistrationType::Pointer registration = RegistrationType::New();
  registration->ParseCommandLine( argc, argv );
  std::cout << " Run Reg " << std::endl;
  try
    {
    registration->RunRegistration();
    }
  catch( ... )
    {
    std::cerr << "Exception thrown: ANTS" << std::endl;
    return EXIT_FAILURE;
    }
  registration->GetTransformationModel()->SetWriteComponentImages(true);
  registration->GetTransformationModel()->Write();

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[] )
{
  int dim = 2;

//  while(argc--) printf("%s\n", *argv++);
  if( argc < 2 )
    {
    std::cout <<  "  Call  to get detailed usage :  \n    ANTS -h  \n or  \n ANTS --help \n \n  " << std::endl;
    std::cout
    <<
    " typical use, for a 3D image matching with mutual information : \n  ANTS  3  -m MI[fixedimage.nii,movingimage.nii,1,32] -o MutualInformationMapping.nii -i 30x20x0 -r Gauss[3,1] -t Elast[3] "
    << std::endl;
    std::cout << " the value 32 in the metric call ( -m MI[ ...] ) is the joint histogram #of bins " << std::endl;
    std::cout << " the transformation model is elastic with gradient descent step-length of 3 " << std::endl;
    std::cout
    <<
    "  3 levels of optimization with 30 iterations at the coarsest (fastest) level, 20 at the next and 0 at full resolution "
    << std::endl;
    std::cout
    <<
    " \n typical use, for a 3D image matching with correlation : \n  ANTS  3  -m PR[fixedimage.nii,movingimage.nii,1,2] -o CC.nii -i 30x20x0 -r Gauss[3,0] -t SyN[1] "
    << std::endl;
    std::cout << " the value 2 in the metric call ( -m PR[ ...] ) is the radius of the correlation window "
              << std::endl;
    std::cout
    << " the transformation is the greedy symmetric diffeomorphic mapping  with xgradient descent step-length of  1"
    << std::endl;
    std::cout
    <<
    " \n typical use, for a 3D image matching with correlation implementation 2 : \n  ANTS  3  -m CC[fixedimage.nii,movingimage.nii,1,2] -o CC2.nii -i 30x20x0 -r Gauss[3,0.5] -t SyN[1,2,0.1] "
    << std::endl;
    std::cout << " the value 2 in the metric call ( -m CC[ ...] ) is the radius of the correlation window "
              << std::endl;
    std::cout
    <<
    " the transformation is the symmetric diffeomorphic mapping with time-varying velocity fied - 10 integration points in time and gradient step length of 1 "
    << std::endl;
    std::cout <<  " \n Call  to get detailed usage :  \n    ANTS -h  \n or  \n ANTS --help \n \n  " << std::endl;
    return 1;
    }
  else
    {
    dim = atoi( argv[1] );
    }

  if( dim <= 1 || dim > 3 )
    {
    std::cout << " you passed image dimension (first argument) as " << dim
              << " ANTS does not function with images of this dimension " << std::endl;
    exit(1);
    }
  /**
   * Try the simple case of the call "ANTS fixedImage movingImage"
   */
  if( argc == 3 && ( atoi( argv[1] ) != 2 || atoi( argv[1] ) != 3 ) )
    {
    itk::ImageIOBase::Pointer fixedImageIO
      = itk::ImageIOFactory::CreateImageIO(
          argv[1], itk::ImageIOFactory::ReadMode );
    if( fixedImageIO.IsNull() )
      {
      std::cerr << "Invalid fixed image: " << argv[1] << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer movingImageIO
      = itk::ImageIOFactory::CreateImageIO(
          argv[2], itk::ImageIOFactory::ReadMode );
    if( movingImageIO.IsNull() )
      {
      std::cerr << "Invalid moving image: " << argv[2] << std::endl;
      return EXIT_FAILURE;
      }
    fixedImageIO->SetFileName( argv[1] );
    fixedImageIO->ReadImageInformation();
    movingImageIO->SetFileName( argv[2] );
    movingImageIO->ReadImageInformation();

    unsigned int fdim = fixedImageIO->GetNumberOfDimensions();
    unsigned int mdim = movingImageIO->GetNumberOfDimensions();

    if( fdim != mdim )
      {
      std::cerr << "Fixed image dimension does not equal "
                << "the moving image dimension (" << fdim << " != " << mdim << ")"
                << std::endl;
      return EXIT_FAILURE;
      }
    if( fdim != 2 && fdim != 3 )
      {
      std::cerr << "Unsupported image dimension" << std::endl;
      return EXIT_FAILURE;
      }

    /**
     * After the checking, we can add the default parameters;
     */
    std::string arguments;

    std::string transformation( " -t SyN[1.0] ");
    std::string regularization( " -r Gauss[3,0.5] ");
    std::string metric = std::string( " -m PR[" ) + std::string( argv[1] )
      + std::string( "," ) + std::string( argv[2] ) + std::string( ",1,3] " );
    std::string outputNaming( " -o ANTS.nii.gz " );

    long maxSize = vnl_math_min( fixedImageIO->GetDimensions( 0 ),
                                 movingImageIO->GetDimensions( 0 ) );
    for( unsigned int d = 1; d < fdim; d++ )
      {
      long tmpMax = vnl_math_max( fixedImageIO->GetDimensions( d ),
                                  movingImageIO->GetDimensions( d ) );
      if( maxSize < tmpMax )
        {
        maxSize = tmpMax;
        }
      }

    unsigned int numberOfLevels = static_cast<unsigned int>(
        vcl_log( maxSize / 32 ) / vcl_log( 2 ) ) + 1;
    std::string iterations( " -i " );
    for( int n = numberOfLevels; n > 0; n-- )
      {
      itk::OStringStream buf;
      unsigned long      numberOfIterations = ( 5 << (n - 1) );
      buf << numberOfIterations << "x";
      iterations += buf.str();
      }
    iterations = iterations.substr( 0, iterations.length() - 1 );

    itk::OStringStream dimBuf;
    dimBuf << fdim;
    arguments.clear();
    arguments = std::string( " ANTS " ) + dimBuf.str() + iterations
      + transformation + regularization + metric + outputNaming;

    dim = fdim;

    std::cout << arguments << std::endl;

    unsigned int           my_argc = 0;
    std::string::size_type delimPos = 0;
    std::string::size_type tokenPos = 0;
    std::string::size_type pos = 0;
    while( true )
      {
      delimPos = arguments.find_first_of( " ", pos );
      tokenPos = arguments.find_first_not_of( " ", pos );
      if( std::string::npos != delimPos &&
          std::string::npos != tokenPos && tokenPos < delimPos )
        {
        my_argc++;
        }
      else if( std::string::npos == delimPos )
        {
        if( std::string::npos != tokenPos )
          {
          my_argc++;
          }
        else
          {
          break;
          }
        }
      pos = delimPos + 1;
      }

    unsigned int arg_count = 0;
    char *       my_argv[my_argc];
    delimPos = 0;
    tokenPos = 0;
    pos = 0;
    while( true )
      {
      delimPos = arguments.find_first_of( " ", pos );
      tokenPos = arguments.find_first_not_of( " ", pos );
      if( std::string::npos != delimPos &&
          std::string::npos != tokenPos && tokenPos < delimPos )
        {
        std::string arg = arguments.substr( pos, delimPos - pos );
        my_argv[arg_count] = new char[arg.size() + 1];
        strcpy( my_argv[arg_count], arg.c_str() );
        arg_count++;
        }
      else if( std::string::npos == delimPos )
        {
        if( std::string::npos != tokenPos )
          {
          std::string arg = arguments.substr( pos );
          my_argv[arg_count] = new char[arg.size() + 1];
          strcpy( my_argv[arg_count], arg.c_str() );
          arg_count++;
          }
        else
          {
          break;
          }
        }
      pos = delimPos + 1;
      }

    switch( dim )
      {
      case 3:
        {
        ANTSex<3>( my_argc, my_argv );
        }
        break;
      default:
        ANTSex<2>( my_argc, my_argv );
      }
    }
  else
    {
    switch( dim )
      {
      case 3:
        {
        ANTSex<3>( argc, argv );
        }
        break;
      default:
        ANTSex<2>( argc, argv );
      }
    }

  return EXIT_SUCCESS;
}
