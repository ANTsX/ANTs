/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ANTS.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:05 $
  Version:   $Revision: 1.19 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include "itkPICSLAdvancedNormalizationToolKit.h"
#include "itkANTSImageTransformation.h"
#include "itkANTSImageRegistrationOptimizer.h"
#include <algorithm>
#include <string>

namespace ants
{
void PrintCommandLineHelp( const std::string & progName )
{
  antscout <<  " \n " << std::endl;
  antscout <<  "Example usage: \n " << std::endl;
  antscout << progName
           <<
    " ImageDimension -m MI[fixedimage.nii.gz,movingimage.nii.gz,1,32] -o Outputfname.nii.gz -i 30x20x0 -r Gauss[3,1] -t Elast[3] \n \n "
           << std::endl;
  antscout << " Compulsory arguments:\n " << std::endl;
  antscout << " ImageDimension: 2 or 3 (for 2 or 3 Dimensional registration)\n " << std::endl;
  antscout << " -m:    Type of similarity model used for registration. \n " << std::endl;
  antscout << "    For intramodal image registration, use: " << std::endl;
  antscout << "        CC = cross-correlation " << std::endl;
  antscout << "        MI = mutual information " << std::endl;
  antscout << "        PR = probability mapping " << std::endl;
  antscout << "        MSQ = mean square difference " << std::endl;
  antscout << " \n " << std::endl;
  antscout << "    For intermodal image registration, use: " << std::endl;
  antscout << "        MI = mutual information " << std::endl;
  antscout << "        PR = probability mapping " << std::endl;
  antscout << " \n " << std::endl;
  antscout << " -o     Outputfname.nii.gz: the name of the resulting image.\n " << std::endl;
  antscout << " -i     Max-iterations in format: JxKxL, where: " << std::endl;
  antscout << "        J = max iterations at coarsest resolution (here, reduce by power of 2^2) " << std::endl;
  antscout << "        K = middle resolution iterations (here,reduce by power of 2) " << std::endl;
  antscout
    <<
    "        L = fine resolution iterations (here, full resolution). This level takes much more time per iteration!\n "
    << std::endl;
  antscout
    << "        Adding an extra value before JxKxL (i.e. resulting in IxJxKxL) would add another iteration level.\n "
    << std::endl;
  antscout << " -r     Regularization \n" << std::endl;
  antscout << " -t     Type of transformation model used for registration \n" << std::endl;
  antscout << "    For elastic image registration, use: " << std::endl;
  antscout << "        Elast = elastic transformation model (less deformation possible)\n " << std::endl;
  antscout << "    For diffeomorphic image registration, use: " << std::endl;
  antscout
    <<
    "        Syn[GradStep,TimePoints,IntegrationStep] --geodesic 2 = SyN with time with arbitrary number of time points in time discretization  "
    << std::endl;
  antscout
    <<
    "        SyN[GradStep,2,IntegrationStep] = SyN with time optimized specifically for 2 time points in the time discretization "
    << std::endl;
  antscout << "        SyN[GradStep] = Greedy SyN, typicall GradStep=0.25  " << std::endl;
  antscout << "        Exp[GradStep,TimePoints] = Exponential " << std::endl;
  antscout << "        GreedyExp = Diffeomorphic Demons style exponential mapping " << std::endl;
  antscout << " \n " << std::endl;
  antscout
    <<
    " Please use the `ANTS -h ` call or refer to the ANTS.pdf manual or antsIntroduction.sh script for additional information and typical values for transformation models\n "
    << std::endl;
}

template <unsigned int ImageDimension>
int ANTSex(int argc, char *argv[])
{
  typedef itk::PICSLAdvancedNormalizationToolKit<ImageDimension, float> RegistrationType;
  typename RegistrationType::Pointer registration = RegistrationType::New();
  registration->ParseCommandLine( argc, argv );
  antscout << " Run Reg " << std::endl;
  try
    {
    registration->RunRegistration();
    }
  catch( std::exception const& e )
    {
    antscout << "Exception caught in ANTS: " << std::endl << e.what() << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    antscout << "Non-standard exception caught in ANTS. No more information available." << std::endl;
    return EXIT_FAILURE;
    }

  registration->GetTransformationModel()->SetWriteComponentImages(true);
  registration->GetTransformationModel()->Write();
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTS( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTS" );
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

//   if( argc >= 2 && ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") )
// )
//     {
//     ants::PrintCommandLineHelp(argv[0]);
//     if( argc < 2 )
//       {
//       return EXIT_FAILURE;
//       }
//     return EXIT_SUCCESS;
//     }

  int dim = 0;
  if( argc > 1 )
    {
    dim = atoi( argv[1] );
    }

//   if( dim <= 1 || dim > 3 )
//     {
//     antscout << " You passed ImageDimension: " << dim
//              << " . Please use only 2 or 3 (for 2 or 3 Dimensional registration)  " << std::endl;
//     ants::PrintCommandLineHelp(argv[0]);
//     return EXIT_FAILURE;
//     }

  /**
   * Try the simple case of the call "ANTS fixedImage movingImage"
   */
  if( argc == 3 && ( atoi( argv[1] ) != '-' || atoi( argv[1] ) != 2 || atoi( argv[1] ) != 3 ) )
    {
    itk::ImageIOBase::Pointer fixedImageIO
      = itk::ImageIOFactory::CreateImageIO( argv[1], itk::ImageIOFactory::ReadMode );
    if( fixedImageIO.IsNull() )
      {
      antscout << "Invalid fixed image: " << argv[1] << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer movingImageIO
      = itk::ImageIOFactory::CreateImageIO( argv[2], itk::ImageIOFactory::ReadMode );
    if( movingImageIO.IsNull() )
      {
      antscout << "Invalid moving image: " << argv[2] << std::endl;
      return EXIT_FAILURE;
      }
    fixedImageIO->SetFileName( argv[1] );
    fixedImageIO->ReadImageInformation();
    movingImageIO->SetFileName( argv[2] );
    movingImageIO->ReadImageInformation();

    const unsigned int fdim = fixedImageIO->GetNumberOfDimensions();
    const unsigned int mdim = movingImageIO->GetNumberOfDimensions();
    if( fdim != mdim )
      {
      antscout << "Fixed image dimension does not equal "
               << "the moving image dimension (" << fdim << " != " << mdim << ")"
               << std::endl;
      return EXIT_FAILURE;
      }
    if( fdim != 2 && fdim != 3 )
      {
      antscout << "Unsupported image dimension" << std::endl;
      return EXIT_FAILURE;
      }

    /**
     * After the checking, we can add the default parameters;
     */

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
        vcl_log( (double)maxSize / 32 ) / vcl_log( (double) 2 ) ) + 1;
    std::string iterations( " -i " );
    for( int n = numberOfLevels; n > 0; n-- )
      {
      std::ostringstream buf;
      unsigned long      numberOfIterations = ( 5 << (n - 1) );
      buf << numberOfIterations << "x";
      iterations += buf.str();
      }
    iterations = iterations.substr( 0, iterations.length() - 1 );

    std::ostringstream dimBuf;
    dimBuf << fdim;
    std::string arguments;
    arguments.clear();
    arguments = std::string( " ANTS " ) + dimBuf.str() + iterations
      + transformation + regularization + metric + outputNaming;

    dim = fdim;

    antscout << arguments << std::endl;

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
    char * *     my_argv = new char *[my_argc];
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
        return ANTSex<3>( my_argc, my_argv );
        }
        break;
      default:
        return ANTSex<2>( my_argc, my_argv );
      }
    }
  else
    {
    switch( dim )
      {
      case 3:
        {
        return ANTSex<3>( argc, argv );
        }
        break;
      default:
        return ANTSex<2>( argc, argv );
      }
    }

  antscout << "Shoudln't have gotten here." << std::endl;
  return EXIT_FAILURE;
}
} // namespace ants
