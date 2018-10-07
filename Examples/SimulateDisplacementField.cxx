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
#include "ReadWriteData.h"

#include "itkSimulatedBSplineDisplacementFieldSource.h"
#include "itkSimulatedExponentialDisplacementFieldSource.h"

namespace ants
{

template <unsigned int ImageDimension>
int SimulateDisplacementField( int argc, char *argv[] )
{
  typedef double RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typename RealImageType::Pointer domainImage;
  ReadImage<RealImageType>( domainImage, argv[2] );

  typedef itk::SimulatedBSplineDisplacementFieldSource<DisplacementFieldType> BSplineSimulatorType;
  typename BSplineSimulatorType::ArrayType numberOfControlPoints;
  numberOfControlPoints.Fill( 8 );

  typename BSplineSimulatorType::Pointer bsplineSimulator = BSplineSimulatorType::New();
  bsplineSimulator->SetDisplacementFieldDomainFromImage( domainImage );
  bsplineSimulator->SetNumberOfRandomPoints( 10 );
  bsplineSimulator->SetNumberOfFittingLevels( 4 );
  bsplineSimulator->SetNumberOfControlPoints( numberOfControlPoints );
  bsplineSimulator->SetDisplacementNoiseStandardDeviation( 10.0 );
  bsplineSimulator->Update();

  WriteImage<DisplacementFieldType>( bsplineSimulator->GetOutput(), argv[3] );

  typedef itk::SimulatedExponentialDisplacementFieldSource<DisplacementFieldType> ExponentialSimulatorType;
  typename ExponentialSimulatorType::Pointer exponentialSimulator = ExponentialSimulatorType::New();
  exponentialSimulator->SetDisplacementFieldDomainFromImage( domainImage );
  exponentialSimulator->SetNumberOfRandomPoints( 1000 );
  exponentialSimulator->SetDisplacementNoiseStandardDeviation( 10.0 );
  exponentialSimulator->SetSmoothingStandardDeviation( 10.0 );
  exponentialSimulator->Update();

  WriteImage<DisplacementFieldType>( exponentialSimulator->GetOutput(), argv[3] );



  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SimulateDisplacementField( std::vector<std::string> args, std::ostream* /*out_stream = ITK_NULLPTR */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SimulateDisplacementField" );

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
  argv[argc] = ITK_NULLPTR;
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

  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension domainImage outputField" << std::endl;

    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     return SimulateDisplacementField<2>( argc, argv );
     break;
   case 3:
     return SimulateDisplacementField<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
   }
  return EXIT_SUCCESS;
}

} // namespace ants
