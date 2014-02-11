#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkDeformationFieldGradientTensorImageFilter.h"
#include "itkDeterminantTensorImageFilter.h"
#include "itkGeometricJacobianDeterminantImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkMaximumImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension>
int CreateJacobianDeterminantImage( int argc, char *argv[] )
{
  typedef double RealType;
  typedef itk::Image<RealType, ImageDimension> ImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::Pointer jacobian = NULL;

  typename ImageType::Pointer minimumConstantImage = ImageType::New();
  minimumConstantImage->CopyInformation( reader->GetOutput() );
  minimumConstantImage->SetRegions( reader->GetOutput()->GetRequestedRegion() );
  minimumConstantImage->Allocate();
  minimumConstantImage->FillBuffer( 0.001 );

  bool calculateLogJacobian = false;
  if ( argc > 4 )
    {
    calculateLogJacobian = static_cast<bool>( atoi( argv[4] ) );
    }

  bool calculateGeometricJacobian = false;
  if ( argc > 5 )
    {
    calculateGeometricJacobian = static_cast<bool>( atoi( argv[5] ) );
    }

  if( calculateGeometricJacobian )
    {
    typedef itk::GeometricJacobianDeterminantImageFilter
      <VectorImageType, RealType, ImageType> JacobianFilterType;
    typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
    jacobianFilter->SetInput( reader->GetOutput() );

    jacobian = jacobianFilter->GetOutput();
    jacobian->Update();
    jacobian->DisconnectPipeline();
    }
  else
    {
    typedef itk::DeformationFieldGradientTensorImageFilter<VectorImageType, RealType> JacobianFilterType;
    typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
    jacobianFilter->SetInput( reader->GetOutput() );
    jacobianFilter->SetCalculateJacobian( true );
    jacobianFilter->SetUseImageSpacing( true );
    jacobianFilter->SetOrder( 2 );
    jacobianFilter->SetUseCenteredDifference( true );

    typedef itk::DeterminantTensorImageFilter<typename JacobianFilterType::OutputImageType, RealType>
      DeterminantFilterType;
    typename DeterminantFilterType::Pointer determinantFilter = DeterminantFilterType::New();
    determinantFilter->SetInput( jacobianFilter->GetOutput() );
    determinantFilter->Update();

    minimumConstantImage->FillBuffer( 0.0 );

    typedef itk::MaximumImageFilter<ImageType, ImageType, ImageType> MaxFilterType;
    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1( determinantFilter->GetOutput() );
    maxFilter->SetInput2( minimumConstantImage );

    jacobian = maxFilter->GetOutput();
    jacobian->Update();
    jacobian->DisconnectPipeline();
    }

  if( calculateLogJacobian )
    {
    minimumConstantImage->FillBuffer( 0.001 );

    typedef itk::MaximumImageFilter<ImageType, ImageType, ImageType> MaxFilterType;
    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1( jacobian );
    maxFilter->SetInput2( minimumConstantImage );

    typedef itk::LogImageFilter<ImageType, ImageType> LogFilterType;
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput( maxFilter->GetOutput() );
    logFilter->Update();

    typedef itk::ImageFileWriter<ImageType> ImageWriterType;
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( logFilter->GetOutput() );
    writer->Update();
    }
  else
    {
    typedef itk::ImageFileWriter<ImageType> ImageWriterType;
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( jacobian );
    writer->Update();
    }

  return EXIT_SUCCESS;
}





// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int CreateJacobianDeterminantImage( std::vector<std::string> args, std::ostream* itkNotUsed( out_stream ) )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CreateJacobianDeterminantImage" );

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

  // antscout->set_stream( out_stream );

  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension deformationField outputImage [doLogJacobian=0] [useGeometric=0]" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      {
      CreateJacobianDeterminantImage<2>( argc, argv );
      }
      break;
    case 3:
      {
      CreateJacobianDeterminantImage<3>( argc, argv );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
