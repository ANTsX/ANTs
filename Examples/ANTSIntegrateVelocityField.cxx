

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkANTSImageRegistrationOptimizer.h"
#include "itkTimeVaryingVelocityFieldIntegrationImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkTimeVaryingVelocityFieldTransform.h"

#include "itkImageFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int IntegrateVelocityField(int argc, char *argv[])
{
  int         argct = 1;
  std::string imgfn = std::string(argv[argct]); argct++;

  std::string vectorfn = std::string(argv[argct]); argct++;
  std::string outname = std::string(argv[argct]); argct++;

  typedef float PixelType;
  PixelType timezero = 0;
  PixelType timeone = 1;
  PixelType dT = 0.01;
  if( argc > argct )
    {
    timezero = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    timeone = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    dT = atof(argv[argct]);
    }
  argct++;
  std::cout << " time-0 " << timezero << " dt " << dT << " time-1 " << timeone << std::endl;
  PixelType starttime = timezero;
  PixelType finishtime = timeone;
  typedef float                                       PixelType;
  typedef itk::Vector<PixelType, ImageDimension>      VectorType;
  typedef itk::Image<VectorType, ImageDimension>      DisplacementFieldType;
  typedef itk::Image<VectorType, ImageDimension + 1>  TimeVaryingVelocityFieldType;
  typedef itk::Image<PixelType, ImageDimension>       ImageType;

  typename ImageType::Pointer image;
  ReadImage<ImageType>(image, imgfn.c_str() );
  typedef TimeVaryingVelocityFieldType                tvt;
  typename tvt::Pointer timeVaryingVelocity;
  ReadImage<tvt>(timeVaryingVelocity, vectorfn.c_str() );

  VectorType zero;
  zero.Fill(0);
  typename DisplacementFieldType::Pointer deformation =
    AllocImage<DisplacementFieldType>(image, zero);

  if( !timeVaryingVelocity )
    {
    std::cerr << " No TV Field " << std::endl;  return EXIT_FAILURE;
    }

  if( starttime < 0 )
    {
    starttime = 0;
    }
  if( starttime > 1 )
    {
    starttime = 1;
    }
  if( finishtime < 0 )
    {
    finishtime = 0;
    }
  if( finishtime > 1 )
    {
    finishtime = 1;
    }

  typedef itk::TimeVaryingVelocityFieldIntegrationImageFilter
    <TimeVaryingVelocityFieldType, DisplacementFieldType> IntegratorType;
  typename IntegratorType::Pointer integrator = IntegratorType::New();
  integrator->SetInput( timeVaryingVelocity );
  integrator->SetLowerTimeBound( starttime );
  integrator->SetUpperTimeBound( finishtime );
  integrator->SetNumberOfIntegrationSteps( (unsigned int ) 1 / dT );
  integrator->Update();

  WriteImage<DisplacementFieldType>( integrator->GetOutput(), outname.c_str() );
  return 0;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSIntegrateVelocityField( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTSIntegrateVelocityField" );
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

  if( argc < 4 )
    {
    std::cerr << "Usage:   " << argv[0]
             << " reference_image  VelocityIn.mhd DeformationOut.nii.gz  time0 time1 dT  " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }
  std::cout << " start " << std::endl;
  std::string               ifn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(ifn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(ifn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim =  imageIO->GetNumberOfDimensions();
  std::cout << " dim " << dim << std::endl;

  switch( dim )
    {
    case 2:
      {
      IntegrateVelocityField<2>(argc, argv);
      }
      break;
    case 3:
      {
      IntegrateVelocityField<3>(argc, argv);
      }
      break;
    case 4:
      {
      IntegrateVelocityField<4>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // namespace ants
