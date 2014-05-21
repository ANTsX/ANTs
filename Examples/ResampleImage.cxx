
#include "antsUtilities.h"
#include <algorithm>

#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"

#include "itkConstantBoundaryCondition.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include <string>
#include <vector>

namespace ants
{

template <unsigned int ImageDimension>
int ResampleImage( int argc, char *argv[] )
{
  typedef double                                RealType;
  typedef double                                PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::IdentityTransform<RealType, ImageDimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::LinearInterpolateImageFunction<ImageType, RealType>
    LinearInterpolatorType;
  typename LinearInterpolatorType::Pointer interpolator
    = LinearInterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, RealType>
    NearestNeighborInterpolatorType;
  typename NearestNeighborInterpolatorType::Pointer nn_interpolator
    = NearestNeighborInterpolatorType::New();
  nn_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::BSplineInterpolateImageFunction<ImageType, RealType>
    BSplineInterpolatorType;
  typename BSplineInterpolatorType::Pointer bs_interpolator
    = BSplineInterpolatorType::New();
  bs_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::GaussianInterpolateImageFunction<ImageType, RealType>
    GaussianInterpolatorType;
  typename GaussianInterpolatorType::Pointer g_interpolator
    = GaussianInterpolatorType::New();
  g_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3> HammingInterpolatorType;
  typename HammingInterpolatorType::Pointer sh_interpolator = HammingInterpolatorType::New();
  sh_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::CosineWindowFunction<3> > Sinc1InterpolatorType;
  typename Sinc1InterpolatorType::Pointer sc_interpolator = Sinc1InterpolatorType::New();
  sc_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::WelchWindowFunction<3> > Sinc2InterpolatorType;
  typename Sinc2InterpolatorType::Pointer sw_interpolator = Sinc2InterpolatorType::New();
  sw_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::LanczosWindowFunction<3> > Sinc3InterpolatorType;
  typename Sinc3InterpolatorType::Pointer sl_interpolator = Sinc3InterpolatorType::New();
  sl_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::BlackmanWindowFunction<3> > Sinc4InterpolatorType;
  typename Sinc3InterpolatorType::Pointer sb_interpolator = Sinc3InterpolatorType::New();
  sb_interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  typename ResamplerType::Pointer resampler = ResamplerType::New();
  typename ResamplerType::SpacingType spacing;
  typename ResamplerType::SizeType size;

  std::vector<RealType> sp = ConvertVector<RealType>( std::string( argv[4] ) );

  if( argc <= 5 || atoi( argv[5] ) == 0 )
    {
    if( sp.size() == 1 )
      {
      spacing.Fill( sp[0] );
      }
    else if( sp.size() == ImageDimension )
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        spacing[d] = sp[d];
        }
      }
    else
      {
      std::cout << "Invalid spacing." << std::endl;
      }
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      RealType spacing_old = reader->GetOutput()->GetSpacing()[i];
      RealType size_old = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
      size[i] = static_cast<int>( ( spacing_old * size_old ) / spacing[i] + 0.5 );
      }
    }
  else
    {
    if( sp.size() == 1 )
      {
      size.Fill( static_cast<unsigned int>( sp[0] ) );
      }
    else if( sp.size() == ImageDimension )
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        size[d] = static_cast<unsigned int>( sp[d] );
        }
      }
    else
      {
      std::cout << "Invalid size." << std::endl;
      }
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      RealType spacing_old = reader->GetOutput()->GetSpacing()[i];
      RealType size_old = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
      spacing[i] = spacing_old * static_cast<float>( size_old - 1.0 )
        / static_cast<float>( size[i] - 1.0 );
      }
    }

  char arg7 = '\0';
  if( argc > 7 )
    {
    arg7 = *argv[7];
    }

  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  if( argc > 6 && atoi( argv[6] ) )
    {
    switch( atoi( argv[6] ) )
      {
      case 0: default:
        {
        resampler->SetInterpolator( interpolator );
        }
        break;
      case 1:
        {
        resampler->SetInterpolator( nn_interpolator );
        }
        break;
      case 2:
        {
        double sigma[ImageDimension];
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          sigma[d] = reader->GetOutput()->GetSpacing()[d];
          }
        double alpha = 1.0;

        if( argc > 7 )
          {
          std::vector<RealType> sg = ConvertVector<RealType>( std::string( argv[7] ) );
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            sigma[d] = sg[d];
            }
          }
        if( argc > 8 )
          {
          alpha = static_cast<double>( atof( argv[8] ) );
          }
        g_interpolator->SetParameters( sigma, alpha );

        resampler->SetInterpolator( g_interpolator );
        }
        break;
      case 3:
        {
        switch( arg7 )
          {
          case 'h': default:
            {
            resampler->SetInterpolator( sh_interpolator );
            }
            break;
          case 'c':
            {
            resampler->SetInterpolator( sc_interpolator );
            }
            break;
          case 'l':
            {
            resampler->SetInterpolator( sl_interpolator );
            }
            break;
          case 'w':
            {
            resampler->SetInterpolator( sw_interpolator );
            }
            break;
          case 'b':
            {
            resampler->SetInterpolator( sb_interpolator );
            }
            break;
          }
        }
      case 4:
        {
        if( argc > 7 && atoi( argv[7] ) >= 0 && atoi( argv[7] ) <= 5 )
          {
          bs_interpolator->SetSplineOrder( atoi( argv[7] ) );
          }
        else
          {
          bs_interpolator->SetSplineOrder( 3 );
          }
        resampler->SetInterpolator( bs_interpolator );
        }
        break;
      }
    }
  resampler->SetInput( reader->GetOutput() );
  resampler->SetSize( size );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetOutputDirection( reader->GetOutput()->GetDirection() );
  resampler->SetOutputSpacing( spacing );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ResampleImage( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ResampleImage" );

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

  if( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
             << "outputImage MxNxO [size=1,spacing=0] [interpolate type]" << std::endl;
    std::cout << "  Interpolation type: " << std::endl;
    std::cout << "    0. linear (default)" << std::endl;
    std::cout << "    1. nn " << std::endl;
    std::cout << "    2. gaussian [sigma=imageSpacing] [alpha=1.0]" << std::endl;
    std::cout << "    3. windowedSinc [type = 'c'osine, 'w'elch, 'b'lackman, 'l'anczos, 'h'amming]" << std::endl;
    std::cout << "    4. B-Spline [order=3]" << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      {
      ResampleImage<2>( argc, argv );
      }
      break;
    case 3:
      {
      ResampleImage<3>( argc, argv );
      }
      break;
    case 4:
      {
      ResampleImage<4>( argc, argv );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
