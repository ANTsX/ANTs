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

template <class TValue>
TValue Convert( std::string optionString )
{
  TValue             value;
  std::istringstream iss( optionString );

  iss >> value;
  return value;
}

template <class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string        element = optionString.substr( 0, crosspos );
    TValue             value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos );
        }
      std::istringstream iss( element );
      iss >> value;
      values.push_back( value );
      }
    }
  return values;
}

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
      std::cerr << "Invalid spacing." << std::endl;
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
      std::cerr << "Invalid size." << std::endl;
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
        resampler->SetInterpolator( interpolator );
        break;
      case 1:
        resampler->SetInterpolator( nn_interpolator );
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
            resampler->SetInterpolator( sh_interpolator );
            break;
          case 'c':
            resampler->SetInterpolator( sc_interpolator );
            break;
          case 'l':
            resampler->SetInterpolator( sl_interpolator );
            break;
          case 'w':
            resampler->SetInterpolator( sw_interpolator );
            break;
          case 'b':
            resampler->SetInterpolator( sb_interpolator );
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
        break;
        }
      }
    }
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetSize( size );
  resampler->SetOutputDirection( reader->GetOutput()->GetDirection() );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
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
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      ResampleImage<2>( argc, argv );
      break;
    case 3:
      ResampleImage<3>( argc, argv );
      break;
    case 4:
      ResampleImage<4>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
