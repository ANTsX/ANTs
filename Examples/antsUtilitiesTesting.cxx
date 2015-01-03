#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include <stdio.h>

#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"

#include <string>
#include <vector>

namespace ants
{
int antsUtilitiesTesting( std::vector<std::string> args, std::ostream* itkNotUsed( out_stream ) )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsUtilitiesTesting" );

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

  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " testImage whichMetric[CCorMI] lowerScalexUpperScalexNumberOfScaleSamples numberOfRotationSamples doPrintMatchingFileName transformName setOfTrainingImages" << std::endl;
    std::cerr << "Notes:  if transformName='none', no transform is printed" << std::endl;
    std::cerr << "Example call:  " << std::endl;
    std::cerr << "   " << argv[0] << "test.nii.gz CC 0.25x3x10 10 testTransform.txt 0 training*.nii.gz" << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

  const unsigned int ImageDimension = 2;

  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType> MIMetricType;
  typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType, ImageType> CCMetricType;

  CCMetricType::RadiusType ccRadius;
  ccRadius.Fill( 4 );

  unsigned int metricType = 0;
  if( strcmp( argv[2], "CC" ) == 0 )
    {
    metricType = 1;
    }
  else if( strcmp( argv[2], "MI" ) == 0 )
    {
    metricType = 2;
    }
  else
    {
    std::cerr << "Unrecognized metric " << argv[2] << std::endl;
    return EXIT_FAILURE;
    }

  // read in test image

  ImageType::Pointer testImage = ImageType::New();
  ReadImage<ImageType>( testImage, argv[1] );

  // read in training image file names

  std::vector<std::string> trainingImageFileNames;
  for( int n = 7; n < argc; n++ )
    {
    trainingImageFileNames.push_back( std::string( argv[n] ) );
    }

//   // Get scale parameters

  std::vector<float> scaleParameters = ConvertVector<float>( std::string( argv[3] ) );
  if( scaleParameters.size() != 3 )
    {
    std::cerr << "The scale parameters were improperly specified.  See usage." << std::endl;
    return EXIT_FAILURE;
    }
  float scaleLowerBoundLog = vcl_log( scaleParameters[0] );
  float scaleUpperBoundLog = vcl_log( scaleParameters[1] );
  unsigned int scaleNumberOfSamples = static_cast<unsigned int>( scaleParameters[2] );
  float scaleDelta = ( scaleUpperBoundLog - scaleLowerBoundLog ) / static_cast<float>( scaleNumberOfSamples - 1 );

//   // Get rotation sampling resolution

  unsigned int rotationNumberOfSamples = static_cast<unsigned int>( atoi( argv[4] ) );
  float rotationDelta = ( 2.0 * vnl_math::pi - 0.0 ) / static_cast<float>( rotationNumberOfSamples - 1 );

//   // Now go through the rotations + scalings to find the optimal pose.

  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;

  float optimalMetricValue = itk::NumericTraits<float>::max();
  AffineTransformType::Pointer optimalTransform = ITK_NULLPTR;
  unsigned int optimalMetricIndex = 0;

  for( unsigned int n = 0; n < trainingImageFileNames.size(); n++ )
    {
    ImageType::Pointer trainingImage = ImageType::New();
    ReadImage<ImageType>( trainingImage, trainingImageFileNames[n].c_str() );

    // Initialize centered transform (based on the center of the image)

    AffineTransformType::Pointer initialTransform = AffineTransformType::New();

    typedef itk::CenteredTransformInitializer<AffineTransformType, ImageType, ImageType> TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform( initialTransform );
    initializer->SetFixedImage( trainingImage );
    initializer->SetMovingImage( testImage );
    initializer->GeometryOn();

    initializer->InitializeTransform();

    for( float angle = 0.0; angle <= 2.0 * vnl_math::pi; angle += rotationDelta )
      {
      AffineTransformType::MatrixType rotationMatrix;
      rotationMatrix( 0, 0 ) = rotationMatrix( 1, 1 ) = vcl_cos( angle );
      rotationMatrix( 1, 0 ) = vcl_sin( angle );
      rotationMatrix( 0, 1 ) = -rotationMatrix( 1, 0 );

      for( float scaleLog = scaleLowerBoundLog; scaleLog <= scaleUpperBoundLog; scaleLog += scaleDelta )
        {
        float scale = vcl_exp( scaleLog );

        AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        affineTransform->SetCenter( initialTransform->GetCenter() );
        affineTransform->SetTranslation( initialTransform->GetTranslation() );
        affineTransform->SetMatrix( initialTransform->GetMatrix() * rotationMatrix * scale );

//         typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResamplerType;
//         ResamplerType::Pointer resampleFilter = ResamplerType::New();
//         resampleFilter->SetInput( testImage );
//         resampleFilter->SetOutputParametersFromImage( trainingImage );
//         resampleFilter->SetTransform( affineTransform );
//         resampleFilter->SetDefaultPixelValue( 255 );
//         resampleFilter->Update();
//         WriteImage<ImageType>( resampleFilter->GetOutput(), "test.nii.gz" );

        float metricValue = 0;

        if( metricType == 1 )
          {
          CCMetricType::Pointer metric = CCMetricType::New();
          metric->SetFixedImage( trainingImage );
          metric->SetMovingImage( testImage );
          metric->SetRadius( ccRadius );
          metric->SetMovingTransform( affineTransform );
          metric->SetVirtualDomainFromImage( trainingImage );
          metric->Initialize();

          try
            {
            metricValue = metric->GetValue();
            }
          catch(...)
            {
            continue;
            }
          }
        else if( metricType == 2 )
          {
          MIMetricType::Pointer metric = MIMetricType::New();
          metric->SetFixedImage( trainingImage );
          metric->SetMovingImage( testImage );
          metric->SetNumberOfHistogramBins( 20 );
          metric->SetMovingTransform( affineTransform );
          metric->SetVirtualDomainFromImage( trainingImage );
          metric->Initialize();

          try
            {
            metricValue = metric->GetValue();
            }
          catch(...)
            {
            continue;
            }
          }
        if( metricValue < optimalMetricValue )
          {
          optimalMetricValue = metricValue;
          optimalTransform = affineTransform;
          optimalMetricIndex = n;
          }
        }
      }
    }

  if( strcmp( argv[6], "none" ) != 0 )
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetInput( optimalTransform );
    transformWriter->SetFileName( argv[6] );
    transformWriter->Update();
    }

  if( static_cast<bool>( atoi( argv[5] ) ) )
    {
    std::cout << trainingImageFileNames[optimalMetricIndex] << std::endl;
    }
  std::cout << optimalMetricValue << std::endl;

  return EXIT_SUCCESS;
}

} // namespace ants
