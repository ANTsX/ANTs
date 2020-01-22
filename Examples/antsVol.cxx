#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkImageToVTKImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "vtkCamera.h"
#include "vtkColorTransferFunction.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkImageShiftScale.h"
#include "vtkMultiThreader.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPointData.h"
#include "vtkPNGWriter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSampleFunction.h"
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkSphere.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkWindowToImageFilter.h"

#include "itkMath.h"

#include <vector>
#include <string>

namespace ants
{

int antsVolumetricRendering( itk::ants::CommandLineParser *parser )
{
  constexpr unsigned int ImageDimension = 3;

  typedef float                                      RealType;
  typedef itk::Image<RealType, ImageDimension>       ImageType;
  typedef itk::Image<unsigned int, ImageDimension>   MaskImageType;

  typedef unsigned char                              RgbComponentType;

  typedef itk::RGBPixel<RgbComponentType>            RgbPixelType;
  typedef itk::Image<RgbPixelType, ImageDimension>   RgbImageType;

  typedef itk::RGBAPixel<RgbComponentType>            RgbaPixelType;
  typedef itk::Image<RgbaPixelType, ImageDimension>   RgbaImageType;

  // Read in input image

  ImageType::Pointer inputImage = nullptr;

  itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "input-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = std::string( "" );

    std::vector<RealType> clipPercentage;
    clipPercentage.push_back( 0.0 );
    clipPercentage.push_back( 1.0 );

    if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      inputFile = inputImageOption->GetFunction( 0 )->GetName();
      }
    else if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      inputFile = inputImageOption->GetFunction( 0 )->GetParameter( 0 );
      if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        clipPercentage = parser->ConvertVector<RealType>(
          inputImageOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      }

    ImageType::Pointer readImage = nullptr;
    ReadImage<ImageType>( readImage, inputFile.c_str() );

    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum( 0.0 );
    rescaler->SetOutputMaximum( 1.0 );
    rescaler->SetInput( readImage );

    typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;
    IntensityWindowingImageFilterType::Pointer windower = IntensityWindowingImageFilterType::New();
    windower->SetInput( rescaler->GetOutput() );
    windower->SetWindowMinimum( clipPercentage[0] );
    windower->SetWindowMaximum( clipPercentage[1] );
    windower->SetOutputMinimum( 0.0 );
    windower->SetOutputMaximum( 255.0 );

    inputImage = windower->GetOutput();
    inputImage->Update();
    inputImage->DisconnectPipeline();
    }
  else
    {
    std::cerr << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }

  // Read in input image

  MaskImageType::Pointer maskImage = nullptr;

  itk::ants::CommandLineParser::OptionType::Pointer maskImageOption =
    parser->GetOption( "mask-image" );
  if( maskImageOption && maskImageOption->GetNumberOfFunctions() )
    {
    std::string maskFile = maskImageOption->GetFunction( 0 )->GetName();
    ReadImage<MaskImageType>( maskImage, maskFile.c_str() );
    }

  // Read in the functional overlays and alpha values

  std::vector<RgbImageType::Pointer>    functionalRgbImages;
  std::vector<MaskImageType::Pointer>   functionalMaskImages;

  itk::ants::CommandLineParser::OptionType::Pointer functionalOverlayOption =
    parser->GetOption( "functional-overlay" );
  if( functionalOverlayOption && functionalOverlayOption->GetNumberOfFunctions() )
    {
    for( unsigned int n = 0; n < functionalOverlayOption->GetNumberOfFunctions(); n++ )
      {
      // read RGB and mask image
      std::string rgbFileName = std::string( "" );
      std::string maskFileName = std::string( "" );

      if( functionalOverlayOption->GetFunction( n )->GetNumberOfParameters() == 0 )
        {
        rgbFileName = functionalOverlayOption->GetFunction( n )->GetName();
        }
      else
        {
        rgbFileName = functionalOverlayOption->GetFunction( n )->GetParameter( 0 );
        if( functionalOverlayOption->GetFunction( n )->GetNumberOfParameters() > 1 )
          {
          maskFileName = functionalOverlayOption->GetFunction( n )->GetParameter( 1 );
          }
        }

      typedef itk::ImageFileReader<RgbImageType> RgbReaderType;
      RgbReaderType::Pointer rgbReader = RgbReaderType::New();
      rgbReader->SetFileName( rgbFileName.c_str() );
      try
        {
        rgbReader->Update();
        }
      catch( ... )
        {
        std::cerr << "Error reading RGB file " << rgbFileName << std::endl;
        return EXIT_FAILURE;
        }
      functionalRgbImages.emplace_back(rgbReader->GetOutput() );

      if( ! maskFileName.empty() )
        {
        typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
        MaskReaderType::Pointer maskReader = MaskReaderType::New();
        maskReader->SetFileName( maskFileName.c_str() );
        maskReader->Update();

        functionalMaskImages.emplace_back(maskReader->GetOutput() );
        }
      else
        {
        functionalMaskImages.emplace_back(nullptr );
        }
      }
    }

  // Combine the functional overlays and alpha values

  RgbaImageType::Pointer rgbaImage = RgbaImageType::New();
  rgbaImage->CopyInformation( inputImage );
  rgbaImage->SetRegions( inputImage->GetRequestedRegion() );
  rgbaImage->Allocate();

  itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImage, inputImage->GetRequestedRegion() );

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ImageType::IndexType index = It.GetIndex();
    ImageType::PixelType pixel = It.Get();

    if( maskImage.IsNotNull() && maskImage->GetPixel( index ) == 0 )
      {
      pixel = 0.0;
      }

    // The rgb values are in the range 0,255 but we manipulate alpha values in the range
    // [0,1] since that is what is specified on the command line and simply renormalize
    // to the range [0,255] when setting the voxel.

    RealType currentRed   = static_cast<RealType>( pixel ) / static_cast<RealType>( 255.0 );
    RealType currentGreen = static_cast<RealType>( pixel ) / static_cast<RealType>( 255.0 );
    RealType currentBlue  = static_cast<RealType>( pixel ) / static_cast<RealType>( 255.0 );

    RealType currentAlpha = static_cast<RealType>( pixel ) / static_cast<RealType>( 255.0 );

    for( int i = functionalRgbImages.size() - 1; i >= 0; i-- )
      {
      if( functionalMaskImages[i].IsNotNull() && functionalMaskImages[i]->GetPixel( index ) == 0 )
        {
        continue;
        }

      RgbPixelType rgbPixel = functionalRgbImages[i]->GetPixel( index );

      RealType functionalRed = rgbPixel.GetRed() / static_cast<RealType>( 255.0 );
      RealType functionalGreen = rgbPixel.GetGreen() / static_cast<RealType>( 255.0 );
      RealType functionalBlue = rgbPixel.GetBlue() / static_cast<RealType>( 255.0 );

      if( functionalRed + functionalGreen + functionalBlue > itk::NumericTraits<RealType>::ZeroValue() )
        {
        currentRed   = functionalRed;
        currentGreen = functionalGreen;
        currentBlue  = functionalBlue;
        }
      }

    RgbaPixelType currentColor;
    currentColor.SetRed( static_cast<unsigned char>( currentRed * static_cast<RealType>( 255.0 ) ) );
    currentColor.SetGreen( static_cast<unsigned char>( currentGreen * static_cast<RealType>( 255.0 ) ) );
    currentColor.SetBlue( static_cast<unsigned char>( currentBlue * static_cast<RealType>( 255.0 ) ) );
    currentColor.SetAlpha( static_cast<unsigned char>( currentAlpha * static_cast<RealType>( 255.0 ) ) );

    rgbaImage->SetPixel( index, currentColor );
    }

  // Get display options

  float magnificationFactor = 3.0;

  std::vector<float> rotationAnglesInDegrees;
  rotationAnglesInDegrees.push_back( 0.0 );
  rotationAnglesInDegrees.push_back( 0.0 );
  rotationAnglesInDegrees.push_back( 0.0 );

  std::vector<float> backgroundColor;
  backgroundColor.push_back( 255.0 );
  backgroundColor.push_back( 255.0 );
  backgroundColor.push_back( 255.0 );

  std::string screenCaptureFileName = std::string( "" );
  bool writeScreenCaptureToFile = false;

  itk::ants::CommandLineParser::OptionType::Pointer displayOption = parser->GetOption( "display" );
  if( displayOption && displayOption->GetNumberOfFunctions() )
    {
    screenCaptureFileName = displayOption->GetFunction( 0 )->GetName();

    if( strcmp( screenCaptureFileName.c_str(), "false" ) == 0 ||
        strcmp( screenCaptureFileName.c_str(), "0" ) == 0 )
      {
      // do not render and exit
      return EXIT_SUCCESS;
      }

    std::size_t position = screenCaptureFileName.find( "png" );
    if( position == std::string::npos )
      {
      screenCaptureFileName.clear();
      }
    else
      {
      writeScreenCaptureToFile = true;
      std::cout << "Writing screenshot to image file " << screenCaptureFileName << "." << std::endl;
      }

    if( displayOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      magnificationFactor = parser->Convert<float>(
        displayOption->GetFunction( 0 )->GetParameter( 0 ) );
      }

    if( displayOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      rotationAnglesInDegrees = parser->ConvertVector<float>(
        displayOption->GetFunction( 0 )->GetParameter( 1 ) );
      }

    if( displayOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      backgroundColor = parser->ConvertVector<float>(
        displayOption->GetFunction( 0 )->GetParameter( 2 ) );
      if( backgroundColor.size() == 1 )
        {
        backgroundColor.push_back( backgroundColor[0] );
        backgroundColor.push_back( backgroundColor[0] );
        }
      }
    }

  // Set up rendering window

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( backgroundColor[0] / static_cast<RealType>( 255.0 ),
                           backgroundColor[1] / static_cast<RealType>( 255.0 ),
                           backgroundColor[2] / static_cast<RealType>( 255.0 ) );

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 301, 300 );

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderWindow->Render();

  // Turn off multi-threading?

  vtkSmartPointer<vtkMultiThreader> multiThreader = vtkSmartPointer<vtkMultiThreader>::New();
  multiThreader->SetGlobalMaximumNumberOfThreads( 1 );

  // Do volumetric rendering

  typedef itk::ImageToVTKImageFilter<RgbaImageType> ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput( rgbaImage );
  connector->Update();

  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  imageData->ShallowCopy( connector->GetOutput() );

  vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
  volumeMapper->SetInputData( imageData );

  vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
  compositeOpacity->AddPoint(   0.0, 0.0  );
  compositeOpacity->AddPoint( 124.0, 0.25  );
  compositeOpacity->AddPoint( 125.0, 0.5 );
  compositeOpacity->AddPoint( 255.0, 1.0 );

  vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationType( VTK_LINEAR_INTERPOLATION );
  volumeProperty->SetScalarOpacity( 0, compositeOpacity );
  volumeProperty->IndependentComponentsOff();

  vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
  volume->SetMapper( volumeMapper );
  volume->SetProperty( volumeProperty );

  renderer->AddViewProp( volume );
  renderer->ResetCamera();

  renderer->GetActiveCamera()->Azimuth( rotationAnglesInDegrees[0] );
  renderer->GetActiveCamera()->Elevation( rotationAnglesInDegrees[1] );
  renderer->GetActiveCamera()->Roll( rotationAnglesInDegrees[2] );
  renderer->GetActiveCamera()->Zoom( magnificationFactor );

  renderWindow->Render();

  // Screenshot

  if( writeScreenCaptureToFile )
    {
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput( renderWindow );
    windowToImageFilter->SetScale( magnificationFactor );
    windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName( screenCaptureFileName.c_str() );
    writer->SetInputConnection( windowToImageFilter->GetOutputPort() );
    writer->Write();
    }
  else
    {
    renderWindowInteractor->Start();
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Main input image for 3-D rendering." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input-image" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "imageFilename" );
    option->SetUsageOption( 1, "[imageFilename,<lowerClipPercentagexupperClipPercentage=0.0x1.0>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Mask for determining volumetric rendering of main input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "A functional overlay can be specified using an rgb image." )
      + std::string( "Note that more than one functional overlays can " )
      + std::string( "be rendered, the order in which they are specified " )
      + std::string( "on the command line matters, and rgb images are " )
      + std::string( "assumed to be unsigned char [0,255].  An optional mask " )
      + std::string( "can also be provided." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "functional-overlay" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "rgbImageFileName" );
    option->SetUsageOption( 1, "[rgbImageFileName,rgbMaskFileName]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Display output volume rendering in VTK window.  Rotation " )
      + std::string( "angles are in degrees and the default background color " )
      + std::string( "is white (255x255x255).  Note that the filename, to be " )
      + std::string( "considered such, must have a \"png\" extension.  If the " )
      + std::string( "filename is omitted in the third usage option, then the " )
      + std::string( "window is displayed." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "display" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "doWindowDisplay" );
    option->SetUsageOption( 1, "filename" );
    option->SetUsageOption( 2, "<filename>[<magnificationFactor=1>,<AzimuthxElevationxRoll=0x0x0>,<backgroundColor=255x255x255>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    parser->AddOption( option );
    }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsVol( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsVol" );

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

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "Produce a 3-D volume rendering with optional RGB overlay. " );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc == 1 )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_FAILURE;
    }
  else if( parser->GetOption( "help" )->GetFunction() && parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' )->GetFunction() && parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Get dimensionality

  itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "input-image" );

  if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
    {
    std::string inputFile;
    if( imageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      inputFile = imageOption->GetFunction( 0 )->GetName();
      }
    else if( imageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      inputFile = imageOption->GetFunction( 0 )->GetParameter( 0 );
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        inputFile.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
    unsigned int dimension = imageIO->GetNumberOfDimensions();

    if( dimension == 3 )
      {
      antsVolumetricRendering( parser );
      }
    else
      {
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "Input image not specified.  See help menu." << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
} // namespace ants
