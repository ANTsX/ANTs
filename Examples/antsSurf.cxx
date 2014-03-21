#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteImage.h"

#include "itkAffineTransform.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkActor.h"
#include "vtkCallbackCommand.h"
#include "vtkExtractEdges.h"
#include "vtkGraphicsFactory.h"
#include "vtkMarchingCubes.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkTriangleFilter.h"
#include "vtkUnsignedCharArray.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkPNGWriter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkWindowToImageFilter.h"

#include "vnl/vnl_math.h"

#include <vector>
#include <string>

namespace ants
{

float CalculateGenus( vtkPolyData *mesh, bool verbose )
{
  vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputData( mesh );
  extractEdges->Update();

  float numberOfEdges = static_cast<float>( extractEdges->GetOutput()->GetNumberOfLines() );
  float numberOfVertices = static_cast<float>( mesh->GetNumberOfPoints() );
  float numberOfFaces = static_cast<float>( mesh->GetNumberOfPolys() );

  float genus = 0.5 * ( 2.0 - numberOfVertices + numberOfEdges - numberOfFaces );

  if( verbose )
    {
    std::cout << "Genus = "  << genus << std::endl;
    std::cout << "  number of vertices = "  << numberOfVertices << std::endl;
    std::cout << "  number of edges = "  << numberOfEdges << std::endl;
    std::cout << "  number of faces = "  << numberOfFaces << std::endl;
    }

  return genus;
}

void Display( vtkPolyData *vtkMesh, const std::vector<float> rotationAngleInDegrees,
              const std::vector<float> backgroundColor,
              const std::string screenCaptureFileName )
{
  vtkSmartPointer<vtkGraphicsFactory> graphicsFactory =
    vtkSmartPointer<vtkGraphicsFactory>::New();
  graphicsFactory->SetOffScreenOnlyMode( false );
  graphicsFactory->SetUseMesaClasses( 1 );

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData( vtkMesh );
  mapper->ScalarVisibilityOn();

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper( mapper );
  actor->GetProperty()->SetInterpolationToFlat();
  actor->GetProperty()->ShadingOff();
  actor->GetProperty()->SetSpecular( 0.0 );
  actor->GetProperty()->SetSpecularPower( 0 );
  actor->RotateX( rotationAngleInDegrees[0] );
  actor->RotateY( rotationAngleInDegrees[1] );
  actor->RotateZ( rotationAngleInDegrees[2] );

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( backgroundColor[0] / 255.0, backgroundColor[1] / 255.0, backgroundColor[2] / 255.0 );

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );

  vtkSmartPointer<vtkCallbackCommand> callback = vtkSmartPointer<vtkCallbackCommand>::New();
  renderer->AddObserver( vtkCommand::KeyPressEvent, callback );

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderer->AddActor( actor );
  renderWindow->Render();

  if( screenCaptureFileName.empty() )
    {
    renderWindowInteractor->Start();
    }
  else
    {
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput( renderWindow );
    windowToImageFilter->SetMagnification( 4 );
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName( screenCaptureFileName.c_str()  );
    writer->SetInputConnection( windowToImageFilter->GetOutputPort() );
    writer->Write();
    }
}

int antsSurf( itk::ants::CommandLineParser *parser )
{
  const unsigned int ImageDimension = 3;

  typedef float                                      RealType;
  typedef itk::Image<RealType, ImageDimension>       ImageType;
  typedef itk::Image<int, ImageDimension>            MaskImageType;

  typedef unsigned char                              RgbComponentType;
  typedef itk::RGBPixel<RgbComponentType>            RgbPixelType;
  typedef itk::Image<RgbPixelType, ImageDimension>   RgbImageType;

  // Read in input surface image

  ImageType::Pointer inputImage = NULL;

  RealType defaultColorRed = 255.0;
  RealType defaultColorGreen = 255.0;
  RealType defaultColorBlue = 255.0;

  itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "surface-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      std::string inputFile = inputImageOption->GetFunction( 0 )->GetName();
      ReadImage<ImageType>( inputImage, inputFile.c_str() );
      }
    else
      {
      std::string inputFile = inputImageOption->GetFunction( 0 )->GetParameter( 0 );
      ReadImage<ImageType>( inputImage, inputFile.c_str() );

      if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        std::vector<RealType> defaultColors = parser->ConvertVector<RealType>(
          inputImageOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( defaultColors.size() == 1 )
          {
          defaultColorRed = defaultColors[0];
          defaultColorGreen = defaultColors[0];
          defaultColorBlue = defaultColors[0];
          }
        else if( defaultColors.size() == 3 )
          {
          defaultColorRed = defaultColors[0];
          defaultColorGreen = defaultColors[1];
          defaultColorBlue = defaultColors[2];
          }
        else
          {
          std::cerr << "Incorrect color format specified." << std::endl;
          return EXIT_FAILURE;
          }
        }
      }
    }
  else
    {
    std::cerr << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }


  // There's a reorientation issue between itk image physical space and the mesh space
  // for which we have to account.  See
  // http://www.vtk.org/pipermail/vtkusers/2011-July/068595.html
  //   and
  // http://www.vtk.org/Wiki/VTK/ExamplesBoneYard/Cxx/VolumeRendering/itkVtkImageConvert

  typedef itk::AffineTransform<RealType> RigidTransformType;
  RigidTransformType::Pointer meshToItkImageTransform = RigidTransformType::New();
  RigidTransformType::OutputVectorType offset;
  offset[0] = -inputImage->GetOrigin()[0];
  offset[1] = -inputImage->GetOrigin()[1];
  offset[2] = -inputImage->GetOrigin()[2];
  RigidTransformType::MatrixType matrix;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      matrix( i, j ) = inputImage->GetDirection()( i, j );
      }
    }
  meshToItkImageTransform->SetMatrix( matrix );
  meshToItkImageTransform->SetOffset( offset );

  // Get anti-alias RMSE parameter

  RealType antiAliasRmseParameter = 0.03;

  itk::ants::CommandLineParser::OptionType::Pointer antiAliasRmseOption =
    parser->GetOption( "anti-alias-rmse" );
  if( antiAliasRmseOption && antiAliasRmseOption->GetNumberOfFunctions() )
    {
    antiAliasRmseParameter = parser->Convert<RealType>( antiAliasRmseOption->GetFunction( 0 )->GetName() );
    }

  typedef itk::AntiAliasBinaryImageFilter<ImageType, ImageType> AntiAliasFilterType;
  AntiAliasFilterType::Pointer antiAlias = AntiAliasFilterType::New();
  antiAlias->SetMaximumRMSError( antiAliasRmseParameter );
  antiAlias->SetInput( inputImage );
  antiAlias->Update();

  // Reconstruct binary surface.

  typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput( antiAlias->GetOutput() );
  connector->Update();

  vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
  marchingCubes->SetInputData( connector->GetOutput() );
  marchingCubes->ComputeScalarsOff();
  marchingCubes->ComputeGradientsOff();
  marchingCubes->SetNumberOfContours( 1 );
  marchingCubes->SetValue( 0, 0.0 );
  marchingCubes->Update();

  vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter =
    vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  connectivityFilter->SetExtractionModeToLargestRegion();
  connectivityFilter->SetInputData( marchingCubes->GetOutput() );
  connectivityFilter->Update();

  vtkSmartPointer<vtkTriangleFilter> triangularizer = vtkSmartPointer<vtkTriangleFilter>::New();
  triangularizer->SetInputData( connectivityFilter->GetOutput() );
  triangularizer->Update();

  // Smooth mesh

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> meshSmoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  meshSmoother->SetInputData( triangularizer->GetOutput() );
  meshSmoother->SetNumberOfIterations( 25 );
  meshSmoother->BoundarySmoothingOff();
  meshSmoother->FeatureEdgeSmoothingOff();
  meshSmoother->SetFeatureAngle( 120.0 );
  meshSmoother->SetPassBand( 0.01 );
  meshSmoother->NonManifoldSmoothingOn();
  meshSmoother->NormalizeCoordinatesOn();
  meshSmoother->Update();

  vtkPolyData *vtkMesh = meshSmoother->GetOutput();

  CalculateGenus( vtkMesh, true );

  // Add the functional overlays

  std::vector<RgbImageType::Pointer>    functionalRgbImages;
  std::vector<MaskImageType::Pointer>   functionalMaskImages;
  std::vector<RealType>                 functionalAlphaValues;

  itk::ants::CommandLineParser::OptionType::Pointer functionalOverlayOption =
    parser->GetOption( "functional-overlay" );
  if( functionalOverlayOption && functionalOverlayOption->GetNumberOfFunctions() )
    {
    for( unsigned int n = 0; n < functionalOverlayOption->GetNumberOfFunctions(); n++ )
      {
      if( functionalOverlayOption->GetFunction( n )->GetNumberOfParameters() < 2 )
        {
        std::cerr << "Error:  each functional overlay must have an RGB image and mask."
          << "See help menu." << std::endl;
        return EXIT_FAILURE;
        }

      // read RGB image

      std::string rgbFileName = functionalOverlayOption->GetFunction( n )->GetParameter( 0 );

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
      functionalRgbImages.push_back( rgbReader->GetOutput() );

      // read mask

      std::string maskFileName = functionalOverlayOption->GetFunction( n )->GetParameter( 1 );

      typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( maskFileName.c_str() );
      try
        {
        maskReader->Update();
        }
      catch( ... )
        {
        std::cerr << "Error reading mask file " << maskFileName << std::endl;
        return EXIT_FAILURE;
        }
      functionalMaskImages.push_back( maskReader->GetOutput() );

      if( functionalOverlayOption->GetFunction( n )->GetNumberOfParameters() > 2 )
        {
        RealType alpha = parser->Convert<RealType>(
          functionalOverlayOption->GetFunction( n )->GetParameter( 2 ) );
        functionalAlphaValues.push_back( alpha );
        }
      else
        {
        functionalAlphaValues.push_back( 1.0 );
        }
      }
    }

  // Reset mesh points to physical space of ITK images

  vtkPoints* meshPoints = vtkMesh->GetPoints();
  int        numberOfPoints = meshPoints->GetNumberOfPoints();

  for( int n = 0; n < numberOfPoints; n++ )
    {
    RigidTransformType::InputPointType inputTransformPoint;
    RigidTransformType::OutputPointType outputTransformPoint;

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      inputTransformPoint[d] = meshPoints->GetPoint( n )[d];
      }
    outputTransformPoint = meshToItkImageTransform->TransformPoint( inputTransformPoint );

    meshPoints->SetPoint( n, outputTransformPoint[0], outputTransformPoint[1], outputTransformPoint[2] );
    }


  // Do the painting
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents( 3 );   // R, G, B, and alpha components
  colors->SetName( "Colors" );

  for( int n = 0; n < numberOfPoints; n++ )
    {
    ImageType::IndexType index;
    ImageType::PointType imagePoint;

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      imagePoint[d] = meshPoints->GetPoint( n )[d];
      }

    RealType currentRed   = defaultColorRed / 255.0;
    RealType currentGreen = defaultColorGreen / 255.0;
    RealType currentBlue  = defaultColorBlue / 255.0;
    RealType currentAlpha = 1.0;

    for( unsigned int i = 0; i < functionalAlphaValues.size(); i++ )
      {
      bool isInsideImage = functionalMaskImages[i]->TransformPhysicalPointToIndex( imagePoint, index );

      if( isInsideImage && functionalMaskImages[i]->GetPixel( index ) != 0 )
        {
        // http://stackoverflow.com/questions/726549/algorithm-for-additive-color-mixing-for-rgb-values
        // or
        // http://en.wikipedia.org/wiki/Alpha_compositing

        RgbPixelType rgbPixel = functionalRgbImages[i]->GetPixel( index );

        RealType functionalRed = rgbPixel.GetRed() / 255.0;
        RealType functionalGreen = rgbPixel.GetGreen() / 255.0;
        RealType functionalBlue = rgbPixel.GetBlue() / 255.0;
        RealType functionalAlpha = functionalAlphaValues[i];

        RealType backgroundRed   = currentRed;
        RealType backgroundGreen = currentGreen;
        RealType backgroundBlue  = currentBlue;
        RealType backgroundAlpha = currentAlpha;

        currentAlpha = 1.0 - ( 1.0 - functionalAlpha ) * ( 1.0 - backgroundAlpha );
        currentRed   = functionalRed   * functionalAlpha / currentAlpha + backgroundRed   * backgroundAlpha * ( 1.0 - functionalAlpha ) / currentAlpha;
        currentGreen = functionalGreen * functionalAlpha / currentAlpha + backgroundGreen * backgroundAlpha * ( 1.0 - functionalAlpha ) / currentAlpha;
        currentBlue  = functionalBlue  * functionalAlpha / currentAlpha + backgroundBlue  * backgroundAlpha * ( 1.0 - functionalAlpha ) / currentAlpha;
        }
      }

    unsigned char currentColor[3];
    currentColor[0] = static_cast<unsigned char>( currentRed   * 255.0 );
    currentColor[1] = static_cast<unsigned char>( currentGreen * 255.0 );
    currentColor[2] = static_cast<unsigned char>( currentBlue  * 255.0 );

    colors->InsertNextTupleValue( currentColor );
    }
  vtkMesh->GetPointData()->SetScalars( colors );

  // Inflation

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> inflater =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();

  itk::ants::CommandLineParser::OptionType::Pointer inflationOption = parser->GetOption( "inflation" );
  if( inflationOption && inflationOption->GetNumberOfFunctions() )
    {
    unsigned int numberOfIterations = 0;
    if( inflationOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      numberOfIterations = parser->Convert<unsigned int>( inflationOption->GetFunction( 0 )->GetName() );
      }
    else
      {
      numberOfIterations = parser->Convert<unsigned int>( inflationOption->GetFunction( 0 )->GetParameter( 0 ) );
      }

    inflater->SetInputData( vtkMesh );
    inflater->SetNumberOfIterations( numberOfIterations );
    inflater->BoundarySmoothingOn();
    inflater->FeatureEdgeSmoothingOff();
    inflater->SetFeatureAngle( 180.0 );
    inflater->SetEdgeAngle( 180.0 );
    inflater->SetPassBand( 0.001 );
    inflater->NonManifoldSmoothingOn();
    inflater->NormalizeCoordinatesOff();
    inflater->Update();

    vtkMesh = inflater->GetOutput();
    }

  // Write the vtk mesh to file.

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string outputFile = outputOption->GetFunction( 0 )->GetName();

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInputData( vtkMesh );
    writer->SetFileName( outputFile.c_str() );
    writer->Write();
    }
  else
    {
    std::cerr << "No output filename specified." << std::endl;
    return EXIT_FAILURE;
    }

  // Display vtk mesh

  itk::ants::CommandLineParser::OptionType::Pointer displayOption = parser->GetOption( "display" );
  if( displayOption && displayOption->GetNumberOfFunctions() )
    {
    std::vector<float> rotationAnglesInDegrees;
    rotationAnglesInDegrees.push_back( 0.0 );
    rotationAnglesInDegrees.push_back( 0.0 );
    rotationAnglesInDegrees.push_back( 0.0 );

    std::vector<float> backgroundColor;
    backgroundColor.push_back( 255.0 );
    backgroundColor.push_back( 255.0 );
    backgroundColor.push_back( 255.0 );

    std::string screenCaptureFileName = std::string( "" );
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
      std::cout << "Writing surface to image file " << screenCaptureFileName << "." << std::endl;
      }

    if( displayOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      Display( vtkMesh, rotationAnglesInDegrees, backgroundColor, screenCaptureFileName );
      }
    else
      {
      if( displayOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        rotationAnglesInDegrees = parser->ConvertVector<float>(
          displayOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      if( displayOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        backgroundColor = parser->ConvertVector<float>(
          displayOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( backgroundColor.size() == 1 )
          {
          backgroundColor.push_back( backgroundColor[0] );
          backgroundColor.push_back( backgroundColor[0] );
          }
        }

      Display( vtkMesh, rotationAnglesInDegrees, backgroundColor, screenCaptureFileName );
      }
    }

  // Clean up
  return EXIT_SUCCESS;
}



void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Main input binary image for 3-D rendering.  One can also " )
      + std::string( "set a default color value in the range [0,255]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "surface-image" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "surfaceImageFilename" );
    option->SetUsageOption( 1, "[surfaceImageFilename,<defaultColor=255x255x255>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "A functional overlay can be specified using both " )
      + std::string( "and rgb image and a mask specifying where that " )
      + std::string( "rgb image should be applied.  Both images must " )
      + std::string( "have the same image geometry as the input image. " )
      + std::string( "Optionally, an alpha parameter can be specified." )
      + std::string( "Note that more than one functional overlays can " )
      + std::string( "be rendered, the order in which they are specified " )
      + std::string( "on the command line matters, and rgb images are " )
      + std::string( "assumed to be unsigned char [0,255]." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "functional-overlay" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "[rgbImageFileName,maskImageFileName,<alpha=1>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Anti-alias maximum RMSE parameter for surface reconstruction " )
      + std::string( "(default = 0.03)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "anti-alias-rmse" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "value" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Perform inflation of the mesh." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "inflation" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "numberOfIterations" );
    option->SetUsageOption( 1, "[numberOfIterations]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Display output surface function in VTK window.  Rotation " )
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
    option->SetUsageOption( 2, "<filename>[rotateXxrotateYxrotateZ,<backgroundColor=255x255x255>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output is a vtk polydata file which can be viewed in programs " )
      + std::string( "such as Paraview." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "surface.vtk" );
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
int antsSurf( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsSurf" );

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

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "Produce a 3-D surface rendering with optional RGB overlay." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

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

  std::string filename;

  itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "surface-image" );
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
        inputFile.c_str(), itk::ImageIOFactory::ReadMode );
    unsigned int dimension = imageIO->GetNumberOfDimensions();

    if( dimension == 3 )
      {
      antsSurf( parser );
      }
    else
      {
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cout << "Input surface image not specified." << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
} // namespace ants
