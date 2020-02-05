#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkAffineTransform.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkSTLReader.h"
#include "vtkSTLWriter.h"
#include "vtkPLYReader.h"
#include "vtkPLYWriter.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkActor.h"
#include "vtkCallbackCommand.h"
#include "vtkExtractEdges.h"
#include "vtkGraphicsFactory.h"
#include "vtkImageData.h"
#include "vtkImageStencil.h"
#include "vtkLookupTable.h"
#include "vtkMarchingCubes.h"
#include "vtkMetaImageWriter.h"
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
#include "vtkPolyDataReader.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkPNGWriter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkScalarBarActor.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTextProperty.h"
#include "vtkWindowToImageFilter.h"

#include "itkMath.h"

#include <vector>
#include <string>

namespace ants
{

float CalculateGenus( vtkPolyData *mesh, bool verbose )
{
  vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputData( mesh );
  extractEdges->Update();

  auto numberOfEdges = static_cast<float>( extractEdges->GetOutput()->GetNumberOfLines() );
  auto numberOfVertices = static_cast<float>( mesh->GetNumberOfPoints() );
  auto numberOfFaces = static_cast<float>( mesh->GetNumberOfPolys() );

  float genus = 0.5f * ( 2.0f - numberOfVertices + numberOfEdges - numberOfFaces );

  if( verbose )
    {
    std::cout << "Genus = "  << genus << std::endl;
    std::cout << "  number of vertices = "  << numberOfVertices << std::endl;
    std::cout << "  number of edges = "  << numberOfEdges << std::endl;
    std::cout << "  number of faces = "  << numberOfFaces << std::endl;
    }

  return genus;
}

void Display( vtkPolyData *vtkMesh,
              const std::vector<float> rotationAngleInDegrees,
              const std::vector<float> backgroundColor,
              const std::string screenCaptureFileName,
              const bool renderScalarBar = false,
              vtkLookupTable *scalarBarLookupTable = nullptr,
              const std::string scalarBarTitle = std::string( "" ),
              unsigned int scalarBarNumberOfLabels = 5,
              unsigned int scalarBarWidthInPixels = 0,
              unsigned int scalarBarHeightInPixels = 0
              )
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
  renderer->SetBackground( backgroundColor[0] / 255.0f, backgroundColor[1] / 255.0f, backgroundColor[2] / 255.0f );

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );

  vtkSmartPointer<vtkCallbackCommand> callback = vtkSmartPointer<vtkCallbackCommand>::New();
  renderer->AddObserver( vtkCommand::KeyPressEvent, callback );

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderer->AddActor( actor );

  if( renderScalarBar )
    {
    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable( scalarBarLookupTable );
    scalarBar->SetTitle( scalarBarTitle.c_str() );
    scalarBar->SetMaximumNumberOfColors( 256 );
    scalarBar->SetNumberOfLabels( scalarBarNumberOfLabels );
    scalarBar->SetLabelFormat( "%.2g" );

    if( scalarBarWidthInPixels > 0 && scalarBarHeightInPixels > 0 )
      {
      if( scalarBarWidthInPixels > scalarBarHeightInPixels )
        {
        scalarBar->SetOrientationToHorizontal();
        }
      else
        {
        scalarBar->SetOrientationToVertical();
        }
      scalarBar->SetMaximumWidthInPixels( scalarBarWidthInPixels );
      scalarBar->SetMaximumHeightInPixels( scalarBarHeightInPixels );
      }

    vtkSmartPointer<vtkTextProperty> titleTextProperty = vtkSmartPointer<vtkTextProperty>::New();
    titleTextProperty->ItalicOff();
    titleTextProperty->BoldOn();
    titleTextProperty->SetColor( 0.0, 0.0, 0.0 );
    titleTextProperty->SetJustificationToCentered();
    // titleTextProperty->SetFontSize( 50 );

    scalarBar->SetTitleTextProperty( titleTextProperty );

    vtkSmartPointer<vtkTextProperty> labelTextProperty = vtkSmartPointer<vtkTextProperty>::New();
    labelTextProperty->ItalicOff();
    labelTextProperty->BoldOff();
    labelTextProperty->SetColor( 0.0, 0.0, 0.0 );
    // labelTextProperty->SetFontSize( 5 );

    scalarBar->SetLabelTextProperty( labelTextProperty );

    scalarBar->VisibilityOn();

    renderer->AddActor2D( scalarBar );
    }

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
    windowToImageFilter->SetScale( 5 );
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName( screenCaptureFileName.c_str()  );
    writer->SetInputConnection( windowToImageFilter->GetOutputPort() );
    writer->Write();
    }

}

int antsImageToSurface( itk::ants::CommandLineParser *parser )
{
  constexpr unsigned int ImageDimension = 3;

  typedef float                                      RealType;
  typedef itk::Image<RealType, ImageDimension>       ImageType;
  typedef itk::Image<int, ImageDimension>            MaskImageType;

  typedef unsigned char                              RgbComponentType;
  typedef itk::RGBPixel<RgbComponentType>            RgbPixelType;
  typedef itk::Image<RgbPixelType, ImageDimension>   RgbImageType;

  ImageType::PointType zeroOrigin;
  zeroOrigin.Fill( 0.0 );

  // Read in input surface image

  ImageType::Pointer inputImage = nullptr;

  RealType defaultColorRed = 255.0;
  RealType defaultColorGreen = 255.0;
  RealType defaultColorBlue = 255.0;
  RealType defaultAlpha = 1.0;

  itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "surface-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      std::string inputFile = inputImageOption->GetFunction( 0 )->GetName();
      ReadImage<ImageType>( inputImage, inputFile.c_str() );
      inputImage->SetOrigin( zeroOrigin );
      }
    else
      {
      std::string inputFile = inputImageOption->GetFunction( 0 )->GetParameter( 0 );
      ReadImage<ImageType>( inputImage, inputFile.c_str() );
      inputImage->SetOrigin( zeroOrigin );

      if( inputImageOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        std::vector<RealType> defaultColors = parser->ConvertVector<RealType>(
          inputImageOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( defaultColors.size() == 1 )
          {
          defaultColorRed = defaultColors[0];
          defaultColorGreen = defaultColors[0];
          defaultColorBlue = defaultColors[0];
          defaultAlpha = 1.0;
          }
        else if( defaultColors.size() == 3 )
          {
          defaultColorRed = defaultColors[0];
          defaultColorGreen = defaultColors[1];
          defaultColorBlue = defaultColors[2];
          defaultAlpha = 1.0;
          }
        else if( defaultColors.size() == 4 )
          {
          defaultColorRed = defaultColors[0];
          defaultColorGreen = defaultColors[1];
          defaultColorBlue = defaultColors[2];
          defaultAlpha = defaultColors[3];
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

  vtkPolyData *vtkMesh = triangularizer->GetOutput();
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
        rgbReader->GetOutput()->SetOrigin( zeroOrigin );
        }
      catch( ... )
        {
        std::cerr << "Error reading RGB file " << rgbFileName << std::endl;
        return EXIT_FAILURE;
        }
      functionalRgbImages.emplace_back(rgbReader->GetOutput() );

      // read mask

      std::string maskFileName = functionalOverlayOption->GetFunction( n )->GetParameter( 1 );

      typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( maskFileName.c_str() );
      try
        {
        maskReader->Update();
        maskReader->GetOutput()->SetOrigin( zeroOrigin );
        }
      catch( ... )
        {
        std::cerr << "Error reading mask file " << maskFileName << std::endl;
        return EXIT_FAILURE;
        }
      functionalMaskImages.emplace_back(maskReader->GetOutput() );

      if( functionalOverlayOption->GetFunction( n )->GetNumberOfParameters() > 2 )
        {
        auto alpha = parser->Convert<RealType>(
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

  vtkSmartPointer<vtkPoints> meshPoints = vtkMesh->GetPoints();
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
  colors->SetNumberOfComponents( 4 );   // R, G, B, and alpha components
  colors->SetName( "Colors" );

  for( int n = 0; n < numberOfPoints; n++ )
    {
    ImageType::IndexType index;
    ImageType::PointType imagePoint;

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      imagePoint[d] = meshPoints->GetPoint( n )[d];
      }

    RealType currentRed   = defaultColorRed / static_cast<RealType>( 255.0 );
    RealType currentGreen = defaultColorGreen / static_cast<RealType>( 255.0 );
    RealType currentBlue  = defaultColorBlue / static_cast<RealType>( 255.0 );
    RealType currentAlpha = defaultAlpha;

    for( int i = functionalAlphaValues.size() - 1; i >= 0; i-- )
      {
      bool isInsideImage = functionalMaskImages[i]->TransformPhysicalPointToIndex( imagePoint, index );

      if( isInsideImage && functionalMaskImages[i]->GetPixel( index ) != 0 )
        {
        // http://stackoverflow.com/questions/726549/algorithm-for-additive-color-mixing-for-rgb-values
        // or
        // http://en.wikipedia.org/wiki/Alpha_compositing

        RgbPixelType rgbPixel = functionalRgbImages[i]->GetPixel( index );

        RealType functionalRed = rgbPixel.GetRed() / static_cast<RealType>( 255.0 );
        RealType functionalGreen = rgbPixel.GetGreen() / static_cast<RealType>( 255.0 );
        RealType functionalBlue = rgbPixel.GetBlue() / static_cast<RealType>( 255.0 );
        RealType functionalAlpha = functionalAlphaValues[i];

        RealType backgroundRed   = currentRed;
        RealType backgroundGreen = currentGreen;
        RealType backgroundBlue  = currentBlue;
        RealType backgroundAlpha = currentAlpha;

        currentAlpha = 1.0f - ( 1.0f - functionalAlpha ) * ( 1.0f - backgroundAlpha );
        currentRed   = functionalRed   * functionalAlpha / currentAlpha + backgroundRed   * backgroundAlpha * ( 1.0f - functionalAlpha ) / currentAlpha;
        currentGreen = functionalGreen * functionalAlpha / currentAlpha + backgroundGreen * backgroundAlpha * ( 1.0f - functionalAlpha ) / currentAlpha;
        currentBlue  = functionalBlue  * functionalAlpha / currentAlpha + backgroundBlue  * backgroundAlpha * ( 1.0f - functionalAlpha ) / currentAlpha;
        }
      }

    unsigned char currentColor[4];
    currentColor[0] = static_cast<unsigned char>( currentRed   * 255.0f );
    currentColor[1] = static_cast<unsigned char>( currentGreen * 255.0f );
    currentColor[2] = static_cast<unsigned char>( currentBlue  * 255.0f );
    currentColor[3] = static_cast<unsigned char>( currentAlpha * 255.0f );

    colors->InsertNextTypedTuple( currentColor );
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

    if( numberOfIterations > 0 )
      {
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
    }

  // Write the vtk mesh to file.

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string outputFile = outputOption->GetFunction( 0 )->GetName();
    std::string ext = itksys::SystemTools::GetFilenameExtension( outputFile );
    if( strcmp( ext.c_str(), ".stl" ) == 0 )
      {
      vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
      writer->SetInputData( vtkMesh );
      writer->SetFileName( outputFile.c_str() );
      writer->Write();
      }
    if( strcmp( ext.c_str(), ".ply" ) == 0 )
      {
      vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
      writer->SetInputData( vtkMesh );
      writer->SetFileName( outputFile.c_str() );
      writer->Write();
      }
    if( strcmp( ext.c_str(), ".vtk" ) == 0 )
      {
      vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
      writer->SetInputData( vtkMesh );
      writer->SetFileName( outputFile.c_str() );
      writer->Write();
      }
    }

  vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
  std::string scalarBarTitle( "antsSurf" );
  unsigned int scalarBarNumberOfLabels = 5;

  unsigned int scalarBarWidthInPixels = 0;
  unsigned int scalarBarHeightInPixels = 0;

  bool renderScalarBar = false;

  itk::ants::CommandLineParser::OptionType::Pointer scalarBarOption = parser->GetOption( "scalar-bar" );
  if( scalarBarOption && scalarBarOption->GetNumberOfFunctions() )
    {
    renderScalarBar = true;

    std::string lookupTableFile;
    if( scalarBarOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      lookupTableFile = scalarBarOption->GetFunction( 0 )->GetName();
      }
    else
      {
      lookupTableFile = scalarBarOption->GetFunction( 0 )->GetParameter( 0 );
      if( scalarBarOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        scalarBarTitle = scalarBarOption->GetFunction( 0 )->GetParameter( 1 );
        }
      if( scalarBarOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        scalarBarNumberOfLabels = parser->Convert<unsigned int>( scalarBarOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      if( scalarBarOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        std::vector<unsigned int> dimensions = parser->ConvertVector<unsigned int>( scalarBarOption->GetFunction( 0 )->GetParameter( 3 ) );
        scalarBarWidthInPixels = dimensions[0];
        scalarBarHeightInPixels = dimensions[1];
       }
      }

    // Read in color table

    std::ifstream fileStr( lookupTableFile.c_str() );
    if( !fileStr.is_open() )
      {
      std::cerr << " Could not open file " << lookupTableFile << '\n';
      renderScalarBar = false;
      }

    int tableSize = std::count( std::istreambuf_iterator<char>( fileStr ), std::istreambuf_iterator<char>(), '\n' );
    fileStr.clear();
    fileStr.seekg( 0, std::ios::beg );

    lookupTable->SetNumberOfTableValues( tableSize );
    lookupTable->Build();

    RealType value;
    RealType redComponent;
    RealType greenComponent;
    RealType blueComponent;
    RealType alphaComponent;

    char comma;

    RealType minValue = itk::NumericTraits<RealType>::max();
    RealType maxValue = itk::NumericTraits<RealType>::min();

    unsigned int index = 0;
    while( fileStr >> value >> comma >> redComponent >> comma >> greenComponent >> comma >> blueComponent >> comma >> alphaComponent )
      {
      lookupTable->SetTableValue( index++, redComponent / static_cast<RealType>( 255.0 ),
                                           greenComponent / static_cast<RealType>( 255.0 ),
                                           blueComponent / static_cast<RealType>( 255.0 ),
                                           alphaComponent );

      if( value < minValue )
        {
        minValue = value;
        }
      if( value > maxValue )
        {
        maxValue = value;
        }
      }
    lookupTable->SetTableRange( minValue, maxValue );

    fileStr.close();
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
      Display( vtkMesh, rotationAnglesInDegrees, backgroundColor, screenCaptureFileName,
               renderScalarBar, lookupTable, scalarBarTitle, scalarBarNumberOfLabels,
               scalarBarWidthInPixels, scalarBarHeightInPixels );
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

      Display( vtkMesh, rotationAnglesInDegrees, backgroundColor, screenCaptureFileName,
               renderScalarBar, lookupTable, scalarBarTitle, scalarBarNumberOfLabels,
               scalarBarWidthInPixels, scalarBarHeightInPixels );
      }
    }

  return EXIT_SUCCESS;
}

int antsSurfaceToImage( itk::ants::CommandLineParser *parser )
{
  itk::ants::CommandLineParser::OptionType::Pointer surfaceOption =
    parser->GetOption( "mesh" );

  vtkSmartPointer<vtkPolyData> vtkMesh;

  std::string inputFile;
  if( surfaceOption && surfaceOption->GetNumberOfFunctions() > 0 )
    {
    inputFile = surfaceOption->GetFunction( 0 )->GetName();
    std::string ext = itksys::SystemTools::GetFilenameExtension( inputFile );
    try
      {
      if( strcmp( ext.c_str(), ".stl" ) == 0 )
        {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName( inputFile.c_str() );
        reader->Update();
        vtkMesh = reader->GetOutput();
        }
      if( strcmp( ext.c_str(), ".ply" ) == 0 )
        {
        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        reader->SetFileName( inputFile.c_str() );
        reader->Update();
        vtkMesh = reader->GetOutput();
        }
      if( strcmp( ext.c_str(), ".vtk" ) == 0 )
        {
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName( inputFile.c_str() );
        reader->Update();
        vtkMesh = reader->GetOutput();
        }
      }
    catch( ... )
      {
      std::cerr << "Error.  Unable to read mesh input file." << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "No mesh file specified." << std::endl;
    return EXIT_FAILURE;
    }

  double bounds[6];
  vtkMesh->GetBounds( bounds );

  std::string outputFile;
  std::vector<double> spacing;

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    outputFile = outputOption->GetFunction( 0 )->GetName();

    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      spacing.push_back( 1.0 );
      std::cout << "Warning.  No spacing is specified---defaulting to 1.0." << std::endl;
      }
    else
      {
      spacing = parser->ConvertVector<double>(
        outputOption->GetFunction( 0 )->GetParameter( 0 ) );
      }
    }
  else
    {
    std::cerr << "Error.  No output specified." << std::endl;
    return EXIT_FAILURE;
    }

  vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

  double spacing2[3]; // desired volume spacing
  if( spacing.size() == 1 )
    {
    spacing2[0] = spacing[0];
    spacing2[1] = spacing[0];
    spacing2[2] = spacing[0];
    }
  else if( spacing.size() == 3 )
    {
    spacing2[0] = spacing[0];
    spacing2[1] = spacing[1];
    spacing2[2] = spacing[2];
    }
  else
    {
    std::cerr << "Error. Incorrect spacing specified." << std::endl;
    return EXIT_FAILURE;
    }
  whiteImage->SetSpacing( spacing2 );

  // compute dimensions
  int dim[3];
  for( unsigned int i = 0; i < 3; i++ )
    {
    dim[i] = static_cast<int>( std::ceil( ( bounds[i * 2 + 1] - bounds[i * 2] ) / spacing2[i] ) );
    }
  whiteImage->SetDimensions( dim );
  whiteImage->SetExtent( 0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1 );

  double origin[3];
  origin[0] = bounds[0] + spacing2[0] / 2;
  origin[1] = bounds[2] + spacing2[1] / 2;
  origin[2] = bounds[4] + spacing2[2] / 2;
  whiteImage->SetOrigin( origin );

  whiteImage->AllocateScalars( VTK_UNSIGNED_CHAR, 1 );

  // fill the image with foreground voxels:
  unsigned char inval = 1;
  unsigned char outval = 0;
  vtkIdType count = whiteImage->GetNumberOfPoints();
  for( vtkIdType i = 0; i < count; ++i )
    {
    whiteImage->GetPointData()->GetScalars()->SetTuple1( i, inval );
    }

  // polygonal data --> image stencil:
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  pol2stenc->SetInputData( vtkMesh );
  pol2stenc->SetOutputOrigin( origin );
  pol2stenc->SetOutputSpacing( spacing2 );
  pol2stenc->SetOutputWholeExtent( whiteImage->GetExtent() );
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
  imgstenc->SetInputData( whiteImage );
  imgstenc->SetStencilConnection( pol2stenc->GetOutputPort() );
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue( outval );
  imgstenc->Update();

  // Write the vtk mesh to image file.

  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetFileName( outputFile.c_str() );
    writer->SetInputData( imgstenc->GetOutput() );
    writer->Write();
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Main input binary image for 3-D rendering.  One can also " )
      + std::string( "set a default color value in the range [0,255].  The " )
      + std::string( "fourth default color element is the alpha value in " )
      + std::string( "the range [0,1]." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "surface-image" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "surfaceImageFilename" );
    option->SetUsageOption( 1, "[surfaceImageFilename,<defaultColor=255x255x255x1>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The user can also specify a vtk polydata file to be converted " )
      + std::string( "to a binary image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mesh" );
    option->SetShortName( 'm' );
    option->SetUsageOption( 0, "meshFilename" );
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
      std::string( "Given a binary image input, the output is a vtk polydata file (possible " )
      + std::string( "extensions include .stl, .ply, and .vtk). " )
      + std::string( "Alternatively, if a mesh file is specified as input, the output  " )
      + std::string( "is an itk binary image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "surfaceFilename" );
    option->SetUsageOption( 1, "imageFilename[spacing]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
    std::string( "Add a scalar bar to the rendering for the final overlay.  One can tailor " )
    + std::string( "the aesthetic by changing the number of labels and/or the orientation and ")
    + std::string( "size of the scalar bar.  If the \'width\' > \'height\' (in pixels) then the ")
    + std::string( "orientation is horizontal.  Otherwise it is vertical (default)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "scalar-bar" );
    option->SetShortName( 'b' );
    option->SetUsageOption( 0, "lookupTable" );
    option->SetUsageOption( 1, "[lookupTable,<title=antsSurf>,<numberOfLabels=5>,<widthxheight>]" );
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
int antsSurf( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
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
    std::string( "Produce a 3-D surface rendering with optional RGB overlay.  Alternatively, " )
     + std::string( "one can input a mesh which can then be converted to a binary image. " );

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
    parser->GetOption( "surface-image" );

  itk::ants::CommandLineParser::OptionType::Pointer surfaceOption =
    parser->GetOption( "mesh" );

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
      antsImageToSurface( parser );
      }
    else
      {
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( surfaceOption && surfaceOption->GetNumberOfFunctions() > 0 )
    {
    antsSurfaceToImage( parser );
    }
  else
    {
    std::cerr << "Input not specified.  See help menu." << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
} // namespace ants
