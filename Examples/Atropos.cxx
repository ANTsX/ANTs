#include "itkAddImageFilter.h"
#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkDivideImageFilter.h"
#include "itkImage.h"
#include "itkMaskImageFilter.h"
#include "itkAtroposSegmentationImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include <string>
#include <algorithm>
#include <vector>

template <class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate()
  {
  };
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter =
      dynamic_cast<const TFilter *>( object );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    std::cout << "Iteration " << filter->GetElapsedIterations()
              << " (of " << filter->GetMaximumNumberOfIterations() << "): ";
    std::cout << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
  }
};

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
}

template <unsigned int ImageDimension>
int AtroposSegmentation( itk::CommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;

  typedef unsigned char                         LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef  itk::AtroposSegmentationImageFilter
    <InputImageType, LabelImageType> SegmentationFilterType;
  typename SegmentationFilterType::Pointer segmenter
    = SegmentationFilterType::New();

  typedef CommandIterationUpdate<SegmentationFilterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  segmenter->AddObserver( itk::IterationEvent(), observer );
  // segmenter->DebugOn();

  /**
   * Initialization
   */
  typename itk::CommandLineParser::OptionType::Pointer initializationOption =
    parser->GetOption( "initialization" );
  if( initializationOption
      && initializationOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << "Incorrect initialization option specification." << std::endl;
    std::cerr << "   " << initializationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(
      ( initializationOption->GetParameter( 0 ) ).c_str() );
    reader->Update();

    segmenter->SetInput( reader->GetOutput() );
    segmenter->SetNumberOfClasses( parser->Convert<unsigned int>(
                                     initializationOption->GetParameter( 1 ) ) );

    std::string initializationStrategy = initializationOption->GetValue();
    ConvertToLowerCase( initializationStrategy );
    if( !initializationStrategy.compare( std::string( "random" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::Random );
      }
    else if( !initializationStrategy.compare( std::string( "otsu" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::Otsu );
      }
    else if( !initializationStrategy.compare( std::string( "kmeans" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::KMeans );
      }
    else if( !initializationStrategy.compare(
               std::string( "priorprobabilityimages" ) ) )
      {
      segmenter->SetInitializationStrategy(
        SegmentationFilterType::PriorProbabilityImages );
      if( initializationOption->GetNumberOfParameters() < 4 )
        {
        std::cerr << "Incorrect initialization option specification."
                  << std::endl;
        std::cerr << "   " << initializationOption->GetDescription()
                  << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetParameter( 3 ) ) );
      segmenter->SetAdaptiveSmoothingWeight( 0, parser->Convert<float>(
                                               initializationOption->GetParameter( 4 ) ) );

      std::string filename = initializationOption->GetParameter( 2 );

      if( filename.find( std::string( "%" ) ) != std::string::npos )
        {
        itk::NumericSeriesFileNames::Pointer fileNamesCreator =
          itk::NumericSeriesFileNames::New();
        fileNamesCreator->SetStartIndex( 1 );
        fileNamesCreator->SetEndIndex( segmenter->GetNumberOfClasses() );
        fileNamesCreator->SetSeriesFormat( filename.c_str() );
        const std::vector<std::string> & imageNames
          = fileNamesCreator->GetFileNames();
        for( unsigned int k = 0; k < imageNames.size(); k++ )
          {
          typedef itk::ImageFileReader<InputImageType> ReaderType;
          typename ReaderType::Pointer reader = ReaderType::New();
          reader->SetFileName( imageNames[k].c_str() );
          reader->Update();

          segmenter->SetPriorProbabilityImage( k + 1, reader->GetOutput() );
          }
        }
      else
        {
        typedef itk::VectorImage<PixelType, ImageDimension> VectorImageType;
        typedef itk::ImageFileReader<VectorImageType>       ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( filename.c_str() );
        reader->Update();

        if(  reader->GetOutput()->GetNumberOfComponentsPerPixel()
             != segmenter->GetNumberOfClasses() )
          {
          std::cerr << "The number of components does not match the number of "
                    << "classes." << std::endl;
          return EXIT_FAILURE;
          }

        typedef itk::VectorIndexSelectionCastImageFilter
          <VectorImageType, InputImageType> CasterType;
        typename CasterType::Pointer caster = CasterType::New();
        caster->SetInput( reader->GetOutput() );
        for( unsigned int k = 0; k < segmenter->GetNumberOfClasses(); k++ )
          {
          caster->SetIndex( k );
          caster->Update();
          segmenter->SetPriorProbabilityImage( k + 1, caster->GetOutput() );
          }
        }
      }
    else if( !initializationStrategy.compare( std::string( "priorlabelimage" ) ) )
      {
      segmenter->SetInitializationStrategy(
        SegmentationFilterType::PriorLabelImage );

      if( initializationOption->GetNumberOfParameters() < 4 )
        {
        std::cerr << "Incorrect initialization option specification." << std::endl;
        std::cerr << "   " << initializationOption->GetDescription() << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetParameter( 3 ) ) );
      segmenter->SetAdaptiveSmoothingWeight( 0, parser->Convert<float>(
                                               initializationOption->GetParameter( 4 ) ) );

      std::string filename = initializationOption->GetParameter( 2 );
      typedef itk::ImageFileReader<LabelImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filename.c_str() );
      reader->Update();

      segmenter->SetPriorLabelImage( reader->GetOutput() );
      }
    else
      {
      std::cerr << "Unrecognized initialization strategy request." << std::endl;
      return EXIT_FAILURE;
      }
    }

  /**
   * number of iterations
   */
  typename itk::CommandLineParser::OptionType::Pointer iterationsOption =
    parser->GetOption( "number-of-iterations" );
  if( iterationsOption )
    {
    segmenter->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
                                               iterationsOption->GetValue() ) );
    }

  /**
   * convergence threshold
   */
  typename itk::CommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence-threshold" );
  if( convergenceOption )
    {
    segmenter->SetConvergenceThreshold( parser->Convert<float>(
                                          convergenceOption->GetValue() ) );
    }

  /**
   * Mask image
   */
  typename itk::CommandLineParser::OptionType::Pointer maskOption =
    parser->GetOption( "mask-image" );
  if( maskOption && maskOption->GetNumberOfValues() )
    {
    try
      {
      typedef  itk::ImageFileReader<LabelImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( ( maskOption->GetValue() ).c_str() );
      reader->Update();

      segmenter->SetMaskImage( reader->GetOutput() );
      }
    catch( ... )
      {
      }
    }

  /**
   * BSpline options
   */

  typename itk::CommandLineParser::OptionType::Pointer bsplineOption =
    parser->GetOption( "bspline" );
  if( bsplineOption && bsplineOption->GetNumberOfValues() )
    {
    if( bsplineOption->GetNumberOfParameters() > 0 )
      {
      std::vector<unsigned int> numLevels = parser->ConvertVector<unsigned int>(
          bsplineOption->GetParameter( 0 ) );
      typename SegmentationFilterType::ArrayType numberOfFittingLevels;

      if( numLevels.size() == 1 )
        {
        numberOfFittingLevels.Fill( numLevels[0] );
        }
      else if( numLevels.size() == ImageDimension )
        {
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          numberOfFittingLevels[d] = numLevels[d];
          }
        }
      else
        {
        std::cerr << "Incorrect number of levels" << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetNumberOfLevels( numberOfFittingLevels );
      }
    if( bsplineOption->GetNumberOfParameters() > 2 )
      {
      segmenter->SetSplineOrder( parser->Convert<unsigned int>(
                                   bsplineOption->GetParameter( 2 ) ) );
      }
    if( bsplineOption->GetNumberOfParameters() > 1 )
      {
      std::vector<unsigned int> array = parser->ConvertVector<unsigned int>(
          bsplineOption->GetParameter( 1 ) );
      typename SegmentationFilterType::ArrayType numberOfControlPoints;
      if( array.size() == 1 )
        {
        numberOfControlPoints.Fill( array[0] + segmenter->GetSplineOrder() );
        }
      else if( array.size() == ImageDimension )
        {
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          numberOfControlPoints[d] = array[d] + segmenter->GetSplineOrder();
          }
        }
      else
        {
        std::cerr << "Incorrect mesh resolution" << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetNumberOfControlPoints( numberOfControlPoints );
      }
    }

  /**
   * labels
   */
  typename itk::CommandLineParser::OptionType::Pointer labelOption =
    parser->GetOption( "labels" );
  if( labelOption && labelOption->GetNumberOfValues() > 0 )
    {
    typename SegmentationFilterType::LabelParameterMapType labelMap;
    for( unsigned int n = 0; n < labelOption->GetNumberOfValues(); n++ )
      {
      typename SegmentationFilterType::LabelParametersType labelPair;

      float labelSigma = parser->Convert<float>(
          labelOption->GetParameter( n, 0 ) );
      float labelBoundaryProbability = 0.75;
      if( labelOption->GetNumberOfParameters( n ) > 1 )
        {
        labelBoundaryProbability = parser->Convert<float>(
            labelOption->GetParameter( n, 1 ) );
        if( labelBoundaryProbability < 0.0 )
          {
          labelBoundaryProbability = 0.0;
          }
        if( labelBoundaryProbability > 1.0 )
          {
          labelBoundaryProbability = 1.0;
          }
        }
      labelPair.first = labelSigma;
      labelPair.second = labelBoundaryProbability;

      unsigned int whichClass = parser->Convert<unsigned int>(
          labelOption->GetValue( n ) );

      labelMap[whichClass] = labelPair;
      }
    segmenter->SetPriorLabelParameterMap( labelMap );
    }

  /**
   * auxiliary images
   */
  typename itk::CommandLineParser::OptionType::Pointer auxOption =
    parser->GetOption( "auxiliary-image" );
  if( auxOption && auxOption->GetNumberOfValues() > 0 )
    {
    for( unsigned int n = 0; n < auxOption->GetNumberOfValues(); n++ )
      {
      typedef itk::ImageFileReader<InputImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( auxOption->GetParameter( n, 0 ) );
      reader->Update();

      segmenter->SetAuxiliaryImage( n + 1, reader->GetOutput() );
      if( initializationOption->GetNumberOfParameters( n ) > 1 )
        {
        segmenter->SetAdaptiveSmoothingWeight( n + 1, parser->Convert<float>(
                                                 initializationOption->GetParameter( n, 1 ) ) );
        }
      }
    }

  /**
   * MRF options
   */
  typename itk::CommandLineParser::OptionType::Pointer mrfOption =
    parser->GetOption( "mrf" );
  if( mrfOption && mrfOption->GetNumberOfValues() )
    {
    if( mrfOption->GetNumberOfParameters() > 0 )
      {
      segmenter->SetMRFSmoothingFactor( parser->Convert<float>(
                                          mrfOption->GetParameter( 0 ) ) );
      }
    if( mrfOption->GetNumberOfParameters() > 1 )
      {
      std::vector<unsigned int> array = parser->ConvertVector<unsigned int>(
          mrfOption->GetParameter( 1 ) );
      typename SegmentationFilterType::ArrayType radius;
      if( array.size() == 1 )
        {
        radius.Fill( array[0] );
        }
      else if( array.size() == ImageDimension )
        {
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          radius[d] = array[d];
          }
        }
      else
        {
        std::cerr << "MRF radius size needs to be equal to the image dimension."
                  << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetMRFRadius( radius );
      }
    if( mrfOption->GetNumberOfParameters() > 2 )
      {
      segmenter->SetMRFSigmoidAlpha( parser->Convert<float>(
                                       mrfOption->GetParameter( 2 ) ) );
      }
    if( mrfOption->GetNumberOfParameters() > 3 )
      {
      segmenter->SetMRFSigmoidBeta( parser->Convert<float>(
                                      mrfOption->GetParameter( 3 ) ) );
      }
    }

  /**
   * euclidean distance
   */
  typename itk::CommandLineParser::OptionType::Pointer distanceOption =
    parser->GetOption( "use-euclidean-distance" );
  if( distanceOption )
    {
    segmenter->SetUseEuclideanDistanceForPriorLabels(
      parser->Convert<bool>( distanceOption->GetValue() ) );
    }

  /**
   * memory-usage
   */
  typename itk::CommandLineParser::OptionType::Pointer memoryOption =
    parser->GetOption( "minimize-memory-usage" );
  if( memoryOption )
    {
    segmenter->SetMinimizeMemoryUsage( parser->Convert<bool>(
                                         memoryOption->GetValue() ) );
    }

  try
    {
    segmenter->Update();
    }
  catch( itk::ExceptionObject exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * output
   */
  typename itk::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption )
    {
    if( outputOption->GetNumberOfParameters() == 0 )
      {
      typedef  itk::ImageFileWriter<LabelImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( segmenter->GetOutput() );
      writer->SetFileName( ( outputOption->GetValue() ).c_str() );
      writer->Update();
      }
    if( outputOption->GetNumberOfParameters() > 0 )
      {
      typedef  itk::ImageFileWriter<LabelImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( segmenter->GetOutput() );
      writer->SetFileName( ( outputOption->GetParameter( 0 ) ).c_str() );
      writer->Update();
      }
    if( outputOption->GetNumberOfParameters() > 1 )
      {
      std::string filename = outputOption->GetParameter( 1 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator =
        itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames
        = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < imageNames.size(); i++ )
        {
        std::cout << "Writing posterior image (class " << i + 1 << ")"
                  << std::endl;
        typename InputImageType::Pointer probabilityImage
          = segmenter->GetPosteriorProbabilityImage( i + 1 );

        if( segmenter->GetMaskImage() )
          {
          typedef itk::MaskImageFilter<InputImageType,
                                       LabelImageType, InputImageType> MaskerType;
          typename MaskerType::Pointer masker = MaskerType::New();
          masker->SetInput1( probabilityImage );
          masker->SetInput2( segmenter->GetMaskImage() );
          masker->SetOutsideValue( 0 );
          masker->Update();

          probabilityImage = masker->GetOutput();
          }

        typedef  itk::ImageFileWriter<InputImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( probabilityImage );
        writer->SetFileName( imageNames[i].c_str() );
        writer->Update();
        }
      }
    if( outputOption->GetNumberOfParameters() > 2 )
      {
      std::string filename = outputOption->GetParameter( 2 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator =
        itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames
        = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < segmenter->GetNumberOfClasses(); i++ )
        {
        if( segmenter->GetPriorProbabilityImage( i + 1 ) ||
            segmenter->GetPriorLabelImage() )
          {
          std::cout << "Writing B-spline image (class " << i + 1 << ")"
                    << std::endl;

          typename InputImageType::Pointer bsplineImage = segmenter->
            CalculateSmoothIntensityImageFromPriorProbabilityImage( 0, i + 1 );

          typedef  itk::ImageFileWriter<InputImageType> WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( bsplineImage );
          writer->SetFileName( imageNames[i].c_str() );
          writer->Update();
          }
        }
      }
    if( outputOption->GetNumberOfParameters() > 3 )
      {
      std::string filename = outputOption->GetParameter( 3 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator =
        itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames
        = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < segmenter->GetNumberOfClasses(); i++ )
        {
        if( segmenter->GetPriorProbabilityImage( i + 1 ) ||
            segmenter->GetPriorLabelImage() )
          {
          std::cout << "Writing distance image (class " << i + 1 << ")"
                    << std::endl;

          typename InputImageType::Pointer distanceImage = segmenter->
            GetDistancePriorProbabilityImageFromPriorLabelImage( i + 1 );

          typedef  itk::ImageFileWriter<InputImageType> WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( distanceImage );
          writer->SetFileName( imageNames[i].c_str() );
          writer->Update();
          }
        }
      }
    }

  segmenter->Print( std::cout, 5 );

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::CommandLineParser *parser )
{
  typedef itk::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Option 1:  Random[inputImage,numberOfClasses]\n" )
      + std::string( "\t  Option 2:  Kmeans[inputImage,numberOfClasses]\n" )
      + std::string( "\t  Option 3:  Otsu[inputImage,numberOfClasses]\n" )
      + std::string( "\t  Option 4:  PriorProbabilityImages[inputImage,numberOfClasses," )
      + std::string(
        "fileSeriesFormat(index=1 to numberOfClasses-1) or vectorImage,priorWeighting,<adaptiveSmoothingWeight>]\n" )
      + std::string(
        "\t  Option 5:  PriorLabelImage[inputImage,numberOfClasses,labelImage,priorWeighting,<adaptiveSmoothingWeight>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initialization" );
    option->SetShortName( 'i' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "maskImage" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "[image,<adaptiveSmoothingWeight>]" )
      + std::string( " -- adaptive smoothing only applies to initialization with prior image(s)" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "auxiliary-image" );
    option->SetShortName( 'a' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Number of iterations" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "number-of-iterations" );
    option->SetShortName( 'n' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Convergence threshold" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "convergence-threshold" );
    option->SetShortName( 'c' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "[<smoothingFactor>,<radius>,<sigmoidAlpha>,<sigmoidBeta>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mrf" );
    option->SetShortName( 'm' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "[classifiedImage," )
        + std::string( "<posteriorProbabilityImageFileNameFormat>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "minimize-memory-usage=true/false" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "minimize-memory-usage" );
    option->SetShortName( 'u' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "[<numberOfLevels>,<initialMeshResolution>,<splineOrder>]" )
      + std::string( " -- only applies to initialization with prior image(s)" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bspline" );
    option->SetShortName( 'b' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "use-euclidean-distance-for-prior-label-image=true/false" )
      + std::string( " -- only applies to initialization with a prior label image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "use-euclidean-distance" );
    option->SetShortName( 'e' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "whichLabel[sigma,<boundaryProbability>]" )
      + std::string( " -- only applies to initialization with a prior label image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "labels" );
    option->SetShortName( 'l' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cout << "Usage: " << argv[0]
              << " imageDimension args" << std::endl;
    exit( 1 );
    }

  itk::CommandLineParser::Pointer parser = itk::CommandLineParser::New();
  parser->SetCommand( argv[0] );

  parser->SetCommandDescription( "Atropos:  A priori classification with registration initialized template assistance." );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 3 || parser->Convert<bool>(
        parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5 );
    exit( EXIT_FAILURE );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      AtroposSegmentation<2>( parser );
      break;
    case 3:
      AtroposSegmentation<3>( parser );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
