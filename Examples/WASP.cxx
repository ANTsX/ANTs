#include "itkAddImageFilter.h"
#include "itkWASPCommandLineOption.h"
#include "itkWASPCommandLineParser.h"
#include "itkDivideImageFilter.h"
#include "itkImage.h"
#include "itkMaskImageFilter.h"
#include "itkWASPSegmentationImageFilter.h"
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
int WASPSegmentation( itk::WASPCommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;

  typedef unsigned char                         LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef  itk::WASPSegmentationImageFilter
  <InputImageType, LabelImageType> SegmentationFilterType;
  typename SegmentationFilterType::Pointer segmenter
    = SegmentationFilterType::New();

  typedef CommandIterationUpdate<SegmentationFilterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  segmenter->AddObserver( itk::IterationEvent(), observer );

  /**
   * Initialization
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer initializationOption =
    parser->GetOption( "initialization" );
  if( initializationOption
      && initializationOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << "Incorrect point set option specification." << std::endl;
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
    if( !initializationStrategy.compare( std::string( "otsu" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::Otsu );
      }
    else if( !initializationStrategy.compare( std::string( "kmeans" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::KMeans );
      }
    else if( !initializationStrategy.compare( std::string( "priorprobabilityimages" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::PriorProbabilityImages );

      if( initializationOption->GetNumberOfParameters() > 3 )
        {
        segmenter->SetPriorProbabilityWeighting( parser->Convert<float>(
                                                   initializationOption->GetParameter( 3 ) ) );
        }
      if( initializationOption->GetNumberOfParameters() > 2 )
        {
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
            std::cerr << "The number of components does not match the number of classes." << std::endl;
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
      }
    else if( !initializationStrategy.compare( std::string( "priorlabelimage" ) ) )
      {
      segmenter->SetInitializationStrategy( SegmentationFilterType::PriorLabelImage );

      if( initializationOption->GetNumberOfParameters() > 3 )
        {
        segmenter->SetPriorProbabilityWeighting( parser->Convert<float>(
                                                   initializationOption->GetParameter( 3 ) ) );
        }
      if( initializationOption->GetNumberOfParameters() > 2 )
        {
        std::string filename = initializationOption->GetParameter( 2 );

        typedef itk::ImageFileReader<LabelImageType> ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( filename.c_str() );
        reader->Update();

        segmenter->SetPriorLabelImage( reader->GetOutput() );
        }
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
  typename itk::WASPCommandLineParser::OptionType::Pointer iterationsOption =
    parser->GetOption( "number-of-iterations" );
  if( iterationsOption )
    {
    segmenter->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
                                               iterationsOption->GetValue() ) );
    }

  /**
   * convergence threshold
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence-threshold" );
  if( convergenceOption )
    {
    segmenter->SetConvergenceThreshold( parser->Convert<float>(
                                          convergenceOption->GetValue() ) );
    }
  std::cout << segmenter->GetConvergenceThreshold() << std::endl;

  /**
   * Mask image
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer maskOption =
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

  typename itk::WASPCommandLineParser::OptionType::Pointer bsplineOption =
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
   * label-sigma
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer labelSigmasOption =
    parser->GetOption( "label-sigma" );
  if( labelSigmasOption && labelSigmasOption->GetNumberOfValues() )
    {
    std::vector<float> labelSigmas;
    for( unsigned int n = 0; n < segmenter->GetNumberOfClasses(); n++ )
      {
      labelSigmas.push_back( 0.0 );
      }
    for( unsigned int n = 0; n < labelSigmasOption->GetNumberOfValues(); n++ )
      {
      unsigned int whichClass = parser->Convert<unsigned int>(
          labelSigmasOption->GetValue( n ) );
      if( whichClass <= labelSigmas.size() )
        {
        labelSigmas[whichClass - 1] = parser->Convert<float>(
            labelSigmasOption->GetParameter( n, 0 ) );
        }
      }
    segmenter->SetPriorLabelSigmas( labelSigmas );
    }

  /**
   * MRF options
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer mrfOption =
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

  segmenter->Update();

  /**
   * output
   */
  typename itk::WASPCommandLineParser::OptionType::Pointer outputOption =
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
        std::cout << "Writing posterior image (class " << i + 1 << ")" << std::endl;
        typename InputImageType::Pointer probabilityImage
          = segmenter->GetPosteriorImage( i ); // CalculatePosteriorProbabilityImage( i + 1 );

        if( segmenter->GetMaskImage() )
          {
          typedef itk::MaskImageFilter<InputImageType, LabelImageType, InputImageType> MaskerType;
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
        if( segmenter->GetPriorProbabilityImage( i + 1 ) || segmenter->GetPriorLabelImage() )
          {
          std::cout << "Writing B-spline image (class " << i + 1 << ")" << std::endl;

          typename InputImageType::Pointer bsplineImage =
            segmenter->CalculateSmoothIntensityImageFromPriorProbabilityImage( i + 1 );

          typedef  itk::ImageFileWriter<InputImageType> WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( bsplineImage );
          writer->SetFileName( imageNames[i].c_str() );
          writer->Update();
          }
        }
      }
    }

  segmenter->Print( std::cout, 5 );

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::WASPCommandLineParser *parser )
{
  typedef itk::WASPCommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "\n\tOption 1:  Kmeans[inputImage,numberOfClasses]\n" )
      + std::string( "\tOption 2:  Otsu[inputImage,numberOfClasses]\n" )
      + std::string( "\tOption 3:  PriorProbabilities[inputImage,numberOfClasses," )
      + std::string( "fileSeriesFormat(index=1 to numberOfClasses) or vectorImage,priorWeighting]\n" )
      + std::string( "\tOption 4:  PriorLabelImage[inputImage,numberOfClasses,labelImage,priorWeighting]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initialization" );
    option->SetShortName( 'i' );
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
    std::string description = std::string( "maskImage" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "[<smoothingFactor>,<radius>,<sigmoidAlpha>,<sigmoidBeta>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mrf" );
    option->SetShortName( 'm' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "[<numberOfLevels>,<initialMeshResolution>,<splineOrder>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bspline" );
    option->SetShortName( 'b' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "whichLabel[sigma]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "label-sigma" );
    option->SetShortName( 'l' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "[classifiedImage," )
        + std::string( "<posteriorProbabilityImageFileNameFormat>," )
        + std::string( "<bSplineImageFileNameFormat>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description
      = std::string( "Print menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetShortName( 'h' );
    option->SetDescription( description );
    parser->AddOption( option );
    }
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cout
    << " WASP : Wieldy Atlas-based Segmentation Protocol -- stinging the segmentation problem in the behind. "
    << std::endl;
    std::cout << " wieldy [weel-dee] " << std::endl;
    std::cout <<  " –adjective: readily wielded or managed, as in use or action. " << std::endl;
    std::cout << "Usage: " << argv[0] << " imageDimension args" << std::endl;
    std::cout << "call: " << argv[0] << " -h " << std::endl;
    std::cout << " this program has not been evaluated " << std::endl;
    exit( 1 );
    }

  itk::WASPCommandLineParser::Pointer parser = itk::WASPCommandLineParser::New();
  parser->SetCommand( argv[0] );

  parser->SetCommandDescription( "WASP segmentation." );
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
      WASPSegmentation<2>( parser );
      break;
    case 3:
      WASPSegmentation<3>( parser );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
