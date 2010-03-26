#include "itkAtroposSegmentationImageFilter.h"
#include "itkBoxPlotQuantileListSampleFilter.h"
#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkGaussianListSampleFunction.h"
#include "itkGrubbsRosnerListSampleFilter.h"
#include "itkHistogramParzenWindowsListSampleFunction.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkListSampleToListSampleFilter.h"
#include "itkManifoldParzenWindowsListSampleFunction.h"
#include "itkMaskImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkPassThroughListSampleFilter.h"
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
      && initializationOption->GetNumberOfParameters() < 1 )
    {
    std::cerr << "Incorrect initialization option specification." << std::endl;
    std::cerr << "   " << initializationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    segmenter->SetNumberOfClasses( parser->Convert<unsigned int>(
                                     initializationOption->GetParameter( 0 ) ) );

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
      if( initializationOption->GetNumberOfParameters() < 3 )
        {
        std::cerr << "Incorrect initialization option specification."
                  << std::endl;
        std::cerr << "   " << initializationOption->GetDescription()
                  << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetParameter( 2 ) ) );

      std::string filename = initializationOption->GetParameter( 1 );

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

      if( initializationOption->GetNumberOfParameters() < 3 )
        {
        std::cerr << "Incorrect initialization option specification." << std::endl;
        std::cerr << "   " << initializationOption->GetDescription() << std::endl;
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetParameter( 2 ) ) );

      std::string filename = initializationOption->GetParameter( 1 );
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
   * convergence options
   */
  typename itk::CommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence" );
  if( convergenceOption )
    {
    if( convergenceOption->GetNumberOfParameters() > 0 )
      {
      segmenter->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
                                                 convergenceOption->GetParameter( 0 ) ) );
      }
    if( convergenceOption->GetNumberOfParameters() > 1 )
      {
      segmenter->SetConvergenceThreshold( parser->Convert<float>(
                                            convergenceOption->GetParameter( 1 ) ) );
      }
    }

  /**
   * Mask image
   */
  typename itk::CommandLineParser::OptionType::Pointer maskOption =
    parser->GetOption( "mask-image" );
  if( maskOption && maskOption->GetNumberOfValues() > 0 )
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
   * intensity images
   */
  typename itk::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "intensity-image" );
  if( imageOption && imageOption->GetNumberOfValues() > 0 )
    {
    for( unsigned int n = 0; n < imageOption->GetNumberOfValues(); n++ )
      {
      typedef itk::ImageFileReader<InputImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      if( imageOption->GetNumberOfParameters( n ) > 0 )
        {
        reader->SetFileName( imageOption->GetParameter( n, 0 ) );
        }
      else
        {
        reader->SetFileName( imageOption->GetValue( n ) );
        }
      reader->Update();

      if( n == imageOption->GetNumberOfValues() - 1 )
        {
        segmenter->SetInput( reader->GetOutput() );
        if( imageOption->GetNumberOfParameters( 0 ) > 1 )
          {
          segmenter->SetAdaptiveSmoothingWeight( 0, parser->Convert<float>(
                                                   imageOption->GetParameter( 0, 1 ) ) );
          }
        else
          {
          segmenter->SetAdaptiveSmoothingWeight( 0, 0.0 );
          }
        }
      else
        {
        segmenter->SetIntensityImage( n + 1, reader->GetOutput() );
        if( imageOption->GetNumberOfParameters( n + 1 ) > 1 )
          {
          segmenter->SetAdaptiveSmoothingWeight( n + 1, parser->Convert<float>(
                                                   imageOption->GetParameter( n + 1, 1 ) ) );
          }
        else
          {
          segmenter->SetAdaptiveSmoothingWeight( n + 1, 0.0 );
          }
        }
      }
    }
  else
    {
    std::cerr << "No input images were specified.  Specify an input image"
              << " with the -a option" << std::endl;
    }

  /**
   * MRF options
   */
  typename itk::CommandLineParser::OptionType::Pointer mrfOption =
    parser->GetOption( "mrf" );
  if( mrfOption && mrfOption->GetNumberOfValues() > 0 )
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
    }

  /**
   * euclidean distance
   */
  typename itk::CommandLineParser::OptionType::Pointer distanceOption =
    parser->GetOption( "use-euclidean-distance" );
  if( distanceOption && distanceOption->GetNumberOfValues() > 0 )
    {
    segmenter->SetUseEuclideanDistanceForPriorLabels(
      parser->Convert<bool>( distanceOption->GetValue() ) );
    }

  /**
   * memory-usage
   */
  typename itk::CommandLineParser::OptionType::Pointer memoryOption =
    parser->GetOption( "minimize-memory-usage" );
  if( memoryOption && memoryOption->GetNumberOfValues() > 0 )
    {
    segmenter->SetMinimizeMemoryUsage( parser->Convert<bool>(
                                         memoryOption->GetValue() ) );
    }

  /**
   * likelihood
   */
  typename itk::CommandLineParser::OptionType::Pointer likelihoodOption =
    parser->GetOption( "likelihood-model" );
  if( likelihoodOption && likelihoodOption->GetNumberOfValues() > 0 )
    {
    std::string likelihoodModel = likelihoodOption->GetValue();
    ConvertToLowerCase( likelihoodModel );
    if( !likelihoodModel.compare( std::string( "gaussian" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::Statistics::GaussianListSampleFunction
        <SampleType, float, float> LikelihoodType;
      for( unsigned int n = 0; n < segmenter->GetNumberOfClasses(); n++ )
        {
        typename LikelihoodType::Pointer gaussianLikelihood =
          LikelihoodType::New();
        segmenter->SetLikelihoodFunction( n, gaussianLikelihood );
        }
      }
    else if( !likelihoodModel.compare( std::string( "manifoldparzenwindows" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::Statistics::ManifoldParzenWindowsListSampleFunction
        <SampleType, float, float> LikelihoodType;

      float regularizationSigma = 1.0;
      if( likelihoodOption->GetNumberOfParameters() > 0 )
        {
        regularizationSigma = parser->Convert<float>(
            likelihoodOption->GetParameter( 0 ) );
        }
      unsigned int evalNeighborhood = 50;
      if( likelihoodOption->GetNumberOfParameters() > 1 )
        {
        evalNeighborhood = parser->Convert<unsigned int>(
            likelihoodOption->GetParameter( 1 ) );
        }
      unsigned int covNeighborhood = 0;
      if( likelihoodOption->GetNumberOfParameters() > 2 )
        {
        covNeighborhood = parser->Convert<unsigned int>(
            likelihoodOption->GetParameter( 2 ) );
        }
      float covSigma = 1.0;
      if( likelihoodOption->GetNumberOfParameters() > 3 )
        {
        covSigma = parser->Convert<float>(
            likelihoodOption->GetParameter( 3 ) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfClasses(); n++ )
        {
        typename LikelihoodType::Pointer mpwLikelihood =
          LikelihoodType::New();
        mpwLikelihood->SetRegularizationSigma( regularizationSigma );
        mpwLikelihood->SetEvaluationKNeighborhood( evalNeighborhood );
        mpwLikelihood->SetCovarianceKNeighborhood( covNeighborhood );
        mpwLikelihood->SetKernelSigma( covSigma );
        segmenter->SetLikelihoodFunction( n, mpwLikelihood );
        }
      }
    else if( !likelihoodModel.compare( std::string( "histogramparzenwindows" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::Statistics::HistogramParzenWindowsListSampleFunction
        <SampleType, float, float> LikelihoodType;

      float sigma = 1.0;
      if( likelihoodOption->GetNumberOfParameters() > 0 )
        {
        sigma = parser->Convert<float>(
            likelihoodOption->GetParameter( 0 ) );
        }
      unsigned int numberOfBins = 32;
      if( likelihoodOption->GetNumberOfParameters() > 1 )
        {
        numberOfBins = parser->Convert<unsigned int>(
            likelihoodOption->GetParameter( 1 ) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfClasses(); n++ )
        {
        typename LikelihoodType::Pointer hpwLikelihood =
          LikelihoodType::New();
        hpwLikelihood->SetSigma( sigma );
        hpwLikelihood->SetNumberOfHistogramBins( numberOfBins );
        segmenter->SetLikelihoodFunction( n, hpwLikelihood );
        }
      }
    else
      {
      std::cerr << "Unrecognized likelihood model request." << std::endl;
      return EXIT_FAILURE;
      }
    }

  /**
   * outliers?
   */
  typename itk::CommandLineParser::OptionType::Pointer outlierOption =
    parser->GetOption( "winsorize-outliers" );
  if( outlierOption && outlierOption->GetNumberOfValues() > 0 )
    {
    std::string outlierStrategy = outlierOption->GetValue();
    ConvertToLowerCase( outlierStrategy );
    if( !outlierStrategy.compare( std::string( "boxplot" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::Statistics::BoxPlotQuantileListSampleFilter<SampleType>
        SampleFilterType;
      typename SampleFilterType::Pointer boxplotFilter =
        SampleFilterType::New();

      if( outlierOption->GetNumberOfParameters( 0 ) > 0 )
        {
        boxplotFilter->SetLowerPercentile( parser->Convert<float>(
                                             outlierOption->GetParameter( 0 ) ) );
        }
      if( outlierOption->GetNumberOfParameters( 0 ) > 1 )
        {
        boxplotFilter->SetUpperPercentile( parser->Convert<float>(
                                             outlierOption->GetParameter( 1 ) ) );
        }
      if( outlierOption->GetNumberOfParameters( 0 ) > 2 )
        {
        boxplotFilter->SetWhiskerScalingFactor( parser->Convert<float>(
                                                  outlierOption->GetParameter( 2 ) ) );
        }
      segmenter->SetOutlierHandlingFilter( boxplotFilter );
      }
    else if( !outlierStrategy.compare( std::string( "grubbsrosner" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::Statistics::GrubbsRosnerListSampleFilter<SampleType>
        SampleFilterType;
      typename SampleFilterType::Pointer grubbsFilter =
        SampleFilterType::New();

      if( outlierOption->GetNumberOfParameters( 0 ) > 0 )
        {
        grubbsFilter->SetSignificanceLevel( parser->Convert<float>(
                                              outlierOption->GetParameter( 0 ) ) );
        }
      if( outlierOption->GetNumberOfParameters( 0 ) > 1 )
        {
        grubbsFilter->SetWinsorizingLevel( parser->Convert<float>(
                                             outlierOption->GetParameter( 1 ) ) );
        }
      segmenter->SetOutlierHandlingFilter( grubbsFilter );
      }
    else
      {
      std::cerr << "Unrecognized outlier handling strategy request." << std::endl;
      return EXIT_FAILURE;
      }
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
  if( outputOption && outputOption->GetNumberOfValues() > 0 )
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
      std::string( "Option 1:  Random[numberOfClasses]\n" )
      + std::string( "\t  Option 2:  Kmeans[numberOfClasses]\n" )
      + std::string( "\t  Option 3:  Otsu[numberOfClasses]\n" )
      + std::string( "\t  Option 4:  PriorProbabilityImages[numberOfClasses," )
      + std::string( "fileSeriesFormat(index=1 to numberOfClasses) or vectorImage,priorWeighting]\n" )
      + std::string( "\t  Option 5:  PriorLabelImage[numberOfClasses,labelImage,priorWeighting]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initialization" );
    option->SetShortName( 'i' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "[intensityImage,<adaptiveSmoothingWeight>]" )
      + std::string( " -- adaptive smoothing only applies to initialization with prior image(s)" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "intensity-image" );
    option->SetShortName( 'a' );
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
    std::string description =
      std::string( "[<numberOfIterations>,<convergenceThreshold>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "convergence" );
    option->SetShortName( 'c' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Option 1: Gaussian\n" )
      + std::string(
        "\t  Option 2: ManifoldParzenWindows[<pointSetSigma=1.0>,<evaluationKNeighborhood=50>,<CovarianceKNeighborhood=0>,<kernelSigma=0>] \n" )
      + std::string( "\t  Option 3: HistogramParzenWindows[<Sigma=1.0>,<numberOfBins=32>]" );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "likelihood-model" );
    option->SetShortName( 'k' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "[<smoothingFactor>,<radius>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mrf" );
    option->SetShortName( 'm' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

//   {
//   std::string description = std::string( "[image,<adaptiveSmoothingWeight>]" ) +
//     std::string( " -- adaptive smoothing only applies to initialization with prior image(s)" );
//
//   OptionType::Pointer option = OptionType::New();
//   option->SetLongName( "vector-image" );
//   option->SetShortName( 'v' );
//   option->SetDescription( description );
//   parser->AddOption( option );
//   }

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
    std::string description = std::string( "minimize-memory-usage=1/(0)" );

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
      + std::string( " -- only applies to initialization with prior image(s) and non-zero adaptive smoothing parameter." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bspline" );
    option->SetShortName( 'b' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "use euclidean or geodesic distance for prior label distance maps = 1/(0)" )
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
      std::string( "Option 1: BoxPlot[<lowerPercentile=0.25>,<upperPercentile=0.75>,<whiskerLength=1.5>]\n" )
      + std::string( "\t  Option 2: GrubbsRosner[<significanceLevel=0.05>,<winsorizingLevel=0.10>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "winsorize-outliers" );
    option->SetShortName( 'w' );
    option->SetDescription( description );
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
