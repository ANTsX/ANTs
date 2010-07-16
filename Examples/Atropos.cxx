#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMaskImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "antsAtroposSegmentationImageFilter.h"
#include "antsBoxPlotQuantileListSampleFilter.h"
#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"
#include "antsGaussianListSampleFunction.h"
#include "antsGrubbsRosnerListSampleFilter.h"
#include "antsHistogramParzenWindowsListSampleFunction.h"
#include "antsListSampleToListSampleFilter.h"
#include "antsManifoldParzenWindowsListSampleFunction.h"
#include "antsPassThroughListSampleFilter.h"

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
int AtroposSegmentation( itk::ants::CommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;

  typedef unsigned char                         LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef  itk::ants::AtroposSegmentationImageFilter
    <InputImageType, LabelImageType> SegmentationFilterType;
  typename SegmentationFilterType::Pointer segmenter
    = SegmentationFilterType::New();

  typedef CommandIterationUpdate<SegmentationFilterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  segmenter->AddObserver( itk::IterationEvent(), observer );
  // segmenter->DebugOn();

  /**
   * memory-usage -- need to set before setting the prior probability images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer memoryOption =
    parser->GetOption( "minimize-memory-usage" );
  if( memoryOption && memoryOption->GetNumberOfValues() > 0 )
    {
    segmenter->SetMinimizeMemoryUsage( parser->Convert<bool>(
                                         memoryOption->GetValue() ) );
    }

  /**
   * Initialization
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer initializationOption =
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
      if( initializationOption->GetNumberOfParameters() > 3 )
        {
        segmenter->SetPriorProbabilityThreshold( parser->Convert<float>(
                                                   initializationOption->GetParameter( 3 ) ) );
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
  typename itk::ants::CommandLineParser::OptionType::Pointer convergenceOption =
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
   * mask image
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer maskOption =
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
  else
    {
    std::cerr << "No mask images was specified.  Specify a mask image"
              << " with the -x option." << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * BSpline options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer bsplineOption =
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
  typename itk::ants::CommandLineParser::OptionType::Pointer labelOption =
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
  typename itk::ants::CommandLineParser::OptionType::Pointer imageOption =
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
              << " with the -a option." << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * MRF options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer mrfOption =
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
  typename itk::ants::CommandLineParser::OptionType::Pointer distanceOption =
    parser->GetOption( "use-euclidean-distance" );
  if( distanceOption && distanceOption->GetNumberOfValues() > 0 )
    {
    segmenter->SetUseEuclideanDistanceForPriorLabels(
      parser->Convert<bool>( distanceOption->GetValue() ) );
    }

  /**
   * likelihood
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer likelihoodOption =
    parser->GetOption( "likelihood-model" );
  if( likelihoodOption && likelihoodOption->GetNumberOfValues() > 0 )
    {
    std::string likelihoodModel = likelihoodOption->GetValue();
    ConvertToLowerCase( likelihoodModel );
    if( !likelihoodModel.compare( std::string( "gaussian" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::GaussianListSampleFunction
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
      typedef itk::ants::Statistics::ManifoldParzenWindowsListSampleFunction
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
      typedef itk::ants::Statistics::HistogramParzenWindowsListSampleFunction
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
  typename itk::ants::CommandLineParser::OptionType::Pointer outlierOption =
    parser->GetOption( "winsorize-outliers" );
  if( outlierOption && outlierOption->GetNumberOfValues() > 0 )
    {
    std::string outlierStrategy = outlierOption->GetValue();
    ConvertToLowerCase( outlierStrategy );
    if( !outlierStrategy.compare( std::string( "boxplot" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::BoxPlotQuantileListSampleFilter<SampleType>
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
      typedef itk::ants::Statistics::GrubbsRosnerListSampleFilter<SampleType>
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
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
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

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "One or more scalar images is specified for segmentation " )
      + std::string( "using the -a/--intensity-image option.  For segmentation " )
      + std::string( "scenarios with no prior information, the first scalar " )
      + std::string( "image encountered on the command line is used to order " )
      + std::string( "labelings such that the class with the smallest intensity " )
      + std::string( "signature is class \'1\' through class \'N\' which represents " )
      + std::string( "the voxels with the largest intensity values.  The " )
      + std::string( "optional adaptive smoothing weight parameter is applicable " )
      + std::string( "only when using prior label or probability images.  This " )
      + std::string( "scalar parameter is to be specified between [0,1] which " )
      + std::string( "smooths each labeled region separately and modulates the " )
      + std::string( "intensity measurement at each voxel in each intensity image " )
      + std::string( "between the original intensity and its smoothed " )
      + std::string( "counterpart.  The smoothness parameters are governed by the " )
      + std::string( "-b/--bspline option." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "intensity-image" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "[intensityImage,<adaptiveSmoothingWeight>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "If the adaptive smoothing weights are > 0, the intensity " )
      + std::string( "images are smoothed in calculating the likelihood values. " )
      + std::string( "This is to account for subtle intensity differences " )
      + std::string( "across the same tissue regions." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bspline" );
    option->SetShortName( 'b' );
    option->SetUsageOption( 0,
                            "[<numberOfLevels=6>,<initialMeshResolution=1x1x...>,<splineOrder=3>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "To initialize the FMM parameters, one of the following " )
      + std::string( "options must be specified.  If one does not have " )
      + std::string( "prior label or probability images we recommend " )
      + std::string( "using kmeans as it is typically faster than otsu and can " )
      + std::string( "be used with multivariate initialization. " )
      + std::string( "Random initialization is meant purely for intellectual " )
      + std::string( "curiosity. The prior weighting (specified in the range " )
      + std::string( "[0,1]) is used to modulate the calculation of the " )
      + std::string( "posterior probabilities between the likelihood*mrfprior " )
      + std::string( "and the likelihood*mrfprior*prior.  For specifying many " )
      + std::string( "prior probability images for a multi-label segmentation, " )
      + std::string( "we offer a minimize usage option (see -m).  With that option " )
      + std::string( "one can specify a prior probability threshold in which only " )
      + std::string( "those pixels exceeding that threshold are stored in memory. ");

//     std::string( "\t  Usage: \n" ) +
//     std::string( "\t    Option 1:  Random[numberOfClasses]\n" ) +
//     std::string( "\t    Option 2:  Kmeans[numberOfClasses]\n" ) +
//     std::string( "\t    Option 3:  Otsu[numberOfClasses]\n" ) +
//     std::string( "\t    Option 4:  PriorProbabilityImages[numberOfClasses," ) +
//     std::string( "fileSeriesFormat(index=1 to numberOfClasses) or vectorImage,priorWeighting]\n" ) +
//     std::string( "\t    Option 5:  PriorLabelImage[numberOfClasses,labelImage,priorWeighting]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initialization" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "Random[numberOfClasses]" );
    option->SetUsageOption( 1, "KMeans[numberOfClasses]" );
    option->SetUsageOption( 2, "Otsu[numberOfClasses]" );
    option->SetUsageOption( 3,
                            "PriorProbabilityImages[numberOfClasses,fileSeriesFormat(index=1 to numberOfClasses) or vectorImage,priorWeighting,<priorProbabilityThreshold>]" );
    option->SetUsageOption( 4, "PriorLabelImage[numberOfClasses,labelImage,priorWeighting]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The image mask defines the region which is to be labeled by " )
      + std::string( "the Atropos algorithm." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Convergence is determined by calculating the mean maximum " )
      + std::string( "posterior probability over the region of interest at " )
      + std::string( "each iteration. When this value decreases or increases " )
      + std::string( "less than the specified threshold from the previous " )
      + std::string( "iteration the program terminates.");

//     std::string( "\t  Usage: \n" ) +
//     std::string( "\t    [<numberOfIterations=5>,<convergenceThreshold=0.001>]" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "convergence" );
    option->SetShortName( 'c' );
    option->SetUsageOption( 0, "[<numberOfIterations=5>,<convergenceThreshold=0.001>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Both parametric and non-parametric options exist in Atropos. " )
      + std::string( "The Gaussian parametric option is commonly used " )
      + std::string( "(e.g. SPM & FAST) where the mean and standard deviation " )
      + std::string( "for the Gaussian of each class is calculated at each " )
      + std::string( "iteration.  Other groups use non-parametric approaches " )
      + std::string( "exemplified by option 2.  We recommend using options 1 " )
      + std::string( "or 2 as they are fairly standard and the " )
      + std::string( "default parameters work adequately." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "likelihood-model" );
    option->SetShortName( 'k' );
    option->SetUsageOption( 0, "Gaussian" );
    option->SetUsageOption( 1, "HistogramParzenWindows[<Sigma=1.0>,<numberOfBins=32>]" );
    option->SetUsageOption( 2,
                            "ManifoldParzenWindows[<pointSetSigma=1.0>,<evaluationKNeighborhood=50>,<CovarianceKNeighborhood=0>,<kernelSigma=0>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Markov random field (MRF) theory provides a general " )
      + std::string( "framework for enforcing spatially contextual constraints " )
      + std::string( "on the segmentation solution.  The default smoothing " )
      + std::string( "factor of 0.3 provides a moderate amount of smoothing. " )
      + std::string( "Increasing this number causes more smoothing whereas " )
      + std::string( "decreasing the number lessens the smoothing. The radius " )
      + std::string( "parameter specifies the mrf neighborhood." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mrf" );
    option->SetShortName( 'm' );
    option->SetUsageOption( 0, "[<smoothingFactor=0.3>,<radius=1x1x...>]" );
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
    std::string description =
      std::string( "The output consists of a labeled image where each voxel " )
      + std::string( "in the masked region is assigned a label from 1, 2, " )
      + std::string( "..., N.  Optionally, one can also output the posterior " )
      + std::string( "probability images specified in the same format as the " )
      + std::string( "prior probability images, e.g. posterior%02d.nii.gz " )
      + std::string( "(C-style file name formatting)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "[classifiedImage,<posteriorProbabilityImageFileNameFormat>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "By default, memory usage is not minimized, however, if " )
      + std::string( "this is needed, the various probability and distance " )
      + std::string( "images are calculated on the fly instead of being " )
      + std::string( "stored in memory at each iteration. Also, if prior " )
      + std::string( "probability images are used, only the non-negligible " )
      + std::string( "pixel values are stored in memory. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "minimize-memory-usage" );
    option->SetShortName( 'u' );
    option->SetUsageOption( 0, "(0)/1" );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "To remove the effects of outliers in calculating the " )
      + std::string( "weighted mean and weighted covariance, the user can " )
      + std::string( "opt to remove the outliers through the options " )
      + std::string( "specified below." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "winsorize-outliers" );
    option->SetShortName( 'w' );
    option->SetUsageOption( 0, "BoxPlot[<lowerPercentile=0.25>,<upperPercentile=0.75>,<whiskerLength=1.5>]" );
    option->SetUsageOption( 1, "GrubbsRosner[<significanceLevel=0.05>,<winsorizingLevel=0.10>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Given prior label or probability images, the labels are " )
      + std::string( "propagated throughout the masked region so that every " )
      + std::string( "voxel in the mask is labeled.  Propagation is done " )
      + std::string( "by using a signed distance transform of the label. " )
      + std::string( "Alternatively, propagation of the labels with the " )
      + std::string( "fast marching filter respects the distance along the " )
      + std::string( "shape of the mask (e.g. the sinuous sulci and gyri " )
      + std::string( "of the cortex." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "use-euclidean-distance" );
    option->SetShortName( 'e' );
    option->SetUsageOption( 0, "(0)/1" );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The propagation of each prior label can be controlled " )
      + std::string( "by the sigma and boundary probability parameters.  The " )
      + std::string( "latter parameter is the probability (in the range " )
      + std::string( "[0,1]) of the label on the boundary which increases linearly " )
      + std::string( "to a maximum value of 1.0 in the interior of the labeled " )
      + std::string( "region.  The former parameter dictates the exponential " )
      + std::string( "decay of probability propagation outside the labeled " )
      + std::string( "region from the boundary probability, i.e. " )
      + std::string( "boundaryProbability*exp( -distance / sigma^2 )." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "label-propagation" );
    option->SetShortName( 'l' );
    option->SetUsageOption( 0, "whichLabel[sigma=0.0,<boundaryProbability=1.0>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, Atropos tries to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "image-dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3/4" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }
}

int main( int argc, char *argv[] )
{
  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "A finite mixture modeling (FMM) segmentation approach " )
    + std::string( "with possibilities for specifying prior constraints. " )
    + std::string( "These prior constraints include the specification " )
    + std::string( "of a prior label image, prior probability images " )
    + std::string( "(one for each class), and/or an MRF prior to " )
    + std::string( "enforce spatial smoothing of the labels.  Similar algorithms " )
    + std::string( "include FAST and SPM.  " );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>(
        parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    exit( EXIT_FAILURE );
    }
  else if( parser->Convert<bool>(
             parser->GetOption( 'h' )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    exit( EXIT_FAILURE );
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetValue() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "intensity-image" );
    if( imageOption && imageOption->GetNumberOfValues() > 0 )
      {
      if( imageOption->GetNumberOfParameters( 0 ) > 0 )
        {
        filename = imageOption->GetParameter( 0, 0 );
        }
      else
        {
        filename = imageOption->GetValue( 0 );
        }
      }
    else
      {
      std::cerr << "No input images were specified.  Specify an input image"
                << " with the -a option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  switch( dimension )
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
