#include <string>
#include <algorithm>
#include <vector>
#include <algorithm>
#include "antsUtilities.h"
#include "ReadWriteData.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMaskImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "antsAtroposSegmentationImageFilter.h"
#include "antsBoxPlotQuantileListSampleFilter.h"
#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"
#include "antsGaussianListSampleFunction.h"
#include "antsLogEuclideanGaussianListSampleFunction.h"
#include "antsGrubbsRosnerListSampleFilter.h"
#include "antsHistogramParzenWindowsListSampleFunction.h"
#include "antsJointHistogramParzenShapeAndOrientationListSampleFunction.h"
#include "antsListSampleToListSampleFilter.h"
#include "antsManifoldParzenWindowsListSampleFunction.h"
#include "antsPassThroughListSampleFilter.h"
#include "antsPartialVolumeGaussianListSampleFunction.h"
#include "itkTimeProbe.h"

namespace ants
{
template <typename TFilter>
class CommandIterationUpdate final : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() = default;
public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    const auto * filter =
      dynamic_cast<const TFilter *>( object );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    std::cout << "  Iteration " << filter->GetElapsedIterations()
             << " (of " << filter->GetMaximumNumberOfIterations() << "): ";
    std::cout << "posterior probability = "
             << filter->GetCurrentPosteriorProbability();

    typedef typename TFilter::RealType RealType;

    RealType annealingTemperature = filter->GetInitialAnnealingTemperature()
      * std::pow( filter->GetAnnealingRate(), static_cast<RealType>(
                   filter->GetElapsedIterations() ) );

    annealingTemperature = std::max( annealingTemperature,
                                         filter->GetMinimumAnnealingTemperature() );

    std::cout << " (annealing temperature = "
             << annealingTemperature << ")" << std::endl;
  }
};

template <unsigned int ImageDimension>
int AtroposSegmentation( itk::ants::CommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;

  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  bool verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( verbose )
    {
    std::cout << std::endl << "Running Atropos for "
             << ImageDimension << "-dimensional images." << std::endl;
    }

  typedef  itk::ants::AtroposSegmentationImageFilter
    <InputImageType, LabelImageType> SegmentationFilterType;
  typename SegmentationFilterType::Pointer segmenter
    = SegmentationFilterType::New();

  if( verbose )
    {
    typedef CommandIterationUpdate<SegmentationFilterType> CommandType;
    typename CommandType::Pointer observer = CommandType::New();
    segmenter->AddObserver( itk::IterationEvent(), observer );
    }

  /**
   * memory-usage -- need to set before setting the prior probability images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer memoryOption =
    parser->GetOption( "minimize-memory-usage" );
  if( memoryOption && memoryOption->GetNumberOfFunctions() )
    {
    segmenter->SetMinimizeMemoryUsage( parser->Convert<bool>(
                                         memoryOption->GetFunction( 0 )->GetName() ) );
    }

  /**
   * Initialization
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer initializationOption =
    parser->GetOption( "initialization" );
  if( initializationOption && initializationOption->GetNumberOfFunctions() &&
      initializationOption->GetFunction( 0 )->GetNumberOfParameters() < 1 )
    {
    if( verbose )
      {
      std::cerr << "Incorrect initialization option specification." << std::endl;
      std::cerr << "   " << initializationOption->GetDescription() << std::endl;
      }
    return EXIT_FAILURE;
    }
  else
    {
    segmenter->SetNumberOfTissueClasses( parser->Convert<unsigned int>(
                                           initializationOption->GetFunction( 0 )->GetParameter( 0 ) ) );

    std::string initializationStrategy = initializationOption->GetFunction( 0 )->GetName();
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
      if( initializationOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        std::vector<float> clusterCenters = parser->ConvertVector<float>(
            initializationOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( clusterCenters.size() != segmenter->GetNumberOfTissueClasses() )
          {
          if( verbose )
            {
            std::cerr << "The cluster center vector size does not equal the "
                     << "specified number of classes." << std::endl;
            }
          return EXIT_FAILURE;
          }
        else
          {
          typename SegmentationFilterType::ParametersType parameters;
          parameters.SetSize( segmenter->GetNumberOfTissueClasses() );
          for( unsigned int n = 0; n < parameters.GetSize(); n++ )
            {
            parameters[n] = clusterCenters[n];
            }
          segmenter->SetInitialKMeansParameters( parameters );
          }
        }
      }
    else if( !initializationStrategy.compare(
               std::string( "priorprobabilityimages" ) ) )
      {
      segmenter->SetInitializationStrategy(
        SegmentationFilterType::PriorProbabilityImages );
      if( initializationOption->GetFunction( 0 )->GetNumberOfParameters() < 3 )
        {
        if( verbose )
          {
          std::cerr << "Incorrect initialization option specification."
                   << std::endl;
          std::cerr << "   " << initializationOption->GetDescription()
                   << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetFunction( 0 )->GetParameter( 2 ) ) );
      if( initializationOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        segmenter->SetProbabilityThreshold( parser->Convert<float>(
                                              initializationOption->GetFunction( 0 )->GetParameter( 3 ) ) );
        }

      std::string filename = initializationOption->GetFunction( 0 )->GetParameter( 1 );

      if( filename.find( std::string( "%" ) ) != std::string::npos )
        {
        itk::NumericSeriesFileNames::Pointer fileNamesCreator =
          itk::NumericSeriesFileNames::New();
        fileNamesCreator->SetStartIndex( 1 );
        fileNamesCreator->SetEndIndex( segmenter->GetNumberOfTissueClasses() );
        fileNamesCreator->SetSeriesFormat( filename.c_str() );
        const std::vector<std::string> & imageNames
          = fileNamesCreator->GetFileNames();
        for( unsigned int k = 0; k < imageNames.size(); k++ )
          {
          typename InputImageType::Pointer image;
          ReadImage<InputImageType>( image, imageNames[k].c_str() );
          segmenter->SetPriorProbabilityImage( k + 1, image );
          }
        }
      else
        {
        typedef itk::VectorImage<PixelType, ImageDimension> VectorImageType;
        typename VectorImageType::Pointer image;
        ReadImage<VectorImageType>( image, filename.c_str() );

        if(  image->GetNumberOfComponentsPerPixel()
             != segmenter->GetNumberOfTissueClasses() )
          {
          if( verbose )
            {
            std::cerr << "The number of components does not match the number of "
                     << "classes." << std::endl;
            }
          return EXIT_FAILURE;
          }

        typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, InputImageType> CasterType;
        typename CasterType::Pointer caster = CasterType::New();
        caster->SetInput( image );
        for( unsigned int k = 0; k < segmenter->GetNumberOfTissueClasses(); k++ )
          {
          caster->SetIndex( k );
          caster->Update();
          segmenter->SetPriorProbabilityImage( k + 1, caster->GetOutput() );
          }
        }
      if( initializationOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        segmenter->SetProbabilityThreshold( parser->Convert<float>(
                                              initializationOption->GetFunction( 0 )->GetParameter( 3 ) ) );
        }
      }
    else if( !initializationStrategy.compare( std::string( "priorlabelimage" ) ) )
      {
      segmenter->SetInitializationStrategy(
        SegmentationFilterType::PriorLabelImage );

      if( initializationOption->GetFunction( 0 )->GetNumberOfParameters() < 3 )
        {
        if( verbose )
          {
          std::cerr << "Incorrect initialization option specification." << std::endl;
          std::cerr << "   " << initializationOption->GetDescription() << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->SetPriorProbabilityWeight( parser->Convert<float>(
                                              initializationOption->GetFunction( 0 )->GetParameter( 2 ) ) );

      std::string filename = initializationOption->GetFunction( 0 )->GetParameter( 1 );
      typename LabelImageType::Pointer image;
      ReadImage<LabelImageType>( image, filename.c_str() );
      segmenter->SetPriorLabelImage( image );
      }
    else
      {
      if( verbose )
        {
        std::cerr << "Unrecognized initialization strategy request." << std::endl;
        }
      return EXIT_FAILURE;
      }
    }

  /**
   * Posterior probability formulation
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer posteriorOption =
    parser->GetOption( "posterior-formulation" );
  if( posteriorOption && posteriorOption->GetNumberOfFunctions() )
    {
    if( posteriorOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      segmenter->SetUseMixtureModelProportions( parser->Convert<bool>(
                                                  posteriorOption->GetFunction( 0 )->GetParameter( 0 ) ) );

      RealType annealingTemperature = 1.0;
      if( posteriorOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        annealingTemperature =
          parser->Convert<RealType>( posteriorOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( annealingTemperature <= itk::NumericTraits<RealType>::ZeroValue() )
          {
          if( verbose )
            {
            std::cerr << "Annealing temperature must be positive." << std::endl;
            }
          return EXIT_FAILURE;
          }
        }
      segmenter->SetInitialAnnealingTemperature( annealingTemperature );

      RealType annealingRate = 1.0;
      if( posteriorOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        annealingRate =
          parser->Convert<RealType>( posteriorOption->GetFunction( 0 )->GetParameter( 2 ) );
        if( annealingRate < itk::NumericTraits<RealType>::ZeroValue() ||
            annealingRate > itk::NumericTraits<RealType>::OneValue() )
          {
          if( verbose )
            {
            std::cerr << "Annealing rate must be in the range [0, 1]." << std::endl;
            }
          return EXIT_FAILURE;
          }
        }
      segmenter->SetAnnealingRate( annealingRate );

      if( posteriorOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        auto minimumAnnealingTemperature =
          parser->Convert<RealType>( posteriorOption->GetFunction( 0 )->GetParameter( 3 ) );
        segmenter->SetMinimumAnnealingTemperature( minimumAnnealingTemperature );
        }
      }
    std::string posteriorStrategy = posteriorOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( posteriorStrategy );

    if( !posteriorStrategy.compare( std::string( "socrates" ) ) )
      {
      segmenter->SetPosteriorProbabilityFormulation(
        SegmentationFilterType::Socrates );
      }
    else if( !posteriorStrategy.compare( std::string( "plato" ) ) )
      {
      segmenter->SetPosteriorProbabilityFormulation(
        SegmentationFilterType::Plato );
      }
    else if( !posteriorStrategy.compare( std::string( "aristotle" ) ) )
      {
      segmenter->SetPosteriorProbabilityFormulation(
        SegmentationFilterType::Aristotle );
      }
    else if( !posteriorStrategy.compare( std::string( "sigmoid" ) ) )
      {
      segmenter->SetPosteriorProbabilityFormulation(
        SegmentationFilterType::Sigmoid );
      }
    }

  /**
   * convergence options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence" );
  if( convergenceOption && convergenceOption->GetNumberOfFunctions() )
    {
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      segmenter->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
                                                 convergenceOption->GetFunction( 0 )->GetName() ) );
      }
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      segmenter->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
                                                 convergenceOption->GetFunction( 0 )->GetParameter( 0 ) ) );
      }
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      segmenter->SetConvergenceThreshold( parser->Convert<float>(
                                            convergenceOption->GetFunction( 0 )->GetParameter( 1 ) ) );
      }
    }

  /**
   * mask image
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer maskOption =
    parser->GetOption( "mask-image" );
  if( maskOption && maskOption->GetNumberOfFunctions() )
    {
    try
      {
      typename LabelImageType::Pointer image;
      ReadImage<LabelImageType>( image, maskOption->GetFunction( 0 )->GetName().c_str()  );
      segmenter->SetMaskImage( image );

      // Check to see that the labels in the prior label image or the non-zero
      // probability voxels in the prior probability images encompass the entire
      // mask region.

      if( segmenter->GetInitializationStrategy() ==
          SegmentationFilterType::PriorLabelImage )
        {
        itk::ImageRegionConstIterator<LabelImageType> ItM( segmenter->GetMaskImage(),
                                                           segmenter->GetMaskImage()->GetLargestPossibleRegion() );
        itk::ImageRegionConstIterator<LabelImageType> ItP( segmenter->GetPriorLabelImage(),
                                                           segmenter->GetPriorLabelImage()->GetLargestPossibleRegion() );
        for( ItM.GoToBegin(), ItP.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItP )
          {
          if( ItM.Get() != itk::NumericTraits<LabelType>::ZeroValue() && ItP.Get() == 0 )
            {
            if( verbose )
              {
              std::cout << std::endl;
              std::cout << "Warning: the labels in the the prior label image do "
                       << "not encompass the entire mask region.  As a result each unlabeled voxel will be "
                       << "initially assigned a random label.  The user might want to consider "
                       << "various alternative strategies like assigning an additional "
                       << "\"background\" label to the unlabeled voxels or propagating "
                       << "the labels within the mask region."
                       << std::endl;
              std::cout << std::endl;
              }
            break;
            }
          }
        }
      else if( segmenter->GetInitializationStrategy() ==
               SegmentationFilterType::PriorProbabilityImages )
        {
        itk::ImageRegionConstIteratorWithIndex<LabelImageType> ItM( segmenter->GetMaskImage(),
                                                                    segmenter->GetMaskImage()->GetLargestPossibleRegion() );
        for( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
          {
          if( ItM.Get() != itk::NumericTraits<LabelType>::ZeroValue() )
            {
            RealType sumPriorProbability = 0.0;
            for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
              {
              sumPriorProbability +=
                segmenter->GetPriorProbabilityImage( n + 1 )->GetPixel( ItM.GetIndex() );
              }
            if( sumPriorProbability < segmenter->GetProbabilityThreshold() )
              {
              if( verbose )
                {
                std::cout << std::endl;
                std::cout << "Warning: the sum of the priors from the the prior probability images are "
                         << "less than the probability threshold within the mask region.  As a result "
                         << "each zero probability voxel will be "
                         << "initially assigned a random label.  The user might want to consider "
                         << "various alternative strategies like assigning an additional "
                         << "\"background\" label to the zero probability voxels or propagating "
                         << "the probabilities within the mask region."
                         << std::endl;
                std::cout << std::endl;
                }
              break;
              }
            }
          }
        }
      }
    catch( ... )
      {
      }
    }
  else
    {
    if( verbose )
      {
      std::cerr << "An image mask is required.  Specify a mask image"
               << " with the -x option." << std::endl;
      }
    return EXIT_FAILURE;
    }

  /**
   * BSpline options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer bsplineOption =
    parser->GetOption( "bspline" );
  if( bsplineOption && bsplineOption->GetNumberOfFunctions() )
    {
    if( bsplineOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      std::vector<unsigned int> numLevels = parser->ConvertVector<unsigned int>(
          bsplineOption->GetFunction( 0 )->GetParameter( 0 ) );
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
        if( verbose )
          {
          std::cerr << "Incorrect number of levels" << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->SetNumberOfLevels( numberOfFittingLevels );
      }
    if( bsplineOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      segmenter->SetSplineOrder( parser->Convert<unsigned int>(
                                   bsplineOption->GetFunction( 0 )->GetParameter( 2 ) ) );
      }
    if( bsplineOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      std::vector<unsigned int> array = parser->ConvertVector<unsigned int>(
          bsplineOption->GetFunction( 0 )->GetParameter( 1 ) );
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
        if( verbose )
          {
          std::cerr << "Incorrect mesh resolution" << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->SetNumberOfControlPoints( numberOfControlPoints );
      }
    }

  /**
   * labels
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer labelOption =
    parser->GetOption( "label-propagation" );
  if( labelOption && labelOption->GetNumberOfFunctions() )
    {
    if( labelOption->GetNumberOfFunctions() == 1 &&
        ( labelOption->GetFunction( 0 )->GetName() ).empty() )
      {
      typename SegmentationFilterType::LabelParameterMapType labelMap;

      auto labelLambda = parser->Convert<float>(
          labelOption->GetFunction( 0 )->GetParameter( 0 ) );
      float labelBoundaryProbability = 1.0;
      if( labelOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        labelBoundaryProbability = parser->Convert<float>(
            labelOption->GetFunction( 0 )->GetParameter( 1 ) );
        if( labelBoundaryProbability < itk::NumericTraits<float>::ZeroValue() )
          {
          labelBoundaryProbability = itk::NumericTraits<float>::ZeroValue();
          }
        if( labelBoundaryProbability > itk::NumericTraits<float>::OneValue() )
          {
          labelBoundaryProbability = itk::NumericTraits<float>::OneValue();
          }
        }
      for( unsigned int n = 1; n <= segmenter->GetNumberOfTissueClasses(); n++ )
        {
        typename SegmentationFilterType::LabelParametersType labelPair;
        labelPair.first = labelLambda;
        labelPair.second = labelBoundaryProbability;
        labelMap[n] = labelPair;
        }
      segmenter->SetPriorLabelParameterMap( labelMap );
      }
    else
      {
      typename SegmentationFilterType::LabelParameterMapType labelMap;
      for( unsigned int n = 0; n < labelOption->GetNumberOfFunctions(); n++ )
        {
        typename SegmentationFilterType::LabelParametersType labelPair;

        auto labelLambda = parser->Convert<float>(
            labelOption->GetFunction( n )->GetParameter( 0 ) );
        float labelBoundaryProbability = 1.0;
        if( labelOption->GetFunction( n )->GetNumberOfParameters() > 1 )
          {
          labelBoundaryProbability = parser->Convert<float>( labelOption->GetFunction( n )->GetParameter( 1 ) );
          if( labelBoundaryProbability < itk::NumericTraits<float>::ZeroValue() )
            {
            labelBoundaryProbability = itk::NumericTraits<float>::ZeroValue();
            }
          if( labelBoundaryProbability > itk::NumericTraits<float>::OneValue() )
            {
            labelBoundaryProbability = itk::NumericTraits<float>::OneValue();
            }
          }
        labelPair.first = labelLambda;
        labelPair.second = labelBoundaryProbability;

        auto whichClass = parser->Convert<unsigned int>( labelOption->GetFunction( n )->GetName() );

        labelMap[whichClass] = labelPair;
        }
      segmenter->SetPriorLabelParameterMap( labelMap );
      }
    }

  /**
   * intensity images
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "intensity-image" );
  if( imageOption && imageOption->GetNumberOfFunctions() )
    {
    unsigned int count = 0;
    for( int n = imageOption->GetNumberOfFunctions() - 1; n >= 0; n-- )
      {
      std::string imagename;
      if( imageOption->GetFunction( n )->GetNumberOfParameters() > 0 )
        {
        imagename = imageOption->GetFunction( n )->GetParameter( 0 );
        }
      else
        {
        imagename = imageOption->GetFunction( n )->GetName();
        }
      typename InputImageType::Pointer image;
      ReadImage<InputImageType>( image, imagename.c_str() );
      segmenter->SetIntensityImage( count, image );
      if( imageOption->GetFunction( count )->GetNumberOfParameters() > 1 )
        {
        segmenter->SetAdaptiveSmoothingWeight( count,
                                               parser->Convert<float>( imageOption->GetFunction( count )->GetParameter(
                                                                         1 ) ) );
        }
      else
        {
        segmenter->SetAdaptiveSmoothingWeight( count, 0.0 );
        }
      count++;
      }
    }
  else
    {
    if( verbose )
      {
      std::cerr << "No input images were specified.  Specify an input image"
               << " with the -a option." << std::endl;
      }
    return EXIT_FAILURE;
    }

  /**
   * MRF options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer mrfOption = parser->GetOption( "mrf" );
  if( mrfOption && mrfOption->GetNumberOfFunctions() )
    {
    if( mrfOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      try
        {
        typedef typename SegmentationFilterType::RealImageType
          MRFCoefficientImageType;
        typedef itk::ImageFileReader<MRFCoefficientImageType>
          MRFNeighborhoodImageReaderType;
        typename MRFNeighborhoodImageReaderType::Pointer mrfNeighborhoodReader =
          MRFNeighborhoodImageReaderType::New();
        mrfNeighborhoodReader->SetFileName( mrfOption->GetFunction( 0 )->GetParameter( 0 ) );

        typename MRFCoefficientImageType::Pointer mrfCoefficientImage =
          mrfNeighborhoodReader->GetOutput();
        mrfCoefficientImage->Update();
        mrfCoefficientImage->DisconnectPipeline();

        segmenter->SetMRFCoefficientImage( mrfCoefficientImage );
        }
      catch( ... )
        {
        segmenter->SetMRFSmoothingFactor( parser->Convert<float>( mrfOption->GetFunction( 0 )->GetParameter( 0 ) ) );
        }
      }
    if( mrfOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      std::vector<unsigned int> array =
        parser->ConvertVector<unsigned int>( mrfOption->GetFunction( 0 )->GetParameter( 1 ) );
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
        if( verbose )
          {
          std::cerr << "MRF radius size needs to be equal to the image dimension."
                   << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->SetMRFRadius( radius );
      }
    }

  /**
   * ICM options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer icmOption =
    parser->GetOption( "icm" );
  if( icmOption && icmOption->GetNumberOfFunctions() == 1 )
    {
    segmenter->SetUseAsynchronousUpdating( parser->Convert<bool>(
                                             icmOption->GetFunction( 0 )->GetName() ) );
    }
  if( icmOption && icmOption->GetNumberOfFunctions() )
    {
    if( icmOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      segmenter->SetUseAsynchronousUpdating( parser->Convert<bool>( icmOption->GetFunction( 0 )->GetParameter( 0 ) ) );
      }
    if( icmOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      segmenter->SetMaximumNumberOfICMIterations( parser->Convert<unsigned int>( icmOption->GetFunction( 0 )->
                                                                                 GetParameter( 1 ) ) );
      }
    }

  /**
   * random seed
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer seedOption =
    parser->GetOption( "use-random-seed" );
  if( seedOption && seedOption->GetNumberOfFunctions() )
    {
    bool useRandomSeed = parser->Convert<bool>( seedOption->GetFunction( 0 )->GetName() );
    if( !useRandomSeed )
      {
      // assign seed from itkMersenneTwisterRandomVariateGenerator.h (line 347)
      segmenter->SetRandomizerInitializationSeed( 19650218UL );
      }
    }

  /**
   * euclidean distance
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer distanceOption =
    parser->GetOption( "use-euclidean-distance" );
  if( distanceOption && distanceOption->GetNumberOfFunctions() )
    {
    segmenter->SetUseEuclideanDistanceForPriorLabels(
      parser->Convert<bool>( distanceOption->GetFunction( 0 )->GetName() ) );
    }

  /**
   * likelihood
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer likelihoodOption =
    parser->GetOption( "likelihood-model" );
  if( likelihoodOption && likelihoodOption->GetNumberOfFunctions() )
    {
    std::string likelihoodModel = likelihoodOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( likelihoodModel );
    if( !likelihoodModel.compare( std::string( "gaussian" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::GaussianListSampleFunction
        <SampleType, float, float> LikelihoodType;
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
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
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        regularizationSigma = parser->Convert<float>( likelihoodOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      unsigned int evalNeighborhood = 50;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        evalNeighborhood = parser->Convert<unsigned int>( likelihoodOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      unsigned int covNeighborhood = 0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        covNeighborhood = parser->Convert<unsigned int>( likelihoodOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      float covSigma = 1.0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        covSigma = parser->Convert<float>( likelihoodOption->GetFunction( 0 )->GetParameter( 3 ) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
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
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        sigma = parser->Convert<float>( likelihoodOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      unsigned int numberOfBins = 32;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        numberOfBins = parser->Convert<unsigned int>( likelihoodOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
        {
        typename LikelihoodType::Pointer hpwLikelihood =
          LikelihoodType::New();
        hpwLikelihood->SetSigma( sigma );
        hpwLikelihood->SetNumberOfHistogramBins( numberOfBins );
        segmenter->SetLikelihoodFunction( n, hpwLikelihood );
        }
      }
    else if( !likelihoodModel.compare( std::string( "jointshapeandorientationprobability" ) ) )
      {
      if( segmenter->GetNumberOfIntensityImages() !=
          static_cast<unsigned int>( ImageDimension * ( ImageDimension + 1 ) / 2 ) )
        {
        if( verbose )
          {
          std::cerr << " Expect images in upper triangular order " << std::endl;
          std::cerr << " xx xy xz yy yz zz " << std::endl;
          std::cerr << "Incorrect number of intensity images specified." << std::endl;
          }
        return EXIT_FAILURE;
        }
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::
        JointHistogramParzenShapeAndOrientationListSampleFunction
        <SampleType, float, float> LikelihoodType;

      float shapeSigma = 2.0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        shapeSigma = parser->Convert<float>( likelihoodOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      unsigned int numberOfShapeBins = 64;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        numberOfShapeBins = parser->Convert<unsigned int>( likelihoodOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      float orientationSigma = 1.0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        orientationSigma = parser->Convert<float>( likelihoodOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      unsigned int numberOfOrientationBins = 32;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        numberOfOrientationBins =
          parser->Convert<unsigned int>( likelihoodOption->GetFunction( 0 )->GetParameter( 3 ) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
        {
        typename LikelihoodType::Pointer hpwLikelihood =
          LikelihoodType::New();
        hpwLikelihood->SetShapeSigma( shapeSigma );
        hpwLikelihood->SetOrientationSigma( orientationSigma);
        hpwLikelihood->SetNumberOfShapeJointHistogramBins( numberOfShapeBins );
        hpwLikelihood->SetNumberOfOrientationJointHistogramBins( numberOfOrientationBins);
        segmenter->SetLikelihoodFunction( n, hpwLikelihood );
        }
      }
    else if( !likelihoodModel.compare( std::string( "logeuclideangaussian" ) ) )
      {
      if( segmenter->GetNumberOfIntensityImages() !=
          static_cast<unsigned int>( ImageDimension * ( ImageDimension + 1 ) / 2 ) )
        {
        if( verbose )
          {
          std::cerr << " Expect images in upper triangular order " << std::endl;
          std::cerr << " xx xy xz yy yz zz " << std::endl;
          std::cerr << "Incorrect number of intensity images specified." << std::endl;
          }
        return EXIT_FAILURE;
        }
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::LogEuclideanGaussianListSampleFunction
        <SampleType, float, float> LikelihoodType;
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
        {
        typename LikelihoodType::Pointer gaussianLikelihood =
          LikelihoodType::New();
        segmenter->SetLikelihoodFunction( n, gaussianLikelihood );
        }
      }
    else
      {
      if( verbose )
        {
        std::cerr << "Unrecognized likelihood model request." << std::endl;
        }
      return EXIT_FAILURE;
      }
    }

  /**
   * partial volume
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer pvOption =
    parser->GetOption( "partial-volume-label-set" );

  if( pvOption && pvOption->GetNumberOfFunctions() )
    {
    unsigned int labelSetCount = 0;
    for( int n = pvOption->GetNumberOfFunctions() - 1; n >= 0; n-- )
      {
      typename SegmentationFilterType::PartialVolumeLabelSetType labelSet =
        parser->ConvertVector<LabelType>( pvOption->GetFunction( n )->GetName() );
      if( labelSet.size() != 2 )
        {
        if( verbose )
          {
          std::cerr << "Error:  Currently Atropos only supports partial "
                   << "volume label sets of size equal to 2." << std::endl;
          }
        return EXIT_FAILURE;
        }
      segmenter->AddPartialVolumeLabelSet( labelSet );

      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::PartialVolumeGaussianListSampleFunction
        <SampleType, float, float> LikelihoodType;

      typename LikelihoodType::Pointer partialVolumeLikelihood =
        LikelihoodType::New();
      segmenter->SetLikelihoodFunction( labelSetCount
                                        + segmenter->GetNumberOfTissueClasses(), partialVolumeLikelihood );
      labelSetCount++;
      }

    typename itk::ants::CommandLineParser::OptionType::Pointer pvlOption =
      parser->GetOption( "use-partial-volume-likelihoods" );

    bool useLikelihoods = false;
    if( pvlOption && pvlOption->GetNumberOfFunctions() )
      {
      std::string value = pvlOption->GetFunction( 0 )->GetName();
      ConvertToLowerCase( value );
      if( !value.compare( "true" ) || !value.compare( "1" ) )
        {
        useLikelihoods = true;
        }
      else
        {
        useLikelihoods = false;
        }
      }
    segmenter->SetUsePartialVolumeLikelihoods( useLikelihoods );
    }

  /**
   * outliers handling
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outlierOption =
    parser->GetOption( "winsorize-outliers" );
  if( outlierOption && outlierOption->GetNumberOfFunctions() )
    {
    std::string outlierStrategy = outlierOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( outlierStrategy );
    if( !outlierStrategy.compare( std::string( "boxplot" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::BoxPlotQuantileListSampleFilter<SampleType>
        SampleFilterType;
      typename SampleFilterType::Pointer boxplotFilter =
        SampleFilterType::New();

      if( outlierOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        boxplotFilter->SetLowerPercentile( parser->Convert<float>( outlierOption->GetFunction( 0 )->GetParameter( 0 ) ) );
        }
      if( outlierOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        boxplotFilter->SetUpperPercentile( parser->Convert<float>( outlierOption->GetFunction( 0 )->GetParameter( 1 ) ) );
        }
      if( outlierOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        boxplotFilter->SetWhiskerScalingFactor( parser->Convert<float>( outlierOption->GetFunction( 0 )->GetParameter(
                                                                          2 ) ) );
        }
      segmenter->SetOutlierHandlingFilter( boxplotFilter );
      }
    else if( !outlierStrategy.compare( std::string( "grubbsrosner" ) ) )
      {
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::GrubbsRosnerListSampleFilter<SampleType>
        SampleFilterType;
      typename SampleFilterType::Pointer grubbsFilter = SampleFilterType::New();

      if( outlierOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        grubbsFilter->SetSignificanceLevel( parser->Convert<float>( outlierOption->GetFunction( 0 )->GetParameter( 0 ) ) );
        }
      if( outlierOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        grubbsFilter->SetWinsorizingLevel( parser->Convert<float>( outlierOption->GetFunction( 0 )->GetParameter( 1 ) ) );
        }
      segmenter->SetOutlierHandlingFilter( grubbsFilter );
      }
    else
      {
      if( verbose )
        {
        std::cerr << "Unrecognized outlier handling strategy request." << std::endl;
        }
      return EXIT_FAILURE;
      }
    }

  itk::TimeProbe timer;
  timer.Start();

  try
    {
    if( verbose )
      {
      std::cout << std::endl << "Progress: " << std::endl;
      }

//    segmenter->DebugOn();
    segmenter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    if( verbose )
      {
      std::cerr << exp << std::endl;
      }
    return EXIT_FAILURE;
    }

  timer.Stop();

  /**
   * output
   */
  if( icmOption && icmOption->GetNumberOfFunctions() && icmOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
    {
    if( segmenter->GetUseAsynchronousUpdating() && segmenter->GetICMCodeImage() )
      {
      typedef  itk::ImageFileWriter<LabelImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( segmenter->GetICMCodeImage() );
      writer->SetFileName( ( icmOption->GetFunction( 0 )->GetParameter( 2 ) ).c_str() );
      writer->Update();
      }
    }

  if( verbose )
    {
    std::cout << std::endl << "Writing output:" << std::endl;
    }
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      WriteImage<LabelImageType>(segmenter->GetOutput(), ( outputOption->GetFunction( 0 )->GetName() ).c_str() );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      WriteImage<LabelImageType>(segmenter->GetOutput(), ( outputOption->GetFunction( 0 )->GetParameter( 0 ) ).c_str() );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      std::string                          filename = outputOption->GetFunction( 0 )->GetParameter( 1 );
      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfTissueClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < imageNames.size(); i++ )
        {
        if( verbose )
          {
          std::cout << "  Writing posterior image (class " << i + 1 << ")" << std::endl;
          }
        typename InputImageType::Pointer probabilityImage
          = segmenter->GetPosteriorProbabilityImage( i + 1 );

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
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      std::string filename = outputOption->GetFunction( 0 )->GetParameter( 2 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfTissueClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < segmenter->GetNumberOfTissueClasses(); i++ )
        {
        if( verbose )
          {
          std::cout << "  Writing likelihood image (class " << i + 1 << ")" << std::endl;
          }
        typename InputImageType::Pointer likelihoodImage = segmenter->GetLikelihoodImage( i + 1 );
        typedef  itk::ImageFileWriter<InputImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( likelihoodImage );
        writer->SetFileName( imageNames[i].c_str() );
        writer->Update();
        }
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
      {
      std::string filename = outputOption->GetFunction( 0 )->GetParameter( 3 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfTissueClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < segmenter->GetNumberOfTissueClasses(); i++ )
        {
        if( segmenter->GetPriorProbabilityImage( i + 1 ) || segmenter->GetPriorLabelImage() )
          {
          if( verbose )
            {
            std::cout << "  Writing distance image (class " << i + 1 << ")" << std::endl;
            }

          typename InputImageType::Pointer distanceImage = segmenter->GetDistancePriorProbabilityImage( i + 1 );

          typedef  itk::ImageFileWriter<InputImageType> WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( distanceImage );
          writer->SetFileName( imageNames[i].c_str() );
          writer->Update();
          }
        }
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 4 )
      {
      std::string filename = outputOption->GetFunction( 0 )->GetParameter( 4 );

      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( segmenter->GetNumberOfTissueClasses() );
      fileNamesCreator->SetSeriesFormat( filename.c_str() );
      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();

      if( segmenter->GetAdaptiveSmoothingWeight( 0 ) > itk::NumericTraits<RealType>::ZeroValue() )
        {
        for( unsigned int i = 0; i < segmenter->GetNumberOfTissueClasses(); i++ )
          {
          if( segmenter->GetPriorProbabilityImage( i + 1 ) ||
              segmenter->GetPriorLabelImage() )
            {
            if( verbose )
              {
              std::cout << "  Writing B-spline image (class " << i + 1 << ")" << std::endl;
              }

            typename InputImageType::Pointer bsplineImage =
              segmenter->GetSmoothIntensityImageFromPriorImage( 0, i + 1 );

            typedef  itk::ImageFileWriter<InputImageType> WriterType;
            typename WriterType::Pointer writer = WriterType::New();
            writer->SetInput( bsplineImage );
            writer->SetFileName( imageNames[i].c_str() );
            writer->Update();
            }
          }
        }
      }
    }

  if( verbose )
    {
    std::cout << std::endl;
    segmenter->Print( std::cout, 2 );
    std::cout << "Elapsed time: " << timer.GetMean() << std::endl;
    }

  return EXIT_SUCCESS;
}

void AtroposInitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "This option forces the image to be treated as a specified-" )
    + std::string( "dimensional image.  If not specified, Atropos tries to " )
    + std::string( "infer the dimensionality from the first input image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3/4" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

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
    + std::string( "be used with multivariate initialization. However, since a " )
    + std::string( "Euclidean distance on the inter cluster distances is used, one " )
    + std::string( "might have to appropriately scale the additional input images. " )
    + std::string( "Random initialization is meant purely for intellectual " )
    + std::string( "curiosity. The prior weighting (specified in the range " )
    + std::string( "[0,1]) is used to modulate the calculation of the " )
    + std::string( "posterior probabilities between the likelihood*mrfprior " )
    + std::string( "and the likelihood*mrfprior*prior.  For specifying many " )
    + std::string( "prior probability images for a multi-label segmentation, " )
    + std::string( "we offer a minimize usage option (see -m).  With that option " )
    + std::string( "one can specify a prior probability threshold in which only " )
    + std::string( "those pixels exceeding that threshold are stored in memory. ");

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "initialization" );
  option->SetShortName( 'i' );
  option->SetUsageOption( 0, "Random[numberOfClasses]" );
  option->SetUsageOption( 1, "Otsu[numberOfTissueClasses]" );
  option->SetUsageOption(
    2, "KMeans[numberOfTissueClasses,<clusterCenters(in ascending order and for first intensity image only)>]" );
  option->SetUsageOption(
    3,
    "PriorProbabilityImages[numberOfTissueClasses,fileSeriesFormat(index=1 to numberOfClasses) or vectorImage,priorWeighting,<priorProbabilityThreshold>]" );
  option->SetUsageOption( 4, "PriorLabelImage[numberOfTissueClasses,labelImage,priorWeighting]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The partial volume estimation option allows one to model" )
    + std::string( "mixtures of classes within single voxels.  Atropos " )
    + std::string( "currently allows the user to model two class mixtures " )
    + std::string( "per partial volume class. The user specifies a set of " )
    + std::string( "class labels per partial volume class requested.  For " )
    + std::string( "example, suppose the user was performing a classic 3-" )
    + std::string( "tissue segmentation (csf, gm, wm) using kmeans " )
    + std::string( "initialization.  Suppose the user also wanted to model the " )
    + std::string( "partial voluming effects between csf/gm and gm/wm. " )
    + std::string( "The user would specify it using -i kmeans[3] " )
    + std::string( "and -s 1x2 -s 2x3.  So, for this example, there would be 3 " )
    + std::string( "tissue classes and 2 partial volume classes.  Optionally," )
    + std::string( "the user can limit partial volume handling to mrf considerations " )
    + std::string( "only whereby the output would only be the three tissues." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "partial-volume-label-set" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "label1xlabel2" );
  option->SetUsageOption( 0, "label1xlabel2xlabel3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The user can specify whether or not to use the partial " )
    + std::string( "volume likelihoods, in which case the partial volume class " )
    + std::string( "is considered separate from the tissue classes.  " )
    + std::string( "Alternatively, one can use the MRF only to handle partial " )
    + std::string( "volume in which case, partial volume voxels are not " )
    + std::string( "considered as separate classes." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-partial-volume-likelihoods" );
  option->SetUsageOption( 0, "1/(0)" );
  option->SetUsageOption( 1, "true/(false)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Different posterior probability formulations are possible as ")
    + std::string( "are different update options.  To guarantee theoretical " )
    + std::string( "convergence properties, a proper formulation of the well-known " )
    + std::string( "iterated conditional modes (ICM) uses an asynchronous update step " )
    + std::string( "modulated by a specified annealing temperature.  If one sets " )
    + std::string( "the AnnealingTemperature > 1 in the posterior formulation " )
    + std::string( "a traditional code set for a proper ICM update will be created. ")
    + std::string( "Otherwise, a synchronous update step will take place at each iteration. ")
    + std::string( "The annealing temperature, T, converts the posteriorProbability " )
    + std::string( "to posteriorProbability^(1/T) over the course of optimization. ");
  std::string( "Options include the following:  " )
  + std::string( " Socrates: posteriorProbability = (spatialPrior)^priorWeight" )
  + std::string( "*(likelihood*mrfPrior)^(1-priorWeight), " )
  + std::string( " Plato: posteriorProbability = 1.0, " )
  + std::string( " Aristotle: posteriorProbability = 1.0, " )
  + std::string( " Sigmoid: posteriorProbability = 1.0, " )/* +
  std::string( " Zeno: posteriorProbability = 1.0\n" ) +
  std::string( " Diogenes: posteriorProbability = 1.0\n" ) +
  std::string( " Thales: posteriorProbability = 1.0\n" ) +
  std::string( " Democritus: posteriorProbability = 1.0.\n" ) */;

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "posterior-formulation" );
  option->SetShortName( 'p' );
  option->SetUsageOption(
    0,
    "Socrates[<useMixtureModelProportions=1>,<initialAnnealingTemperature=1>,<annealingRate=1>,<minimumTemperature=0.1>]" );
  option->SetUsageOption(
    1,
    "Plato[<useMixtureModelProportions=1>,<initialAnnealingTemperature=1>,<annealingRate=1>,<minimumTemperature=0.1>]" );
  option->SetUsageOption(
    2,
    "Aristotle[<useMixtureModelProportions=1>,<initialAnnealingTemperature=1>,<annealingRate=1>,<minimumTemperature=0.1>]" );
  option->SetUsageOption(
    3,
    "Sigmoid[<useMixtureModelProportions=1>,<initialAnnealingTemperature=1>,<annealingRate=1>,<minimumTemperature=0.1>]]" );
//  option->SetUsageOption( 5, "Thales[<useMixtureModelProportions=1>]" );
//  option->SetUsageOption( 6, "Democritus" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The image mask (which is required) defines the region which " )
    + std::string( "is to be labeled by the Atropos algorithm." );

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
    + std::string( "iteration or the maximum number of iterations is exceeded " )
    + std::string( "the program terminates.");

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "convergence" );
  option->SetShortName( 'c' );
  option->SetUsageOption( 0, "numberOfIterations" );
  option->SetUsageOption( 1, "[<numberOfIterations=5>,<convergenceThreshold=0.001>]" );
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
  option->SetUsageOption( 1, "HistogramParzenWindows[<sigma=1.0>,<numberOfBins=32>]" );
  option->SetUsageOption(
    2,
    "ManifoldParzenWindows[<pointSetSigma=1.0>,<evaluationKNeighborhood=50>,<CovarianceKNeighborhood=0>,<kernelSigma=0>]" );
  option->SetUsageOption(
    3,
    "JointShapeAndOrientationProbability[<shapeSigma=1.0>,<numberOfShapeBins=64>, <orientationSigma=1.0>, <numberOfOrientationBins=32>]" );
  option->SetUsageOption( 4, "LogEuclideanGaussian" );
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
    + std::string( "parameter specifies the mrf neighborhood.  Different " )
    + std::string( "update schemes are possible but only the asynchronous " )
    + std::string( "updating has theoretical convergence properties. " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "mrf" );
  option->SetShortName( 'm' );
  option->SetUsageOption( 0, "[<smoothingFactor=0.3>,<radius=1x1x...>]" );
  option->SetUsageOption( 1, "[<mrfCoefficientImage>,<radius=1x1x...>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Asynchronous updating requires the construction of an " )
    + std::string( "ICM code image which is a label image (with labels in the " )
    + std::string( "range {1,..,MaximumICMCode}) constructed such that no MRF " )
    + std::string( "neighborhood has duplicate ICM code labels.  Thus, to update " )
    + std::string( "the voxel class labels we iterate through the code labels " )
    + std::string( "and, for each code label, we iterate through the image " )
    + std::string( "and update the voxel class label that has the corresponding " )
    + std::string( "ICM code label.  One can print out the ICM code image by " )
    + std::string( "specifying an ITK-compatible image filename." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "icm" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "[<useAsynchronousUpdate=1>,<maximumNumberOfICMIterations=1>,<icmCodeImage=''>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Initialize internal random number generator with a random seed. " ) +
    std::string( "Otherwise, initialize with a constant seed number." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-random-seed" );
  option->SetShortName( 'r' );
  option->SetUsageOption( 0, "0/(1)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

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
  option->AddFunction( std::string( "0" ) );
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
  option->AddFunction( std::string( "0" ) );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The propagation of each prior label can be controlled " )
    + std::string( "by the lambda and boundary probability parameters.  The " )
    + std::string( "latter parameter is the probability (in the range " )
    + std::string( "[0,1]) of the label on the boundary which increases linearly " )
    + std::string( "to a maximum value of 1.0 in the interior of the labeled " )
    + std::string( "region.  The former parameter dictates the exponential " )
    + std::string( "decay of probability propagation outside the labeled " )
    + std::string( "region from the boundary probability, i.e. " )
    + std::string( "boundaryProbability*exp( -lambda * distance )." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "label-propagation" );
  option->SetShortName( 'l' );
  option->SetUsageOption( 0, "whichLabel[lambda=0.0,<boundaryProbability=1.0>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Verbose output." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'v' );
  option->SetLongName( "verbose" );
  option->SetUsageOption( 0, "(0)/1" );
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
int Atropos( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "Atropos" );
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
    std::string( "A finite mixture modeling (FMM) segmentation approach " )
    + std::string( "with possibilities for specifying prior constraints. " )
    + std::string( "These prior constraints include the specification " )
    + std::string( "of a prior label image, prior probability images " )
    + std::string( "(one for each class), and/or an MRF prior to " )
    + std::string( "enforce spatial smoothing of the labels.  All segmentation " )
    + std::string( "images including priors and masks must be in the same " )
    + std::string( "voxel and physical space.  Similar algorithms include FAST " )
    + std::string( "and SPM.  Reference: " )
    + std::string( "Avants BB, Tustison NJ, Wu J, Cook PA, Gee JC. An open " )
    + std::string( "source multivariate framework for n-tissue segmentation " )
    + std::string( "with evaluation on public data. Neuroinformatics. " )
    + std::string( "2011 Dec;9(4):381-400.");

  parser->SetCommandDescription( commandDescription );
  AtroposInitializeCommandLineOptions( parser );

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
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "intensity-image" );
    if( imageOption && imageOption->GetNumberOfFunctions() )
      {
      if( imageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        filename = imageOption->GetFunction( 0 )->GetParameter( 0 );
        }
      else
        {
        filename = imageOption->GetFunction( 0 )->GetName();
        }
      }
    else
      {
      std::cerr << "No input images were specified.  Specify an input image"
               << " with the -a option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  switch( dimension )
    {
    case 2:
      return AtroposSegmentation<2>( parser );
      break;
    case 3:
      return AtroposSegmentation<3>( parser );
      break;
    case 4:
      return AtroposSegmentation<4>( parser );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
