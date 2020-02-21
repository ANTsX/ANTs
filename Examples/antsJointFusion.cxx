#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "ReadWriteData.h"

#include "itkNumericSeriesFileNames.h"
#include "itkTimeProbe.h"
#include "itkWeightedVotingFusionImageFilter.h"

#include "stdio.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "ANTsVersion.h"

namespace ants
{

template <typename TFilter>
class CommandProgressUpdate : public itk::Command
{
public:
  typedef  CommandProgressUpdate                      Self;
  typedef  itk::Command                               Superclass;
  typedef  itk::SmartPointer<CommandProgressUpdate>  Pointer;
  itkNewMacro( CommandProgressUpdate );
protected:

  CommandProgressUpdate() : m_CurrentProgress( 0 ), m_StartNewLine( true ) {};

  typedef TFilter FilterType;

  unsigned int m_CurrentProgress;
  bool         m_StartNewLine;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    const auto * filter = dynamic_cast<const TFilter *>( caller );

    if( this->m_CurrentProgress == 0 && ! filter->GetIsWeightedAveragingComplete() )
      {
      std::cout << "Weighted averaging: " << std::flush;
      }
    else if( this->m_StartNewLine && filter->GetIsWeightedAveragingComplete() )
      {
      std::cout << std::endl << "Reconstruction: " << std::flush;
      this->m_StartNewLine = false;
      this->m_CurrentProgress = 0;
      }

    auto *po = dynamic_cast<itk::ProcessObject *>( caller );
    if (! po) return;
//    std::cout << po->GetProgress() << std::endl;
    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    const auto * filter = dynamic_cast<const TFilter *>( object );

    if( this->m_CurrentProgress == 0 && ! filter->GetIsWeightedAveragingComplete() )
      {
      std::cout << "Weighted averaging: " << std::flush;
      }
    else if( this->m_StartNewLine && filter->GetIsWeightedAveragingComplete() )
      {
      std::cout << std::endl << "Reconstruction: " << std::flush;
      this->m_StartNewLine = false;
      this->m_CurrentProgress = 0;
      }

    auto *po = dynamic_cast<itk::ProcessObject *>(
      const_cast<itk::Object *>( object ) );
    if (! po) return;

    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {

      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }
};

template <unsigned int ImageDimension>
int antsJointFusion( itk::ants::CommandLineParser *parser )
{
  typedef float                                               RealType;
  typedef itk::Image<RealType, ImageDimension>                ImageType;
  typedef itk::Image<unsigned int, ImageDimension>            LabelImageType;
  typedef LabelImageType                                      MaskImageType;

  typedef typename itk::ants::CommandLineParser::OptionType   OptionType;

  // Determine verbosity of output

  bool verbose = false;
  typename OptionType::Pointer verboseOption = parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( verbose )
    {
    std::cout << std::endl << "Running antsJointFusion for "
             << ImageDimension << "-dimensional images." << std::endl << std::endl;
    }

  // Instantiate the joint fusion filter

  typedef itk::WeightedVotingFusionImageFilter<ImageType, LabelImageType> FusionFilterType;
  typename FusionFilterType::Pointer fusionFilter = FusionFilterType::New();
  typedef typename LabelImageType::PixelType                   LabelType;

  // Get the alpha and beta parameters

  RealType alpha = 0.1;
  typename OptionType::Pointer alphaOption = parser->GetOption( "alpha" );
  if( alphaOption && alphaOption->GetNumberOfFunctions() )
    {
    alpha = parser->Convert<RealType>( alphaOption->GetFunction( 0 )->GetName() );
    }

  RealType beta = 2.0;
  typename OptionType::Pointer betaOption = parser->GetOption( "beta" );
  if( betaOption && betaOption->GetNumberOfFunctions() )
    {
    beta = parser->Convert<RealType>( betaOption->GetFunction( 0 )->GetName() );
    }

  fusionFilter->SetAlpha( alpha );
  fusionFilter->SetBeta( beta );

  // Get the search and patch radii
  typename OptionType::Pointer searchRadiusOption = parser->GetOption( "search-radius" );

  if( searchRadiusOption && searchRadiusOption->GetNumberOfFunctions() )
    {
    // try reading the search radius as an image first.
    std::string searchRadiusString = searchRadiusOption->GetFunction( 0 )->GetName();
    if( itksys::SystemTools::FileExists( searchRadiusString.c_str() ) )
      {
      typedef typename FusionFilterType::RadiusImageType  RadiusImageType;
      typename RadiusImageType::Pointer searchRadiusImage;
      bool fileReadSuccessfully = ReadImage<RadiusImageType>( searchRadiusImage, searchRadiusString.c_str() );
      if( fileReadSuccessfully )
        {
        fusionFilter->SetNeighborhoodSearchRadiusImage( searchRadiusImage );
        }
      else
        {
        if( verbose )
          {
          std::cerr << "Search radius image exists but was not read successfully.." << std::endl;
          }
        return EXIT_FAILURE;
        }
      }
    else
      {
      std::vector<unsigned int> searchRadius;
      searchRadius.push_back( 3 );
      if( searchRadiusOption && searchRadiusOption->GetNumberOfFunctions() )
        {
        searchRadius = parser->ConvertVector<unsigned int>( searchRadiusString );
        }
      if( searchRadius.size() == 1 )
        {
        for( unsigned int d = 1; d < ImageDimension; d++ )
          {
          searchRadius.push_back( searchRadius[0] );
          }
        }
      if( searchRadius.size() != ImageDimension )
        {
        if( verbose )
          {
          std::cerr << "Search radius specified incorrectly.  Please see usage options." << std::endl;
          }
        return EXIT_FAILURE;
        }
      typename FusionFilterType::NeighborhoodRadiusType searchNeighborhoodRadius;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        searchNeighborhoodRadius[d] = searchRadius[d];
        }
      fusionFilter->SetNeighborhoodSearchRadius( searchNeighborhoodRadius );
      }
    }

  std::vector<unsigned int> patchRadius;
  patchRadius.push_back( 2 );
  typename OptionType::Pointer patchRadiusOption = parser->GetOption( "patch-radius" );
  if( patchRadiusOption && patchRadiusOption->GetNumberOfFunctions() )
    {
    patchRadius = parser->ConvertVector<unsigned int>( patchRadiusOption->GetFunction( 0 )->GetName() );
    }
  if( patchRadius.size() == 1 )
    {
    for( unsigned int d = 1; d < ImageDimension; d++ )
      {
      patchRadius.push_back( patchRadius[0] );
      }
    }
  if( patchRadius.size() != ImageDimension )
    {
    if( verbose )
      {
      std::cerr << "Patch radius specified incorrectly.  Please see usage options." << std::endl;
      }
    return EXIT_FAILURE;
    }
  typename FusionFilterType::NeighborhoodRadiusType patchNeighborhoodRadius;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    patchNeighborhoodRadius[d] = patchRadius[d];
    }

  fusionFilter->SetNeighborhoodPatchRadius( patchNeighborhoodRadius );

  // Check if the user wants to retain atlas voting and/or label posterior images

  bool retainAtlasVotingImages = false;
  bool retainLabelPosteriorImages = false;

  typename OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      retainLabelPosteriorImages = true;
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
      {
      retainAtlasVotingImages = true;
      }
    }

//   typename OptionType::Pointer retainLabelPosteriorOption = parser->GetOption( "retain-label-posterior-images" );
//   if( retainLabelPosteriorOption && retainLabelPosteriorOption->GetNumberOfFunctions() > 0 )
//     {
//     retainLabelPosteriorImages = parser->Convert<bool>( retainLabelPosteriorOption->GetFunction()->GetName() );
//     }
//
//   typename OptionType::Pointer retainAtlasVotingOption = parser->GetOption( "retain-atlas-voting-images" );
//   if( retainAtlasVotingOption && retainAtlasVotingOption->GetNumberOfFunctions() > 0 )
//     {
//     retainAtlasVotingImages = parser->Convert<bool>( retainAtlasVotingOption->GetFunction()->GetName() );
//     }

  bool constrainSolutionToNonnegativeWeights = false;

  typename OptionType::Pointer constrainWeightsOption = parser->GetOption( "constrain-nonnegative" );
  if( constrainWeightsOption && constrainWeightsOption->GetNumberOfFunctions() > 0 )
    {
    constrainSolutionToNonnegativeWeights = parser->Convert<bool>( constrainWeightsOption->GetFunction()->GetName() );
    }

  typename OptionType::Pointer metricOption = parser->GetOption( "patch-metric" );
  if( metricOption && metricOption->GetNumberOfFunctions() > 0 )
    {
    std::string metricString = metricOption->GetFunction()->GetName();
    ConvertToLowerCase( metricString );
    if( metricString.compare( "pc" ) == 0 )
      {
      fusionFilter->SetSimilarityMetric( FusionFilterType::PEARSON_CORRELATION );
      }
    else if( metricString.compare( "msq" ) == 0 )
      {
      fusionFilter->SetSimilarityMetric( FusionFilterType::MEAN_SQUARES );
      }
    else
      {
      std::cerr << "Unrecognized metric option. See help menu." << std::endl;
      return EXIT_FAILURE;
      }
    }

  fusionFilter->SetRetainAtlasVotingWeightImages( retainAtlasVotingImages );
  fusionFilter->SetRetainLabelPosteriorProbabilityImages( retainLabelPosteriorImages );
  fusionFilter->SetConstrainSolutionToNonnegativeWeights( constrainSolutionToNonnegativeWeights );

  // Get the target image

  unsigned int numberOfTargetModalities = 0;

  typename FusionFilterType::InputImageList targetImageList;

  typename OptionType::Pointer targetImageOption = parser->GetOption( "target-image" );
  if( targetImageOption && targetImageOption->GetNumberOfFunctions() )
    {
    if( targetImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      typename ImageType::Pointer targetImage = nullptr;

      std::string targetFile = targetImageOption->GetFunction( 0 )->GetName();
      ReadImage<ImageType>( targetImage, targetFile.c_str() );

      targetImageList.push_back( targetImage );

      numberOfTargetModalities = 1;
      }
    else
      {
      numberOfTargetModalities = targetImageOption->GetFunction( 0 )->GetNumberOfParameters();
      for( unsigned int n = 0; n < numberOfTargetModalities; n++ )
        {
        typename ImageType::Pointer targetImage = nullptr;

        std::string targetFile = targetImageOption->GetFunction( 0 )->GetParameter( n );
        ReadImage<ImageType>( targetImage, targetFile.c_str() );

        targetImageList.push_back( targetImage );
        }
      }
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Target image(s) not specified." << std::endl;
      }
    return EXIT_FAILURE;
    }

  fusionFilter->SetTargetImage( targetImageList );

  // Get the atlas images and segmentations

  typename OptionType::Pointer atlasImageOption = parser->GetOption( "atlas-image" );
  typename OptionType::Pointer atlasSegmentationOption = parser->GetOption( "atlas-segmentation" );

  unsigned int numberOfAtlases = 0;
  unsigned int numberOfAtlasSegmentations = 0;
  unsigned int numberOfAtlasModalities = 0;

  if( atlasImageOption && atlasImageOption->GetNumberOfFunctions() )
    {
    numberOfAtlases = atlasImageOption->GetNumberOfFunctions();
    }
  if( atlasSegmentationOption && atlasSegmentationOption->GetNumberOfFunctions() )
    {
    numberOfAtlasSegmentations = atlasSegmentationOption->GetNumberOfFunctions();
    }

  if( numberOfAtlases < 2 )
    {
    if( verbose )
      {
      std::cerr << "At least 2 atlases are required." << std::endl;
      }
    return EXIT_FAILURE;
    }
  if( numberOfAtlasSegmentations != 0 && numberOfAtlasSegmentations != numberOfAtlases )
    {
    if( verbose )
      {
      std::cout << "Warning:  the number of atlases does not match the number of "
        << "segmentations.  Only performing joint intensity fusion."  << std::endl;
      }
    numberOfAtlasSegmentations = 0;
    }

  for( unsigned int m = 0; m < numberOfAtlases; m++ )
    {
    typename FusionFilterType::InputImageList atlasImageList;
    typename LabelImageType::Pointer atlasSegmentation = nullptr;

    if( atlasImageOption->GetFunction( m )->GetNumberOfParameters() == 0 )
      {
      numberOfAtlasModalities = 1;

      if( numberOfTargetModalities != 1 )
        {
        if( verbose )
          {
          std::cerr << "The number of atlas modalities does not match the number of target modalities." << std::endl;
          }
        return EXIT_FAILURE;
        }
      typename ImageType::Pointer atlasImage = nullptr;

      std::string atlasFile = atlasImageOption->GetFunction( m )->GetName();
      ReadImage<ImageType>( atlasImage, atlasFile.c_str() );
      atlasImageList.push_back( atlasImage );
      }
    else
      {
      if( m == 0 )
        {
        numberOfAtlasModalities = atlasImageOption->GetFunction( m )->GetNumberOfParameters();
        }

      if( numberOfAtlasModalities != atlasImageOption->GetFunction( m )->GetNumberOfParameters() )
        {
        if( verbose )
          {
          std::cerr << "The number of atlas modalities does not match the number of target modalities." << std::endl;
          }
        return EXIT_FAILURE;
        }
      for( unsigned int n = 0; n < numberOfAtlasModalities; n++ )
        {
        typename ImageType::Pointer atlasImage = nullptr;

        std::string atlasFile = atlasImageOption->GetFunction( m )->GetParameter( n );
        ReadImage<ImageType>( atlasImage, atlasFile.c_str() );

        atlasImageList.push_back( atlasImage );
        }
      }
    if( numberOfAtlasSegmentations > 0 )
      {
      std::string atlasSegmentationFile = atlasSegmentationOption->GetFunction( m )->GetName();
      ReadImage<LabelImageType>( atlasSegmentation, atlasSegmentationFile.c_str() );
      }
    fusionFilter->AddAtlas( atlasImageList, atlasSegmentation );
    }

  // Get the exclusion images

  typename OptionType::Pointer exclusionImageOption = parser->GetOption( "exclusion-image" );
  if( exclusionImageOption && exclusionImageOption->GetNumberOfFunctions() )
    {
    for( unsigned int n = 0; n < exclusionImageOption->GetNumberOfFunctions(); n++ )
      {
      auto label = parser->Convert<LabelType>( exclusionImageOption->GetFunction( n )->GetName() );

      typename LabelImageType::Pointer exclusionImage = nullptr;
      std::string exclusionFile = exclusionImageOption->GetFunction( n )->GetParameter( 0 );
      ReadImage<LabelImageType>( exclusionImage, exclusionFile.c_str() );
      fusionFilter->AddLabelExclusionImage( label, exclusionImage );
      }
    }

  // Get the mask

  typename itk::ants::CommandLineParser::OptionType::Pointer maskImageOption =
    parser->GetOption( "mask-image" );
  if( maskImageOption && maskImageOption->GetNumberOfFunctions() )
    {
    typename MaskImageType::Pointer maskImage = nullptr;

    std::string inputFile = maskImageOption->GetFunction( 0 )->GetName();
    ReadImage<MaskImageType>( maskImage, inputFile.c_str() );

    fusionFilter->SetMaskImage( maskImage );
    }

  // Run the fusion program

  itk::TimeProbe timer;
  timer.Start();

  if( verbose )
    {
    typedef CommandProgressUpdate<FusionFilterType> CommandType;
    typename CommandType::Pointer observer = CommandType::New();
    fusionFilter->AddObserver( itk::ProgressEvent(), observer );
    }

  try
    {
    fusionFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    if( verbose )
      {
      std::cerr << "Exception caught: " << e << std::endl;
      }
    return EXIT_FAILURE;
    }

  timer.Stop();

  if( verbose )
    {
    std::cout << std::endl << std::endl;
    fusionFilter->Print( std::cout, 3 );
    }

  // write the output

  if( verbose )
    {
    std::cout << std::endl << "Writing output:" << std::endl;
    }
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string labelFusionName;
    std::string intensityFusionName;
    std::string labelPosteriorName;
    std::string atlasVotingName;

    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      if( numberOfAtlasSegmentations != 0 )
        {
        labelFusionName = outputOption->GetFunction( 0 )->GetName();
        }
      else
        {
        intensityFusionName = outputOption->GetFunction( 0 )->GetName();
        }
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      if( numberOfAtlasSegmentations != 0 )
        {
        labelFusionName = outputOption->GetFunction( 0 )->GetParameter( 0 );
        }
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      intensityFusionName = outputOption->GetFunction( 0 )->GetParameter( 1 );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      if( numberOfAtlasSegmentations != 0 )
        {
        labelPosteriorName = outputOption->GetFunction( 0 )->GetParameter( 2 );
        }
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
      {
      atlasVotingName = outputOption->GetFunction( 0 )->GetParameter( 3 );
      }

    if( !labelFusionName.empty() )
      {
      WriteImage<LabelImageType>( fusionFilter->GetOutput(), labelFusionName.c_str() );
      }
    if( !intensityFusionName.empty() )
      {
      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( numberOfAtlasModalities );
      fileNamesCreator->SetSeriesFormat( intensityFusionName.c_str() );

      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < imageNames.size(); i++ )
        {
        if( verbose )
          {
          std::cout << "  Writing intensity fusion image (modality " << i + 1 << ")" << std::endl;
          }
        typename ImageType::Pointer jointIntensityFusionImage
          = fusionFilter->GetJointIntensityFusionImage( i );
        WriteImage<ImageType>( jointIntensityFusionImage, imageNames[i].c_str() );
        }
      }
    if( !labelPosteriorName.empty() && fusionFilter->GetRetainLabelPosteriorProbabilityImages() )
      {
      typename FusionFilterType::LabelSetType labelSet = fusionFilter->GetLabelSet();

      typename FusionFilterType::LabelSetType::const_iterator labelIt;
      for( labelIt = labelSet.begin(); labelIt != labelSet.end(); ++labelIt )
        {
        if( *labelIt == 0 )
          {
          continue;
          }
        if( verbose )
          {
          std::cout << "  Writing label probability image (label " << *labelIt << ")" << std::endl;
          }

        char buffer[256];
        std::snprintf( buffer, sizeof( buffer ), labelPosteriorName.c_str(), *labelIt );
        WriteImage<typename FusionFilterType::ProbabilityImageType>( fusionFilter->GetLabelPosteriorProbabilityImage( *labelIt ), buffer );
        }
      }
    if( !atlasVotingName.empty() && fusionFilter->GetRetainAtlasVotingWeightImages() )
      {
      itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
      fileNamesCreator->SetStartIndex( 1 );
      fileNamesCreator->SetEndIndex( numberOfAtlases );
      fileNamesCreator->SetSeriesFormat( atlasVotingName.c_str() );

      const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();
      for( unsigned int i = 0; i < imageNames.size(); i++ )
        {
        if( verbose )
          {
          std::cout << "  Writing atlas voting image (atlas " << i+1 << ")" << std::endl;
          }
        WriteImage<typename FusionFilterType::ProbabilityImageType>( fusionFilter->GetAtlasVotingWeightImage( i ), imageNames[i].c_str() );
        }
      }
    }

  if( verbose )
    {
    std::cout << "Elapsed time: " << timer.GetMean() << std::endl;
    }

  return EXIT_SUCCESS;
}

void ajfInitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, the program tries to " )
      + std::string( "infer the dimensionality from the input image." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3/4" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The target image (or multimodal target images) assumed to be " )
    + std::string( "aligned to a common image domain." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "target-image" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "targetImage" );
  option->SetUsageOption( 1, "[targetImageModality0,targetImageModality1,...,targetImageModalityN]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The atlas image (or multimodal atlas images) assumed to be " )
    + std::string( "aligned to a common image domain." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "atlas-image" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "atlasImage" );
  option->SetUsageOption( 1, "[atlasImageModality0,atlasImageModality1,...,atlasImageModalityN]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The atlas segmentation images.  For performing label fusion the number of " )
    + std::string( "specified segmentations should be identical to the number of atlas image sets." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "atlas-segmentation" );
  option->SetShortName( 'l' );
  option->SetUsageOption( 0, "atlasSegmentation" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Regularization term added to matrix Mx for calculating the inverse.  Default = 0.1" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "alpha" );
  option->SetShortName( 'a' );
  option->SetUsageOption( 0, "0.1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Exponent for mapping intensity difference to the joint error.  Default = 2.0" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "beta" );
  option->SetShortName( 'b' );
  option->SetUsageOption( 0, "2.0" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

//   {
//   std::string description =
//     std::string( "Retain label posterior probability images.  Requires atlas segmentations " )
//     + std::string( "to be specified.  Default = false" );
//
//   OptionType::Pointer option = OptionType::New();
//   option->SetLongName( "retain-label-posterior-images" );
//   option->SetShortName( 'r' );
//   option->SetUsageOption( 0, "(0)/1" );
//   option->SetDescription( description );
//   parser->AddOption( option );
//   }
//
//   {
//   std::string description =
//     std::string( "Retain atlas voting images.  Default = false" );
//
//   OptionType::Pointer option = OptionType::New();
//   option->SetLongName( "retain-atlas-voting-images" );
//   option->SetShortName( 'f' );
//   option->SetUsageOption( 0, "(0)/1" );
//   option->SetDescription( description );
//   parser->AddOption( option );
//   }

  {
  std::string description = std::string( "Constrain solution to non-negative weights." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'c' );
  option->SetLongName( "constrain-nonnegative" );
  option->SetUsageOption( 0, "(0)/1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Patch radius for similarity measures.  Default = 2x2x2" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "patch-radius" );
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "2" );
  option->SetUsageOption( 1, "2x2x2" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Metric to be used in determining the most similar neighborhood patch.  " )
    + std::string( "Options include Pearson's correlation (PC) and mean squares (MSQ). " )
    + std::string( "Default = PC (Pearson correlation)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "patch-metric" );
  option->SetShortName( 'm' );
  option->SetUsageOption( 0, "(PC)/MSQ" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Search radius for similarity measures.  Default = 3x3x3.  One " )
    + std::string( "can also specify an image where the value at the voxel specifies " )
    + std::string( "the isotropic search radius at that voxel." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "search-radius" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "3" );
  option->SetUsageOption( 1, "3x3x3" );
  option->SetUsageOption( 2, "searchRadiusMap.nii.gz" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Specify an exclusion region for the given label." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "exclusion-image" );
  option->SetShortName( 'e' );
  option->SetUsageOption( 0, "label[exclusionImage]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "If a mask image is specified, fusion is only performed in the mask region." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "mask-image" );
  option->SetShortName( 'x' );
  option->SetUsageOption( 0, "maskImageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The output is the intensity and/or label fusion image.  Additional " )
    + std::string( "optional outputs include the label posterior probability images " )
    + std::string( "and the atlas voting weight images." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "labelFusionImage" );
  option->SetUsageOption( 1, "intensityFusionImageFileNameFormat" );
  option->SetUsageOption( 2, "[labelFusionImage,intensityFusionImageFileNameFormat,<labelPosteriorProbabilityImageFileNameFormat>,<atlasVotingWeightImageFileNameFormat>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }


  {
  std::string description = std::string( "Get version information." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "version" );
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
int antsJointFusion( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsJointFusion" );

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
    std::string( "antsJointFusion is an image fusion algorithm developed by Hongzhi Wang and " )
    + std::string( "Paul Yushkevich which won segmentation challenges at MICCAI 2012 and MICCAI 2013. " )
    + std::string( "The original label fusion framework was extended to accommodate intensities by " )
    + std::string( "Brian Avants.  This implementation is based on Paul's original ITK-style " )
    + std::string( "implementation and Brian's ANTsR implementation.  References include  1) H. Wang, " )
    + std::string( "J. W. Suh, S. Das, J. Pluta, C. Craige, P. Yushkevich, Multi-atlas " )
    + std::string( "segmentation with joint label fusion IEEE Trans. on Pattern " )
    + std::string( "Analysis and Machine Intelligence, 35(3), 611-623, 2013. and 2) " )
    + std::string( "H. Wang and P. A. Yushkevich, Multi-atlas segmentation with joint " )
    + std::string( "label fusion and corrective learning--an open source implementation, " )
    + std::string( "Front. Neuroinform., 2013. " );

  parser->SetCommandDescription( commandDescription );
  ajfInitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc == 1 )
    {
    parser->PrintMenu( std::cerr, 5, false );
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
  // Show automatic version
  itk::ants::CommandLineParser::OptionType::Pointer versionOption = parser->GetOption( "version" );
  if( versionOption && versionOption->GetNumberOfFunctions() )
    {
    std::string versionFunction = versionOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( versionFunction );
    if( versionFunction.compare( "1" ) == 0 || versionFunction.compare( "true" ) == 0 )
      {
      //Print Version Information
      std::cout << ANTs::Version::ExtendedVersionString() << std::endl;
      return EXIT_SUCCESS;
      }
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
      parser->GetOption( "target-image" );
    if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
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
               << " with the -t option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  switch( dimension )
    {
    case 2:
      {
      return antsJointFusion<2>( parser );
      }
      break;
    case 3:
      {
      return antsJointFusion<3>( parser );
      }
      break;
    case 4:
      {
      return antsJointFusion<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
