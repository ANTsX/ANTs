#include "itkBSplineControlPointImageFilter.h"
#include "antsCommandLineParser.h"
#include "itkConstantPadImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "itkTimeProbe.h"

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
    if( filter->GetElapsedIterations() == 1 )
      {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1
                << std::endl;
      }
    std::cout << "  Iteration " << filter->GetElapsedIterations()
              << " (of "
              << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()]
              << ").  ";
    std::cout << " Current convergence value = "
              << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
  }
};

template <unsigned int ImageDimension>
int N4( itk::ants::CommandLineParser *parser )
{
  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typename ImageType::Pointer inputImage = NULL;

  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage = NULL;

  typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, MaskImageType,
                                                ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  typename itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "input-image" );
  if( inputImageOption )
    {
    std::string inputFile = inputImageOption->GetValue();
    reader->SetFileName( inputFile.c_str() );

    inputImage = reader->GetOutput();
    inputImage->Update();
    inputImage->DisconnectPipeline();
    }
  else
    {
    std::cerr << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * handle the mask image
   */

  typename itk::ants::CommandLineParser::OptionType::Pointer maskImageOption =
    parser->GetOption( "mask-image" );
  if( maskImageOption && maskImageOption->GetNumberOfValues() )
    {
    std::string inputFile = maskImageOption->GetValue();
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( inputFile.c_str() );
    try
      {
      maskImage = maskreader->GetOutput();
      maskImage->Update();
      maskImage->DisconnectPipeline();
      }
    catch( ... )
      {
      }
    }
  if( !maskImage )
    {
    std::cout << "Mask not read.  Creating Otsu mask." << std::endl;
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType>
      ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput( inputImage );
    otsu->SetNumberOfHistogramBins( 200 );
    otsu->SetInsideValue( 0 );
    otsu->SetOutsideValue( 1 );

    maskImage = otsu->GetOutput();
    maskImage->Update();
    maskImage->DisconnectPipeline();
    }

  typename ImageType::Pointer weightImage = NULL;

  typename itk::ants::CommandLineParser::OptionType::Pointer weightImageOption =
    parser->GetOption( "weight-image" );
  if( weightImageOption  && weightImageOption->GetNumberOfValues() )
    {
    std::string inputFile = weightImageOption->GetValue();
    typedef itk::ImageFileReader<ImageType> WeightReaderType;
    typename WeightReaderType::Pointer weightreader = WeightReaderType::New();
    weightreader->SetFileName( inputFile.c_str() );
    weightImage = weightreader->GetOutput();
    weightImage->Update();
    weightImage->DisconnectPipeline();
    }

  /**
   * convergence options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence" );
  if( convergenceOption && convergenceOption->GetNumberOfValues() )
    {
    if( convergenceOption->GetNumberOfParameters() > 0 )
      {
      std::vector<unsigned int> numIters = parser->ConvertVector<unsigned int>(
          convergenceOption->GetParameter( 0 ) );
      typename CorrecterType::VariableSizeArrayType
        maximumNumberOfIterations( numIters.size() );
      for( unsigned int d = 0; d < numIters.size(); d++ )
        {
        maximumNumberOfIterations[d] = numIters[d];
        }
      correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

      typename CorrecterType::ArrayType numberOfFittingLevels;
      numberOfFittingLevels.Fill( numIters.size() );
      correcter->SetNumberOfFittingLevels( numberOfFittingLevels );
      }
    if( convergenceOption->GetNumberOfParameters() > 1 )
      {
      correcter->SetConvergenceThreshold( parser->Convert<float>(
                                            convergenceOption->GetParameter( 1 ) ) );
      }
    }
  else // set default values
    {
    typename CorrecterType::VariableSizeArrayType
      maximumNumberOfIterations( 4 );
    maximumNumberOfIterations.Fill( 50 );
    correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
    correcter->SetNumberOfFittingLevels( 4 );
    correcter->SetConvergenceThreshold( 0.000001 );
    }

  /**
   * B-spline options -- we place this here to take care of the case where
   * the user wants to specify things in terms of the spline distance.
   */

  typename ImageType::IndexType inputImageIndex =
    inputImage->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType inputImageSize =
    inputImage->GetLargestPossibleRegion().GetSize();

  typename ImageType::PointType newOrigin = inputImage->GetOrigin();

  typename itk::ants::CommandLineParser::OptionType::Pointer bsplineOption =
    parser->GetOption( "bspline-fitting" );
  if( bsplineOption )
    {
    if( bsplineOption->GetNumberOfParameters() > 1 )
      {
      correcter->SetSplineOrder( parser->Convert<unsigned int>(
                                   bsplineOption->GetParameter( 1 ) ) );
      }
    if( bsplineOption->GetNumberOfParameters() > 0 )
      {
      std::vector<float> array = parser->ConvertVector<float>(
          bsplineOption->GetParameter( 0 ) );
      typename CorrecterType::ArrayType numberOfControlPoints;
      if( array.size() == 1 )
        {
        // the user wants to specify things in terms of spline distance.
        //  1. need to pad the images to get as close to possible to the
        //     requested domain size.
        float splineDistance = array[0];

        unsigned long lowerBound[ImageDimension];
        unsigned long upperBound[ImageDimension];
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          float domain = static_cast<RealType>( inputImage->
                                                GetLargestPossibleRegion().GetSize()[d]
                                                - 1 ) * inputImage->GetSpacing()[d];
          unsigned int numberOfSpans = static_cast<unsigned int>(
              vcl_ceil( domain / splineDistance ) );
          unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans
                                                                     * splineDistance
                                                                     - domain ) / inputImage->GetSpacing()[d] + 0.5 );
          lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
          upperBound[d] = extraPadding - lowerBound[d];
          newOrigin[d] -= ( static_cast<RealType>( lowerBound[d] )
                            * inputImage->GetSpacing()[d] );
          numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
          }

        typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
        typename PadderType::Pointer padder = PadderType::New();
        padder->SetInput( inputImage );
        padder->SetPadLowerBound( lowerBound );
        padder->SetPadUpperBound( upperBound );
        padder->SetConstant( 0 );
        padder->Update();

        inputImage = padder->GetOutput();
        inputImage->DisconnectPipeline();

        typedef itk::ConstantPadImageFilter<MaskImageType, MaskImageType> MaskPadderType;
        typename MaskPadderType::Pointer maskPadder = MaskPadderType::New();
        maskPadder->SetInput( maskImage );
        maskPadder->SetPadLowerBound( lowerBound );
        maskPadder->SetPadUpperBound( upperBound );
        maskPadder->SetConstant( 0 );
        maskPadder->Update();

        maskImage = maskPadder->GetOutput();
        maskImage->DisconnectPipeline();

        if( weightImage )
          {
          typename PadderType::Pointer weightPadder = PadderType::New();
          weightPadder->SetInput( weightImage );
          weightPadder->SetPadLowerBound( lowerBound );
          weightPadder->SetPadUpperBound( upperBound );
          weightPadder->SetConstant( 0 );
          weightPadder->Update();

          weightImage = weightPadder->GetOutput();
          weightImage->DisconnectPipeline();
          }
        }
      else if( array.size() == ImageDimension )
        {
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          numberOfControlPoints[d] = static_cast<unsigned int>( array[d] )
            + correcter->GetSplineOrder();
          }
        }
      else
        {
        std::cerr << "Incorrect mesh resolution" << std::endl;
        return EXIT_FAILURE;
        }
      correcter->SetNumberOfControlPoints( numberOfControlPoints );
      }
    }

  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( inputImage );
  shrinker->SetShrinkFactors( 1 );

  typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType> MaskShrinkerType;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput( maskImage );
  maskshrinker->SetShrinkFactors( 1 );

  typename itk::ants::CommandLineParser::OptionType::Pointer shrinkFactorOption =
    parser->GetOption( "shrink-factor" );
  int shrinkFactor = 4;
  if( shrinkFactorOption && shrinkFactorOption->GetNumberOfValues() )
    {
    shrinkFactor = parser->Convert<int>( shrinkFactorOption->GetValue() );
    }
  shrinker->SetShrinkFactors( shrinkFactor );
  maskshrinker->SetShrinkFactors( shrinkFactor );
  shrinker->Update();
  maskshrinker->Update();

  itk::TimeProbe timer;
  timer.Start();

  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  typedef itk::ShrinkImageFilter<ImageType, ImageType> WeightShrinkerType;
  typename WeightShrinkerType::Pointer weightshrinker = WeightShrinkerType::New();
  if( weightImage )
    {
    weightshrinker->SetInput( weightImage );
    weightshrinker->SetShrinkFactors(shrinkFactor);
    weightshrinker->Update();

    correcter->SetConfidenceImage( weightshrinker->GetOutput() );
    }

  typedef CommandIterationUpdate<CorrecterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  correcter->AddObserver( itk::IterationEvent(), observer );

  /**
   * histogram sharpening options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer histOption =
    parser->GetOption( "histogram-sharpening" );
  if( histOption )
    {
    if( histOption->GetNumberOfParameters() > 0 )
      {
      correcter->SetBiasFieldFullWidthAtHalfMaximum( parser->Convert<float>(
                                                       histOption->GetParameter( 0 ) ) );
      }
    if( histOption->GetNumberOfParameters() > 1 )
      {
      correcter->SetWienerFilterNoise( parser->Convert<float>(
                                         histOption->GetParameter( 1 ) ) );
      }
    if( histOption->GetNumberOfParameters() > 2 )
      {
      correcter->SetNumberOfHistogramBins( parser->Convert<unsigned int>(
                                             histOption->GetParameter( 2 ) ) );
      }
    }

  try
    {
    // correcter->DebugOn();
    correcter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }

  correcter->Print( std::cout, 3 );

  timer.Stop();
  std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;

  /**
   * output
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption )
    {
    /**
                    * Reconstruct the bias field at full image resolution.  Divide
                    * the original input image by the bias field to get the final
                    * corrected image.
                    */
    typedef itk::BSplineControlPointImageFilter<typename
                                                CorrecterType::BiasFieldControlPointLatticeType, typename
                                                CorrecterType::ScalarImageType> BSplinerType;
    typename BSplinerType::Pointer bspliner = BSplinerType::New();
    bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
    bspliner->SetSplineOrder( correcter->GetSplineOrder() );
    bspliner->SetSize( inputImage->GetLargestPossibleRegion().GetSize() );
    bspliner->SetOrigin( newOrigin );
    bspliner->SetDirection( inputImage->GetDirection() );
    bspliner->SetSpacing( inputImage->GetSpacing() );
    bspliner->Update();

    typename ImageType::Pointer logField = ImageType::New();
    logField->SetOrigin( inputImage->GetOrigin() );
    logField->SetSpacing( inputImage->GetSpacing() );
    logField->SetRegions( inputImage->GetLargestPossibleRegion() );
    logField->SetDirection( inputImage->GetDirection() );
    logField->Allocate();

    itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
      bspliner->GetOutput(),
      bspliner->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> ItF( logField,
                                             logField->GetLargestPossibleRegion() );
    for( ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF )
      {
      ItF.Set( ItB.Get()[0] );
      }

    typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput( logField );
    expFilter->Update();

    typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput1( inputImage );
    divider->SetInput2( expFilter->GetOutput() );
    divider->Update();

    if( weightImage &&
        ( maskImageOption && maskImageOption->GetNumberOfValues() > 0 ) )
      {
      itk::ImageRegionIteratorWithIndex<ImageType> ItD( divider->GetOutput(),
                                                        divider->GetOutput()->GetLargestPossibleRegion() );
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        if( maskImage->GetPixel( ItD.GetIndex() ) != correcter->GetMaskLabel() )
          {
          ItD.Set( inputImage->GetPixel( ItD.GetIndex() ) );
          }
        }
      }

    typename ImageType::RegionType inputRegion;
    inputRegion.SetIndex( inputImageIndex );
    inputRegion.SetSize( inputImageSize );

    typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
    typename CropperType::Pointer cropper = CropperType::New();
    cropper->SetInput( divider->GetOutput() );
    cropper->SetExtractionRegion( inputRegion );
    cropper->SetDirectionCollapseToSubmatrix();
    cropper->Update();

    typename CropperType::Pointer biasFieldCropper = CropperType::New();
    biasFieldCropper->SetInput( expFilter->GetOutput() );
    biasFieldCropper->SetExtractionRegion( inputRegion );
    biasFieldCropper->SetDirectionCollapseToSubmatrix();
    biasFieldCropper->Update();

    if( outputOption->GetNumberOfParameters() == 0 )
      {
      typedef  itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( cropper->GetOutput() );
      writer->SetFileName( ( outputOption->GetValue() ).c_str() );
      writer->Update();
      }
    if( outputOption->GetNumberOfParameters() > 0 )
      {
      typedef  itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( cropper->GetOutput() );
      writer->SetFileName( ( outputOption->GetParameter( 0 ) ).c_str() );
      writer->Update();
      }
    if( outputOption->GetNumberOfParameters() > 1 )
      {
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( ( outputOption->GetParameter( 1 ) ).c_str() );
      writer->SetInput( biasFieldCropper->GetOutput() );
      writer->Update();
      }
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, N4 tries to " )
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
      std::string( "A scalar image is expected as input for bias correction.  " )
      + std::string( "Since N4 log transforms the intensities, negative values " )
      + std::string( "or values close to zero should be processed prior to " )
      + std::string( "correction." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input-image" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "inputImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "If a mask image is specified, the final bias correction is " )
      + std::string( "only performed in the mask region.  If a weight image is not " )
      + std::string( "specified, only intensity values inside the masked region are " )
      + std::string( "used during the execution of the algorithm.  If a weight " )
      + std::string( "image is specified, only the non-zero weights are used in the " )
      + std::string( "execution of the algorithm although the mask region defines " )
      + std::string( "where bias correction is performed in the final output. " )
      + std::string( "Otherwise bias correction occurs over the entire image domain. " )
      + std::string( "See also the option description for the weight image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The weight image allows the user to perform a relative " )
      + std::string( "weighting of specific voxels during the B-spline fitting. " )
      + std::string( "For example, some studies have shown that N3 performed on " )
      + std::string( "white matter segmentations improves performance.  If one " )
      + std::string( "has a spatial probability map of the white matter, one can " )
      + std::string( "use this map to weight the b-spline fitting towards those " )
      + std::string( "voxels which are more probabilistically classified as white " )
      + std::string( "matter.  See also the option description for the mask image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "weight-image" );
    option->SetUsageOption( 0, "weightImageFilename" );
    option->SetShortName( 'w' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Running N4 on large images can be time consuming. " )
      + std::string( "To lessen computation time, the input image can be resampled. " )
      + std::string( "The shrink factor, specified as a single integer, describes " )
      + std::string( "this resampling.  Shrink factors <= 4 are commonly used." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "shrink-factor" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "1/2/3/4/..." );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Convergence is determined by calculating the coefficient of " )
      + std::string( "variation between subsequent iterations. When this value " )
      + std::string( "is less than the specified threshold " )
      + std::string( "from the previous iteration or the maximum number of " )
      + std::string( "iterations is exceeded the program terminates.  Multiple " )
      + std::string( "resolutions can be specified by using 'x' between the number " )
      + std::string( "of iterations at each resolution, e.g. 100x50x50." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "convergence" );
    option->SetShortName( 'c' );
    option->SetUsageOption( 0, "[<numberOfIterations=50x50x50x50>,<convergenceThreshold=0.000001>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "These options describe the b-spline fitting parameters. " )
      + std::string( "The initial b-spline mesh at the coarsest resolution is " )
      + std::string( "specified either as the number of elements in each dimension, " )
      + std::string( "e.g. 2x2x3 for 3-D images, or it can be specified as a " )
      + std::string( "single scalar parameter which describes the isotropic sizing " )
      + std::string( "of the mesh elements. The latter option is typically preferred. " )
      + std::string( "For each subsequent level, the spline distance decreases in " )
      + std::string( "half, or equivalently, the number of mesh elements doubles " )
      + std::string( "Cubic splines (order = 3) are typically used." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bspline-fitting" );
    option->SetShortName( 'b' );
    option->SetUsageOption( 0, "[splineDistance,<splineOrder=3>]" );
    option->SetUsageOption( 1, "[initialMeshResolution,<splineOrder=3>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "These options describe the histogram sharpening parameters, " )
      + std::string( "i.e. the deconvolution step parameters described in the " )
      + std::string( "original N3 algorithm.  The default values have been shown " )
      + std::string( "to work fairly well." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "histogram-sharpening" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "[<FWHM=0.15>,<wienerNoise=0.01>,<numberOfHistogramBins=200>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of the bias corrected version of the " )
      + std::string( "input image.  Optionally, one can also output the estimated " )
      + std::string( "bias field." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "[correctedImage,<biasField>]" );
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
    std::string( "N4 is a variant of the popular N3 (nonparameteric nonuniform " )
    + std::string( "normalization) retrospective bias correction algorithm. Based " )
    + std::string( "on the assumption that the corruption of the low frequency bias " )
    + std::string( "field can be modeled as a convolution of the intensity histogram " )
    + std::string( "by a Gaussian, the basic algorithmic protocol is to iterate " )
    + std::string( "between deconvolving the intensity histogram by a Gaussian, " )
    + std::string( "remapping the intensities, and then spatially smoothing this " )
    + std::string( "result by a B-spline modeling of the bias field itself. "  )
    + std::string( "The modifications from and improvements obtained over " )
    + std::string( "the original N3 algorithm are described in the following paper:  " )
    + std::string( "N. Tustison et al., N4ITK:  Improved N3 Bias Correction, " )
    + std::string( "IEEE Transactions on Medical Imaging, 29(6):1310-1320, June 2010." );

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
      parser->GetOption( "input-image" );
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
                << " with the -i option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  std::cout << std::endl << "Running N4 for "
            << dimension << "-dimensional images." << std::endl << std::endl;

  switch( dimension )
    {
    case 2:
      N4<2>( parser );
      break;
    case 3:
      N4<3>( parser );
      break;
    case 4:
      N4<4>( parser );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
