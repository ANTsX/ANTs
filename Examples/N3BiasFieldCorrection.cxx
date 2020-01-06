
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"

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

    std::cout << "Iteration " << filter->GetElapsedIterations()
             << " (of " << filter->GetMaximumNumberOfIterations() << ").  ";
    std::cout << " Current convergence value = "
             << filter->GetCurrentConvergenceMeasurement()
             << " (threshold = " << filter->GetConvergenceThreshold()
             << ")" << std::endl;
  }
};

template <unsigned int ImageDimension>
int N3BiasFieldCorrection( int argc, char *argv[] )
{
  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension>      ImageType;
  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  typename ImageType::Pointer image;
  ReadImage<ImageType>( image, argv[2] );

  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( image );
  shrinker->SetShrinkFactors( 1 );

  typename MaskImageType::Pointer maskImage = nullptr;

  if( argc > 5 )
    {
    ReadImage<MaskImageType>( maskImage, argv[5] );
    }
  if( !maskImage )
    {
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType>
      ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput( image );
    otsu->SetNumberOfHistogramBins( 200 );
    otsu->SetInsideValue( 0 );
    otsu->SetOutsideValue( 1 );
    otsu->Update();

    maskImage = otsu->GetOutput();
    }
  typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType> MaskShrinkerType;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput( maskImage );
  maskshrinker->SetShrinkFactors( 1 );

  // 4 is [shrinkFactor]
  // 5 is [maskImage]
  if( argc > 4 )
    {
    shrinker->SetShrinkFactors( std::stoi( argv[4] ) );
    maskshrinker->SetShrinkFactors( std::stoi( argv[4] ) );
    }
  shrinker->Update();
  maskshrinker->Update();

  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, MaskImageType,
                                                   ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  // 6 is [numberOfIterations]
  if( argc > 6 )
    {
    correcter->SetMaximumNumberOfIterations( std::stoi( argv[6] ) );
    }
  // 7 is [numberOfFittingLevels]
  if( argc > 7 )
    {
    correcter->SetNumberOfFittingLevels( std::stoi( argv[7] ) );
    }

  // 8 is [outputBiasField]
  // 9 is [verbose]
  bool verbose = false;
  if( argc > 9 )
    {
    verbose = static_cast<bool>(std::stoi( argv[9] ));
    }

  typedef CommandIterationUpdate<CorrecterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  if ( verbose )
    {
    correcter->AddObserver( itk::IterationEvent(), observer );
    }
  try
    {
    correcter->Update();
    }
  catch( ... )
    {
    std::cout << "Exception caught." << std::endl;
    return EXIT_FAILURE;
    }

//  correcter->Print( std::cout, 3 );

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
  bspliner->SetSize(
    image->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( image->GetOrigin() );
  bspliner->SetDirection( image->GetDirection() );
  bspliner->SetSpacing( image->GetSpacing() );
  bspliner->Update();

  typename ImageType::Pointer logField =
    AllocImage<ImageType>(bspliner->GetOutput() );

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
  divider->SetInput1( image );
  divider->SetInput2( expFilter->GetOutput() );
  divider->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  if( argc < 4 )
    {
    std::cout << "missing divider image filename" << std::endl;
    throw;
    }
  WriteImage<ImageType>( divider->GetOutput(), argv[3] );

  if( argc > 8 )
    {
    WriteImage<ImageType>( expFilter->GetOutput(), argv[8] );
    }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int N3BiasFieldCorrection( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "N3BiasFieldCorrection" );

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

  if( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
             << "outputImage [shrinkFactor] [maskImage] [numberOfIterations] "
             << "[numberOfFittingLevels] [outputBiasField] [verbose]" << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  switch( std::stoi( argv[1] ) )
    {
    case 2:
      {
      return N3BiasFieldCorrection<2>( argc, argv );
      }
      break;
    case 3:
      {
      return N3BiasFieldCorrection<3>( argc, argv );
      }
      break;
    case 4:
      {
      return N3BiasFieldCorrection<4>( argc, argv );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
