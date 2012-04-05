
#include "antscout.hxx"
#include <algorithm>

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

    antscout << "Iteration " << filter->GetElapsedIterations()
             << " (of " << filter->GetMaximumNumberOfIterations() << ").  ";
    antscout << " Current convergence value = "
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

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  if( argc < 3 )
    {
    antscout << "missing 1st filename" << std::endl;
    throw;
    }
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( reader->GetOutput() );
  shrinker->SetShrinkFactors( 1 );

  typename MaskImageType::Pointer maskImage = NULL;

  if( argc > 5 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( argv[5] );

    try
      {
      maskreader->Update();
      maskImage = maskreader->GetOutput();
      }
    catch( ... )
      {
      antscout << "Mask file not read.  Generating mask file using otsu"
               << " thresholding." << std::endl;
      }
    }
  if( !maskImage )
    {
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType>
      ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput( reader->GetOutput() );
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

  if( argc > 4 )
    {
    shrinker->SetShrinkFactors( atoi( argv[4] ) );
    maskshrinker->SetShrinkFactors( atoi( argv[4] ) );
    }
  shrinker->Update();
  maskshrinker->Update();

  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, MaskImageType,
                                                   ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  if( argc > 6 )
    {
    correcter->SetMaximumNumberOfIterations( atoi( argv[6] ) );
    }
  if( argc > 7 )
    {
    correcter->SetNumberOfFittingLevels( atoi( argv[7] ) );
    }

  typedef CommandIterationUpdate<CorrecterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  correcter->AddObserver( itk::IterationEvent(), observer );

  try
    {
    correcter->Update();
    }
  catch( ... )
    {
    antscout << "Exception caught." << std::endl;
    return EXIT_FAILURE;
    }

//  correcter->Print( antscout, 3 );

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
    reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( reader->GetOutput()->GetOrigin() );
  bspliner->SetDirection( reader->GetOutput()->GetDirection() );
  bspliner->SetSpacing( reader->GetOutput()->GetSpacing() );
  bspliner->Update();

  typename ImageType::Pointer logField = ImageType::New();
  logField->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  logField->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  logField->SetRegions(
    bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
  logField->SetDirection( bspliner->GetOutput()->GetDirection() );
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
  divider->SetInput1( reader->GetOutput() );
  divider->SetInput2( expFilter->GetOutput() );
  divider->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  if( argc < 4 )
    {
    antscout << "missing divider image filename" << std::endl;
    throw;
    }
  writer->SetFileName( argv[3] );
  writer->SetInput( divider->GetOutput() );
  writer->Update();

  if( argc > 8 )
    {
    writer = WriterType::New();
    writer->SetFileName( argv[8] );
    writer->SetInput( expFilter->GetOutput() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int N3BiasFieldCorrection( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "N3BiasFieldCorrection" );

  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  antscout->set_stream( out_stream );

  if( argc < 4 )
    {
    antscout << "Usage: " << argv[0] << " imageDimension inputImage "
             << "outputImage [shrinkFactor] [maskImage] [numberOfIterations] "
             << "[numberOfFittingLevels] [outputBiasField] " << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      {
      N3BiasFieldCorrection<2>( argc, argv );
      }
      break;
    case 3:
      {
      N3BiasFieldCorrection<3>( argc, argv );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
}
} // namespace ants
