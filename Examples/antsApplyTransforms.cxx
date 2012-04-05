

#include "antscout.hxx"
#include <algorithm>

#include "antsCommandLineParser.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformToDisplacementFieldSource.h"
#include "itkVector.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

#include <deque>
#include <string>
#include <vector>
#include <algorithm>

namespace ants
{
void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
}

template <class TPixel, unsigned int VDim>
class VectorPixelCompare
{
public:
  bool operator()( const itk::Vector<TPixel, VDim> & v1,
                   const itk::Vector<TPixel, VDim> & v2 )
  {
    // Ordering of vectors based on 1st component, then second, etc.
    for( size_t i = 0; i < VDim; i++ )
      {
      if( v1[i] < v2[i] )
        {
        return true;
        }
      else if( v1[i] > v2[i] )
        {
        return false;
        }
      }
    return false;
  }
};

template <unsigned int Dimension>
int antsApplyTransforms( itk::ants::CommandLineParser *parser )
{
  typedef double                           RealType;
  typedef double                           PixelType;
  typedef itk::Vector<RealType, Dimension> VectorType;

  typedef itk::Image<PixelType, Dimension>  ImageType;
  typedef itk::Image<VectorType, Dimension> DisplacementFieldType;

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  typename ResamplerType::Pointer resampleFilter = ResamplerType::New();

  /**
   * Input object option - for now, we're limiting this to images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption =
    parser->GetOption( "input" );
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( inputOption && inputOption->GetNumberOfValues() > 0 )
    {
    antscout << "Input object: " << inputOption->GetValue() << std::endl;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( inputOption->GetValue() ).c_str() );
    reader->Update();

    resampleFilter->SetInput( reader->GetOutput() );
    }
  else if( outputOption && outputOption->GetNumberOfValues() > 0 )
    {
    if( outputOption->GetNumberOfParameters( 0 ) > 1 &&
        parser->Convert<unsigned int>( outputOption->GetParameter( 0, 1 ) ) == 0 )
      {
      antscout << "An input image is required." << std::endl;
      return EXIT_FAILURE;
      }
    }

  /**
   * Reference image option
   */

  // read in the image as char since we only need the header information.
  typedef itk::Image<char, Dimension> ReferenceImageType;
  typename ReferenceImageType::Pointer referenceImage = ReferenceImageType::New();

  typename itk::ants::CommandLineParser::OptionType::Pointer referenceOption =
    parser->GetOption( "reference-image" );
  if( referenceOption && referenceOption->GetNumberOfValues() > 0 )
    {
    antscout << "Reference image: " << referenceOption->GetValue() << std::endl;

    // read in the image as char since we only need the header information.
    typedef itk::Image<char, Dimension>              ReferenceImageType;
    typedef itk::ImageFileReader<ReferenceImageType> ReferenceReaderType;
    typename ReferenceReaderType::Pointer referenceReader =
      ReferenceReaderType::New();
    referenceReader->SetFileName( ( referenceOption->GetValue() ).c_str() );

    referenceImage = referenceReader->GetOutput();
    referenceImage->Update();
    referenceImage->DisconnectPipeline();
    }
  else
    {
    antscout << "Error:  No reference image specified." << std::endl;
    return EXIT_FAILURE;
    }
  resampleFilter->SetOutputParametersFromImage( referenceImage );

  /**
   * Transform option
   */
  // Register the matrix offset transform base class to the
  // transform factory for compatibility with the current ANTs.
  typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
  typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

  /**
   * Load an identity transform in case no transforms are loaded.
   */
  typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
  typename IdentityTransformType::Pointer identityTransform =
    IdentityTransformType::New();
  identityTransform->SetIdentity();

  typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform =
    CompositeTransformType::New();
  compositeTransform->AddTransform( identityTransform );

  typename itk::ants::CommandLineParser::OptionType::Pointer transformOption =
    parser->GetOption( "transform" );
  if( transformOption && transformOption->GetNumberOfValues() > 0 )
    {
    std::deque<std::string> transformNames;
    std::deque<std::string> transformTypes;
    for( unsigned int n = 0; n < transformOption->GetNumberOfValues(); n++ )
      {
      std::string transformName;
      std::string transformType;

      typedef itk::Transform<double, Dimension, Dimension> TransformType;
      typename TransformType::Pointer transform;

      bool hasTransformBeenRead = false;
      try
        {
        transformName = transformOption->GetValue( n );

        typedef itk::DisplacementFieldTransform<PixelType, Dimension>
          DisplacementFieldTransformType;

        typedef typename DisplacementFieldTransformType::DisplacementFieldType
          DisplacementFieldType;

        typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
        typename DisplacementFieldReaderType::Pointer fieldReader =
          DisplacementFieldReaderType::New();
        fieldReader->SetFileName( transformName.c_str() );
        fieldReader->Update();

        typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
          DisplacementFieldTransformType::New();
        displacementFieldTransform->SetDisplacementField( fieldReader->GetOutput() );
        transform = dynamic_cast<TransformType *>( displacementFieldTransform.GetPointer() );

        hasTransformBeenRead = true;
        }
      catch( ... )
        {
        hasTransformBeenRead = false;
        }

      if( !hasTransformBeenRead )
        {
        try
          {
          typedef itk::TransformFileReader TransformReaderType;
          typename TransformReaderType::Pointer transformReader
            = TransformReaderType::New();

          if( transformOption->GetNumberOfParameters( n ) == 0 )
            {
            transformName = transformOption->GetValue( n );
            transformReader->SetFileName( transformName.c_str() );
            try
              {
              transformReader->Update();
              }
            catch( itk::ExceptionObject & excp )
              {
              antscout << excp << std::endl;
              return EXIT_FAILURE;
              }
            transform = dynamic_cast<TransformType *>(
                ( ( transformReader->GetTransformList() )->front() ).GetPointer() );
            }
          else
            {
            transformName = transformOption->GetParameter( n, 0 );

            transformReader->SetFileName( transformName.c_str() );
            try
              {
              transformReader->Update();
              }
            catch( itk::ExceptionObject & excp )
              {
              antscout << excp << std::endl;
              return EXIT_FAILURE;
              }

            transform = dynamic_cast<TransformType *>(
                ( ( transformReader->GetTransformList() )->front() ).GetPointer() );
            if( ( transformOption->GetNumberOfParameters( n ) > 1 ) &&
                parser->Convert<bool>( transformOption->GetParameter( n, 1 ) ) )
              {
              transform = dynamic_cast<TransformType *>(
                  transform->GetInverseTransform().GetPointer() );
              if( !transform )
                {
                antscout << "Inverse does not exist for " << transformName
                         << std::endl;
                return EXIT_FAILURE;
                }
              transformName = std::string( "inverse of " ) + transformName;
              }
            }
          }
        catch( const itk::ExceptionObject & e )
          {
          antscout << "Transform reader for "
                   << transformName << " caught an ITK exception:\n";
          e.Print( antscout );
          return EXIT_FAILURE;
          }
        catch( const std::exception & e )
          {
          antscout << "Transform reader for "
                   << transformName << " caught an exception:\n";
          antscout << e.what() << std::endl;
          return EXIT_FAILURE;
          }
        catch( ... )
          {
          antscout << "Transform reader for "
                   << transformName << " caught an unknown exception!!!\n";
          return EXIT_FAILURE;
          }
        }
      compositeTransform->AddTransform( transform );

      transformNames.push_back( transformName );
      transformTypes.push_back( transform->GetNameOfClass() );
      }
    antscout << "The composite transform is comprised of the following transforms "
             << "(in order): " << std::endl;
    for( unsigned int n = 0; n < transformNames.size(); n++ )
      {
      antscout << "  " << n + 1 << ". " << transformNames[n] << " (type = "
               << transformTypes[n] << ")" << std::endl;
      }
    }
  resampleFilter->SetTransform( compositeTransform );

  /**
   * Interpolation option
   */
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType>
    LinearInterpolatorType;
  typename LinearInterpolatorType::Pointer linearInterpolator
    = LinearInterpolatorType::New();

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, RealType>
    NearestNeighborInterpolatorType;
  typename NearestNeighborInterpolatorType::Pointer nearestNeighborInterpolator
    = NearestNeighborInterpolatorType::New();

  typedef itk::BSplineInterpolateImageFunction<ImageType, RealType>
    BSplineInterpolatorType;
  typename BSplineInterpolatorType::Pointer bSplineInterpolator
    = BSplineInterpolatorType::New();

  typedef itk::GaussianInterpolateImageFunction<ImageType, RealType>
    GaussianInterpolatorType;
  typename GaussianInterpolatorType::Pointer gaussianInterpolator
    = GaussianInterpolatorType::New();

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3>
    HammingInterpolatorType;
  typename HammingInterpolatorType::Pointer hammingInterpolator =
    HammingInterpolatorType::New();

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::CosineWindowFunction<3> > CosineInterpolatorType;
  typename CosineInterpolatorType::Pointer cosineInterpolator =
    CosineInterpolatorType::New();

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::WelchWindowFunction<3> > WelchInterpolatorType;
  typename WelchInterpolatorType::Pointer welchInterpolator =
    WelchInterpolatorType::New();

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::LanczosWindowFunction<3> > LanczosInterpolatorType;
  typename LanczosInterpolatorType::Pointer lanczosInterpolator =
    LanczosInterpolatorType::New();

  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3,
                                                    itk::Function::BlackmanWindowFunction<3> > BlackmanInterpolatorType;
  typename BlackmanInterpolatorType::Pointer blackmanInterpolator =
    BlackmanInterpolatorType::New();

  const unsigned int NVectorComponents = 1;
  typedef VectorPixelCompare<RealType, NVectorComponents> CompareType;
  typedef typename itk::LabelImageGaussianInterpolateImageFunction<ImageType,
                                                                   RealType, CompareType> MultiLabelInterpolatorType;
  typename MultiLabelInterpolatorType::Pointer multiLabelInterpolator =
    MultiLabelInterpolatorType::New();

  std::string whichInterpolator( "linear" );

  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption =
    parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfValues() > 0 )
    {
    whichInterpolator = interpolationOption->GetValue();
    ConvertToLowerCase( whichInterpolator );

    if( !std::strcmp( whichInterpolator.c_str(), "linear" ) )
      {
      linearInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( linearInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "nearestneighbor" ) )
      {
      nearestNeighborInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( nearestNeighborInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "bspline" ) )
      {
      bSplineInterpolator->SetInputImage( resampleFilter->GetInput() );
      if( interpolationOption->GetNumberOfParameters() > 0 )
        {
        unsigned int bsplineOrder = parser->Convert<unsigned int>(
            interpolationOption->GetParameter( 0, 0 ) );
        bSplineInterpolator->SetSplineOrder( bsplineOrder );
        }
      resampleFilter->SetInterpolator( bSplineInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "gaussian" ) )
      {
      gaussianInterpolator->SetInputImage( resampleFilter->GetInput() );
      double sigma[Dimension];
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        sigma[d] = resampleFilter->GetInput()->GetSpacing()[d];
        }
      double alpha = 1.0;

      if( interpolationOption->GetNumberOfParameters() > 0 )
        {
        std::vector<double> s = parser->ConvertVector<double>(
            interpolationOption->GetParameter( 0 ) );
        if( s.size() == Dimension )
          {
          for( unsigned int d = 0; d < Dimension; d++ )
            {
            sigma[d] = s[d];
            }
          }
        else
          {
          for( unsigned int d = 0; d < Dimension; d++ )
            {
            sigma[d] = s[0];
            }
          }
        }
      if( interpolationOption->GetNumberOfParameters() > 1 )
        {
        alpha = parser->Convert<double>(
            interpolationOption->GetParameter( 1 ) );
        }
      gaussianInterpolator->SetParameters( sigma, alpha );
      resampleFilter->SetInterpolator( gaussianInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "cosinewindowedsinc" ) )
      {
      cosineInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( cosineInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "hammingwindowedsinc" ) )
      {
      hammingInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( hammingInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "lanczoswindowedsinc" ) )
      {
      lanczosInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( lanczosInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "blackmanwindowedsinc" ) )
      {
      blackmanInterpolator->SetInputImage( resampleFilter->GetInput() );
      resampleFilter->SetInterpolator( blackmanInterpolator );
      }
    else if( !std::strcmp( whichInterpolator.c_str(), "multilabel" ) )
      {
      multiLabelInterpolator->SetInputImage( resampleFilter->GetInput() );

      double sigma[Dimension];
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        sigma[d] = resampleFilter->GetInput()->GetSpacing()[d];
        }
      double alpha = 4.0;

      if( interpolationOption->GetNumberOfParameters() > 0 )
        {
        std::vector<double> s = parser->ConvertVector<double>(
            interpolationOption->GetParameter( 0 ) );
        if( s.size() == Dimension )
          {
          for( unsigned int d = 0; d < Dimension; d++ )
            {
            sigma[d] = s[d];
            }
          }
        else
          {
          for( unsigned int d = 0; d < Dimension; d++ )
            {
            sigma[d] = s[0];
            }
          }
        }

      multiLabelInterpolator->SetParameters( sigma, alpha );
      resampleFilter->SetInterpolator( multiLabelInterpolator );
      }
    else
      {
      antscout << "Error:  Unrecognized interpolation option." << std::endl;
      return EXIT_FAILURE;
      }
    }
  antscout << "Interpolation type: "
           << resampleFilter->GetInterpolator()->GetNameOfClass() << std::endl;

  /**
   * Default voxel value
   */
  resampleFilter->SetDefaultPixelValue( 0 );
  typename itk::ants::CommandLineParser::OptionType::Pointer defaultOption =
    parser->GetOption( "default-value" );
  if( defaultOption && defaultOption->GetNumberOfValues() > 0 )
    {
    PixelType defaultValue =
      parser->Convert<PixelType>( defaultOption->GetValue() );
    resampleFilter->SetDefaultPixelValue( defaultValue );
    }
  antscout << "Default pixel value: "
           << resampleFilter->GetDefaultPixelValue() << std::endl;

  /**
   * output
   */
  if( outputOption && outputOption->GetNumberOfValues() > 0 )
    {
    if( outputOption->GetNumberOfParameters( 0 ) > 1 &&
        parser->Convert<unsigned int>( outputOption->GetParameter( 0, 1 ) ) != 0 )
      {
      antscout << "Output composite transform displacement field: " << outputOption->GetParameter( 0, 0 ) << std::endl;

      typedef typename itk::TransformToDisplacementFieldSource<DisplacementFieldType> ConverterType;
      typename ConverterType::Pointer converter = ConverterType::New();
      converter->SetOutputParametersFromImage( referenceImage );
      converter->SetTransform( resampleFilter->GetTransform() );

      typedef  itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;
      typename DisplacementFieldWriterType::Pointer displacementFieldWriter = DisplacementFieldWriterType::New();
      displacementFieldWriter->SetInput( converter->GetOutput() );
      displacementFieldWriter->SetFileName( ( outputOption->GetParameter( 0, 0 ) ).c_str() );
      displacementFieldWriter->Update();
      }
    else
      {
      std::string outputFileName = "";
      if( outputOption->GetNumberOfParameters( 0 ) > 1 &&
          parser->Convert<unsigned int>( outputOption->GetParameter( 0, 1 ) ) == 0 )
        {
        outputFileName = outputOption->GetParameter( 0, 0 );
        }
      else
        {
        outputFileName = outputOption->GetValue();
        }
      antscout << "Output warped image: " << outputFileName << std::endl;

      typedef  itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( resampleFilter->GetOutput() );
      writer->SetFileName( ( outputFileName ).c_str() );
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
      + std::string( "dimensional image.  If not specified, antsWarp tries to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Currently, the only input objects supported are image " )
      + std::string( "objects.  However, the current framework allows for " )
      + std::string( "warping of other objects such as meshes and point sets. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "inputFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "For warping input images, the reference image defines the " )
      + std::string( "spacing, origin, size, and direction of the output warped " )
      + std::string( "image. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "reference-image" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "imageFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "One can either output the warped image or, if the boolean " )
      + std::string( "is set, one can print out the displacement field based on the" )
      + std::string( "composite transform and the reference image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "warpedOutputFileName" );
    option->SetUsageOption( 1, "[compositeDisplacementField,<printOutCompositeWarpFile=0>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Several interpolation options are available in ITK. " )
      + std::string( "These have all been made available." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "interpolation" );
    option->SetShortName( 'n' );
    option->SetUsageOption( 0, "Linear" );
    option->SetUsageOption( 1, "NearestNeighbor" );
    option->SetUsageOption( 2, "MultiLabel[<sigma=imageSpacing>,<alpha=4.0>]" );
    option->SetUsageOption( 3, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]" );
    option->SetUsageOption( 4, "BSpline[<order=3>]" );
    option->SetUsageOption( 5, "CosineWindowedSinc" );
    option->SetUsageOption( 6, "WelchWindowedSinc" );
    option->SetUsageOption( 7, "HammingWindowedSinc" );
    option->SetUsageOption( 8, "LanczosWindowedSinc" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Several transform options are supported including all " )
      + std::string( "those defined in the ITK library in addition to " )
      + std::string( "a deformation field transform.  The ordering of " )
      + std::string( "the transformations follows the ordering specified " )
      + std::string( "on the command line.  An identity transform is pushed " )
      + std::string( "onto the transformation stack. Each new transform " )
      + std::string( "encountered on the command line is also pushed onto " )
      + std::string( "the transformation stack. Then, to warp the input object, " )
      + std::string( "each point comprising the input object is warped first " )
      + std::string( "according to the last transform pushed onto the stack " )
      + std::string( "followed by the second to last transform, etc. until " )
      + std::string( "the last transform encountered which is the identity " )
      + std::string( "transform. " )
      + std::string( "Also, it should be noted that the inverse transform can " )
      + std::string( "be accommodated with the usual caveat that such an inverse " )
      + std::string( "must be defined by the specified transform class " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "transform" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "transformFileName" );
    option->SetUsageOption( 1, "[transformFileName,useInverse]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Default voxel value to be used with input images only. " )
      + std::string( "Specifies the voxel value when the input point maps outside " )
      + std::string( "the output domain" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "default-value" );
    option->SetShortName( 'v' );
    option->SetUsageOption( 0, "value" );
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

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsApplyTransforms( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsApplyTransforms" );
  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "antsApplyTransforms, applied to an input image, transforms it " )
    + std::string( "according to a reference image and a transform " )
    + std::string( "(or a set of transforms)." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || ( parser->GetOption( "help" ) &&
                    ( parser->Convert<bool>( parser->GetOption( "help" )->GetValue() ) ) ) )
    {
    parser->PrintMenu( antscout, 5, false );
    return EXIT_FAILURE;
    }
  else if( parser->GetOption( 'h' ) &&
           ( parser->Convert<bool>( parser->GetOption( 'h' )->GetValue() ) ) )
    {
    parser->PrintMenu( antscout, 5, true );
    return EXIT_FAILURE;
    }

  // Read in the first intensity image to get the image dimension.
  std::string filename;

  itk::ants::CommandLineParser::OptionType::Pointer inputOption =
    parser->GetOption( "reference-image" );
  if( inputOption && inputOption->GetNumberOfValues() > 0 )
    {
    if( inputOption->GetNumberOfParameters( 0 ) > 0 )
      {
      filename = inputOption->GetParameter( 0, 0 );
      }
    else
      {
      filename = inputOption->GetValue( 0 );
      }
    }
  else
    {
    antscout << "No reference image was specified." << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int              dimension = 3;
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
      filename.c_str(), itk::ImageIOFactory::ReadMode );
  dimension = imageIO->GetNumberOfDimensions();

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetValue() );
    }

  switch( dimension )
    {
    case 2:
      {
      antsApplyTransforms<2>( parser );
      }
      break;
    case 3:
      {
      antsApplyTransforms<3>( parser );
      }
      break;
    case 4:
      {
      antsApplyTransforms<4>( parser );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
