#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkantsRegistrationHelper.h"
#include "ReadWriteImage.h"
#include "TensorFunctions.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformToDisplacementFieldSource.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
template <typename TensorImageType, typename ImageType>
void
CorrectImageTensorDirection( TensorImageType * movingTensorImage, ImageType * referenceImage )
{
  typedef typename TensorImageType::DirectionType DirectionType;
  typedef typename DirectionType::InternalMatrixType MatrixType;
  MatrixType direction =
    movingTensorImage->GetDirection().GetTranspose() * referenceImage->GetDirection().GetVnlMatrix();

  if( !direction.is_identity( 0.00001 ) )
    {
    itk::ImageRegionIterator<TensorImageType> It( movingTensorImage, movingTensorImage->GetBufferedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typedef typename TensorImageType::PixelType TensorType;
      typedef typename TensorImageType::DirectionType::InternalMatrixType TensorMatrixType;

      TensorType tensor = It.Get();
      TensorMatrixType dt;

      Vector2Matrix<TensorType,TensorMatrixType>(tensor,dt);

      dt = direction * dt * direction.transpose();

      tensor = Matrix2Vector<TensorType,TensorMatrixType>(dt);

      It.Set( tensor );
      }
    }
}

template <typename DisplacementFieldType, typename ImageType>
void
CorrectImageVectorDirection( DisplacementFieldType * movingVectorImage, ImageType * referenceImage )
{
  typedef typename DisplacementFieldType::DirectionType DirectionType;

  typename DirectionType::InternalMatrixType direction =
    movingVectorImage->GetDirection().GetTranspose() * referenceImage->GetDirection().GetVnlMatrix();

  typedef typename DisplacementFieldType::PixelType VectorType;
  typedef typename VectorType::ComponentType        ComponentType;

  const unsigned int dimension = ImageType::ImageDimension;

  if( !direction.is_identity( 0.00001 ) )
    {
    itk::ImageRegionIterator<DisplacementFieldType> It( movingVectorImage, movingVectorImage->GetBufferedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      VectorType vector = It.Get();

      vnl_vector<double> internalVector( dimension );
      for( unsigned int d = 0; d < dimension; d++ )
        {
        internalVector[d] = vector[d];
        }

      internalVector.pre_multiply( direction );;
      for( unsigned int d = 0; d < dimension; d++ )
        {
        vector[d] = internalVector[d];
        }

      It.Set( vector );
      }
    }
}

template <unsigned int NDim>
unsigned int numTensorElements()
{
  return NDim + numTensorElements<NDim-1>();
}

template <>
unsigned int numTensorElements<0>()
{
  return 0;
}

template <class T, unsigned int Dimension>
int antsApplyTransforms( itk::ants::CommandLineParser::Pointer & parser, unsigned int inputImageType = 0 )
{
  typedef T                                RealType;
  typedef T                                PixelType;
  typedef itk::Vector<RealType, Dimension> VectorType;

  // typedef unsigned int                     LabelPixelType;
  // typedef itk::Image<PixelType, Dimension> LabelImageType;

  typedef itk::Image<PixelType, Dimension>     ImageType;
  typedef itk::Image<PixelType, Dimension + 1> TimeSeriesImageType;
  typedef itk::Image<VectorType, Dimension>    DisplacementFieldType;
  typedef ImageType                            ReferenceImageType;

  typedef itk::SymmetricSecondRankTensor<RealType, Dimension> TensorPixelType;
  typedef itk::Image<TensorPixelType, Dimension>              TensorImageType;

  const unsigned int NumberOfTensorElements = numTensorElements<Dimension>();

  typename TimeSeriesImageType::Pointer timeSeriesImage = NULL;
  typename TensorImageType::Pointer tensorImage = NULL;
  typename DisplacementFieldType::Pointer vectorImage = NULL;

  std::vector<typename ImageType::Pointer> inputImages;
  inputImages.clear();

  std::vector<typename ImageType::Pointer> outputImages;
  outputImages.clear();

  /**
   * Input object option - for now, we're limiting this to images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption( "input" );
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );

  if( inputImageType == 3 && inputOption && inputOption->GetNumberOfFunctions() )
    {
    antscout << "Input time-series image: " << inputOption->GetFunction( 0 )->GetName() << std::endl;
    ReadImage<TimeSeriesImageType>( timeSeriesImage, ( inputOption->GetFunction( 0 )->GetName() ).c_str() );
    }
  else if( inputImageType == 2 && inputOption && inputOption->GetNumberOfFunctions() )
    {
    antscout << "Input tensor image: " << inputOption->GetFunction( 0 )->GetName() << std::endl;
    ReadTensorImage<TensorImageType>( tensorImage, ( inputOption->GetFunction( 0 )->GetName() ).c_str(), true );
    }
  else if( inputImageType == 0 && inputOption && inputOption->GetNumberOfFunctions() )
    {
    antscout << "Input scalar image: " << inputOption->GetFunction( 0 )->GetName() << std::endl;
    typename ImageType::Pointer image;
    ReadImage<ImageType>( image, ( inputOption->GetFunction( 0 )->GetName() ).c_str()  );
    inputImages.push_back( image );
    }
  else if( inputImageType == 1 && inputOption && inputOption->GetNumberOfFunctions() )
    {
    antscout << "Input vector image: " << inputOption->GetFunction( 0 )->GetName() << std::endl;

    typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( inputOption->GetFunction( 0 )->GetName() ).c_str() );

    try
      {
      vectorImage = reader->GetOutput();
      vectorImage->Update();
      vectorImage->DisconnectPipeline();
      }
    catch( ... )
      {
      std::cerr << "Unable to read vector image " << reader->GetFileName() << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 &&
        parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) == 0 )
      {
      antscout << "An input image is required." << std::endl;
      return EXIT_FAILURE;
      }
    }

  /**
   * Reference image option
   */
  typedef ImageType ReferenceImageType;
  typename ReferenceImageType::Pointer referenceImage;

  typename itk::ants::CommandLineParser::OptionType::Pointer referenceOption =
    parser->GetOption( "reference-image" );
  if( referenceOption && referenceOption->GetNumberOfFunctions() )
    {
    antscout << "Reference image: " << referenceOption->GetFunction( 0 )->GetName() << std::endl;
    ReadImage<ReferenceImageType>( referenceImage,  ( referenceOption->GetFunction( 0 )->GetName() ).c_str() );
    }
  else
    {
    antscout << "Error:  No reference image specified." << std::endl;
    return EXIT_FAILURE;
    }

  if( inputImageType == 1 )
    {
    CorrectImageVectorDirection<DisplacementFieldType, ReferenceImageType>( vectorImage, referenceImage );
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      typedef itk::VectorIndexSelectionCastImageFilter<DisplacementFieldType, ImageType> SelectorType;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput( vectorImage );
      selector->SetIndex( i );
      selector->Update();

      inputImages.push_back( selector->GetOutput() );
      }
    }
  else if( inputImageType == 2 )
    {
    CorrectImageTensorDirection<TensorImageType, ReferenceImageType>( tensorImage, referenceImage );
    for( unsigned int i = 0; i < NumberOfTensorElements; i++ )
      {
      typedef itk::VectorIndexSelectionCastImageFilter<TensorImageType, ImageType> SelectorType;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput( tensorImage );
      selector->SetIndex( i );
      selector->Update();

      inputImages.push_back( selector->GetOutput() );
      }
    }
  else if( inputImageType == 3 )
    {
    typename TimeSeriesImageType::RegionType extractRegion = timeSeriesImage->GetLargestPossibleRegion();
    unsigned int numberOfTimePoints = extractRegion.GetSize()[Dimension];
    int          startTimeIndex = extractRegion.GetIndex()[Dimension];

    extractRegion.SetSize( Dimension, 0 );
    for( unsigned int i = 0; i < numberOfTimePoints; i++ )
      {
      extractRegion.SetIndex( Dimension, startTimeIndex + i );

      typedef itk::ExtractImageFilter<TimeSeriesImageType, ImageType> ExtracterType;
      typename ExtracterType::Pointer extracter = ExtracterType::New();
      extracter->SetInput( timeSeriesImage );
      extracter->SetExtractionRegion( extractRegion );
      extracter->SetDirectionCollapseToSubmatrix();
      extracter->Update();

      inputImages.push_back( extracter->GetOutput() );
      }
    }

  /**
   * Transform option
   */
  // Register the matrix offset transform base class to the
  // transform factory for compatibility with the current ANTs.
  typedef itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
  typedef itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

  typedef itk::CompositeTransform<RealType, Dimension> CompositeTransformType;
  typename itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption( "transform" );

  bool useStaticCastForR = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer rOption =
    parser->GetOption( "static-cast-for-R" );
  if( rOption && rOption->GetNumberOfFunctions() )
    {
    useStaticCastForR = parser->Convert<bool>(  rOption->GetFunction( 0 )->GetName() );
    }

  std::vector<bool> isDerivedTransform;
  typename CompositeTransformType::Pointer compositeTransform =
    GetCompositeTransformFromParserOption<RealType, Dimension>( parser, transformOption, isDerivedTransform, useStaticCastForR );
  if( compositeTransform.IsNull() )
    {
    return EXIT_FAILURE;
    }

  if( !compositeTransform->GetNumberOfParameters() )
    {
    antscout << "WARNING: No transforms found, using identify transform" << std::endl;
    typename MatrixOffsetTransformType::Pointer idTransform = MatrixOffsetTransformType::New();
    idTransform->SetIdentity();
    compositeTransform->AddTransform( idTransform );
    }

  std::string whichInterpolator( "linear" );
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfFunctions() )
    {
    whichInterpolator = interpolationOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( whichInterpolator );
    }

  const size_t VImageDimension = Dimension;
  typename ImageType::SpacingType
    cache_spacing_for_smoothing_sigmas(itk::NumericTraits<typename ImageType::SpacingType::ValueType>::Zero);
  if( !std::strcmp( whichInterpolator.c_str(), "gaussian" )
      ||   !std::strcmp( whichInterpolator.c_str(), "multilabel" )
      )
    {
    cache_spacing_for_smoothing_sigmas = inputImages[0]->GetSpacing();
    }

#include "make_interpolator_snip.tmpl"

  /**
   * Default voxel value
   */
  PixelType defaultValue = 0;
  typename itk::ants::CommandLineParser::OptionType::Pointer defaultOption =
    parser->GetOption( "default-value" );
  if( defaultOption && defaultOption->GetNumberOfFunctions() )
    {
    defaultValue = parser->Convert<PixelType>( defaultOption->GetFunction( 0 )->GetName() );
    }
  antscout << "Default pixel value: " << defaultValue << std::endl;
  for( unsigned int n = 0; n < inputImages.size(); n++ )
    {
    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
    typename ResamplerType::Pointer resampleFilter = ResamplerType::New();
    resampleFilter->SetInput( inputImages[n] );
    resampleFilter->SetOutputParametersFromImage( referenceImage );
    resampleFilter->SetTransform( compositeTransform );
    resampleFilter->SetDefaultPixelValue( defaultValue );

    interpolator->SetInputImage( inputImages[n] );
    resampleFilter->SetInterpolator( interpolator );
    if( n == 0 )
      {
      antscout << "Interpolation type: " << resampleFilter->GetInterpolator()->GetNameOfClass() << std::endl;
      }
    if( inputImageType == 3 )
      {
      antscout << "  Applying transform(s) to time point " << n << " (out of " << inputImages.size() << ")." << std::endl;
      }
    resampleFilter->Update();
    outputImages.push_back( resampleFilter->GetOutput() );
    }

  /**
   * output
   */
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 &&
        parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) != 0 )
      {
      antscout << "Output composite transform displacement field: "
               << outputOption->GetFunction( 0 )->GetParameter( 0 ) << std::endl;

      typedef typename itk::TransformToDisplacementFieldSource<DisplacementFieldType, RealType> ConverterType;
      typename ConverterType::Pointer converter = ConverterType::New();
      converter->SetOutputParametersFromImage( referenceImage );
      converter->SetTransform( compositeTransform );

      typedef  itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;
      typename DisplacementFieldWriterType::Pointer displacementFieldWriter = DisplacementFieldWriterType::New();
      displacementFieldWriter->SetInput( converter->GetOutput() );
      displacementFieldWriter->SetFileName( ( outputOption->GetFunction( 0 )->GetParameter( 0 ) ).c_str() );
      displacementFieldWriter->Update();
      }
    else
      {
      std::string outputFileName = "";
      if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) == 0 )
        {
        outputFileName = outputOption->GetFunction( 0 )->GetParameter( 0 );
        }
      else
        {
        outputFileName = outputOption->GetFunction( 0 )->GetName();
        }
      antscout << "Output warped image: " << outputFileName << std::endl;

      if( inputImageType == 1 )
        {
        if( outputImages.size() != Dimension )
          {
          antscout << "The number of output images does not match the number of vector components." << std::endl;
          return EXIT_FAILURE;
          }

        VectorType zeroVector( 0.0 );

        typename DisplacementFieldType::Pointer outputVectorImage = DisplacementFieldType::New();
        outputVectorImage->CopyInformation( referenceImage );
        outputVectorImage->SetRegions( referenceImage->GetRequestedRegion() );
        outputVectorImage->Allocate();
        outputVectorImage->FillBuffer( zeroVector );

        itk::ImageRegionIteratorWithIndex<DisplacementFieldType> It( outputVectorImage,
                                                                     outputVectorImage->GetRequestedRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          VectorType vector = It.Get();
          typename DisplacementFieldType::IndexType index = It.GetIndex();
          for( unsigned int n = 0; n < Dimension; n++ )
            {
            vector.SetNthComponent( n, outputImages[n]->GetPixel( index ) );
            }
          It.Set( vector );
          }
        typedef  itk::ImageFileWriter<DisplacementFieldType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( outputVectorImage );
        writer->SetFileName( ( outputFileName ).c_str() );
        writer->Update();
        }
      else if( inputImageType == 2 )
        {
        if( outputImages.size() != NumberOfTensorElements )
          {
          antscout << "The number of output images does not match the number of tensor elements." << std::endl;
          return EXIT_FAILURE;
          }

        TensorPixelType zeroTensor( 0.0 );

        typename TensorImageType::Pointer outputTensorImage = TensorImageType::New();
        outputTensorImage->CopyInformation( referenceImage );
        outputTensorImage->SetRegions( referenceImage->GetRequestedRegion() );
        outputTensorImage->Allocate();
        outputTensorImage->FillBuffer( zeroTensor );

        itk::ImageRegionIteratorWithIndex<TensorImageType> It( outputTensorImage,
                                                               outputTensorImage->GetRequestedRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          TensorPixelType tensor = It.Get();
          typename TensorImageType::IndexType index = It.GetIndex();
          for( unsigned int n = 0; n < NumberOfTensorElements; n++ )
            {
            tensor.SetNthComponent( n, outputImages[n]->GetPixel( index ) );
            }
          It.Set( tensor );
          }
        WriteTensorImage<TensorImageType>( outputTensorImage, ( outputFileName ).c_str(), true );
        }
      else if( inputImageType == 3 )
        {
        unsigned int numberOfTimePoints = timeSeriesImage->GetLargestPossibleRegion().GetSize()[Dimension];

        if( outputImages.size() != numberOfTimePoints )
          {
          antscout << "The number of output images does not match the number of image time points." << std::endl;
          return EXIT_FAILURE;
          }

        typename TimeSeriesImageType::Pointer outputTimeSeriesImage = TimeSeriesImageType::New();

        typename TimeSeriesImageType::PointType origin = timeSeriesImage->GetOrigin();
        typename TimeSeriesImageType::SizeType size = timeSeriesImage->GetLargestPossibleRegion().GetSize();
        typename TimeSeriesImageType::DirectionType direction = timeSeriesImage->GetDirection();
        typename TimeSeriesImageType::IndexType index = timeSeriesImage->GetLargestPossibleRegion().GetIndex();
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          origin[i] = referenceImage->GetOrigin()[i];
          size[i] = referenceImage->GetRequestedRegion().GetSize()[i];
          index[i] = referenceImage->GetRequestedRegion().GetIndex()[i];
          for( unsigned int j = 0; j < Dimension; j++ )
            {
            direction[i][j] = referenceImage->GetDirection()[i][j];
            }
          }

        typename TimeSeriesImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( index );

        int startTimeIndex = timeSeriesImage->GetLargestPossibleRegion().GetIndex()[Dimension];

        outputTimeSeriesImage->CopyInformation( timeSeriesImage );
        outputTimeSeriesImage->SetOrigin( origin );
        outputTimeSeriesImage->SetDirection( direction );
        outputTimeSeriesImage->SetRegions( region );
        outputTimeSeriesImage->Allocate();
        outputTimeSeriesImage->FillBuffer( 0 );

        typename ImageType::IndexType referenceIndex;

        itk::ImageRegionIteratorWithIndex<TimeSeriesImageType> It( outputTimeSeriesImage,
                                                                   outputTimeSeriesImage->GetRequestedRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          typename TimeSeriesImageType::IndexType timeImageIndex = It.GetIndex();
          for( unsigned int i = 0; i < Dimension; i++ )
            {
            referenceIndex[i] = timeImageIndex[i];
            }
          It.Set( outputImages[timeImageIndex[Dimension] - startTimeIndex]->GetPixel( referenceIndex ) );
          }
        WriteImage<TimeSeriesImageType>( outputTimeSeriesImage, ( outputFileName ).c_str() );
        }
      else
        {
        try
          {
          WriteImage<ImageType>( outputImages[0], ( outputFileName ).c_str()  );
          }
        catch( itk::ExceptionObject & err )
          {
          std::cout << "Caught an ITK exception: " << std::endl;
          std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
          throw & err;
          }
        catch( ... )
          {
          std::cout << "Error while writing in image: " << outputFileName << std::endl;
          throw;
          }
        }
      }
    }

  return EXIT_SUCCESS;
}

static void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
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
    option->SetUsageOption( 0, "2/3/4" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Option specifying the input image type of scalar (default), " )
      + std::string( "vector, tensor, or time series." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input-image-type" );
    option->SetShortName( 'e' );
    option->SetUsageOption( 0, "0/1/2/3 " );
    option->SetUsageOption( 1, "scalar/vector/tensor/time-series " );
    option->AddFunction( std::string( "0" ) );
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
    std::string description = std::string( "forces static cast in ReadTransform (for R)" );
    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'z' );
    option->SetDescription( description );
    option->SetLongName( "static-cast-for-R" );
    option->SetUsageOption( 0, "value" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Use 'float' instead of 'double' for computations." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "float" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "antsApplyTransforms, applied to an input image, transforms it " )
    + std::string( "according to a reference image and a transform " )
    + std::string( "(or a set of transforms)." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || ( parser->GetOption( "help" ) &&
                    ( parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) ) ) )
    {
    parser->PrintMenu( antscout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' ) &&
           ( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) ) )
    {
    parser->PrintMenu( antscout, 5, true );
    return EXIT_SUCCESS;
    }

#if 0 // HACK This makes no sense here, filename is never used.
  //Perhaps the "input" option is not needed in this program
  //but is a copy/paste error from another program.
  // Read in the first intensity image to get the image dimension.
  std::string filename;
  itk::ants::CommandLineParser::OptionType::Pointer inputOption =
    parser->GetOption( "reference-image" );
  if( inputOption && inputOption->GetNumberOfFunctions() )
    {
    if( inputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      filename = inputOption->GetFunction( 0 )->GetParameter( 0 );
      }
    else
      {
      filename = inputOption->GetFunction( 0 )->GetName();
      }
    }
  else
    {
    antscout << "No reference image was specified." << std::endl;
    return EXIT_FAILURE;
    }
#endif

  itk::ants::CommandLineParser::OptionType::Pointer inputImageTypeOption =
    parser->GetOption( "input-image-type" );

  std::string inputImageType;
  if( inputImageTypeOption )
    {
    inputImageType = inputImageTypeOption->GetFunction( 0 )->GetName();
    }

  unsigned int dimension = 3;

  // BA - code below creates problems in ANTsR
  //  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
  //                                                            filename.c_str(), itk::ImageIOFactory::ReadMode );
  //  dimension = imageIO->GetNumberOfDimensions();

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }

  bool useDoublePrecision = true;

  std::string precisionType;
  OptionType::Pointer typeOption = parser->GetOption( "float" );
  if( typeOption && parser->Convert<bool>( typeOption->GetFunction( 0 )->GetName() ) )
    {
    antscout << "Using single precision for computations." << std::endl;
    precisionType = "float";
    useDoublePrecision = false;
    }
  else
    {
    antscout << "Using double precision for computations." << std::endl;
    precisionType = "double";
    useDoublePrecision = true;
    }

  enum InputImageType
    {
    SCALAR = 0,
    VECTOR,
    TENSOR,
    TIME_SERIES
    };

  InputImageType imageType = SCALAR;

  if( inputImageTypeOption )
    {
    if( !std::strcmp( inputImageType.c_str(), "scalar" ) || !std::strcmp( inputImageType.c_str(), "0" ) )
      {
      imageType = SCALAR;
      }
    else if( !std::strcmp( inputImageType.c_str(), "vector" ) || !std::strcmp( inputImageType.c_str(), "1" ) )
      {
      imageType = VECTOR;
      }
    else if( !std::strcmp( inputImageType.c_str(), "tensor" ) || !std::strcmp( inputImageType.c_str(), "2" ) )
      {
      imageType = TENSOR;
      }
    else if( !std::strcmp( inputImageType.c_str(), "time-series" ) || !std::strcmp( inputImageType.c_str(), "3" ) )
      {
      imageType = TIME_SERIES;
      }
    else
      {
      antscout << "Unrecognized input image type (cf --input-image-type option)." << std::endl;
      return EXIT_FAILURE;
      }
    }

  switch( dimension )
    {
    case 2:
      {
      if( imageType == TENSOR )
        {
        antscout << "antsApplyTransforms is not implemented for 2-D tensor images." << std::endl;
        return EXIT_FAILURE;
        }
      else
        {
        if( useDoublePrecision )
          {
          antsApplyTransforms<double, 2>( parser, imageType );
          }
        else
          {
          antsApplyTransforms<float, 2>( parser, imageType );
          }
        }
      }
      break;
    case 3:
      {
      if( useDoublePrecision )
        {
        antsApplyTransforms<double, 3>( parser, imageType );
        }
      else
        {
        antsApplyTransforms<float, 3>( parser, imageType );
        }
      }
      break;
    case 4:
      {
      if( imageType == TENSOR )
        {
        antscout << "antsApplyTransforms is not implemented for 4-D tensor images." << std::endl;
        }
      else if( imageType == TIME_SERIES )
        {
        antscout << "antsApplyTransforms is not implemented for 4-D + time images." << std::endl;
        }
      else
        {
        if( useDoublePrecision )
          {
          antsApplyTransforms<double, 4>( parser, imageType );
          }
        else
          {
          antsApplyTransforms<float, 4>( parser, imageType );
          }
        }
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
