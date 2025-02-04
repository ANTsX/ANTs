#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkantsRegistrationHelper.h"
#include "ReadWriteData.h"
#include "TensorFunctions.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformToDisplacementFieldFilter.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkLabelImageGenericInterpolateImageFunction.h"

#include "itksys/SystemTools.hxx"

namespace ants
{
template <typename TensorImageType, typename ImageType>
void
CorrectImageTensorDirection(TensorImageType * movingTensorImage, ImageType * referenceImage)
{
  using DirectionType = typename TensorImageType::DirectionType;
  using MatrixType = typename DirectionType::InternalMatrixType;

  // Assume tensors start in moving voxel space, we want to put them in fixed voxel space

  MatrixType direction =
    referenceImage->GetDirection().GetTranspose() * movingTensorImage->GetDirection().GetVnlMatrix();

  if (!direction.is_identity(0.00001))
  {
    itk::ImageRegionIterator<TensorImageType> It(movingTensorImage, movingTensorImage->GetBufferedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      using TensorType = typename TensorImageType::PixelType;
      using TensorMatrixType = typename TensorImageType::DirectionType::InternalMatrixType;

      TensorType       tensor = It.Get();
      TensorMatrixType dt;

      Vector2Matrix<TensorType, TensorMatrixType>(tensor, dt);

      dt = direction * dt * direction.transpose();

      tensor = Matrix2Vector<TensorType, TensorMatrixType>(dt);

      It.Set(tensor);
    }
  }
}

template <typename DisplacementFieldType, typename ImageType>
void
CorrectImageVectorDirection(DisplacementFieldType * movingVectorImage, ImageType * referenceImage)
{
  using DirectionType = typename DisplacementFieldType::DirectionType;

  // Assume vectors are initially described in the voxel space of the moving image, like tensors
  typename DirectionType::InternalMatrixType direction =
    referenceImage->GetDirection().GetTranspose() * movingVectorImage->GetDirection().GetVnlMatrix();

  using VectorType = typename DisplacementFieldType::PixelType;

  const unsigned int dimension = ImageType::ImageDimension;

  if (!direction.is_identity(0.00001))
  {
    itk::ImageRegionIterator<DisplacementFieldType> It(movingVectorImage, movingVectorImage->GetBufferedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      VectorType vector = It.Get();

      vnl_vector<double> internalVector(dimension);
      for (unsigned int d = 0; d < dimension; d++)
      {
        internalVector[d] = vector[d];
      }

      internalVector.pre_multiply(direction.as_matrix());
      for (unsigned int d = 0; d < dimension; d++)
      {
        vector[d] = internalVector[d];
      }

      It.Set(vector);
    }
  }
}

template <unsigned int NDim>
unsigned int
numTensorElements()
{
  return NDim + numTensorElements<NDim - 1>();
}

template <>
unsigned int
numTensorElements<0>()
{
  return 0;
}

template <unsigned int NDim>
std::vector<unsigned int>
tensorDiagonalArrayIndices()
{
  std::vector<unsigned int> diagElements;
  for (unsigned int d = 0; d < NDim; d++)
  { // I think this is correct for upper-triangular ordering but only tested on 3D tensors
    diagElements.push_back(d * NDim - d * (d - 1) / 2);
  }
  return diagElements;
}

template <unsigned int NDim>
bool
isDiagonalElement(std::vector<unsigned int> diagElements, unsigned int ind)
{
  for (unsigned int i = 0; i < NDim; i++)
  {
    if (diagElements[i] == ind)
    {
      return true;
    }
  }

  return false;
}

template <typename T, unsigned int Dimension, typename OT>
int
antsApplyTransforms(itk::ants::CommandLineParser::Pointer & parser, unsigned int inputImageType = 0)
{
  using RealType = T;
  using PixelType = T;
  using VectorType = itk::Vector<RealType, Dimension>;

  using TensorPixelType = itk::SymmetricSecondRankTensor<RealType, Dimension>;

  // typedef unsigned int                     LabelPixelType;
  // typedef itk::Image<PixelType, Dimension> LabelImageType;

  using ImageType = itk::Image<PixelType, Dimension>;
  using ReferenceImageType = ImageType;

  using TimeSeriesImageType = itk::Image<PixelType, Dimension + 1>;
  using DisplacementFieldType = itk::Image<VectorType, Dimension>;
  using MultiChannelImageType = itk::VectorImage<PixelType, Dimension>;
  using FiveDimensionalImageType = itk::Image<PixelType, 5>;
  using TensorImageType = itk::Image<TensorPixelType, Dimension>;

  using OutputPixelType = OT;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  using OutputTimeSeriesImageType = itk::Image<OutputPixelType, Dimension + 1>;
  using OutputMultiChannelImageType = itk::VectorImage<OutputPixelType, Dimension>;
  using OutputFiveDimensionalImageType = itk::Image<OutputPixelType, 5>;
  using OutputVectorType = itk::Vector<OutputPixelType, Dimension>;
  using OutputDisplacementFieldType = itk::Image<OutputVectorType, Dimension>;

  using RegistrationHelperType = typename ants::RegistrationHelper<T, Dimension>;
  using AffineTransformType = typename RegistrationHelperType::AffineTransformType;
  using CompositeTransformType = typename RegistrationHelperType::CompositeTransformType;

  const unsigned int NumberOfTensorElements = numTensorElements<Dimension>();

  std::vector<unsigned int> tensorDiagIndices = tensorDiagonalArrayIndices<Dimension>();

  typename TimeSeriesImageType::Pointer      timeSeriesImage = nullptr;
  typename TensorImageType::Pointer          tensorImage = nullptr;
  typename MultiChannelImageType::Pointer    multiChannelImage = nullptr;
  typename FiveDimensionalImageType::Pointer fiveDimensionalImage = nullptr;
  typename DisplacementFieldType::Pointer    vectorImage = nullptr;

  std::vector<typename ImageType::Pointer> inputImages;
  inputImages.clear();

  std::vector<typename ImageType::Pointer> outputImages;
  outputImages.clear();

  bool                                                       verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  /**
   * Input object option - for now, we're limiting this to images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption("input");
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");

  // Time-index to extract from time-series image
  unsigned long extractTimeIndex = 0;

  bool extractTimeIndexSet = false;

  typename itk::ants::CommandLineParser::OptionType::Pointer timeIndexOption = parser->GetOption("time-index");
  if (timeIndexOption && timeIndexOption->GetNumberOfFunctions())
  {
    extractTimeIndexSet = true;
    extractTimeIndex = parser->Convert<unsigned int>(timeIndexOption->GetFunction(0)->GetName());
  }

  if (inputImageType == 5 && inputOption && inputOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Input 5-dimensional image: " << inputOption->GetFunction(0)->GetName() << std::endl;
    }
    ReadImage<FiveDimensionalImageType>(fiveDimensionalImage, (inputOption->GetFunction(0)->GetName()).c_str());
    timeSeriesImage =
      ConvertFiveDimensionalImageToTimeSeriesImage<FiveDimensionalImageType, TimeSeriesImageType>(fiveDimensionalImage);
  }
  else if (inputImageType == 4 && inputOption && inputOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Input multichannel image: " << inputOption->GetFunction(0)->GetName() << std::endl;
    }
    ReadImage<MultiChannelImageType>(multiChannelImage, (inputOption->GetFunction(0)->GetName()).c_str());
    timeSeriesImage =
      ConvertMultiChannelImageToTimeSeriesImage<MultiChannelImageType, TimeSeriesImageType>(multiChannelImage);
  }
  else if (inputImageType == 3 && inputOption && inputOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Input time-series image: " << inputOption->GetFunction(0)->GetName() << std::endl;
      if (extractTimeIndexSet)
      {
        std::cout << "Extracting time index: " << extractTimeIndex << std::endl;
      }
    }
    if (!extractTimeIndexSet)
    {
      ReadImage<TimeSeriesImageType>(timeSeriesImage, (inputOption->GetFunction(0)->GetName()).c_str());
    }
    else
    {
      // Modifying inputImageType, because we're going to extract a single time point
      // if we don't do this, the code will try to read from timeSeriesImage below and output a time series
      inputImageType = 0;

      // Set up the image reader with streaming support
      using ReaderType = itk::ImageFileReader<TimeSeriesImageType>;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName((inputOption->GetFunction(0)->GetName()).c_str());

      // Update the reader's output information without loading the entire image into memory
      reader->UpdateOutputInformation();
      typename TimeSeriesImageType::RegionType extractRegion = reader->GetOutput()->GetLargestPossibleRegion();
      typename TimeSeriesImageType::SizeType size = extractRegion.GetSize();

      // Check if the extractTimeIndex is within range
      if (extractTimeIndex >= size[3])
      {
        std::cerr << "Error: time index to extract is out of range [0, " << (size[3] - 1) << "]" << std::endl;
        return EXIT_FAILURE;
      }

      extractRegion.SetIndex(Dimension, extractTimeIndex);
      extractRegion.SetSize(Dimension, 0);
      using ExtracterType = itk::ExtractImageFilter<TimeSeriesImageType, ImageType>;
      typename ExtracterType::Pointer extracter = ExtracterType::New();
      extracter->SetInput(reader->GetOutput());
      extracter->SetExtractionRegion(extractRegion);
      extracter->SetDirectionCollapseToSubmatrix();
      extracter->Update();
      inputImages.push_back(extracter->GetOutput());
    }
  }
  else if (inputImageType == 2 && inputOption && inputOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Input tensor image: " << inputOption->GetFunction(0)->GetName() << std::endl;
    }
    ReadTensorImage<TensorImageType>(tensorImage, (inputOption->GetFunction(0)->GetName()).c_str(), true);
  }
  else if (inputImageType == 0 && inputOption && inputOption->GetNumberOfFunctions())
  {
    const std::string inputFN = inputOption->GetFunction(0)->GetName();

    if (verbose)
    {
      std::cout << "Input scalar image: " << inputFN << std::endl;
    }

    // As this is the default input type, check that input is correct otherwise instruct user to use
    // the -e option

    itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(inputFN.c_str(), itk::IOFileModeEnum::ReadMode);

    if (!imageIO)
    {
      // can only check files, so skip checks if using a pointer in ANTsR / ANTsPy
      if (verbose)
      {
        std::cout << "Input is not a file, assuming dimension = " << Dimension <<
            " and scalar pixel type" << std::endl;
      }
    }
    else
    {
      imageIO->SetFileName(inputFN.c_str());
      imageIO->ReadImageInformation();

      const size_t inputDimension = imageIO->GetNumberOfDimensions();

      // Check if the dimension matches
      if (inputDimension != Dimension)
      {
        if (verbose)
          {
            std::cout << "Input image dimension does not match. Expected: " << Dimension
                    << ", but got: " << inputDimension << std::endl << "See -e option for available input types."
                    << std::endl;
          }
          return EXIT_FAILURE;
      }

      // Check if the pixel type is scalar
      if (imageIO->GetPixelType() != itk::IOPixelEnum::SCALAR)
      {
        if (verbose)
        {
          std::cout << "Image pixel type is not scalar." << std::endl << "See -e option for available input types."
                  << std::endl;
        }
        return EXIT_FAILURE;
      }
    }
    typename ImageType::Pointer image;
    ReadImage<ImageType>(image, inputFN.c_str());
    inputImages.push_back(image);
  }
  else if (inputImageType == 1 && inputOption && inputOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Input vector image: " << inputOption->GetFunction(0)->GetName() << std::endl;
    }

    using ReaderType = itk::ImageFileReader<DisplacementFieldType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName((inputOption->GetFunction(0)->GetName()).c_str());

    try
    {
      vectorImage = reader->GetOutput();
      vectorImage->Update();
      vectorImage->DisconnectPipeline();
    }
    catch (...)
    {
      if (verbose)
      {
        std::cerr << "Unable to read vector image " << reader->GetFileName() << std::endl;
      }
      return EXIT_FAILURE;
    }
  }
  else if (outputOption && outputOption->GetNumberOfFunctions())
  {
    if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
        parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) == 0)
    {
      // Here we have no input image, valid usage is
      // -o [warp.nii.gz, 1] or -o Linear[affine.mat, 0]
      // Exit failure with -o [outputImage.nii.gz,0] but no input to warp
      std::string outputOptionName = outputOption->GetFunction(0)->GetName();
      ConvertToLowerCase(outputOptionName);
      if (std::strcmp(outputOptionName.c_str(), "linear"))
      {
        if (verbose)
        {
          std::cerr << "An input image is required." << std::endl;
        }
        return EXIT_FAILURE;
      }
    }
  }
  /**
   * Reference image option
   */
  bool needReferenceImage = true;

  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    std::string outputOptionName = outputOption->GetFunction(0)->GetName();
    ConvertToLowerCase(outputOptionName);

    // Check if the function name matches keywords requiring no reference image
    if (outputOptionName == "linear" || outputOptionName == "compositetransform")
    {
        if (outputOption->GetFunction(0)->GetNumberOfParameters() > 0)
        {
            needReferenceImage = false;
        }
        else
        {
          // no output parameter where we need one
          if (verbose)
          {
            std::cerr << "An output parameter is required with -o Linear or -o CompositeTransform, see usage" << std::endl;
          }
          return EXIT_FAILURE;
        }
    }
  }

  using ReferenceImageType = ImageType;
  typename ReferenceImageType::Pointer referenceImage;

  typename itk::ants::CommandLineParser::OptionType::Pointer referenceOption = parser->GetOption("reference-image");

  if (referenceOption && referenceOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "Reference image: " << referenceOption->GetFunction(0)->GetName() << std::endl;
    }
    ReadImage<ReferenceImageType>(referenceImage, (referenceOption->GetFunction(0)->GetName()).c_str());
  }
  else if (needReferenceImage == true)
  {
    if (verbose)
    {
      std::cerr << "A reference image is required." << std::endl;
    }
    return EXIT_FAILURE;
  }

  if (inputImageType == 1)
  {
    CorrectImageVectorDirection<DisplacementFieldType, ReferenceImageType>(vectorImage, referenceImage);
    for (unsigned int i = 0; i < Dimension; i++)
    {
      using SelectorType = itk::VectorIndexSelectionCastImageFilter<DisplacementFieldType, ImageType>;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput(vectorImage);
      selector->SetIndex(i);
      selector->Update();

      inputImages.push_back(selector->GetOutput());
    }
  }
  else if (inputImageType == 2)
  {
    CorrectImageTensorDirection<TensorImageType, ReferenceImageType>(tensorImage, referenceImage);
    for (unsigned int i = 0; i < NumberOfTensorElements; i++)
    {
      using SelectorType = itk::VectorIndexSelectionCastImageFilter<TensorImageType, ImageType>;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput(tensorImage);
      selector->SetIndex(i);
      selector->Update();

      inputImages.push_back(selector->GetOutput());
    }
  }
  else if (inputImageType == 3 || inputImageType == 4 || inputImageType == 5)
  {
    typename TimeSeriesImageType::RegionType extractRegion = timeSeriesImage->GetLargestPossibleRegion();
    unsigned int                             numberOfTimePoints = extractRegion.GetSize()[Dimension];
    int                                      startTimeIndex = extractRegion.GetIndex()[Dimension];

    extractRegion.SetSize(Dimension, 0);
    for (unsigned int i = 0; i < numberOfTimePoints; i++)
    {
      extractRegion.SetIndex(Dimension, startTimeIndex + i);

      using ExtracterType = itk::ExtractImageFilter<TimeSeriesImageType, ImageType>;
      typename ExtracterType::Pointer extracter = ExtracterType::New();
      extracter->SetInput(timeSeriesImage);
      extracter->SetExtractionRegion(extractRegion);
      extracter->SetDirectionCollapseToSubmatrix();
      extracter->Update();

      inputImages.push_back(extracter->GetOutput());
    }
  }

  /**
   * Transform option
   */
  // Register the matrix offset transform base class to the
  // transform factory for compatibility with the current ANTs.
  using MatrixOffsetTransformType = itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension>;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

  typename itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption("transform");

  bool                                                       useStaticCastForR = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer rOption = parser->GetOption("static-cast-for-R");
  if (rOption && rOption->GetNumberOfFunctions())
  {
    useStaticCastForR = parser->Convert<bool>(rOption->GetFunction(0)->GetName());
  }

  std::vector<bool>                        isDerivedTransform;
  typename CompositeTransformType::Pointer compositeTransform =
    GetCompositeTransformFromParserOption<RealType, Dimension>(
      parser, transformOption, isDerivedTransform, useStaticCastForR);
  if (compositeTransform.IsNull())
  {
    return EXIT_FAILURE;
  }

  if (!compositeTransform->GetNumberOfParameters())
  {
    if (verbose)
    {
      std::cout << "WARNING: No transforms found, using identify transform" << std::endl;
    }
    typename MatrixOffsetTransformType::Pointer idTransform = MatrixOffsetTransformType::New();
    idTransform->SetIdentity();
    compositeTransform->AddTransform(idTransform);
  }

  std::string                                                whichInterpolator("linear");
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption("interpolation");
  if (interpolationOption && interpolationOption->GetNumberOfFunctions())
  {
    whichInterpolator = interpolationOption->GetFunction(0)->GetName();
    ConvertToLowerCase(whichInterpolator);
  }

  const size_t                    VImageDimension = Dimension;
  typename ImageType::SpacingType cache_spacing_for_smoothing_sigmas(
    itk::NumericTraits<typename ImageType::SpacingType::ValueType>::ZeroValue());
  if (!std::strcmp(whichInterpolator.c_str(), "gaussian") || !std::strcmp(whichInterpolator.c_str(), "multilabel"))
  {
    cache_spacing_for_smoothing_sigmas = inputImages[0]->GetSpacing();
  }

#include "make_interpolator_snip.tmpl"

  /**
   * Default voxel value
   */
  PixelType                                                  defaultValue = 0;
  typename itk::ants::CommandLineParser::OptionType::Pointer defaultOption = parser->GetOption("default-value");
  if (defaultOption && defaultOption->GetNumberOfFunctions())
  {
    defaultValue = parser->Convert<PixelType>(defaultOption->GetFunction(0)->GetName());
  }
  if (verbose)
  {
    if (inputImageType == 2)
      std::cout << "Default pixel mean diffusivity: " << defaultValue << std::endl;
    else
      std::cout << "Default pixel value: " << defaultValue << std::endl;
  }
  for (unsigned int n = 0; n < inputImages.size(); n++)
  {
    using ResamplerType = itk::ResampleImageFilter<ImageType, ImageType, RealType>;
    typename ResamplerType::Pointer resampleFilter = ResamplerType::New();
    resampleFilter->SetInput(inputImages[n]);
    resampleFilter->SetOutputParametersFromImage(referenceImage);
    resampleFilter->SetTransform(compositeTransform);
    if (inputImageType == 2)
    {
      // Set background pixel values in tensor images to produce an isotropic tensor
      if (defaultValue > 0 && isDiagonalElement<Dimension>(tensorDiagIndices, n))
      {
        // defaultValue == MD of isotropic tensor. Resampling is done in log space
        resampleFilter->SetDefaultPixelValue(log(defaultValue));
      }
      else
      {
        resampleFilter->SetDefaultPixelValue(0);
      }
    }
    else
    {
      // for non-tensor images, set the same background value for each component
      resampleFilter->SetDefaultPixelValue(defaultValue);
    }

    interpolator->SetInputImage(inputImages[n]);
    resampleFilter->SetInterpolator(interpolator);
    if (n == 0)
    {
      if (verbose)
      {
        std::cout << "Interpolation type: " << resampleFilter->GetInterpolator()->GetNameOfClass() << std::endl;
      }
    }
    if (inputImageType == 3 || inputImageType == 4 || inputImageType == 5)
    {
      if (verbose)
      {
        std::cout << "  Applying transform(s) to timePoint/channel/dimension5 " << n << " (out of "
                  << inputImages.size() << ")." << std::endl;
      }
    }
    resampleFilter->Update();
    outputImages.push_back(resampleFilter->GetOutput());
  }

  /**
   * output
   */
  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    std::string outputOptionName = outputOption->GetFunction(0)->GetName();
    ConvertToLowerCase(outputOptionName);
    if (!std::strcmp(outputOptionName.c_str(), "linear"))
    {
      if (verbose)
      {
        std::cout << "Output collapsed affine transform: " << outputOption->GetFunction(0)->GetParameter(0) << std::endl;
      }

      if (!compositeTransform->IsLinear())
      {
        if (verbose)
        {
          std::cerr << "The transform or set of transforms is not linear." << std::endl;
        }
        return EXIT_FAILURE;
      }
      else
      {
        typename RegistrationHelperType::Pointer helper = RegistrationHelperType::New();

        typename AffineTransformType::Pointer transform = helper->CollapseLinearTransforms(compositeTransform);

        using TransformWriterType = itk::TransformFileWriterTemplate<T>;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetFileName((outputOption->GetFunction(0)->GetParameter(0)).c_str());

        if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
            parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) != 0)
        {
          typename AffineTransformType::Pointer inverseTransform = AffineTransformType::New();
          inverseTransform->SetMatrix(
            dynamic_cast<MatrixOffsetTransformType *>(transform->GetInverseTransform().GetPointer())->GetMatrix());
          inverseTransform->SetOffset(-(inverseTransform->GetMatrix() * transform->GetOffset()));
          transformWriter->SetInput(inverseTransform);
        }
        else
        {
          transformWriter->SetInput(transform);
        }
        // Use compression just to be consistent
        transformWriter->SetUseCompression(true);
        transformWriter->Update();
      }
    }
    else if (!std::strcmp(outputOptionName.c_str(), "compositetransform"))
    {
      // above should match -o CompositeTransform[output.ext]
      if (verbose)
      {
        std::cout << "Output composite transform: " << outputOption->GetFunction(0)->GetParameter(0) << std::endl;
      }

      using TransformWriterType = itk::TransformFileWriterTemplate<T>;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetFileName((outputOption->GetFunction(0)->GetParameter(0)).c_str());
      transformWriter->SetInput(compositeTransform);
      transformWriter->SetUseCompression(true);
      transformWriter->Update();
    }
    else if ((!std::strcmp(outputOptionName.c_str(), "displacementfield") ||
                (!std::strcmp(outputOptionName.c_str(), "") && outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
                    parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) != 0)))
    {
      // above should match -o [ compositeDisplacement.ext, 1 ] or the preferred -o DisplacementField[output.ext]
      if (verbose)
      {
        std::cout << "Output collapsed transform displacement field: " << outputOption->GetFunction(0)->GetParameter(0)
                  << std::endl;
      }

      using ConverterType = typename itk::TransformToDisplacementFieldFilter<DisplacementFieldType, RealType>;
      typename ConverterType::Pointer converter = ConverterType::New();
      converter->SetOutputOrigin(referenceImage->GetOrigin());
      converter->SetOutputStartIndex(referenceImage->GetBufferedRegion().GetIndex());
      converter->SetSize(referenceImage->GetBufferedRegion().GetSize());
      converter->SetOutputSpacing(referenceImage->GetSpacing());
      converter->SetOutputDirection(referenceImage->GetDirection());
      converter->SetTransform(compositeTransform);
      converter->Update();

      if ((itksys::SystemTools::GetFilenameLastExtension(outputOption->GetFunction(0)->GetParameter(0)) == ".h5") ||
          (itksys::SystemTools::GetFilenameLastExtension(outputOption->GetFunction(0)->GetParameter(0)) == ".hdf5") ||
          (itksys::SystemTools::GetFilenameLastExtension(outputOption->GetFunction(0)->GetParameter(0)) == ".xfm"))
      {
        // unlike NIFTI images, xfm or hdf5 files must be written with a transform writer
        // Note here we are just outputing a single combined displacement field as an hdf5 file, not the full composite
        // transform with its components stored separately
        using DisplacementFieldTransformType = itk::DisplacementFieldTransform<RealType, Dimension>;
        typename DisplacementFieldTransformType::Pointer compositeDisplacementTransform = DisplacementFieldTransformType::New();
        compositeDisplacementTransform->SetDisplacementField(converter->GetOutput());

        using TransformWriterType = itk::TransformFileWriterTemplate<RealType>;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetFileName((outputOption->GetFunction(0)->GetParameter(0)).c_str());
        transformWriter->SetUseCompression(true);
        transformWriter->SetInput(compositeDisplacementTransform);
        transformWriter->Update();
      }
      else
      {
        // write as a vector image
        using DisplacementFieldWriterType = itk::ImageFileWriter<DisplacementFieldType>;
        typename DisplacementFieldWriterType::Pointer displacementFieldWriter = DisplacementFieldWriterType::New();
        displacementFieldWriter->SetInput(converter->GetOutput());
        displacementFieldWriter->SetFileName((outputOption->GetFunction(0)->GetParameter(0)).c_str());
        displacementFieldWriter->Update();
      }
    }
    else
    {

      // if we get here, we are applying the transform to an image, input is required
      if (!(inputOption && inputOption->GetNumberOfFunctions()))
      {
        if (verbose)
        {
          std::cerr << "Output is not a transform, and there is no input image. Exiting" << std::endl;
        }
       return EXIT_FAILURE;
      }

      std::string outputFileName = "";
      if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) == 0)
      {
        outputFileName = outputOption->GetFunction(0)->GetParameter(0);
      }
      else
      {
        outputFileName = outputOption->GetFunction(0)->GetName();
      }
      if (verbose)
      {
        std::cout << "Output warped image: " << outputFileName << std::endl;
      }

      if (inputImageType == 1)
      {
        if (outputImages.size() != Dimension)
        {
          if (verbose)
          {
            std::cerr << "The number of output images does not match the number of vector components." << std::endl;
          }
          return EXIT_FAILURE;
        }

        typename DisplacementFieldType::Pointer outputVectorImage = DisplacementFieldType::New();
        outputVectorImage->CopyInformation(referenceImage);
        outputVectorImage->SetRegions(referenceImage->GetRequestedRegion());
        outputVectorImage->AllocateInitialized();

        itk::ImageRegionIteratorWithIndex<DisplacementFieldType> It(outputVectorImage,
                                                                    outputVectorImage->GetRequestedRegion());
        for (It.GoToBegin(); !It.IsAtEnd(); ++It)
        {
          VectorType                                vector = It.Get();
          typename DisplacementFieldType::IndexType index = It.GetIndex();
          for (unsigned int n = 0; n < Dimension; n++)
          {
            vector.SetNthComponent(n, outputImages[n]->GetPixel(index));
          }
          It.Set(vector);
        }

        using CastFilterType = itk::CastImageFilter<DisplacementFieldType, OutputDisplacementFieldType>;
        typename CastFilterType::Pointer caster = CastFilterType::New();
        caster->SetInput(outputVectorImage);
        caster->Update();

        using WriterType = itk::ImageFileWriter<OutputDisplacementFieldType>;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput(caster->GetOutput());
        writer->SetFileName((outputFileName).c_str());
        writer->Update();
      }
      else if (inputImageType == 2)
      {
        if (outputImages.size() != NumberOfTensorElements)
        {
          if (verbose)
          {
            std::cerr << "The number of output images does not match the number of tensor elements." << std::endl;
          }
          return EXIT_FAILURE;
        }

        typename TensorImageType::Pointer outputTensorImage = TensorImageType::New();
        outputTensorImage->CopyInformation(referenceImage);
        outputTensorImage->SetRegions(referenceImage->GetRequestedRegion());
        outputTensorImage->AllocateInitialized();

        itk::ImageRegionIteratorWithIndex<TensorImageType> It(outputTensorImage,
                                                              outputTensorImage->GetRequestedRegion());
        for (It.GoToBegin(); !It.IsAtEnd(); ++It)
        {
          TensorPixelType                     tensor = It.Get();
          typename TensorImageType::IndexType index = It.GetIndex();
          for (unsigned int n = 0; n < NumberOfTensorElements; n++)
          {
            tensor.SetNthComponent(n, outputImages[n]->GetPixel(index));
          }
          It.Set(tensor);
        }
        WriteTensorImage<TensorImageType>(outputTensorImage, (outputFileName).c_str(), true);
      }
      else if (inputImageType == 3 || inputImageType == 4 || inputImageType == 5)
      {
        unsigned int numberOfTimePoints = timeSeriesImage->GetLargestPossibleRegion().GetSize()[Dimension];

        if (outputImages.size() != numberOfTimePoints)
        {
          if (verbose)
          {
            std::cerr << "The number of output images does not match the number of image time/channel points."
                      << std::endl;
          }
          return EXIT_FAILURE;
        }

        typename TimeSeriesImageType::Pointer outputTimeSeriesImage = TimeSeriesImageType::New();

        typename TimeSeriesImageType::PointType     origin = timeSeriesImage->GetOrigin();
        typename TimeSeriesImageType::SizeType      size = timeSeriesImage->GetLargestPossibleRegion().GetSize();
        typename TimeSeriesImageType::DirectionType direction = timeSeriesImage->GetDirection();
        typename TimeSeriesImageType::IndexType     index = timeSeriesImage->GetLargestPossibleRegion().GetIndex();
        typename TimeSeriesImageType::SpacingType   spacing = timeSeriesImage->GetSpacing();

        for (unsigned int i = 0; i < Dimension; i++)
        {
          origin[i] = referenceImage->GetOrigin()[i];
          size[i] = referenceImage->GetRequestedRegion().GetSize()[i];
          index[i] = referenceImage->GetRequestedRegion().GetIndex()[i];
          spacing[i] = referenceImage->GetSpacing()[i];
          for (unsigned int j = 0; j < Dimension; j++)
          {
            direction[i][j] = referenceImage->GetDirection()[i][j];
          }
        }

        typename TimeSeriesImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(index);

        int startTimeIndex = timeSeriesImage->GetLargestPossibleRegion().GetIndex()[Dimension];

        outputTimeSeriesImage->CopyInformation(timeSeriesImage);
        outputTimeSeriesImage->SetOrigin(origin);
        outputTimeSeriesImage->SetDirection(direction);
        outputTimeSeriesImage->SetSpacing(spacing);
        outputTimeSeriesImage->SetRegions(region);
        outputTimeSeriesImage->AllocateInitialized();

        typename ImageType::IndexType referenceIndex;

        itk::ImageRegionIteratorWithIndex<TimeSeriesImageType> It(outputTimeSeriesImage,
                                                                  outputTimeSeriesImage->GetRequestedRegion());
        for (It.GoToBegin(); !It.IsAtEnd(); ++It)
        {
          typename TimeSeriesImageType::IndexType timeImageIndex = It.GetIndex();
          for (unsigned int i = 0; i < Dimension; i++)
          {
            referenceIndex[i] = timeImageIndex[i];
          }
          It.Set(outputImages[timeImageIndex[Dimension] - startTimeIndex]->GetPixel(referenceIndex));
        }
        if (inputImageType == 3)
        {
          using CastFilterType = itk::CastImageFilter<TimeSeriesImageType, OutputTimeSeriesImageType>;
          typename CastFilterType::Pointer caster = CastFilterType::New();
          caster->SetInput(outputTimeSeriesImage);
          caster->Update();

          ANTs::WriteImage<OutputTimeSeriesImageType>(caster->GetOutput(), (outputFileName).c_str());
        }
        else if (inputImageType == 4)
        {
          multiChannelImage = ConvertTimeSeriesImageToMultiChannelImage<TimeSeriesImageType, MultiChannelImageType>(
            outputTimeSeriesImage);

          using CastFilterType = itk::CastImageFilter<MultiChannelImageType, OutputMultiChannelImageType>;
          typename CastFilterType::Pointer caster = CastFilterType::New();
          caster->SetInput(multiChannelImage);
          caster->Update();

          ANTs::WriteImage<OutputMultiChannelImageType>(caster->GetOutput(), (outputFileName).c_str());
        }
        else //  if( inputImageType == 5 )
        {
          fiveDimensionalImage =
            ConvertTimeSeriesImageToFiveDimensionalImage<TimeSeriesImageType, FiveDimensionalImageType>(
              outputTimeSeriesImage);

          using CastFilterType = itk::CastImageFilter<FiveDimensionalImageType, OutputFiveDimensionalImageType>;
          typename CastFilterType::Pointer caster = CastFilterType::New();
          caster->SetInput(fiveDimensionalImage);
          caster->Update();

          ANTs::WriteImage<OutputFiveDimensionalImageType>(caster->GetOutput(), (outputFileName).c_str());
        }
      }
      else
      {
        try
        {
          using CastFilterType = itk::CastImageFilter<ImageType, OutputImageType>;
          typename CastFilterType::Pointer caster = CastFilterType::New();
          caster->SetInput(outputImages[0]);
          caster->Update();

          ANTs::WriteImage<OutputImageType>(caster->GetOutput(), (outputFileName).c_str());
        }
        catch (const itk::ExceptionObject & err)
        {
          if (verbose)
          {
            std::cerr << "Caught an ITK exception: " << std::endl;
            std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
          }
          return EXIT_FAILURE;
          // throw &err;
        }
        catch (...)
        {
          if (verbose)
          {
            std::cerr << "Error while writing in image: " << outputFileName << std::endl;
          }
          return EXIT_FAILURE;
          // throw;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

static void
antsApplyTransformsInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  {
    std::string description = std::string("Sets the dimensionality of transforms and scalar inputs.") +
                              std::string("dimensional image.  This will be the same dimension as used in ") +
                              std::string("antsRegistration. This does not change for multi-valued inputs, use ") +
                              std::string("the -e option for time series and other multi-component images.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3/4");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Option specifying the input image type of scalar (default), ") +
                              std::string("vector, tensor, time series, or multi-channel.  A time series ") +
                              std::string("image is a scalar image defined by an additional dimension ") +
                              std::string("for the time component whereas a multi-channel image is a ") +
                              std::string("vector image with only spatial dimensions.  Five-dimensional") +
                              std::string("images are e.g., AFNI stats image.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input-image-type");
    option->SetShortName('e');
    option->SetUsageOption(0, "0/1/2/3/4/5");
    option->SetUsageOption(1, "scalar/vector/tensor/time-series/multichannel/five-dimensional");
    option->AddFunction(std::string("0"));
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Time index to extract from time series input (-e 3). ") +
                              std::string("This selects a single slice from time series input without reading ") +
                              std::string("the entire dataset into memory.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("time-index");
    option->SetUsageOption(0, "<timeIndex>");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Currently, the only input objects supported are image ") +
                              std::string("objects.  However, the current framework allows for ") +
                              std::string("warping of other objects such as meshes and point sets. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("For warping input images, the reference image defines the ") +
                              std::string("spacing, origin, size, and direction of the output warped ") +
                              std::string("image. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("reference-image");
    option->SetShortName('r');
    option->SetUsageOption(0, "imageFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
        "The default output is the warped input image, resliced into the space of the reference image. Alternatively, one "
        "can output the composite transform containing all input transforms, or a collapsed affine or displacement field "
        "created from the input transforms. "
        "If all input transforms are linear, they can be collapsed and (if boolean is set) inverted, then written to "
        "an ITK transform file. "
        "If the input transforms contain displacement fields, they can be collapsed into a single displacement field with "
        "DisplacementField[collapsedDFFileName]. This replaces the usage '-o [collapsedDFFileName, 1]', which is "
        "deprecated but still allowed for backwards compatibility. "
        "To write a non-collapsed composite transform, use CompositeTransform[compositeTransform.h5].";

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "warpedOutputFileName");
    option->SetUsageOption(1, "DisplacementField[compositeDFFileName]");
    option->SetUsageOption(2, "CompositeTransform[compositeTransformFileName]");
    option->SetUsageOption(3, "Linear[collapsedAffine.mat, invert=0]");
    option->SetUsageOption(4, "[collapsedDisplacementField,1]");

    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Several interpolation options are available in ITK. ") +
                              std::string("These have all been made available.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("interpolation");
    option->SetShortName('n');
    option->SetUsageOption(0, "Linear");
    option->SetUsageOption(1, "NearestNeighbor");
    option->SetUsageOption(2, "MultiLabel[<sigma=imageSpacing>,<alpha=4.0>]");
    option->SetUsageOption(3, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]");
    option->SetUsageOption(4, "BSpline[<order=3>]");
    option->SetUsageOption(5, "CosineWindowedSinc");
    option->SetUsageOption(6, "WelchWindowedSinc");
    option->SetUsageOption(7, "HammingWindowedSinc");
    option->SetUsageOption(8, "LanczosWindowedSinc");
    option->SetUsageOption(9, "GenericLabel[<interpolator=Linear>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Output image data type. ") +
                              std::string("This is a direct typecast; output values are not rescaled. ") +
                              std::string("Default is to use the internal data type (float or double). ") +
                              std::string("uchar is unsigned char; others are signed. ") +
                              std::string("WARNING: Outputs will be incorrect (overflowed/reinterpreted) ") +
                              std::string("if values exceed the range allowed by your choice. ") +
                              std::string("Note that some pixel types are not supported by some image ") +
                              std::string("formats. e.g. int is not supported by jpg. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output-data-type");
    option->SetShortName('u');
    option->SetUsageOption(0, "char");
    option->SetUsageOption(1, "uchar");
    option->SetUsageOption(2, "short");
    option->SetUsageOption(3, "int");
    option->SetUsageOption(4, "float");
    option->SetUsageOption(5, "double");
    option->SetUsageOption(6, "default");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Several transform options are supported including all ") +
                              std::string("those defined in the ITK library in addition to ") +
                              std::string("a deformation field transform.  The ordering of ") +
                              std::string("the transformations follows the ordering specified ") +
                              std::string("on the command line.  An identity transform is pushed ") +
                              std::string("onto the transformation stack. Each new transform ") +
                              std::string("encountered on the command line is also pushed onto ") +
                              std::string("the transformation stack. Then, to warp the input object, ") +
                              std::string("each point comprising the input object is warped first ") +
                              std::string("according to the last transform pushed onto the stack ") +
                              std::string("followed by the second to last transform, etc. until ") +
                              std::string("the last transform encountered which is the identity ") +
                              std::string("transform. ") +
                              std::string("Also, it should be noted that the inverse transform can ") +
                              std::string("be accommodated with the usual caveat that such an inverse ") +
                              std::string("must be defined by the specified transform class ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transform");
    option->SetShortName('t');
    option->SetUsageOption(0, "transformFileName");
    option->SetUsageOption(1, "[transformFileName,useInverse]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Default voxel value to be used with input images only. ") +
                              std::string("Specifies the voxel value when the input point maps outside ") +
                              std::string("the output domain. With tensor input images, specifies the ") +
                              std::string("default voxel eigenvalues. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("default-value");
    option->SetShortName('f');
    option->SetUsageOption(0, "value");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("forces static cast in ReadTransform (for R)");
    OptionType::Pointer option = OptionType::New();
    option->SetShortName('z');
    option->SetDescription(description);
    option->SetLongName("static-cast-for-R");
    option->SetUsageOption(0, "value");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Use 'float' instead of 'double' for computations.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("float");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Verbose output.");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('v');
    option->SetLongName("verbose");
    option->SetUsageOption(0, "(0)/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsApplyTransforms(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsApplyTransforms");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand(argv[0]);

  std::string commandDescription = std::string("antsApplyTransforms, applied to an input image, transforms it ") +
                                   std::string("according to a reference image and a transform ") +
                                   std::string("(or a set of transforms).");

  parser->SetCommandDescription(commandDescription);
  antsApplyTransformsInitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  bool                                              verbose = false;
  itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  if (argc == 1)
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_FAILURE;
  }
  else if (parser->GetOption("help")->GetFunction() &&
           parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_SUCCESS;
  }
  else if (parser->GetOption('h')->GetFunction() &&
           parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  itk::ants::CommandLineParser::OptionType::Pointer inputImageTypeOption = parser->GetOption("input-image-type");

  std::string inputImageType;
  if (inputImageTypeOption)
  {
    inputImageType = inputImageTypeOption->GetFunction(0)->GetName();
  }

  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption("dimensionality");
  if (dimOption && dimOption->GetNumberOfFunctions())
  {
    dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
  }

  bool useDoublePrecision = true;

  std::string         precisionType;
  OptionType::Pointer typeOption = parser->GetOption("float");
  if (typeOption && parser->Convert<bool>(typeOption->GetFunction(0)->GetName()))
  {
    if (verbose)
    {
      std::cout << "Using single precision for computations." << std::endl;
    }
    precisionType = "float";
    useDoublePrecision = false;
  }
  else
  {
    if (verbose)
    {
      std::cout << "Using double precision for computations." << std::endl;
    }
    precisionType = "double";
    useDoublePrecision = true;
  }

  enum InputImageType
  {
    SCALAR = 0,
    VECTOR,
    TENSOR,
    TIME_SERIES,
    MULTICHANNEL,
    FIVEDIMENSIONAL
  };

  InputImageType imageType = SCALAR;

  if (inputImageTypeOption)
  {
    if (!std::strcmp(inputImageType.c_str(), "scalar") || !std::strcmp(inputImageType.c_str(), "0"))
    {
      imageType = SCALAR;
    }
    else if (!std::strcmp(inputImageType.c_str(), "vector") || !std::strcmp(inputImageType.c_str(), "1"))
    {
      imageType = VECTOR;
    }
    else if (!std::strcmp(inputImageType.c_str(), "tensor") || !std::strcmp(inputImageType.c_str(), "2"))
    {
      imageType = TENSOR;
    }
    else if (!std::strcmp(inputImageType.c_str(), "time-series") || !std::strcmp(inputImageType.c_str(), "3"))
    {
      imageType = TIME_SERIES;
    }
    else if (!std::strcmp(inputImageType.c_str(), "multichannel") || !std::strcmp(inputImageType.c_str(), "4"))
    {
      imageType = MULTICHANNEL;
    }
    else if (!std::strcmp(inputImageType.c_str(), "five-dimensional") || !std::strcmp(inputImageType.c_str(), "5"))
    {
      imageType = FIVEDIMENSIONAL;
    }
    else
    {
      if (verbose)
      {
        std::cerr << "Unrecognized input image type (cf --input-image-type option)." << std::endl;
      }
      return EXIT_FAILURE;
    }
  }

  std::string                                                outputDataType("default");
  typename itk::ants::CommandLineParser::OptionType::Pointer outputDataTypeOption =
    parser->GetOption("output-data-type");
  if (outputDataTypeOption && outputDataTypeOption->GetNumberOfFunctions())
  {
    outputDataType = outputDataTypeOption->GetFunction(0)->GetName();
    ConvertToLowerCase(outputDataType);
  }

  switch (dimension)
  {
    case 2:
    {
      if (imageType == TENSOR)
      {
        if (verbose)
        {
          std::cerr << "antsApplyTransforms is not implemented for 2-D tensor images." << std::endl;
        }
        return EXIT_FAILURE;
      }
      else
      {
        if (useDoublePrecision)
        {
          if (!std::strcmp(outputDataType.c_str(), "char"))
            return antsApplyTransforms<double, 2, char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "uchar"))
            return antsApplyTransforms<double, 2, unsigned char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "short"))
            return antsApplyTransforms<double, 2, short>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "int"))
            return antsApplyTransforms<double, 2, int>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "float"))
            return antsApplyTransforms<double, 2, float>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "double"))
            return antsApplyTransforms<double, 2, double>(parser, imageType);
          else
            return antsApplyTransforms<double, 2, double>(parser, imageType);
        }
        else
        {
          if (!std::strcmp(outputDataType.c_str(), "char"))
            return antsApplyTransforms<float, 2, char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "uchar"))
            return antsApplyTransforms<float, 2, unsigned char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "short"))
            return antsApplyTransforms<float, 2, short>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "int"))
            return antsApplyTransforms<float, 2, int>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "float"))
            return antsApplyTransforms<float, 2, float>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "double"))
            return antsApplyTransforms<float, 2, double>(parser, imageType);
          else
            return antsApplyTransforms<float, 2, float>(parser, imageType);
        }
      }
    }
    break;
    case 3:
    {
      if (useDoublePrecision)
      {
        if (!std::strcmp(outputDataType.c_str(), "char"))
          return antsApplyTransforms<double, 3, char>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "uchar"))
          return antsApplyTransforms<double, 3, unsigned char>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "short"))
          return antsApplyTransforms<double, 3, short>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "int"))
          return antsApplyTransforms<double, 3, int>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "float"))
          return antsApplyTransforms<double, 3, float>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "double"))
          return antsApplyTransforms<double, 3, double>(parser, imageType);
        else
          return antsApplyTransforms<double, 3, double>(parser, imageType);
      }
      else
      {
        if (!std::strcmp(outputDataType.c_str(), "char"))
          return antsApplyTransforms<float, 3, char>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "uchar"))
          return antsApplyTransforms<float, 3, unsigned char>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "short"))
          return antsApplyTransforms<float, 3, short>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "int"))
          return antsApplyTransforms<float, 3, int>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "float"))
          return antsApplyTransforms<float, 3, float>(parser, imageType);
        else if (!std::strcmp(outputDataType.c_str(), "double"))
          return antsApplyTransforms<float, 3, double>(parser, imageType);
        else
          return antsApplyTransforms<float, 3, float>(parser, imageType);
      }
    }
    break;
    case 4:
    {
      if (imageType == TENSOR)
      {
        if (verbose)
        {
          std::cerr << "antsApplyTransforms is not implemented for 4-D tensor images." << std::endl;
        }
      }
      else if (imageType == TIME_SERIES)
      {
        if (verbose)
        {
          std::cerr << "antsApplyTransforms is not implemented for 4-D + time images." << std::endl;
        }
      }
      else
      {
        if (useDoublePrecision)
        {
          if (!std::strcmp(outputDataType.c_str(), "char"))
            return antsApplyTransforms<double, 4, char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "uchar"))
            return antsApplyTransforms<double, 4, unsigned char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "short"))
            return antsApplyTransforms<double, 4, short>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "int"))
            return antsApplyTransforms<double, 4, int>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "float"))
            return antsApplyTransforms<double, 4, float>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "double"))
            return antsApplyTransforms<double, 4, double>(parser, imageType);
          else
            return antsApplyTransforms<double, 4, double>(parser, imageType);
        }
        else
        {
          if (!std::strcmp(outputDataType.c_str(), "char"))
            return antsApplyTransforms<float, 4, char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "uchar"))
            return antsApplyTransforms<float, 4, unsigned char>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "short"))
            return antsApplyTransforms<float, 4, short>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "int"))
            return antsApplyTransforms<float, 4, int>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "float"))
            return antsApplyTransforms<float, 4, float>(parser, imageType);
          else if (!std::strcmp(outputDataType.c_str(), "double"))
            return antsApplyTransforms<float, 4, double>(parser, imageType);
          else
            return antsApplyTransforms<float, 4, float>(parser, imageType);
        }
      }
    }
    break;
    default:
    {
      if (verbose)
      {
        std::cerr << "Unsupported dimension" << std::endl;
      }
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

} // namespace ants
