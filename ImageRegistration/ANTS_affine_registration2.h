#ifndef ANTS_AFFINE_REGISTRATION2_H_
#define ANTS_AFFINE_REGISTRATION2_H_

#include <vector>
#include <cstdlib>
#include <ctime>
#include "itkImage.h"
#include "itkPoint.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkCorrelationCoefficientHistogramImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkANTSAffine3DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"
#include "itkWarpImageWAffineFilter.h"
#include "itkImageMomentsCalculator.h"
#include <vector>
#include "ReadWriteData.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkGradientDifferenceImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkImageRegionIterator.h"
#include "itkRandomImageSource.h"
#include "itkAddImageFilter.h"

typedef enum
{
  AffineWithMutualInformation = 1,
  AffineWithMeanSquareDifference,
  AffineWithHistogramCorrelation,
  AffineWithNormalizedCorrelation,
  AffineWithGradientDifference
} AffineMetricType;

template <typename TAffineTransform, typename TMaskImage>
class OptAffine
{
public:
  typedef typename TAffineTransform::Pointer AffineTransformPointerType;
  typedef typename TMaskImage::Pointer       MaskImagePointerType;
  typedef TAffineTransform                   AffineTransformType;
  typedef TMaskImage                         MaskImageType;
  typedef typename MaskImageType::Pointer    MaskObjectPointerType;
  OptAffine()
    : metric_type(AffineWithMutualInformation)
  {
    MI_bins = 32;
    MI_samples = 6000;
    number_of_seeds = 0;
    time_seed = (unsigned int)time(nullptr);
    number_of_levels = 3;
    number_of_iteration_list.resize(number_of_levels, 10000);
    const int kParaDim = AffineTransformType::ParametersDimension;
    gradient_scales.resize(kParaDim, 1.0);
    is_rigid = false;

    maximum_step_length = 0.1;
    relaxation_factor = 0.5;
    minimum_step_length = 1.e-5;
    translation_scales = 1.e-4;

    use_rotation_header = false;
    ignore_void_orgin = true;
  };

  ~OptAffine() = default;

  AffineTransformPointerType transform_initial;
  MaskImagePointerType       mask_fixed;

  int                 MI_bins;
  int                 MI_samples;
  int                 number_of_seeds;
  unsigned int        time_seed;
  int                 number_of_levels;
  std::vector<int>    number_of_iteration_list;
  std::vector<double> gradient_scales;
  AffineMetricType    metric_type;
  bool                is_rigid;

  double maximum_step_length;
  double relaxation_factor;
  double minimum_step_length;
  double translation_scales;

  bool use_rotation_header;
  bool ignore_void_orgin;
};

template <typename TAffineTransform, typename TMaskImage>
std::ostream &
operator<<(std::ostream & os, const OptAffine<TAffineTransform, TMaskImage> & p)
{
  os << "OptAffine: ";
  os << "metric_type=";

  switch (p.metric_type)
  {
    case AffineWithMutualInformation:
      os << "AffineWithMutualInformation" << std::endl;
      break;
    case AffineWithMeanSquareDifference:
      os << "AffineWithMeanSquareDifference" << std::endl;
      break;
    case AffineWithHistogramCorrelation:
      os << "AffineWithHistogramCorrelation" << std::endl;
      break;
    case AffineWithNormalizedCorrelation:
      os << "AffineWithNormalizedCorrelation" << std::endl;
      break;
    case AffineWithGradientDifference:
      os << "AffineWithGradientDifference" << std::endl;
      break;
  }
  os << "MI_bins=" << p.MI_bins << " "
     << "MI_samples=" << p.MI_samples << std::endl;
  os << "number_of_seeds=" << p.number_of_seeds << " "
     << "time_seed=" << p.time_seed << std::endl;
  os << "number_of_levels=" << p.number_of_levels << std::endl;
  os << "number_of_iteration_list="
     << "[";
  for (unsigned int i = 0; i < p.number_of_iteration_list.size() - 1; i++)
  {
    os << p.number_of_iteration_list[i] << ",";
  }
  if (p.number_of_iteration_list.size() > 0)
  {
    os << p.number_of_iteration_list[p.number_of_iteration_list.size() - 1];
  }
  os << "]" << std::endl;
  os << "graident_scales="
     << "[";
  for (unsigned int i = 0; i < p.gradient_scales.size() - 1; i++)
  {
    os << p.gradient_scales[i] << ",";
  }
  if (p.gradient_scales.size() > 0)
  {
    os << p.gradient_scales[p.gradient_scales.size() - 1];
  }
  os << "]" << std::endl;
  os << "is_rigid = " << p.is_rigid << std::endl;
  os << "mask null: " << p.mask_fixed.IsNull() << std::endl;

  os << "maximum_step_length=" << p.maximum_step_length << std::endl;
  ;
  os << "relaxation_factor=" << p.relaxation_factor << std::endl;
  os << "minimum_step_length=" << p.minimum_step_length << std::endl;
  os << "translation_scales=" << p.translation_scales << std::endl;

  return os;
}

template <typename TransformType>
void
WriteAffineTransformFile(typename TransformType::Pointer & transform, const std::string & filename)
{
  itk::TransformFileWriter::Pointer transform_writer;

  transform_writer = itk::TransformFileWriter::New();
  transform_writer->SetFileName(filename);
  transform_writer->SetInput(transform);
#if ITK_VERSION_MAJOR >= 5
  transform_writer->SetUseCompression(true);
#endif

  try
  {
    transform_writer->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << err << std::endl
              << "Exception in writing transform file: " << std::endl
              << filename << std::endl;
    return;
  }

  return;
}

template <typename CastTransformType>
void
ReadAffineTransformFile(const std::string & filename, typename CastTransformType::Pointer & transform)
{
  //    const unsigned int InputSpaceDimension = CastTransformType::InputSpaceDimension;
  //    const unsigned int OutputSpaceDimension = CastTransformType::OutputSpaceDimension;

  itk::TransformFactory<CastTransformType>::RegisterTransform();
  itk::TransformFactory<itk::ANTSAffine3DTransform<double>>::RegisterTransform();

  typedef typename itk::TransformFileReader TranReaderType;
  TranReaderType::Pointer                   tran_reader = TranReaderType::New();
  tran_reader->SetFileName(filename);

  try
  {
    tran_reader->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << err << std::endl;
    std::cout << "Exception caught in reading tran para file: " << filename << std::endl;
    return;
  }

  transform = dynamic_cast<CastTransformType *>((tran_reader->GetTransformList())->front().GetPointer());

  return;
}

template <typename OptAffine>
void
InitializeAffineOptmizationParameters(OptAffine & opt, double translationScale)
{
  const int kImageDim = OptAffine::MaskImageType::ImageDimension;

  switch (kImageDim)
  {
    case 2:
    {
      //        const double translationScale = 1.0 / 1000.0;
      opt.gradient_scales[0] = 1.0;
      opt.gradient_scales[1] = 1.0;
      opt.gradient_scales[2] = 1.0;
      opt.gradient_scales[3] = 1.0;
      opt.gradient_scales[4] = translationScale;
      opt.gradient_scales[5] = translationScale;
      opt.gradient_scales[6] = translationScale;
      opt.gradient_scales[7] = translationScale;
    }
    break;
    case 3:
    {
      //        const double translationScale = 1.0/1.e4;
      opt.gradient_scales[0] = 1.0; // quaternion
      opt.gradient_scales[1] = 1.0; // quaternion
      opt.gradient_scales[2] = 1.0; // quaternion
      opt.gradient_scales[3] = 1.0; // quaternion
      opt.gradient_scales[4] = 1.0; // s1
      opt.gradient_scales[5] = 1.0; // s2
      opt.gradient_scales[6] = 1.0; // s3
      opt.gradient_scales[7] = 1.0; // k1
      opt.gradient_scales[8] = 1.0; // k2
      opt.gradient_scales[9] = 1.0; // k3
      opt.gradient_scales[10] = translationScale;
      opt.gradient_scales[11] = translationScale;
      opt.gradient_scales[12] = translationScale;
    }
    break;
  }

  std::cout << opt;
}

template <typename TMaskObjectType, typename TImagePyramid, typename TMetricType, typename TInterpolatorType>
class RunningAffineCache
{
public:
  typedef TImagePyramid                         ImagePyramidType;
  typedef typename ImagePyramidType::value_type ImagePointerType;
  typedef typename ImagePointerType::ObjectType ImageType;
  typedef TMetricType                           MetricType;
  typedef typename MetricType::Pointer          MetricPointerType;
  typedef TInterpolatorType                     InterpolatorType;
  typedef typename TInterpolatorType::Pointer   InterpolatorPointerType;
  typedef TMaskObjectType                       MaskObjectType;
  typedef typename MaskObjectType::Pointer      MaskObjectPointerType;

  RunningAffineCache() = default;
  ~RunningAffineCache() = default;

  MaskObjectPointerType   mask_fixed_object;
  ImagePyramidType        fixed_image_pyramid;
  ImagePyramidType        moving_image_pyramid;
  MetricPointerType       metric;
  MetricPointerType       invmetric;
  InterpolatorPointerType interpolator;
};

template <typename ImageTypePointer, typename AffineTransformPointer>
void
GetAffineTransformFromImage(const ImageTypePointer & img, AffineTransformPointer & aff)
{
  typedef typename ImageTypePointer::ObjectType                        ImageType;
  typedef typename ImageType::DirectionType                            DirectionType;
  typedef typename ImageType::PointType                                PointType;
  typedef typename AffineTransformPointer::ObjectType::TranslationType VectorType;

  DirectionType direction = img->GetDirection();
  PointType     pt = img->GetOrigin();
  VectorType    translation;
  translation.Fill(0);
  aff->SetMatrix(direction);
  aff->SetCenter(pt);
  aff->SetTranslation(translation);
}

// //////////////////////////////////////////////////////////////////////
template <typename ImagePointerType,
          typename RunningImagePointerType,
          typename OptAffineType,
          typename RunningOptAffineType>
inline void
PreConversionInAffine(ImagePointerType &        fixedImage,
                      RunningImagePointerType & R_fixedImage,
                      ImagePointerType &        movingImage,
                      RunningImagePointerType & R_movingImage,
                      OptAffineType &           opt,
                      RunningOptAffineType &    R_opt)
{
  typedef typename OptAffineType::AffineTransformType        AffineTransformType;
  typedef typename RunningOptAffineType::AffineTransformType RunningAffineTransformType;

  if (opt.use_rotation_header)
  {
    std::cout << "===================>initialize from rotation header ... " << std::endl;
    // use the rotation header to initialize the affine: inv(Tm) * Tf
    typename AffineTransformType::Pointer aff_Im = AffineTransformType::New();
    GetAffineTransformFromImage(movingImage, aff_Im);
    typename AffineTransformType::Pointer aff_If = AffineTransformType::New();
    GetAffineTransformFromImage(fixedImage, aff_If);
    typename AffineTransformType::Pointer aff_combined = AffineTransformType::New();
    aff_combined->SetFixedParameters(aff_If->GetFixedParameters());
    aff_combined->SetParameters(aff_If->GetParameters());
    typename AffineTransformType::Pointer aff_Im_inv = AffineTransformType::New();
    aff_Im->GetInverse(aff_Im_inv);
    aff_combined->Compose(aff_Im_inv, 0);
    opt.transform_initial = aff_combined;

    //            std::cout << "aff_If: " << aff_If << std::endl;
    //            std::cout << "aff_Im: " << aff_Im << std::endl;
    //            std::cout << "aff_combined: " << aff_combined << std::endl;
  }

  if (!opt.use_rotation_header && opt.ignore_void_orgin)
  {
    std::cout
      << "===================> ignore void origins which are too far away to be possible alignments: use 0 instead."
      << std::endl;
    typename AffineTransformType::Pointer aff_Im = AffineTransformType::New();
    GetAffineTransformFromImage(movingImage, aff_Im);
    typename AffineTransformType::Pointer aff_If = AffineTransformType::New();
    GetAffineTransformFromImage(fixedImage, aff_If);

    bool b_far_origin_without_rotation = false;
    // bool b_far_origin_without_rotation = HaveFarOriginWithoutRotation(aff_If, aff_Im);

    if (b_far_origin_without_rotation)
    {
      typename AffineTransformType::Pointer aff_combined = AffineTransformType::New();
      aff_combined->SetFixedParameters(aff_If->GetFixedParameters());
      aff_combined->SetParameters(aff_If->GetParameters());
      typename AffineTransformType::Pointer aff_Im_inv = AffineTransformType::New();
      aff_Im->GetInverse(aff_Im_inv);
      aff_combined->Compose(aff_Im_inv, 0);
      opt.transform_initial = aff_combined;
    }
  }

  if (opt.transform_initial.IsNotNull())
  {
    R_opt.transform_initial = RunningAffineTransformType::New();

    R_opt.transform_initial->SetCenter(*(reinterpret_cast<typename RunningAffineTransformType::InputPointType *>(
      const_cast<typename AffineTransformType::InputPointType *>(&(opt.transform_initial->GetCenter())))));
    R_opt.transform_initial->SetMatrix(*(reinterpret_cast<typename RunningAffineTransformType::MatrixType *>(
      const_cast<typename AffineTransformType::MatrixType *>(&(opt.transform_initial->GetMatrix())))));
    R_opt.transform_initial->SetTranslation(*(reinterpret_cast<typename RunningAffineTransformType::OutputVectorType *>(
      const_cast<typename AffineTransformType::OutputVectorType *>(&(opt.transform_initial->GetTranslation())))));
  }

  // std::cout << "R_opt.transform_initial" << R_opt.transform_initial << std::endl;

  if (opt.mask_fixed.IsNotNull())
  {
    R_opt.mask_fixed = dynamic_cast<typename RunningOptAffineType::MaskImageType *>(opt.mask_fixed.GetPointer());
    if (R_opt.mask_fixed.IsNull())
    {
      itkGenericExceptionMacro(<< "Can't convert optimizer mask to proper mask type.");
    }
    // have to set " -fno-strict-aliasing " in gcc to remove the following compilation warning:
    //  warning: dereferencing type-punned pointer will break strict-aliasing rules
  }

  R_fixedImage = reinterpret_cast<RunningImagePointerType &>(fixedImage);
  R_movingImage = reinterpret_cast<RunningImagePointerType &>(movingImage);

  R_opt.MI_bins = opt.MI_bins;
  R_opt.MI_samples = opt.MI_samples;
  R_opt.number_of_seeds = opt.number_of_seeds;
  R_opt.time_seed = opt.time_seed;
  R_opt.number_of_levels = opt.number_of_levels;
  R_opt.number_of_iteration_list = opt.number_of_iteration_list;
  // R_opt.gradient_scales = opt.gradient_scales; // does not need, will assign value later.
  R_opt.metric_type = opt.metric_type;
  R_opt.is_rigid = opt.is_rigid;

  R_opt.maximum_step_length = opt.maximum_step_length;
  R_opt.relaxation_factor = opt.relaxation_factor;
  R_opt.minimum_step_length = opt.minimum_step_length;
  R_opt.translation_scales = opt.translation_scales;

  R_opt.use_rotation_header = opt.use_rotation_header;
  R_opt.ignore_void_orgin = opt.ignore_void_orgin;
}

// //////////////////////////////////////////////////////////////////////
template <typename RunningAffineTransformPointerType, typename AffineTransformPointerType>
inline void
PostConversionInAffine(RunningAffineTransformPointerType & transform_running, AffineTransformPointerType & transform)
{
  typedef typename RunningAffineTransformPointerType::ObjectType RunningAffineTransformType;
  typedef typename AffineTransformPointerType::ObjectType        AffineTransformType;

  transform->SetCenter(*(reinterpret_cast<typename AffineTransformType::InputPointType *>(
    const_cast<typename RunningAffineTransformType::InputPointType *>(&(transform_running->GetCenter())))));
  transform->SetTranslation(*(reinterpret_cast<typename AffineTransformType::OutputVectorType *>(
    const_cast<typename RunningAffineTransformType::OutputVectorType *>(&(transform_running->GetTranslation())))));
  transform->SetMatrix(*(reinterpret_cast<typename AffineTransformType::MatrixType *>(
    const_cast<typename RunningAffineTransformType::MatrixType *>(&(transform_running->GetMatrix())))));

  // std::cout << "transform_running" << transform_running << std::endl;
  // std::cout << "transform" << transform << std::endl;
}

// /////////////////////////////////////////////////////////////////////////////
template <typename ImageType, typename TransformType, typename OptAffineType>
void
ComputeSingleAffineTransform2D3D(typename ImageType::Pointer &     fixed_image,
                                 typename ImageType::Pointer &     moving_image,
                                 OptAffineType &                   opt,
                                 typename TransformType::Pointer & transform)
{
  const int ImageDimension = ImageType::ImageDimension;

  typedef std::vector<typename ImageType::Pointer>               ImagePyramidType;
  typedef itk::ImageMaskSpatialObject<ImageDimension>            ImageMaskSpatialObjectType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;

  typedef typename TransformType::ParametersType ParaType;

  InitializeAffineOptmizationParameters(opt, opt.translation_scales);

  // std::cout << "DEBUG: opt.gradient_scales.size() = " << opt.gradient_scales.size() << std::endl;

  InitializeAffineTransform(fixed_image, moving_image, opt);

  std::cout << "input affine center: " << opt.transform_initial->GetCenter() << std::endl;
  std::cout << "input affine para: " << opt.transform_initial->GetParameters() << std::endl;

  transform = TransformType::New();
  ParaType para_final(TransformType::ParametersDimension);

  switch (opt.metric_type)
  {
    case AffineWithMeanSquareDifference:
    {
      typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
      typedef RunningAffineCache<ImageMaskSpatialObjectType, ImagePyramidType, MetricType, InterpolatorType>
        RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);
      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
    }
    break;
    case AffineWithHistogramCorrelation:
    {
      typedef itk::CorrelationCoefficientHistogramImageToImageMetric<ImageType, ImageType> MetricType;
      typedef RunningAffineCache<ImageMaskSpatialObjectType, ImagePyramidType, MetricType, InterpolatorType>
        RunningAffineCacheType;

      RunningAffineCacheType running_cache;

      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);
      unsigned int                                 nBins = 32;
      typename MetricType::HistogramType::SizeType histSize;
      histSize[0] = nBins;
      histSize[1] = nBins;
      running_cache.metric->SetHistogramSize(histSize);
      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
    }
    break;
    case AffineWithNormalizedCorrelation:
    {
      typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
      typedef RunningAffineCache<ImageMaskSpatialObjectType, ImagePyramidType, MetricType, InterpolatorType>
        RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);
      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
    }
    break;
    case AffineWithGradientDifference:
    {
      typedef itk::GradientDifferenceImageToImageMetric<ImageType, ImageType> MetricType;
      typedef RunningAffineCache<ImageMaskSpatialObjectType, ImagePyramidType, MetricType, InterpolatorType>
        RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);
      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
    }
    break;
    case AffineWithMutualInformation:
    {
      typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> MetricType;
      typedef RunningAffineCache<ImageMaskSpatialObjectType, ImagePyramidType, MetricType, InterpolatorType>
        RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);

      running_cache.metric->SetNumberOfHistogramBins(opt.MI_bins);
      running_cache.metric->SetNumberOfSpatialSamples(opt.MI_samples);

      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
    }
    break;
    default:
      break;
  }

  bool noaffine = true;
  for (int i = 0; i < opt.number_of_levels; i++)
  {
    if (opt.number_of_iteration_list[i] > 0)
    {
      noaffine = false;
    }
  }
  if (noaffine)
  {
    for (unsigned int i = TransformType::ParametersDimension - ImageDimension; i < TransformType::ParametersDimension;
         i++)
    {
      para_final[i] = 0;
    }
  }

  transform->SetParameters(para_final);
  transform->SetCenter(opt.transform_initial->GetCenter());

  double rval_init = TestCostValueMMI(
    fixed_image, moving_image, opt.transform_initial->GetParameters(), opt.transform_initial->GetCenter(), transform);

  // std::cout << "ABCDABCD: " << transform << std::endl;

  double rval_final =
    TestCostValueMMI(fixed_image, moving_image, para_final, opt.transform_initial->GetCenter(), transform);

  std::cout << "outputput affine center: " << transform->GetCenter() << std::endl;
  std::cout << "output affine para: " << transform->GetParameters() << std::endl;
  std::cout << "initial measure value (MMI): rval = " << rval_init << std::endl;
  std::cout << "final measure value (MMI): rval = " << rval_final << std::endl;
  std::cout << "finish affine registeration..." << std::endl;
}

// /////////////////////////////////////////////////////////////////////////
// the initial transform maybe any derivative class type from MatrixOffsetTransformBase,
// it will be automatically converted to the my 2D/3D affine type
template <typename ImageType, typename TransformType, typename OptAffineType>
void
ComputeSingleAffineTransform(typename ImageType::Pointer &     fixedImage,
                             typename ImageType::Pointer &     movingImage,
                             OptAffineType &                   opt,
                             typename TransformType::Pointer & transform)
{
  const int ImageDimension = ImageType::ImageDimension;

  typedef typename ImageType::IOPixelType PixelType;

  std::cout << "transform_initial: IsNotNull():" << opt.transform_initial.IsNotNull() << std::endl;

  if (ImageDimension == 2)
  {
    typedef itk::ANTSCenteredAffine2DTransform<double>   RunningAffineTransformType;
    typedef typename RunningAffineTransformType::Pointer RunningAffineTransformPointerType;
    constexpr unsigned int                               RunningImageDimension = 2;

    typedef typename itk::Image<PixelType, RunningImageDimension>   RunningImageType;
    typedef typename RunningImageType::Pointer                      RunningImagePointerType;
    typedef OptAffine<RunningAffineTransformType, RunningImageType> RunningOptAffineType;

    RunningImagePointerType R_fixedImage, R_movingImage;
    RunningOptAffineType    R_opt;

    PreConversionInAffine(fixedImage, R_fixedImage, movingImage, R_movingImage, opt, R_opt);

    RunningAffineTransformPointerType transform_running = nullptr;
    ComputeSingleAffineTransform2D3D<RunningImageType, RunningAffineTransformType, RunningOptAffineType>(
      R_fixedImage, R_movingImage, R_opt, transform_running);

    PostConversionInAffine(transform_running, transform);
  }
  else if (ImageDimension == 3)
  {
    typedef itk::ANTSAffine3DTransform<double>           RunningAffineTransformType;
    typedef typename RunningAffineTransformType::Pointer RunningAffineTransformPointerType;
    constexpr unsigned int                               RunningImageDimension = 3;

    typedef typename itk::Image<PixelType, RunningImageDimension>   RunningImageType;
    typedef typename RunningImageType::Pointer                      RunningImagePointerType;
    typedef OptAffine<RunningAffineTransformType, RunningImageType> RunningOptAffineType;

    RunningImagePointerType R_fixedImage, R_movingImage;
    RunningOptAffineType    R_opt;

    PreConversionInAffine(fixedImage, R_fixedImage, movingImage, R_movingImage, opt, R_opt);

    RunningAffineTransformPointerType transform_running = nullptr;
    ComputeSingleAffineTransform2D3D<RunningImageType, RunningAffineTransformType, RunningOptAffineType>(
      R_fixedImage, R_movingImage, R_opt, transform_running);

    PostConversionInAffine(transform_running, transform);
  }
  else
  {
    std::cout << "Unsupported, not 2D/ 3D" << std::endl;
    return;
  }
}

// /////////////////////////////////////////////////////////////////////////////
template <typename MaskImagePointerType, typename ImageMaskSpatialObjectPointerType>
void
InitialzeImageMask(MaskImagePointerType & mask_fixed, ImageMaskSpatialObjectPointerType & mask_fixed_object)
{
  if (mask_fixed.IsNull())
  {
    return;
  }

  const unsigned int                                ImageDimension = MaskImagePointerType::ObjectType::ImageDimension;
  typedef typename MaskImagePointerType::ObjectType MaskImageType;
  typedef typename ImageMaskSpatialObjectPointerType::ObjectType ImageMaskSpatialObjectType;

  typedef itk::Image<unsigned char, ImageDimension>              CharMaskImageType;
  typedef itk::CastImageFilter<MaskImageType, CharMaskImageType> CastFilterType;
  typename CastFilterType::Pointer                               cast_filter = CastFilterType::New();
  cast_filter->SetInput(mask_fixed);
  cast_filter->Update();
  typename CharMaskImageType::Pointer mask_fixed_char = cast_filter->GetOutput();

  mask_fixed_object = ImageMaskSpatialObjectType::New();
  mask_fixed_object->SetImage(mask_fixed_char);
}

// /////////////////////////////////////////////////////////////////////////////
template <typename ImagePointerType>
ImagePointerType
AddRandomNoise(ImagePointerType & I)
{
  typedef typename ImagePointerType::ObjectType ImageType;

  typename itk::RandomImageSource<ImageType>::Pointer randomImageSource = itk::RandomImageSource<ImageType>::New();

  randomImageSource->SetOrigin(I->GetOrigin());
  randomImageSource->SetSpacing(I->GetSpacing());
  randomImageSource->SetSize(I->GetBufferedRegion().GetSize());
  randomImageSource->SetMin(0);
  randomImageSource->SetMax(0.001);
  randomImageSource->Update();
  randomImageSource->GetOutput()->SetDirection(I->GetDirection());

  typedef itk::AddImageFilter<ImageType, ImageType> AddImageFilterType;

  typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
  addFilter->SetInput1(I);
  addFilter->SetInput2(randomImageSource->GetOutput());
  addFilter->Update();

  return addFilter->GetOutput();
}

// /////////////////////////////////////////////////////////////////////////////
template <typename ImagePointerType, typename PointType, typename VectorType>
void
ComputeInitialPosition(ImagePointerType & I_fixed,
                       ImagePointerType & I_moving,
                       PointType &        center,
                       VectorType &       translation_vec)
{
  typedef typename ImagePointerType::ObjectType           ImageType;
  typedef typename itk::ImageMomentsCalculator<ImageType> ImageCalculatorType;

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typename ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();

  typename ImageCalculatorType::VectorType fixed_center;
  typename ImageCalculatorType::VectorType moving_center;

  // a dirty fix to handle the constant/blank images after preprocessing
  try
  {
    calculator->SetImage(I_fixed);
    calculator->Compute();
    fixed_center = calculator->GetCenterOfGravity();

    calculator->SetImage(I_moving);
    calculator->Compute();
    moving_center = calculator->GetCenterOfGravity();
  }
  catch (...)
  {
    // try to add a small amount of noise to avoid exception from computing moments
    std::cout << "try to add a small amount of noise to avoid exception"
                 " from computing moments"
              << std::endl;
    ImagePointerType If1 = AddRandomNoise(I_fixed);
    ImagePointerType Im1 = AddRandomNoise(I_moving);

    //  calculator->SetImage(  I_fixed );
    calculator->SetImage(If1);
    calculator->Compute();
    fixed_center = calculator->GetCenterOfGravity();

    //  calculator->SetImage(  I_moving );
    calculator->SetImage(Im1);
    calculator->Compute();
    moving_center = calculator->GetCenterOfGravity();
  }
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    center[i] = fixed_center[i];
    translation_vec[i] = moving_center[i] - fixed_center[i];
  }
}

// fake a all-zero vector
template <typename ImagePointerType, typename PointType, typename VectorType>
void
ComputeInitialPosition_tmp(ImagePointerType & I_fixed,
                           ImagePointerType & I_moving,
                           PointType &        center,
                           VectorType &       translation_vec)
{
  typedef typename ImagePointerType::ObjectType           ImageType;
  typedef typename itk::ImageMomentsCalculator<ImageType> ImageCalculatorType;

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typename ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();

  calculator->SetImage(I_fixed);
  calculator->Compute();
  typename ImageCalculatorType::VectorType fixed_center = calculator->GetCenterOfGravity();

  calculator->SetImage(I_moving);
  calculator->Compute();
  typename ImageCalculatorType::VectorType moving_center = calculator->GetCenterOfGravity();
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    center[i] = fixed_center[i];
    translation_vec[i] = fixed_center[i] - fixed_center[i];
  }
}

// /////////////////////////////////////////////////////////////////////////////
template <typename PointType, typename VectorType, typename TransformPointerType>
void
InjectInitialPara(PointType & center, VectorType & translation_vec, TransformPointerType & transform)
{
  typedef typename TransformPointerType::ObjectType::ParametersType ParaType;
  ParaType para0(TransformPointerType::ObjectType::ParametersDimension);

  switch ((unsigned int)PointType::PointDimension)
  {
    case 2:
      para0[0] = 0;                  // para1[0]; // theta
      para0[1] = 1.0;                // s1
      para0[2] = 1.0;                // s2
      para0[3] = 0.0;                // k
      para0[4] = center[0];          // para1[1]; //c1
      para0[5] = center[1];          // para1[2]; //c2
      para0[6] = translation_vec[0]; // 0;//para1[3]; //t1
      para0[7] = translation_vec[1]; // 0; //para1[4]; //t2

      transform->SetParameters(para0);
      transform->SetCenter(center);

      break;
    case 3:
      para0[0] = 0.0;
      para0[1] = 0.0;
      para0[2] = 0.0;
      para0[3] = 1.0;
      para0[4] = 1.0;
      para0[5] = 1.0;
      para0[6] = 1.0;
      para0[7] = 0.0;
      para0[8] = 0.0;
      para0[9] = 0.0;
      para0[10] = translation_vec[0];
      para0[11] = translation_vec[1];
      para0[12] = translation_vec[2];
      // para0[10] = 0.0; para0[11] = 0.0;   para0[12] = 0.0;

      transform->SetParameters(para0);
      transform->SetCenter(center);

      break;
  }
}

// ////////////////////////////////////////////////////////////////////////////////////////
template <typename ImagePointerType, typename ParaType, typename PointType, typename TransformTypePointer>
double
TestCostValueMMI(ImagePointerType fixedImage,
                 ImagePointerType movingImage,
                 ParaType         para,
                 PointType        center,
                 TransformTypePointer /* null_transform */)
{
  typedef typename ImagePointerType::ObjectType     ImageType;
  typedef typename TransformTypePointer::ObjectType TransformType;

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetCenter(center);
  // transform->SetParameters(para);

  typedef typename itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> mattesMutualInfoMetricType;
  typename mattesMutualInfoMetricType::Pointer mattesMutualInfo = mattesMutualInfoMetricType::New();

  typedef typename itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer                                      interpolator = InterpolatorType::New();

  interpolator->SetInputImage(movingImage);

  mattesMutualInfo->SetFixedImage(fixedImage);
  mattesMutualInfo->SetMovingImage(movingImage);
  mattesMutualInfo->SetFixedImageRegion(fixedImage->GetBufferedRegion());
  mattesMutualInfo->SetTransform(transform);
  mattesMutualInfo->SetInterpolator(interpolator);
  mattesMutualInfo->SetNumberOfHistogramBins(32);
  mattesMutualInfo->SetNumberOfSpatialSamples(5000);
  mattesMutualInfo->SetTransformParameters(para);
  mattesMutualInfo->Initialize();
  double rval = 0;
  try
  {
    rval = mattesMutualInfo->GetValue(para);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << err << std::endl
              << "Exception caught in computing mattesMutualInfo after registration" << std::endl
              << "Maybe: Too many samples map outside moving image buffer" << std::endl
              << "Set the cost value = 0 (max for MutualInfo) " << std::endl;
    rval = 0;
  }

  return rval;
}

template <typename ImagePointer>
ImagePointer
ShrinkImageToScale(ImagePointer image, float scalingFactor)
{
  typedef typename ImagePointer::ObjectType ImageType;
  typedef typename ImageType::PixelType     RealType;

  typedef typename ImagePointer::ObjectType ImageType;
  typename ImageType::SpacingType           inputSpacing = image->GetSpacing();
  typename ImageType::RegionType::SizeType  inputSize = image->GetRequestedRegion().GetSize();

  typename ImageType::SpacingType          outputSpacing;
  typename ImageType::RegionType::SizeType outputSize;

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer                   resampler = ResampleFilterType::New();

  RealType minimumSpacing = inputSpacing.GetVnlVector().min_value();
  //    RealType maximumSpacing = inputSpacing.GetVnlVector().max_value();

  ImagePointer current_image = image;
  for (unsigned int d = 0; d < ImageType::ImageDimension; d++)
  {
    RealType scaling = static_cast<RealType>(
      std::min(scalingFactor * static_cast<float>(minimumSpacing) / static_cast<float>(inputSpacing[d]),
               static_cast<float>(inputSize[d]) / 32.0f));
    outputSpacing[d] = inputSpacing[d] * static_cast<double>(scaling);
    outputSize[d] =
      static_cast<unsigned long>(static_cast<RealType>(inputSpacing[d]) * static_cast<RealType>(inputSize[d]) /
                                   static_cast<RealType>(outputSpacing[d]) +
                                 static_cast<RealType>(0.5));

    typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
    typename GaussianFilterType::Pointer                            smoother = GaussianFilterType::New();
    smoother->SetInputImage(current_image);
    smoother->SetDirection(d);
    smoother->SetNormalizeAcrossScale(false);
    smoother->SetSigma(0.25 * (outputSpacing[d] / inputSpacing[d]));

    if (smoother->GetSigma() > 0.0)
    {
      smoother->Update();
      current_image = smoother->GetOutput();
    }
  }

  resampler->SetInput(current_image);
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(image->GetOrigin());
  resampler->SetOutputDirection(image->GetDirection());

  resampler->Update();

  image = resampler->GetOutput();

  // std::cout << "DEBUG: " << outputSize << std::endl;

  return image;
}

template <typename ImagePointerType, typename ImagePyramidType>
void
BuildImagePyramid(const ImagePointerType & image, int number_of_levels, ImagePyramidType & image_pyramid)
{
  image_pyramid.resize(number_of_levels);

  image_pyramid[number_of_levels - 1] = image;
  double scale_factor = 2;
  for (int i = 0; i < number_of_levels - 1; i++)
  {
    image_pyramid[number_of_levels - 2 - i] = ShrinkImageToScale(image, scale_factor);
    scale_factor *= 2;
  }

  //    for(int i=0; i < number_of_levels; i++)
  //        std::cout << "level " << i << ": size: " << image_pyramid[i]->GetLargestPossibleRegion().GetSize() <<
  // std::endl;
}

template <typename ImagePointerType, typename OptAffineType, typename RunningAffineCacheType>
void
InitializeRunningAffineCache(ImagePointerType &       fixed_image,
                             ImagePointerType &       moving_image,
                             OptAffineType &          opt,
                             RunningAffineCacheType & running_cache)
{
  typedef typename RunningAffineCacheType::InterpolatorType InterpolatorType;
  typedef typename RunningAffineCacheType::MetricType       MetricType;

  BuildImagePyramid(fixed_image, opt.number_of_levels, running_cache.fixed_image_pyramid);
  BuildImagePyramid(moving_image, opt.number_of_levels, running_cache.moving_image_pyramid);
  InitialzeImageMask(opt.mask_fixed, running_cache.mask_fixed_object);

  running_cache.interpolator = InterpolatorType::New();
  running_cache.metric = MetricType::New();
  running_cache.invmetric = MetricType::New();
}

template <typename ImagePointerType, typename OptAffineType>
void
InitializeAffineTransform(ImagePointerType & fixed_image, ImagePointerType & moving_image, OptAffineType & opt)
{
  typedef typename OptAffineType::AffineTransformType TransformType;
  typedef typename TransformType::InputPointType      PointType;
  typedef typename TransformType::OutputVectorType    VectorType;

  std::cout << "opt.transform_initial.IsNull(): " << opt.transform_initial.IsNull() << std::endl;
  std::cout << " opt.use_rotation_header: " << opt.use_rotation_header << std::endl;
  std::cout << " opt.ignore_void_orgin: " << opt.ignore_void_orgin << std::endl;

  if (opt.transform_initial.IsNull())
  {
    PointType  center;
    VectorType translation_vec;
    // std::cout << "GS: debug: fake a all zero translation_vec" << std::endl;
    // ComputeInitialPosition_tmp(fixed_image, moving_image, center, translation_vec);
    ComputeInitialPosition(fixed_image, moving_image, center, translation_vec);
    opt.transform_initial = TransformType::New();
    InjectInitialPara(center, translation_vec, opt.transform_initial);
  }
}

template <typename ParaType>
ParaType
NormalizeGradientForRigidTransform(ParaType & original_gradient, int kImageDim)
{
  ParaType new_gradient(original_gradient.Size());

  new_gradient = original_gradient;

  switch (kImageDim)
  {
    case 2: // theta, s1, s2, k
      for (int j = 1; j <= 3; j++)
      {
        new_gradient[j] = 0.;
      }
      break;
    case 3: // q1,q2,q3,q4,s1,s2,s3,k1,k2,k3
      for (int j = 4; j <= 9; j++)
      {
        new_gradient[j] = 0.;
      }
      break;
  }
  return new_gradient;
}

// /////////////////////////////////////////////////////////////////////////////
// template<typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
template <typename RunningAffineCacheType, typename OptAffine, typename ParaType>
bool
SymmRegisterImageAffineMutualInformationMultiResolution(RunningAffineCacheType & running_cache,
                                                        OptAffine &              opt,
                                                        ParaType &               para_final)
{
  typedef typename RunningAffineCacheType::ImagePyramidType      ImagePyramidType;
  typedef typename RunningAffineCacheType::ImageType             ImageType;
  typedef typename RunningAffineCacheType::ImagePointerType      ImagePointerType;
  typedef typename RunningAffineCacheType::MetricType            MetricType;
  typedef typename RunningAffineCacheType::InterpolatorType      InterpolatorType;
  typedef typename OptAffine::AffineTransformType                TransformType;
  typedef typename RunningAffineCacheType::MaskObjectPointerType MaskObjectPointerType;

  const unsigned int kImageDim = ImageType::ImageDimension;

  ImagePyramidType &      fixed_image_pyramid = running_cache.fixed_image_pyramid;
  ImagePyramidType &      moving_image_pyramid = running_cache.moving_image_pyramid;
  MaskObjectPointerType & mask_fixed_object = running_cache.mask_fixed_object;

  int                   number_of_levels = opt.number_of_levels;
  std::vector<int> &    number_of_iteration_list = opt.number_of_iteration_list;
  std::vector<double> & gradient_scales = opt.gradient_scales;
  bool                  is_rigid = opt.is_rigid;

  // use my own's registration routine of image pyramid and gradient descent , only use ITK's implementation of mutual
  // information
  // try to use my own transform class together with image pyramid when transform are required (not much for MI though)

  typename TransformType::Pointer transform = TransformType::New();
  typename TransformType::Pointer invtransform = TransformType::New();

  typename InterpolatorType::Pointer & interpolator = running_cache.interpolator;

  typename InterpolatorType::Pointer invinterpolator = InterpolatorType::New();

  typename MetricType::Pointer & metric = running_cache.metric;
  typename MetricType::Pointer & invmetric = running_cache.invmetric;

  const int kParaDim = TransformType::ParametersDimension;

  ParaType current_para(kParaDim);
  current_para = opt.transform_initial->GetParameters();

  double maximum_step_length = opt.maximum_step_length;
  double relaxation_factor = opt.relaxation_factor;
  double minimum_step_length = opt.minimum_step_length;
  double current_step_length;
  for (int i = 0; i < number_of_levels; i++)
  {
    transform->SetParameters(current_para);
    transform->SetCenter(opt.transform_initial->GetCenter());

    ImagePointerType fixed_image = fixed_image_pyramid[i];
    ImagePointerType moving_image = moving_image_pyramid[i];
    int              number_of_iteration_current_level = number_of_iteration_list[i];
    interpolator->SetInputImage(moving_image);
    metric->SetMovingImage(moving_image);
    metric->SetFixedImage(fixed_image);
    metric->SetTransform(transform);
    metric->SetInterpolator(interpolator);
    metric->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion());

    if (mask_fixed_object.IsNotNull())
    {
      metric->SetFixedImageMask(mask_fixed_object);
    }
    metric->Initialize();

    ParaType last_gradient(kParaDim);
    ParaType invlast_gradient(kParaDim);
    for (int j = 0; j < kParaDim; j++)
    {
      last_gradient[j] = 0;
      invlast_gradient[j] = 0;
    }
    current_step_length = maximum_step_length;

    bool is_converged = false;
    int  used_iterations = 0;
    for (used_iterations = 0; used_iterations < number_of_iteration_current_level; used_iterations++)
    {
      transform->GetInverse(invtransform);
      invinterpolator->SetInputImage(fixed_image);
      invmetric->SetMovingImage(fixed_image);
      invmetric->SetFixedImage(moving_image);
      invmetric->SetTransform(invtransform);
      invmetric->SetInterpolator(invinterpolator);
      invmetric->SetFixedImageRegion(moving_image->GetLargestPossibleRegion());
      invmetric->Initialize();

      ParaType original_gradient(kParaDim);
      ParaType current_gradient(kParaDim);
      ParaType invoriginal_gradient(kParaDim);
      ParaType invcurrent_gradient(kParaDim);

      double value, invvalue;
      try
      {
        metric->GetValueAndDerivative(current_para, value, original_gradient);
        invmetric->GetValueAndDerivative(invtransform->GetParameters(), invvalue, invoriginal_gradient);
      }
      catch (const itk::ExceptionObject & err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return false;
      }

      // use the similar routine as RegularStepGradientDescentBaseOptimizer::AdvanceOneStep
      // to use oscillation as the minimization convergence
      // notice this is always a minimization procedure
      if (is_rigid)
      {
        original_gradient = NormalizeGradientForRigidTransform(original_gradient, kImageDim);
        invoriginal_gradient = NormalizeGradientForRigidTransform(invoriginal_gradient, kImageDim);
      }
      for (int j = 0; j < kParaDim; j++)
      {
        current_gradient[j] = original_gradient[j] / gradient_scales[j];
      }
      double gradient_magnitude = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        gradient_magnitude += current_gradient[j] * current_gradient[j];
      }
      gradient_magnitude = sqrt(gradient_magnitude);
      double inner_product_last_current_gradient = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        inner_product_last_current_gradient += current_gradient[j] * last_gradient[j];
      }
      if (inner_product_last_current_gradient < 0)
      {
        current_step_length *= relaxation_factor;
      }

      if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
      {
        is_converged = true;
        break;
      }
      // for inverse
      for (int j = 0; j < kParaDim; j++)
      {
        invcurrent_gradient[j] = invoriginal_gradient[j] / gradient_scales[j];
      }
      double invgradient_magnitude = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        invgradient_magnitude += invcurrent_gradient[j] * invcurrent_gradient[j];
      }
      invgradient_magnitude = sqrt(invgradient_magnitude);
      inner_product_last_current_gradient = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        inner_product_last_current_gradient += invcurrent_gradient[j] * invlast_gradient[j];
      }
      if (inner_product_last_current_gradient < 0)
      {
        current_step_length *= relaxation_factor;
      }

      if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
      {
        is_converged = true;
        break;
      }
      for (int j = 0; j < kParaDim; j++)
      {
        current_para[j] += (-1.0) * current_gradient[j] * current_step_length / gradient_magnitude;
        current_para[j] += (1.0) * invcurrent_gradient[j] * current_step_length / invgradient_magnitude;
      }
      if (kImageDim == 3) // normalize quaternion
      {
        double quat_mag = 0.0;
        for (int j = 0; j < 4; j++)
        {
          quat_mag += current_para[j] * current_para[j];
        }
        quat_mag = sqrt(quat_mag);
        for (int j = 0; j < 4; j++)
        {
          current_para[j] /= quat_mag;
        }
        if (!is_rigid)
        {
          for (int j = 4; j < 7; j++)
          {
            current_para[j] *= quat_mag;
          }
        }
      }

      last_gradient = current_gradient;
      invlast_gradient = invcurrent_gradient;
    }

    std::cout << "level " << i << ", iter " << used_iterations << ", size: fix"
              << fixed_image->GetRequestedRegion().GetSize() << "-mov" << moving_image->GetRequestedRegion().GetSize();

    std::cout << ", affine para: " << current_para << std::endl;

    if (is_converged)
    {
      std::cout << "    reach oscillation, current step: " << current_step_length << "<" << minimum_step_length
                << std::endl;
    }
    else
    {
      std::cout << "    does not reach oscillation, current step: " << current_step_length << ">" << minimum_step_length
                << std::endl;
    }
  }
  para_final = current_para;
  return true;
}

// /////////////////////////////////////////////////////////////////////////////
// template<typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
template <typename RunningAffineCacheType, typename OptAffine, typename ParaType>
bool
RegisterImageAffineMutualInformationMultiResolution(RunningAffineCacheType & running_cache,
                                                    OptAffine &              opt,
                                                    ParaType &               para_final)
{
  typedef typename RunningAffineCacheType::ImagePyramidType      ImagePyramidType;
  typedef typename RunningAffineCacheType::ImagePointerType      ImagePointerType;
  typedef typename RunningAffineCacheType::ImageType             ImageType;
  typedef typename RunningAffineCacheType::MetricType            MetricType;
  typedef typename RunningAffineCacheType::InterpolatorType      InterpolatorType;
  typedef typename OptAffine::AffineTransformType                TransformType;
  typedef typename RunningAffineCacheType::MaskObjectPointerType MaskObjectPointerType;

  const unsigned int kImageDim = ImageType::ImageDimension;

  ImagePyramidType &      fixed_image_pyramid = running_cache.fixed_image_pyramid;
  ImagePyramidType &      moving_image_pyramid = running_cache.moving_image_pyramid;
  MaskObjectPointerType & mask_fixed_object = running_cache.mask_fixed_object;

  int                   number_of_levels = opt.number_of_levels;
  std::vector<int> &    number_of_iteration_list = opt.number_of_iteration_list;
  std::vector<double> & gradient_scales = opt.gradient_scales;
  bool                  is_rigid = opt.is_rigid;

  // use my own's registration routine of image pyramid and gradient descent , only use ITK's implementation of mutual
  // information
  // try to use my own transform class together with image pyramid when transform are required (not much for MI though)

  typename TransformType::Pointer transform = TransformType::New();

  typename InterpolatorType::Pointer & interpolator = running_cache.interpolator;
  typename MetricType::Pointer &       metric = running_cache.metric;
  const int                            kParaDim = TransformType::ParametersDimension;

  ParaType current_para(kParaDim);
  current_para = opt.transform_initial->GetParameters();

  double maximum_step_length = opt.maximum_step_length;
  double relaxation_factor = opt.relaxation_factor;
  double minimum_step_length = opt.minimum_step_length;
  double current_step_length;
  double value = 0;
  for (int i = 0; i < number_of_levels; i++)
  {
    transform->SetParameters(current_para);
    transform->SetCenter(opt.transform_initial->GetCenter());

    ImagePointerType fixed_image = fixed_image_pyramid[i];
    ImagePointerType moving_image = moving_image_pyramid[i];
    int              number_of_iteration_current_level = number_of_iteration_list[i];
    interpolator->SetInputImage(moving_image);
    metric->SetMovingImage(moving_image);
    metric->SetFixedImage(fixed_image);
    metric->SetTransform(transform);
    metric->SetInterpolator(interpolator);
    metric->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion());

    if (mask_fixed_object.IsNotNull())
    {
      metric->SetFixedImageMask(mask_fixed_object);
    }
    metric->Initialize();

    ParaType last_gradient(kParaDim);
    ParaType invlast_gradient(kParaDim);
    for (int j = 0; j < kParaDim; j++)
    {
      last_gradient[j] = 0;
      invlast_gradient[j] = 0;
    }
    current_step_length = maximum_step_length;

    bool is_converged = false;
    int  used_iterations = 0;
    for (used_iterations = 0; used_iterations < number_of_iteration_current_level; used_iterations++)
    {
      ParaType original_gradient(kParaDim);
      ParaType current_gradient(kParaDim);

      value = 0;
      try
      {
        metric->GetValueAndDerivative(current_para, value, original_gradient);
      }
      catch (const itk::ExceptionObject & err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return false;
      }

      // use the similar routine as RegularStepGradientDescentBaseOptimizer::AdvanceOneStep
      // to use oscillation as the minimization convergence
      // notice this is always a minimization procedure
      if (is_rigid)
      {
        original_gradient = NormalizeGradientForRigidTransform(original_gradient, kImageDim);
      }
      for (int j = 0; j < kParaDim; j++)
      {
        current_gradient[j] = original_gradient[j] / gradient_scales[j];
      }
      double gradient_magnitude = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        gradient_magnitude += current_gradient[j] * current_gradient[j];
      }
      gradient_magnitude = sqrt(gradient_magnitude);
      double inner_product_last_current_gradient = 0.0;
      for (int j = 0; j < kParaDim; j++)
      {
        inner_product_last_current_gradient += current_gradient[j] * last_gradient[j];
      }
      if (inner_product_last_current_gradient < 0)
      {
        current_step_length *= relaxation_factor;
      }

      if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
      {
        is_converged = true;
        break;
      }

      if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
      {
        is_converged = true;
        break;
      }
      for (int j = 0; j < kParaDim; j++)
      {
        current_para[j] += (-1.0) * current_gradient[j] * current_step_length / gradient_magnitude;
      }
      if (kImageDim == 3) // normalize quaternion
      {
        double quat_mag = 0.0;
        for (int j = 0; j < 4; j++)
        {
          quat_mag += current_para[j] * current_para[j];
        }
        quat_mag = sqrt(quat_mag);
        for (int j = 0; j < 4; j++)
        {
          current_para[j] /= quat_mag;
        }
        if (!is_rigid)
        {
          for (int j = 4; j < 7; j++)
          {
            current_para[j] *= quat_mag;
          }
        }
      }

      last_gradient = current_gradient;
    }

    std::cout << "level " << i << ", iter " << used_iterations << ", size: fix"
              << fixed_image->GetRequestedRegion().GetSize() << "-mov" << moving_image->GetRequestedRegion().GetSize();

    std::cout << ", affine para: " << current_para << std::endl;

    if (is_converged)
    {
      std::cout << "    reach oscillation, current step: " << current_step_length << "<" << minimum_step_length
                << std::endl;
    }
    else
    {
      std::cout << "    does not reach oscillation, current step: " << current_step_length << ">" << minimum_step_length
                << std::endl;
    }
  }

  double value1 = value;

  if (!mask_fixed_object.IsNotNull())
  {
    typename TransformType::Pointer transform2 = TransformType::New();
    ParaType                        current_para2(kParaDim);

    // GS: should use the inverse of initial transform here, removed the next line:
    // current_para2 = opt.transform_initial->GetParameters();
    typename TransformType::Pointer transform2_initial = TransformType::New();
    opt.transform_initial->GetInverse(transform2_initial);
    current_para2 = transform2_initial->GetParameters();

    maximum_step_length = opt.maximum_step_length;
    relaxation_factor = opt.relaxation_factor;
    minimum_step_length = opt.minimum_step_length;
    value = 0;
    for (int i = 0; i < number_of_levels; i++)
    {
      transform2->SetParameters(current_para2);
      transform2->SetCenter(transform2_initial->GetCenter());

      /** see below -- we switch fixed and moving!!  */
      ImagePointerType fixed_image = moving_image_pyramid[i];
      ImagePointerType moving_image = fixed_image_pyramid[i];
      int              number_of_iteration_current_level = number_of_iteration_list[i];
      interpolator->SetInputImage(moving_image);
      metric->SetMovingImage(moving_image);
      metric->SetFixedImage(fixed_image);
      metric->SetTransform(transform2);
      metric->SetInterpolator(interpolator);
      metric->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion());

      /** FIXME --- need a moving mask ... */
      //        if (mask_fixed_object.IsNotNull()) metric->SetFixedImageMask(mask_fixed_object);
      metric->Initialize();

      ParaType last_gradient(kParaDim);
      for (int j = 0; j < kParaDim; j++)
      {
        last_gradient[j] = 0;
      }
      current_step_length = maximum_step_length;

      bool is_converged = false;
      int  used_iterations = 0;
      for (used_iterations = 0; used_iterations < number_of_iteration_current_level; used_iterations++)
      {
        ParaType original_gradient(kParaDim);
        ParaType current_gradient(kParaDim);

        value = 0;
        try
        {
          metric->GetValueAndDerivative(current_para2, value, original_gradient);
        }
        catch (const itk::ExceptionObject & err)
        {
          std::cout << "ExceptionObject caught !" << std::endl;
          std::cout << err << std::endl;
          // don't have to return here if got anything from the previous forward direction
          //          return false;
          break;
        }

        // use the similar routine as RegularStepGradientDescentBaseOptimizer::AdvanceOneStep
        // to use oscillation as the minimization convergence
        // notice this is always a minimization procedure
        if (is_rigid)
        {
          original_gradient = NormalizeGradientForRigidTransform(original_gradient, kImageDim);
        }
        for (int j = 0; j < kParaDim; j++)
        {
          current_gradient[j] = original_gradient[j] / gradient_scales[j];
        }
        double gradient_magnitude = 0.0;
        for (int j = 0; j < kParaDim; j++)
        {
          gradient_magnitude += current_gradient[j] * current_gradient[j];
        }
        gradient_magnitude = sqrt(gradient_magnitude);
        double inner_product_last_current_gradient = 0.0;
        for (int j = 0; j < kParaDim; j++)
        {
          inner_product_last_current_gradient += current_gradient[j] * last_gradient[j];
        }
        if (inner_product_last_current_gradient < 0)
        {
          current_step_length *= relaxation_factor;
        }

        if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
        {
          is_converged = true;
          break;
        }

        if (current_step_length < minimum_step_length || itk::Math::FloatAlmostEqual(gradient_magnitude, 0.0))
        {
          is_converged = true;
          break;
        }
        for (int j = 0; j < kParaDim; j++)
        {
          current_para2[j] += (-1.0) * current_gradient[j] * current_step_length / gradient_magnitude;
        }
        if (kImageDim == 3) // normalize quaternion
        {
          double quat_mag = 0.0;
          for (int j = 0; j < 4; j++)
          {
            quat_mag += current_para2[j] * current_para2[j];
          }
          quat_mag = sqrt(quat_mag);
          for (int j = 0; j < 4; j++)
          {
            current_para2[j] /= quat_mag;
          }
          if (!is_rigid)
          {
            for (int j = 4; j < 7; j++)
            {
              current_para2[j] *= quat_mag;
            }
          }
        }

        last_gradient = current_gradient;
      }

      std::cout << "level " << i << ", iter " << used_iterations << ", size: fix"
                << fixed_image->GetRequestedRegion().GetSize() << "-mov"
                << moving_image->GetRequestedRegion().GetSize();

      std::cout << ", affine para: " << current_para2 << std::endl;

      if (is_converged)
      {
        std::cout << "    reach oscillation, current step: " << current_step_length << "<" << minimum_step_length
                  << std::endl;
      }
      else
      {
        std::cout << "    does not reach oscillation, current step: " << current_step_length << ">"
                  << minimum_step_length << std::endl;
      }
    }
    std::cout << " v1 " << value1 << " v2 " << value << std::endl;
    if (value < value1)
    {
      std::cout << " last params " << transform->GetParameters() << std::endl;
      std::cout << " my params " << transform2->GetParameters() << std::endl;
      transform2->GetInverse(transform);
      std::cout << " new params " << transform->GetParameters() << std::endl;
      para_final = transform->GetParameters();
      return true;
    }
  }

  para_final = current_para;
  std::cout << "final " << para_final << std::endl;
  return true;
}

#endif /*ANTS_AFFINE_REGISTRATION2_H_*/
