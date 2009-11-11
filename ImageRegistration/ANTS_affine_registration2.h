#ifndef ANTS_AFFINE_REGISTRATION2_H_
#define ANTS_AFFINE_REGISTRATION2_H_

#include <vector>
#include <stdlib.h>
#include <time.h>
#include "itkImage.h"
#include "itkPoint.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkANTSAffine3DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include  "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"
#include "itkWarpImageWAffineFilter.h"
#include "itkImageMomentsCalculator.h"
#include <vector>
#include "ReadWriteImage.h"
#include "itkMeanSquaresImageToImageMetric.h"

typedef enum { AffineWithMutualInformation = 1, AffineWithMeanSquareDifference } AffineMetricType;

template <class TAffineTransformPointer, class TMaskImagePointer>
class OptAffine
{
public:
  typedef TAffineTransformPointer                         AffineTransformPointerType;
  typedef TMaskImagePointer                               MaskImagePointerType;
  typedef typename AffineTransformPointerType::ObjectType AffineTransformType;

  OptAffine()
  {
    MI_bins = 32;
    MI_samples = 6000;
    number_of_seeds = 0;
    time_seed = (unsigned int) time(NULL);
    number_of_levels = 3;
    number_of_iteration_list.resize(number_of_levels);
    for( int i = 0; i < number_of_levels; i++ )
      {
      number_of_iteration_list[i] = 10000;
      }
    const int kParaDim = AffineTransformType::ParametersDimension;
    gradient_scales.resize(kParaDim);
    for( int i = 0; i < kParaDim; i++ )
      {
      gradient_scales[i] = 1.0;
      }
    metric_type = AffineWithMutualInformation;
    is_rigid = false;

    maximum_step_length = 0.1;
    relaxation_factor = 0.5;
    minimum_step_length = 1.e-5;
    translation_scales = 1.e-4;

    use_rotation_header = false;
    ignore_void_orgin = true;
  };

  ~OptAffine()
  {
  };

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

template <class TAffineTransformPointer, class TMaskImagePointer>
std::ostream & operator<<(std::ostream& os, const OptAffine<TAffineTransformPointer, TMaskImagePointer>& p)
{
  typedef OptAffine<TAffineTransformPointer, TMaskImagePointer> OptAffineType;
  os << "OptAffine: ";
  os << "metric_type=";

  switch( p.metric_type )
    {
    case AffineWithMutualInformation:
      os << "AffineWithMutualInformation" << std::endl; break;
    case AffineWithMeanSquareDifference:
      os << "AffineWithMeanSquareDifference" << std::endl; break;
    }
  os << "MI_bins=" << p.MI_bins << " " << "MI_samples=" << p.MI_samples << std::endl;
  os << "number_of_seeds=" << p.number_of_seeds << " " << "time_seed=" << p.time_seed << std::endl;
  os << "number_of_levels=" << p.number_of_levels << std::endl;
  os << "number_of_iteration_list=" << "[";
  for( unsigned int i = 0; i < p.number_of_iteration_list.size() - 1; i++ )
    {
    os << p.number_of_iteration_list[i] << ",";
    }
  if( p.number_of_iteration_list.size() > 0 )
    {
    os << p.number_of_iteration_list[p.number_of_iteration_list.size() - 1];
    }
  os << "]" << std::endl;
  os << "graident_scales=" << "[";
  for( unsigned int i = 0; i < p.gradient_scales.size() - 1; i++ )
    {
    os << p.gradient_scales[i] << ",";
    }
  if( p.gradient_scales.size() > 0 )
    {
    os << p.gradient_scales[p.gradient_scales.size() - 1];
    }
  os << "]" << std::endl;
  os << "is_rigid = " << p.is_rigid << std::endl;
  os << "mask null: " << p.mask_fixed.IsNull() << std::endl;

  os << "maximum_step_length=" << p.maximum_step_length << std::endl;;
  os << "relaxation_factor=" << p.relaxation_factor << std::endl;
  os << "minimum_step_length=" << p.minimum_step_length << std::endl;
  os << "translation_scales=" << p.translation_scales << std::endl;

  return os;
};

template <class TransformPointerType, class StringType>
void WriteAffineTransformFile(TransformPointerType & transform, StringType filename)
{
  itk::TransformFileWriter::Pointer transform_writer;

  transform_writer = itk::TransformFileWriter::New();
  transform_writer->SetFileName(filename);
  transform_writer->SetInput(transform);

  try
    {
    transform_writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
              << "Exception in writing tranform file: " << std::endl
              << filename << std::endl;
    return;
    }

  return;
}

template <class StringType, class CastTransformPointerType>
void ReadAffineTransformFile(StringType filename, CastTransformPointerType & transform)
{
  typedef typename CastTransformPointerType::ObjectType CastTransformType;
//    const unsigned int InputSpaceDimension = CastTransformType::InputSpaceDimension;
//    const unsigned int OutputSpaceDimension = CastTransformType::OutputSpaceDimension;

  itk::TransformFactory<CastTransformType>::RegisterTransform();
  itk::TransformFactory<itk::ANTSAffine3DTransform<double> >::RegisterTransform();

  typedef typename itk::TransformFileReader TranReaderType;
  TranReaderType::Pointer tran_reader = TranReaderType::New();
  tran_reader->SetFileName(filename);

  try
    {
    tran_reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    std::cerr << "Exception caught in reading tran para file: "
              << filename << std::endl;
    return;
    }

  transform = dynamic_cast<CastTransformType *>( (tran_reader->GetTransformList() )->front().GetPointer() );

  return;
}

template <class OptAffine>
void InitializeAffineOptmizationParameters(OptAffine & opt, double translationScale)
{
  const int kImageDim = OptAffine::MaskImagePointerType::ObjectType::ImageDimension;

  switch( kImageDim )
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
      opt.gradient_scales[0] = 1.0;   // quaternion
      opt.gradient_scales[1] = 1.0;   // quaternion
      opt.gradient_scales[2] = 1.0;   // quaternion
      opt.gradient_scales[3] = 1.0;   // quaternion
      opt.gradient_scales[4] = 1.0;   // s1
      opt.gradient_scales[5] = 1.0;   // s2
      opt.gradient_scales[6] = 1.0;   // s3
      opt.gradient_scales[7] = 1.0;   // k1
      opt.gradient_scales[8] = 1.0;   // k2
      opt.gradient_scales[9] = 1.0;   // k3
      opt.gradient_scales[10] = translationScale;
      opt.gradient_scales[11] = translationScale;
      opt.gradient_scales[12] = translationScale;
      }
      break;
    }

  std::cout << opt;
}

template <class TMaskObjectPointerType, class TImagePyramid, class TMetricPointerType, class TInterpolatorPointerType>
class RunningAffineCache
{
public:
  typedef TImagePyramid                                ImagePyramidType;
  typedef typename ImagePyramidType::value_type        ImagePointerType;
  typedef typename ImagePointerType::ObjectType        ImageType;
  typedef TMetricPointerType                           MetricPointerType;
  typedef typename MetricPointerType::ObjectType       MetricType;
  typedef TInterpolatorPointerType                     InterpolatorPointerType;
  typedef typename InterpolatorPointerType::ObjectType InterpolatorType;
  typedef TMaskObjectPointerType                       MaskObjectPointerType;

  RunningAffineCache()
  {
  };
  ~RunningAffineCache()
  {
  };

  TMaskObjectPointerType   mask_fixed_object;
  TImagePyramid            fixed_image_pyramid;
  TImagePyramid            moving_image_pyramid;
  TMetricPointerType       metric;
  TInterpolatorPointerType interpolator;
};

template <class ImagePointerType, class TransformPointerType, class OptAffineType>
void ComputeSingleAffineTransform(ImagePointerType fixedImage, ImagePointerType movingImage, OptAffineType & opt,
                                  TransformPointerType & transform);

template <class ImagePointerType, class TransformPointerType, class OptAffineType>
void ComputeSingleAffineTransform2D3D(ImagePointerType I_fixed, ImagePointerType I_moving, OptAffineType & opt,
                                      TransformPointerType & transform);

// //////////////////////////////////////////////////////////////////////
template <class ImagePointerType, class RunningImagePointerType, class OptAffineType, class RunningOptAffineType>
inline void PreConversionInAffine(ImagePointerType & fixedImage, RunningImagePointerType& R_fixedImage,
                                  ImagePointerType & movingImage, RunningImagePointerType& R_movingImage,
                                  OptAffineType & opt, RunningOptAffineType & R_opt)
{
  typedef typename OptAffineType::AffineTransformPointerType::ObjectType        AffineTransformType;
  typedef typename RunningOptAffineType::AffineTransformPointerType::ObjectType RunningAffineTransformType;

  if( opt.use_rotation_header )
    {
    std::cout << "===================>initialize from rotation header ... " << std::endl;
    // use the rotation header to initialize the affine: inv(Tm) * Tf
    typename AffineTransformType::Pointer aff_Im = AffineTransformType::New();
    GetAffineTransformFromImage(movingImage, aff_Im);
    typename AffineTransformType::Pointer aff_If = AffineTransformType::New();
    GetAffineTransformFromImage(fixedImage, aff_If);
    typename AffineTransformType::Pointer aff_combined = AffineTransformType::New();
    aff_combined->SetFixedParameters(aff_If->GetFixedParameters() );
    aff_combined->SetParameters(aff_If->GetParameters() );
    typename AffineTransformType::Pointer aff_Im_inv = AffineTransformType::New();
    aff_Im->GetInverse(aff_Im_inv);
    aff_combined->Compose(aff_Im_inv, 0);
    opt.transform_initial = aff_combined;

//            std::cout << "aff_If: " << aff_If << std::endl;
//            std::cout << "aff_Im: " << aff_Im << std::endl;
//            std::cout << "aff_combined: " << aff_combined << std::endl;
    }

  if( !opt.use_rotation_header && opt.ignore_void_orgin )
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

    if( b_far_origin_without_rotation )
      {
      typename AffineTransformType::Pointer aff_combined = AffineTransformType::New();
      aff_combined->SetFixedParameters(aff_If->GetFixedParameters() );
      aff_combined->SetParameters(aff_If->GetParameters() );
      typename AffineTransformType::Pointer aff_Im_inv = AffineTransformType::New();
      aff_Im->GetInverse(aff_Im_inv);
      aff_combined->Compose(aff_Im_inv, 0);
      opt.transform_initial = aff_combined;
      }
    }

  if( opt.transform_initial.IsNotNull() )
    {
    R_opt.transform_initial = RunningAffineTransformType::New();

    R_opt.transform_initial->SetCenter(*(reinterpret_cast<typename RunningAffineTransformType::InputPointType *>
                                         (const_cast<typename AffineTransformType::InputPointType *>(&(opt.
                                                                                                       transform_initial
                                                                                                       ->GetCenter() ) ) ) ) );
    R_opt.transform_initial->SetMatrix(*(reinterpret_cast<typename RunningAffineTransformType::MatrixType *>
                                         (const_cast<typename AffineTransformType::MatrixType *>(&(opt.
                                                                                                   transform_initial->
                                                                                                   GetMatrix() ) ) ) ) );
    R_opt.transform_initial->SetTranslation(*(reinterpret_cast<typename RunningAffineTransformType::OutputVectorType *>
                                              (const_cast<typename AffineTransformType::OutputVectorType *>(&(opt.
                                                                                                              transform_initial
                                                                                                              ->
                                                                                                              GetTranslation() ) ) ) ) );
    }

  // std::cout << "R_opt.transform_initial" << R_opt.transform_initial << std::endl;

  if( opt.mask_fixed.IsNotNull() )
    {
    R_opt.mask_fixed = RunningOptAffineType::MaskImagePointerType::ObjectType::New();
    R_opt.mask_fixed = reinterpret_cast<typename RunningOptAffineType::MaskImagePointerType &>(opt.mask_fixed);
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
template <class RunningAffineTransformPointerType, class AffineTransformPointerType>
inline void PostConversionInAffine(RunningAffineTransformPointerType& transform_running,
                                   AffineTransformPointerType & transform)
{
  typedef typename RunningAffineTransformPointerType::ObjectType RunningAffineTransformType;
  typedef typename AffineTransformPointerType::ObjectType        AffineTransformType;

  transform->SetCenter(*(reinterpret_cast<typename AffineTransformType::InputPointType *>
                         (const_cast<typename RunningAffineTransformType::InputPointType *>(&(transform_running->
                                                                                              GetCenter() ) ) ) ) );
  transform->SetTranslation(*(reinterpret_cast<typename AffineTransformType::OutputVectorType *>
                              (const_cast<typename RunningAffineTransformType::OutputVectorType *>(&(transform_running
                                                                                                     ->GetTranslation() ) ) ) ) );
  transform->SetMatrix(*(reinterpret_cast<typename AffineTransformType::MatrixType *>
                         (const_cast<typename RunningAffineTransformType::MatrixType *>(&(transform_running->GetMatrix() ) ) ) ) );

  // std::cout << "transform_running" << transform_running << std::endl;
  // std::cout << "transform" << transform << std::endl;
}

// /////////////////////////////////////////////////////////////////////////
// the initial transform maybe any derivative class type from MatrixOffsetTransformBase,
// it will be automatically converted to the my 2D/3D affine type
template <class ImagePointerType, class TransformPointerType, class OptAffineType>
void ComputeSingleAffineTransform(ImagePointerType fixedImage, ImagePointerType movingImage, OptAffineType & opt,
                                  TransformPointerType & transform)
{
  typedef typename ImagePointerType::ObjectType ImageType;
  const int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::IOPixelType PixelType;

  std::cout << "transform_initial: IsNotNull():" << opt.transform_initial.IsNotNull() << std::endl;

  if( ImageDimension == 2 )
    {
    typedef itk::ANTSCenteredAffine2DTransform<double>::Pointer RunningAffineTransformPointerType;
    const unsigned int RunningImageDimension = 2;

    typedef typename itk::Image<PixelType, RunningImageDimension>::Pointer        RunningImagePointerType;
    typedef OptAffine<RunningAffineTransformPointerType, RunningImagePointerType> RunningOptAffineType;
    RunningImagePointerType           R_fixedImage, R_movingImage;
    RunningOptAffineType              R_opt;
    RunningAffineTransformPointerType transform_running;

    PreConversionInAffine(fixedImage, R_fixedImage, movingImage, R_movingImage, opt, R_opt);

    ComputeSingleAffineTransform2D3D(R_fixedImage, R_movingImage, R_opt, transform_running);

    PostConversionInAffine(transform_running, transform);
    }
  else if( ImageDimension == 3 )
    {
    typedef itk::ANTSAffine3DTransform<double>::Pointer RunningAffineTransformPointerType;
    const unsigned int RunningImageDimension = 3;

    typedef typename itk::Image<PixelType, RunningImageDimension>::Pointer        RunningImagePointerType;
    typedef OptAffine<RunningAffineTransformPointerType, RunningImagePointerType> RunningOptAffineType;

    RunningImagePointerType           R_fixedImage, R_movingImage;
    RunningOptAffineType              R_opt;
    RunningAffineTransformPointerType transform_running;

    PreConversionInAffine(fixedImage, R_fixedImage, movingImage, R_movingImage, opt, R_opt);

    ComputeSingleAffineTransform2D3D(R_fixedImage, R_movingImage, R_opt, transform_running);

    PostConversionInAffine(transform_running, transform);
    }
  else
    {
    std::cout << "Unsupported, not 2D/ 3D" << std::endl;
    return;
    }
}

// /////////////////////////////////////////////////////////////////////////////
template <class MaskImagePointerType, class ImageMaskSpatialObjectPointerType>
void InitialzeImageMask(MaskImagePointerType & mask_fixed, ImageMaskSpatialObjectPointerType & mask_fixed_object)
{
  if( mask_fixed.IsNull() )
    {
    return;
    }

  const unsigned int ImageDimension = MaskImagePointerType::ObjectType::ImageDimension;
  typedef typename MaskImagePointerType::ObjectType              MaskImageType;
  typedef typename ImageMaskSpatialObjectPointerType::ObjectType ImageMaskSpatialObjectType;

  typedef itk::Image<unsigned char, ImageDimension>              CharMaskImageType;
  typedef itk::CastImageFilter<MaskImageType, CharMaskImageType> CastFilterType;
  typename CastFilterType::Pointer cast_filter = CastFilterType::New();
  cast_filter->SetInput(mask_fixed);
  cast_filter->Update();
  typename CharMaskImageType::Pointer mask_fixed_char = cast_filter->GetOutput();

  mask_fixed_object = ImageMaskSpatialObjectType::New();
  mask_fixed_object->SetImage(mask_fixed_char);
}

// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
template <class ImagePointerType, class PointType, class VectorType>
void ComputeInitialPosition(ImagePointerType & I_fixed, ImagePointerType & I_moving, PointType & center,
                            VectorType & translation_vec)
{
  typedef typename ImagePointerType::ObjectType           ImageType;
  typedef typename itk::ImageMomentsCalculator<ImageType> ImageCalculatorType;

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typename ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();

  calculator->SetImage(  I_fixed );
  calculator->Compute();
  typename ImageCalculatorType::VectorType fixed_center = calculator->GetCenterOfGravity();

  calculator->SetImage(  I_moving );
  calculator->Compute();
  typename ImageCalculatorType::VectorType moving_center = calculator->GetCenterOfGravity();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    center[i] = fixed_center[i];
    translation_vec[i] = moving_center[i] - fixed_center[i];
    }
}

// /////////////////////////////////////////////////////////////////////////////
template <class PointType, class VectorType, class TransformPointerType>
void InjectInitialPara(PointType & center, VectorType & translation_vec, TransformPointerType & transform)
{
  typedef typename TransformPointerType::ObjectType::ParametersType ParaType;
  ParaType para0(TransformPointerType::ObjectType::ParametersDimension);

  switch( (unsigned int) PointType::PointDimension )
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
      para0[0] = 0.0; para0[1] = 0.0;      para0[2] = 0.0;        para0[3] = 1.0;
      para0[4] = 1.0; para0[5] = 1.0;   para0[6] = 1.0;
      para0[7] = 0.0;  para0[8] = 0.0;     para0[9] = 0.0;
      para0[10] = translation_vec[0]; para0[11] = translation_vec[1];   para0[12] = translation_vec[2];
      // para0[10] = 0.0; para0[11] = 0.0;   para0[12] = 0.0;

      transform->SetParameters(para0);
      transform->SetCenter(center);

      break;
    }
}

// ////////////////////////////////////////////////////////////////////////////////////////
template <class ImagePointerType, class ParaType, class PointType, class TransformTypePointer>
double TestCostValueMMI(ImagePointerType fixedImage, ImagePointerType movingImage, ParaType para, PointType center,
                        TransformTypePointer null_transform)
{
  typedef typename ImagePointerType::ObjectType     ImageType;
  typedef typename TransformTypePointer::ObjectType TransformType;

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetCenter(center);
  // transform->SetParameters(para);

  typedef typename itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> mattesMutualInfoMetricType;
  typename mattesMutualInfoMetricType::Pointer mattesMutualInfo = mattesMutualInfoMetricType::New();

  typedef typename itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputImage(movingImage);

  mattesMutualInfo->SetFixedImage(fixedImage);
  mattesMutualInfo->SetMovingImage(movingImage);
  mattesMutualInfo->SetFixedImageRegion(fixedImage->GetBufferedRegion() );
  mattesMutualInfo->SetTransform(transform);
  mattesMutualInfo->SetInterpolator(interpolator);
  mattesMutualInfo->SetNumberOfHistogramBins( 32 );
  mattesMutualInfo->SetNumberOfSpatialSamples( 5000 );
  mattesMutualInfo->SetTransformParameters(para);
  mattesMutualInfo->Initialize();
  double rval = 0;
  try
    {
    rval = mattesMutualInfo->GetValue(para);
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
              << "Exception caught in computing mattesMutualInfo after registration" << std::endl
              << "Maybe: Too many samples map outside moving image buffer" << std::endl
              << "Set the cost value = 0 (max for MutualInfo) " << std::endl;
    rval = 0;
    }

  return rval;
}

template <class ImagePointerType, class OptAffineType, class RunningAffineCacheType>
void InitializeRunningAffineCache(ImagePointerType & fixed_image, ImagePointerType & moving_image, OptAffineType & opt,
                                  RunningAffineCacheType & running_cache)
{
  typedef typename ImagePointerType::ObjectType ImageType;
  // typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typedef typename RunningAffineCacheType::InterpolatorType InterpolatorType;
  typedef typename RunningAffineCacheType::MetricType       MetricType;
  // typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> MetricType;

  BuildImagePyramid(fixed_image, opt.number_of_levels, running_cache.fixed_image_pyramid);
  BuildImagePyramid(moving_image, opt.number_of_levels, running_cache.moving_image_pyramid);
  InitialzeImageMask(opt.mask_fixed, running_cache.mask_fixed_object);

  running_cache.interpolator = InterpolatorType::New();
  running_cache.metric = MetricType::New();
}

template <class ImageTypePointer, class AffineTransformPointer>
void GetAffineTransformFromImage(const ImageTypePointer& img, AffineTransformPointer & aff)
{
  typedef typename ImageTypePointer::ObjectType                        ImageType;
  typedef typename ImageType::DirectionType                            DirectionType;
  typedef typename ImageType::PointType                                PointType;
  typedef typename ImageType::SpacingType                              SpacingType;
  typedef typename AffineTransformPointer::ObjectType::TranslationType VectorType;

  DirectionType direction = img->GetDirection();
  PointType     pt = img->GetOrigin();
  SpacingType   spacing = img->GetSpacing();
  VectorType    translation;
  translation.Fill(0);
  aff->SetMatrix(direction);
  aff->SetCenter(pt);
  aff->SetTranslation(translation);
}

template <class ImagePointerType, class OptAffineType>
void  InitializeAffineTransform(ImagePointerType & fixed_image, ImagePointerType & moving_image, OptAffineType& opt)
{
  typedef typename OptAffineType::AffineTransformType TransformType;
  typedef typename TransformType::ParametersType      ParaType;
  typedef typename TransformType::InputPointType      PointType;
  typedef typename TransformType::OutputVectorType    VectorType;

  std::cout << "opt.transform_initial.IsNull(): " << opt.transform_initial.IsNull() << std::endl;
  std::cout << " opt.use_rotation_header: " << opt.use_rotation_header << std::endl;
  std::cout << " opt.ignore_void_orgin: " << opt.ignore_void_orgin << std::endl;

  if( opt.transform_initial.IsNull() )
    {
    PointType  center;
    VectorType translation_vec;
    ComputeInitialPosition(fixed_image, moving_image, center, translation_vec);
    opt.transform_initial = TransformType::New();
    InjectInitialPara(center, translation_vec, opt.transform_initial);
    }
}

template <class ImagePointer>
ImagePointer  ShrinkImageToScale(ImagePointer image,  float scalingFactor )
{
  typedef float RealType;

  typedef typename ImagePointer::ObjectType ImageType;
  typename ImageType::SpacingType inputSpacing = image->GetSpacing();
  typename ImageType::RegionType::SizeType inputSize = image->GetRequestedRegion().GetSize();

  typename ImageType::SpacingType outputSpacing;
  typename ImageType::RegionType::SizeType outputSize;

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  RealType minimumSpacing = inputSpacing.GetVnlVector().min_value();
//    RealType maximumSpacing = inputSpacing.GetVnlVector().max_value();

  ImagePointer current_image = image;
  for( unsigned int d = 0; d < ImageType::ImageDimension; d++ )
    {
    RealType scaling = vnl_math_min( scalingFactor * minimumSpacing / inputSpacing[d],
                                     static_cast<RealType>( inputSize[d] ) / 32.0 );
    outputSpacing[d] = inputSpacing[d] * scaling;
    outputSize[d] = static_cast<unsigned long>( inputSpacing[d]
                                                * static_cast<RealType>( inputSize[d] ) / outputSpacing[d] + 0.5 );

    typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
    typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
    smoother->SetInputImage( current_image  );
    smoother->SetDirection( d );
    smoother->SetNormalizeAcrossScale( false );
    smoother->SetSigma( 0.25 * ( outputSpacing[d] / inputSpacing[d]  ) );

    if( smoother->GetSigma() > 0.0 )
      {
      smoother->Update();
      current_image  = smoother->GetOutput();
      }
    }

  resampler->SetInput(current_image );
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(image->GetOrigin() );
  resampler->SetOutputDirection(image->GetDirection() );

  resampler->Update();

  image = resampler->GetOutput();

  // std::cout << "DEBUG: " << outputSize << std::endl;

  return image;
}

template <class ImagePointerType, class ImagePyramidType>
void BuildImagePyramid(const ImagePointerType & image, int number_of_levels, ImagePyramidType & image_pyramid)
{
  image_pyramid.resize(number_of_levels);

  image_pyramid[number_of_levels - 1] = image;
  double scale_factor = 2;
  for( int i = 0; i < number_of_levels - 1; i++ )
    {
    image_pyramid[number_of_levels - 2 - i] = ShrinkImageToScale(image, scale_factor);
    scale_factor *= 2;
    }

//    for(int i=0; i < number_of_levels; i++)
//      std::cout << "level " << i << ": size: " << image_pyramid[i]->GetLargestPossibleRegion().GetSize() << std::endl;
}

template <class ParaType>
ParaType NormalizeGradientForRigidTransform(ParaType & original_gradient, int kImageDim)
{
  ParaType new_gradient(original_gradient.Size() );

  new_gradient = original_gradient;

  switch( kImageDim )
    {
    case 2: // theta, s1, s2, k
      for( int j = 1; j <= 3; j++ )
        {
        new_gradient[j] = 0.;
        }
      break;
    case 3: // q1,q2,q3,q4,s1,s2,s3,k1,k2,k3
      for( int j = 4; j <= 9; j++ )
        {
        new_gradient[j] = 0.;
        }
      break;
    }
  return new_gradient;
}

// /////////////////////////////////////////////////////////////////////////////
// template<class ImagePointerType, class ImageMaskSpatialObjectPointerType, class ParaType>
template <class RunningAffineCacheType, class OptAffine, class ParaType>
bool RegisterImageAffineMutualInformationMultiResolution(RunningAffineCacheType & running_cache, OptAffine & opt,
                                                         ParaType & para_final)
{
  typedef typename RunningAffineCacheType::ImagePyramidType      ImagePyramidType;
  typedef typename RunningAffineCacheType::ImagePointerType      ImagePointerType;
  typedef typename ImagePointerType::ObjectType                  ImageType;
  typedef typename RunningAffineCacheType::MetricType            MetricType;
  typedef typename RunningAffineCacheType::InterpolatorType      InterpolatorType;
  typedef typename OptAffine::AffineTransformType                TransformType;
  typedef typename RunningAffineCacheType::MaskObjectPointerType MaskObjectPointerType;

  const unsigned int kImageDim = ImageType::ImageDimension;

  ImagePyramidType&      fixed_image_pyramid = running_cache.fixed_image_pyramid;
  ImagePyramidType&      moving_image_pyramid = running_cache.moving_image_pyramid;
  MaskObjectPointerType& mask_fixed_object = running_cache.mask_fixed_object;

  int                  number_of_levels = opt.number_of_levels;
  std::vector<int>&    number_of_iteration_list = opt.number_of_iteration_list;
  std::vector<double>& gradient_scales = opt.gradient_scales;
  bool                 is_rigid = opt.is_rigid;

  // use my own's registration routine of image pyramid and gradient descent , only use ITK's implementation of mutual
  // information
  // try to use my own transform class together with image pyramid when transform are required (not much for MI though)

  typename TransformType::Pointer transform = TransformType::New();

  typename InterpolatorType::Pointer& interpolator = running_cache.interpolator;
  typename MetricType::Pointer& metric = running_cache.metric;

  const int kParaDim = TransformType::ParametersDimension;

  ParaType current_para(kParaDim);
  current_para = opt.transform_initial->GetParameters();

  double maximum_step_length = opt.maximum_step_length;
  double relaxation_factor = opt.relaxation_factor;
  double minimum_step_length = opt.minimum_step_length;
  double current_step_length;
  for( int i = 0; i < number_of_levels; i++ )
    {
    transform->SetParameters(current_para);
    transform->SetCenter(opt.transform_initial->GetCenter() );

    ImagePointerType fixed_image = fixed_image_pyramid[i];
    ImagePointerType moving_image = moving_image_pyramid[i];
    int              number_of_iteration_current_level = number_of_iteration_list[i];
    interpolator->SetInputImage( moving_image );
    metric->SetMovingImage( moving_image );
    metric->SetFixedImage( fixed_image );
    metric->SetTransform( transform );
    metric->SetInterpolator( interpolator );
    metric->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion() );

    if( mask_fixed_object.IsNotNull() )
      {
      metric->SetFixedImageMask(mask_fixed_object);
      }
    metric->Initialize();

    ParaType last_gradient(kParaDim);
    for( int j = 0; j < kParaDim; j++ )
      {
      last_gradient[j] = 0;
      }
    current_step_length = maximum_step_length;

    bool is_converged = false;
    int  used_iterations = 0;
    for( used_iterations = 0; used_iterations < number_of_iteration_current_level; used_iterations++ )
      {
      ParaType original_gradient(kParaDim);
      ParaType current_gradient(kParaDim);

      double value;
      try
        {
        metric->GetValueAndDerivative(current_para, value, original_gradient);
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return false;
        break;
        }

      // use the similar routine as RegularStepGradientDescentBaseOptimizer::AdvanceOneStep
      // to use oscillation as the minimization convergence
      // notice this is always a minimization procedure
      if( is_rigid )
        {
        original_gradient = NormalizeGradientForRigidTransform(original_gradient, kImageDim);
        }
      for( int j = 0; j < kParaDim; j++ )
        {
        current_gradient[j] = original_gradient[j] / gradient_scales[j];
        }
      double gradient_magnitude = 0.0;
      for( int j = 0; j < kParaDim; j++ )
        {
        gradient_magnitude += current_gradient[j] * current_gradient[j];
        }
      gradient_magnitude = sqrt(gradient_magnitude);
      double inner_product_last_current_gradient = 0.0;
      for( int j = 0; j < kParaDim; j++ )
        {
        inner_product_last_current_gradient += current_gradient[j] * last_gradient[j];
        }
      if( inner_product_last_current_gradient < 0 )
        {
        current_step_length *= relaxation_factor;
        }
      if( current_step_length < minimum_step_length || gradient_magnitude == 0.0 )
        {
        is_converged = true;
        break;
        }
//            for(int j = 0; j < kParaDim; j++)  std::cerr<< current_gradient[i] << ", ";
//            std::cerr<< std::endl;
      for( int j = 0; j < kParaDim; j++ )
        {
        current_para[j] += (-1.0) * current_gradient[j] * current_step_length / gradient_magnitude;
        }

      if( kImageDim == 3 )     // normalize quaternion
        {
        double quat_mag = 0.0;
        for( int j = 0; j < 4; j++ )
          {
          quat_mag += current_para[j] * current_para[j];
          }
        quat_mag = sqrt(quat_mag);
        for( int j = 0; j < 4; j++ )
          {
          current_para[j] /= quat_mag;
          }
        if( !is_rigid )
          {
          for( int j = 4; j < 7; j++ )
            {
            current_para[j] *= quat_mag;
            }
          }
        }

//            if (used_iterations%100 == 1){
//                std::cout << "["<<used_iterations<<"]:" << value <<  ", current_para=" << current_para << std::endl;
//                std::cout << "["<<used_iterations<<"]:" << "step_length=" << current_step_length <<  ", factor="
//                    << current_step_length / gradient_magnitude <<  " original_gradient=" << original_gradient <<
// std::endl;
//            }

      last_gradient = current_gradient;
      }

    std::cout << "level " << i << ", iter " << used_iterations
              << ", size: fix" << fixed_image->GetRequestedRegion().GetSize()
              << "-mov" << moving_image->GetRequestedRegion().GetSize();

    std::cout << ", affine para: " << current_para << std::endl;

    if( is_converged )
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
  para_final = current_para;
  return true;
}

// /////////////////////////////////////////////////////////////////////////////
template <class ImagePointerType, class TransformPointerType, class OptAffineType>
void ComputeSingleAffineTransform2D3D(ImagePointerType fixed_image, ImagePointerType moving_image, OptAffineType & opt,
                                      TransformPointerType & transform)
{
  typedef typename ImagePointerType::ObjectType ImageType;
  const int ImageDimension = ImageType::ImageDimension;
  typedef std::vector<ImagePointerType>                          ImagePyramidType;
  typedef itk::ImageMaskSpatialObject<ImageDimension>            ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::Pointer           MaskObjectPointerType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typedef typename InterpolatorType::Pointer                     InterpolatorPointerType;

  typedef typename TransformPointerType::ObjectType TransformType;
  typedef typename TransformType::ParametersType    ParaType;

  InitializeAffineOptmizationParameters(opt, opt.translation_scales);

  // std::cout << "DEBUG: opt.gradient_scales.size() = " << opt.gradient_scales.size() << std::endl;

  InitializeAffineTransform(fixed_image, moving_image, opt);

  std::cout << "input affine center: " << opt.transform_initial->GetCenter() << std::endl;
  std::cout << "input affine para: " << opt.transform_initial->GetParameters() << std::endl;

  transform = TransformType::New();
  ParaType para_final(TransformType::ParametersDimension);

  switch( opt.metric_type )
    {
    case AffineWithMeanSquareDifference:
      {
      typedef itk::MeanSquaresImageToImageMetric<ImageType,
                                                 ImageType> MetricType;
      typedef typename MetricType::Pointer
                                                            MetricPointerType;
      typedef RunningAffineCache<MaskObjectPointerType, ImagePyramidType, MetricPointerType,
                                 InterpolatorPointerType>   RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);
      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
      }
      break;
    case AffineWithMutualInformation:
      {
      typedef itk::MattesMutualInformationImageToImageMetric<ImageType,
                                                             ImageType> MetricType;
      typedef typename MetricType::Pointer
                                                                        MetricPointerType;
      typedef RunningAffineCache<MaskObjectPointerType, ImagePyramidType, MetricPointerType,
                                 InterpolatorPointerType>               RunningAffineCacheType;

      RunningAffineCacheType running_cache;
      InitializeRunningAffineCache(fixed_image, moving_image, opt, running_cache);

      running_cache.metric->SetNumberOfHistogramBins( opt.MI_bins );
      running_cache.metric->SetNumberOfSpatialSamples( opt.MI_samples );

      RegisterImageAffineMutualInformationMultiResolution(running_cache, opt, para_final);
      }
      break;
    default:
      break;
    }

  bool noaffine = true;
  for( int i = 0; i < opt.number_of_levels; i++ )
    {
    if( opt.number_of_iteration_list[i] > 0 )
      {
      noaffine = false;
      }
    }
  if( noaffine )
    {
    for( int i = TransformType::ParametersDimension - ImageDimension;  i < TransformType::ParametersDimension; i++ )
      {
      para_final[i] = 0;
      }
    }

  transform->SetParameters(para_final);
  transform->SetCenter(opt.transform_initial->GetCenter() );

  double rval_init = TestCostValueMMI(fixed_image, moving_image,
                                      opt.transform_initial->GetParameters(),
                                      opt.transform_initial->GetCenter(), transform);

  // std::cerr << "ABCDABCD: " << transform << std::endl;

  double rval_final = TestCostValueMMI(fixed_image, moving_image, para_final,
                                       opt.transform_initial->GetCenter(), transform);

  std::cout << "outputput affine center: " << transform->GetCenter() << std::endl;
  std::cout << "output affine para: " << transform->GetParameters() << std::endl;
  std::cout << "initial measure value (MMI): rval = " << rval_init << std::endl;
  std::cout << "final measure value (MMI): rval = " << rval_final << std::endl;
  std::cout << "finish affine registeration..."  <<  std::endl;
};

#endif /*ANTS_AFFINE_REGISTRATION2_H_*/
