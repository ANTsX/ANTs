#ifndef ANTS_AFFINE_REGISTRATION_H_
#define ANTS_AFFINE_REGISTRATION_H_

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

template <typename StringType, typename CastTransformPointerType>
void
read_transform_file(StringType filename, CastTransformPointerType & transform);

template <typename TransformPointerType, typename StringType>
void
write_transform_file(TransformPointerType & transform, StringType str);

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform_3d(ImagePointerType       I_fixed,
                                   ImagePointerType       I_moving,
                                   MaskImagePointerType   mask_fixed,
                                   TransformPointerType & transform,
                                   TransformPointerType & transform_initial);

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform_2d(ImagePointerType       I_fixed,
                                   ImagePointerType       I_moving,
                                   MaskImagePointerType   mask_fixed,
                                   TransformPointerType & transform,
                                   TransformPointerType & transform_initial);

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform(ImagePointerType       fixedImage,
                                ImagePointerType       movingImage,
                                MaskImagePointerType   maskImage,
                                TransformPointerType & transform,
                                TransformPointerType & transform_initial);

// void compute_single_affine_transform_2d(itk::Image<float, 2>::Pointer I_fixed,
//        itk::Image<float, 2>::Pointer I_moving,
//        itk::Image<float, 2>::Pointer mask_fixed,
//        itk::CenteredAffine2DTransform<double>::Pointer &transform);
template <typename DisplacementFieldPointerType>
void
create_deformation_field_byref(const DisplacementFieldPointerType & ref, DisplacementFieldPointerType & field);

// this is obsolet, use itkWarpImageWAffineFilter
template <typename TransformPointerType, typename DisplacementFieldPointerType>
void
compose_affine_with_field(const TransformPointerType &         aff,
                          const DisplacementFieldPointerType & field,
                          DisplacementFieldPointerType &       field_output);

template <typename ImagePointerType, typename DisplacementFieldPointerType>
void
warp_image_field(const ImagePointerType &             img_input,
                 const DisplacementFieldPointerType & field,
                 ImagePointerType &                   img_output);

template <typename ImagePointerType, typename TransformPointerType, typename DisplacementFieldPointerType>
void
warp_image_field_waffine(const ImagePointerType &             img_input,
                         const TransformPointerType &         aff,
                         const DisplacementFieldPointerType & field,
                         ImagePointerType &                   img_output);

template <typename ImageTypePointer, typename RefImageTypePointer, typename TransformTypePointer>
void
affine_image(const ImageTypePointer &     input_image,
             const TransformTypePointer & transform,
             const RefImageTypePointer &  ref_image,
             ImageTypePointer &           img_aff);

template <typename TransformPointerType, typename StringType>
void
write_transform_file(TransformPointerType & transform, StringType filename)
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
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
              << "Exception in writing transform file: " << std::endl
              << filename << std::endl;
    return;
  }

  return;
}

template <typename StringType, typename CastTransformPointerType>
void
read_transform_file(StringType filename, CastTransformPointerType & transform)
{
  typedef typename CastTransformPointerType::ObjectType CastTransformType;
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

template <typename TransformPointerType, typename DisplacementFieldType>
void
convert_affine_para_to_deformation_field(TransformPointerType & transform, DisplacementFieldType & def)
{
  return;
}

template <typename ParaType, unsigned Dimension = 3>
typename SEARCH_POINT_TYPE
{
public:
  ParaType                      para0;               // seed parameter
  ParaType                      para1;               // local optimal parameter
  itk::Point<double, Dimension> center;              // transformation center
  double                        rval;                // registered value
  unsigned int                  index;               // the index in the list, //for sorting
  int                           number_of_iteration; // the number of iteration used for gradient descent
};

template <typename SEARCH_LIST, typename ParaType>
void
generate_search_seed_3d(SEARCH_LIST & search_list, ParaType & ret);

template <typename SEARCH_LIST, typename ParaType>
void
generate_search_seed_2d(SEARCH_LIST & search_list, ParaType & ret);

template <typename SEARCH_LIST, typename SEARCH_POINT>
void
get_best_search_point(SEARCH_LIST & search_list, SEARCH_POINT & spt);

template <typename SEARCH_LIST, typename SEARCH_POINT>
void
add_search_point(SEARCH_LIST & search_list, SEARCH_POINT spt);

// use moment of the image  initializer to get cx and cy only
template <typename ImagePointerType, typename ParaType>
bool
register_image_cxyz(ImagePointerType fixed_image, ImagePointerType moving_image, ParaType & para1, double & rval);

// use an optional mask for the fixed image
template <typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
bool
register_image_affine3d_mres_mask(ImagePointerType                  fixed_image,
                                  ImagePointerType                  moving_image,
                                  ImageMaskSpatialObjectPointerType mask_fixed_object,
                                  ParaType                          para0,
                                  itk::Point<double, 3>             center,
                                  int                               number_of_iteration,
                                  int                               MI_bins,
                                  int                               MI_samples,
                                  ParaType &                        para1,
                                  double &                          rval);

template <typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
bool
register_image_affine2d_mres_mask(ImagePointerType                  fixed_image,
                                  ImagePointerType                  moving_image,
                                  ImageMaskSpatialObjectPointerType mask_fixed_object,
                                  ParaType                          para0,
                                  itk::Point<double, 2>             center,
                                  int                               number_of_iteration,
                                  int                               MI_bins,
                                  int                               MI_samples,
                                  ParaType &                        para1,
                                  double &                          rval);

template <typename ImagePointerType, typename ParaType, typename PointType, typename TransformTypePointer>
double
get_cost_value_mmi(ImagePointerType     fixedImage,
                   ImagePointerType     movingImage,
                   ParaType             para,
                   PointType            center,
                   TransformTypePointer null_transform);

template <typename SEARCH_POINT>
bool
compare_search_point(const SEARCH_POINT & lhs, const SEARCH_POINT & rhs)
{
  return lhs.rval < rhs.rval;
}

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform_3d(ImagePointerType       I_fixed,
                                   ImagePointerType       I_moving,
                                   MaskImagePointerType   mask_fixed,
                                   TransformPointerType & transform,
                                   TransformPointerType & transform_initial)
{
  typedef typename ImagePointerType::ObjectType     ImageType;
  const int                                         ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::IOPixelType           PixelType;
  typedef typename MaskImagePointerType::ObjectType MaskImageType;
  typedef typename MaskImageType::IOPixelType       MaskPixelType;

  typedef typename TransformPointerType::ObjectType TransformType;
  typedef typename TransformType::ParametersType    ParaType;

  // option parameters
  int          number_of_seeds = 0;
  int          number_of_iteration = 10000;
  int          MI_bins = 32;
  int          MI_samples = 6000;
  unsigned int time_seed = (unsigned int)time(nullptr);
  srand(time_seed);
  // TODO: need to fix here
  bool b_use_mask = 0; // (mask_fixed == nullptr);

  std::cout << "number_of_seeds: " << number_of_seeds << std::endl;
  std::cout << "rand_time_seed: " << time_seed << std::endl;
  std::cout << "number_of_iteration: " << number_of_iteration << std::endl;
  std::cout << "MI_bins: " << MI_bins << std::endl;
  std::cout << "MI_samples: " << MI_samples << std::endl;
  std::cout << "use mask: " << b_use_mask << std::endl;

  // memory of searched results
  typedef SEARCH_POINT_TYPE<ParaType, ImageDimension> SEARCH_POINT;
  typedef std::vector<SEARCH_POINT>                   SEARCH_LIST;

  SEARCH_LIST search_list;

  // typedef itk::ImageMaskSpatialObject<ImageDimension, MaskPixelType> ImageMaskSpatialObject;
  typedef itk::ImageMaskSpatialObject<ImageDimension> ImageMaskSpatialObject;
  typename ImageMaskSpatialObject::Pointer            mask_fixed_object = 0;
  if (b_use_mask)
  {
    typedef itk::Image<unsigned char, ImageDimension>              CharMaskImageType;
    typedef itk::CastImageFilter<MaskImageType, CharMaskImageType> CastFilterType;
    typename CastFilterType::Pointer                               cast_filter = CastFilterType::New();
    cast_filter->SetInput(mask_fixed);
    cast_filter->Update();
    typename CharMaskImageType::Pointer mask_fixed_char = cast_filter->GetOutput();

    mask_fixed_object = ImageMaskSpatialObject::New();
    mask_fixed_object->SetImage(mask_fixed_char);
  }

  typename ImageType::PointType origin;
  origin.Fill(0);
  I_moving->SetOrigin(origin);
  I_fixed->SetOrigin(origin);

  typename TransformType::Pointer trans = TransformType::New();
  ParaType                        para0(trans->GetNumberOfParameters()), para1(trans->GetNumberOfParameters());

  double       rval;
  SEARCH_POINT spt;

  typedef itk::CenteredEuler3DTransform<double> TransformType_Euler3D;
  ParaType para_cxy(TransformType_Euler3D::New()->GetNumberOfParameters()); // translated from pre registration

  // find initial center
  bool                  is_ok = false;
  itk::Point<double, 3> center;

  if (transform_initial.IsNull())
  {
    is_ok = register_image_cxyz(I_fixed, I_moving, para_cxy, rval);
    if (!is_ok)
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl << "initial affine registeration falied" << std::endl;
      std::exception();
    }
    center[0] = para_cxy[3];
    center[1] = para_cxy[4];
    center[2] = para_cxy[5];
  }
  else
  {
    center[0] = transform_initial->GetCenter()[0];
    center[1] = transform_initial->GetCenter()[1];
    center[2] = transform_initial->GetCenter()[2];
  }

  std::cout << std::endl;
  for (int n = 0; (number_of_seeds > 0) ? (n < number_of_seeds) : (n <= number_of_seeds); n++)
  {
    if (n == 0)
    {
      if (transform_initial.IsNull())
      {
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
        para0[10] = 0.0;
        para0[11] = 0.0;
        para0[12] = 0.0;
      }
      else // use input as intial
      {
        for (unsigned int i = 0; i < transform_initial->GetParameters().Size(); i++)
        {
          para0[i] = transform_initial->GetParameters()[i];
        }
      }
    }
    else
    {
      generate_search_seed_3d(search_list, para0);
    }

    // main registration using affine transform
    is_ok = register_image_affine3d_mres_mask(
      I_fixed, I_moving, mask_fixed_object, para0, center, number_of_iteration, MI_bins, MI_samples, para1, rval);
    if (!is_ok)
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                << "affine registration failed!" << std::endl
                << "use the initial parameters" << std::endl;
      // return -1;
    }

    spt.para0 = para0;
    spt.para1 = para1;
    spt.center = center;
    spt.rval = rval;
    spt.index = n;
    spt.number_of_iteration = number_of_iteration;

    std::cout << "para0: " << para0 << std::endl;
    std::cout << "para1: " << para1 << std::endl;
    std::cout << "center: " << center << std::endl;
    std::cout << "rval: " << rval << std::endl;
    std::cout << "add the search result to the list ... seed [" << n << "]" << std::endl << std::endl;

    add_search_point(search_list, spt);
  }
  get_best_search_point(search_list, spt);

  std::cout << std::endl << "History: " << std::endl;
  for (int ii = 0; ii < search_list.size(); ii++)
  {
    std::cout << "[" << ii << "]: " << search_list[ii].rval << std::endl;
    std::cout << "center: " << search_list[ii].center << std::endl;
    std::cout << "para0: " << search_list[ii].para0 << std::endl;
    std::cout << "para1: " << search_list[ii].para1 << std::endl;
  }

  typename TransformType::Pointer transform_final = TransformType::New();
  transform_final->SetParameters(spt.para1);
  transform_final->SetCenter(center);

  std::cout << "final transform  parameters = " << transform_final->GetParameters() << std::endl;

  transform = transform_final;
}

double
myrand(double l, double u)
{
  //  std::cout<<"in myrand" << std::endl;
  double r = l + (u - l) * rand() / (RAND_MAX + 1.0);

  return r;
  // std::cout<<"out myrand" << std::endl;
}

template <typename ParaType>
double
dist2_search_point(ParaType para0, ParaType para1)
{
  // use theta / scale / skew to compare two nodes
  double d2 = 0;

  //  double scale[4] = {1.0, 1.0, 1.0, 1.0};
  for (int ii = 0; ii < 12; ii++)
  {
    double a = para0[ii] - para1[ii];
    d2 += a * a; //  * scale[ii];
  }
  return d2;
}

template <typename SEARCH_LIST, typename ParaType>
void
generate_search_seed_3d(SEARCH_LIST & search_list, ParaType & para)
{
  bool   b_found = 0;
  double s1, s2, s3, k1, k2, k3 = 0;
  double a = 0;               // rotation alpha
  double u = 0, v = 0, w = 0; // rotation axis

  // ParaType para(13);

  constexpr double scale_upper = 1.5;
  const double     scale_lower = 1 / 1.5;
  const double     skew_lower = -1.0;
  constexpr double skew_upper = 1.0;
  constexpr double dist2_thres = 0.3;
  constexpr double mypi = 3.1415926536;
  const double     rot_angle_upper = mypi * (0.5);
  const double     rot_angle_lower = mypi * (-0.5);

  unsigned int iteration = 0;
  unsigned int maxiteration = 50;

  while (b_found == 0 && iteration < maxiteration)
  {
    iteration++;
    s1 = myrand(scale_lower, scale_upper);
    s2 = myrand(scale_lower, scale_upper);
    s3 = myrand(scale_lower, scale_upper);
    k1 = myrand(skew_lower, skew_upper);
    k2 = myrand(skew_lower, skew_upper);
    k3 = myrand(skew_lower, skew_upper);
    a = myrand(rot_angle_lower, rot_angle_upper);
    u = myrand(-1, 1);
    v = myrand(-1, 1);
    w = myrand(-1, 1);

    double n1 = 1.0 / sqrt(u * u + v * v + w * w);
    double sin_half_a = sin(a * 0.5);

    para.Fill(0.0);

    para[0] = sin_half_a * u * n1;
    para[1] = sin_half_a * v * n1;
    para[2] = sin_half_a * w * n1;
    para[3] = cos(a * 0.5);
    para[4] = s1;
    para[5] = s2;
    para[6] = s3;
    para[7] = k1;
    para[8] = k2;
    para[9] = k3;

    //
    std::cout << "test rand: " << para << " iteration " << iteration << std::endl;

    // search nearby search points
    bool bfar = 1;
    for (int i = 0; (bfar) & (i < search_list.size()); i++)
    {
      ParaType para0 = search_list[i].para0;
      ParaType para1 = search_list[i].para1;
      double   d0 = dist2_search_point(para0, para);
      double   d1 = dist2_search_point(para1, para);

      // std::cout << "compare with para0: " << d0 << para0 << std::endl;
      // std::cout << "compare with para1: " << d1 << para1 << std::endl;

      bfar = bfar & (d0 > dist2_thres) & (d1 > dist2_thres);
    }

    b_found = bfar;
    // std::cout << "b_found = " << b_found << " bfar = " << bfar << std::endl;
  }

  // ret = para;
}

template <typename SEARCH_LIST, typename ParaType>
void
generate_search_seed_2d(SEARCH_LIST & search_list, ParaType & para)
{
  bool   b_found = 0;
  double r1, s1, s2, k = 0;

  constexpr double pi = 3.1415927;
  const double     theta_upper = pi / 4;
  const double     theta_lower = -pi / 4;
  constexpr double scale_upper = 1.5;
  const double     scale_lower = 1 / 1.5;
  const double     skew_lower = -1.0;
  constexpr double skew_upper = 1.0;
  constexpr double dist2_thres = 0.1;

  unsigned int iteration = 0;
  unsigned int maxiteration = 50;

  while (b_found == 0 && iteration < maxiteration)
  {
    //  for(;~b_found; ){
    // std::cout << "b_found = " << b_found << std::endl;
    iteration++;
    r1 = myrand(theta_lower, theta_upper);
    s1 = myrand(scale_lower, scale_upper);
    s2 = myrand(scale_lower, scale_upper);
    k = myrand(skew_lower, skew_upper);
    para.Fill(0.0);
    para[0] = r1;
    para[1] = s1;
    para[2] = s2;
    para[3] = k;

    //
    // std::cout << "test rand: " << para << " iteration " << iteration <<  std::endl;

    // search nearby search points
    bool bfar = 1;
    for (int i = 0; (bfar) & (i < search_list.size()); i++)
    {
      ParaType para0 = search_list[i].para0;
      ParaType para1 = search_list[i].para1;
      double   d0 = dist2_search_point(para0, para);
      double   d1 = dist2_search_point(para1, para);

      // std::cout << "compare with para0: " << d0 << para0 << std::endl;
      // std::cout << "compare with para1: " << d1 << para1 << std::endl;

      bfar = bfar & (d0 > dist2_thres) & (d1 > dist2_thres);
    }

    b_found = bfar;
    // std::cout << "b_found = " << b_found << " bfar = " << bfar << std::endl;
  }
}

template <typename SEARCH_LIST, typename SEARCH_POINT>
void
get_best_search_point(SEARCH_LIST & search_list, SEARCH_POINT & spt)
{
  std::sort(search_list.begin(), search_list.end(), compare_search_point<SEARCH_POINT>);

  spt = search_list[0];
}

template <typename SEARCH_LIST, typename SEARCH_POINT>
void
add_search_point(SEARCH_LIST & search_list, SEARCH_POINT spt)
{
  search_list.push_back(spt);
}

template <typename ImagePointerType, typename ParaType, typename PointType, typename TransformTypePointer>
double
get_cost_value_mmi(ImagePointerType     fixedImage,
                   ImagePointerType     movingImage,
                   ParaType             para,
                   PointType            center,
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
  typename InterpolatorType::Pointer                                      interpolator = InterpolatorType::New();

  interpolator->SetInputImage(movingImage);

  mattesMutualInfo->SetFixedImage(fixedImage);
  mattesMutualInfo->SetMovingImage(movingImage);
  mattesMutualInfo->SetFixedImageRegion(fixedImage->GetBufferedRegion());
  //  mattesMutualInfo->SetMovingImage(outputImage);
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
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
              << "Exception caught in computing mattesMutualInfo after registration" << std::endl
              << "Maybe: Too many samples map outside moving image buffer" << std::endl
              << "Set the cost value = 0 (max for MutualInfo) " << std::endl;
    rval = 0;
  }

  // test the cost before registration
  transform->SetIdentity();
  transform->SetCenter(center);
  ParaType para0 = transform->GetParameters();
  mattesMutualInfo->SetTransformParameters(para0);
  mattesMutualInfo->Initialize();
  double rval0 = 0;
  try
  {
    rval0 = mattesMutualInfo->GetValue(para0);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
              << "Exception caught in computing mattesMutualInfo before registration" << std::endl
              << "Maybe: Too many samples map outside moving image buffer" << std::endl
              << "Set the cost value = 0 (max for MutualInfo) " << std::endl;
    rval0 = 0;
  }

  std::cout << "in cost: before register: cost = " << rval0 << std::endl;
  std::cout << "in cost: after register: cost = " << rval << std::endl;

  return rval;
}

template <typename ImagePointerType, typename ParaType>
bool
register_image_cxy(ImagePointerType fixed_image, ImagePointerType moving_image, ParaType & para1, double & rval)
{
  // for 3d image use CenteredRigid2DTransform
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::CenteredRigid2DTransform<double>    TransformType_Rigid2D;
  typedef typename ImagePointerType::ObjectType    ImageType;
  typename TransformType_Rigid2D::Pointer          transform = TransformType_Rigid2D::New();

  typedef itk::CenteredTransformInitializer<TransformType_Rigid2D, ImageType, ImageType> TransformInitializerType;

  typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  try
  {
    transform->SetIdentity();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixed_image);
    initializer->SetMovingImage(moving_image);
    initializer->MomentsOn();
    initializer->InitializeTransform();
    transform->SetAngle(0.0);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!1" << std::endl << "Exception in InitializeTransform" << std::endl;
    return false;
  }

  para1 = transform->GetParameters();

  std::cout << "finish initialize cx/cy/cz ..." << std::endl;
  std::cout << "cx/cy parameters (Euler3D): " << para1 << std::endl;

  return true;
}

template <typename ImagePointerType, typename ParaType>
bool
register_image_cxyz(ImagePointerType fixed_image, ImagePointerType moving_image, ParaType & para1, double & rval)
{
  // for 3d image use CenteredRigid2DTransform
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::CenteredEuler3DTransform<double>    TransformType_Euler3D;
  typedef typename ImagePointerType::ObjectType    ImageType;

  TransformType_Euler3D::Pointer transform = TransformType_Euler3D::New();

  typedef itk::CenteredTransformInitializer<TransformType_Euler3D, ImageType, ImageType> TransformInitializerType;

  typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  try
  {
    transform->SetIdentity();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixed_image);
    initializer->SetMovingImage(moving_image);
    initializer->MomentsOn();
    initializer->InitializeTransform();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!1" << std::endl << "Exception in InitializeTransform" << std::endl;
    return false;
  }

  para1 = transform->GetParameters();

  std::cout << "finish initialize cx/cy/cz ..." << std::endl;
  std::cout << "cx/cy parameters (Euler3D): " << para1 << std::endl;

  return true;
}

template <typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
bool
register_image_affine3d_mres_mask(ImagePointerType                  fixed_image,
                                  ImagePointerType                  moving_image,
                                  ImageMaskSpatialObjectPointerType mask_fixed_object,
                                  ParaType                          para0,
                                  itk::Point<double, 3>             center,
                                  int                               number_of_iteration,
                                  int                               MI_bins,
                                  int                               MI_samples,
                                  ParaType &                        para1,
                                  double &                          rval)
{
  typedef typename ImagePointerType::ObjectType                                ImageType;
  typedef itk::RegularStepGradientDescentOptimizer                             OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> MetricType;
  typedef typename ImagePointerType::ObjectType                                ImageType;

  typename MetricType::Pointer metric = MetricType::New();

  // for mattesMutualInfo
  metric->SetNumberOfHistogramBins(MI_bins);     /** 32 BA */
  metric->SetNumberOfSpatialSamples(MI_samples); /** 6000 BA */
  if (mask_fixed_object)
  {
    metric->SetFixedImageMask(mask_fixed_object);
    std::cout << mask_fixed_object << std::endl;
    std::cout << mask_fixed_object->GetImage() << std::endl;
  }

  // typedef TransformType_Rigid2D TransformType_Pre;

  typedef itk::LinearInterpolateImageFunction<ImageType, double>            InterpolatorType;
  typedef itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType> RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> fPyramidType;
  typename fPyramidType::Pointer                                       fixedImagePyramid = fPyramidType::New();
  typename fPyramidType::Pointer                                       movingImagePyramid = fPyramidType::New();

  typename ImageType::PointType origin;
  origin.Fill(0);
  moving_image->SetOrigin(origin);
  fixed_image->SetOrigin(origin);

  typename OptimizerType::Pointer    optimizer = OptimizerType::New();
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename RegistrationType::Pointer registration = RegistrationType::New();

  typedef OptimizerType::ScalesType OptimizerScalesType;

  /*******************************************/
  /* translate to para0 here */

  std::cout << "pre registration para :" << para0 << std::endl;

  typedef itk::ANTSAffine3DTransform<double>   TransformType_ANTSAffine3D;
  typename TransformType_ANTSAffine3D::Pointer transform_a = TransformType_ANTSAffine3D::New();
  transform_a->SetCenter(center);
  std::cout << "initial center: " << transform_a->GetCenter() << std::endl;

  //  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales_a(transform_a->GetNumberOfParameters());
  const double        translationScale = 1.0 / 1.e4; /** BA */ // /1.e6; // 1.e6; //1.e6; //1.0 / 1000.0;
  optimizerScales_a[0] = 1.0;                                  // quaternion
  optimizerScales_a[1] = 1.0;                                  // quaternion
  optimizerScales_a[2] = 1.0;                                  // quaternion
  optimizerScales_a[3] = 1.0;                                  // quaternion
  optimizerScales_a[4] = 1.0;                                  // s1
  optimizerScales_a[5] = 1.0;                                  // s2
  optimizerScales_a[6] = 1.0;                                  // s3
  optimizerScales_a[7] = 1.0;                                  // k1
  optimizerScales_a[8] = 1.0;                                  // k2
  optimizerScales_a[9] = 1.0;                                  // k3
  optimizerScales_a[10] = translationScale;
  optimizerScales_a[11] = translationScale;
  optimizerScales_a[12] = translationScale;
  optimizer->SetScales(optimizerScales_a);

  optimizer->SetMaximumStepLength(0.1);
  optimizer->SetMinimumStepLength(1.e-5);
  //  optimizer->SetNumberOfIterations( 1000 );
  //  optimizer->SetNumberOfIterations( 500 );
  optimizer->SetNumberOfIterations(number_of_iteration); // (500); /** BA */
  optimizer->MinimizeOn();

  // Create the Command observer and register it with the optimizer.
  //
  // CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  // optimizer->AddObserver( itk::IterationEvent(), observer );

  registration = RegistrationType::New();

  registration->SetNumberOfLevels(3);
  registration->SetFixedImagePyramid(fixedImagePyramid);
  registration->SetMovingImagePyramid(movingImagePyramid);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  registration->SetFixedImage(fixed_image);
  registration->SetMovingImage(moving_image);
  registration->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion());
  registration->SetTransform(transform_a); // reset the parameter
  //  registration->SetInitialTransformParameters(para_pre );
  registration->SetInitialTransformParameters(para0);

  std::cout << "reset initial transform parameters" << std::endl;
  std::cout << "para_pre: " << para0 << std::endl;
  // std::cout << "transform_a: " << transform_a << std::endl;

  rval = get_cost_value_mmi(fixed_image, moving_image, para0, center, transform_a);
  std::cout << "init measure value: rval = " << rval << std::endl;

  // rval = optimizer->GetValue(para0);
  // rval = metric->GetValue(para0);

  bool bsuc = 1;
  try
  {
    registration->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    bsuc = 0;
    // std::exception();
  }

  if (bsuc)
  {
    // OptimizerType::ParametersType
    // finalParameters = registration->GetLastTransformParameters();
    para1 = registration->GetLastTransformParameters();
    //  para1 = finalParameters;
    //  rval = optimizer->GetValue(para1);
    // rval = metric->GetValue(para1);
    // rval = registration->GetMetric()->GetValue(para1);

    rval = get_cost_value_mmi(fixed_image, moving_image, para1, center, transform_a);
    // double rval2 = optimizer->GetValue();
    // std::cout << "measure value: rval2 = " << rval2 << std::endl;
  }
  else
  {
    // register failed
    para1 = para0; // registration->GetLastTransformParameters();
    rval = get_cost_value_mmi(fixed_image, moving_image, para1, center, transform_a);
  }

  std::cout << "final affine3d registration para :" << para1 << std::endl;
  std::cout << "use iteration: " << optimizer->GetNumberOfIterations() << std::endl;
  std::cout << "measure value: rval = " << rval << std::endl;
  std::cout << "finish register..." << std::endl;

  return bsuc;
}

template <typename ImagePointerType, typename ImageMaskSpatialObjectPointerType, typename ParaType>
bool
register_image_affine2d_mres_mask(ImagePointerType                  fixed_image,
                                  ImagePointerType                  moving_image,
                                  ImageMaskSpatialObjectPointerType mask_fixed_object,
                                  ParaType                          para0,
                                  itk::Point<double, 2>             center,
                                  int                               number_of_iteration,
                                  int                               MI_bins,
                                  int                               MI_samples,
                                  ParaType &                        para1,
                                  double &                          rval)
{
  typedef typename ImagePointerType::ObjectType                                ImageType;
  typedef itk::RegularStepGradientDescentOptimizer                             OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> MetricType;
  typedef typename ImagePointerType::ObjectType                                ImageType;

  typename MetricType::Pointer metric = MetricType::New();

  // for mattesMutualInfo
  metric->SetNumberOfHistogramBins(MI_bins);     /** 32 BA */
  metric->SetNumberOfSpatialSamples(MI_samples); /** 6000 BA */
  if (mask_fixed_object)
  {
    metric->SetFixedImageMask(mask_fixed_object);
    std::cout << mask_fixed_object << std::endl;
    std::cout << mask_fixed_object->GetImage() << std::endl;
  }

  // typedef TransformType_Rigid2D TransformType_Pre;

  typedef itk::LinearInterpolateImageFunction<ImageType, double>            InterpolatorType;
  typedef itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType> RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> fPyramidType;
  typename fPyramidType::Pointer                                       fixedImagePyramid = fPyramidType::New();
  typename fPyramidType::Pointer                                       movingImagePyramid = fPyramidType::New();

  typename ImageType::PointType origin;
  origin.Fill(0);
  moving_image->SetOrigin(origin);
  fixed_image->SetOrigin(origin);

  typename OptimizerType::Pointer    optimizer = OptimizerType::New();
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename RegistrationType::Pointer registration = RegistrationType::New();

  typedef OptimizerType::ScalesType OptimizerScalesType;

  /*******************************************/
  /* translate to para0 here */

  std::cout << "pre registration para :" << para0 << std::endl;

  // typedef itk::CenteredAffine2DTransform<double> TransformType;
  //    typedef itk::CenteredAffine2DTransform<double> TransformType_GSAffine2D;
  //    typename TransformType_GSAffine2D::Pointer  transform_a = TransformType_GSAffine2D::New();
  typedef itk::ANTSCenteredAffine2DTransform<double> TransformType_ANTSAffine2D;
  typename TransformType_ANTSAffine2D::Pointer       transform_a = TransformType_ANTSAffine2D::New();
  // transform_a->SetCenter(center);
  // std::cout<<"initial center: " << transform_a->GetCenter() << std::endl;

  OptimizerScalesType optimizerScales_a(transform_a->GetNumberOfParameters());
  const double        translationScale = 1.0 / 1000.0;
  optimizerScales_a[0] = 1.0;
  optimizerScales_a[1] = 1.0;
  optimizerScales_a[2] = 1.0;
  optimizerScales_a[3] = 1.0;
  optimizerScales_a[4] = translationScale;
  optimizerScales_a[5] = translationScale;
  optimizerScales_a[6] = translationScale;
  optimizerScales_a[7] = translationScale;
  optimizer->SetScales(optimizerScales_a);

  optimizer->SetMaximumStepLength(0.1);
  optimizer->SetMinimumStepLength(0.01);
  optimizer->SetNumberOfIterations(number_of_iteration);
  optimizer->MinimizeOn();

  // Create the Command observer and register it with the optimizer.
  //
  // CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  // optimizer->AddObserver( itk::IterationEvent(), observer );

  registration = RegistrationType::New();

  registration->SetNumberOfLevels(3);
  registration->SetFixedImagePyramid(fixedImagePyramid);
  registration->SetMovingImagePyramid(movingImagePyramid);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  registration->SetFixedImage(fixed_image);
  registration->SetMovingImage(moving_image);
  registration->SetFixedImageRegion(fixed_image->GetLargestPossibleRegion());
  registration->SetTransform(transform_a); // reset the parameter
  //  registration->SetInitialTransformParameters(para_pre );
  registration->SetInitialTransformParameters(para0);

  std::cout << "reset initial transform parameters" << std::endl;
  std::cout << "para_pre: " << para0 << std::endl;
  // std::cout << "transform_a: " << transform_a << std::endl;

  rval = get_cost_value_mmi(fixed_image, moving_image, para0, center, transform_a);
  std::cout << "init measure value: rval = " << rval << std::endl;

  // rval = optimizer->GetValue(para0);
  // rval = metric->GetValue(para0);

  bool bsuc = 1;
  try
  {
    registration->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    bsuc = 0;
    // std::exception();
  }

  if (bsuc)
  {
    // OptimizerType::ParametersType
    // finalParameters = registration->GetLastTransformParameters();
    para1 = registration->GetLastTransformParameters();
    //  para1 = finalParameters;
    //  rval = optimizer->GetValue(para1);
    // rval = metric->GetValue(para1);
    // rval = registration->GetMetric()->GetValue(para1);

    rval = get_cost_value_mmi(fixed_image, moving_image, para1, center, transform_a);
    // double rval2 = optimizer->GetValue();
    // std::cout << "measure value: rval2 = " << rval2 << std::endl;
  }
  else
  {
    // register failed
    para1 = para0; // registration->GetLastTransformParameters();
    rval = get_cost_value_mmi(fixed_image, moving_image, para1, center, transform_a);
  }

  std::cout << "final affine2d registration para :" << para1 << std::endl;
  std::cout << "use iteration: " << optimizer->GetNumberOfIterations() << std::endl;
  std::cout << "measure value: rval = " << rval << std::endl;
  std::cout << "finish register..." << std::endl;

  return bsuc;
}

// template<typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
// void compute_single_affine_transform_2d(ImagePointerType I_fixed, ImagePointerType I_moving, MaskImagePointerType
// mask_fixed, TransformPointerType &transform){

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform_2d(ImagePointerType       I_fixed,
                                   ImagePointerType       I_moving,
                                   MaskImagePointerType   mask_fixed,
                                   TransformPointerType & transform,
                                   TransformPointerType & transform_initial)
{
  typedef typename ImagePointerType::ObjectType     ImageType;
  const int                                         ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::IOPixelType           PixelType;
  typedef typename MaskImagePointerType::ObjectType MaskImageType;
  typedef typename MaskImageType::IOPixelType       MaskPixelType;

  typedef typename TransformPointerType::ObjectType TransformType;
  typedef typename TransformType::ParametersType    ParaType;

  // option parameters
  int          number_of_seeds = 0;
  int          number_of_iteration = 500;
  int          MI_bins = 32;
  int          MI_samples = 6000;
  unsigned int time_seed = (unsigned int)time(nullptr);
  srand(time_seed);
  // TODO: need to fix here
  bool b_use_mask = 0; // (mask_fixed == nullptr);

  std::cout << "number_of_seeds: " << number_of_seeds << std::endl;
  std::cout << "rand_time_seed: " << time_seed << std::endl;
  std::cout << "number_of_iteration: " << number_of_iteration << std::endl;
  std::cout << "MI_bins: " << MI_bins << std::endl;
  std::cout << "MI_samples: " << MI_samples << std::endl;
  std::cout << "use mask: " << b_use_mask << std::endl;

  // memory of searched results
  typedef SEARCH_POINT_TYPE<ParaType, ImageDimension> SEARCH_POINT;
  typedef std::vector<SEARCH_POINT>                   SEARCH_LIST;

  SEARCH_LIST search_list;

  // typedef itk::ImageMaskSpatialObject<ImageDimension, MaskPixelType> ImageMaskSpatialObject;
  typedef typename itk::ImageMaskSpatialObject<ImageDimension> ImageMaskSpatialObject;
  typename ImageMaskSpatialObject::Pointer                     mask_fixed_object = 0;
  if (b_use_mask)
  {
    typedef typename itk::Image<unsigned char, ImageDimension>              CharMaskImageType;
    typedef typename itk::CastImageFilter<MaskImageType, CharMaskImageType> CastFilterType;
    typename CastFilterType::Pointer                                        cast_filter = CastFilterType::New();
    cast_filter->SetInput(mask_fixed);
    cast_filter->Update();
    typename CharMaskImageType::Pointer mask_fixed_char = cast_filter->GetOutput();

    mask_fixed_object = ImageMaskSpatialObject::New();
    mask_fixed_object->SetImage(mask_fixed_char);
  }

  typename ImageType::PointType origin;
  origin.Fill(0);
  I_moving->SetOrigin(origin);
  I_fixed->SetOrigin(origin);

  typename TransformType::Pointer trans = TransformType::New();
  ParaType                        para0(trans->GetNumberOfParameters()), para1(trans->GetNumberOfParameters());

  double       rval;
  SEARCH_POINT spt;

  typedef typename itk::CenteredRigid2DTransform<double> TransformType_Rigid2D;
  ParaType para_cxy(TransformType_Rigid2D::New()->GetNumberOfParameters()); // translated from pre registration

  // find initial center
  bool                  is_ok = false;
  itk::Point<double, 2> center;

  if (transform_initial.IsNull())
  {
    //    if (1) {
    is_ok = register_image_cxy(I_fixed, I_moving, para_cxy, rval);
    if (!is_ok)
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl << "initial affine registeration falied" << std::endl;
      std::exception();
    }

    center[0] = para_cxy[1];
    center[1] = para_cxy[2];

    std::cout << "is_ok=" << is_ok << "para_cxy:" << para_cxy << std::endl;
  }
  else
  {
    center[0] = transform_initial->GetCenter()[0];
    center[1] = transform_initial->GetCenter()[1];
    std::cout << "input transform: " << transform_initial << std::endl;
  }

  std::cout << "initial center: (" << center[0] << "," << center[1] << ")" << std::endl;
  for (int n = 0; (number_of_seeds > 0) ? (n < number_of_seeds) : (n <= number_of_seeds); n++)
  {
    if (n == 0)
    {
      // if (1) {
      if (transform_initial.IsNull())
      {
        para0[0] = 0;           // para1[0]; // theta
        para0[1] = 1.0;         // s1
        para0[2] = 1.0;         // s2
        para0[3] = 0.0;         // k
        para0[4] = center[0];   // para1[1]; //c1
        para0[5] = center[1];   // para1[2]; //c2
        para0[6] = para_cxy[3]; // 0;//para1[3]; //t1
        para0[7] = para_cxy[4]; // 0; //para1[4]; //t2

        std::cout << "ABC: " << std::endl;
      }
      else
      {
        for (unsigned int i = 0; i < transform_initial->GetParameters().Size(); i++)
        {
          para0[i] = transform_initial->GetParameters()[i];
        }
        std::cout << "DEF: " << std::endl;
      }
    }
    else
    {
      generate_search_seed_2d(search_list, para0);
    }

    // main registration using affine transform
    is_ok = register_image_affine2d_mres_mask(
      I_fixed, I_moving, mask_fixed_object, para0, center, number_of_iteration, MI_bins, MI_samples, para1, rval);
    if (!is_ok)
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                << "affine registration failed!" << std::endl
                << "use the initial parameters" << std::endl;
      // return -1;
    }

    spt.para0 = para0;
    spt.para1 = para1;
    spt.center = center;
    spt.rval = rval;
    spt.index = n;
    spt.number_of_iteration = number_of_iteration;

    std::cout << "para0: " << para0 << std::endl;
    std::cout << "para1: " << para1 << std::endl;
    std::cout << "center: " << center << std::endl;
    std::cout << "rval: " << rval << std::endl;
    std::cout << "add the search result to the list ... seed [" << n << "]" << std::endl << std::endl;

    add_search_point(search_list, spt);
  }
  get_best_search_point(search_list, spt);

  std::cout << std::endl << "History: " << std::endl;
  for (int ii = 0; ii < search_list.size(); ii++)
  {
    std::cout << "[" << ii << "]: " << search_list[ii].rval << std::endl;
    std::cout << "center: " << search_list[ii].center << std::endl;
    std::cout << "para0: " << search_list[ii].para0 << std::endl;
    std::cout << "para1: " << search_list[ii].para1 << std::endl;
  }

  typename TransformType::Pointer transform_final = TransformType::New();
  transform_final->SetParameters(spt.para1);
  transform_final->SetCenter(center);

  std::cout << "final transform  parameters = " << transform_final->GetParameters() << std::endl;

  transform = transform_final;
}

template <typename ImagePointerType, typename MaskImagePointerType, typename TransformPointerType>
void
compute_single_affine_transform(ImagePointerType       fixedImage,
                                ImagePointerType       movingImage,
                                MaskImagePointerType   maskImage,
                                TransformPointerType & transform,
                                TransformPointerType & transform_initial)
{
  typedef typename ImagePointerType::ObjectType     ImageType;
  const int                                         ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::IOPixelType           PixelType;
  typedef typename MaskImagePointerType::ObjectType MaskImageType;
  typedef typename MaskImageType::IOPixelType       MaskPixelType;
  typedef typename TransformPointerType::ObjectType TransformType;

  std::cout << "transform_initial: IsNotNull():" << transform_initial.IsNotNull() << std::endl;

  if (ImageDimension == 2)
  {
    typedef itk::ANTSCenteredAffine2DTransform<double> RunningAffineTransformType;
    RunningAffineTransformType::Pointer                transform_running = RunningAffineTransformType::New();
    RunningAffineTransformType::Pointer transform_running_initial; // = RunningAffineTransformType::New();

    std::cout << "1: transform_running_initial: IsNotNull():" << transform_running_initial.IsNotNull() << std::endl;

    if (transform_initial.IsNotNull())
    {
      transform_running_initial->SetCenter(*(reinterpret_cast<typename RunningAffineTransformType::InputPointType *>(
        const_cast<typename TransformType::InputPointType *>(&(transform_initial->GetCenter())))));
      transform_running_initial->SetMatrix(*(reinterpret_cast<typename RunningAffineTransformType::MatrixType *>(
        const_cast<typename TransformType::MatrixType *>(&(transform_initial->GetMatrix())))));
      transform_running_initial->SetTranslation(
        *(reinterpret_cast<typename RunningAffineTransformType::OutputVectorType *>(
          const_cast<typename TransformType::OutputVectorType *>(&(transform_initial->GetTranslation())))));
    }

    // Use type casting
    typedef typename itk::Image<PixelType, 2>::Pointer R_ImagePointerType;
    R_ImagePointerType & R_fixedImage = reinterpret_cast<R_ImagePointerType &>(fixedImage);
    R_ImagePointerType & R_movingImage = reinterpret_cast<R_ImagePointerType &>(movingImage);
    R_ImagePointerType & R_maskImage = reinterpret_cast<R_ImagePointerType &>(maskImage);

    std::cout << "2: transform_running_initial: IsNotNull():" << transform_running_initial.IsNotNull() << std::endl;

    compute_single_affine_transform_2d(
      R_fixedImage, R_movingImage, R_maskImage, transform_running, transform_running_initial);

    // TODO:
    transform->SetCenter(*(reinterpret_cast<typename TransformType::InputPointType *>(
      const_cast<typename RunningAffineTransformType::InputPointType *>(&(transform_running->GetCenter())))));
    transform->SetTranslation(*(reinterpret_cast<typename TransformType::OutputVectorType *>(
      const_cast<typename RunningAffineTransformType::OutputVectorType *>(&(transform_running->GetTranslation())))));
    transform->SetMatrix(*(reinterpret_cast<typename TransformType::MatrixType *>(
      const_cast<typename RunningAffineTransformType::MatrixType *>(&(transform_running->GetMatrix())))));

    // transform->SetFixedParameters(transform_running->GetFixedParameters());
    // transform->SetParameters(transform_running->GetParameters());
  }
  else if (ImageDimension == 3)
  {
    typedef itk::ANTSAffine3DTransform<double> RunningAffineTransformType;
    RunningAffineTransformType::Pointer        transform_running = RunningAffineTransformType::New();
    RunningAffineTransformType::Pointer        transform_running_initial = RunningAffineTransformType::New();
    ;
    // compute_single_affine_transform_3d(fixedImage, movingImage, maskImage, transform_running);

    if (transform_initial.IsNotNull())
    {
      transform_running_initial->SetCenter(*(reinterpret_cast<typename RunningAffineTransformType::InputPointType *>(
        const_cast<typename TransformType::InputPointType *>(&(transform_initial->GetCenter())))));
      transform_running_initial->SetMatrix(*(reinterpret_cast<typename RunningAffineTransformType::MatrixType *>(
        const_cast<typename TransformType::MatrixType *>(&(transform_initial->GetMatrix())))));
      transform_running_initial->SetTranslation(
        *(reinterpret_cast<typename RunningAffineTransformType::OutputVectorType *>(
          const_cast<typename TransformType::OutputVectorType *>(&(transform_initial->GetTranslation())))));
    }

    // Use type casting
    typedef typename itk::Image<PixelType, 3>::Pointer R_ImagePointerType;
    R_ImagePointerType & R_fixedImage = reinterpret_cast<R_ImagePointerType &>(fixedImage);
    R_ImagePointerType & R_movingImage = reinterpret_cast<R_ImagePointerType &>(movingImage);
    R_ImagePointerType & R_maskImage = reinterpret_cast<R_ImagePointerType &>(maskImage);

    compute_single_affine_transform_3d(
      R_fixedImage, R_movingImage, R_maskImage, transform_running, transform_running_initial);

    // TODO:
    transform->SetCenter(*(reinterpret_cast<typename TransformType::InputPointType *>(
      const_cast<typename RunningAffineTransformType::InputPointType *>(&(transform_running->GetCenter())))));
    transform->SetTranslation(*(reinterpret_cast<typename TransformType::OutputVectorType *>(
      const_cast<typename RunningAffineTransformType::OutputVectorType *>(&(transform_running->GetTranslation())))));
    transform->SetMatrix(*(reinterpret_cast<typename TransformType::MatrixType *>(
      const_cast<typename RunningAffineTransformType::MatrixType *>(&(transform_running->GetMatrix())))));

    //        transform->SetFixedParameters(transform_running->GetFixedParameters());
    //        transform->SetParameters(transform_running->GetParameters());
  }
  else
  {
    std::cout << "Unsupported, not 2D/ 3D" << std::endl;
    return;
  }
}

template <typename DisplacementFieldPointerType>
void
create_deformation_field_byref(const DisplacementFieldPointerType & ref, DisplacementFieldPointerType & field)
{
  typedef typename DisplacementFieldPointerType::ObjectType DisplacementFieldType;
  // field = DisplacementFieldType::New();

  typename DisplacementFieldType::RegionType region;
  region.SetSize(ref->GetLargestPossibleRegion().GetSize());
  region.SetIndex(ref->GetLargestPossibleRegion().GetIndex());
  field->SetRegions(region);
  field->SetSpacing(ref->GetSpacing());
  field->SetOrigin(ref->GetOrigin());
  field->AllocateInitialized();
}

// compose affine transform (in a matrix format A: (Ax+b)) with a deformation field F:
// the new field is: F_new (x)  = F ( A (x) )
// output should be allocated outside
template <typename TransformPointerType, typename DisplacementFieldPointerType>
void
compose_affine_with_field(const TransformPointerType &         aff,
                          const DisplacementFieldPointerType & field,
                          DisplacementFieldPointerType &       field_output)
{
  typedef typename DisplacementFieldPointerType::ObjectType        DisplacementFieldType;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef typename DisplacementFieldType::IndexType                IndexType;
  typedef typename DisplacementFieldType::PointType                PointType;
  typedef typename DisplacementFieldType::PixelType                VectorType;

  const unsigned int ImageDimension = DisplacementFieldType::ImageDimension;

  //    PointType zeroorigin;
  //    zeroorigin.Fill(0);
  //    field->SetOrigin(zeroorigin);
  //    field_output->SetOrigin(zeroorigin);

  PointType pointIn1;
  PointType pointIn2;
  PointType pointIn3;

  // iterate through field_output finding the points that it maps to via field.
  // then take the difference from the original point and put it in the output field.
  // std::cout << " begin iteration " << std::endl;
  FieldIterator iter_field(field, field->GetLargestPossibleRegion());

  // std::cout << field_output->GetLargestPossibleRegion() << std::endl;

  int cnt = 0;
  for (iter_field.GoToBegin(); !iter_field.IsAtEnd(); ++iter_field)
  {
    IndexType index = iter_field.GetIndex();
    field_output->TransformIndexToPhysicalPoint(index, pointIn1);
    VectorType disp = iter_field.Get();
    for (int jj = 0; jj < ImageDimension; jj++)
    {
      pointIn2[jj] = disp[jj] + pointIn1[jj];
    }

    // use affine transform
    pointIn3 = aff->TransformPoint(pointIn2);

    VectorType out;
    for (int jj = 0; jj < ImageDimension; jj++)
    {
      out[jj] = pointIn3[jj] - pointIn1[jj];
    }

    field_output->SetPixel(iter_field.GetIndex(), out);
  }

  // std::cout << " end iteration " << std::endl;
}

// this is obsolet, use itkWarpImageWAffineFilter
template <typename ImagePointerType, typename DisplacementFieldPointerType>
void
warp_image_field(const ImagePointerType &             img_input,
                 const DisplacementFieldPointerType & field,
                 ImagePointerType &                   img_output)
{
  typedef typename ImagePointerType::ObjectType             ImageType;
  typedef typename DisplacementFieldPointerType::ObjectType DisplacementFieldType;

  typedef typename itk::WarpImageFilter<ImageType, ImageType, DisplacementFieldType> WarperType;
  typename WarperType::Pointer                                                       warper = WarperType::New();

  warper->SetInput(img_input);
  warper->SetDisplacementField(field);
  warper->SetEdgePaddingValue(0);
  warper->SetOutputSpacing(field->GetSpacing());
  warper->SetOutputOrigin(field->GetOrigin());
  warper->Update();

  img_output = warper->GetOutput();
}

template <typename ImageTypePointer, typename RefImageTypePointer, typename TransformPointerType>
void
affine_image(const ImageTypePointer &     input_image,
             const TransformPointerType & transform,
             const RefImageTypePointer &  ref_image,
             ImageTypePointer &           img_aff)
{
  typedef typename ImageTypePointer::ObjectType     ImageType;
  typedef typename TransformPointerType::ObjectType TransformType;

  // apply the transform
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer                   resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer                             interpolator = InterpolatorType::New();
  resampler->SetInterpolator(interpolator);
  resampler->SetInput(input_image);
  resampler->SetSize(ref_image->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(ref_image->GetOrigin());
  resampler->SetOutputSpacing(ref_image->GetSpacing());
  resampler->SetDefaultPixelValue(0);
  resampler->SetTransform(transform);

  resampler->Update();
  img_aff = resampler->GetOutput();
}

template <typename ImagePointerType, typename TransformPointerType, typename DisplacementFieldPointerType>
void
warp_image_field_waffine(const ImagePointerType &             img_input,
                         const TransformPointerType &         aff,
                         const DisplacementFieldPointerType & field,
                         ImagePointerType &                   img_output)
{
  // TODO: add a new WarpImageFilter to support affine as an input
  // temporary solution:
  typedef typename DisplacementFieldPointerType::ObjectType DisplacementFieldType;
  typename DisplacementFieldType::Pointer                   field_comp = DisplacementFieldType::New();

  //    create_deformation_field_byref(field, field_comp);
  //    compose_affine_with_field(aff, field, field_comp);
  //    warp_image_field(img_input, field_comp, img_output);

  typedef typename TransformPointerType::ObjectType                                               TransformType;
  typedef typename ImagePointerType::ObjectType                                                   ImageType;
  typedef itk::WarpImageWAffineFilter<ImageType, ImageType, DisplacementFieldType, TransformType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  warper->SetInput(img_input);
  warper->SetDisplacementField(field);
  warper->SetAffineTransform(aff);
  warper->SetEdgePaddingValue(0);
  warper->SetOutputSpacing(field->GetSpacing());
  warper->SetOutputOrigin(field->GetOrigin());
  warper->SetOutputSize(field->GetLargestPossibleRegion().GetSize());
  warper->Update();

  img_output = warper->GetOutput();
  return;
}

// TODO: use my own code to implement all the optimization codes

#endif /*ANTS_AFFINE_REGISTRATION_H_*/
