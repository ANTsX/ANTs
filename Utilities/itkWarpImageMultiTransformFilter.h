/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef ITKWARPIMAGEMULTITRANSFORMFILTER_H_
#define ITKWARPIMAGEMULTITRANSFORMFILTER_H_

#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorGaussianInterpolateImageFunction.h"
#include "itkPoint.h"
#include "itkFixedArray.h"
#include "itkRecursiveGaussianImageFilter.h"
#include <list>

namespace itk
{
/** \class WarpImageMultiTransformFilter
 * \brief Warps an image using an input deformation field.
 *
 * WarpImageMultiTransformFilter warps an existing image with respect to
 * a given deformation field.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the input image. The vector type must support element access via operator
 * [].
 *
 * The output image is produced by inverse mapping: the output pixels
 * are mapped back onto the input image. This scheme avoids the creation of
 * any holes and overlaps in the output image.
 *
 * Each vector in the deformation field represent the distance between
 * a geometric point in the input space and a point in the output space such
 * that:
 *
 * \f[ p_{in} = p_{out} + d \f]
 *
 * Typically the mapped position does not correspond to an integer pixel
 * position in the input image. Interpolation via an image function
 * is used to compute values at non-integer positions. The default
 * interpolation typed used is the LinearInterpolateImageFunction.
 * The user can specify a particular interpolation function via
 * SetInterpolator(). Note that the input interpolator must derive
 * from base class InterpolateImageFunction.
 *
 * Position mapped to outside of the input image buffer are assigned
 * a edge padding value.
 *
 * The LargetPossibleRegion for the output is inherited from the
 * input deformation field. The output image spacing and origin may be set
 * via SetOutputSpacing, SetOutputOrigin. The default are respectively a
 * vector of 1's and a vector of 0's.
 *
 * This class is templated over the type of the input image, the
 * type of the output image and the type of the deformation field.
 *
 * The input image is set via SetInput. The input deformation field
 * is set via SetDisplacementField.
 *
 * This filter is implemented as a multithreaded filter.
 *
 * \warning This filter assumes that the input type, output type
 * and deformation field type all have the same number of dimensions.
 *
 * \ingroup GeometricTransforms MultiThreaded
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
class WarpImageMultiTransformFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** transform order type **/
  // typedef enum _TransformOrderType {AffineFirst=0, AffineLast} TransformOrderType;
  /** transform type **/
  typedef enum _SingleTransformType
  {
    EnumAffineType = 0,
    EnumDisplacementFieldType
  } SingleTransformType;

  /** Standard class typedefs. */
  typedef WarpImageMultiTransformFilter                 Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(WarpImageMultiTransformFilter);

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::PixelType         PixelType;
  typedef typename OutputImageType::SpacingType       SpacingType;
  typedef typename OutputImageType::DirectionType     DirectionType;

  /** Determine the image dimension. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int DisplacementFieldDimension = TDisplacementField::ImageDimension;

  /** base type for images of the current ImageDimension */
  typedef ImageBase<Self::ImageDimension> ImageBaseType;

  /** Deformation field typedef support. */
  typedef TDisplacementField                        DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename DisplacementFieldType::PixelType DisplacementType;
  typedef typename DisplacementType::ValueType      DisplacementScalarValueType;

  /** songgang: Affine transform typedef support. */
  typedef TTransform                      TransformType;
  typedef typename TransformType::Pointer TransformTypePointer;

  /** Interpolator typedef support. */
  typedef double                                                                      CoordRepType;
  typedef InterpolateImageFunction<InputImageType, CoordRepType>                      InterpolatorType;
  typedef typename InterpolatorType::Pointer                                          InterpolatorPointer;
  typedef LinearInterpolateImageFunction<InputImageType, CoordRepType>                DefaultInterpolatorType;
  typedef VectorLinearInterpolateImageFunction<DisplacementFieldType, CoordRepType>   DefaultVectorInterpolatorType;
  typedef VectorGaussianInterpolateImageFunction<DisplacementFieldType, CoordRepType> DefaultVectorInterpolatorType2;
  typedef typename DefaultVectorInterpolatorType::Pointer                             VectorInterpolatorPointer;

  /** Point type */
  typedef Point<CoordRepType, Self::ImageDimension> PointType;

  typedef struct _DeformationTypeEx
  {
    DisplacementFieldPointer  field;
    VectorInterpolatorPointer vinterp;
  } DeformationTypeEx;

  typedef struct _AffineTypeEx
  {
    TransformTypePointer aff;
  } AffineTypeEx;

  typedef struct _VarTransformType
  {
    AffineTypeEx      aex;
    DeformationTypeEx dex;
  } VarTransformType;

  typedef std::pair<SingleTransformType, VarTransformType> SingleTransformItemType;
  typedef std::list<SingleTransformItemType>               TransformListType;

  /** Set the interpolator function. */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetModifiableObjectMacro(Interpolator, InterpolatorType);

  /** Convenience funciton to set all the necessary items */
  virtual void
  SetOutputParametersFromImage(const ImageBaseType * image);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image spacing. */
  itkSetMacro(OutputDirection, DirectionType);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, PointType);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, PointType);

  /** Set the output image size. */
  itkSetMacro(OutputSize, SizeType);

  // virtual void SetOutputSize( const double *values);
  itkGetConstReferenceMacro(OutputSize, SizeType);

  /** Set the start index of the output largest possible region.
   * The default is an index of all zeros. */
  itkSetMacro(OutputStartIndex, IndexType);

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro(OutputStartIndex, IndexType);

  /** Set the edge padding value */
  itkSetMacro(EdgePaddingValue, PixelType);

  /** Get the edge padding value */
  itkGetMacro(EdgePaddingValue, PixelType);

  TransformListType &
  GetTransformList()
  {
    return m_TransformList;
  }

  void
  PrintTransformList();

  /** WarpImageMultiTransformFilter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  void
  GenerateOutputInformation() override;

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  void
  GenerateInputRequestedRegion() override;

  /** This method is used to set the state of the filter before
   * multi-threading. */
  void
  BeforeThreadedGenerateData() override;

  /** This method is used to set the state of the filter after
   * multi-threading. */
  void
  AfterThreadedGenerateData() override;

  /** precompute the smoothed image if necessary **/
  void
  SetSmoothScale(double scale);

  double
  GetSmoothScale()
  {
    return m_SmoothScale;
  };

  // void UpdateSizeByScale();

  void
  PushBackAffineTransform(const TransformType * t);

  void
  PushBackDisplacementFieldTransform(const DisplacementFieldType * t);

  void
  ComposeAffineOnlySequence(const PointType & center_output, TransformTypePointer & affine_output);

  bool
  MultiInverseAffineOnlySinglePoint(const PointType & point1, PointType & point2);

  bool
  MultiTransformSinglePoint(const PointType & point1, PointType & point2);

  bool
  MultiTransformPoint(const PointType & point1, PointType & point2, bool bFisrtDeformNoInterp, const IndexType & index);

  void
  DetermineFirstDeformNoInterp();

  inline bool
  IsOutOfNumericBoundary(const PointType & p);

  // set interpolator from outside
  // virtual void SetInterpolator(InterpolatorPointer interp);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1, (Concept::SameDimension<ImageDimension, InputImageDimension>));
  itkConceptMacro(SameDimensionCheck2, (Concept::SameDimension<ImageDimension, DisplacementFieldDimension>));

  // removed to be compatible with vector form input image
  //    itkConceptMacro(InputHasNumericTraitsCheck,
  //            (Concept::HasNumericTraits<typename TInputImage::PixelType>));
  itkConceptMacro(DisplacementFieldHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<typename TDisplacementField::PixelType::ValueType>));
  /** End concept checking */
#endif

  bool m_bFirstDeformNoInterp;

protected:
  WarpImageMultiTransformFilter();
  ~WarpImageMultiTransformFilter() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** WarpImageMultiTransformFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  PointType     m_OutputOrigin;
  SpacingType   m_OutputSpacing;
  DirectionType m_OutputDirection;
  SizeType      m_OutputSize;
  IndexType     m_OutputStartIndex; // output image start index

  InterpolatorPointer m_Interpolator;

  PixelType         m_EdgePaddingValue;
  TransformListType m_TransformList;

  double m_SmoothScale;

  InputImagePointer m_CachedSmoothImage;

private:
  WarpImageMultiTransformFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkWarpImageMultiTransformFilter.hxx"
#endif

#endif /*ITKWARPIMAGEMULTITRANSFORMFILTER_H_*/
