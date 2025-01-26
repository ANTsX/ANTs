/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWarpImageWAffineFilter_h
#define __itkWarpImageWAffineFilter_h

#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPoint.h"
#include "itkFixedArray.h"
#include "itkRecursiveGaussianImageFilter.h"

namespace itk
{
/** \class WarpImageWAffineFilter
 * \brief Warps an image using an input deformation field.
 *
 * WarpImageWAffineFilter warps an existing image with respect to
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
class WarpImageWAffineFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** transform order type **/
  typedef enum _TransformOrderType
  {
    AffineFirst = 0,
    AffineLast
  } TransformOrderType;

  /** Standard class typedefs. */
  typedef WarpImageWAffineFilter                        Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(WarpImageWAffineFilter);

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

  /** Determine the image dimension. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int DisplacementFieldDimension = TDisplacementField::ImageDimension;

  /** Deformation field typedef support. */
  typedef TDisplacementField                        DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename DisplacementFieldType::PixelType DisplacementType;

  /** songgang: Affine transform typedef support. */
  typedef TTransform                      TransformType;
  typedef typename TransformType::Pointer TransformTypePointer;

  /** Interpolator typedef support. */
  typedef double                                                       CoordRepType;
  typedef InterpolateImageFunction<InputImageType, CoordRepType>       InterpolatorType;
  typedef typename InterpolatorType::Pointer                           InterpolatorPointer;
  typedef LinearInterpolateImageFunction<InputImageType, CoordRepType> DefaultInterpolatorType;

  /** Point type */
  typedef Point<CoordRepType, Self::ImageDimension> PointType;

  /** Set the deformation field. */
  void
  SetDisplacementField(const DisplacementFieldType * field);

  /** Get a pointer the deformation field. */
  DisplacementFieldType *
  GetDisplacementField();

  /** songgang: Set / Get  the affine transform. */
  // void SetAffineTransform( const TransformType * aff );
  // TransformType * GetAffineTransform();

  itkSetObjectMacro(AffineTransform, TransformType);

  itkGetModifiableObjectMacro(AffineTransform, TransformType);

  itkSetEnumMacro(TransformOrder, TransformOrderType);
  itkGetEnumMacro(TransformOrder, TransformOrderType);

  /** Set the interpolator function. */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetModifiableObjectMacro(Interpolator, InterpolatorType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void
  SetOutputSpacing(const double * values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, PointType);
  virtual void
  SetOutputOrigin(const double * values);

  /** Set the output image size. */
  itkSetMacro(OutputSize, SizeType);
  // virtual void SetOutputSize( const double *values);
  itkGetConstReferenceMacro(OutputSize, SizeType);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, PointType);

  /** Set the edge padding value */
  itkSetMacro(EdgePaddingValue, PixelType);

  /** Get the edge padding value */
  itkGetMacro(EdgePaddingValue, PixelType);

  /** WarpImageWAffineFilter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  virtual void
  GenerateOutputInformation();

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void
  GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before
   * multi-threading. */
  virtual void
  BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after
   * multi-threading. */
  virtual void
  AfterThreadedGenerateData();

  /** precompute the smoothed image if necessary **/
  void
  SetSmoothScale(double scale);

  double
  GetSmoothScale()
  {
    return m_SmoothScale;
  };

  void
  UpdateSizeByScale();

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1, (Concept::SameDimension<ImageDimension, InputImageDimension>));
  itkConceptMacro(SameDimensionCheck2, (Concept::SameDimension<ImageDimension, DisplacementFieldDimension>));
  itkConceptMacro(InputHasNumericTraitsCheck, (Concept::HasNumericTraits<typename TInputImage::PixelType>));
  itkConceptMacro(DisplacementFieldHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<typename TDisplacementField::PixelType::ValueType>));
  /** End concept checking */
#endif
protected:
  WarpImageWAffineFilter();
  ~WarpImageWAffineFilter() = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** WarpImageWAffineFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

private:
  WarpImageWAffineFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  PixelType          m_EdgePaddingValue;
  SpacingType        m_OutputSpacing;
  PointType          m_OutputOrigin;
  SizeType           m_OutputSize;
  TransformOrderType m_TransformOrder;

  // songgang: affine transform
  TransformTypePointer m_AffineTransform;

  InterpolatorPointer m_Interpolator;

  /** cached smoothed image (for downsampling) **/
  InputImagePointer m_CachedSmoothImage;
  double            m_SmoothScale;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkWarpImageWAffineFilter.hxx"
#endif

#endif
