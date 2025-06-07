/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkComposeDiffeomorphismsImageFilter_h
#define __itkComposeDiffeomorphismsImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorInterpolateImageFunction.h"

namespace itk
{
/**
 * \class ComposeDiffeomorphismsImageFilter
 *
 * \brief
 *
 * \par
 *
 * \author Nicholas J. Tustison
 */

template <typename TInputImage, typename TOutputImage = TInputImage>
class ComposeDiffeomorphismsImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef ComposeDiffeomorphismsImageFilter             Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  typedef TInputImage  InputFieldType;
  typedef TOutputImage OutputFieldType;

  /** Image typedef support. */
  typedef typename OutputFieldType::PixelType  PixelType;
  typedef typename OutputFieldType::PixelType  VectorType;
  typedef typename OutputFieldType::RegionType RegionType;
  typedef typename OutputFieldType::IndexType  IndexType;

  typedef typename OutputFieldType::PointType     PointType;
  typedef typename OutputFieldType::SpacingType   SpacingType;
  typedef typename OutputFieldType::PointType     OriginType;
  typedef typename OutputFieldType::SizeType      SizeType;
  typedef typename OutputFieldType::DirectionType DirectionType;

  /** Other typedef */
  typedef typename VectorType::ComponentType                       RealType;
  typedef VectorInterpolateImageFunction<InputFieldType, RealType> InterpolatorType;

  /** Get the interpolator. */
  itkGetModifiableObjectMacro(Interpolator, InterpolatorType);

  /** Set the deformation field */
  void
  SetDisplacementField(const InputFieldType * field)
  {
    itkDebugMacro("setting deformation field to " << field);
    if (field != this->GetInput(0))
    {
      this->SetInput(0, field);
      this->Modified();
      if (!this->m_Interpolator.IsNull())
      {
        this->m_Interpolator->SetInputImage(field);
      }
    }
  }

  /**
   * Get the deformation field.
   */
  const InputFieldType *
  GetDisplacementField() const
  {
    return this->GetInput(0);
  }

  /** Set the warping field */
  void
  SetWarpingField(const InputFieldType * field)
  {
    itkDebugMacro("setting warping field to " << field);
    if (field != this->GetInput(1))
    {
      this->SetInput(1, field);
    }
  }

  /**
   * Get the warping field.
   */
  const InputFieldType *
  GetWarpingField() const
  {
    return this->GetInput(1);
  }

  /* Set the interpolator. */
  virtual void
  SetInterpolator(InterpolatorType * interpolator);

protected:
  /** Constructor */
  ComposeDiffeomorphismsImageFilter();

  /** Deconstructor */
  virtual ~ComposeDiffeomorphismsImageFilter();

  /** Standard print self function **/
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** preprocessing function */
  void
  BeforeThreadedGenerateData();

  /** Multithreaded function which generates the output field. */
  void
  ThreadedGenerateData(const RegionType &, ThreadIdType);

private:
  ComposeDiffeomorphismsImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  /** The interpolator. */
  typename InterpolatorType::Pointer m_Interpolator;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkComposeDiffeomorphismsImageFilter.hxx"
#endif

#endif
