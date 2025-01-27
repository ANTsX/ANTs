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
#ifndef itkSimulatedBSplineDisplacementFieldSource_h
#define itkSimulatedBSplineDisplacementFieldSource_h

#include "itkSimulatedDisplacementFieldSource.h"
#include "itkDisplacementFieldToBSplineImageFilter.h"
#include "itkFixedArray.h"

namespace itk
{
/** \class SimulatedBSplineDisplacementFieldSource
 * \brief Computes a randomly SimulatedBSpline displacement field.
 *
 * SimulatedBSplineDisplacementFieldSource produces a random B-spline
 * displacement field.
 *
 * This source object expects the image to be of pixel type Vector.
 *
 * \ingroup ImageSource
 * \ingroup ITKDisplacementField
 */
template <typename TOutputImage>
class SimulatedBSplineDisplacementFieldSource final : public SimulatedDisplacementFieldSource<TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SimulatedBSplineDisplacementFieldSource);

  /** Standard class type aliases. */
  using Self = SimulatedBSplineDisplacementFieldSource;
  using Superclass = SimulatedDisplacementFieldSource<TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SimulatedBSplineDisplacementFieldSource);

  /** Number of dimensions */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  using OutputPixelType = typename Superclass::OutputPixelType;
  using VectorType = typename Superclass::VectorType;
  using RealType = typename Superclass::RealType;
  using RealImageType = typename Superclass::RealImageType;

  /** B-sline filter type alias */
  using BSplineFilterType = DisplacementFieldToBSplineImageFilter<OutputImageType>;
  using ArrayType = typename BSplineFilterType::ArrayType;
  using StandardDeviationArrayType = FixedArray<RealType, ImageDimension>;

  /** Point set type alias support */
  using PointType = typename BSplineFilterType::PointType;
  using PointDataType = typename BSplineFilterType::PixelType;
  using PointSetType = typename BSplineFilterType::InputPointSetType;

  /** Image spacing type alias */
  using SpacingType = typename Superclass::SpacingType;
  using OriginType = typename Superclass::OriginType;
  using DirectionType = typename Superclass::DirectionType;
  using RegionType = typename Superclass::RegionType;
  using SizeType = typename Superclass::SizeType;

  /**
   * Set/Get the spline order defining the bias field estimate.  Default = 3.
   */
  itkSetMacro(SplineOrder, unsigned int);
  itkGetConstMacro(SplineOrder, unsigned int);

  /**
   * Set/Get the control point grid size definining the B-spline estimate of the
   * scalar bias field.  In each dimension, the B-spline mesh size is equal
   * to the number of control points in that dimension minus the spline order.
   * Default = 4 control points in each dimension for a mesh size of 1 in each
   * dimension.
   */
  itkSetMacro(NumberOfControlPoints, ArrayType);
  itkGetConstMacro(NumberOfControlPoints, ArrayType);

  /**
   * Set/Get the number of fitting levels.  Default = 1 level.
   */
  itkSetMacro(NumberOfFittingLevels, ArrayType);
  itkGetConstMacro(NumberOfFittingLevels, ArrayType);

  /**
   * Set the number of fitting levels.  Default = 1 level.
   */
  void
  SetNumberOfFittingLevels(unsigned int n)
  {
    ArrayType nlevels;

    nlevels.Fill(n);
    this->SetNumberOfFittingLevels(nlevels);
  }

  /**
   * Set/Get the standard deviation of the noise for the displacements at the
   * random points.  Default = c( 1, 1, ... ) mm.
   */
  itkSetMacro(DisplacementNoiseStandardDeviation, StandardDeviationArrayType);
  itkGetConstMacro(DisplacementNoiseStandardDeviation, StandardDeviationArrayType);

  /**
   * Set the standard deviation of the displacement noise.  Default = 1 (mm).
   */
  void
  SetDisplacementNoiseStandardDeviation(RealType sd)
  {
    StandardDeviationArrayType standardDeviation;

    standardDeviation.Fill(sd);
    this->SetDisplacementNoiseStandardDeviation(standardDeviation);
  }

protected:
  SimulatedBSplineDisplacementFieldSource();
  ~SimulatedBSplineDisplacementFieldSource() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  unsigned int               m_SplineOrder;
  ArrayType                  m_NumberOfControlPoints;
  ArrayType                  m_NumberOfFittingLevels;
  StandardDeviationArrayType m_DisplacementNoiseStandardDeviation;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSimulatedBSplineDisplacementFieldSource.hxx"
#endif

#endif
