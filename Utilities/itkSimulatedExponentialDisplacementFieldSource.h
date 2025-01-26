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
#ifndef itkSimulatedExponentialDisplacementFieldSource_h
#define itkSimulatedExponentialDisplacementFieldSource_h

#include "itkSimulatedDisplacementFieldSource.h"

namespace itk
{
/** \class SimulatedExponentialDisplacementFieldSource
 * \brief Computes a randomly SimulatedExponential displacement field.
 *
 * SimulatedExponentialDisplacementFieldSource produces a random exponential
 * displacement field.
 *
 * This source object expects the image to be of pixel type Vector.
 *
 * \ingroup ImageSource
 * \ingroup ITKDisplacementField
 */
template <typename TOutputImage>
class SimulatedExponentialDisplacementFieldSource final : public SimulatedDisplacementFieldSource<TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SimulatedExponentialDisplacementFieldSource);

  /** Standard class type aliases. */
  using Self = SimulatedExponentialDisplacementFieldSource;
  using Superclass = SimulatedDisplacementFieldSource<TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelComponentType = typename Superclass::OutputPixelComponentType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SimulatedExponentialDisplacementFieldSource);

  /** Number of dimensions */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  using OutputPixelType = typename Superclass::OutputPixelType;
  using VectorType = typename Superclass::VectorType;
  using RealType = typename Superclass::RealType;
  using RealImageType = typename Superclass::RealImageType;

  using StandardDeviationArrayType = FixedArray<RealType, ImageDimension>;

  /** Image spacing type alias */
  using SpacingType = typename Superclass::SpacingType;
  using OriginType = typename Superclass::OriginType;
  using DirectionType = typename Superclass::DirectionType;
  using RegionType = typename Superclass::RegionType;
  using SizeType = typename Superclass::SizeType;

  /**
   * Set/Get the number of integration steps.  If 0, this number is computed
   * automatically.  Default = 0.
   */
  itkSetMacro(NumberOfIntegrationSteps, unsigned int);
  itkGetConstMacro(NumberOfIntegrationSteps, unsigned int);

  /**
   * Set/Get the dilation radius.  Default = 5.
   */
  itkSetMacro(DilationRadius, unsigned int);
  itkGetConstMacro(DilationRadius, unsigned int);

  /**
   * Set/Get the standard deviation of the noise for the displacements at the
   * random points.  Default = c( 1, 1, ... ) mm.
   */
  itkSetMacro(DisplacementNoiseStandardDeviation, StandardDeviationArrayType);
  itkGetConstMacro(DisplacementNoiseStandardDeviation, StandardDeviationArrayType);

  /**
   * Set the standard deviation of the displacement noise.  Default = 1 mm.
   */
  void
  SetDisplacementNoiseStandardDeviation(RealType sd)
  {
    StandardDeviationArrayType standardDeviation;

    standardDeviation.Fill(sd);
    this->SetDisplacementNoiseStandardDeviation(standardDeviation);
  }

  /**
   * Set/Get the standard deviation of the noise for the smoothing.  Default = c( 1, 1, ... )
   */
  itkSetMacro(SmoothingStandardDeviation, RealType);
  itkGetConstMacro(SmoothingStandardDeviation, RealType);

protected:
  SimulatedExponentialDisplacementFieldSource();
  ~SimulatedExponentialDisplacementFieldSource() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  unsigned int               m_NumberOfIntegrationSteps;
  RealType                   m_GradientStep;
  unsigned int               m_DilationRadius;
  StandardDeviationArrayType m_DisplacementNoiseStandardDeviation;
  RealType                   m_SmoothingStandardDeviation;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSimulatedExponentialDisplacementFieldSource.hxx"
#endif

#endif
