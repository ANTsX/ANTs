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
#ifndef itkSimulatedDisplacementFieldSource_h
#define itkSimulatedDisplacementFieldSource_h

#include "itkImageSource.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

namespace itk
{
/** \class SimulatedDisplacementFieldSource
 * \brief Computes a randomly simulated displacement field.
 *
 * SimulatedDisplacementFieldSource is an abstract class for producing a random
 * displacement field.  It is used to define the domain of the displacement
 * field.
 *
 * This source object expects the image to be of pixel type Vector.
 *
 * \ingroup ImageSource
 * \ingroup ITKDisplacementField
 */
template <typename TOutputImage>
class SimulatedDisplacementFieldSource : public ImageSource<TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SimulatedDisplacementFieldSource);

  /** Standard class type aliases. */
  using Self = SimulatedDisplacementFieldSource;
  using Superclass = ImageSource<TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SimulatedDisplacementFieldSource);

  /** Number of dimensions. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Image size type alias. */
  using OutputSizeType = typename OutputImageType::SizeType;

  /** Image index type alias. */
  using OutputIndexType = typename OutputImageType::IndexType;

  /** Image pixel value type alias. */
  using OutputPixelType = typename TOutputImage::PixelType;
  using OutputPixelComponentType = typename OutputPixelType::ValueType;
  using VectorType = typename OutputImageType::PixelType;
  using RealType = typename VectorType::RealValueType;
  using RealImageType = Image<RealType, ImageDimension>;

  /** Image spacing type alias */
  using SpacingType = typename TOutputImage::SpacingType;
  using OriginType = typename TOutputImage::PointType;
  using DirectionType = typename TOutputImage::DirectionType;
  using RegionType = typename TOutputImage::RegionType;
  using SizeType = typename RegionType::SizeType;

  /** Randomizer typedefs */
  using RandomizerType = Statistics::MersenneTwisterRandomVariateGenerator;
  using RandomizerPointer = typename RandomizerType::Pointer;
  using RandomizerSeedType = typename RandomizerType::IntegerType;

  /** Define the output displacement field domain from an image */
  void
  SetDisplacementFieldDomainFromImage(RealImageType *);

  /** Define the output displacement field domain from a const image */
  void
  SetDisplacementFieldDomainFromImage(const RealImageType * image)
  {
    this->SetDisplacementFieldDomainFromImage(const_cast<RealImageType *>(image));
  }

  /** Define the output displacement field domain from a displacement field */
  void
  SetDisplacementFieldDomainFromField(OutputImageType *);

  /** Define the output displacement field domain from a const displacement field */
  void
  SetDisplacementFieldDomainFromField(const OutputImageType * field)
  {
    this->SetDisplacementFieldDomainFromField(const_cast<OutputImageType *>(field));
  }

  /** Define the displacement field domain explicitly. */
  void SetDisplacementFieldDomain(OriginType, SpacingType, SizeType, DirectionType);

  /** Set/Get the size of the output image. */
  itkSetMacro(OutputSize, SizeType);
  itkGetConstReferenceMacro(OutputSize, SizeType);

  /** Set/Get the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set/Get the output image origin. */
  itkSetMacro(OutputOrigin, OriginType);
  itkGetConstReferenceMacro(OutputOrigin, OriginType);

  /** Set the output direction cosine matrix. */
  itkSetMacro(OutputDirection, DirectionType);
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Set the number of random points.  Default = 100. */
  itkSetMacro(NumberOfRandomPoints, SizeValueType);
  itkGetConstMacro(NumberOfRandomPoints, SizeValueType);

  /**
   * Enforce stationary boundary conditions.  Default = true.
   */
  itkBooleanMacro(EnforceStationaryBoundary);
  itkSetMacro(EnforceStationaryBoundary, bool);
  itkGetConstMacro(EnforceStationaryBoundary, bool);

  /**
   * Set/Get initialization random see for random number
   * generator.  Default is to initialize randomly using the system
   * clock.
   */
  void
  SetRandomizerInitializationSeed(const RandomizerSeedType);
  itkGetConstMacro(RandomizerInitializationSeed, RandomizerSeedType);

  /**
   * Get randomizer.
   */
  RandomizerType *
  GetRandomizer() const
  {
    return static_cast<RandomizerType *>(this->m_Randomizer);
  }

  /**
   * SimulatedDisplacementFieldSource produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  void
  GenerateOutputInformation() override;

protected:
  SimulatedDisplacementFieldSource();
  ~SimulatedDisplacementFieldSource() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /**
   * GenerateData() computes the simulated displacement field.  Is overwritten
   * in derived classes.
   */
  void
  GenerateData() override{};

private:
  SizeType      m_OutputSize;
  SpacingType   m_OutputSpacing;
  OriginType    m_OutputOrigin;
  DirectionType m_OutputDirection;

  RandomizerPointer  m_Randomizer;
  RandomizerSeedType m_RandomizerInitializationSeed;

  bool m_EnforceStationaryBoundary;

  SizeValueType m_NumberOfRandomPoints;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSimulatedDisplacementFieldSource.hxx"
#endif

#endif
