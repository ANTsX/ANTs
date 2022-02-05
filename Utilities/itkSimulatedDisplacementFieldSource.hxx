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
#ifndef itkSimulatedDisplacementFieldSource_hxx
#define itkSimulatedDisplacementFieldSource_hxx


namespace itk
{
template <typename TOutputImage>
SimulatedDisplacementFieldSource<TOutputImage>::SimulatedDisplacementFieldSource()
  : m_EnforceStationaryBoundary(true)
  , m_NumberOfRandomPoints(100)
{
  this->m_OutputSpacing.Fill(1.0);
  this->m_OutputOrigin.Fill(0.0);
  this->m_OutputSize.Fill(0);
  this->m_OutputDirection.SetIdentity();

  this->m_RandomizerInitializationSeed = std::numeric_limits<RandomizerSeedType>::quiet_NaN();
  this->m_Randomizer = RandomizerType::New();
  this->m_Randomizer->Initialize();
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::SetDisplacementFieldDomainFromImage(RealImageType * image)
{
  this->SetDisplacementFieldDomain(
    image->GetOrigin(), image->GetSpacing(), image->GetRequestedRegion().GetSize(), image->GetDirection());
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::SetDisplacementFieldDomainFromField(OutputImageType * field)
{
  this->SetDisplacementFieldDomain(
    field->GetOrigin(), field->GetSpacing(), field->GetRequestedRegion().GetSize(), field->GetDirection());
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::SetDisplacementFieldDomain(OriginType    origin,
                                                                           SpacingType   spacing,
                                                                           SizeType      size,
                                                                           DirectionType direction)
{
  if (this->m_OutputOrigin != origin || this->m_OutputSpacing != spacing || this->m_OutputSize != size ||
      this->m_OutputDirection != direction)
  {
    this->m_OutputOrigin = origin;
    this->m_OutputSpacing = spacing;
    this->m_OutputSize = size;
    this->m_OutputDirection = direction;

    this->Modified();
  }
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::SetRandomizerInitializationSeed(const RandomizerSeedType seed)
{
  if (seed != this->m_RandomizerInitializationSeed)
  {
    this->m_RandomizerInitializationSeed = seed;
    this->m_Randomizer->Initialize(this->m_RandomizerInitializationSeed);
    this->Modified();
  }
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  OutputImagePointer output = this->GetOutput();
  if (!output)
  {
    return;
  }

  // Define displacement field domain
  output->SetRegions(this->m_OutputSize);
  output->SetSpacing(this->m_OutputSpacing);
  output->SetOrigin(this->m_OutputOrigin);
  output->SetDirection(this->m_OutputDirection);
}

template <typename TOutputImage>
void
SimulatedDisplacementFieldSource<TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Enforce stationary boundary: ";
  if (this->m_EnforceStationaryBoundary)
  {
    os << "true" << std::endl;
  }
  else
  {
    os << "false" << std::endl;
  }

  os << indent << "Displacement field domain:" << std::endl;
  os << indent << "  Origin: " << this->m_OutputOrigin << std::endl;
  os << indent << "  Spacing: " << this->m_OutputSpacing << std::endl;
  os << indent << "  Size: " << this->m_OutputSize << std::endl;
  os << indent << "  Direction: " << this->m_OutputDirection << std::endl;
}


} // end namespace itk

#endif
