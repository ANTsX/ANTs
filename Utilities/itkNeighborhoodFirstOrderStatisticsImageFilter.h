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
#ifndef __itkNeighborhoodFirstOrderStatisticsImageFilter_h
#define __itkNeighborhoodFirstOrderStatisticsImageFilter_h

#include "itkMovingHistogramImageFilter.h"
#include "itkTextureHistogram.h"

namespace itk
{
/**
 * \class NeighborhoodFirstOrderStatisticsImageFilter
 * \brief Compute first order statistics in a neighborhood at each pixel
 *
 * \ingroup ITK-TextureAnalysis
 */

template <typename TInputImage, typename TOutputImage, typename TKernel>
class ITK_EXPORT NeighborhoodFirstOrderStatisticsImageFilter
  : public MovingHistogramImageFilter<
      TInputImage,
      TOutputImage,
      TKernel,
      typename Function::TextureHistogram<typename TInputImage::PixelType, typename TOutputImage::PixelType>>
{
public:
  /** Standard class typedefs. */
  typedef NeighborhoodFirstOrderStatisticsImageFilter Self;
  typedef MovingHistogramImageFilter<
    TInputImage,
    TOutputImage,
    TKernel,
    typename Function::TextureHistogram<typename TInputImage::PixelType, typename TOutputImage::PixelType>>
                                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(NeighborhoodFirstOrderStatisticsImageFilter);

  /** Image related typedefs. */
  typedef TInputImage                                InputImageType;
  typedef TOutputImage                               OutputImageType;
  typedef typename TInputImage::RegionType           RegionType;
  typedef typename TInputImage::SizeType             SizeType;
  typedef typename TInputImage::IndexType            IndexType;
  typedef typename TInputImage::PixelType            PixelType;
  typedef typename TInputImage::OffsetType           OffsetType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename TOutputImage::PixelType           OutputPixelType;

  /** Image related typedefs. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

protected:
  unsigned int
  GetNumberOfOutputComponents()
  {
    return 8;
  }

  NeighborhoodFirstOrderStatisticsImageFilter()
  {
    // this->m_Boundary = NumericTraits< PixelType >::max();
  }

  void
  GenerateOutputInformation() override
  {
    // this methods is overloaded so that if the output image is a
    // VectorImage then the correct number of components are set.

    Superclass::GenerateOutputInformation();
    OutputImageType * output = this->GetOutput();

    if (!output)
    {
      return;
    }
    if (output->GetNumberOfComponentsPerPixel() != this->GetNumberOfOutputComponents())
    {
      output->SetNumberOfComponentsPerPixel(this->GetNumberOfOutputComponents());
    }
  }


  ~NeighborhoodFirstOrderStatisticsImageFilter() override = default;

private:
  NeighborhoodFirstOrderStatisticsImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;
}; // end of class
} // end namespace itk

#endif
