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
#ifndef __itkSliceTimingCorrectionImageFilter_hxx
#define __itkSliceTimingCorrectionImageFilter_hxx

#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::SliceTimingCorrectionImageFilter()
{
  m_Interpolator = dynamic_cast<InterpolatorType *>(DefaultInterpolatorType::New().GetPointer());
  m_IndexPadding = 1;
  m_TimeDimension = InputImageDimension - 1;
  m_SliceDimension = InputImageDimension - 2;
  m_SliceTiming = 0.0;
  m_ExtrapolateEdges = true;
}

template <typename TInputImage, typename TOutputImage>
void
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "TimeDimension: " << m_TimeDimension << std::endl;
}

template <typename TInputImage, typename TOutputImage>
void
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::VerifyInputInformation() const
{
  Superclass::VerifyInputInformation();

  itkDebugMacro(<< "Verify Input Information");

  // Verify that all input have the same number of components

  typename InputImageType::ConstPointer image = this->GetInput();

  if (image.IsNull())
  {
    itkExceptionMacro(<< "Input not set as expected!");
  }

  const unsigned int numComponents = image->GetNumberOfComponentsPerPixel();
  for (IndexValueType idx = 1; idx < this->GetNumberOfInputs(); ++idx)
  {
    image = this->GetInput(idx);

    // if the input was not set it could still be null
    if (image.IsNull())
    {
      // an invalid requested region exception will be generated later.
      continue;
    }

    if (numComponents != image->GetNumberOfComponentsPerPixel())
    {
      itkExceptionMacro(<< "Primary input has " << numComponents << " numberOfComponents "
                        << "but input " << idx << " has " << image->GetNumberOfComponentsPerPixel() << "!");
    }
  }
}

/**
 * \sa UnaryFunctorImageFilter::GenerateOutputInformation()
 */
template <typename TInputImage, typename TOutputImage>
void
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{

  // std::cout << "GenerateOutputInformation" << std::endl;
  itkDebugMacro(<< "GenerateOutputInformation");

  // do not call the superclass' implementation of this method since
  // this filter allows the input and output to be of different size

  // get pointers to the input and output
  typename Superclass::OutputImagePointer     outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer inputPtr = this->GetInput();

  if (!outputPtr || !inputPtr)
  {
    return;
  }

  // Check for valid time dimension
  if (!(this->m_TimeDimension < InputImageDimension))
  {
    itkExceptionMacro(<< "Time dimension must be within dimensions of the input image");
  }

  // Set the output image largest possible region.  Use a RegionCopier
  // so that the input and output images can be different dimensions.
  OutputImageRegionType outputLargestPossibleRegion;
  this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion, inputPtr->GetLargestPossibleRegion());

  TimingType timingCoverage = inputPtr->GetLargestPossibleRegion().GetSize()[m_SliceDimension] * m_SliceTiming;

  itkDebugMacro("Timing coverage = " << timingCoverage);
  double timingDiff = inputPtr->GetSpacing()[m_TimeDimension] - timingCoverage;


  if (timingDiff < 0)
  {
    std::cout << "Timing diff = " << timingDiff << std::endl;
    std::cout << "Timing coverage = " << timingCoverage << std::endl;
    std::cout << "Volume timing = " << inputPtr->GetSpacing()[m_TimeDimension] << std::endl;
    std::cout << "Slice dimension = " << m_SliceDimension << std::endl;
    std::cout << "Slice Dimension Size = " << inputPtr->GetLargestPossibleRegion().GetSize()[m_SliceDimension]
              << std::endl;
    itkExceptionMacro(<< "SliceTiming * Dim[SliceDimension] should not be greater than Spacing[TimeDimension]");
  }

  // for the time dimension
  if ((2 * m_IndexPadding) >= this->GetInput()->GetLargestPossibleRegion().GetSize()[m_TimeDimension])
  {
    itkExceptionMacro(<< "IndexPadding is too large for this input image");
  }

  // unsigned int nTimePointsOut = this->GetInput()->GetLargestPossibleRegion().GetSize()[m_TimeDimension]
  //  - 2 * m_IndexPadding;
  unsigned int nTimePointsOut = this->GetInput()->GetLargestPossibleRegion().GetSize()[m_TimeDimension];

  outputLargestPossibleRegion.SetSize(m_TimeDimension, nTimePointsOut);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);

  // Set the output spacing and origin
  const ImageBase<InputImageDimension> * phyData;

  phyData = dynamic_cast<const ImageBase<InputImageDimension> *>(this->GetInput());

  if (phyData)
  {
    // Copy what we can from the image from spacing and origin of the input
    // This logic needs to be augmented with logic that select which
    // dimensions to copy
    unsigned int                                 ii;
    const typename InputImageType::SpacingType & inputSpacing = inputPtr->GetSpacing();
    const typename InputImageType::PointType &   inputOrigin = inputPtr->GetOrigin();

    typename OutputImageType::SpacingType outputSpacing;
    typename OutputImageType::PointType   outputOrigin;
    // copy the input to the output and fill the rest of the
    // output with zeros.
    for (ii = 0; ii < OutputImageDimension; ++ii)
    {
      outputSpacing[ii] = inputSpacing[ii];
      outputOrigin[ii] = inputOrigin[ii];
    }
    // outputOrigin[ m_TimeDimension ] += outputSpacing[ m_TimeDimension ]*m_IndexPadding;

    // set the spacing and origin
    outputPtr->SetSpacing(outputSpacing);
    outputPtr->SetOrigin(outputOrigin);

    //
    // Copy the direction cosines from the input to the output.
    // On join, the output dim is always >= input dim
    typedef typename InputImageType::DirectionType  InputDirectionType;
    typedef typename OutputImageType::DirectionType OutputDirectionType;
    InputDirectionType                              inputDir = inputPtr->GetDirection();
    unsigned int                                    outputdim = OutputImageType::GetImageDimension();
    OutputDirectionType                             outputDir = outputPtr->GetDirection();
    for (unsigned int i = 0; i < outputdim; i++)
    {
      for (unsigned int j = 0; j < outputdim; j++)
      {
        outputDir[i][j] = inputDir[i][j];
      }
    }
    outputPtr->SetDirection(outputDir);
  }
  else
  {
    // pointer could not be cast back down
    itkExceptionMacro(<< "itk::SliceTimingCorrectionImageFilter::GenerateOutputInformation "
                      << "cannot cast input to " << typeid(ImageBase<InputImageDimension> *).name());
  }

  // Support VectorImages by setting number of components on output.
  const unsigned int numComponents = inputPtr->GetNumberOfComponentsPerPixel();
  if (numComponents != outputPtr->GetNumberOfComponentsPerPixel())
  {
    outputPtr->SetNumberOfComponentsPerPixel(numComponents);
  }
}


/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage>
void
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  itkDebugMacro(<< "Before Threaded Generate Data");
  if (!m_Interpolator)
  {
    itkExceptionMacro(<< "Interpolator not set");
  }

  // Connect input image to interpolator
  m_Interpolator->SetInputImage(this->GetInput());

  // Connect input image to extrapolator
  // if( !m_Extrapolator.IsNull() )
  //  {
  //  m_Extrapolator->SetInputImage( this->GetInput() );
  //  }
}

template <typename TInputImage, typename TOutputImage>
void
SliceTimingCorrectionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  itkDebugMacro(<< "Actually executing");
  // std::cout << "Actually executing" << std::endl;

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  OutputImageRegionType outputRegion = outputRegionForThread;

  InputImageRegionType inputRegion;
  this->CallCopyOutputRegionToInputRegion(inputRegion, outputRegionForThread);

  ImageRegionIteratorWithIndex<OutputImageType> outIt(this->GetOutput(), outputRegion);

  typename InputImageType::PointType pt;

  while (!outIt.IsAtEnd())
  {

    typename InputImageType::IndexType idx = outIt.GetIndex();
    unsigned int maxDim = this->GetOutput()->GetLargestPossibleRegion().GetSize()[m_TimeDimension] - 1;
    // std::cout << pt << std::endl;

    if (this->m_ExtrapolateEdges)
    {
      if (idx[m_TimeDimension] < m_IndexPadding)
      {
        idx[m_TimeDimension] = m_IndexPadding;
      }
      else if (idx[m_TimeDimension] > (maxDim - m_IndexPadding))
      {
        idx[m_TimeDimension] = maxDim - m_IndexPadding;
      }
    }

    this->GetOutput()->TransformIndexToPhysicalPoint(idx, pt);
    pt[m_TimeDimension] -= idx[m_SliceDimension] * m_SliceTiming;

    if (this->m_Interpolator->IsInsideBuffer(pt))
    {
      float oValue = this->m_Interpolator->Evaluate(pt);
      outIt.Set(oValue);
    }

    else
    {
      outIt.Set(0.0);
    }

    ++outIt;
  }
}
} // end namespace itk

#endif
