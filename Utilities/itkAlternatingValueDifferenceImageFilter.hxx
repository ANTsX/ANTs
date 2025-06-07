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
#ifndef __itkAlternatingValueDifferenceImageFilter_hxx
#define __itkAlternatingValueDifferenceImageFilter_hxx

#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::AlternatingValueDifferenceImageFilter()
{
  m_ControlInterpolator = dynamic_cast<InterpolatorType *>(DefaultInterpolatorType::New().GetPointer());
  m_LabelInterpolator = dynamic_cast<InterpolatorType *>(DefaultInterpolatorType::New().GetPointer());
  m_IndexPadding = 1;
  this->DynamicMultiThreadingOff();
}

template <typename TInputImage, typename TOutputImage>
void
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SubtractionDimension: " << m_SubtractionDimension << std::endl;
}

template <typename TInputImage, typename TOutputImage>
void
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::VerifyInputInformation() const
{
  Superclass::VerifyInputInformation();

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
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  // do not call the superclass' implementation of this method since
  // this filter allows the input and output to be of different dimensions

  // get pointers to the input and output
  typename Superclass::OutputImagePointer     outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer inputPtr = this->GetInput();

  if (!outputPtr || !inputPtr)
  {
    return;
  }

  // Set the output image largest possible region.  Use a RegionCopier
  // so that the input and output images can be different dimensions.
  OutputImageRegionType outputLargestPossibleRegion;
  this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion, inputPtr->GetLargestPossibleRegion());

  // for the subtraction dimension
  unsigned int nPairs =
    this->GetInput()->GetLargestPossibleRegion().GetSize()[m_SubtractionDimension] - 2 * m_IndexPadding;
  outputLargestPossibleRegion.SetSize(m_SubtractionDimension, nPairs);

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
    outputOrigin[m_SubtractionDimension] += outputSpacing[m_SubtractionDimension] * m_IndexPadding;

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
    itkExceptionMacro(<< "itk::DimensionSimpleSubtractionImageFilter::GenerateOutputInformation "
                      << "cannot cast input to " << typeid(ImageBase<InputImageDimension> *).name());
  }

  // Support VectorImages by setting number of components on output.
  const unsigned int numComponents = inputPtr->GetNumberOfComponentsPerPixel();
  if (numComponents != outputPtr->GetNumberOfComponentsPerPixel())
  {
    outputPtr->SetNumberOfComponentsPerPixel(numComponents);
  }
}

/*
template <typename TInputImage, typename TOutputImage>
void
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  if( !this->GetOutput() )
    {
    return;
    }
  OutputImageRegionType outputRegion = this->GetOutput()->GetRequestedRegion();
  IndexValueType        begin = outputRegion.GetIndex(InputImageDimension);
  IndexValueType        end = begin + outputRegion.GetSize(InputImageDimension);
  for( IndexValueType idx = 0; idx < this->GetNumberOfIndexedInputs(); ++idx )
    {
    InputImagePointer inputPtr =
      const_cast<InputImageType *>( this->GetInput(idx) );
    if( !inputPtr )
      {
      // Because DataObject::PropagateRequestedRegion() allows only
      // InvalidRequestedRegionError, it's impossible to write simply:
      // itkExceptionMacro(<< "Missing input " << idx);
      InvalidRequestedRegionError e(__FILE__, __LINE__);
      e.SetLocation(ITK_LOCATION);
      e.SetDescription("Missing input.");
      e.SetDataObject( this->GetOutput() );
      throw e;
      }

    InputImageRegionType inputRegion;
    if( begin <= idx && idx < end )
      {
      this->CallCopyOutputRegionToInputRegion(inputRegion, outputRegion);
      }
    else
      {
      // to tell the pipeline that updating this input is unncesseary
      inputRegion = inputPtr->GetBufferedRegion();
      }
    inputPtr->SetRequestedRegion(inputRegion);
    }
}
*/

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage>
void
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  if (!m_ControlInterpolator)
  {
    itkExceptionMacro(<< "Control interpolator not set");
  }
  if (!m_LabelInterpolator)
  {
    itkExceptionMacro(<< "Label interpolator not set");
  }

  // Split image into control and label images
  typename InputImageType::RegionType           region = this->GetInput()->GetLargestPossibleRegion();
  typename InputImageType::RegionType::SizeType size = region.GetSize();
  size[InputImageDimension - 1] = size[InputImageDimension - 1] / 2;
  region.SetSize(size);
  typename InputImageType::SpacingType spacing = this->GetInput()->GetSpacing();
  spacing[InputImageDimension - 1] = spacing[InputImageDimension - 1] * 2.0;
  typename InputImageType::PointType origin = this->GetInput()->GetOrigin();

  this->m_ControlImage = InputImageType::New();
  this->m_ControlImage->SetRegions(region);
  this->m_ControlImage->SetSpacing(spacing);
  this->m_ControlImage->SetOrigin(origin);
  this->m_ControlImage->SetDirection(this->GetInput()->GetDirection());
  this->m_ControlImage->AllocateInitialized();

  origin[InputImageDimension - 1] = origin[InputImageDimension - 1] + spacing[InputImageDimension - 1];

  this->m_LabelImage = InputImageType::New();
  this->m_LabelImage->SetRegions(region);
  this->m_LabelImage->SetSpacing(spacing);
  this->m_LabelImage->SetOrigin(origin);
  this->m_LabelImage->SetDirection(this->GetInput()->GetDirection());
  this->m_LabelImage->AllocateInitialized();


  this->m_ControlOutputImage = InputImageType::New();
  this->m_ControlOutputImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());
  this->m_ControlOutputImage->SetSpacing(this->GetOutput()->GetSpacing());
  this->m_ControlOutputImage->SetOrigin(this->GetOutput()->GetOrigin());
  this->m_ControlOutputImage->SetDirection(this->GetOutput()->GetDirection());
  this->m_ControlOutputImage->AllocateInitialized();

  this->m_LabelOutputImage = InputImageType::New();
  this->m_LabelOutputImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());
  this->m_LabelOutputImage->SetSpacing(this->GetOutput()->GetSpacing());
  this->m_LabelOutputImage->SetOrigin(this->GetOutput()->GetOrigin());
  this->m_LabelOutputImage->SetDirection(this->GetOutput()->GetDirection());
  this->m_LabelOutputImage->AllocateInitialized();

  ImageRegionConstIterator<InputImageType> it(this->GetInput(), this->GetInput()->GetLargestPossibleRegion());

  while (!it.IsAtEnd())
  {
    typename InputImageType::IndexType idx = it.GetIndex();

    unsigned int offset = idx[InputImageDimension - 1] % 2;
    unsigned int tp = idx[InputImageDimension - 1] / 2;
    idx[InputImageDimension - 1] = tp;

    if (offset)
    {
      this->m_LabelImage->SetPixel(idx, it.Value());
    }
    else
    {
      this->m_ControlImage->SetPixel(idx, it.Value());
    }


    ++it;
  }

  // Connect input image to interpolator
  m_ControlInterpolator->SetInputImage(this->m_ControlImage);
  m_LabelInterpolator->SetInputImage(this->m_LabelImage);

  // Connect input image to extrapolator
  // if( !m_Extrapolator.IsNull() )
  //  {
  //  m_Extrapolator->SetInputImage( this->GetInput() );
  //  }
}

template <typename TInputImage, typename TOutputImage>
void
AlternatingValueDifferenceImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  itkDebugMacro(<< "Actually executing");

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  OutputImageRegionType outputRegion = outputRegionForThread;

  InputImageRegionType inputRegion;
  this->CallCopyOutputRegionToInputRegion(inputRegion, outputRegionForThread);

  ImageRegionIterator<OutputImageType> outIt(this->GetOutput(), outputRegion);

  typename InputImageType::PointType pt;
  while (!outIt.IsAtEnd())
  {
    typename InputImageType::IndexType idx = outIt.GetIndex();
    this->GetOutput()->TransformIndexToPhysicalPoint(idx, pt);

    float cValue = this->m_ControlInterpolator->Evaluate(pt);
    float lValue = this->m_LabelInterpolator->Evaluate(pt);

    this->m_ControlOutputImage->SetPixel(outIt.GetIndex(), cValue);
    this->m_LabelOutputImage->SetPixel(outIt.GetIndex(), lValue);

    outIt.Set(cValue - lValue);
    ++outIt;
  }
}
} // end namespace itk

#endif
