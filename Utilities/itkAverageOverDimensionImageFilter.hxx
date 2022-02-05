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
#ifndef __itkAverageOverDimensionImageFilter_hxx
#define __itkAverageOverDimensionImageFilter_hxx

#include "itkImageAlgorithm.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
/**
 *
 */
template <typename TInputImage, typename TOutputImage>
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::AverageOverDimensionImageFilter()
  :
#ifdef ITKV3_COMPATIBILITY
  m_DirectionCollapseStrategy(DIRECTIONCOLLAPSETOGUESS)
#else
  m_DirectionCollapseStrategy(DIRECTIONCOLLAPSETOUNKOWN)
#endif
{
  Superclass::InPlaceOff();
  this->DynamicMultiThreadingOff();
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OutputImageRegion: " << m_OutputImageRegion << std::endl;
  os << indent << "DirectionCollapseStrategy: " << m_DirectionCollapseStrategy << std::endl;
}

template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::CallCopyOutputRegionToInputRegion(
  InputImageRegionType &        destRegion,
  const OutputImageRegionType & srcRegion)
{
  typename InputImageRegionType::SizeType  inputSize;
  typename InputImageRegionType::IndexType inputIndex;

  unsigned int countDim = 0;
  for (unsigned int i = 0; i < InputImageDimension; i++)
  {
    if (i != this->m_AveragingDimension)
    {
      inputSize[i] = srcRegion.GetSize()[countDim];
      inputIndex[i] = srcRegion.GetIndex()[countDim];
      countDim++;
    }
  }

  inputSize[this->m_AveragingDimension] = this->GetInput()->GetRequestedRegion().GetSize()[this->m_AveragingDimension];
  inputIndex[this->m_AveragingDimension] =
    this->GetInput()->GetRequestedRegion().GetIndex()[this->m_AveragingDimension];

  destRegion.SetSize(inputSize);
  destRegion.SetIndex(inputIndex);
}

template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::SetAveragingDimension(unsigned int averagingDimension)
{
  this->m_AveragingDimension = averagingDimension;

  if (this->m_AveragingDimension > InputImageDimension)
  {
    itkExceptionMacro("Averaging dimension is larger than input image dimension");
  }

  InputImageSizeType  inputSize = this->GetInput()->GetRequestedRegion().GetSize();
  OutputImageSizeType outputSize;
  outputSize.Fill(0);
  OutputImageIndexType outputIndex;
  outputIndex.Fill(0);

  /**
   * check to see if the number of non-zero entries in the extraction region
   * matches the number of dimensions in the output image.
   */
  unsigned int dimCount = 0;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    if (i != this->m_AveragingDimension)
    {
      outputSize[dimCount] = inputSize[i];
      outputIndex[dimCount] = this->GetInput()->GetRequestedRegion().GetIndex()[i];
      dimCount++;
    }
  }

  m_OutputImageRegion.SetSize(outputSize);
  m_OutputImageRegion.SetIndex(outputIndex);
  this->Modified();
}

/**
 * AverageOverDimensionImageFilter can produce an image which is a different resolution
 * than its input image.  As such, AverageOverDimensionImageFilter needs to provide an
 * implementation for GenerateOutputInformation() in order to inform
 * the pipeline execution model.  The original documentation of this
 * method is below.
 *
 * \sa ProcessObject::GenerateOutputInformaton()
 */
template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  // do not call the superclass' implementation of this method since
  // this filter allows the input and the output to be of different dimensions

  // get pointers to the input and output
  typename Superclass::OutputImagePointer     outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer inputPtr = this->GetInput();

  if (!outputPtr || !inputPtr)
  {
    return;
  }

  // Set the output image size to the same value as the extraction region.
  outputPtr->SetLargestPossibleRegion(m_OutputImageRegion);

  // Set the output spacing and origin
  const ImageBase<InputImageDimension> * phyData;

  phyData = dynamic_cast<const ImageBase<InputImageDimension> *>(this->GetInput());

  if (phyData)
  {
    // Copy what we can from the image from spacing and origin of the input
    // This logic needs to be augmented with logic that select which
    // dimensions to copy

    const typename InputImageType::SpacingType &   inputSpacing = inputPtr->GetSpacing();
    const typename InputImageType::DirectionType & inputDirection = inputPtr->GetDirection();
    const typename InputImageType::PointType &     inputOrigin = inputPtr->GetOrigin();

    typename OutputImageType::SpacingType   outputSpacing;
    typename OutputImageType::DirectionType outputDirection;
    typename OutputImageType::PointType     outputOrigin;
    outputOrigin.Fill(0.0);
    outputDirection.SetIdentity();

    unsigned int countDim = 0;
    for (unsigned int i = 0; i < InputImageDimension; ++i)
    {
      if (i != this->m_AveragingDimension)
      {
        outputSpacing[countDim] = inputSpacing[i];
        outputOrigin[countDim] = inputOrigin[i];
        unsigned int countDim2 = 0;
        for (unsigned int dim = 0; dim < InputImageDimension; ++dim)
        {
          if (dim != this->m_AveragingDimension)
          {
            outputDirection[countDim][countDim2] = inputDirection[i][dim];
            ++countDim2;
          }
        }
        countDim++;
      }
    }

    // if the filter changes from a higher to a lower dimension, or
    // if, after rebuilding the direction cosines, there's a zero
    // length cosine vector, reset the directions to identity.
    switch (m_DirectionCollapseStrategy)
    {
      case DIRECTIONCOLLAPSETOIDENTITY:
      {
        outputDirection.SetIdentity();
      }
      break;
      case DIRECTIONCOLLAPSETOSUBMATRIX:
      {
        if (itk::Math::FloatAlmostEqual(vnl_determinant(outputDirection.GetVnlMatrix()), 0.0))
        {
          itkExceptionMacro(<< "Invalid submatrix extracted for collapsed direction.");
        }
      }
      break;
      case DIRECTIONCOLLAPSETOGUESS:
      {
        if (itk::Math::FloatAlmostEqual(vnl_determinant(outputDirection.GetVnlMatrix()), 0.0))
        {
          outputDirection.SetIdentity();
        }
      }
      break;
      case DIRECTIONCOLLAPSETOUNKOWN:
      default:
      {
        itkExceptionMacro(
          << "It is required that the strategy for collapsing the direction matrix be explicitly specified. "
          << "Set with either myfilter->SetDirectionCollapseToIdentity() or "
             "myfilter->SetDirectionCollapseToSubmatrix() "
          << typeid(ImageBase<InputImageDimension> *).name());
      }
    }

    // set the spacing and origin
    outputPtr->SetSpacing(outputSpacing);
    outputPtr->SetDirection(outputDirection);
    outputPtr->SetOrigin(outputOrigin);
    outputPtr->SetNumberOfComponentsPerPixel(inputPtr->GetNumberOfComponentsPerPixel());
  }
  else
  {
    // pointer could not be cast back down
    itkExceptionMacro(<< "itk::AverageOverDimensionImageFilter::GenerateOutputInformation "
                      << "cannot cast input to " << typeid(ImageBase<InputImageDimension> *).name());
  }
}

template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::GenerateData()
{

  // InPlace::AllocateOutputs set the running in place ivar.
  // This method will be called again, by GenerateData, but there is
  // no harm done.
  this->AllocateOutputs();

  // The input matched the output, nothing to do.
  if (this->GetRunningInPlace())
  {
    OutputImageType * outputPtr = this->GetOutput();

    // the in-place grafting copies the meta data, this needs to be
    // set back.
    outputPtr->SetLargestPossibleRegion(m_OutputImageRegion);

    this->UpdateProgress(1.0);
    return;
  }

  this->Superclass::GenerateData();
}

/**
 * AverageOverDimensionImageFilter can be implemented as a multithreaded filter.
 * Therefore, this implementation provides a ThreadedGenerateData()
 * routine which is called for each processing thread. The output
 * image data is allocated automatically by the superclass prior to
 * calling ThreadedGenerateData().  ThreadedGenerateData can only
 * write to the portion of the output image specified by the
 * parameter "outputRegionForThread"
 *
 * \sa ImageToImageFilter::ThreadedGenerateData(),
 *     ImageToImageFilter::GenerateData()
 */
template <typename TInputImage, typename TOutputImage>
void
AverageOverDimensionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  itkDebugMacro(<< "Actually executing");

  // Get the input and output pointers
  const InputImageType * inputPtr = this->GetInput();
  OutputImageType *      outputPtr = this->GetOutput();

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, 1);

  // Define the portion of the input to walk for this thread
  InputImageRegionType inputRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

  unsigned int nValues = inputRegionForThread.GetSize()[this->m_AveragingDimension];
  unsigned int offset = inputRegionForThread.GetIndex()[this->m_AveragingDimension];

  ImageRegionIteratorWithIndex<OutputImageType> it(outputPtr, outputRegionForThread);

  while (!it.IsAtEnd())
  {
    typename InputImageType::IndexType idx;
    unsigned int                       dimCount = 0;
    for (unsigned int i = 0; i < InputImageDimension; i++)
    {
      if (i != this->m_AveragingDimension)
      {
        idx[i] = it.GetIndex()[dimCount];
        ++dimCount;
      }
    }

    typename OutputImageType::PixelType value = 0.0;
    typename OutputImageType::PixelType nTimes = (nValues - offset + 1);
    for (unsigned int i = offset; i < nValues; i++)
    {
      idx[this->m_AveragingDimension] = i;
      value += inputPtr->GetPixel(idx);
    }
    if (nTimes > 0)
    {
      value /= nTimes;
    }
    it.Set(value);

    ++it;
  }

  progress.CompletedPixel();
}
} // end namespace itk

#endif
