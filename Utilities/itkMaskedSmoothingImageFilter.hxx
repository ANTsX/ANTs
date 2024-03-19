/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskedSmoothingImageFilter_hxx
#define __itkMaskedSmoothingImageFilter_hxx


#include "itkCastImageFilter.h"

namespace itk
{
template <typename TInputImage, typename TMaskImage, typename TOutputImage>
MaskedSmoothingImageFilter<TInputImage, TMaskImage, TOutputImage>::MaskedSmoothingImageFilter()
  : m_SparseImageNeighborhoodRadius(2)
  , m_SmoothingVariance(1.0)
  , m_TimeSmoothingVariance(1.0)
{
  this->SetNumberOfRequiredInputs(2);

  this->m_SparseMatrixIndexImage = nullptr;
  this->m_TimePoints.clear();
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
MaskedSmoothingImageFilter<TInputImage, TMaskImage, TOutputImage>::~MaskedSmoothingImageFilter() = default;

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MaskedSmoothingImageFilter<TInputImage, TMaskImage, TOutputImage>::GenerateData()
{
  this->AllocateOutputs();

  const InputImageType * inputImage = this->GetInput();
  const MaskImageType *  maskImage = this->GetMaskImage();

  // Check if this is longitudinal data
  if (this->m_TimePoints.size() > 0 &&
      inputImage->GetLargestPossibleRegion().GetSize()[ImageDimension - 1] != this->m_TimePoints.size())
  {
    itkExceptionMacro("Longitudinal data does not match expected size.");
  }

  // Initialize the sparse matrix smoother.   Iterate over the input
  // mask and count the entries.  Label an image with the numeric index to
  // make a sparse matrix of that size and fill it with neighborhood data.

  this->m_SparseMatrixIndexImage = RealImageType::New();
  this->m_SparseMatrixIndexImage->CopyInformation(inputImage);
  this->m_SparseMatrixIndexImage->SetRegions(inputImage->GetRequestedRegion());
  this->m_SparseMatrixIndexImage->AllocateInitialized();

  ImageRegionIteratorWithIndex<RealImageType> ItSparseImage(this->m_SparseMatrixIndexImage,
                                                            inputImage->GetRequestedRegion());

  SizeValueType count = 0;
  for (ItSparseImage.GoToBegin(); !ItSparseImage.IsAtEnd(); ++ItSparseImage)
  {
    if (maskImage->GetPixel(ItSparseImage.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue())
    {
      ItSparseImage.Set(count + 1);
      ++count;
    }
  }
  this->m_SparseMatrix.set_size(count, count);

  typename RealImageType::SpacingType spacing = inputImage->GetSpacing();
  typename RealImageType::SizeType    sparseRadius;
  sparseRadius.Fill(this->m_SparseImageNeighborhoodRadius);

  SizeValueType numberOfNeighborhoodIndices =
    static_cast<unsigned int>(std::pow(this->m_SparseImageNeighborhoodRadius * 2 + 1, ImageDimension));

  NeighborhoodIterator<RealImageType> ItSparseImageNeighborhood(
    sparseRadius, this->m_SparseMatrixIndexImage, inputImage->GetRequestedRegion());
  for (ItSparseImageNeighborhood.GoToBegin(); !ItSparseImageNeighborhood.IsAtEnd(); ++ItSparseImageNeighborhood)
  {
    typename RealImageType::IndexType centerIndex = ItSparseImageNeighborhood.GetIndex();

    if (maskImage->GetPixel(centerIndex) != NumericTraits<MaskPixelType>::ZeroValue())
    {
      long location =
        static_cast<long>(this->m_SparseMatrixIndexImage->GetPixel(centerIndex) + static_cast<RealType>(0.5)) - 1;

      for (SizeValueType i = 0; i < numberOfNeighborhoodIndices; i++)
      {
        bool isInBounds;
        ItSparseImageNeighborhood.GetPixel(i, isInBounds);
        if (isInBounds)
        {
          IndexType index = ItSparseImageNeighborhood.GetIndex(i);

          if (maskImage->GetPixel(index) != NumericTraits<MaskPixelType>::ZeroValue())
          {
            long next =
              static_cast<long>(ItSparseImageNeighborhood.GetPixel(i, isInBounds) + static_cast<RealType>(0.5)) - 1;

            if (next >= 0 && location >= 0 && isInBounds)
            {
              RealType spaceValue = NumericTraits<RealType>::ZeroValue();
              RealType timeValue = NumericTraits<RealType>::ZeroValue();
              if (this->m_TimePoints.size() == 0)
              {
                for (SizeValueType k = 0; k < ImageDimension; k++)
                {
                  spaceValue += static_cast<RealType>(std::pow((centerIndex[k] - index[k]) * spacing[k], 2.0));
                }
              }
              else // handle temporal regularization
              {
                spaceValue = NumericTraits<RealType>::ZeroValue();
                for (unsigned int k = 0; k < ImageDimension - 1; k++)
                {
                  spaceValue += static_cast<RealType>(std::pow((centerIndex[k] - index[k]) * spacing[k], 2.0));
                }
                RealType timeDistance =
                  this->m_TimePoints[centerIndex[ImageDimension - 1]] - this->m_TimePoints[index[ImageDimension - 1]];

                timeValue = static_cast<RealType>(std::pow(timeDistance, 2.0));
              }
              RealType smoothValue =
                std::exp(static_cast<RealType>(-1.0) *
                         (spaceValue / this->m_SmoothingVariance + timeValue / this->m_TimeSmoothingVariance));
              this->m_SparseMatrix(location, next) = smoothValue;
            }
          }
        }
      }
    }
  }

  // We use this for smoothing so we force rows to sum to one.
  for (unsigned int k = 0; k < this->m_SparseMatrix.rows(); k++)
  {
    RealType rowSum = this->m_SparseMatrix.sum_row(k);
    if (rowSum > 0)
    {
      this->m_SparseMatrix = this->m_SparseMatrix.scale_row(k, NumericTraits<RealType>::OneValue() / rowSum);
    }
  }

  using CasterType = CastImageFilter<InputImageType, OutputImageType>;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput(inputImage);
  caster->Update();

  OutputImagePointer outputImage = caster->GetOutput();
  outputImage->Update();
  outputImage->DisconnectPipeline();

  // 0. use the fact that we already have a properly sized sparse matrix
  // 1. build a vnl matrix that holds the current input image components
  // 2. sparseMatrix * displacementMatrix
  // 3. put the result in the outputImage

  SizeValueType        numberOfMaskedVoxels = this->m_SparseMatrix.rows();
  vnl_matrix<RealType> smoothImageMatrix(numberOfMaskedVoxels, outputImage->GetNumberOfComponentsPerPixel());
  smoothImageMatrix.fill(0);

  ImageRegionIteratorWithIndex<OutputImageType> It(outputImage, outputImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    IndexType index = It.GetIndex();
    if (maskImage->GetPixel(index) != NumericTraits<MaskPixelType>::ZeroValue())
    {
      long location =
        static_cast<long>(this->m_SparseMatrixIndexImage->GetPixel(index) + static_cast<RealType>(0.5)) - 1;

      InputPixelType inputValue = It.Get();
      for (SizeValueType d = 0; d < smoothImageMatrix.columns(); d++)
      {
        smoothImageMatrix(location, d) = inputValue[d];
      }
    }
  }

  // perform the smoothing operation
  for (SizeValueType d = 0; d < ImageDimension; d++)
  {
    vnl_vector<RealType> smoothColumn;
    this->m_SparseMatrix.mult(smoothImageMatrix.get_column(d), smoothColumn);
    smoothImageMatrix.set_column(d, smoothColumn);
  }

  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (maskImage->GetPixel(It.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue())
    {
      long location =
        static_cast<long>(this->m_SparseMatrixIndexImage->GetPixel(It.GetIndex()) + static_cast<RealType>(0.5)) - 1;
      OutputPixelType outputValue = It.Get();
      for (SizeValueType d = 0; d < smoothImageMatrix.columns(); d++)
      {
        outputValue[d] = smoothImageMatrix(location, d);
      }
      It.Set(outputValue);
    }
  }

  this->SetNthOutput(0, outputImage);
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MaskedSmoothingImageFilter<TInputImage, TMaskImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Smoothing variance = " << this->m_SmoothingVariance << std::endl;
  os << indent << "Time smoothing variance = " << this->m_TimeSmoothingVariance << std::endl;
  os << indent << "Sparse neighborhood radius == " << this->m_SparseImageNeighborhoodRadius << std::endl;
}

} // end namespace itk

#endif
