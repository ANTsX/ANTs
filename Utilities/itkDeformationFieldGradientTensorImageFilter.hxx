/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDeformationFieldGradientTensorImageFilter_hxx
#define _itkDeformationFieldGradientTensorImageFilter_hxx


#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"
#include "itkCastImageFilter.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

#include "vnl/vnl_det.h"

namespace itk
{

template <typename TInputImage, typename TRealType, typename TOutputImage>
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::
  DeformationFieldGradientTensorImageFilter()
{
  this->m_CalculateJacobian = false;
  this->m_UseImageSpacing = true;
  this->m_Order = 2;
  this->m_DerivativeWeights.Fill(1.0);
  this->DynamicMultiThreadingOff();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  this->m_NeighborhoodRadius.Fill(this->m_Order);

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius(this->m_NeighborhoodRadius);

  // crop the input requested region at the input's largest possible region
  if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
  {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
  }
  else
  {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::BeforeThreadedGenerateData()
{
  if (this->m_UseImageSpacing)
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      if (itk::Math::FloatAlmostEqual(this->GetInput()->GetSpacing()[i], 0.0))
      {
        itkExceptionMacro(<< "Image spacing in dimension " << i << " is zero.");
      }
      this->m_DerivativeWeights[i] = 1.0 / static_cast<RealType>(this->GetInput()->GetSpacing()[i]);
    }
  }

  /** If the input needs casting to a real-valued vector type, create the
      appropriate image and set the m_RealValuedInputImage pointer to this
      image.  Otherwise just point to the input image. */
  if (typeid(typename InputImageType::PixelType) != typeid(RealVectorType))
  {
    typename CastImageFilter<TInputImage, RealVectorImageType>::Pointer caster =
      CastImageFilter<TInputImage, RealVectorImageType>::New();
    caster->SetInput(this->GetInput());
    caster->Update();
    this->m_RealValuedInputImage = caster->GetOutput();
  }
  else
  {
    this->m_RealValuedInputImage = dynamic_cast<const RealVectorImageType *>(this->GetInput());
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  ZeroFluxNeumannBoundaryCondition<RealVectorImageType> nbc;
  ConstNeighborhoodIteratorType                         bit;
  ImageRegionIterator<TOutputImage>                     it;

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>                        bC;
  faceList = bC(dynamic_cast<const RealVectorImageType *>(this->m_RealValuedInputImage.GetPointer()),
                outputRegionForThread,
                this->m_NeighborhoodRadius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>::FaceListType::iterator fit;

  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  const typename InputImageType::DirectionType R = this->m_RealValuedInputImage->GetDirection();
  // const typename InputImageType::DirectionType RI = this->m_RealValuedInputImage->GetInverseDirection();

  // Process each of the data set faces.  The iterator is reinitialized on each
  // face so that it can determine whether or not to check for boundary
  // conditions.
  for (fit = faceList.begin(); fit != faceList.end(); ++fit)
  {
    bit = ConstNeighborhoodIteratorType(
      this->m_NeighborhoodRadius, dynamic_cast<const RealVectorImageType *>(m_RealValuedInputImage.GetPointer()), *fit);
    it = ImageRegionIterator<TOutputImage>(this->GetOutput(), *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    while (!bit.IsAtEnd())
    {
      RealMatrixType F = this->EvaluateAtNeighborhood(bit);

      // (PAC 2025-05-22): F is the gradient of the deformation field. The displacements are in physical space,
      // but the gradient is estimated in the index coordinates (ie, using a neighborhood on the voxel grid)
      // This next step is from itkDisplacementFieldTransform and transforms the rows from index to physical space.
      // I think this is correct because it produces more accurate volume change estimation in my tests.
      for (unsigned int i = 0; i < ImageDimension; ++i)
      {
        itk::Vector<TRealType, ImageDimension> row_i;
        for (unsigned int j = 0; j < ImageDimension; ++j)
        {
          row_i[j] = F[i][j];
        }
        itk::Vector<TRealType, ImageDimension> rotated_row = R * row_i;
        for (unsigned int j = 0; j < ImageDimension; ++j)
        {
          F[i][j] = rotated_row[j];
        }
      }

      if (this->m_CalculateJacobian)
      {
        for (unsigned int i = 0; i < ImageDimension; i++)
        {
          F[i][i] += 1.0;
        }
      }

      it.Set(F);
      ++bit;
      ++it;
      progress.CompletedPixel();
    }
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::RealMatrixType
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::EvaluateAtNeighborhood(
  const ConstNeighborhoodIteratorType & it) const
{
  unsigned       i, j;
  RealMatrixType F;

  for (j = 0; j < ImageDimension; ++j)
  {
    RealVectorType physicalVectorNext1 = it.GetNext(j, 1);
    RealVectorType physicalVectorNext2 = it.GetNext(j, 2);
    RealVectorType physicalVectorPrevious1 = it.GetPrevious(j, 1);
    RealVectorType physicalVectorPrevious2 = it.GetPrevious(j, 2);

    RealType weight = this->m_DerivativeWeights[j];
    for (i = 0; i < ImageDimension; ++i)
    {
       F[i][j] = weight *
                      (-physicalVectorNext2[i] + 8.0 * physicalVectorNext1[i] - 8.0 * physicalVectorPrevious1[i] +
                       physicalVectorPrevious2[i]) /
                      12.0;
    }
  }
  return F;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeformationFieldGradientTensorImageFilter<TInputImage, TRealType, TOutputImage>::PrintSelf(std::ostream & os,
                                                                                           Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "m_CalculateJacobian = " << this->m_CalculateJacobian << std::endl;
  os << indent << "m_UseImageSpacing = " << this->m_UseImageSpacing << std::endl;
  os << indent << "m_Order = " << this->m_Order << std::endl;
  os << indent << "m_NeighborhoodRadius = " << this->m_NeighborhoodRadius << std::endl;
  os << indent << "m_RealValuedInputImage = " << m_RealValuedInputImage.GetPointer() << std::endl;
}

} // end namespace itk

#endif
