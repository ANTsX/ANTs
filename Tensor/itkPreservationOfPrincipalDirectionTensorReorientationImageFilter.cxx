/*==========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_cxx
#define _itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_cxx

#include "antsAllocImage.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"
#include "itkVector.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "TensorFunctions.h"

namespace itk
{
template <typename TTensorImage>
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage>::
PreservationOfPrincipalDirectionTensorReorientationImageFilter()
{
  m_CompositeTransform = nullptr;
}

template <typename TTensorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage>::GenerateData()
{
  // get input and output images
  // FIXME - use buffered region, etc
  auto input = this->GetInput();
  auto output = this->GetOutput();

  output->SetRegions(input->GetLargestPossibleRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetDirection(input->GetDirection());
  output->AllocateInitialized();

  ImageRegionIteratorWithIndex<OutputImageType> outputIt(output, output->GetLargestPossibleRegion());

  using DirectionType = typename TTensorImage::DirectionType;
  const auto directionTransformMatrix = input->GetDirection().GetVnlMatrix();
  const auto directionTransformMatrixTranspose = input->GetDirection().GetTranspose();

  using TensorType = typename TTensorImage::PixelType;
  using TensorMatrixType = typename DirectionType::InternalMatrixType;

  // for all voxels
  for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
  {
    const auto & index = outputIt.GetIndex();
    const TensorType inTensor = input->GetPixel(index);

    TensorMatrixType tmpMat; // matrix to hold values for rebasing
    TensorType inTensorPhysical; // inTensor in physical space
    TensorType inTensorReoriented; // inTensor reoriented by transform
    TensorType outTensor; // reoriented tensor rebased to voxel space

    // valid values?
    bool hasNans = false;
    for (unsigned int jj = 0; jj < 6; jj++)
    {
      if (std::isnan(inTensor[jj]) || std::isinf(inTensor[jj]))
      {
        hasNans = true;
        break;
      }
    }

    bool isNull = false;
    RealType trace = inTensor[0] + inTensor[3] + inTensor[5];
    if (trace <= 0.0)
    {
      isNull = true;
    }

    if (hasNans || isNull)
    {
      outTensor = inTensor;
    }
    else
    {
      // rebase tensor to physical space using directionTransformMatrix and directionTransformMatrixTranspose;
      tmpMat = Vector2Matrix<TensorType, TensorMatrixType>(inTensor);
      tmpMat = directionTransformMatrix * tmpMat * directionTransformMatrixTranspose;
      inTensorPhysical = Matrix2Vector<TensorType, TensorMatrixType>(tmpMat);

      typename InputImageType::PointType pt;
      input->TransformIndexToPhysicalPoint(index, pt);
      inTensorReoriented = this->m_CompositeTransform->TransformDiffusionTensor3D(inTensorPhysical, pt);
    }

    // valid values?
    for (unsigned int jj = 0; jj < 6; jj++)
    {
      if (std::isnan(inTensorReoriented[jj]) || std::isinf(inTensorReoriented[jj]))
      {
        inTensorReoriented[jj] = 0;
      }
    }

    // rebase tensor to voxel space
    tmpMat = Vector2Matrix<TensorType, TensorMatrixType>(inTensorReoriented);
    tmpMat = directionTransformMatrixTranspose * tmpMat * directionTransformMatrix;
    outTensor = Matrix2Vector<TensorType, TensorMatrixType>(tmpMat);

    outputIt.Set(outTensor);
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TTensorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage>::PrintSelf(
  std::ostream & os,
  Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
