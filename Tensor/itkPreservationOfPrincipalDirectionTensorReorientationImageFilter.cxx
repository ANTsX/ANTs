/*=========================================================================

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
template <typename TTensorImage, typename TVectorImage>
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>::
  PreservationOfPrincipalDirectionTensorReorientationImageFilter()
{
  m_DisplacementField = nullptr;
  m_AffineTransform = nullptr;
  m_UseAffine = false;
}


template <typename TTensorImage, typename TVectorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>::GenerateData()
{
  // get input and output images
  // FIXME - use buffered region, etc
  InputImagePointer  input = this->GetInput();
  OutputImagePointer output = this->GetOutput();

  if (this->m_UseAffine)
  {
    output->SetRegions(input->GetLargestPossibleRegion());
    output->SetSpacing(input->GetSpacing());
    output->SetOrigin(input->GetOrigin());
    output->SetDirection(input->GetDirection());
    output->AllocateInitialized();
  }
  else
  {
    output->SetRegions(input->GetLargestPossibleRegion());
    output->SetSpacing(input->GetSpacing());
    output->SetOrigin(input->GetOrigin());
    output->SetDirection(input->GetDirection());
    output->AllocateInitialized();

    this->m_DisplacementTransform = DisplacementFieldTransformType::New();
    this->m_DisplacementTransform->SetDisplacementField(m_DisplacementField);
  }

  ImageRegionIteratorWithIndex<OutputImageType> outputIt(output, output->GetLargestPossibleRegion());

  std::cout << "Iterating over image" << std::endl;

  using DirectionType = typename TTensorImage::DirectionType;

  auto directionTransformMatrix = input->GetDirection().GetVnlMatrix();
  auto directionTransformMatrixTranspose = input->GetDirection().GetTranspose();

  using TensorType = typename TTensorImage::PixelType;
  using TensorMatrixType = typename DirectionType::InternalMatrixType;

  // for all voxels
  for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
  {

    TensorType inTensor = input->GetPixel(outputIt.GetIndex());

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
      }
    }

    bool     isNull = false;
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

      if (this->m_UseAffine)
      {
        Vector2Matrix<TensorType, TensorMatrixType>(inTensor, tmpMat);
        inTensorReoriented = this->m_AffineTransform->TransformDiffusionTensor3D(inTensorPhysical);
      }
      else
      {
        typename DisplacementFieldType::PointType pt;
        this->m_DisplacementField->TransformIndexToPhysicalPoint(outputIt.GetIndex(), pt);
        inTensorReoriented = this->m_DisplacementTransform->TransformDiffusionTensor3D(inTensorPhysical, pt);
      }
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
template <typename TTensorImage, typename TVectorImage>
void
PreservationOfPrincipalDirectionTensorReorientationImageFilter<TTensorImage, TVectorImage>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
