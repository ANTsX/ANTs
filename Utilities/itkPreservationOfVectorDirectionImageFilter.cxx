/*==========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPreservationOfVectorDirectionImageFilter_cxx
#define _itkPreservationOfVectorDirectionImageFilter_cxx

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
#include "itkPreservationOfVectorDirectionImageFilter.h"

namespace itk
{
template <typename TVectorImage>
PreservationOfVectorDirectionImageFilter<TVectorImage>::
PreservationOfVectorDirectionImageFilter()
{
  m_CompositeTransform = nullptr;
}

template <typename TVectorImage>
void
PreservationOfVectorDirectionImageFilter<TVectorImage>::GenerateData()
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

  const MatrixType directionTransformMatrix = input->GetDirection().GetVnlMatrix();
  const MatrixType directionTransformMatrixTranspose = input->GetDirection().GetTranspose();

  using VectorType = typename TVectorImage::PixelType;

  auto transformQueue = m_CompositeTransform->GetTransformQueue();

  // for all voxels
  for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
  {
    const auto & index = outputIt.GetIndex();
    const VectorType inVector = input->GetPixel(index);

    // initialize output vector to new VectorType
    VectorType outVector;

    // valid values?
    bool hasNans = false;
    bool isNull = false;
    RealType inVectorNorm = 0.0;
    for (unsigned int jj = 0; jj < ImageDimension; jj++)
    {
      if (std::isnan(inVector[jj]) || std::isinf(inVector[jj]))
      {
        hasNans = true;
        break;
      }
      inVectorNorm += inVector[jj] * inVector[jj];
    }
    inVectorNorm = std::sqrt(inVectorNorm);
    // isNull is true if norm is within machine epsilon of zero
    isNull = (inVectorNorm < NumericTraits<RealType>::epsilon());

    if (hasNans || isNull)
    {
      outVector = inVector;
    }
    else
    {
      // rebase Vector to physical space using directionTransformMatrix
      VectorMatrixType inVectorPhysical;

      for (unsigned int jj = 0; jj < ImageDimension; jj++)
      {
        inVectorPhysical[jj][0] = inVector[jj];
      }

      if (!m_InputVectorsInPhysicalSpace)
      {
        // reorient the vector to physical space
        inVectorPhysical = directionTransformMatrix * inVectorPhysical;
      }

      typename InputImageType::PointType pt;
      input->TransformIndexToPhysicalPoint(index, pt);

      JacobianMatrixType compositeInverseJacobianMatrix;
      compositeInverseJacobianMatrix.fill(0.0);
      for (unsigned int ii = 0; ii < ImageDimension; ++ii)
      {
        compositeInverseJacobianMatrix[ii][ii] = 1.0;
      }

      for (auto it = transformQueue.rbegin(); it != transformQueue.rend(); ++it)
      {
        JacobianMatrixType localJacobian;
        (*it)->ComputeInverseJacobianWithRespectToPosition(pt, localJacobian);
        compositeInverseJacobianMatrix = localJacobian * compositeInverseJacobianMatrix;
        pt = (*it)->TransformPoint(pt);
      }

      // reorient the vector to the new space
      // Get the Jacobian as a local matrix type
      vnl_matrix<double> tmpMat(ImageDimension, ImageDimension);
      for (unsigned int ii = 0; ii < ImageDimension; ++ii)
      {
        for (unsigned int jj = 0; jj < ImageDimension; ++jj)
        {
          tmpMat[ii][jj] = compositeInverseJacobianMatrix[ii][jj];
        }
      }

      // Get a rotation from the composite Jacobian matrix
      vnl_svd<double> svd(tmpMat);
      MatrixType R = svd.U() * svd.V().transpose();

      if (vnl_determinant(R) < 0)
      {
        MatrixType V = svd.V();
        for (unsigned int i = 0; i < ImageDimension; ++i)
        {
          V(i, ImageDimension - 1) *= -1;
        }
        R = svd.U() * V.transpose();
      }

      VectorMatrixType inVectorReoriented = R * inVectorPhysical;

      // rebase Vector to voxel space if neded
      if (!m_InputVectorsInPhysicalSpace)
      {
        // reorient the vector to voxel space
        inVectorReoriented = directionTransformMatrixTranspose * inVectorReoriented;
      }

      for (unsigned int jj = 0; jj < ImageDimension; jj++)
      {
        outVector[jj] = inVectorReoriented[jj][0];
      }

    } // end if hasNans or isNull

    outputIt.Set(outVector);
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TVectorImage>
void
PreservationOfVectorDirectionImageFilter<TVectorImage>::PrintSelf(
  std::ostream & os,
  Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
