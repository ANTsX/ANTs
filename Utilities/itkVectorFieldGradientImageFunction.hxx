/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkVectorFieldGradientImageFunction_hxx
#define _itkVectorFieldGradientImageFunction_hxx


#include "itkDecomposeTensorFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_vector.h"

namespace itk
{
template <typename TInputImage, typename TRealType, typename TOutputImage>
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::VectorFieldGradientImageFunction()
{}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateDeformationGradientTensor(
  const PointType & point) const
{
  if (!this->IsInsideBuffer(point))
  {
    unsigned int minDimension = ImageDimension;

    MatrixType F;
    F.SetSize(ImageDimension, VectorDimension);
    F.Fill(0.0);
    for (unsigned int i = 0; i < minDimension; i++)
    {
      F[i][i] = 1.0;
    }
    itkWarningMacro("The specified point, " << point << ", is outside the image boundaries.");
    return F;
  }

  typename InputImageType::SpacingType spacing = this->GetInputImage()->GetSpacing();

  typedef VectorLinearInterpolateImageFunction<InputImageType> InterpolatorType;
  typename InterpolatorType::Pointer                           interpolator = InterpolatorType::New();
  interpolator->SetInputImage(this->GetInputImage());

  MatrixType F;
  F.SetSize(ImageDimension, VectorDimension);

  typename InterpolatorType::OutputType x;

  typename InterpolatorType::PointType ipoint;
  ipoint.CastFrom(point);
  x = interpolator->Evaluate(ipoint);
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    typename PointType::VectorType delta;
    delta.Fill(0.0);
    delta[i] = spacing[i];

    typename InterpolatorType::PointType::VectorType idelta;
    idelta.Fill(0.0);
    idelta[i] = spacing[i];

    typename InterpolatorType::OutputType xp1;
    typename InterpolatorType::OutputType xp2;
    typename InterpolatorType::OutputType xm1;
    typename InterpolatorType::OutputType xm2;

    if (!this->IsInsideBuffer(point + delta))
    {
      xp2 = xp1 = x;
    }
    else
    {
      xp1 = interpolator->Evaluate(ipoint + idelta);
      if (this->IsInsideBuffer(point + delta * 2.0))
      {
        xp2 = interpolator->Evaluate(ipoint + idelta * 2.0);
      }
      else
      {
        xp2 = xp1;
      }
    }

    if (!this->IsInsideBuffer(point - delta))
    {
      xm2 = xm1 = x;
    }
    else
    {
      xm1 = interpolator->Evaluate(ipoint - idelta);
      if (this->IsInsideBuffer(point - delta * 2.0))
      {
        xm2 = interpolator->Evaluate(ipoint - idelta * 2.0);
      }
      else
      {
        xm2 = xm1;
      }
    }

    RealType weight = itk::NumericTraits<RealType>::OneValue() / (static_cast<RealType>(12.0) * delta[i]);
    RealType eightValue = static_cast<RealType>(8.0);
    for (unsigned int j = 0; j < VectorDimension; j++)
    {
      F[i][j] = static_cast<typename MatrixType::ComponentType>(
        weight * (static_cast<RealType>(-xp2[j]) + eightValue * static_cast<RealType>(xp1[j]) -
                  eightValue * static_cast<RealType>(xm1[j]) + static_cast<RealType>(xm2[j])));
    }
  }

  unsigned int minDimension = ImageDimension;
  if (static_cast<unsigned int>(VectorDimension) < static_cast<unsigned int>(ImageDimension))
  {
    minDimension = VectorDimension;
  }
  for (unsigned int i = 0; i < minDimension; i++)
  {
    F[i][i] += itk::NumericTraits<RealType>::OneValue();
  }

  return F;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateDeformationGradientTensorAtIndex(
  const IndexType & index) const
{
  if (!this->IsInsideBuffer(index))
  {
    unsigned int minDimension = ImageDimension;

    MatrixType F;
    F.SetSize(ImageDimension, VectorDimension);
    F.Fill(0.0);
    for (unsigned int i = 0; i < minDimension; i++)
    {
      F[i][i] = 1.0;
    }
    itkWarningMacro("The specified index, " << index << ", is outside the image boundaries.");
    return F;
  }

  MatrixType F;
  F.SetSize(ImageDimension, VectorDimension);

  typename InputImageType::PixelType x;
  x = this->GetInputImage()->GetPixel(index);

  typename InputImageType::SpacingType spacing = this->GetInputImage()->GetSpacing();
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    typename InputImageType::OffsetType offset1;
    offset1.Fill(0);
    offset1[i] = 1;
    typename InputImageType::OffsetType offset2;
    offset2.Fill(0);
    offset2[i] = 2;

    typename InputImageType::PixelType xp1;
    typename InputImageType::PixelType xp2;
    typename InputImageType::PixelType xm1;
    typename InputImageType::PixelType xm2;

    if (!this->IsInsideBuffer(index + offset1))
    {
      xp2 = xp1 = x;
    }
    else
    {
      xp1 = this->GetInputImage()->GetPixel(index + offset1);
      if (this->IsInsideBuffer(index + offset2))
      {
        xp2 = this->GetInputImage()->GetPixel(index + offset2);
      }
      else
      {
        xp2 = xp1;
      }
    }

    if (!this->IsInsideBuffer(index - offset1))
    {
      xm2 = xm1 = x;
    }
    else
    {
      xm1 = this->GetInputImage()->GetPixel(index - offset1);
      if (this->IsInsideBuffer(index - offset2))
      {
        xm2 = this->GetInputImage()->GetPixel(index - offset2);
      }
      else
      {
        xm2 = xm1;
      }
    }

    RealType weight =
      NumericTraits<RealType>::OneValue() /
      (static_cast<RealType>(12.0) * static_cast<RealType>(offset1[i]) * static_cast<RealType>(spacing[i]));
    for (unsigned int j = 0; j < VectorDimension; j++)
    {
      F[i][j] = weight * (-xp2[j] + static_cast<RealType>(8.0) * xp1[j] - static_cast<RealType>(8.0) * xm1[j] + xm2[j]);
    }
  }

  unsigned int minDimension = ImageDimension;
  if (static_cast<unsigned int>(VectorDimension) < static_cast<unsigned int>(ImageDimension))
  {
    minDimension = VectorDimension;
  }
  for (unsigned int i = 0; i < minDimension; i++)
  {
    F[i][i] += itk::NumericTraits<RealType>::OneValue();
  }

  return F;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateJacobian(const PointType & point) const
{
  return this->EvaluateDeformationGradientTensor(point);
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateJacobianAtIndex(
  const IndexType & index) const
{
  return this->EvaluateDeformationGradientTensorAtIndex(index);
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateJacobianDeterminant(
  const PointType & point) const
{
  MatrixType J = this->EvaluateJacobian(point);

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  return decomposer->EvaluateDeterminant(J);
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateJacobianDeterminantAtIndex(
  const IndexType & index) const
{
  MatrixType J = this->EvaluateJacobianAtIndex(index);

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  return decomposer->EvaluateDeterminant(J);
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLagrangianStrainTensor(
  const PointType & point) const
{
  MatrixType F = this->EvaluateDeformationGradientTensor(point);
  MatrixType E(ImageDimension, ImageDimension);

  typename MatrixType::InternalMatrixType ff = F.GetTranspose() * F.GetVnlMatrix();
  for (unsigned int i = 0; i < ff.rows(); i++)
  {
    for (unsigned int j = 0; j < ff.columns(); j++)
    {
      E[i][j] = ff.get(i, j);
      if (i == j)
      {
        E[i][j] -= 1.0;
      }
      E[i][j] *= 0.5;
    }
  }
  return E;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLagrangianStrainTensorAtIndex(
  const IndexType & index) const
{
  MatrixType F = this->EvaluateDeformationGradientTensorAtIndex(index);
  MatrixType E(ImageDimension, ImageDimension);

  typename MatrixType::InternalMatrixType ff = F.GetTranspose() * F.GetVnlMatrix();
  for (unsigned int i = 0; i < ff.rows(); i++)
  {
    for (unsigned int j = 0; j < ff.columns(); j++)
    {
      E[i][j] = ff.get(i, j);
      if (i == j)
      {
        E[i][j] -= 1.0;
      }
      E[i][j] *= 0.5;
    }
  }
  return E;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLagrangianDirectionalStrain(
  const PointType &  point,
  const VectorType & V) const
{
  MatrixType E = this->EvaluateLagrangianStrainTensor(point);

  vnl_vector<RealType> v = E * V.GetVnlVector();
  RealType             s = 0.0;
  for (unsigned int i = 0; i < v.size(); i++)
  {
    s += v[i] * V[i];
  }
  return s;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLagrangianDirectionalStrainAtIndex(
  const IndexType &  index,
  const VectorType & V) const
{
  MatrixType E = this->EvaluateLagrangianStrainTensorAtIndex(index);

  vnl_vector<RealType> v = E * V.GetVnlVector();
  RealType             s = 0.0;
  for (unsigned int i = 0; i < v.size(); i++)
  {
    s += v[i] * V[i];
  }
  return s;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateEulerianStrainTensor(
  const PointType & point) const
{
  MatrixType F = this->EvaluateDeformationGradientTensor(point);
  MatrixType E(ImageDimension, ImageDimension);

  typename MatrixType::InternalMatrixType ff = vnl_matrix_inverse<RealType>(F.GetVnlMatrix() * F.GetTranspose());
  for (unsigned int i = 0; i < ff.rows(); i++)
  {
    for (unsigned int j = 0; j < ff.columns(); j++)
    {
      E[i][j] = -ff.get(i, j);
      if (i == j)
      {
        E[i][j] += 1.0;
      }
      E[i][j] *= 0.5;
    }
  }
  return E;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateEulerianStrainTensorAtIndex(
  const IndexType & index) const
{
  MatrixType F = this->EvaluateDeformationGradientTensorAtIndex(index);
  MatrixType E(ImageDimension, ImageDimension);

  typename MatrixType::InternalMatrixType ff =
    vnl_matrix_inverse<RealType>(F.GetVnlMatrix() * F.GetTranspose()).as_matrix();
  for (unsigned int i = 0; i < ff.rows(); i++)
  {
    for (unsigned int j = 0; j < ff.columns(); j++)
    {
      E[i][j] = -ff.get(i, j);
      if (i == j)
      {
        E[i][j] += 1.0;
      }
      E[i][j] *= 0.5;
    }
  }
  return E;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateEulerianDirectionalStrain(
  const PointType &  point,
  const VectorType & V) const
{
  MatrixType E = this->EvaluateEulerianStrainTensor(point);

  vnl_vector<RealType> v = E * V.GetVnlVector();
  RealType             s = 0.0;
  for (unsigned int i = 0; i < v.size(); i++)
  {
    s += v[i] * V[i];
  }
  return s;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::RealType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateEulerianDirectionalStrainAtIndex(
  const IndexType &  index,
  const VectorType & V) const
{
  MatrixType E = this->EvaluateEulerianStrainTensorAtIndex(index);

  vnl_vector<RealType> v = E * V.GetVnlVector();
  RealType             s = 0.0;
  for (unsigned int i = 0; i < v.size(); i++)
  {
    s += v[i] * V[i];
  }
  return s;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLeftCauchyGreenDeformationTensor(
  const PointType & point) const
{
  MatrixType F = this->EvaluateDeformationGradientTensor(point);

  typename MatrixType::InternalMatrixType b = F.GetVnlMatrix() * F.GetTranspose();
  MatrixType                              B;
  B.SetSize(F.Rows(), F.Cols());
  for (unsigned int i = 0; i < F.Rows(); i++)
  {
    for (unsigned int j = 0; j < F.Cols(); j++)
    {
      B[i][j] = b[i][j];
    }
  }
  return B;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLeftCauchyGreenDeformationTensorAtIndex(
  const IndexType & index) const
{
  MatrixType F = this->EvaluateDeformationGradientTensorAtIndex(index);

  typename MatrixType::InternalMatrixType b = F.GetVnlMatrix() * F.GetTranspose();
  MatrixType                              B;
  B.SetSize(F.Rows(), F.Cols());
  for (unsigned int i = 0; i < F.Rows(); i++)
  {
    for (unsigned int j = 0; j < F.Cols(); j++)
    {
      B[i][j] = b[i][j];
    }
  }
  return B;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateRightCauchyGreenDeformationTensor(
  const PointType & point) const
{
  MatrixType F = this->EvaluateDeformationGradientTensor(point);

  typename MatrixType::InternalMatrixType c = F.GetVnlMatrix() * F.GetTranspose();
  MatrixType                              C;
  C.SetSize(F.Rows(), F.Cols());
  for (unsigned int i = 0; i < F.Rows(); i++)
  {
    for (unsigned int j = 0; j < F.Cols(); j++)
    {
      C[i][j] = c[i][j];
    }
  }
  return C;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::
  EvaluateRightCauchyGreenDeformationTensorAtIndex(const IndexType & index) const
{
  MatrixType F = this->EvaluateDeformationGradientTensorAtIndex(index);

  typename MatrixType::InternalMatrixType c = F.GetVnlMatrix() * F.GetTranspose();
  MatrixType                              C;
  C.SetSize(F.Rows(), F.Cols());
  for (unsigned int i = 0; i < F.Rows(); i++)
  {
    for (unsigned int j = 0; j < F.Cols(); j++)
    {
      C[i][j] = c[i][j];
    }
  }
  return C;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateRightStretchTensor(
  const PointType & point) const
{
  MatrixType C = this->EvaluateRightCauchyGreenDeformationTensor(point);
  MatrixType D;
  MatrixType V;

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(C, D, V);

  MatrixType U;
  U.SetSize(C.Rows(), C.Cols());
  U.Fill(0);
  for (unsigned int d = 0; d < C.Rows(); d++)
  {
    RealType lambda = sqrt(D[d][d]);
    for (unsigned int i = 0; i < C.Rows(); i++)
    {
      for (unsigned int j = 0; j < C.Cols(); j++)
      {
        U[i][j] += lambda * V[i][d] * V[j][d];
      }
    }
  }
  return U;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateRightStretchTensorAtIndex(
  const IndexType & index) const
{
  MatrixType C = this->EvaluateRightCauchyGreenDeformationTensorAtIndex(index);
  MatrixType D;
  MatrixType V;

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(C, D, V);

  MatrixType U;
  U.SetSize(C.Rows(), C.Cols());
  U.Fill(0);
  for (unsigned int d = 0; d < C.Rows(); d++)
  {
    RealType lambda = sqrt(D[d][d]);
    for (unsigned int i = 0; i < C.Rows(); i++)
    {
      for (unsigned int j = 0; j < C.Cols(); j++)
      {
        U[i][j] += lambda * V[i][d] * V[j][d];
      }
    }
  }
  return U;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLeftStretchTensor(
  const PointType & point) const
{
  MatrixType B = this->EvaluateLeftCauchyGreenDeformationTensor(point);
  MatrixType D;
  MatrixType V;

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(B, D, V);

  MatrixType U;
  U.SetSize(B.Rows(), B.Cols());
  U.Fill(0);
  for (unsigned int d = 0; d < B.Rows(); d++)
  {
    RealType lambda = sqrt(D[d][d]);
    for (unsigned int i = 0; i < B.Rows(); i++)
    {
      for (unsigned int j = 0; j < B.Cols(); j++)
      {
        U[i][j] += lambda * V[i][d] * V[j][d];
      }
    }
  }
  return U;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::MatrixType
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::EvaluateLeftStretchTensorAtIndex(
  const IndexType & index) const
{
  MatrixType B = this->EvaluateLeftCauchyGreenDeformationTensorAtIndex(index);
  MatrixType D;
  MatrixType V;

  typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
  typename DecomposerType::Pointer                      decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(B, D, V);

  MatrixType U;
  U.SetSize(B.Rows(), B.Cols());
  U.Fill(0);
  for (unsigned int d = 0; d < B.Rows(); d++)
  {
    RealType lambda = sqrt(D[d][d]);
    for (unsigned int i = 0; i < B.Rows(); i++)
    {
      for (unsigned int j = 0; j < B.Cols(); j++)
      {
        U[i][j] += lambda * V[i][d] * V[j][d];
      }
    }
  }
  return U;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
VectorFieldGradientImageFunction<TInputImage, TRealType, TOutputImage>::PrintSelf(std::ostream & os,
                                                                                  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
