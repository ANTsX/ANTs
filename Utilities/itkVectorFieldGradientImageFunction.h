/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorFieldGradientImageFunction_h
#define __itkVectorFieldGradientImageFunction_h

#include "itkImageFunction.h"
#include "itkVariableSizeMatrix.h"

namespace itk
{
/** \class VectorFieldGradientImageFunction
 *
 */
template <typename TInputImage, typename TRealType = float, typename TOutput = itk::VariableSizeMatrix<TRealType>>
class VectorFieldGradientImageFunction final : public ImageFunction<TInputImage, TOutput>
{
public:
  /** Standard class typedefs. */
  typedef VectorFieldGradientImageFunction    Self;
  typedef ImageFunction<TInputImage, TOutput> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(VectorFieldGradientImageFunction);

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef TInputImage                        InputImageType;
  typedef typename InputImageType::PixelType VectorType;
  typedef TOutput                            MatrixType;

  /** Declare typedefs of superclass */
  typedef typename Superclass::PointType           PointType;
  typedef typename Superclass::IndexType           IndexType;
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** The dimensionality of the input and output images. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Length of the vector pixel type of the input image. */
  static constexpr unsigned int VectorDimension = VectorType::Dimension;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType RealType;

  /**
   * Evaluate deformation gradient tensor
   */
  MatrixType
  EvaluateDeformationGradientTensor(const PointType &) const;

  MatrixType
  EvaluateDeformationGradientTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateDeformationGradientTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateDeformationGradientTensor(point);
  }

  /**
   * Evaluate Jacobian
   */
  MatrixType
  Evaluate(const PointType & point) const override
  {
    return this->EvaluateJacobian(point);
  }

  MatrixType
  EvaluateAtIndex(const IndexType & idx) const override
  {
    return this->EvaluateJacobianAtIndex(idx);
  }

  MatrixType
  EvaluateAtContinuousIndex(const ContinuousIndexType & idx) const override
  {
    return this->EvaluateJacobianAtContinuousIndex(idx);
  }

  MatrixType
  EvaluateJacobian(const PointType &) const;

  MatrixType
  EvaluateJacobianAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateJacobianAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateJacobian(point);
  }

  /**
   * Evaluate Jacobian determinant
   */
  RealType
  EvaluateJacobianDeterminant(const PointType &) const;

  RealType
  EvaluateJacobianDeterminantAtIndex(const IndexType & idx) const;

  RealType
  EvaluateJacobianDeterminantAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateJacobianDeterminant(point);
  }

  /**
   * Evaluate Lagrangian strain tensor
   */
  MatrixType
  EvaluateLagrangianStrainTensor(const PointType &) const;

  MatrixType
  EvaluateLagrangianStrainTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateLagrangianStrainTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateLagrangianStrainTensor(point);
  }

  /**
   * Evaluate Lagrangian directional strain
   */
  RealType
  EvaluateLagrangianDirectionalStrain(const PointType &, const VectorType &) const;

  RealType
  EvaluateLagrangianDirectionalStrainAtIndex(const IndexType & idx, const VectorType & V) const;

  RealType
  EvaluateLagrangianDirectionalStrainAtContinuousIndex(const ContinuousIndexType & idx, const VectorType & V) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateLagrangianDirectionalStrain(point, V);
  }

  /**
   * Evaluate Eulerian strain tensor
   */
  MatrixType
  EvaluateEulerianStrainTensor(const PointType &) const;

  MatrixType
  EvaluateEulerianStrainTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateEulerianStrainTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateEulerianStrainTensor(point);
  }

  /**
   * Evaluate Eulerian directional strain
   */
  RealType
  EvaluateEulerianDirectionalStrain(const PointType &, const VectorType &) const;

  RealType
  EvaluateEulerianDirectionalStrainAtIndex(const IndexType & idx, const VectorType & V) const;

  RealType
  EvaluateEulerianDirectionalStrainAtContinuousIndex(const ContinuousIndexType & idx, const VectorType & V) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateEulerianDirectionalStrain(point, V);
  }

  /**
   * Evaluate Right Cauchy-Green strain tensor
   */
  MatrixType
  EvaluateRightCauchyGreenDeformationTensor(const PointType &) const;

  MatrixType
  EvaluateRightCauchyGreenDeformationTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateRightCauchyGreenDeformationTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateRightCauchyGreenDeformationTensor(point);
  }

  /**
   * Evaluate Left Cauchy-Green strain tensor
   */
  MatrixType
  EvaluateLeftCauchyGreenDeformationTensor(const PointType &) const;

  MatrixType
  EvaluateLeftCauchyGreenDeformationTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateLeftCauchyGreenDeformationTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateLeftCauchyGreenDeformationTensor(point);
  }

  /**
   * Evaluate left stretch tensor
   */
  MatrixType
  EvaluateLeftStretchTensor(const PointType &) const;

  MatrixType
  EvaluateLeftStretchTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateLeftStretchTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateLeftStretchTensor(point);
  }

  /**
   * Evaluate right stretch tensor
   */
  MatrixType
  EvaluateRightStretchTensor(const PointType &) const;

  MatrixType
  EvaluateRightStretchTensorAtIndex(const IndexType & idx) const;

  MatrixType
  EvaluateRightStretchTensorAtContinuousIndex(const ContinuousIndexType & idx) const
  {
    PointType point;

    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(idx, point);
    return this->EvaluateRightStretchTensor(point);
  }

protected:
  VectorFieldGradientImageFunction();
  virtual ~VectorFieldGradientImageFunction() override {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  VectorFieldGradientImageFunction(const Self &) = delete;
  MatrixType
  operator=(const Self &) = delete;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkVectorFieldGradientImageFunction.hxx"
#endif

#endif
