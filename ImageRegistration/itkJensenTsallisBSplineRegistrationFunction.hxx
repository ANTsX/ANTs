/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkJensenTsallisBSplineRegistrationFunction_hxx_
#define _itkJensenTsallisBSplineRegistrationFunction_hxx_


#include "itkContinuousIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"
// #include "itkBSplineControlPointImageFilter.h"

namespace itk
{
template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet,
                                         TMovingImage,
                                         TMovingPointSet,
                                         TDisplacementField>::JensenTsallisBSplineRegistrationFunction()
{
  this->m_UseRegularizationTerm = false;
  this->m_UseInputAsSamples = true;
  this->m_UseAnisotropicCovariances = false;

  this->m_NumberOfFixedSamples = 100;
  this->m_FixedPointSetSigma = 1.0;
  this->m_FixedKernelSigma = 0.0;
  this->m_FixedEvaluationKNeighborhood = 50;

  this->m_NumberOfMovingSamples = 100;
  this->m_MovingPointSetSigma = 1.0;
  this->m_MovingKernelSigma = 0.0;
  this->m_MovingEvaluationKNeighborhood = 50;

  unsigned int covarianceKNeighborhood =
    static_cast<unsigned int>(std::pow(3.0, static_cast<RealType>(ImageDimension))) - 1;

  this->m_FixedCovarianceKNeighborhood = covarianceKNeighborhood;
  this->m_MovingCovarianceKNeighborhood = covarianceKNeighborhood;

  this->m_Alpha = 2.0;

  //  this->m_FixedControlPointLattice = nullptr;
  //  this->m_MovingControlPointLattice = nullptr;

  this->m_DerivativeFixedField = nullptr;
  this->m_DerivativeMovingField = nullptr;
  this->m_IsPointSetMetric = true;

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels = 1;
  this->m_MeshResolution.Fill(1);
}

template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
void
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet,
                                         TMovingImage,
                                         TMovingPointSet,
                                         TDisplacementField>::InitializeIteration(void)
{
  if (this->m_FixedKernelSigma == 0)
  {
    double maxFixedSpacing = 0.0;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      if (this->GetFixedImage()->GetSpacing()[d])
      {
        maxFixedSpacing = this->GetFixedImage()->GetSpacing()[d];
      }
    }
    this->m_FixedKernelSigma = 2.0 * maxFixedSpacing;
  }

  if (this->m_MovingKernelSigma == 0)
  {
    double maxMovingSpacing = 0.0;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      if (this->GetMovingImage()->GetSpacing()[d])
      {
        maxMovingSpacing = this->GetMovingImage()->GetSpacing()[d];
      }
    }
    this->m_MovingKernelSigma = 2.0 * maxMovingSpacing;
  }

  typename PointSetMetricType::Pointer pointSetMetric = PointSetMetricType::New();
  pointSetMetric->SetFixedPointSet(this->m_FixedPointSet);
  pointSetMetric->SetMovingPointSet(this->m_MovingPointSet);

  pointSetMetric->SetUseRegularizationTerm(this->m_UseRegularizationTerm);
  pointSetMetric->SetUseInputAsSamples(this->m_UseInputAsSamples);
  pointSetMetric->SetUseAnisotropicCovariances(this->m_UseAnisotropicCovariances);
  pointSetMetric->SetNumberOfFixedSamples(this->m_NumberOfFixedSamples);
  pointSetMetric->SetFixedPointSetSigma(this->m_FixedPointSetSigma);
  pointSetMetric->SetFixedKernelSigma(this->m_FixedKernelSigma);
  pointSetMetric->SetFixedCovarianceKNeighborhood(this->m_FixedCovarianceKNeighborhood);
  pointSetMetric->SetFixedEvaluationKNeighborhood(this->m_FixedEvaluationKNeighborhood);
  pointSetMetric->SetNumberOfMovingSamples(this->m_NumberOfMovingSamples);
  pointSetMetric->SetMovingPointSetSigma(this->m_MovingPointSetSigma);
  pointSetMetric->SetMovingKernelSigma(this->m_MovingKernelSigma);
  pointSetMetric->SetMovingCovarianceKNeighborhood(this->m_MovingCovarianceKNeighborhood);
  pointSetMetric->SetMovingEvaluationKNeighborhood(this->m_MovingEvaluationKNeighborhood);
  pointSetMetric->SetAlpha(this->m_Alpha);

  pointSetMetric->Initialize();

  typename PointSetMetricType::DefaultTransformType::ParametersType parameters;
  parameters.Fill(0.0);

  /**
   * Calculate with respect to the moving point set
   */

  pointSetMetric->SetUseWithRespectToTheMovingPointSet(true);
  typename PointSetMetricType::DerivativeType movingGradient;
  typename PointSetMetricType::MeasureType    movingMeasure;

  pointSetMetric->GetValueAndDerivative(parameters, movingMeasure, movingGradient);
  this->m_Energy += movingMeasure[0];

  typename BSplinePointSetType::Pointer movingGradientPoints = BSplinePointSetType::New();
  movingGradientPoints->Initialize();

  typename BSplineWeightsType::Pointer movingWeights = BSplineWeightsType::New();
  movingWeights->Initialize();

  typename MovingImageType::SizeType  movingSize = this->GetMovingImage()->GetLargestPossibleRegion().GetSize();
  typename MovingImageType::IndexType movingIndex = this->GetMovingImage()->GetLargestPossibleRegion().GetIndex();

  unsigned long count = 0;
  for (unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++)
  {
    typename MovingPointSetType::PointType point;
    this->m_MovingPointSet->GetPoint(n, &point);

    typename MovingImageType::PointType imagePoint;
    imagePoint.CastFrom(point);

    ContinuousIndex<typename MovingImageType::PointValueType, PointDimension> cidx;

    bool isInside = this->GetMovingImage()->TransformPhysicalPointToContinuousIndex(imagePoint, cidx);

    typename BSplinePointSetType::PointType bsplinePoint;
    for (unsigned int d = 0; d < PointDimension; d++)
    {
      if (cidx[d] - movingIndex[d] > movingSize[d] - 1.0)
      {
        isInside = false;
        break;
      }
      bsplinePoint[d] = cidx[d];
    }

    VectorType gradient;

    if (isInside)
    {
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        gradient[d] = movingGradient(n, d);
      }
      movingGradientPoints->SetPoint(count, bsplinePoint);
      movingGradientPoints->SetPointData(count, gradient);
      movingWeights->InsertElement(count, 1.0);
      count++;
    }
  }

  VectorType zeroVector;
  zeroVector.Fill(0.0);

  ImageLinearConstIteratorWithIndex<MovingImageType> ItM(this->GetMovingImage(),
                                                         this->GetMovingImage()->GetLargestPossibleRegion());
  for (unsigned int d = 0; d < PointDimension; d++)
  {
    ItM.SetDirection(d);
    ItM.GoToBegin();
    while (!ItM.IsAtEnd())
    {
      typename MovingImageType::PointType     point;
      typename BSplinePointSetType::PointType bsplinePoint;

      ItM.GoToBeginOfLine();
      typename MovingImageType::IndexType index = ItM.GetIndex();
      for (unsigned int m = 0; m < PointDimension; m++)
      {
        bsplinePoint[m] = index[m];
      }

      movingGradientPoints->SetPoint(count, bsplinePoint);
      movingGradientPoints->SetPointData(count, zeroVector);
      count++;

      ItM.GoToEndOfLine();
      --ItM;
      index = ItM.GetIndex();
      for (unsigned int m = 0; m < PointDimension; m++)
      {
        bsplinePoint[m] = index[m];
      }

      movingGradientPoints->SetPoint(count, bsplinePoint);
      movingGradientPoints->SetPointData(count, zeroVector);
      movingWeights->InsertElement(count, 1000.0);
      count++;

      ItM.NextLine();
    }
  }

  ArrayType numberOfMovingControlPoints;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    numberOfMovingControlPoints[d] =
      this->m_SplineOrder +
      static_cast<unsigned int>(
        std::floor(0.5 + static_cast<RealType>(this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[d]) /
                           static_cast<RealType>(this->m_MeshResolution[d])));
  }

  typename MovingImageType::SpacingType movingParametricSpacing;
  movingParametricSpacing.Fill(1);
  typename MovingImageType::PointType movingParametricOrigin;
  for (unsigned int d = 0; d < PointDimension; d++)
  {
    movingParametricOrigin[d] = this->GetMovingImage()->GetLargestPossibleRegion().GetIndex()[d];
  }

  typename BSplineFilterType::Pointer movingBSpliner = BSplineFilterType::New();
  movingBSpliner->SetInput(movingGradientPoints);
  movingBSpliner->SetPointWeights(movingWeights.GetPointer());
  movingBSpliner->SetOrigin(movingParametricOrigin);
  movingBSpliner->SetSpacing(movingParametricSpacing);
  movingBSpliner->SetSize(this->GetMovingImage()->GetLargestPossibleRegion().GetSize());
  movingBSpliner->SetDirection(this->GetMovingImage()->GetDirection());
  movingBSpliner->SetNumberOfLevels(this->m_NumberOfLevels);
  movingBSpliner->SetSplineOrder(this->m_SplineOrder);
  movingBSpliner->SetNumberOfControlPoints(numberOfMovingControlPoints);
  movingBSpliner->SetGenerateOutputImage(true);
  movingBSpliner->Update();

  movingBSpliner->GetOutput()->SetSpacing(this->GetMovingImage()->GetSpacing());
  this->m_DerivativeMovingField = movingBSpliner->GetOutput();
  this->m_DerivativeMovingField->DisconnectPipeline();

  //  this->m_MovingControlPointLattice = movingBSpliner->GetPhiLattice();

  /**
   * Calculate with respect to the fixed point set
   */

  pointSetMetric->SetUseWithRespectToTheMovingPointSet(false);
  typename PointSetMetricType::DerivativeType fixedGradient;
  typename PointSetMetricType::MeasureType    fixedMeasure;
  pointSetMetric->GetValueAndDerivative(parameters, fixedMeasure, fixedGradient);
  this->m_Energy += fixedMeasure[0];

  typename BSplinePointSetType::Pointer fixedGradientPoints = BSplinePointSetType::New();
  fixedGradientPoints->Initialize();

  typename BSplineWeightsType::Pointer fixedWeights = BSplineWeightsType::New();
  fixedWeights->Initialize();

  typename FixedImageType::SizeType  fixedSize = this->GetFixedImage()->GetLargestPossibleRegion().GetSize();
  typename FixedImageType::IndexType fixedIndex = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex();

  count = 0;
  for (unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++)
  {
    typename FixedPointSetType::PointType point;
    this->m_FixedPointSet->GetPoint(n, &point);

    typename FixedImageType::PointType imagePoint;
    imagePoint.CastFrom(point);

    ContinuousIndex<typename FixedImageType::PointValueType, PointDimension> cidx;

    bool isInside = this->GetFixedImage()->TransformPhysicalPointToContinuousIndex(imagePoint, cidx);

    typename BSplinePointSetType::PointType bsplinePoint;
    for (unsigned int d = 0; d < PointDimension; d++)
    {
      if (cidx[d] - fixedIndex[d] > fixedSize[d] - 1.0)
      {
        isInside = false;
        break;
      }
      bsplinePoint[d] = cidx[d];
    }

    VectorType gradient;

    if (isInside)
    {
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        gradient[d] = fixedGradient(n, d);
      }

      fixedGradientPoints->SetPoint(count, bsplinePoint);
      fixedGradientPoints->SetPointData(count, gradient);
      fixedWeights->InsertElement(count, 1.0);
      count++;
    }
  }

  ImageLinearConstIteratorWithIndex<FixedImageType> ItF(this->GetFixedImage(),
                                                        this->GetFixedImage()->GetLargestPossibleRegion());
  for (unsigned int d = 0; d < PointDimension; d++)
  {
    ItF.SetDirection(d);
    ItF.GoToBegin();
    while (!ItF.IsAtEnd())
    {
      typename FixedImageType::PointType      point;
      typename BSplinePointSetType::PointType bsplinePoint;

      ItF.GoToBeginOfLine();
      typename FixedImageType::IndexType index = ItF.GetIndex();
      for (unsigned int m = 0; m < PointDimension; m++)
      {
        bsplinePoint[m] = index[m];
      }

      fixedGradientPoints->SetPoint(count, bsplinePoint);
      fixedGradientPoints->SetPointData(count, zeroVector);
      fixedWeights->InsertElement(count, 1.0);
      count++;

      ItF.GoToEndOfLine();
      --ItF;
      index = ItF.GetIndex();
      for (unsigned int m = 0; m < PointDimension; m++)
      {
        bsplinePoint[m] = index[m];
      }

      fixedGradientPoints->SetPoint(count, bsplinePoint);
      fixedGradientPoints->SetPointData(count, zeroVector);
      fixedWeights->InsertElement(count, 1000.0);
      count++;

      ItF.NextLine();
    }
  }

  ArrayType numberOfFixedControlPoints;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    numberOfFixedControlPoints[d] =
      this->m_SplineOrder +
      static_cast<unsigned int>(
        std::floor(0.5 + static_cast<RealType>(this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[d]) /
                           static_cast<RealType>(this->m_MeshResolution[d])));
  }

  typename FixedImageType::SpacingType fixedParametricSpacing;
  fixedParametricSpacing.Fill(1.0);
  typename FixedImageType::PointType fixedParametricOrigin;
  for (unsigned int d = 0; d < PointDimension; d++)
  {
    fixedParametricOrigin[d] = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex()[d];
  }

  typename BSplineFilterType::Pointer fixedBSpliner = BSplineFilterType::New();
  fixedBSpliner->SetInput(fixedGradientPoints);
  fixedBSpliner->SetPointWeights(fixedWeights.GetPointer());
  fixedBSpliner->SetOrigin(fixedParametricOrigin);
  fixedBSpliner->SetSpacing(fixedParametricSpacing);
  fixedBSpliner->SetDirection(this->GetFixedImage()->GetDirection());
  fixedBSpliner->SetSize(this->GetFixedImage()->GetLargestPossibleRegion().GetSize());
  fixedBSpliner->SetNumberOfLevels(this->m_NumberOfLevels);
  fixedBSpliner->SetSplineOrder(this->m_SplineOrder);
  fixedBSpliner->SetNumberOfControlPoints(numberOfFixedControlPoints);
  fixedBSpliner->SetGenerateOutputImage(true);
  fixedBSpliner->Update();

  fixedBSpliner->GetOutput()->SetSpacing(this->GetFixedImage()->GetSpacing());
  this->m_DerivativeFixedField = fixedBSpliner->GetOutput();
  this->m_DerivativeFixedField->DisconnectPipeline();

  //  this->m_FixedControlPointLattice = fixedBSpliner->GetPhiLattice();
}

template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
typename JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                                  TFixedPointSet,
                                                  TMovingImage,
                                                  TMovingPointSet,
                                                  TDisplacementField>::VectorType
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet,
                                         TMovingImage,
                                         TMovingPointSet,
                                         TDisplacementField>::ComputeUpdate(const NeighborhoodType & neighborhood,
                                                                            void *                   globalData,
                                                                            const FloatOffsetType &  offset)
{
  if (this->m_DerivativeFixedField)
  {
    return -this->m_DerivativeFixedField->GetPixel(neighborhood.GetIndex());
  }
  else
  {
    itkExceptionMacro("Initialize() has not been called.");
  }

  /*
    typedef BSplineControlPointImageFilter<ControlPointLatticeType,
      DisplacementFieldType> BSplineControlPointImageFilterType;

    typename BSplineControlPointImageFilterType::Pointer movingBSpliner
      = BSplineControlPointImageFilterType::New();
    movingBSpliner->SetInput( this->m_MovingControlPointLattice );
    movingBSpliner->SetOrigin( this->GetMovingImage()->GetOrigin() );
    movingBSpliner->SetSpacing( this->GetMovingImage()->GetSpacing() );
    movingBSpliner->SetSize(
      this->GetMovingImage()->GetLargestPossibleRegion().GetSize() );
    movingBSpliner->SetSplineOrder( this->m_SplineOrder );

    VectorType gradient;
    movingBSpliner->EvaluateAtIndex( neighborhood.GetIndex(),
      gradient );

    return gradient;
  */
}

template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
typename JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                                  TFixedPointSet,
                                                  TMovingImage,
                                                  TMovingPointSet,
                                                  TDisplacementField>::VectorType
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet,
                                         TMovingImage,
                                         TMovingPointSet,
                                         TDisplacementField>::ComputeUpdateInv(const NeighborhoodType & neighborhood,
                                                                               void *                   globalData,
                                                                               const FloatOffsetType &  offset)
{
  if (this->m_DerivativeMovingField)
  {
    return -this->m_DerivativeMovingField->GetPixel(neighborhood.GetIndex());
  }
  else
  {
    itkExceptionMacro("Initialize() has not been called.");
  }

  /*
    typedef BSplineControlPointImageFilter<ControlPointLatticeType,
      DisplacementFieldType> BSplineControlPointImageFilterType;

    typename BSplineControlPointImageFilterType::Pointer fixedBSpliner
      = BSplineControlPointImageFilterType::New();
    fixedBSpliner->SetInput( this->m_FixedControlPointLattice );
    fixedBSpliner->SetOrigin( this->GetFixedImage()->GetOrigin() );
    fixedBSpliner->SetSpacing( this->GetFixedImage()->GetSpacing() );
    fixedBSpliner->SetSize(
      this->GetFixedImage()->GetLargestPossibleRegion().GetSize() );
    fixedBSpliner->SetSplineOrder( this->m_SplineOrder );

    VectorType gradient;
    fixedBSpliner->EvaluateAtIndex( neighborhood.GetIndex(),
      gradient );

    return gradient;
  */
}

template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
void
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet,
                                         TMovingImage,
                                         TMovingPointSet,
                                         TDisplacementField>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Use regularization term: " << this->m_UseRegularizationTerm << std::endl;
  if (!this->m_UseInputAsSamples)
  {
    os << indent << "Number of fixed samples: " << this->m_NumberOfFixedSamples << std::endl;
    os << indent << "Number of moving samples: " << this->m_NumberOfMovingSamples << std::endl;
  }
  os << indent << "Alpha: " << this->m_Alpha << std::endl;
  os << indent << "Fixed sigma: " << this->m_FixedPointSetSigma << std::endl;
  os << indent << "Moving sigma: " << this->m_MovingPointSetSigma << std::endl;

  os << indent << "Spline order: " << this->m_SplineOrder << std::endl;
  os << indent << "Number of levels: " << this->m_NumberOfLevels << std::endl;
  os << indent << "Number of control points: " << this->m_MeshResolution << std::endl;
}
} // end namespace itk

#endif
