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
#ifndef itkMeanSquaresPointSetToPointSetIntensityMetricv4_h
#define itkMeanSquaresPointSetToPointSetIntensityMetricv4_h

#include "itkPointSetToPointSetMetricv4.h"

namespace itk
{
/** \class MeanSquaresPointSetToPointSetIntensityMetricv4
 * \brief Computes a mean-squares intensity metric between two point sets.
 *
 *  We only have to handle the individual point case as the parent
 *  class handles the aggregation.
 *
 *  The PointSet::PixelType contains both the intensity and gradient
 *  information at each point and the surrounding neighborhood.  The ordering
 *  is based on the itk::Neighborhood data structure:
 *
 *  http://www.itk.org/Doxygen/html/classitk_1_1Neighborhood.html
 *
 *  pixelData = [intensity_0, gradientX_0, gradientY_0, gradientZ_0,       // GetOffset( 0 )
 *               intensity_1, gradientX_1, gradientY_1, gradientZ_1,       // GetOffset( 1 )
 *               ...        , ...        , ...        , ...        ,         ...
 *               intensity_N, gradientX_N, gradientY_N, gradientZ_N]       // GetOffset( N )
 *
 *
 * \ingroup ITKMetricsv4
 */
template <typename TFixedPointSet,
          typename TMovingPointSet = TFixedPointSet,
          typename TInternalComputationValueType = double>
class MeanSquaresPointSetToPointSetIntensityMetricv4 final
  : public PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  /** Standard class typedefs. */
  typedef MeanSquaresPointSetToPointSetIntensityMetricv4                                             Self;
  typedef PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType> Superclass;
  typedef SmartPointer<Self>                                                                         Pointer;
  typedef SmartPointer<const Self>                                                                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(MeanSquaresPointSetToPointSetIntensityMetricv4);

  /**  Type of the fixed point set. */
  typedef TFixedPointSet                              FixedPointSetType;
  typedef typename TFixedPointSet::PointType          FixedPointType;
  typedef typename TFixedPointSet::PixelType          FixedPixelType;
  typedef typename TFixedPointSet::PointsContainer    FixedPointsContainer;
  typedef typename TFixedPointSet::PointDataContainer FixedPointDataContainer;

  /** Dimension type */
  typedef typename Superclass::DimensionType DimensionType;

  static constexpr DimensionType FixedPointDimension = Superclass::FixedDimension;

  /**  Type of the moving point set. */
  typedef TMovingPointSet                              MovingPointSetType;
  typedef typename TMovingPointSet::PointType          MovingPointType;
  typedef typename TMovingPointSet::PixelType          MovingPixelType;
  typedef typename TMovingPointSet::PointsContainer    MovingPointsContainer;
  typedef typename TMovingPointSet::PointDataContainer MovingPointDataContainer;

  static constexpr DimensionType MovingPointDimension = Superclass::MovingDimension;

  /** Transform types from Superclass*/
  typedef typename Superclass::FixedTransformType  FixedTransformType;
  typedef typename Superclass::MovingTransformType MovingTransformType;

  static constexpr DimensionType PointDimension = Superclass::PointDimension;

  /** Types transferred from the base class */
  typedef typename Superclass::MeasureType         MeasureType;
  typedef typename Superclass::DerivativeType      DerivativeType;
  typedef typename Superclass::LocalDerivativeType LocalDerivativeType;
  typedef typename Superclass::PointType           PointType;
  typedef typename Superclass::PixelType           PixelType;
  typedef typename Superclass::PointIdentifier     PointIdentifier;
  typedef typename Superclass::PointsConstIterator PointsConstIterator;

  typedef CovariantVector<TInternalComputationValueType, PointDimension> CovariantVectorType;

  /**
   * Set/get estimate intensity distance sigma automatically based on reasonable
   * heuristics.
   */
  itkSetMacro(EstimateIntensityDistanceSigmaAutomatically, bool);
  itkGetConstMacro(EstimateIntensityDistanceSigmaAutomatically, bool);
  itkBooleanMacro(EstimateIntensityDistanceSigmaAutomatically);

  /**
   * Set/get estimate Euclidean distance sigma automatically based on reasonable
   * heuristics.
   */
  itkSetMacro(EstimateEuclideanDistanceSigmaAutomatically, bool);
  itkGetConstMacro(EstimateEuclideanDistanceSigmaAutomatically, bool);
  itkBooleanMacro(EstimateEuclideanDistanceSigmaAutomatically);

  /**
   * Set/get intensity sigma -- modulate the intensity distance contribution in the
   * metric formulation.  Default = sqrt( 5.0 ).
   */
  itkSetMacro(IntensityDistanceSigma, TInternalComputationValueType);
  itkGetConstMacro(IntensityDistanceSigma, TInternalComputationValueType);

  /**
   * Set/get distance sigma -- modulate the Euclidean distance contribution in the
   * metric formulation.   Default = sqrt( 5.0 ).
   */
  itkSetMacro(EuclideanDistanceSigma, TInternalComputationValueType);
  itkGetConstMacro(EuclideanDistanceSigma, TInternalComputationValueType);

  /**
   * Initialize the metric by estimating the intensity and distance sigmas
   */
  void
  Initialize(void) override;

  /**
   * Prepare point sets for use.
   */
  void
  InitializePointSets() const override;

  /**
   * Calculates the local metric value for a single point.
   */
  MeasureType
  GetLocalNeighborhoodValue(const PointType &, const PixelType &) const override;

  /** Helper method allows for code reuse while skipping the metric value
   * calculation when appropriate */
  void
  CalculateValueAndDerivative(MeasureType & value, DerivativeType & derivative, bool calculateValue) const;

  /**
   * Calculates the local value and derivative for a single point.
   */
  void
  GetLocalNeighborhoodValueAndDerivative(const PointType &,
                                         MeasureType &,
                                         LocalDerivativeType &,
                                         const PixelType &) const override;

  /** Clone method will clone the existing instance of this type,
   *  including its internal member variables. */
  typename LightObject::Pointer
  InternalClone() const override;

protected:
  MeanSquaresPointSetToPointSetIntensityMetricv4();
  ~MeanSquaresPointSetToPointSetIntensityMetricv4() override;

  /**
   * Estimate the intensity distance sigma based on simple heuristic
   */
  void
  EstimateIntensityDistanceSigma();

  /**
   * Estimate the Euclidean distance sigma based on simple heuristic
   */
  void
  EstimateEuclideanDistanceSigma();

  /**
   * Warp the fixed point set gradients based on the fixed transform.
   */
  void
  TransformFixedPointSetGradients() const;

  /**
   * Warp the moving point set gradients based on the moving transform.
   */
  void
  TransformMovingPointSetGradients() const;


  /** PrintSelf function */
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  MeanSquaresPointSetToPointSetIntensityMetricv4(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  bool                          m_EstimateIntensityDistanceSigmaAutomatically;
  bool                          m_EstimateEuclideanDistanceSigmaAutomatically;
  TInternalComputationValueType m_IntensityDistanceSigma;
  TInternalComputationValueType m_EuclideanDistanceSigma;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMeanSquaresPointSetToPointSetIntensityMetricv4.hxx"
#endif

#endif
