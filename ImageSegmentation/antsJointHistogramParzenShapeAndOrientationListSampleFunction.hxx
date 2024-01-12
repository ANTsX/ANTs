/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsJointHistogramParzenShapeAndOrientationListSampleFunction_hxx
#define __antsJointHistogramParzenShapeAndOrientationListSampleFunction_hxx


#include "itkArray.h"
#include "itkContinuousIndex.h"
#include "itkDecomposeTensorFunction.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkStatisticsImageFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TListSample, typename TOutput, typename TCoordRep>
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::
  JointHistogramParzenShapeAndOrientationListSampleFunction()
{
  this->m_NumberOfShapeJointHistogramBins = 32;
  this->m_NumberOfOrientationJointHistogramBins = 64;
  this->m_ShapeSigma = 1.0;
  this->m_OrientationSigma = 2.0;
  this->m_UseNearestNeighborIncrements = true;
  this->m_MaximumEigenvalue1 = 0;
  this->m_MaximumEigenvalue2 = 0;
  this->m_MinimumEigenvalue1 = 1;
  this->m_MinimumEigenvalue2 = 1;
  this->m_Interpolator = InterpolatorType::New();
  this->m_Interpolator->SetSplineOrder(3);
  this->m_JointHistogramImages[0] = nullptr;
  this->m_JointHistogramImages[1] = nullptr;
  this->m_JointHistogramImages[2] = nullptr;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::
  ~JointHistogramParzenShapeAndOrientationListSampleFunction() = default;

template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::
  IncrementJointHistogramForShape(RealType eigenvalue1, RealType eigenvalue2)
{
  RealType newWeight = 1.0;

  // now define two joint histograms, one for shape, one for orientation.
  // first, the shape histogram --- 0,0 origin and spacing of 1
  if (!this->m_JointHistogramImages[0])
  {
    typename JointHistogramImageType::SpacingType spacing;
    spacing.Fill(1);
    typename JointHistogramImageType::PointType origin;
    origin.Fill(0);
    typename JointHistogramImageType::SizeType size;
    size.Fill(this->m_NumberOfShapeJointHistogramBins);
    this->m_JointHistogramImages[0] = AllocImage<JointHistogramImageType>(size);
    this->m_JointHistogramImages[0]->SetOrigin(origin);
    this->m_JointHistogramImages[0]->SetSpacing(spacing);
    this->m_JointHistogramImages[0]->FillBuffer(0);
  }

  typename JointHistogramImageType::PointType shapePoint;
  if (eigenvalue1 > NumericTraits<RealType>::OneValue())
  {
    eigenvalue1 = NumericTraits<RealType>::OneValue();
  }
  if (eigenvalue2 > NumericTraits<RealType>::OneValue())
  {
    eigenvalue2 = NumericTraits<RealType>::OneValue();
  }
  if (eigenvalue1 < NumericTraits<RealType>::ZeroValue())
  {
    eigenvalue1 = NumericTraits<RealType>::ZeroValue();
  }
  if (eigenvalue2 < NumericTraits<RealType>::ZeroValue())
  {
    eigenvalue2 = NumericTraits<RealType>::ZeroValue();
  }

  shapePoint[0] = eigenvalue1 * (this->m_NumberOfShapeJointHistogramBins - 1);
  shapePoint[1] = eigenvalue2 * (this->m_NumberOfShapeJointHistogramBins - 1);

  ContinuousIndex<double, 2> shapeCidx;
  shapeCidx = this->m_JointHistogramImages[0]->
          template TransformPhysicalPointToContinuousIndex<double, SpacePrecisionType>(shapePoint);

  typedef typename JointHistogramImageType::IndexType JointHistogramImageIndexType;
  JointHistogramImageIndexType                        shapeIdx;

  /** Nearest neighbor increment to JH */
  if (this->m_UseNearestNeighborIncrements)
  {
    shapeIdx[0] = static_cast<typename JointHistogramImageIndexType::IndexValueType>(std::floor(shapeCidx[0] + 0.5));
    shapeIdx[1] = static_cast<typename JointHistogramImageIndexType::IndexValueType>(std::floor(shapeCidx[1] + 0.5));
    if (this->m_JointHistogramImages[0]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[0]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[0]->SetPixel(shapeIdx, 1 + oldWeight);
    }
  }
  else
  {
    /** linear addition */
    shapeIdx[0] = static_cast<IndexValueType>(std::floor(shapeCidx[0]));
    shapeIdx[1] = static_cast<IndexValueType>(std::floor(shapeCidx[1]));
    RealType distance1 = std::sqrt(Math::sqr(shapeCidx[0] - shapeIdx[0]) + Math::sqr(shapeCidx[1] - shapeIdx[1]));
    shapeIdx[0]++;
    RealType distance2 = std::sqrt(Math::sqr(shapeCidx[0] - shapeIdx[0]) + Math::sqr(shapeCidx[1] - shapeIdx[1]));
    shapeIdx[1]++;
    RealType distance3 = std::sqrt(Math::sqr(shapeCidx[0] - shapeIdx[0]) + Math::sqr(shapeCidx[1] - shapeIdx[1]));
    shapeIdx[0]--;
    RealType distance4 = std::sqrt(Math::sqr(shapeCidx[0] - shapeIdx[0]) + Math::sqr(shapeCidx[1] - shapeIdx[1]));
    RealType sumDistance = distance1 + distance2 + distance3 + distance4;
    distance1 /= sumDistance;
    distance2 /= sumDistance;
    distance3 /= sumDistance;
    distance4 /= sumDistance;

    unsigned int whichHistogram = 0;
    shapeIdx[0] = static_cast<IndexValueType>(std::floor(shapeCidx[0]));
    shapeIdx[1] = static_cast<IndexValueType>(std::floor(shapeCidx[1]));
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        shapeIdx, (NumericTraits<RealType>::OneValue() - distance1) * newWeight + oldWeight);
    }
    shapeIdx[0]++;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        shapeIdx, (NumericTraits<RealType>::OneValue() - distance2) * newWeight + oldWeight);
    }
    shapeIdx[1]++;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        shapeIdx, (NumericTraits<RealType>::OneValue() - distance3) * newWeight + oldWeight);
    }
    shapeIdx[0]--;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        shapeIdx, (NumericTraits<RealType>::OneValue() - distance4) * newWeight + oldWeight);
    }
  }
  return;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::
  IncrementJointHistogramForOrientation(RealType x, RealType y, RealType z, unsigned int whichHistogram)
{
  RealType newWeight = 1.0;

  // 2nd, the orientation histogram.  origin 0,0.  spacing of 1,1.
  //  need to be careful for wrap around in the 0 to 2*pi case.
  if (!this->m_JointHistogramImages[whichHistogram])
  {
    typename JointHistogramImageType::SpacingType spacing2;
    spacing2.Fill(1);
    typename JointHistogramImageType::PointType origin2;
    origin2.Fill(0);
    typename JointHistogramImageType::SizeType size2;
    size2.Fill(this->m_NumberOfOrientationJointHistogramBins);
    size2[0] = size2[0] + 2;
    this->m_JointHistogramImages[whichHistogram] = AllocImage<JointHistogramImageType>(size2);

    this->m_JointHistogramImages[whichHistogram]->SetOrigin(origin2);
    this->m_JointHistogramImages[whichHistogram]->SetSpacing(spacing2);
    this->m_JointHistogramImages[whichHistogram]->FillBuffer(0);
  }

  JointHistogramImagePointType orientPoint;
  RealType                     tp[2];
  tp[1] = 0.0;

  // If eigenvector has negative x, we reflect it about the origin.
  // We do this to avoid redundancy in representation of the eigenvectors,
  // because they are all redundant.

  if (x < NumericTraits<RealType>::ZeroValue())
  {
    x *= -1;
    y *= -1;
    z *= -1;
  }

  tp[0] = std::acos(z);

  // phi goes from 0.0 (+x axis) and goes to -pi/2 and pi/2.
  // theta goes from 0.0 (+z axis) and wraps at PI
  // if x and y are 0.0 or very close, return phi == 0
  // we do this to eliminate redundancy in the distribution of orientations.

  if (Math::abs(x) + Math::abs(y) < static_cast<RealType>(1e-9))
  {
    tp[1] = NumericTraits<RealType>::ZeroValue();
  }
  else
  {
    if (Math::FloatAlmostEqual(y, NumericTraits<RealType>::ZeroValue()))
    {
      if (x > NumericTraits<RealType>::ZeroValue())
      {
        tp[1] = NumericTraits<RealType>::ZeroValue();
      }
      else
      {
        tp[1] = Math::pi;
      }
    }
    else if (Math::FloatAlmostEqual(x, NumericTraits<RealType>::ZeroValue()))
    {
      // avoid div by zero
      if (y > NumericTraits<RealType>::ZeroValue())
      {
        tp[1] = Math::pi_over_2;
      }
      else
      {
        tp[1] = -Math::pi_over_2;
      }
    }
    else if (x > NumericTraits<RealType>::ZeroValue() && y > NumericTraits<RealType>::ZeroValue())
    { // first quadrant
      tp[1] = std::atan(y / x);
    }
    else if (x < NumericTraits<RealType>::ZeroValue() && y > NumericTraits<RealType>::ZeroValue())
    { // second quadrant
      tp[1] = static_cast<RealType>(Math::pi) + std::atan(y / x);
    }
    else if (x < NumericTraits<RealType>::ZeroValue() && y < NumericTraits<RealType>::ZeroValue())
    { // third quadrant
      tp[1] = static_cast<RealType>(Math::pi) + std::atan(y / x);
    }
    else
    { // fourth quadrant
      tp[1] = std::atan(y / x);
    }
  }
  RealType psi = tp[0];
  RealType theta = tp[1];

  // note, if a point maps to 0 or 2*pi then it should contribute to both bins -- pretty much only difference between
  // this function and matlab code is the next 15 or so lines, as far as we see
  orientPoint[0] = static_cast<typename JointHistogramImagePointType::CoordRepType>(
    psi / static_cast<RealType>(Math::pi) * static_cast<RealType>(this->m_NumberOfOrientationJointHistogramBins - 1) +
    NumericTraits<RealType>::OneValue());
  orientPoint[1] = static_cast<typename JointHistogramImagePointType::CoordRepType>(
    (theta + static_cast<RealType>(Math::pi_over_2)) / static_cast<RealType>(Math::pi) *
    static_cast<RealType>(this->m_NumberOfOrientationJointHistogramBins - 1));

  ContinuousIndex<double, 2> orientCidx;
  orientCidx = this->m_JointHistogramImages[whichHistogram]->
          template TransformPhysicalPointToContinuousIndex<double, itk::SpacePrecisionType>(orientPoint);

  typedef typename JointHistogramImageType::IndexType JointHistogramImageIndexType;
  JointHistogramImageIndexType                        orientIdx;

  /** Nearest neighbor interpolation */
  if (this->m_UseNearestNeighborIncrements)
  {
    orientIdx[0] = static_cast<typename JointHistogramImageIndexType::IndexValueType>(std::floor(orientCidx[0] + 0.5));
    orientIdx[1] = static_cast<typename JointHistogramImageIndexType::IndexValueType>(std::floor(orientCidx[1] + 0.5));
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(orientIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(orientIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(orientIdx, 1 + oldWeight);
    }
  }
  else
  {
    orientIdx[0] = static_cast<IndexValueType>(std::floor(orientCidx[0]));
    orientIdx[1] = static_cast<IndexValueType>(std::floor(orientCidx[1]));
    RealType distance1 = std::sqrt(Math::sqr(orientCidx[0] - orientIdx[0]) + Math::sqr(orientCidx[1] - orientIdx[1]));
    orientIdx[0]++;
    RealType distance2 = std::sqrt(Math::sqr(orientCidx[0] - orientIdx[0]) + Math::sqr(orientCidx[1] - orientIdx[1]));
    orientIdx[1]++;
    RealType distance3 = std::sqrt(Math::sqr(orientCidx[0] - orientIdx[0]) + Math::sqr(orientCidx[1] - orientIdx[1]));
    orientIdx[0]--;
    RealType distance4 = std::sqrt(Math::sqr(orientCidx[0] - orientIdx[0]) + Math::sqr(orientCidx[1] - orientIdx[1]));
    RealType sumDistance = distance1 + distance2 + distance3 + distance4;
    distance1 /= sumDistance;
    distance2 /= sumDistance;
    distance3 /= sumDistance;
    distance4 /= sumDistance;

    orientIdx[0] = static_cast<IndexValueType>(std::floor(orientCidx[0]));
    orientIdx[1] = static_cast<IndexValueType>(std::floor(orientCidx[1]));
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(orientIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(orientIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        orientIdx, (NumericTraits<RealType>::OneValue() - distance1) * newWeight + oldWeight);
    }
    orientIdx[0]++;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(orientIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(orientIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        orientIdx, (NumericTraits<RealType>::OneValue() - distance2) * newWeight + oldWeight);
    }
    orientIdx[1]++;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(orientIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(orientIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        orientIdx, (NumericTraits<RealType>::OneValue() - distance3) * newWeight + oldWeight);
    }
    orientIdx[0]--;
    if (this->m_JointHistogramImages[whichHistogram]->GetLargestPossibleRegion().IsInside(orientIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[whichHistogram]->GetPixel(orientIdx);
      this->m_JointHistogramImages[whichHistogram]->SetPixel(
        orientIdx, (NumericTraits<RealType>::OneValue() - distance4) * newWeight + oldWeight);
    }
  }

  // The last thing we do is copy the [1,] column to the [NBins+1,] column and
  // the [NBins,] column to the [0,] column --- circular boundary conditions.

  typedef itk::ImageRegionIteratorWithIndex<JointHistogramImageType> Iterator;
  Iterator                                                           tIter(this->m_JointHistogramImages[whichHistogram],
                 this->m_JointHistogramImages[whichHistogram]->GetBufferedRegion());
  for (tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter)
  {
    IndexType index = tIter.GetIndex();
    IndexType index2 = tIter.GetIndex();
    if (index[0] == 0)
    {
      index2[0] = this->m_NumberOfOrientationJointHistogramBins;
      index2[1] = index[1];
      tIter.Set(this->m_JointHistogramImages[whichHistogram]->GetPixel(index2));
    }
    if (index[0] == static_cast<typename IndexType::IndexValueType>(this->m_NumberOfOrientationJointHistogramBins + 1))
    {
      index2[0] = 1;
      index2[1] = index[1];
      tIter.Set(this->m_JointHistogramImages[whichHistogram]->GetPixel(index2));
    }
  }
  return;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::SetInputListSample(
  const InputListSampleType * ptr)
{
  Superclass::SetInputListSample(ptr);

  if (!this->GetInputListSample())
  {
    return;
  }

  if (this->GetInputListSample()->Size() <= 1)
  {
    itkWarningMacro("The input list sample has <= 1 element."
                    << "Function evaluations will be equal to 0.");
    return;
  }

  const unsigned int Dimension = this->GetInputListSample()->GetMeasurementVectorSize();

  /**
   * Find the min/max values to define the histogram domain
   */
  Array<RealType> minValues(Dimension);
  minValues.Fill(NumericTraits<RealType>::max());
  Array<RealType> maxValues(Dimension);
  maxValues.Fill(NumericTraits<RealType>::NonpositiveMin());

  typename InputListSampleType::ConstIterator It = this->GetInputListSample()->Begin();
  while (It != this->GetInputListSample()->End())
  {
    InputMeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      if (inputMeasurement[d] < minValues[d])
      {
        minValues[d] = inputMeasurement[d];
      }
      if (inputMeasurement[d] > maxValues[d])
      {
        maxValues[d] = inputMeasurement[d];
      }
    }
    ++It;
  }
  for (unsigned int d = 0; d < 3; d++)
  {
    this->m_JointHistogramImages[d] = nullptr;
  }
  RealType     L = static_cast<RealType>(this->GetInputListSample()->GetMeasurementVectorSize());
  unsigned int D = static_cast<unsigned int>(
    static_cast<RealType>(0.5) * (-NumericTraits<RealType>::OneValue() +
                                  std::sqrt(NumericTraits<RealType>::OneValue() + static_cast<RealType>(8.0) * L)));

  It = this->GetInputListSample()->Begin();
  while (It != this->GetInputListSample()->End())
  {
    InputMeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    // convert to a tensor then get its shape and primary orientation vector
    typedef VariableSizeMatrix<RealType> TensorType;
    TensorType                           T(D, D);
    T.Fill(0.0);
    unsigned int index = 0;
    for (unsigned int i = 0; i < D; i++)
    {
      for (unsigned int j = i; j < D; j++)
      {
        T(i, j) = inputMeasurement(index++);
        T(j, i) = T(i, j);
      }
    }
    // now decompose T into shape and orientation
    TensorType                                  V;
    TensorType                                  W;
    TensorType                                  Tc = T;
    typedef DecomposeTensorFunction<TensorType> DecomposerType;
    typename DecomposerType::Pointer            decomposer = DecomposerType::New();
    decomposer->EvaluateSymmetricEigenDecomposition(Tc, W, V);
    // now W holds the eigenvalues ( shape )

    // for each tensor sample, we add its content to the relevant histogram.
    RealType eigenvalue1 = W(2, 2);
    RealType eigenvalue2 = (W(1, 1)); // + W(0, 0) ) * 0.5;
    if (eigenvalue1 > this->m_MaximumEigenvalue1)
    {
      this->m_MaximumEigenvalue1 = eigenvalue1;
    }
    if (eigenvalue2 > this->m_MaximumEigenvalue2)
    {
      this->m_MaximumEigenvalue2 = eigenvalue2;
    }
    if (eigenvalue1 < this->m_MinimumEigenvalue1)
    {
      this->m_MinimumEigenvalue1 = eigenvalue1;
    }
    if (eigenvalue2 < this->m_MinimumEigenvalue2)
    {
      this->m_MinimumEigenvalue2 = eigenvalue2;
    }
    ++It;
  }

  It = this->GetInputListSample()->Begin();
  while (It != this->GetInputListSample()->End())
  {
    InputMeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    // convert to a tensor then get its shape and primary orientation vector
    typedef VariableSizeMatrix<RealType> TensorType;
    TensorType                           T(D, D);
    T.Fill(0.0);
    unsigned int index = 0;
    for (unsigned int i = 0; i < D; i++)
    {
      for (unsigned int j = i; j < D; j++)
      {
        T(i, j) = inputMeasurement(index++);
        T(j, i) = T(i, j);
      }
    }
    // now decompose T into shape and orientation
    TensorType                                  V;
    TensorType                                  W;
    TensorType                                  Tc = T;
    typedef DecomposeTensorFunction<TensorType> DecomposerType;
    typename DecomposerType::Pointer            decomposer = DecomposerType::New();
    decomposer->EvaluateSymmetricEigenDecomposition(Tc, W, V);
    // now W holds the eigenvalues ( shape )

    // for each tensor sample, we add its content to the relevant histogram.
    RealType eigenvalue1 = W(2, 2);
    RealType eigenvalue2 = W(1, 1);
    eigenvalue1 /= (this->m_MaximumEigenvalue1 - this->m_MinimumEigenvalue1);
    eigenvalue2 /= (this->m_MaximumEigenvalue2 - this->m_MinimumEigenvalue2);

    //    std::cout << " ev1 " << eigenvalue1 << " oev1 " << W(2,2) << " ev2 " << eigenvalue2 << " oev2 " <<
    // W(1,1) <<
    // std::endl;

    /** joint-hist model for the eigenvalues */
    this->IncrementJointHistogramForShape(eigenvalue1, eigenvalue2);

    RealType x = V(0, 2);
    RealType y = V(1, 2);
    RealType z = V(2, 2);

    /** joint-hist model for the principal eigenvector */
    this->IncrementJointHistogramForOrientation(x, y, z, 1);
    x = V(0, 1);
    y = V(1, 1);
    z = V(2, 1);
    /** joint-hist model for the second eigenvector */
    this->IncrementJointHistogramForOrientation(x, y, z, 2);

    ++It;
  }
  for (unsigned int d = 0; d < 3; d++)
  {
    typedef DiscreteGaussianImageFilter<JointHistogramImageType, JointHistogramImageType> GaussianFilterType;
    typename GaussianFilterType::Pointer gaussian = GaussianFilterType::New();
    gaussian->SetInput(this->m_JointHistogramImages[d]);
    if (d == 0) // Shape
    {
      gaussian->SetVariance(this->m_ShapeSigma * this->m_ShapeSigma);
    }
    else if (d == 1) // Orientation of 1st eigenvector
    {
      gaussian->SetVariance(this->m_OrientationSigma * this->m_OrientationSigma);
    }
    else if (d == 2) // Orientation of 2nd eigenvector
    {
      gaussian->SetVariance(this->m_ShapeSigma * this->m_OrientationSigma);
    }
    gaussian->SetMaximumError(0.01);
    gaussian->SetUseImageSpacing(false);
    gaussian->Update();

    typedef StatisticsImageFilter<JointHistogramImageType> StatsFilterType;
    typename StatsFilterType::Pointer                      stats = StatsFilterType::New();
    stats->SetInput(gaussian->GetOutput());
    stats->Update();

    typedef DivideImageFilter<JointHistogramImageType, JointHistogramImageType, JointHistogramImageType> DividerType;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput(gaussian->GetOutput());
    divider->SetConstant(stats->GetSum());
    divider->Update();
    this->m_JointHistogramImages[d] = divider->GetOutput();
  }
  /*  write out histograms--for debugging
      static int which_class=0;
      which_class++;
      std::string string;
      std::stringstream outstring;
      outstring<<which_class;
      string=outstring.str();
      typedef ImageFileWriter< JointHistogramImageType >  WriterType;
      typename WriterType::Pointer      writer = WriterType::New();
      std::string output( "output_shape"+string+".nii.gz" );
      writer->SetFileName( output.c_str() );
      writer->SetInput(this->m_JointHistogramImages[0] );
      writer->Update();
      typedef ImageFileWriter< JointHistogramImageType >  WriterType2;
      typename WriterType2::Pointer      writer2 = WriterType::New();
      std::string output2( "output_orientation"+string+".nii.gz" );
      writer2->SetFileName( output2.c_str() );
      writer2->SetInput(this->m_JointHistogramImages[1] );
      writer2->Update();
      std::cout << "Writing output of histograms." << std::endl;
  */
}

template <typename TListSample, typename TOutput, typename TCoordRep>
TOutput
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::Evaluate(
  const InputMeasurementVectorType & measurement) const
{
  try
  {
    RealType probability = 1.0;
    for (unsigned int d = 0; d < 2; d++)
    {
      typename JointHistogramImageType::PointType point;
      point[0] = measurement[d];

      this->m_Interpolator->SetInputImage(this->m_JointHistogramImages[d]);
      if (this->m_Interpolator->IsInsideBuffer(point))
      {
        probability *= static_cast<RealType>(this->m_Interpolator->Evaluate(point));
      }
      else
      {
        return 0;
      }
    }
    return probability;
  }
  catch (...)
  {
    return 0;
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenShapeAndOrientationListSampleFunction<TListSample, TOutput, TCoordRep>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  os << indent << "Shape Sigma: " << this->m_ShapeSigma << std::endl;
  os << indent << "Number of shape histogram bins: " << this->m_NumberOfShapeJointHistogramBins << std::endl;
  os << indent << "Orientation Sigma: " << this->m_OrientationSigma << std::endl;
  os << indent << "Number of orientation histogram bins: " << this->m_NumberOfOrientationJointHistogramBins
     << std::endl;
  os << indent << "Minimum eigenvalue 1: " << this->m_MinimumEigenvalue1;
  os << indent << "Minimum eigenvalue 2: " << this->m_MinimumEigenvalue2;
  os << indent << "Maximum eigenvalue 1: " << this->m_MaximumEigenvalue1;
  os << indent << "Maximum eigenvalue 2: " << this->m_MaximumEigenvalue2;

  if (this->m_UseNearestNeighborIncrements)
  {
    os << indent << "Use nearest neighbor increments." << std::endl;
  }
  else
  {
    os << indent << "Use linear interpolation for increments." << std::endl;
  }
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
