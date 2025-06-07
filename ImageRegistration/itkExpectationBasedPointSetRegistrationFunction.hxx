/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkExpectationBasedPointSetRegistrationFunction_hxx_
#define _itkExpectationBasedPointSetRegistrationFunction_hxx_

#include "itkMacro.h"
#include "itkMath.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{
/*
 * Default constructor
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  ExpectationBasedPointSetRegistrationFunction()
{
  RadiusType   r;
  unsigned int j;

  for (j = 0; j < ImageDimension; j++)
  {
    r[j] = 0;
  }
  this->SetRadius(r);

  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_EuclideanDistanceThreshold = 0.01;
  this->SetMovingImage(nullptr);
  this->SetFixedImage(nullptr);
  m_FixedImageSpacing.Fill(1.0);
  m_FixedImageOrigin.Fill(0.0);
  m_Normalizer = 1.0;
  m_FixedImageGradientCalculator = GradientCalculatorType::New();

  m_Metric = NumericTraits<double>::max();
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RMSChange = NumericTraits<double>::max();
  m_SumOfSquaredChange = 0.0;
  this->m_KNeighborhood = 100;

  m_MovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
  m_UseMovingImageGradient = false;

  this->m_FixedPointSet = nullptr;
  this->m_MovingPointSet = nullptr;
  this->m_DerivativeFixedField = nullptr;
  this->m_DerivativeMovingField = nullptr;
  this->m_IsPointSetMetric = true;
  this->m_UseSymmetricMatching = 100000;
  this->m_Iterations = 0;
}

/*
 * Standard "PrintSelf" method.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  SetEuclideanDistanceThreshold(double threshold)
{
  m_EuclideanDistanceThreshold = threshold;
}

/**
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
double
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  GetEuclideanDistanceThreshold() const
{
  return m_EuclideanDistanceThreshold;
}

/*
 * Set the function state values before each iteration
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  ExpectationLandmarkField(float weight, bool whichdirection)
{

  SpacingType   spacing = this->GetFixedImage()->GetSpacing();
  unsigned long sz1 = this->m_FixedPointSet->GetNumberOfPoints();
  unsigned long sz2 = this->m_MovingPointSet->GetNumberOfPoints();
  if (!whichdirection)
  {
    sz2 = this->m_FixedPointSet->GetNumberOfPoints();
    sz1 = this->m_MovingPointSet->GetNumberOfPoints();
  }

  using MatrixType = vnl_matrix<double>;
  MatrixType EucDist(sz1, sz2);
  EucDist.fill(0);

  MatrixType fixedlms(sz1, ImageDimension);
  MatrixType movinglms(sz2, ImageDimension);
  fixedlms.fill(0);
  movinglms.fill(0);

  if (sz1 <= 0 || sz2 <= 0)
  {
    return;
  }

  DisplacementFieldTypePointer lmField = this->m_DerivativeFixedField;
  if (!whichdirection)
  {
    lmField = this->m_DerivativeMovingField;
  }

  float inweight = weight;
  this->m_LandmarkEnergy = 0.0;

  std::cout << " sz1 " << sz1 << " sz2 " << sz2 << std::endl;
  // if whichdirection is true, then the fixed direction, else moving
  for (unsigned long ii = 0; ii < sz1; ii++)
  {
    PointType     fixedpoint;
    PointDataType fixedlabel = 0;
    if (whichdirection)
    {
      this->m_FixedPointSet->GetPoint(ii, &fixedpoint);
    }
    else
    {
      this->m_MovingPointSet->GetPoint(ii, &fixedpoint);
    }
    if (whichdirection)
    {
      this->m_FixedPointSet->GetPointData(ii, &fixedlabel);
    }
    else
    {
      this->m_MovingPointSet->GetPointData(ii, &fixedlabel);
    }

    float          min = 1.e9;
    ImagePointType fpt;
    IndexType      oindex;
    for (int j = 0; j < ImageDimension; j++)
    {
      fpt[j] = fixedpoint[j];
      fixedlms(ii, j) = fpt[j];
    }
    bool convok = false;
    convok = this->GetFixedImage()->TransformPhysicalPointToIndex(fpt, oindex);
    if (!convok)
    {
      std::cout << " fpt " << fpt << std::endl;
    }
    // if whichdirection is true, then the fixed direction, else moving
    for (unsigned long jj = 0; jj < sz2; jj++)
    {
      VectorType distance;
      IndexType  fixedindex;
      IndexType  movingindex;
      PointType  movingpoint;
      if (whichdirection)
      {
        this->m_MovingPointSet->GetPoint(jj, &movingpoint);
      }
      else
      {
        this->m_FixedPointSet->GetPoint(jj, &movingpoint);
      }

      ImagePointType mpt;
      for (int j = 0; j < ImageDimension; j++)
      {
        mpt[j] = movingpoint[j];
        movinglms(jj, j) = movingpoint[j];
      }

      if (ii == sz1 - 2 && jj == sz2 - 2)
      {
        std::cout << " fpt " << fpt << " mpt " << mpt << std::endl;
      }

      this->GetMovingImage()->TransformPhysicalPointToIndex(mpt, movingindex);
      double prob = 0;
      if (convok)
      {
        double     mag = 0.0;
        VectorType force;
        for (int j = 0; j < ImageDimension; j++)
        {
          distance[j] = movingpoint[j] - fixedpoint[j];
          mag += static_cast<double>(itk::Math::sqr(distance[j] / spacing[j]));
          force[j] = distance[j] * inweight;
        }
        auto sigma = static_cast<double>(this->m_FixedPointSetSigma);
        if (!whichdirection)
        {
          sigma = static_cast<double>(this->m_MovingPointSetSigma);
        }
        // KW -- rename 'prob' to '_prob' because of 'shadow variable' warning.
        double _prob = 1.0 / sqrt(3.14186 * 2.0 * sigma * sigma) * exp(-1.0 * mag / (2.0 * sigma * sigma));
        force = force * _prob;
        PointDataType movinglabel = 0;
        if (whichdirection)
        {
          this->m_MovingPointSet->GetPointData(jj, &movinglabel);
        }
        else
        {
          this->m_FixedPointSet->GetPointData(jj, &movinglabel);
        }
        //        if (ii == 2 && jj==2) std::cout << "_prob " << _prob << " sigma " << sigma << "  " << mag << " fl " <<
        // fixedlabel << " ml " << movinglabel << std::endl;
        if (fixedlabel != movinglabel)
        {
          _prob = 0;
        }
        // || fixedlabel !=4) _prob=0;
        mag = sqrt(mag);
      }
      EucDist(ii, jj) = prob;
    }
    if (min < static_cast<float>(1.e5))
    {
      this->m_LandmarkEnergy += min;
    }
  }

  MatrixType sinkhorn = EucDist;
  for (unsigned int iter = 0; iter < 1; iter++)
  {
    for (unsigned int jj = 0; jj < sz2; jj++)
    {
      double total = 0.0;
      for (unsigned int ii = 0; ii < sz1; ii++)
      {
        total += sinkhorn(ii, jj);
      }
      if (total <= itk::NumericTraits<double>::ZeroValue())
      {
        total = itk::NumericTraits<double>::OneValue();
      }
      for (unsigned int ii = 0; ii < sz1; ii++)
      {
        sinkhorn(ii, jj) /= total;
      }
    }
    for (unsigned int ii = 0; ii < sz1; ii++)
    {
      double total = 0;
      for (unsigned int jj = 0; jj < sz2; jj++)
      {
        total += sinkhorn(ii, jj);
      }
      if (total <= itk::NumericTraits<double>::ZeroValue())
      {
        total = itk::NumericTraits<double>::OneValue();
      }
      for (unsigned int jj = 0; jj < sz2; jj++)
      {
        sinkhorn(ii, jj) /= total;
      }
    }
  }

  MatrixType resultlms = sinkhorn * movinglms;
  MatrixType difflms = fixedlms - resultlms;
  VectorType sforce;
  sforce.Fill(0);
  float energy = 0, maxerr = 0;
  for (unsigned long ii = 0; ii < sz1; ii++)
  {
    VectorType distance;
    distance.Fill(0);
    PointType     movingpoint;
    PointType     fixedpoint;
    PointDataType fixedlabel = 0;
    if (whichdirection)
    {
      this->m_FixedPointSet->GetPoint(ii, &fixedpoint);
    }
    else
    {
      this->m_MovingPointSet->GetPoint(ii, &fixedpoint);
    }
    if (whichdirection)
    {
      this->m_FixedPointSet->GetPointData(ii, &fixedlabel);
    }
    else
    {
      this->m_MovingPointSet->GetPointData(ii, &fixedlabel);
    }

    ImagePointType mpt;
    ImagePointType fpt;
    for (int j = 0; j < ImageDimension; j++)
    {
      fpt[j] = fixedpoint[j];
    }
    for (int j = 0; j < ImageDimension; j++)
    {
      mpt[j] = resultlms(ii, j);
    }

    bool      convok = false;
    IndexType fixedindex;
    convok = this->GetFixedImage()->TransformPhysicalPointToIndex(fpt, fixedindex);
    if (convok)
    {
      double     mag = 0.0;
      VectorType force;
      for (int j = 0; j < ImageDimension; j++)
      {
        distance[j] = mpt[j] - fixedpoint[j];
        mag += static_cast<double>(itk::Math::sqr(distance[j] / spacing[j]));
        force[j] = distance[j] * inweight;
      }
      auto sigma = static_cast<double>(this->m_FixedPointSetSigma);
      if (!whichdirection)
      {
        sigma = static_cast<double>(this->m_MovingPointSetSigma);
      }
      double prob = 1.0 / sqrt(3.14186 * 2.0 * sigma * sigma) * exp(-1.0 * mag / (2.0 * sigma * sigma));
      force = force * prob;

      //      if (fixedlabel !=4 ) force.Fill(0);
      /*
            if (mag > 50)
          {
          float tot=0;
          for (int k=0; k < sz2; k++) { tot+=sinkhorn(ii,k); }
          std::cout << "TOT " << tot << std::endl;
          force.Fill(0);
          mag=0;
          }
      */
      if (mag > static_cast<double>(maxerr))
      {
        maxerr = static_cast<float>(mag);
      }
      energy += static_cast<float>(mag);
      std::cout << " ii " << ii << " force " << force << " mag " << std::sqrt(mag) << " mpt " << mpt << " fpt "
                << fixedpoint << " nrg " << energy / (float)ii << std::endl;
      lmField->SetPixel(fixedindex, force + lmField->GetPixel(fixedindex));
    }
    //    lmField->SetPixel(fixedindex,sforce);
  }
  //  std::cout <<  " max " << maxerr << std::endl;
  this->m_LandmarkEnergy = energy / (float)sz1;
  this->m_Energy = this->m_LandmarkEnergy;
}

/*
 * Set the function state values before each iteration
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  InitializeIteration()
{
  //  std::cout << " INIT ITER " << std::endl;
  if (!this->GetMovingImage() || !this->GetFixedImage())
  {
    itkExceptionMacro(<< "MovingImage, FixedImage  not set");
  }
  // cache fixed image information
  m_FixedImageSpacing = this->GetFixedImage()->GetSpacing();
  m_FixedImageOrigin = this->GetFixedImage()->GetOrigin();

  // compute the normalizer
  m_Normalizer = 0.0;
  for (unsigned int k = 0; k < ImageDimension; k++)
  {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
  }
  m_Normalizer /= static_cast<double>(ImageDimension);

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage(this->GetFixedImage());
  m_MovingImageGradientCalculator->SetInputImage(this->GetMovingImage());

  // initialize metric computation variables
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_SumOfSquaredChange = 0.0;

  typename DisplacementFieldType::PixelType zero;
  zero.Fill(0);
  this->m_DerivativeFixedField =
    AllocImage<DisplacementFieldType>(this->GetFixedImage()->GetLargestPossibleRegion(), zero);
  this->m_DerivativeFixedField->SetSpacing(this->GetFixedImage()->GetSpacing());
  this->m_DerivativeFixedField->SetOrigin(this->GetFixedImage()->GetOrigin());
  this->m_DerivativeMovingField =
    AllocImage<DisplacementFieldType>(this->GetMovingImage()->GetLargestPossibleRegion(), zero);
  this->m_DerivativeMovingField->SetSpacing(this->GetMovingImage()->GetSpacing());
  this->m_DerivativeMovingField->SetOrigin(this->GetMovingImage()->GetOrigin());

  // acquire labels
  if (this->m_LabelSet.empty())
  {
    this->m_LabelSet.clear();
    unsigned long sz1 = this->m_FixedPointSet->GetNumberOfPoints();
    std::cout << " NPTS " << sz1 << std::endl;
    for (unsigned long ii = 0; ii < sz1; ii++)
    {
      PointType     fixedpoint;
      PointDataType label = 0;
      this->m_FixedPointSet->GetPoint(ii, &fixedpoint);
      this->m_FixedPointSet->GetPointData(ii, &label);
      if (label > 0)
      {
        if (find(this->m_LabelSet.begin(), this->m_LabelSet.end(), label) == this->m_LabelSet.end())
        {
          this->m_LabelSet.push_back(label);
        }
      }
    }
  }
  else
  {
    std::cout << " #of Label Values to match " << this->m_LabelSet.size() << std::endl;
  }

  this->m_bpoints = BSplinePointSetType::New();
  this->m_bpoints->Initialize();
  this->m_bweights = BSplineWeightsType::New();
  this->m_bweights->Initialize();
  this->m_bcount = 0;

  // unsigned int                          lct = 0;
  typename LabelSetType::const_iterator it;
  for (it = this->m_LabelSet.begin(); it != this->m_LabelSet.end(); ++it)
  {
    // lct++;
    auto label = (PointDataType)*it;
    //     std::cout << " doing label " << label << std::endl;
    this->SetUpKDTrees(label);
    bool dobsp = false;
    //    if (lct ==  this->m_LabelSet.size()  ) dobsp=true;
    this->FastExpectationLandmarkField(1.0, true, label, dobsp);
    this->FastExpectationLandmarkField(1.0, false, label, dobsp);
  }
  // follow up with BSpline if dospb is true .
}

/*
 * Compute update at a specify neighbourhood
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
typename ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  PixelType
  ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::ComputeUpdate(
    const NeighborhoodType & it,
    void * /* gd */,
    const FloatOffsetType & itkNotUsed(offset))
{
  IndexType index = it.GetIndex();
  PixelType update = this->m_DerivativeFixedField->GetPixel(index);

  if (this->m_Iterations > this->m_UseSymmetricMatching)
  {
    update.Fill(0);
  }
  return update;
}

/*
 * Compute update at a specify neighbourhood
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
typename ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  PixelType
  ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
    ComputeUpdateInv(const NeighborhoodType & it, void * /* gd */, const FloatOffsetType & itkNotUsed(offset))
{
  IndexType index = it.GetIndex();

  return this->m_DerivativeMovingField->GetPixel(index);
}

/*
 * Update the metric and release the per-thread-global data.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  ReleaseGlobalDataPointer(void * gd) const
{
  auto * globalData = (GlobalDataStruct *)gd;

  m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
  m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
  m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
  if (m_NumberOfPixelsProcessed)
  {
    m_Metric = m_SumOfSquaredDifference / static_cast<double>(m_NumberOfPixelsProcessed);
    m_RMSChange = std::sqrt(m_SumOfSquaredChange / static_cast<double>(m_NumberOfPixelsProcessed));
  }

  delete globalData;
}

/*
 * Set the function state values before each iteration
 */

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::SetUpKDTrees(
  long whichlabel)
{
  // convert this->m_FixedPointSet to a sample type
  this->m_FixedSamplePoints = SampleType::New();
  this->m_FixedSamplePoints->SetMeasurementVectorSize(MeasurementDimension);
  MeasurementVectorType mv;
  unsigned int          bucketsize = 4;
  unsigned int          npts = this->m_FixedPointSet->GetNumberOfPoints();
  //  std::cout << " NP MOV " << npts << std::endl;
  for (unsigned int i = 0; i < npts; i++)
  {
    PointType  fixedpoint;
    const bool validFixedPoint = this->m_FixedPointSet->GetPoint(i, &fixedpoint);
    if (!validFixedPoint)
    {
      itkExceptionMacro(<< "Invalid FixedPoint Requested at " << i);
    }
    else
    {
      PointDataType fixedlabel = 0;
      this->m_FixedPointSet->GetPointData(i, &fixedlabel);
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        mv[d] = fixedpoint[d];
      }
      // mv[ImageDimension]=(float) fixedlabel*1.e6;
      if (fixedlabel == whichlabel)
      {
        this->m_FixedSamplePoints->PushBack(mv);
      }
    }
  }
  this->m_FixedKdTreeGenerator = TreeGeneratorType::New();
  this->m_FixedKdTreeGenerator->SetSample(this->m_FixedSamplePoints);
  this->m_FixedKdTreeGenerator->SetBucketSize(bucketsize);
  this->m_FixedKdTreeGenerator->Update();

  this->m_MovingSamplePoints = SampleType::New();
  this->m_MovingSamplePoints->SetMeasurementVectorSize(ImageDimension);
  npts = this->m_MovingPointSet->GetNumberOfPoints();
  //  std::cout << " NP MOV " << npts << std::endl;
  for (unsigned int i = 0; i < npts; i++)
  {
    PointType  movingpoint;
    const bool validMovingPoint = this->m_MovingPointSet->GetPoint(i, &movingpoint);
    if (!validMovingPoint)
    {
      itkExceptionMacro(<< "Invalid MovingPoint Requested at " << i);
    }
    else
    {
      PointDataType movinglabel = 0;
      this->m_MovingPointSet->GetPointData(i, &movinglabel);
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        mv[d] = movingpoint[d];
      }
      // mv[ImageDimension]=(float) movinglabel*1.e6;
      if (movinglabel == whichlabel)
      {
        this->m_MovingSamplePoints->PushBack(mv);
      }
    }
  }
  this->m_MovingKdTreeGenerator = TreeGeneratorType::New();
  this->m_MovingKdTreeGenerator->SetSample(this->m_MovingSamplePoints);
  this->m_MovingKdTreeGenerator->SetBucketSize(bucketsize);
  this->m_MovingKdTreeGenerator->Update();
}

/*
 * Set the function state values before each iteration
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
void
ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>::
  FastExpectationLandmarkField(float weight, bool whichdirection, long /* whichlabel */, bool dobspline)
{
  /**
   * BSpline typedefs
   */
  /** Typedefs for B-spline filter */

  unsigned int m_SplineOrder = 3;
  unsigned int m_NumberOfBLevels = 5;
  ArrayType    m_MeshResolution;

  m_MeshResolution.Fill(1);
  unsigned int PointDimension = ImageDimension;

  SpacingType spacing = this->GetFixedImage()->GetSpacing();

  typename TreeGeneratorType::Pointer fkdtree;
  typename TreeGeneratorType::Pointer mkdtree;
  if (whichdirection)
  {
    mkdtree = this->m_MovingKdTreeGenerator;
    fkdtree = this->m_FixedKdTreeGenerator;
  }
  else
  {
    fkdtree = this->m_MovingKdTreeGenerator;
    mkdtree = this->m_FixedKdTreeGenerator;
  }

  unsigned long sz1 = fkdtree->GetOutput()->Size();
  unsigned long sz2 = mkdtree->GetOutput()->Size();

  //  std::cout << " s1 " << sz1 << " s2 " << sz2 << std::endl;

  if (sz1 <= 0 || sz2 <= 0)
  {
    return;
  }

  DisplacementFieldTypePointer lmField = this->m_DerivativeFixedField;
  if (!whichdirection)
  {
    lmField = this->m_DerivativeMovingField;
  }

  unsigned int KNeighbors = this->m_KNeighborhood;
  if (KNeighbors > sz2)
  {
    KNeighbors = sz2;
  }
  float inweight = weight;
  this->m_LandmarkEnergy = 0.0;
  //  float max=0;

  vnl_vector<double> probabilities(KNeighbors);
  probabilities.fill(0);

  float energy = 0, maxerr = 0;
  for (unsigned long ii = 0; ii < sz1; ii++)
  {
    VectorType distance;
    distance.Fill(0);
    MeasurementVectorType fixedpoint = fkdtree->GetOutput()->GetMeasurementVector(ii);
    ImagePointType        mpt;
    mpt.Fill(0);
    ImagePointType fpt;
    fpt.Fill(0);
    for (unsigned int j = 0; j < ImageDimension; j++)
    {
      fpt[j] = fixedpoint[j];
    }
    //    float err=0;
    bool      convok = false;
    IndexType fixedindex;
    convok = this->GetFixedImage()->TransformPhysicalPointToIndex(fpt, fixedindex);
    //    std::cout << " Orig " << this->GetFixedImage()->GetOrigin() << " ind " << fixedindex << " pt " << fpt <<
    // std::endl;
    if (convok)
    {
      VectorType force;
      force.Fill(0);
      typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
      mkdtree->GetOutput()->Search(fixedpoint, static_cast<unsigned int>(KNeighbors), neighbors);
      double probtotal = 0.0;
      for (unsigned int dd = 0; dd < KNeighbors; dd++)
      {
        unsigned long         wpt = neighbors[dd];
        MeasurementVectorType npt = mkdtree->GetOutput()->GetMeasurementVector(wpt);
        double                _mag = 0;
        for (unsigned int qq = 0; qq < ImageDimension; qq++)
        {
          _mag += static_cast<double>(itk::Math::sqr(fixedpoint[qq] - npt[qq]));
        }
        auto sigma = static_cast<double>(this->m_FixedPointSetSigma);
        if (!whichdirection)
        {
          sigma = static_cast<double>(this->m_MovingPointSetSigma);
        }
        double prob = 1.0 / sqrt(3.14186 * 2.0 * sigma * sigma) * exp(-1.0 * _mag / (2.0 * sigma * sigma));
        probtotal += prob;
        probabilities(dd) = prob;
      }
      if (probtotal > 0)
      {
        for (unsigned int dd = 0; dd < KNeighbors; dd++)
        {
          unsigned long         wpt = neighbors[dd];
          MeasurementVectorType npt = mkdtree->GetOutput()->GetMeasurementVector(wpt);
          double                pp = probabilities(dd) / probtotal;
          if (pp > 0)
          {
            for (unsigned int j = 0; j < ImageDimension; j++)
            {
              mpt[j] += static_cast<typename ImagePointType::CoordRepType>(pp) *
                        static_cast<typename ImagePointType::CoordRepType>(npt[j]);
            }
          }
          //
          //    if (ii % 245 && pp > 1.e-3) std::cout << " prob " << pp <<  " mpt " << mpt << " dd " << dd <<" wpt " <<
          //    wpt
          // << "
          // movinpoint " << movingpoint << " ptot " << probtotal <<  std::endl;
        }
      }

      double                                  mag = 0.0;
      typename BSplinePointSetType::PointType bpoint;
      for (unsigned int j = 0; j < ImageDimension; j++)
      {
        distance[j] = static_cast<typename VectorType::ComponentType>(mpt[j]) -
                      static_cast<typename VectorType::ComponentType>(fixedpoint[j]);
        mag += static_cast<double>(itk::Math::sqr(static_cast<double>(distance[j]) / static_cast<double>(spacing[j])));
        force[j] = distance[j] * inweight;
        bpoint[j] = fixedpoint[j];
      }
      auto sigma = static_cast<double>(this->m_FixedPointSetSigma);
      if (!whichdirection)
      {
        sigma = static_cast<double>(this->m_MovingPointSetSigma);
      }
      double prob = 1.0 / sqrt(3.14186 * 2.0 * sigma * sigma) * exp(-1.0 * mag / (2.0 * sigma * sigma));
      force = force * prob;

      //
      this->m_bpoints->SetPoint(this->m_bcount, bpoint);
      this->m_bpoints->SetPointData(this->m_bcount, distance);
      float bwt = 1;
      this->m_bweights->InsertElement(this->m_bcount, static_cast<typename BSplineWeightsType::Element>(bwt));
      this->m_bcount++;

      mag = sqrt(mag);
      if (mag > static_cast<double>(maxerr))
      {
        maxerr = static_cast<float>(mag);
      }
      energy += static_cast<float>(mag);
      lmField->SetPixel(fixedindex, force + lmField->GetPixel(fixedindex));
    }
  }
  //  std::cout <<  " max " << maxerr << std::endl;
  this->m_LandmarkEnergy = energy / (float)sz1;
  this->m_Energy = this->m_LandmarkEnergy;

  if (dobspline)
  {
    /**
     * Calculate derivative field with respect to the moving points
     */
    {
      // std::cout << " start bsp " << std::endl;

      typename BSplineFilterType::ArrayType nlevels;
      typename BSplineFilterType::ArrayType ncps;

      nlevels.Fill(m_NumberOfBLevels);
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        ncps[d] = m_MeshResolution[d] + m_SplineOrder;
      }

      typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
      bspliner->SetInput(this->m_bpoints);
      bspliner->SetOrigin(this->GetMovingImage()->GetOrigin());
      bspliner->SetSpacing(this->GetMovingImage()->GetSpacing());
      bspliner->SetSize(this->GetMovingImage()->GetLargestPossibleRegion().GetSize());
      bspliner->SetNumberOfLevels(nlevels);
      bspliner->SetSplineOrder(m_SplineOrder);
      bspliner->SetNumberOfControlPoints(ncps);
      bspliner->SetGenerateOutputImage(true);
      bspliner->SetPointWeights(this->m_bweights);
      bspliner->Update();

      lmField = bspliner->GetOutput();

      /**
       * Now calculate the distances after matching.

      ItM = this->m_MovingPointSet->GetPoints()->Begin();
      ItMD = this->m_MovingPointSet->GetPointData()->Begin();
      ItF = this->m_FixedPointSet->GetPoints()->Begin();
      ItFD = this->m_FixedPointSet->GetPointData()->Begin();

      typename BSplineWeightsType::ConstIterator ItW = weights->Begin();

      while( ItM != this->m_MovingPointSet->GetPoints()->End() &&
        ItF != this->m_FixedPointSet->GetPoints()->End() )
        {
        typename BSplinePointSetType::PixelType vector;
        bspliner->EvaluateAtPoint( ItM.Value(), vector );

        RealType distance = 0.0;
        for( unsigned int d = 0; d < PointDimension; d++ )
          {
          distance += itk::Math::sqr ( ItM.Value()[d] + vector[d]
            - ItF.Value()[d] );
          }
        this->m_Energy += ItW.Value() * std::sqrt( distance );

        ++ItF;
        ++ItFD;
        ++ItM;
        ++ItMD;
        ++ItW;
        }
       */
    }
  }
}
} // end namespace itk

#endif
