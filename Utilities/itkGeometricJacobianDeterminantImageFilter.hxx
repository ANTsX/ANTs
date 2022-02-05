/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGeometricJacobianDeterminantImageFilter_hxx
#define _itkGeometricJacobianDeterminantImageFilter_hxx


#include "itkContinuousIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"
#include "itkCastImageFilter.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

#include "vnl/vnl_cross.h"

namespace itk
{

template <typename TInputImage, typename TRealType, typename TOutputImage>
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::GeometricJacobianDeterminantImageFilter()
{
  this->m_Interpolator = nullptr;
  this->m_UndisplacedVolume = 0.0;
  this->DynamicMultiThreadingOff();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::GenerateInputRequestedRegion()
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

  this->m_NeighborhoodRadius.Fill(1);

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
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::BeforeThreadedGenerateData()
{
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

  this->m_Interpolator = InterpolatorType::New();
  this->m_Interpolator->SetInputImage(this->m_RealValuedInputImage);

  PointType origin(0.0);

  if (ImageDimension == 2)
  {
    this->InitializeTriangularDeltaPoints();

    PointType pointA = origin + this->m_DeltaTriangularPointA;
    PointType pointB = origin + this->m_DeltaTriangularPointB;
    PointType pointC = origin + this->m_DeltaTriangularPointC;

    this->m_UndisplacedVolume = this->CalculateTriangularArea(pointA, pointB, pointC);
  }
  else if (ImageDimension == 3)
  {
    this->InitializeTetrahedralDeltaPoints();

    PointType pointA = origin + this->m_DeltaTetrahedralPointA;
    PointType pointB = origin + this->m_DeltaTetrahedralPointB;
    PointType pointC = origin + this->m_DeltaTetrahedralPointC;
    PointType pointD = origin + this->m_DeltaTetrahedralPointD;

    this->m_UndisplacedVolume = this->CalculateTetrahedralVolume(pointA, pointB, pointC, pointD);
  }
  else
  {
    itkExceptionMacro("Computations are only valid for ImageDimension = 2 or 3");
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::InitializeTetrahedralDeltaPoints()
{
  // cf http://en.wikipedia.org/wiki/Tetrahedron#Formulas_for_a_regular_tetrahedron

  typename InputImageType::PointType originPoint;

  ContinuousIndex<RealType, ImageDimension> cidx;
  cidx.Fill(0.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, originPoint);

  typename InputImageType::PointType deltaPoint;

  cidx[0] = 0.5;
  cidx[1] = 0.0;
  cidx[2] = -0.5 / std::sqrt(2.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTetrahedralPointA = deltaPoint - originPoint;

  cidx[0] = -0.5;
  cidx[1] = 0.0;
  cidx[2] = -0.5 / std::sqrt(2.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTetrahedralPointB = deltaPoint - originPoint;

  cidx[0] = 0.0;
  cidx[1] = 0.5;
  cidx[2] = 0.5 / std::sqrt(2.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTetrahedralPointC = deltaPoint - originPoint;

  cidx[0] = 0.0;
  cidx[1] = -0.5;
  cidx[2] = 0.5 / std::sqrt(2.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTetrahedralPointD = deltaPoint - originPoint;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::InitializeTriangularDeltaPoints()
{
  typename InputImageType::PointType originPoint;

  ContinuousIndex<RealType, ImageDimension> cidx;
  cidx.Fill(0.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, originPoint);

  typename InputImageType::PointType deltaPoint;

  cidx[0] = 0.0;
  cidx[1] = 0.25 * std::sqrt(3.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTriangularPointA = deltaPoint - originPoint;

  cidx[0] = -0.5;
  cidx[1] = -0.25 * std::sqrt(3.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTriangularPointB = deltaPoint - originPoint;

  cidx[0] = 0.5;
  cidx[1] = -0.25 * std::sqrt(3.0);

  this->m_RealValuedInputImage->TransformContinuousIndexToPhysicalPoint(cidx, deltaPoint);

  this->m_DeltaTriangularPointC = deltaPoint - originPoint;
}


template <typename TInputImage, typename TRealType, typename TOutputImage>
void
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  ZeroFluxNeumannBoundaryCondition<RealVectorImageType> nbc;
  ConstNeighborhoodIteratorType                         bit;
  ImageRegionIteratorWithIndex<TOutputImage>            it;

  typename OutputImageType::Pointer outputPtr = this->GetOutput();

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>                        bC;
  faceList = bC(dynamic_cast<const RealVectorImageType *>(this->m_RealValuedInputImage.GetPointer()),
                outputRegionForThread,
                this->m_NeighborhoodRadius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<RealVectorImageType>::FaceListType::iterator fit;

  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process each of the data set faces.  The iterator is reinitialized on each
  // face so that it can determine whether or not to check for boundary
  // conditions.
  for (fit = faceList.begin(); fit != faceList.end(); ++fit)
  {
    bit = ConstNeighborhoodIteratorType(
      this->m_NeighborhoodRadius, dynamic_cast<const RealVectorImageType *>(m_RealValuedInputImage.GetPointer()), *fit);
    it = ImageRegionIteratorWithIndex<TOutputImage>(outputPtr, *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    while (!bit.IsAtEnd())
    {
      typename InputImageType::IndexType index = it.GetIndex();
      PointType                          imagePoint;
      outputPtr->TransformIndexToPhysicalPoint(index, imagePoint);

      RealType displacedVolume = 0.0;
      if (ImageDimension == 2)
      {
        typename InterpolatorType::PointType pointA;
        pointA.CastFrom(imagePoint + this->m_DeltaTriangularPointA);
        typename InterpolatorType::PointType pointB;
        pointB.CastFrom(imagePoint + this->m_DeltaTriangularPointB);
        typename InterpolatorType::PointType pointC;
        pointC.CastFrom(imagePoint + this->m_DeltaTriangularPointC);

        PointType displacedPointA = pointA + this->m_Interpolator->Evaluate(pointA);
        PointType displacedPointB = pointB + this->m_Interpolator->Evaluate(pointB);
        PointType displacedPointC = pointC + this->m_Interpolator->Evaluate(pointC);

        displacedVolume = this->CalculateTriangularArea(displacedPointA, displacedPointB, displacedPointC);
      }
      else if (ImageDimension == 3)
      {
        typename InterpolatorType::PointType pointA;
        pointA.CastFrom(imagePoint + this->m_DeltaTetrahedralPointA);
        typename InterpolatorType::PointType pointB;
        pointB.CastFrom(imagePoint + this->m_DeltaTetrahedralPointB);
        typename InterpolatorType::PointType pointC;
        pointC.CastFrom(imagePoint + this->m_DeltaTetrahedralPointC);
        typename InterpolatorType::PointType pointD;
        pointD.CastFrom(imagePoint + this->m_DeltaTetrahedralPointD);

        PointType displacedPointA = pointA + this->m_Interpolator->Evaluate(pointA);
        PointType displacedPointB = pointB + this->m_Interpolator->Evaluate(pointB);
        PointType displacedPointC = pointC + this->m_Interpolator->Evaluate(pointC);
        PointType displacedPointD = pointD + this->m_Interpolator->Evaluate(pointD);
        displacedVolume =
          this->CalculateTetrahedralVolume(displacedPointA, displacedPointB, displacedPointC, displacedPointD);
      }
      RealType volumeDifferential = displacedVolume / (this->m_UndisplacedVolume);
      it.Set(volumeDifferential);
      ++bit;
      ++it;
      progress.CompletedPixel();
    }
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::RealType
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::CalculateTetrahedralVolume(PointType a,
                                                                                                          PointType b,
                                                                                                          PointType c,
                                                                                                          PointType d)
{
  vnl_vector<double> ad = (a - d).GetVnlVector();
  vnl_vector<double> bd = (b - d).GetVnlVector();
  vnl_vector<double> cd = (c - d).GetVnlVector();
  vnl_vector<double> bdxcd = vnl_cross_3d(bd, cd);
  RealType           volume = itk::Math::abs(ad[0] * bdxcd[0] + ad[1] * bdxcd[1] + ad[2] * bdxcd[2]) / 6.0;
  return volume;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::RealType
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::CalculateTriangularArea(PointType a,
                                                                                                       PointType b,
                                                                                                       PointType c)
{
  RealVectorType ab = (a - b);
  RealVectorType ac = (a - c);

  RealType area = 0.5 * itk::Math::abs(ab[0] * ac[1] - ac[0] * ab[1]);

  return area;
}


template <typename TInputImage, typename TRealType, typename TOutputImage>
void
GeometricJacobianDeterminantImageFilter<TInputImage, TRealType, TOutputImage>::PrintSelf(std::ostream & os,
                                                                                         Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "m_RealValuedInputImage = " << m_RealValuedInputImage.GetPointer() << std::endl;
}

} // end namespace itk

#endif
