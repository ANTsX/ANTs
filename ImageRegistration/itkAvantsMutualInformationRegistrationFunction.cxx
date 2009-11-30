/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAvantsMutualInformationRegistrationFunction.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/08 15:14:48 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAvantsMutualInformationRegistrationFunction_txx
#define _itkAvantsMutualInformationRegistrationFunction_txx

#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::AvantsMutualInformationRegistrationFunction()
{
  this->Superclass::m_NormalizeGradient = true;
  this->m_NumberOfSpatialSamples = 500;
  this->m_NumberOfHistogramBins = 50;

  //  this->SetComputeGradient(false); // don't use the default gradient for now

  this->m_InterpolatorIsBSpline = false;
  this->m_TransformIsBSpline    = false;

  // Initialize PDFs to NULL
  m_JointPDF = NULL;
  m_JointPDFDerivatives = NULL;

  m_OpticalFlow = false;
  typename TransformType::Pointer transformer = TransformType::New();
  this->SetTransform(transformer);

  typename BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
  this->SetInterpolator(interpolator);

  m_FixedImageMask = NULL;
  m_MovingImageMask = NULL;

  // Initialize memory
  m_MovingImageNormalizedMin = 0.0;
  m_FixedImageNormalizedMin = 0.0;
  m_MovingImageTrueMin = 0.0;
  m_MovingImageTrueMax = 0.0;
  m_FixedImageBinSize = 0.0;
  m_MovingImageBinSize = 0.0;
  m_CubicBSplineDerivativeKernel = NULL;
  m_BSplineInterpolator = NULL;
  m_DerivativeCalculator = NULL;
  m_NumParametersPerDim = 0;
  m_NumBSplineWeights = 0;
  m_BSplineTransform = NULL;
  m_NumberOfParameters = ImageDimension;

  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  m_MovingImageGradientCalculator = GradientCalculatorType::New();
  this->m_Padding = 4;

  typename DefaultInterpolatorType::Pointer interp =  DefaultInterpolatorType::New();
  typename DefaultInterpolatorType::Pointer interp2 = DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType *>(
      interp.GetPointer() );
  m_FixedImageInterpolator = static_cast<InterpolatorType *>(
      interp2.GetPointer() );
  m_Interpolator = static_cast<InterpolatorType *>(
      interp.GetPointer() );

  // std::cout << " done declaring " << std::endl;
  this->m_RobustnessParameter = -1.e19;
}

/**
 * Print out internal information about this class
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfSpatialSamples: ";
  os << m_NumberOfSpatialSamples << std::endl;
  os << indent << "NumberOfHistogramBins: ";
  os << m_NumberOfHistogramBins << std::endl;

  // Debugging information
  os << indent << "NumberOfParameters: ";
  os << m_NumberOfParameters << std::endl;
  os << indent << "FixedImageNormalizedMin: ";
  os << m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << m_MovingImageNormalizedMin << std::endl;
  os << indent << "MovingImageTrueMin: ";
  os << m_MovingImageTrueMin << std::endl;
  os << indent << "MovingImageTrueMax: ";
  os << m_MovingImageTrueMax << std::endl;
  os << indent << "FixedImageBinSize: ";
  os << m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << m_MovingImageBinSize << std::endl;
  os << indent << "InterpolatorIsBSpline: ";
  os << m_InterpolatorIsBSpline << std::endl;
  os << indent << "TransformIsBSpline: ";
  os << m_TransformIsBSpline << std::endl;
}

/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::InitializeIteration()
{
//  std::cout <<" Obust P " << this->m_RobustnessParameter;
//  std::cout <<" bins " << m_NumberOfHistogramBins << std::endl;
  m_CubicBSplineKernel = CubicBSplineFunctionType::New();
  m_CubicBSplineDerivativeKernel = CubicBSplineDerivativeFunctionType::New();
  this->m_Energy = 0;
  pdfinterpolator = pdfintType::New();
  dpdfinterpolator = dpdfintType::New();
  pdfinterpolator2 = pdfintType2::New();
  pdfinterpolator3 = pdfintType2::New();
  m_DerivativeCalculator = DerivativeFunctionType::New();

//  this->ComputeMetricImage();
/*
    bool makenewimage=false;
    typedef ImageRegionIteratorWithIndex<TFixedImage> ittype;
    TFixedImage* img =const_cast<TFixedImage *>(this->m_FixedImage.GetPointer());
    typename TFixedImage::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();

    if (!this->m_MetricImage )makenewimage=true;
    else if (imagesize[0] != this->m_MetricImage->GetLargestPossibleRegion().GetSize()[0])
      makenewimage = true;
    else this->m_MetricImage->FillBuffer(0);
    if (makenewimage)
    {
      this->m_MetricImage = TFixedImage::New();
      this->m_MetricImage->SetLargestPossibleRegion(img->GetLargestPossibleRegion()  );
      this->m_MetricImage->SetBufferedRegion(img->GetLargestPossibleRegion());
      this->m_MetricImage->SetSpacing(img->GetSpacing());
      this->m_MetricImage->SetOrigin(img->GetOrigin());
      this->m_MetricImage->Allocate();
      ittype it(this->m_MetricImage,this->m_MetricImage->GetLargestPossibleRegion().GetSize());
      for( it.GoToBegin(); !it.IsAtEnd(); ++it ) it.Set(0);
    }
*/
  m_FixedImageGradientCalculator->SetInputImage( this->m_FixedImage );
  m_MovingImageGradientCalculator->SetInputImage( this->m_MovingImage );
  m_FixedImageInterpolator->SetInputImage( this->m_FixedImage );
  m_Interpolator->SetInputImage( this->m_MovingImage );

  m_FixedImageSpacing    = this->m_FixedImage->GetSpacing();
  m_FixedImageOrigin     = this->m_FixedImage->GetOrigin();
  m_Normalizer      = 0.0;
  m_NumberOfSpatialSamples = 1;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
    m_NumberOfSpatialSamples *= this->m_FixedImage->GetLargestPossibleRegion().GetSize()[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );

  // std::cout << " 6 ";

  /**
   * Compute binsize for the histograms.
   *
   * The binsize for the image intensities needs to be adjusted so that
   * we can avoid dealing with boundary conditions using the cubic
   * spline as the Parzen window.  We do this by increasing the size
   * of the bins so that the joint histogram becomes "padded" at the
   * borders. Because we are changing the binsize,
   * we also need to shift the minimum by the padded amount in order to
   * avoid minimum values filling in our padded region.
   *
   * Note that there can still be non-zero bin values in the padded region,
   * it's just that these bins will never be a central bin for the Parzen
   * window.
   *
  double fixedImageMax = 1.0;
  double fixedImageMin = 0.0;
  double movingImageMax = 1.0;
  double movingImageMin = 0.0;
  m_MovingImageTrueMin = movingImageMin;
  m_MovingImageTrueMax = movingImageMax;

   */

  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();

  typedef ImageRegionConstIterator<MovingImageType> MovingIteratorType;
  MovingIteratorType movingImageIterator(
    this->m_MovingImage, this->m_MovingImage->GetBufferedRegion() );
  for( movingImageIterator.GoToBegin();
       !movingImageIterator.IsAtEnd(); ++movingImageIterator )
    {
    bool takesample = true;

    if( this->m_FixedImageMask )
      {
      if( this->m_FixedImageMask->GetPixel( movingImageIterator.GetIndex() ) < 1.e-6 )
        {
        takesample = false;
        }
      }

    if( takesample )
      {
      double sample = static_cast<double>( movingImageIterator.Get() );
      double fsample = static_cast<double>( this->m_FixedImage->GetPixel( movingImageIterator.GetIndex() ) );

      if( sample < movingImageMin )
        {
        movingImageMin = sample;
        }

      if( sample > movingImageMax )
        {
        movingImageMax = sample;
        }

      if( fsample < fixedImageMin )
        {
        fixedImageMin = fsample;
        }

      if( fsample > fixedImageMax )
        {
        fixedImageMax = fsample;
        }
      }
    }
/*
  float epsilon=0.01;
  fixedImageMax=1-epsilon;
  fixedImageMin=0.0+epsilon;
  movingImageMax=fixedImageMax;
  movingImageMin=fixedImageMin;
*/
//  std::cout <<"  fmx " << fixedImageMax << " fmn " << fixedImageMin << std::endl;
//  std::cout <<"  mmx " << movingImageMax << " mmn " << movingImageMin << std::endl;

  fixedImageMax = 1. * ( fixedImageMax - fixedImageMin ) + fixedImageMin;
  movingImageMax = 1. * ( movingImageMax - movingImageMin ) + movingImageMin;

  m_FixedImageBinSize = ( fixedImageMax - fixedImageMin )
    / static_cast<double>( m_NumberOfHistogramBins - 2 * this->m_Padding );
  m_FixedImageNormalizedMin = fixedImageMin / m_FixedImageBinSize
    - static_cast<double>( this->m_Padding );

  m_MovingImageBinSize = ( movingImageMax - movingImageMin )
    / static_cast<double>( m_NumberOfHistogramBins - 2 * this->m_Padding );
  m_MovingImageNormalizedMin = movingImageMin / m_MovingImageBinSize
    - static_cast<double>( this->m_Padding );

  //  m_FixedImageMarginalPDF.resize( m_NumberOfHistogramBins, 0.0 );
  //  m_MovingImageMarginalPDF.resize( m_NumberOfHistogramBins, 0.0 );
  m_FixedImageMarginalPDF = MarginalPDFType::New();
  m_MovingImageMarginalPDF = MarginalPDFType::New();
  typename MarginalPDFType::RegionType            mPDFRegion;
  typename MarginalPDFType::SizeType              mPDFSize;
  typename MarginalPDFType::IndexType              mPDFIndex;
  mPDFIndex.Fill( 0 );
  mPDFSize.Fill( m_NumberOfHistogramBins );
  mPDFRegion.SetIndex( mPDFIndex );
  mPDFRegion.SetSize( mPDFSize );
  m_FixedImageMarginalPDF->SetRegions( mPDFRegion );
  m_FixedImageMarginalPDF->Allocate();
  m_MovingImageMarginalPDF->SetRegions( mPDFRegion );
  m_MovingImageMarginalPDF->Allocate();

  /**
   * Allocate memory for the joint PDF and joint PDF derivatives.
   * The joint PDF and joint PDF derivatives are store as itk::Image.
   */
  m_JointPDF = JointPDFType::New();
  m_JointPDFDerivatives = JointPDFDerivativesType::New();

  // Instantiate a region, index, size
  JointPDFRegionType jointPDFRegion;
  JointPDFIndexType  jointPDFIndex;
  JointPDFSizeType   jointPDFSize;

  JointPDFDerivativesRegionType jointPDFDerivativesRegion;
  JointPDFDerivativesIndexType  jointPDFDerivativesIndex;
  JointPDFDerivativesSizeType   jointPDFDerivativesSize;
  // For the joint PDF define a region starting from {0,0}
  // with size {m_NumberOfHistogramBins, m_NumberOfHistogramBins}.
  // The dimension represents fixed image parzen window index
  // and moving image parzen window index, respectively.
  jointPDFIndex.Fill( 0 );
  jointPDFDerivativesIndex.Fill( 0 );
  jointPDFSize.Fill( m_NumberOfHistogramBins );
  jointPDFDerivativesSize.Fill(m_NumberOfHistogramBins );

  jointPDFRegion.SetIndex( jointPDFIndex );
  jointPDFRegion.SetSize( jointPDFSize );

  // Set the regions and allocate
  m_JointPDF->SetRegions( jointPDFRegion );
  m_JointPDF->Allocate();

  m_JointHist = JointPDFType::New();
  m_JointHist->SetRegions(jointPDFRegion );
  m_JointHist->Allocate();

  // For the derivatives of the joint PDF define a region starting from {0,0,0}
  // with size {m_NumberOfParameters,m_NumberOfHistogramBins,
  // m_NumberOfHistogramBins}. The dimension represents transform parameters,
  // fixed image parzen window index and moving image parzen window index,
  // respectively.
  jointPDFDerivativesIndex.Fill( 0 );
  jointPDFDerivativesSize[0] = m_NumberOfParameters;
  jointPDFDerivativesSize[1] = m_NumberOfHistogramBins;
  //  jointPDFDerivativesSize[2] = m_NumberOfHistogramBins;

  jointPDFDerivativesRegion.SetIndex( jointPDFDerivativesIndex );
  jointPDFDerivativesRegion.SetSize( jointPDFDerivativesSize );

  // std::cout << " 9 ";
  // Set the regions and allocate
  m_JointPDFDerivatives->SetRegions( jointPDFDerivativesRegion );
  //  m_JointPDFDerivatives->Allocate();

  m_InterpolatorIsBSpline = true;

  BSplineInterpolatorType * testPtr = dynamic_cast<BSplineInterpolatorType *>(
      this->m_Interpolator.GetPointer() );
  if( !testPtr )
    {
    m_InterpolatorIsBSpline = false;

    m_DerivativeCalculator = DerivativeFunctionType::New();
    m_DerivativeCalculator->SetInputImage( this->m_MovingImage );

    m_BSplineInterpolator = NULL;
    //    itkDebugMacro( "Interpolator is not BSpline" );
    }
  else
    {
    m_BSplineInterpolator = testPtr;
    m_DerivativeCalculator = NULL;
    // itkDebugMacro( "Interpolator is BSpline" );
    }

  m_TransformIsBSpline = false;
  m_BSplineTransform = NULL;

  m_NormalizeMetric = 1.0;
  for( int i = 0; i < ImageDimension; i++ )
    {
    m_NormalizeMetric *= this->m_FixedImage->GetLargestPossibleRegion().GetSize()[i];
    }

  this->GetProbabilities();
  this->ComputeMutualInformation();

  pdfinterpolator->SetInputImage(m_JointPDF);
  //  dpdfinterpolator->SetInputImage(m_JointPDFDerivatives);
  pdfinterpolator2->SetInputImage(m_FixedImageMarginalPDF);
  pdfinterpolator3->SetInputImage(m_MovingImageMarginalPDF);
  pdfinterpolator->SetSplineOrder(3);
  pdfinterpolator2->SetSplineOrder(3);
  pdfinterpolator3->SetSplineOrder(3);
}

/**
 * Get the both Value and Derivative Measure
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetProbabilities()
{
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion() );
  for( unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    MarginalPDFIndexType mind;
    mind[0] = j;
    m_FixedImageMarginalPDF->SetPixel(mind, 0);
    m_MovingImageMarginalPDF->SetPixel(mind, 0);
    }

  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer( 0.0 );
  m_JointHist->FillBuffer( 0.0 );
  //  m_JointPDFDerivatives->FillBuffer( 0.0 );

  unsigned long nSamples = 0;
  // unsigned long nFixedImageSamples=0;
  RandomIterator iter( this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion() );
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    bool takesample = true;
    if( this->m_FixedImageMask )
      {
      if( this->m_FixedImageMask->GetPixel( iter.GetIndex() ) < 1.e-6 )
        {
        takesample = false;
        }
      }

    if( takesample )
      {
      // Get sampled index
      FixedImageIndexType index = iter.GetIndex();
      double              fixedImageValue = iter.Get();
      CovariantVectorType fixedGradient;
      fixedGradient.Fill(0);
      typename FixedImageType::SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
      bool inimage = true;
      for( unsigned int dd = 0; dd < ImageDimension; dd++ )
        {
        if( index[dd] < 1 ||
            index[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1) )
          {
          inimage = false;
          }
        }

      if( inimage )
        {
        fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
        }

      MovingImagePointType mappedPoint;
      double               movingImageValue;
      movingImageValue = this->m_MovingImage->GetPixel( index );
      //      this->m_FixedImage->TransformIndexToPhysicalPoint(index,mappedPoint);

      // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
      double movingImageParzenWindowTerm =
        movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
      unsigned int movingImageParzenWindowIndex =
        static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

      // Make sure the extreme values are in valid bins
      if( movingImageParzenWindowIndex < this->m_Padding  )
        {
        movingImageParzenWindowIndex = this->m_Padding;
        }
      else if( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding  - 1 ) )
        {
        movingImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding  - 1;
        }

      // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
      double fixedImageParzenWindowTerm =
        fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;
      unsigned int fixedImageParzenWindowIndex =
        static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

      // Make sure the extreme values are in valid bins
      if( fixedImageParzenWindowIndex < this->m_Padding  )
        {
        fixedImageParzenWindowIndex = this->m_Padding;
        }
      else if( fixedImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding  - 1 ) )
        {
        fixedImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding  - 1;
        }

      /**
        * The region of support of the parzen window determines which bins
        * of the joint PDF are effected by the pair of image values.
        * Since we are using a cubic spline for the moving image parzen
        * window, four bins are affected.  The fixed image parzen window is
        * a zero-order spline (box car) and thus effects only one bin.
        *
        *  The PDF is arranged so that moving image bins corresponds to the
        * zero-th (column) dimension and the fixed image bins corresponds
        * to the first (row) dimension.
        *
        */
      // typename JointPDFType::IndexType pdfind;
      // pdfind[0]=fixedImageParzenWindowIndex;
      // pdfind[1]=movingImageParzenWindowIndex;

      //      m_JointPDF->SetPixel(pdfind,m_JointPDF->GetPixel(pdfind)+1);

      JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer()
        + ( fixedImageParzenWindowIndex * m_NumberOfHistogramBins );
      // Move the pointer to the first affected bin
      int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex );
      pdfPtr += pdfMovingIndex;
      *(pdfPtr) += static_cast<PDFValueType>( 1 );

      //	  static_cast<PDFValueType>(m_CubicBSplineKernel->Evaluate( fixedImageParzenWindowArg )));

      /*
      // Move the pointer to the first affected bin
      int pdfFixedIndex = static_cast<int>( fixedImageParzenWindowIndex ) - 1;


      for ( ; pdfFixedIndex <= static_cast<int>( fixedImageParzenWindowIndex ) + 2;
            pdfFixedIndex++ )
        {
    typename JointPDFType::IndexType pdfind;
    pdfind[1]=pdfFixedIndex;
    pdfind[0]=movingImageParzenWindowIndex;


        // Update PDF for the current intensity pair
    double fixedImageParzenWindowArg =
      static_cast<double>( pdfFixedIndex ) -
      static_cast<double>( fixedImageParzenWindowTerm );

    m_JointPDF->SetPixel(pdfind,m_JointPDF->GetPixel(pdfind)+
      static_cast<PDFValueType>(m_CubicBSplineKernel->Evaluate( fixedImageParzenWindowArg )));

    //	  m_JointPDFDerivatives->SetPixel(dpdfind,m_JointPDFDerivatives->GetPixel(dpdfind)+
    float bspd = static_cast<PDFValueType>(m_CubicBSplineDerivativeKernel->Evaluate( fixedImageParzenWindowArg ));


      // Compute PDF derivative contribution.
    //      this->ComputePDFDerivatives( pdfFixedIndex,movingImageParzenWindowIndex,fixedGradient, bspd );



        }  //end parzen windowing for loop



      // Move the pointer to the first affected bin
      pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex ) - 1;

      for ( ; pdfMovingIndex <= static_cast<int>( movingImageParzenWindowIndex ) + 2;
            pdfMovingIndex++ )
        {
    if (pdfMovingIndex == movingImageParzenWindowIndex) pdfMovingIndex++;

    // Update PDF for the current intensity pair
    double movingImageParzenWindowArg =
      static_cast<double>( pdfMovingIndex ) -
      static_cast<double>( movingImageParzenWindowTerm );


    typename JointPDFType::IndexType pdfind;
    pdfind[1]=fixedImageParzenWindowIndex;
    pdfind[0]=pdfMovingIndex;
    m_JointPDF->SetPixel(pdfind,m_JointPDF->GetPixel(pdfind)+
      static_cast<PDFValueType>(m_CubicBSplineKernel->Evaluate( movingImageParzenWindowArg ) ));

        }  //end parzen windowing for loop

      */
      ++nSamples;
      }
    }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator( m_JointPDF, m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;
  float  max = 0;
  float  max2 = 0;
  float  max3 = 0;
  float  max4 = 0;
  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
    m_JointHist->SetPixel(jointPDFIterator.GetIndex(), temp);
    if( temp > max )
      {
      max = temp;
      }
    if( temp > max2 && temp < max )
      {
      max2 = temp;
      }
    if( temp > max3 && temp < max2 )
      {
      max3 = temp;
      }
    if( temp > max4 && temp < max3 )
      {
      max4 = temp;
      }
    ++jointPDFIterator;
    }

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
// we do this to prevent background from over-whelming the computation
//      if (temp > max2) temp=(max2*0.9+max*0.1);
    jointPDFIterator.Set(temp);
    jointPDFSum += temp;
    ++jointPDFIterator;
    }

// of derivatives
  if( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  /*
// Normalize the joint PDF derivatives by the test image binsize and nSamples
typedef ImageRegionIterator<JointPDFDerivativesType> JointPDFDerivativesIteratorType;
JointPDFDerivativesIteratorType jointPDFDerivativesIterator (
m_JointPDFDerivatives, m_JointPDFDerivatives->GetBufferedRegion() );

jointPDFDerivativesIterator.GoToBegin();

double nFactor = 1.0 / ( m_MovingImageBinSize  * static_cast<double>( nSamples ) );

while( !jointPDFDerivativesIterator.IsAtEnd() )
{
jointPDFDerivativesIterator.Value() *= nFactor;
++jointPDFDerivativesIterator;
}

  */

  bool smoothjh = true;
  if( smoothjh )
    {
    typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
    typename dgtype::Pointer dg = dgtype::New();
    dg->SetInput(this->m_JointPDF);
    dg->SetVariance(1.);
    dg->SetUseImageSpacingOff();
    dg->SetMaximumError(.01f);
    dg->Update();
    this->m_JointPDF = dg->GetOutput();
    }

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter(
    m_JointPDF, m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    MarginalPDFIndexType mind;
    mind[0] = fixedIndex;
    m_FixedImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum) );
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;

    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    MarginalPDFIndexType mind;
    mind[0] = movingIndex;
    m_MovingImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum) );

    linearIter.NextLine();
    ++movingIndex;
    }
}

/**
 * Uniformly sample the fixed image domain using a random walk
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::NormalizeJointHist()
{
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator( m_JointPDF, m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;
  float  max = 0;
  float  max2 = 0;

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = m_JointHist->GetPixel(jointPDFIterator.GetIndex() );
    m_JointPDF->SetPixel(jointPDFIterator.GetIndex(), temp);
    // if (temp > max) max = temp;
    //      if (temp > max2 && temp < max) max2=temp;
    ++jointPDFIterator;
    }

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
    //      if (temp > max) temp=max;
    jointPDFIterator.Set(temp);
    jointPDFSum += temp;
    ++jointPDFIterator;
    }

  //  std::cout << " sum " << jointPDFSum << std::endl;
  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  bool smoothjh = true;
  if( smoothjh )
    {
    typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
    typename dgtype::Pointer dg = dgtype::New();
    dg->SetInput(this->m_JointPDF);
    dg->SetVariance(1.0);
    dg->SetUseImageSpacingOn();
    dg->SetMaximumError(.01f);
    dg->Update();
    this->m_JointPDF = dg->GetOutput();
    }

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter(
    m_JointPDF, m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    typename MarginalPDFType::IndexType mind;  mind[0] = fixedIndex;
    m_FixedImageMarginalPDF->GetPixel(mind) = static_cast<PDFValueType>(sum);
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;

    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    typename MarginalPDFType::IndexType mind;  mind[0] = movingIndex;
    m_MovingImageMarginalPDF->GetPixel(mind) = static_cast<PDFValueType>(sum);

    linearIter.NextLine();
    ++movingIndex;
    }
}

/**
 * Uniformly sample the fixed image domain using a random walk
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::SampleFixedImageDomain( FixedImageSpatialSampleContainer& samples )
{
  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<FixedImageType> RandomIterator2;
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion() );

  //  randIter.SetNumberOfSamples( m_NumberOfSpatialSamples );
  randIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();

  if( this->m_FixedImageMask )
    {
    }
  else
    {
    for( iter = samples.begin(); iter != end; ++iter )
      {
      // Get sampled index
      FixedImageIndexType index = randIter.GetIndex();
      // Get sampled fixed image value
      (*iter).FixedImageValue = randIter.Get();
      // Translate index to point
      (*iter).FixedImageIndex = index;
      this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                                                         (*iter).FixedImagePointValue );
      // Jump to random position
      ++randIter;
      //      ++randIter; // jump twice b/c it's not random
      }
    }
}

/**
 * Uniformly sample the fixed image domain using a random walk
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::SampleFixedImageDomainLocal( FixedImageSpatialSampleContainer& samples,
                               typename TFixedImage::IndexType oindex )
{
  typename FixedImageType::SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  typename FixedImageType::SizeType hradius = this->GetRadius();

  // this keeps some global samples
  unsigned int radiusval = 2; // (unsigned int)
  //  (0.35*pow((double)m_NumberOfSpatialSamples,(double)1.0/(double)ImageDimension)+0.5);
  //
  for( int i = 0; i < ImageDimension; i++ )
    {
    hradius[i] = radiusval;
    }

  NeighborhoodIterator<TDeformationField>
  asamIt( hradius,
          this->GetDeformationField(),
          this->GetDeformationField()->GetLargestPossibleRegion() );
  asamIt.SetLocation(oindex);
  unsigned int hoodlen = asamIt.Size();
  int          samplestep = 1;

  int indct = 0;
  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();
  iter = samples.begin();
  while( indct < m_NumberOfSpatialSamples && indct < hoodlen )
    {
    IndexType index = asamIt.GetIndex(indct);
    bool      inimage = true;
    for( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      if( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1) )
        {
        inimage = false;
        }
      }
    if( inimage )
      {
      (*iter).FixedImageValue = this->m_FixedImage->GetPixel(index);
      // Translate index to point
      (*iter).FixedImageIndex = index;
      this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                                                         (*iter).FixedImagePointValue );

      ++iter;
      }
    indct++;
    }
}

/**
 * Uniformly sample the fixed image domain using a random walk
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeFixedImageParzenWindowIndices( FixedImageSpatialSampleContainer& samples )
{
  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();
  for( iter = samples.begin(); iter != end; ++iter )
    {
    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).
    double windowTerm =
      static_cast<double>( (*iter).FixedImageValue ) / m_FixedImageBinSize
      - m_FixedImageNormalizedMin;
    unsigned int pindex = static_cast<unsigned int>( floor( windowTerm ) );

    // Make sure the extreme values are in valid bins
    if( pindex < this->m_Padding  )
      {
      pindex = this->m_Padding;
      }
    else if( pindex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      pindex = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    (*iter).FixedImageParzenWindowIndex = pindex;
    }
}

/**
 * Get the both Value and Derivative Measure
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative(IndexType oindex,
                        MeasureType& valuei,
                        DerivativeType& derivative1, DerivativeType& derivative2)
{
  double         value = 0;
  DerivativeType zero(ImageDimension);

  zero.Fill(0);

  double       fixedImageValue = (double)this->m_FixedImage->GetPixel(oindex);
  double       movingImageValue = this->GetMovingImageValue(oindex, zero);
  unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue);
//  unsigned int movingIndex =this->GetMovingValueIndex(movingImageValue);

  double dJPDF = 0, dFmPDF = 0, jointPDFValue = 0, fixedImagePDFValue = 0;

  // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
  double movingImageParzenWindowTerm =
    movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
  double fixedImageParzenWindowTerm =
    fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;
//      unsigned int fixedImageParzenWindowIndex =
//        static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

  if( movingImageParzenWindowTerm < this->m_Padding  )
    {
    movingImageParzenWindowTerm = this->m_Padding;
    }
  else if( movingImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1) )
    {
    movingImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
    }
  // Make sure the extreme values are in valid bins
  if( fixedImageParzenWindowTerm < this->m_Padding )
    {
    fixedImageParzenWindowTerm = this->m_Padding;
    }
  else if( fixedImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1) )
    {
    fixedImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
    }

  typename JointPDFType::PointType pdfind;
  pdfind[1] = fixedImageParzenWindowTerm;
  pdfind[0] = movingImageParzenWindowTerm;
  jointPDFValue = pdfinterpolator->Evaluate(pdfind);
  dJPDF = (1.0) * (pdfinterpolator->EvaluateDerivative( pdfind ) )[1];
  typename MarginalPDFType::PointType mind;
  mind[0] = fixedImageParzenWindowTerm;
  dFmPDF = (1.0) * (pdfinterpolator2->EvaluateDerivative( mind ) )[0];
  fixedImagePDFValue = pdfinterpolator2->Evaluate(mind);
  typename MarginalPDFType::IndexType mind2;
  mind2[0] = fixedIndex;

  double term1 = 0, term2 = 0, eps = 1.e-12;
  if( jointPDFValue > eps &&  (fixedImagePDFValue) > 0 )
    {
    term1 = dJPDF / jointPDFValue;
    term2 = dFmPDF / fixedImagePDFValue;
    value =  (term1 * (-1.0) + term2);
    //      if (fabs(value)>100)value=0;
    //      double vv=10,vv2=vv;
    //      if (value > vv) value=vv2;
    //      if (value < vv*(-1.0)) value=vv2*(-1.0);
    }    // end if-block to check non-zero bin contribution
  else
    {
    value = 0;
    }

  return value;
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivativeInv(IndexType oindex,
                           MeasureType& valuei,
                           DerivativeType& derivative1, DerivativeType& derivative2)
{
  double         value = 0;
  DerivativeType zero(ImageDimension);

  zero.Fill(0);

  double       fixedImageValue = (double)this->m_FixedImage->GetPixel(oindex);
  double       movingImageValue = (double)this->m_MovingImage->GetPixel(oindex);
  unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue);
  unsigned int movingIndex = this->GetMovingValueIndex(movingImageValue);

  double dJPDF = 0, dMmPDF = 0, jointPDFValue = 0, movingImagePDFValue = 0;

  // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
  double movingImageParzenWindowTerm =
    movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
  double fixedImageParzenWindowTerm =
    fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;

  if( movingImageParzenWindowTerm < this->m_Padding )
    {
    movingImageParzenWindowTerm = this->m_Padding;
    }
  else if( movingImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1 ) )
    {
    movingImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
    }

  // Make sure the extreme values are in valid bins
  if( fixedImageParzenWindowTerm < this->m_Padding )
    {
    fixedImageParzenWindowTerm = this->m_Padding;
    }
  else if( fixedImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1) )
    {
    fixedImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
    }

    {
    typename JointPDFType::IndexType pdfind2;
    pdfind2[1] = fixedIndex;
    pdfind2[0] = movingIndex;

    typename JointPDFType::PointType pdfind;
    pdfind[1] = fixedImageParzenWindowTerm;
    pdfind[0] = movingImageParzenWindowTerm;
    jointPDFValue = pdfinterpolator->Evaluate(pdfind);
    dJPDF = (1.0) * (pdfinterpolator->EvaluateDerivative( pdfind ) )[0];
    }

    {
    typename MarginalPDFType::PointType mind;
    mind[0] = movingImageParzenWindowTerm;
    movingImagePDFValue = pdfinterpolator3->Evaluate(mind);
    dMmPDF = (pdfinterpolator3->EvaluateDerivative(mind) )[0];
    }

  double term1 = 0, term2 = 0, eps = 1.e-12;
  if( jointPDFValue > eps &&  (movingImagePDFValue) > 0 )
    {
    term1 = dJPDF / jointPDFValue;
    term2 = dMmPDF / movingImagePDFValue;
    value =  (term1 * (-1.0) + term2);
//      if (fabs(value)>100)value=0;
    // double vv=10,vv2=vv;
    // if (value > vv) value=vv2;
    // if (value < vv*(-1.0)) value=vv2*(-1.0);
    } // end if-block to check non-zero bin contribution
  else
    {
    value = 0;
    }

  return value;
}

/**
 * Get the both Value and Derivative Measure
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative2(IndexType oindex,
                         MeasureType& valuei,
                         DerivativeType& derivative)
{
  double         value = 0;
  DerivativeType zero(ImageDimension);

  zero.Fill(0);

  double fixedImageValue = (double)this->m_FixedImage->GetPixel(oindex);
//  double movingImageValue0 = (double)this->m_MovingImage->GetPixel(oindex);
  double movingImageValue = this->GetMovingImageValue(oindex, derivative);
  // double fixedImageValue = (double)this->GetFixedImageValue(oindex,derivative);
  // double movingImageValue = this->m_MovingImage->GetPixel(oindex);
  //  double movingImageValueOld = this->GetMovingImageValue(oindex,zero);
  unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue);
  unsigned int movingIndex = this->GetMovingValueIndex(movingImageValue);
  // unsigned int movingIndex0 =this->GetMovingValueIndex(movingImageValue0);

  // Initialize sum to zero

  MarginalPDFIndexType mind;
  mind[0] = fixedIndex;
  double fixedImagePDFValue = m_FixedImageMarginalPDF->GetPixel(mind);
  mind[0] = movingIndex;
  double movingImagePDFValue = m_MovingImageMarginalPDF->GetPixel(mind);
  //  double movingImagePDFValueOld = m_MovingImageMarginalPDF->GetPixel(movingIndexOld];

  JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer()
    + ( fixedIndex * m_NumberOfHistogramBins);
  int pdfMovingIndex = static_cast<int>( movingIndex );
  pdfPtr += pdfMovingIndex;
  double jointPDFValue = *(pdfPtr);
  /*
  pdfPtr = m_JointPDF->GetBufferPointer() +
    ( fixedIndex* m_NumberOfHistogramBins);
  int pdfMovingIndex0 = static_cast<int>( movingIndex0 );
  pdfPtr += pdfMovingIndex0;
  double jointPDFValue0 = *(pdfPtr);
  */
  // check for non-zero bin contribution
//  typename JointPDFDerivativesType::IndexType jointPDFDerivIndex;
  double denom =  movingImagePDFValue * fixedImagePDFValue;
  double MIval = 0;
  if( jointPDFValue > 1e-16 &&  denom > 1e-16 )
    {
    // old
//      double entropynorm = 1.0 ;//
    //	movingImagePDFValue*log(movingImagePDFValue) + fixedImagePDFValue*log(fixedImagePDFValue);

    double pRatio = ( 1.0 + log( jointPDFValue / denom ) );

    MIval = pRatio;  // 1.0 + log( jointPDFValue / denom );
    }  // end if-block to check non-zero bin contribution

  value = static_cast<MeasureType>( MIval );
  /*
  double dJPDF = 0;
  {
    typename JointPDFType::PointType dpdfind;
    dpdfind[1]=fixedIndex;
    dpdfind[0]=movingIndex;
    dJPDF =(-1.)* (pdfinterpolator->EvaluateDerivative(dpdfind))[1];
    }*/

  return MIval;
}

/**
 * Get the match measure derivative
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetDerivative( const ParametersType& parameters, DerivativeType & derivative ) const
{
  MeasureType value;

  // call the combined version
  this->GetValueAndDerivative( parameters, value, derivative );
}

/**
 * Compute image derivatives using a central difference function
 * if we are not using a BSplineInterpolator, which includes
 * derivatives.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeImageDerivatives(
  const MovingImagePointType& mappedPoint,
  ImageDerivativesType& gradient ) const
{
  if( m_InterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_BSplineInterpolator->EvaluateDerivative( mappedPoint );
    }
  else
    {
    // For all generic interpolator use central differencing.
    gradient = m_DerivativeCalculator->Evaluate( mappedPoint );
    }
}

/**
 * Transform a point from FixedImage domain to MovingImage domain.
 * This function also checks if mapped point is within support region.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::TransformPoint(
  unsigned int sampleNumber,
  MovingImagePointType& mappedPoint,
  bool& sampleOk,
  double& movingImageValue ) const
{
  sampleOk = true;
  movingImageValue =
    this->m_MovingImage->GetPixel(  m_FixedImageSamples[sampleNumber].FixedImageIndex);
  this->m_FixedImage->TransformIndexToPhysicalPoint
    (m_FixedImageSamples[sampleNumber].FixedImageIndex, mappedPoint);
  return;
}

// Method to reinitialize the seed of the random number generator
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ReinitializeSeed()
{
  // This method should be the same used in the ImageRandomIterator
  //  vnl_sample_reseed();
}

// Method to reinitialize the seed of the random number generator
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ReinitializeSeed(int seed)
{
  // This method should be the same used in the ImageRandomIterator
  //  vnl_sample_reseed(seed);
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative3(IndexType oindex,
                         MeasureType& valuei,
                         DerivativeType& derivative)
{
  // Set output values to zero
  double value = NumericTraits<MeasureType>::Zero;

  derivative.Fill( NumericTraits<MeasureType>::Zero );

  double       fixedImageValue = (double)this->m_FixedImage->GetPixel(oindex);
  double       movingImageValue = (double)this->m_MovingImage->GetPixel(oindex);
  unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue);
  unsigned int movingIndex = this->GetMovingValueIndex(movingImageValue);

  MarginalPDFIndexType mind;
  mind[0] = fixedIndex;
  double fixedImagePDFValue = m_FixedImageMarginalPDF->GetPixel(mind);
  mind[0] = movingIndex;
  double movingImagePDFValue = m_MovingImageMarginalPDF->GetPixel(mind);

  JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer()
    + ( fixedIndex * m_NumberOfHistogramBins);
  int pdfMovingIndex = static_cast<int>( movingIndex );
  pdfPtr += pdfMovingIndex;
  double jointPDFValue = *(pdfPtr);

  double pRatio = 0;
  // check for non-zero bin contribution
  if( jointPDFValue > 1e-16 &&  fixedImagePDFValue > 1e-16 )
    {
    pRatio = log( jointPDFValue / fixedImagePDFValue );
    for( unsigned int mu = 0; mu < m_NumberOfParameters; mu++ )
      {
      typename JointPDFDerivativesType::IndexType dpdfind;
      dpdfind[0] = mu;
      dpdfind[1] = fixedIndex;
      dpdfind[2] = movingIndex;

      JointPDFDerivativesType::PixelType val = this->m_JointPDFDerivatives->GetPixel(dpdfind);
      derivative[mu] = val;
      }
    }  // end if-block to check non-zero bin contribution

  value = pRatio;

  return;
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputePDFDerivatives(
  int pdfFixedIndex,
  int pdfMovingIndex,
  CovariantVectorType& fixedImageGradientValue,
  double cubicBSplineDerivativeValue ) const
{
  for( unsigned int mu = 0; mu < m_NumberOfParameters; mu++ )
    {
    typename JointPDFDerivativesType::IndexType dpdfind;
    dpdfind[0] = mu;
    dpdfind[1] = pdfFixedIndex;
    dpdfind[2] = pdfMovingIndex;

    JointPDFDerivativesType::PixelType val = this->m_JointPDFDerivatives->GetPixel(dpdfind);
    val -= fixedImageGradientValue[mu] * cubicBSplineDerivativeValue;
    this->m_JointPDFDerivatives->SetPixel(dpdfind, val);
    }
}

/**
 * Get the both Value and Derivative Measure
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetProbabilitiesLocal(typename TFixedImage::IndexType centerindex, float radius)
{
  /*


  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();

  typedef ImageRegionConstIterator<MovingImageType> MovingIteratorType;
  MovingIteratorType movingImageIterator(
    this->m_MovingImage, this->m_MovingImage->GetBufferedRegion() );

  for ( movingImageIterator.GoToBegin();
        !movingImageIterator.IsAtEnd(); ++movingImageIterator)
    {

      IndexType index = movingImageIterator.GetIndex();
      float dist=0;
      for (unsigned int i=0; i<ImageDimension; i++) dist+=(index[i]-centerindex[i])*(index[i]-centerindex[i]);
      dist=sqrt(dist);

      if ( dist < radius)
  {

      double sample = static_cast<double>( movingImageIterator.Get() );
      double fsample = static_cast<double>( this->m_FixedImage->GetPixel( movingImageIterator.GetIndex() ));

      if ( sample < movingImageMin )
  {
    movingImageMin = sample;
  }

      if ( sample > movingImageMax )
  {
    movingImageMax = sample;
  }

      if ( fsample < fixedImageMin )
  {
    fixedImageMin = fsample;
  }

      if ( fsample > fixedImageMax )
  {
    fixedImageMax = fsample;
  }


  }
    }


  fixedImageMax=1.*( fixedImageMax - fixedImageMin )+fixedImageMin;
  movingImageMax=1.*( movingImageMax - movingImageMin )+movingImageMin;
  const int padding = 2;  // this will pad by 2 bins

  m_FixedImageBinSize = ( fixedImageMax - fixedImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_FixedImageNormalizedMin = fixedImageMin / m_FixedImageBinSize -
    static_cast<double>( padding );

  m_MovingImageBinSize = ( movingImageMax - movingImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_MovingImageNormalizedMin = movingImageMin / m_MovingImageBinSize -
    static_cast<double>( padding );

  */

  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion() );
  for( unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    MarginalPDFIndexType mind;
    mind[0] = j;
    m_FixedImageMarginalPDF->SetPixel(mind, 0);
    m_MovingImageMarginalPDF->SetPixel(mind, 0);
    }

  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer( 0.0 );
  m_JointHist->FillBuffer( 0.0 );
  //  m_JointPDFDerivatives->FillBuffer( 0.0 );

  unsigned long  nSamples = 0;
  unsigned long  nFixedImageSamples = 0;
  RandomIterator iter( this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion() );
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    // Get sampled index
    FixedImageIndexType index = iter.GetIndex();

    float dist = 0;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      dist += (index[i] - centerindex[i]) * (index[i] - centerindex[i]);
      }
    dist = sqrt(dist);

    if( dist < radius )
      {
      double              fixedImageValue = iter.Get();
      CovariantVectorType fixedGradient;
      fixedGradient.Fill(0);
      typename FixedImageType::SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
      bool inimage = true;
      for( unsigned int dd = 0; dd < ImageDimension; dd++ )
        {
        if( index[dd] < 1 ||
            index[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1) )
          {
          inimage = false;
          }
        }

      if( inimage )
        {
        fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
        }

      MovingImagePointType mappedPoint;
      double               movingImageValue;
      movingImageValue = this->m_MovingImage->GetPixel( index );

      double movingImageParzenWindowTerm =
        movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
      unsigned int movingImageParzenWindowIndex =
        static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

      // Make sure the extreme values are in valid bins
      if( movingImageParzenWindowIndex < this->m_Padding )
        {
        movingImageParzenWindowIndex = this->m_Padding;
        }
      else if( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
        {
        movingImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
        }

      // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
      double fixedImageParzenWindowTerm =
        fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;
      unsigned int fixedImageParzenWindowIndex =
        static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

      // Make sure the extreme values are in valid bins
      if( fixedImageParzenWindowIndex < this->m_Padding )
        {
        fixedImageParzenWindowIndex = this->m_Padding;
        }
      else if( fixedImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
        {
        fixedImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
        }

      typename JointPDFType::IndexType pdfind;

      JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer()
        + ( fixedImageParzenWindowIndex * m_NumberOfHistogramBins );
      // Move the pointer to the first affected bin
      int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex );
      pdfPtr += pdfMovingIndex;
      *(pdfPtr) += static_cast<PDFValueType>( 1 );

      //	  static_cast<PDFValueType>(m_CubicBSplineKernel->Evaluate( fixedImageParzenWindowArg )));

      ++nSamples;
      }
    }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator( m_JointPDF, m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
    //      std::cout <<" pdf ind " << jointPDFIterator.GetIndex() << " val " << temp << std::endl;
    // if (temp > max2) temp=max2;
    jointPDFIterator.Set(temp);
    jointPDFSum += temp;
    ++jointPDFIterator;
    }

  if( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  bool smoothjh = true;
  if( smoothjh )
    {
    typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
    typename dgtype::Pointer dg = dgtype::New();
    dg->SetInput(this->m_JointPDF);
    dg->SetVariance(1.0);
    dg->SetUseImageSpacingOff();
    dg->SetMaximumError(.01f);
    dg->Update();
    this->m_JointPDF = dg->GetOutput();
    }

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter(
    m_JointPDF, m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    MarginalPDFIndexType mind;
    mind[0] = fixedIndex;
    m_FixedImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum) );
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;

    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    MarginalPDFIndexType mind;
    mind[0] = movingIndex;
    m_MovingImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum) );

    linearIter.NextLine();
    ++movingIndex;
    }
}
} // end namespace itk

#endif
