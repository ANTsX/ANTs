/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWarpImageWAffineFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/10 15:39:12 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWarpImageWAffineFilter_txx
#define __itkWarpImageWAffineFilter_txx
#include "itkWarpImageWAffineFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace itk
{
/**
 * Default constructor.
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::WarpImageWAffineFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs( 2 );

  // Setup default values
  m_OutputSpacing.Fill( 1.0 );
  m_OutputOrigin.Fill( 0.0 );

  m_EdgePaddingValue = NumericTraits<PixelType>::Zero;

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_Interpolator = static_cast<InterpolatorType *>( interp.GetPointer() );

  m_SmoothScale = -1;

  m_TransformOrder = AffineFirst;
}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "EdgePaddingValue: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
}

/**
 * Set the output image spacing.
 *
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);

  this->SetOutputSpacing( s );
}

/**
 * Set the output image origin.
 *
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::SetOutputOrigin(
  const double* origin)
{
  PointType p(origin);

  this->SetOutputOrigin(p);
}

/**
 * Set deformation field as Inputs[1] for this ProcessObject.
 *
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::SetDeformationField(
  const DeformationFieldType * field )
{
  // const cast is needed because the pipeline is not const-correct.
  DeformationFieldType * input =
    const_cast<DeformationFieldType *>( field );

  this->ProcessObject::SetNthInput( 1, input );
}

/**
 * Return a pointer to the deformation field.
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
typename WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::DeformationFieldType
* WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::GetDeformationField(void)
  {
  return static_cast<DeformationFieldType *>
         ( this->ProcessObject::GetInput( 1 ) );
  }

/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::BeforeThreadedGenerateData()
{
  if( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }

  // Connect input image to interpolator
  // m_Interpolator->SetInputImage( this->GetInput() );

  if( m_CachedSmoothImage.IsNull() )
    {
    m_CachedSmoothImage = const_cast<InputImageType *>(this->GetInput() );
    }

  m_Interpolator->SetInputImage( m_CachedSmoothImage );
}

/**
 * Setup state of filter after multi-threading.
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage( NULL );
}

// /**
// * Compute the output for the region specified by outputRegionForThread.
// */
// template <class TInputImage,class TOutputImage,class TDeformationField, class TTransform>
// void
// WarpImageWAffineFilter<TInputImage,TOutputImage,TDeformationField, TTransform>
// ::ThreadedGenerateData(
//        const OutputImageRegionType& outputRegionForThread,
//        int threadId )
// {
//
//    InputImageConstPointer inputPtr = this->GetInput();
//    OutputImagePointer outputPtr = this->GetOutput();
//    DeformationFieldPointer fieldPtr = this->GetDeformationField();
//    TransformTypePointer aff = this->GetAffineTransform();
//
//    // std::cout << aff << std::endl;
//
//    // support progress methods/callbacks
//    ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
//
//    // iterator for the output image
//    ImageRegionIteratorWithIndex<OutputImageType> outputIt(
//            outputPtr, outputRegionForThread );
//
//    // iterator for the deformation field
//    ImageRegionIterator<DeformationFieldType> fieldIt(
//            fieldPtr, outputRegionForThread );
//
//    IndexType index;
//    PointType point1, point2, point3;
//    DisplacementType displacement;
//
//
//    int cnt = 0;
//    while( !outputIt.IsAtEnd() )
//    {
//        // get the output image index
//        index = outputIt.GetIndex();
//        outputPtr->TransformIndexToPhysicalPoint( index, point1 );
//
//        // get the required displacement
//        displacement = fieldIt.Get();
//
//        // compute the required input image point
//        for(unsigned int j = 0; j < ImageDimension; j++ )
//        {
//            point2[j] = point1[j] + displacement[j];
//        }
//
//        // songgang: followed by affine transform
//        // use affine transform
//        point3 = aff->TransformPoint(point2);
//
//        // get the interpolated value
//        if( m_Interpolator->IsInsideBuffer( point3 ) )
//        {
//            PixelType value = static_cast<PixelType>(
//                    m_Interpolator->Evaluate( point3 ) );
//            outputIt.Set( value );
//        }
//        else
//        {
//            outputIt.Set( m_EdgePaddingValue );
//        }
//        ++outputIt;
//        ++fieldIt;
//        progress.CompletedPixel();
//    }
//
// }

template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the input image
  InputImagePointer inputPtr = const_cast<InputImageType *>( this->GetInput() );

  if( inputPtr )
    {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  // just propagate up the output requested region for the
  // deformation field.
  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  //    OutputImagePointer outputPtr = this->GetOutput();
  if( fieldPtr )
    {
    //        fieldPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    fieldPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  return;
}

template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize( m_OutputSize );
  // outputLargestPossibleRegion.SetIndex( 0 );
  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
  outputPtr->SetSpacing( m_OutputSpacing );
  outputPtr->SetOrigin( m_OutputOrigin );

  //    DeformationFieldPointer fieldPtr = this->GetDeformationField();
  //    if( fieldPtr )
  //    {
  //        outputPtr->SetLargestPossibleRegion( fieldPtr->
  //                GetLargestPossibleRegion() );
  //    }
}

template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::SetSmoothScale(double scale)
{
  if( m_SmoothScale != scale )
    {
    // compute the new cached
    m_SmoothScale = scale;

    typename InputImageType::SpacingType inputSpacing = this->GetInput()->GetSpacing();
    typename InputImageType::RegionType::SizeType inputSize = this->GetInput()->GetLargestPossibleRegion().GetSize();

    typename InputImageType::SpacingType outputSpacing;
    typename InputImageType::RegionType::SizeType outputSize;

    double minimumSpacing = inputSpacing.GetVnlVector().min_value();
    double maximumSpacing = inputSpacing.GetVnlVector().max_value();

    InputImagePointer image = const_cast<InputImageType *>(this->GetInput() );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      double scaling = vnl_math_min( 1.0 / scale * minimumSpacing / inputSpacing[d],
                                     static_cast<double>( inputSize[d] ) / 32.0 );
      outputSpacing[d] = inputSpacing[d] * scaling;
      outputSize[d] = static_cast<unsigned long>( inputSpacing[d]
                                                  * static_cast<double>( inputSize[d] ) / outputSpacing[d] + 0.5 );

      double sigma = 0.2 * ( outputSpacing[d] / inputSpacing[d]  );

      typedef RecursiveGaussianImageFilter<InputImageType, InputImageType> GaussianFilterType;
      typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
      smoother->SetInputImage( image );
      smoother->SetDirection( d );
      smoother->SetNormalizeAcrossScale( false );

//            std::cout << "scale = " << scale << " => " << "sigma of dim " << d << ": " << sigma << " outsize " <<
// outputSize <<  std::endl;

      smoother->SetSigma( sigma );
      if( smoother->GetSigma() > 0.0 )
        {
        smoother->Update();
        image = smoother->GetOutput();
        }
      }

    m_CachedSmoothImage = image;

    SetOutputSpacing( outputSpacing );
    SetOutputOrigin( this->GetInput()->GetOrigin() );
    SetOutputSize(outputSize);
    }
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  int threadId )
{
  InputImageConstPointer  inputPtr = this->GetInput();
  OutputImagePointer      outputPtr = this->GetOutput();
  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  TransformTypePointer    aff = this->GetAffineTransform();

  // std::cout << aff << std::endl;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );

  // iterator for the output image
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(
    outputPtr, outputRegionForThread );

  // iterator for the deformation field
  ImageRegionIterator<DeformationFieldType> fieldIt(
    fieldPtr, outputRegionForThread );

  IndexType        index;
  PointType        point1, point2, point3;
  DisplacementType displacement;

  std::cout << "m_TransformOrder: " << m_TransformOrder << std::endl;

  int cnt = 0;
  while( !outputIt.IsAtEnd() )
    {
    // get the output image index
    index = outputIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( index, point1 );

    // get the required displacement
    displacement = fieldIt.Get();

    // compute the required input image point
    bool isinside = false;

    switch( m_TransformOrder )
      {
      case AffineFirst:
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          point2[j] = point1[j] + displacement[j];
          }
        point3 = aff->TransformPoint(point2);
        isinside = true;     // affine transform is always valid
        }
        break;
      case AffineLast:
        {
        point2 = aff->TransformPoint(point1);

        typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, float> DefaultInterpolatorType;
        typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
        vinterp->SetInputImage(fieldPtr);

        typename DefaultInterpolatorType::ContinuousIndexType  contind;
        // isinside = fieldPtr->TransformPhysicalPointToContinuousIndex(point2, contind);
        // explicitly written to avoid double / float type dismatching
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          contind[i] = ( (point2[i] - fieldPtr->GetOrigin()[i]) / fieldPtr->GetSpacing()[i] );
          }
        isinside = fieldPtr->GetLargestPossibleRegion().IsInside( contind );

        typename DefaultInterpolatorType::OutputType disp2;
        if( isinside )
          {
          disp2 = vinterp->EvaluateAtContinuousIndex( contind );
          }
        else
          {
          disp2.Fill(0);
          }
        for( int jj = 0; jj < ImageDimension; jj++ )
          {
          point3[jj] = disp2[jj] + point2[jj];
          }
        }
        break;
      default:
        itkExceptionMacro(<< "Affine order not set");
      }

    // get the interpolated value
    if( isinside && (m_Interpolator->IsInsideBuffer( point3 ) ) )
      {
      PixelType value = static_cast<PixelType>(
          m_Interpolator->Evaluate( point3 ) );
      outputIt.Set( value );
      }
    else
      {
      outputIt.Set( m_EdgePaddingValue );
      }
    ++outputIt;
    ++fieldIt;
    progress.CompletedPixel();
    }
}

template <class TInputImage, class TOutputImage, class TDeformationField, class TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDeformationField, TTransform>
::UpdateSizeByScale()
{
  DeformationFieldPointer field = this->GetDeformationField();

  SetOutputSpacing( field->GetSpacing() / m_SmoothScale );
  SetOutputOrigin( field->GetOrigin() );

  typename InputImageType::SizeType imgsz =  field->GetLargestPossibleRegion().GetSize();
  for( int ii = 0; ii < InputImageType::ImageDimension; ii++ )
    {
    imgsz[ii] = (typename InputImageType::SizeType::SizeValueType)(imgsz[ii] * m_SmoothScale + 0.5);
    }

  SetOutputSize(imgsz);
}
} // end namespace itk

#endif
