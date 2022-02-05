/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWarpImageWAffineFilter_hxx
#define __itkWarpImageWAffineFilter_hxx

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
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::WarpImageWAffineFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs(2);

  // Setup default values
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);

  m_EdgePaddingValue = NumericTraits<PixelType>::ZeroValue();

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_Interpolator = static_cast<InterpolatorType *>(interp.GetPointer());

  m_SmoothScale = -1;

  m_TransformOrder = AffineFirst;
}

/**
 * Standard PrintSelf method.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PrintSelf(std::ostream & os,
                                                                                             Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  ;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "EdgePaddingValue: " << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
}

/**
 * Set the output image spacing.
 *
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetOutputSpacing(
  const double * spacing)
{
  SpacingType s(spacing);

  this->SetOutputSpacing(s);
}

/**
 * Set the output image origin.
 *
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetOutputOrigin(
  const double * origin)
{
  PointType p(origin);

  this->SetOutputOrigin(p);
}

/**
 * Set deformation field as Inputs[1] for this ProcessObject.
 *
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetDisplacementField(
  const DisplacementFieldType * field)
{
  // const cast is needed because the pipeline is not const-correct.
  DisplacementFieldType * input = const_cast<DisplacementFieldType *>(field);

  this->ProcessObject::SetNthInput(1, input);
}

/**
 * Return a pointer to the deformation field.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
typename WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::DisplacementFieldType *
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::GetDisplacementField(void)
{
  return static_cast<DisplacementFieldType *>(this->ProcessObject::GetInput(1));
}

/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::BeforeThreadedGenerateData()
{
  if (!m_Interpolator)
  {
    itkExceptionMacro(<< "Interpolator not set");
  }

  // Connect input image to interpolator
  // m_Interpolator->SetInputImage( this->GetInput() );

  if (m_CachedSmoothImage.IsNull())
  {
    m_CachedSmoothImage = const_cast<InputImageType *>(this->GetInput());
  }

  m_Interpolator->SetInputImage(m_CachedSmoothImage);
}

/**
 * Setup state of filter after multi-threading.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage(nullptr);
}

// /**
// * Compute the output for the region specified by outputRegionForThread.
// */
// template <typename TInputImage,typename TOutputImage,typename TDisplacementField, typename TTransform>
// void
// WarpImageWAffineFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::ThreadedGenerateData(
//        const OutputImageRegionType& outputRegionForThread,
//        int threadId )
// {
//
//    InputImageConstPointer inputPtr = this->GetInput();
//    OutputImagePointer outputPtr = this->GetOutput();
//    DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
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
//    ImageRegionIterator<DisplacementFieldType> fieldIt(
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

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::GenerateInputRequestedRegion()
{
  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the input image
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());

  if (inputPtr)
  {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
  }

  // just propagate up the output requested region for the
  // deformation field.
  DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
  //    OutputImagePointer outputPtr = this->GetOutput();
  if (fieldPtr)
  {
    //        fieldPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    fieldPtr->SetRequestedRegionToLargestPossibleRegion();
  }

  return;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer outputPtr = this->GetOutput();

  if (!outputPtr)
  {
    return;
  }

  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(m_OutputSize);
  // outputLargestPossibleRegion.SetIndex( 0 );
  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr->SetSpacing(m_OutputSpacing);
  outputPtr->SetOrigin(m_OutputOrigin);

  //    DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
  //    if( fieldPtr )
  //    {
  //        outputPtr->SetLargestPossibleRegion( fieldPtr->
  //                GetLargestPossibleRegion() );
  //    }
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetSmoothScale(double scale)
{
  if (!itk::Math::FloatAlmostEqual(m_SmoothScale, scale))
  {
    // compute the new cached
    m_SmoothScale = scale;

    typename InputImageType::SpacingType          inputSpacing = this->GetInput()->GetSpacing();
    typename InputImageType::RegionType::SizeType inputSize = this->GetInput()->GetLargestPossibleRegion().GetSize();

    typename InputImageType::SpacingType          outputSpacing;
    typename InputImageType::RegionType::SizeType outputSize;

    double minimumSpacing = inputSpacing.GetVnlVector().min_value();

    InputImagePointer image = const_cast<InputImageType *>(this->GetInput());
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      double scaling =
        std::min(1.0 / scale * minimumSpacing / inputSpacing[d], static_cast<double>(inputSize[d]) / 32.0);
      outputSpacing[d] = inputSpacing[d] * scaling;
      outputSize[d] =
        static_cast<unsigned long>(inputSpacing[d] * static_cast<double>(inputSize[d]) / outputSpacing[d] + 0.5);

      double sigma = 0.2 * (outputSpacing[d] / inputSpacing[d]);

      typedef RecursiveGaussianImageFilter<InputImageType, InputImageType> GaussianFilterType;
      typename GaussianFilterType::Pointer                                 smoother = GaussianFilterType::New();
      smoother->SetInputImage(image);
      smoother->SetDirection(d);
      smoother->SetNormalizeAcrossScale(false);

      //            std::cout << "scale = " << scale << " => " << "sigma of dim " << d << ": " << sigma << " outsize "
      // <<
      // outputSize <<  std::endl;

      smoother->SetSigma(sigma);
      if (smoother->GetSigma() > 0.0)
      {
        smoother->Update();
        image = smoother->GetOutput();
      }
    }

    m_CachedSmoothImage = image;

    SetOutputSpacing(outputSpacing);
    SetOutputOrigin(this->GetInput()->GetOrigin());
    SetOutputSize(outputSize);
  }
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  InputImageConstPointer   inputPtr = this->GetInput();
  OutputImagePointer       outputPtr = this->GetOutput();
  DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
  TransformTypePointer     aff = this->GetAffineTransform();

  // std::cout << aff << std::endl;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // iterator for the output image
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread);

  // iterator for the deformation field
  ImageRegionIterator<DisplacementFieldType> fieldIt(fieldPtr, outputRegionForThread);

  IndexType        index;
  PointType        point1, point2, point3;
  DisplacementType displacement;

  std::cout << "m_TransformOrder: " << m_TransformOrder << std::endl;

  while (!outputIt.IsAtEnd())
  {
    // get the output image index
    index = outputIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(index, point1);

    // get the required displacement
    displacement = fieldIt.Get();

    // compute the required input image point
    bool isinside = false;

    switch (m_TransformOrder)
    {
      case AffineFirst:
      {
        for (unsigned int j = 0; j < ImageDimension; j++)
        {
          point2[j] = point1[j] + displacement[j];
        }
        point3 = aff->TransformPoint(point2);
        isinside = true; // affine transform is always valid
      }
      break;
      case AffineLast:
      {
        point2 = aff->TransformPoint(point1);

        typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, float> VecLinInterpolatorType;
        typename VecLinInterpolatorType::Pointer vinterp = DefaultInterpolatorType::New();
        vinterp->SetInputImage(fieldPtr);

        typename DefaultInterpolatorType::ContinuousIndexType contind;
        // isinside = fieldPtr->TransformPhysicalPointToContinuousIndex(point2, contind);
        // explicitly written to avoid double / float type dismatching
        for (unsigned int i = 0; i < ImageDimension; i++)
        {
          contind[i] = ((point2[i] - fieldPtr->GetOrigin()[i]) / fieldPtr->GetSpacing()[i]);
        }
        isinside = fieldPtr->GetLargestPossibleRegion().IsInside(contind);

        typename DefaultInterpolatorType::OutputType disp2;
        if (isinside)
        {
          disp2 = vinterp->EvaluateAtContinuousIndex(contind);
        }
        else
        {
          disp2.Fill(0);
        }
        for (int jj = 0; jj < ImageDimension; jj++)
        {
          point3[jj] = disp2[jj] + point2[jj];
        }
      }
      break;
      default:
        itkExceptionMacro(<< "Affine order not set");
    }

    // get the interpolated value
    if (isinside && (m_Interpolator->IsInsideBuffer(point3)))
    {
      PixelType value = static_cast<PixelType>(m_Interpolator->Evaluate(point3));
      outputIt.Set(value);
    }
    else
    {
      outputIt.Set(m_EdgePaddingValue);
    }
    ++outputIt;
    ++fieldIt;
    progress.CompletedPixel();
  }
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageWAffineFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::UpdateSizeByScale()
{
  DisplacementFieldPointer field = this->GetDisplacementField();

  SetOutputSpacing(field->GetSpacing() / m_SmoothScale);
  SetOutputOrigin(field->GetOrigin());

  typename InputImageType::SizeType imgsz = field->GetLargestPossibleRegion().GetSize();
  for (int ii = 0; ii < InputImageType::ImageDimension; ii++)
  {
    imgsz[ii] = (typename InputImageType::SizeType::SizeValueType)(imgsz[ii] * m_SmoothScale + 0.5);
  }

  SetOutputSize(imgsz);
}
} // end namespace itk

#endif
