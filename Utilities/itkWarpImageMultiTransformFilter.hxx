/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWarpImageMultiTransformFilter_hxx
#define __itkWarpImageMultiTransformFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include <limits>

namespace itk
{
/**
 * Default constructor.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  WarpImageMultiTransformFilter()
{
  // Setup default values
  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_OutputDirection.SetIdentity();

  m_OutputSize.Fill(0);
  m_OutputStartIndex.Fill(0);

  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs(1);

  //  PixelType zeropix = NumericTraits<PixelType>::ZeroValue();
  //  m_EdgePaddingValue = zeropix;

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_Interpolator = static_cast<InterpolatorType *>(interp.GetPointer());

  m_SmoothScale = -1;

  this->DynamicMultiThreadingOff();
  // m_bOutputDisplacementField = false;

  // m_TransformOrder = AffineFirst;
}

// template <typename TInputImage,typename TOutputImage,typename TDisplacementField, typename TTransform>
// void
// WarpImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::SetInterpolator1(InterpolatorPointer interp)
// {
//    m_Interpolator = static_cast<InterpolatorType*> (interp.GetPointer());
//    std::cout << "set interpolator in WarpImage:" << interp << std::endl;
// }

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PrintTransformList()
{
  std::cout << "transform list: " << std::endl;

  typename TransformListType::iterator it = (m_TransformList.begin());
  for (int ii = 0; it != m_TransformList.end(); it++, ii++)
  {
    switch (it->first)
    {
      case EnumAffineType:
      {
        std::cout << '[' << ii << "]: EnumAffineType" << std::endl;
        std::cout << it->second.aex.aff << std::endl;
      }
      break;
      case EnumDisplacementFieldType:
        std::cout << '[' << ii << "]: EnumDisplacementFieldType: size"
                  << it->second.dex.field->GetLargestPossibleRegion().GetSize() << std::endl;
    }
  }
}

/**
 * Standard PrintSelf method.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PrintSelf(std::ostream & os,
                                                                                                    Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  ;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "EdgePaddingValue: " << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;

  os << indent << "m_bFirstDeformNoInterp = " << m_bFirstDeformNoInterp << std::endl;
}

/**
 * Set the output image spacing.
 *
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetOutputParametersFromImage(
  const ImageBaseType * image)
{
  this->SetOutputOrigin(image->GetOrigin());
  this->SetOutputSpacing(image->GetSpacing());
  this->SetOutputDirection(image->GetDirection());
  this->SetOutputStartIndex(image->GetLargestPossibleRegion().GetIndex());
  this->SetOutputSize(image->GetLargestPossibleRegion().GetSize());
}

/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::BeforeThreadedGenerateData()
{
  if (!m_Interpolator)
  {
    itkExceptionMacro(<< "Interpolator not set");
  }

  // Connect input image to interpolator
  // m_Interpolator->SetInputImage( this->GetInput() );

  if (m_CachedSmoothImage.IsNull() && (this->GetInput()))
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
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage(nullptr);
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::GenerateInputRequestedRegion()
{
  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the input image
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());

  if (inputPtr)
  {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
  }

  return;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer outputPtr = this->GetOutput();

  if (!outputPtr)
  {
    return;
  }

  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(this->m_OutputSize);
  outputLargestPossibleRegion.SetIndex(this->m_OutputStartIndex);
  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr->SetSpacing(this->m_OutputSpacing);
  outputPtr->SetOrigin(this->m_OutputOrigin);
  outputPtr->SetDirection(this->m_OutputDirection);
  // determine if the deformation field is the same dimension as the image
  // so that it does not need interpolation in the first step
  DetermineFirstDeformNoInterp();
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetSmoothScale(
  double /* scale */)
{
  /*
      if (m_SmoothScale != scale){
          // compute the new cached

  //        std::cout << "change smooth scale: " << m_SmoothScale << " ---> " << scale << std::endl;

          m_SmoothScale = scale;



          typename InputImageType::SpacingType inputSpacing = this->GetInput()->GetSpacing();
          typename InputImageType::RegionType::SizeType inputSize =
  this->GetInput()->GetLargestPossibleRegion().GetSize();

          typename InputImageType::SpacingType outputSpacing;
          typename InputImageType::RegionType::SizeType outputSize;


          double minimumSpacing = inputSpacing.GetVnlVector().min_value();
          double maximumSpacing = inputSpacing.GetVnlVector().max_value();


          InputImagePointer image = const_cast<InputImageType *> (this->GetInput());
          for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
              double scaling = std::min( 1.0 / scale * minimumSpacing / inputSpacing[d],
                      static_cast<double>( inputSize[d] ) / 32.0 );
              outputSpacing[d] = inputSpacing[d] * scaling;
              outputSize[d] = static_cast<unsigned long>( inputSpacing[d] *
                      static_cast<double>( inputSize[d] ) / outputSpacing[d] + 0.5 );

              double sigma = 0.25  * ( outputSpacing[d] / inputSpacing[d]  );
              if (sigma < 0) sigma=0;

              typedef RecursiveGaussianImageFilter<InputImageType, InputImageType> GaussianFilterType;
              typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
              smoother->SetInputImage( image );
              smoother->SetDirection( d );
              smoother->SetNormalizeAcrossScale( false );

  //            std::cout << "scale = " << scale << " => " << "sigma of dim " << d << ": " << sigma << " out size " <<
  outputSize <<  " spc1 " << outputSpacing << " in " << inputSpacing << std::endl;

              smoother->SetSigma( sigma );
              if ( smoother->GetSigma() > 0.0 )
              {
                  smoother->Update();
                  image = smoother->GetOutput();
              }
          }



          SetOutputSpacing( outputSpacing );
          SetOutputOrigin( this->GetInput()->GetOrigin() );
          SetOutputSize(outputSize);
      }
  */

  InputImagePointer image = const_cast<InputImageType *>(this->GetInput());

  m_CachedSmoothImage = image;
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  // std::cout << "inputPtr->GetOrigin():" << inputPtr->GetOrigin() << std::endl;
  // std::cout << "outputPtr->GetOrigin():" << outputPtr->GetOrigin() << std::endl;

  // std::exception();

  IndexType index;

  index.Fill(0);
  this->m_EdgePaddingValue = inputPtr->GetPixel(index);

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // iterator for the output image
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread);

  while (!outputIt.IsAtEnd())
  {
    PointType point1, point2;

    // get the output image index
    IndexType _index = outputIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(_index, point1);

    bool isinside = MultiTransformPoint(point1, point2, m_bFirstDeformNoInterp, _index);

    // std::cout << "point1:" << point1 << "  point2:" << point2 << " index:" << index << std::endl;
    // std::exception();

    // warp the image
    // get the interpolated value
    if (isinside && (m_Interpolator->IsInsideBuffer(point2)))
    {
      PixelType value = static_cast<PixelType>(m_Interpolator->Evaluate(point2));
      outputIt.Set(value);
    }
    else
    {
      // std::cout << "OUTSIDE" << " isinside:" << isinside << " m_Interpolator->IsInsideBuffer( point2 ):" <<
      // m_Interpolator->IsInsideBuffer( point2 ) <<  std::endl;
      outputIt.Set(m_EdgePaddingValue);
    }

    ++outputIt;
  }

  progress.CompletedPixel();
}

// template <typename TInputImage,typename TOutputImage,typename TDisplacementField, typename TTransform>
// void
// WarpImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::UpdateSizeByScale()
// {
//
//    //    DisplacementFieldPointer field = this->GetDisplacementField();
//    //
//    //    SetOutputSpacing( field->GetSpacing() / m_SmoothScale );
//    //    SetOutputOrigin( field->GetOrigin() );
//    //
//    //    typename InputImageType::SizeType imgsz =  field->GetLargestPossibleRegion().GetSize();
//    //    for(int ii = 0; ii < InputImageType::ImageDimension; ii++) imgsz[ii] = (typename
// InputImageType::SizeType::SizeValueType) (imgsz[ii] * m_SmoothScale + 0.5);
//    //
//    //    SetOutputSize(imgsz);
//
// }

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PushBackAffineTransform(
  const TransformType * t)
{
  if (t)
  {
    VarTransformType t1;
    t1.aex.aff = const_cast<TransformType *>(t);
    m_TransformList.push_back(SingleTransformItemType(EnumAffineType, t1));
  }
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  PushBackDisplacementFieldTransform(const DisplacementFieldType * t)
{
  if (t)
  {
    VarTransformType t1;
    t1.dex.field = const_cast<DisplacementFieldType *>(t);
    t1.dex.vinterp = DefaultVectorInterpolatorType::New();
    t1.dex.vinterp->SetInputImage(t1.dex.field);
    //    t1.dex.vinterp->SetParameters(NULL,1);
    m_TransformList.push_back(SingleTransformItemType(EnumDisplacementFieldType, t1));
  }
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::MultiTransformSinglePoint(
  const PointType & point1,
  PointType &       point2)
{
  IndexType null_index;
  null_index.Fill(0);

  bool isinside = MultiTransformPoint(point1, point2, false, null_index);

  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::IsOutOfNumericBoundary(
  const PointType & p)
{
  const DisplacementScalarValueType kMaxDisp = itk::NumericTraits<DisplacementScalarValueType>::max();

  bool b = false;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (p[i] >= static_cast<double>(kMaxDisp))
    {
      b = true;
      break;
    }
  }

  return b;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::ComposeAffineOnlySequence(
  const PointType &      center_output,
  TransformTypePointer & affine_output)
{
  affine_output->SetIdentity();
  affine_output->SetCenter(center_output);
  for (typename TransformListType::iterator it = m_TransformList.begin(); it != m_TransformList.end(); it++)
  {
    SingleTransformType ttype = it->first;

    switch (ttype)
    {
      case EnumAffineType:
      {
        TransformTypePointer aff = it->second.aex.aff;
        affine_output->Compose(aff, 0);
        break;
      }
      case EnumDisplacementFieldType:
      {
        std::cout << " Ignore deformation type ... " << std::endl;
        itkExceptionMacro(<< "Affine Only Sequence must only contain Affine Transforms, DisplacementField Found!");
      }
      break;
      default:
        itkExceptionMacro(<< "Single Transform Not Supported!");
    }
  }
  return;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  MultiInverseAffineOnlySinglePoint(const PointType & p1, PointType & point2)
{
  bool      isinside = true;
  PointType point1 = p1;

  for (typename TransformListType::iterator it = m_TransformList.begin(); it != m_TransformList.end(); it++)
  {
    SingleTransformType ttype = it->first;

    switch (ttype)
    {
      case EnumAffineType:
      {
        TransformTypePointer aff = it->second.aex.aff;
        TransformTypePointer aff_inv = TransformTypePointer::ObjectType::New();
        // std::cout << "aff before:" << aff << std::endl;
        aff->GetInverse(aff_inv);
        // std::cout << "aff after:" << aff << std::endl;
        // std::cout << "aff_inv after:" << aff_inv << std::endl;

        point2 = aff_inv->TransformPoint(point1);
        point1 = point2;
        isinside = true;
        break;
      }
      case EnumDisplacementFieldType:
      {
        itkExceptionMacro(<< "Affine Only Sequence must only contain Affine Transforms, DisplacementField Found!");
      }
      break;
      default:
        itkExceptionMacro(<< "Single Transform Not Supported!");
    }
    if (IsOutOfNumericBoundary(point2))
    {
      isinside = false;
      break;
    }
    point1 = point2;
  }
  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::MultiTransformPoint(
  const PointType & p1,
  PointType &       point2,
  bool              bFirstDeformNoInterp,
  const IndexType & index)
{
  bool      isinside = false;
  PointType point1 = p1;

  for (typename TransformListType::iterator it = m_TransformList.begin(); it != m_TransformList.end(); it++)
  {
    SingleTransformType ttype = it->first;

    switch (ttype)
    {
      case EnumAffineType:
      {
        TransformTypePointer aff = it->second.aex.aff;
        point2 = aff->TransformPoint(point1);
        point1 = point2;
        isinside = true;
      }
      break;
      case EnumDisplacementFieldType:
      {
        DisplacementFieldPointer fieldPtr = it->second.dex.field;
        if (bFirstDeformNoInterp && it == m_TransformList.begin())
        {
          // use discrete coordinates
          DisplacementType displacement = fieldPtr->GetPixel(index);
          for (unsigned int j = 0; j < ImageDimension; j++)
          {
            point2[j] = point1[j] + static_cast<double>(displacement[j]);
          }
          isinside = true;
        }
        else
        {
          // use continous coordinates
          typename DefaultVectorInterpolatorType::ContinuousIndexType contind;
          // use ITK implementation to use orientation header
          isinside = fieldPtr->TransformPhysicalPointToContinuousIndex(point1, contind);

          VectorInterpolatorPointer                          vinterp = it->second.dex.vinterp;
          typename DefaultVectorInterpolatorType::OutputType disp2;
          if (isinside)
          {
            disp2 = vinterp->EvaluateAtContinuousIndex(contind);
          }
          else
          {
            disp2.Fill(0);
          }
          for (unsigned int jj = 0; jj < ImageDimension; jj++)
          {
            point2[jj] = disp2[jj] + point1[jj];
          }
        }
      }
      break;
      default:
        itkExceptionMacro(<< "Single Transform Not Supported!");
    }

    if (IsOutOfNumericBoundary(point2))
    {
      isinside = false;
      break;
    }
    point1 = point2;
  }
  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::DetermineFirstDeformNoInterp()
{
  m_bFirstDeformNoInterp = false;
  if (m_TransformList.size() > 0)
  {
    if (m_TransformList.front().first == EnumDisplacementFieldType)
    {
      DisplacementFieldPointer field = m_TransformList.front().second.dex.field;

      m_bFirstDeformNoInterp = (this->GetOutputSize() == field->GetLargestPossibleRegion().GetSize()) &&
                               (this->GetOutputSpacing() == field->GetSpacing()) &&
                               (this->GetOutputOrigin() == field->GetOrigin());

      //            std::cout << "in set: field size: " << field->GetLargestPossibleRegion().GetSize()
      //            << "output spacing: " << this->GetOutputSize() << std::endl;
      //            std::cout << field->GetSpacing() << " | " << this->GetOutputSpacing() << std::endl;
      //            std::cout << field->GetOrigin() << " | " << this->GetOutputOrigin() << std::endl;
    }
  }
}
} // end namespace itk

#endif // __itkWarpImageMultiTransformFilter_hxx
