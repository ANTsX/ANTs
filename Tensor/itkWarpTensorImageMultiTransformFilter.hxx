/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWarpTensorImageMultiTransformFilter_hxx
#define __itkWarpTensorImageMultiTransformFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include <limits>
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "TensorFunctions.h"

namespace itk
{
template <typename TTensorType>
TTensorType
LogExpTensor(TTensorType inTensor, bool doLog)
{
  TTensorType result;

  if ((inTensor[0] != inTensor[0]) || (inTensor[1] != inTensor[1]) || (inTensor[2] != inTensor[2]) ||
      (inTensor[3] != inTensor[3]) || (inTensor[4] != inTensor[4]) || (inTensor[5] != inTensor[5]))
  {
    result[0] = result[1] = result[2] = result[3] = result[4] = result[5] = 0;
    // std::cout << " - " << inTensor << std::endl;
    return result;
  }

  double sigma = 1.0e-15;
  if ((inTensor[0] < sigma) || (inTensor[1] < sigma) || (inTensor[2] < sigma) || (inTensor[3] < sigma) ||
      (inTensor[4] < sigma) || (inTensor[5] < sigma))
  {
    result[0] = result[1] = result[2] = result[3] = result[4] = result[5] = 0;
    // std::cout << " - " << inTensor << std::endl;
    return result;
  }

  vnl_matrix<double> tensor(3, 3);
  tensor[0][0] = inTensor[0];
  tensor[1][0] = tensor[0][1] = inTensor[1];
  tensor[1][1] = inTensor[2];
  tensor[2][0] = tensor[0][2] = inTensor[3];
  tensor[2][1] = tensor[1][2] = inTensor[4];
  tensor[2][2] = inTensor[5];

  // std::cout << " + " << inTensor << std::endl;
  if ((itk::Math::FloatAlmostEqual(tensor[0][0], 0.0) && itk::Math::FloatAlmostEqual(tensor[0][1], 0.0) &&
       itk::Math::FloatAlmostEqual(tensor[0][2], 0.0) && itk::Math::FloatAlmostEqual(tensor[1][1], 0.0) &&
       itk::Math::FloatAlmostEqual(tensor[1][2], 0.0) && itk::Math::FloatAlmostEqual(tensor[2][2], 0.0)) ||
      (tensor.has_nans()) || (!tensor.is_finite()))
  {
    result[0] = result[1] = result[2] = result[3] = result[4] = result[5] = 0;
  }
  else
  {
    vnl_symmetric_eigensystem<double> eSystem(tensor);
    for (unsigned int i = 0; i < 3; i++)
    {
      if (doLog)
      {
        eSystem.D[i] = std::log(std::fabs(eSystem.D[i]));
      }
      else
      {
        eSystem.D[i] = std::exp(eSystem.D[i]);
      }
    }

    vnl_matrix<double> leTensor = eSystem.recompose();

    result[0] = leTensor[0][0];
    result[1] = leTensor[1][0];
    result[2] = leTensor[1][1];
    result[3] = leTensor[2][0];
    result[4] = leTensor[2][1];
    result[5] = leTensor[2][2];
  }

  return result;
}

/**
 * Default constructor.
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  WarpTensorImageMultiTransformFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs(1);

  // Setup default values
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);

  m_EdgePaddingValue = NumericTraits<PixelType>::ZeroValue();

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_Interpolator = static_cast<InterpolatorType *>(interp.GetPointer());

  m_SmoothScale = -1;

  // m_bOutputDisplacementField = false;

  // m_TransformOrder = AffineFirst;
}

// template <typename TInputImage,typename TOutputImage,typename TDisplacementField, typename TTransform>
// void
// WarpTensorImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::SetInterpolator1(InterpolatorPointer interp)
// {
//    m_Interpolator = static_cast<InterpolatorType*> (interp.GetPointer());
//   std::cout << "set interpolator in WarpImage:" << interp << std::endl;
// }

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PrintTransformList()
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetOutputSpacing(
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetOutputOrigin(
  const double * origin)
{
  PointType p(origin);

  this->SetOutputOrigin(p);
}

/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  BeforeThreadedGenerateData()
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage(nullptr);
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  GenerateInputRequestedRegion()
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  GenerateOutputInformation()
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
  // outputLargestPossibleRegion.SetIndex( 0 );
  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr->SetSpacing(this->m_OutputSpacing);
  outputPtr->SetOrigin(this->m_OutputOrigin);
  outputPtr->SetDirection(this->m_OutputDirection);

  // this->m_FullWarp->SetLargestPossibleRegion( outputLargestPossibleRegion );
  // this->m_FullWarp->SetSpacing( this->m_OutputSpacing );
  // this->m_FullWarp->SetOrigin( this->m_OutputOrigin );
  // this->m_FullWarp->SetDirection( this->m_OutputDirection  );
  // this->m_FullWarp->AllocateInitialized();

  // determine if the deformation field is the same dimension as the image
  // so that it does not need interpolation in the first step
  DetermineFirstDeformNoInterp();
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::SetSmoothScale(
  double /* scale */)
{
  /*
      if (m_SmoothScale != scale){
          // compute the new cached

  //       std::cout << "change smooth scale: " << m_SmoothScale << " ---> " << scale << std::endl;

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

  //           std::cout << "scale = " << scale << " => " << "sigma of dim " << d << ": " << sigma << " out size " <<
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  int                           threadId)
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  // std::cout << "inputPtr->GetOrigin():" << inputPtr->GetOrigin() << std::endl;
  // std::cout << "outputPtr->GetOrigin():" << outputPtr->GetOrigin() << std::endl;

  // std::exception();

  // Need to get full warp for reorientation of tensors

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
    IndexType index2 = outputIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(index2, point1);

    bool isinside = MultiTransformPoint(point1, point2, m_bFirstDeformNoInterp, index2);

    // std::cout << "point1:" << point1 << "  point2:" << point2 << " index:" << index2 << std::endl;
    // std::exception();
    // DisplacementFieldType::PixelType diff;
    // diff[0] = point2[0]-point1[0];
    // diff[1] = point2[1]-point1[1];
    // diff[2] = point2[2]-point1[2];

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
// WarpTensorImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::PushBackAffineTransform(
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
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  PushBackDisplacementFieldTransform(const DisplacementFieldType * t)
{
  if (t)
  {
    VarTransformType t1;
    t1.dex.field = const_cast<DisplacementFieldType *>(t);
    t1.dex.vinterp = DefaultVectorInterpolatorType::New();
    t1.dex.vinterp->SetInputImage(t1.dex.field);

    m_TransformList.push_back(SingleTransformItemType(EnumDisplacementFieldType, t1));
  }
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  MultiTransformSinglePoint(const PointType & point1, PointType & point2)
{
  IndexType null_index;

  bool isinside = MultiTransformPoint(point1, point2, false, null_index);

  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::IsOutOfNumericBoundary(
  const PointType & p)
{
  const DisplacementScalarValueType kMaxDisp = itk::NumericTraits<DisplacementScalarValueType>::max();

  bool b = false;

  for (int i = 0; i < ImageDimension; i++)
  {
    if (p[i] >= kMaxDisp)
    {
      b = true;
      break;
    }
  }

  return b;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  MultiInverseAffineOnlySinglePoint(const PointType & p1, PointType & p2)
{
  bool      isinside = true;
  PointType point1 = p1, point2;

  typename TransformListType::iterator it = m_TransformList.begin();
  for (; it != m_TransformList.end(); it++)
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
        // aff->GetInverse(aff);
        // std::cout << "aff after:" << aff << std::endl;
        // std::cout << "aff_inv after:" << aff_inv << std::endl;

        point2 = aff_inv->TransformPoint(point1);
        point1 = point2;
        isinside = true;
        break;
      }

      case EnumDisplacementFieldType:
      {
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

  p2 = point2;

  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
bool
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::MultiTransformPoint(
  const PointType & p1,
  PointType &       p2,
  bool              bFisrtDeformNoInterp,
  const IndexType & index)
{
  bool      isinside = false;
  PointType point1 = p1, point2;

  typename TransformListType::iterator it = m_TransformList.begin();
  for (; it != m_TransformList.end(); it++)
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
        if (bFisrtDeformNoInterp && it == m_TransformList.begin())
        {
          // use discrete coordinates
          DisplacementType displacement = fieldPtr->GetPixel(index);
          for (int j = 0; j < ImageDimension; j++)
          {
            point2[j] = point1[j] + displacement[j];
          }
          isinside = true;
        }
        else
        {
          // use continous coordinates
          typename DefaultVectorInterpolatorType::ContinuousIndexType contind;

          // use ITK implementation to use orientation header
          fieldPtr->TransformPhysicalPointToContinuousIndex(point1, contind);

          isinside = fieldPtr->GetLargestPossibleRegion().IsInside(contind);

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
          for (int jj = 0; jj < ImageDimension; jj++)
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

  p2 = point2;

  return isinside;
}

template <typename TInputImage, typename TOutputImage, typename TDisplacementField, typename TTransform>
void
WarpTensorImageMultiTransformFilter<TInputImage, TOutputImage, TDisplacementField, TTransform>::
  DetermineFirstDeformNoInterp()
{
  m_bFirstDeformNoInterp = false;
  if (m_TransformList.size() > 0)
  {
    if ((m_TransformList.front().first == EnumDisplacementFieldType))
    {
      DisplacementFieldPointer field = m_TransformList.front().second.dex.field;

      m_bFirstDeformNoInterp = (this->GetOutputSize() == field->GetLargestPossibleRegion().GetSize()) &
                               (this->GetOutputSpacing() == field->GetSpacing()) &
                               (this->GetOutputOrigin() == field->GetOrigin());

      //           std::cout << "in set: field size: " << field->GetLargestPossibleRegion().GetSize()
      //            << "output spacing: " << this->GetOutputSize() << std::endl;
      //           std::cout << field->GetSpacing() << " | " << this->GetOutputSpacing() << std::endl;
      //           std::cout << field->GetOrigin() << " | " << this->GetOutputOrigin() << std::endl;
    }
  }
}
} // end namespace itk

#endif // __itkWarpTensorImageMultiTransformFilter_hxx
