#ifndef __itkGeneralToBSplineDisplacementFieldFilter_hxx
#define __itkGeneralToBSplineDisplacementFieldFilter_hxx


#include "itkBSplineControlPointImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
GeneralToBSplineDisplacementFieldFilter<TInputImage, TOutputImage>::GeneralToBSplineDisplacementFieldFilter()
{
  this->m_IgnorePixelValue.Fill(NumericTraits<InputPixelComponentType>::max());
  this->m_SplineOrder = 3;
  this->m_NumberOfControlPoints.Fill(this->m_SplineOrder + 1);

  //  this->m_ConfidenceImage = nullptr;
}

template <typename TInputImage, typename TOutputImage>
GeneralToBSplineDisplacementFieldFilter<TInputImage, TOutputImage>::~GeneralToBSplineDisplacementFieldFilter()
{}

template <typename TInputImage, typename TOutputImage>
void
GeneralToBSplineDisplacementFieldFilter<TInputImage, TOutputImage>::GenerateData()
{
  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  //  typename BSplineFilterType::WeightsContainerType confidenceValues;
  //  confidenceValues->Initialize();

  ImageRegionConstIteratorWithIndex<InputImageType> It(this->GetInput(), this->GetInput()->GetRequestedRegion());

  itkDebugMacro(<< "Extracting points from input deformation field. ");

  unsigned int N = 0;
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    InputPixelType data = It.Get();

    if (data != this->m_IgnorePixelValue)
    {
      typename PointSetType::PointType point;
      this->GetInput()->TransformIndexToPhysicalPoint(It.GetIndex(), point);

      fieldPoints->SetPointData(N, data);
      fieldPoints->SetPoint(N, point);

      //      if ( this->m_ConfidenceImage )
      //        {
      //        confidenceValues->InsertElement
      //          ( N, this->m_ConfidenceImage->GetPixel( It.GetIndex() ) );
      //        }
      N++;
    }
  }

  itkDebugMacro("Calculating the B-spline deformation field. ");

  typename OutputImageType::PointType   origin;
  typename OutputImageType::SpacingType spacing;
  typename OutputImageType::SizeType    size;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    origin[i] = this->GetInput(0)->GetOrigin()[i];
    spacing[i] = this->GetInput(0)->GetSpacing()[i];
    size[i] = this->GetInput(0)->GetRequestedRegion().GetSize()[i];
  }

  typename BSplineFilterType::ArrayType close;
  close.Fill(false);

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetOrigin(origin);
  bspliner->SetSpacing(spacing);
  bspliner->SetSize(size);
  bspliner->SetNumberOfLevels(this->m_NumberOfLevels);
  bspliner->SetSplineOrder(this->m_SplineOrder);
  bspliner->SetNumberOfControlPoints(this->m_NumberOfControlPoints);
  bspliner->SetCloseDimension(close);
  bspliner->SetInput(fieldPoints);
  bspliner->SetGenerateOutputImage(true);
  bspliner->Update();
  this->SetNthOutput(0, bspliner->GetOutput());
}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputImage, typename TOutputImage>
void
GeneralToBSplineDisplacementFieldFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Ignore pixel value: " << this->m_IgnorePixelValue << std::endl;
  os << indent << "Number of control points: " << this->m_NumberOfControlPoints << std::endl;
  os << indent << "Number of levels: " << this->m_NumberOfLevels << std::endl;
  os << indent << "Spline order: " << this->m_SplineOrder << std::endl;
}
} // end namespace itk

#endif
