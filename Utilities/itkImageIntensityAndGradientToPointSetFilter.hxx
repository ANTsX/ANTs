/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkImageIntensityAndGradientToPointSetFilter_hxx
#define itkImageIntensityAndGradientToPointSetFilter_hxx


#include "itkCentralDifferenceImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

namespace itk
{
//
// Constructor
//
template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>::
  ImageIntensityAndGradientToPointSetFilter()
  : m_Sigma(1.5)
  , m_NeighborhoodRadius()
  , m_UseCentralDifferenceFunction(true)
{
  this->m_NeighborhoodRadius.Fill(1), this->ProcessObject::SetNumberOfRequiredInputs(2);

  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());
}

template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
void
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>::Update()
{
  this->GenerateData();
}

template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
void
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>::GenerateData()
{
  const InputImageType * inputImage = this->GetInputImage();
  const MaskImageType *  maskImage = this->GetMaskImage();

  typename OutputMeshType::Pointer output = this->GetOutput();

  // Calculate gradient image

  typename GradientImageType::Pointer gradientImage = nullptr;
  if (this->m_UseCentralDifferenceFunction)
  {
    GradientPixelType zeroVector;
    zeroVector.Fill(0);

    gradientImage = GradientImageType::New();
    gradientImage->CopyInformation(inputImage);
    gradientImage->SetRegions(inputImage->GetRequestedRegion());
    gradientImage->AllocateInitialized();

    typedef CentralDifferenceImageFunction<InputImageType, InputImagePixelType, GradientPixelType>
                                             GradientCalculatorType;
    typename GradientCalculatorType::Pointer gradientCalculator = GradientCalculatorType::New();
    gradientCalculator->SetInputImage(inputImage);
    gradientCalculator->SetUseImageDirection(true);

    ImageRegionIteratorWithIndex<GradientImageType> It(gradientImage, gradientImage->GetRequestedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      It.Set(gradientCalculator->EvaluateAtIndex(It.GetIndex()));
    }
  }
  else
  {
    typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
    gradientFilter->SetInput(inputImage);
    gradientFilter->SetSigma(this->m_Sigma);
    gradientFilter->SetUseImageDirection(true);

    gradientImage = gradientFilter->GetOutput();
    gradientImage->Update();
    gradientImage->DisconnectPipeline();
  }

  // Set up the point set pixel type characteristics
  //    size of pixel type array = intensities + gradientXs + gradientYs + ...
  //                             = ( numberOfNeighborhoodVoxels + 1 ) * Dimension

  SizeValueType numberOfNeighborhoodVoxels = 1;
  for (SizeValueType d = 0; d < Dimension; d++)
  {
    numberOfNeighborhoodVoxels *= (2 * this->m_NeighborhoodRadius[d] + 1);
  }
  const SizeValueType sizeOfPixelTypeArray = numberOfNeighborhoodVoxels * (1 + Dimension);

  // Initialize output point set

  SizeValueType count = 0;

  ConstNeighborhoodIteratorType ItN(this->m_NeighborhoodRadius, gradientImage, gradientImage->GetRequestedRegion());
  for (ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN)
  {
    typename InputImageType::IndexType index = ItN.GetIndex();

    if (maskImage->GetPixel(index) && ItN.InBounds())
    {
      typename InputImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint(index, imagePoint);

      PointType point;
      point.CastFrom(imagePoint);

      PointSetPixelType array(sizeOfPixelTypeArray);

      unsigned int arrayIndex = 0;
      for (SizeValueType n = 0; n < numberOfNeighborhoodVoxels; n++)
      {
        typename InputImageType::IndexType offsetIndex = ItN.GetIndex(n);

        InputImagePixelType intensity = inputImage->GetPixel(offsetIndex);
        GradientPixelType   gradient = gradientImage->GetPixel(offsetIndex);

        array[arrayIndex++] = intensity;
        for (SizeValueType d = 0; d < Dimension; d++)
        {
          array[arrayIndex++] = gradient[d];
        }
      }
      output->SetPoint(count, point);
      output->SetPointData(count++, array);
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
void
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>::PrintSelf(std::ostream & os,
                                                                                           Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Sigma = " << this->m_Sigma << std::endl;
  os << "Neighborhood radius = " << this->m_NeighborhoodRadius << std::endl;
}
} // end of namespace itk

#endif
