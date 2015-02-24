/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageIntensityAndGradientToPointSetFilter_hxx
#define __itkImageIntensityAndGradientToPointSetFilter_hxx

#include "itkImageIntensityAndGradientToPointSetFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

namespace itk
{
//
// Constructor
//
template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>
::ImageIntensityAndGradientToPointSetFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 2 );

  this->m_NeighborhoodRadius.Fill( 1 );
  this->m_Sigma = 1.5;

  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}

template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
void
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>
::GenerateData()
{
  const InputImageType * inputImage = this->GetInputImage();
  const MaskImageType * maskImage = this->GetMaskImage();

  // Calculate gradient image

  typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
  gradientFilter->SetInput( inputImage );
  gradientFilter->SetSigma( this->m_Sigma );

  typename GradientImageType::Pointer gradientImage = gradientFilter->GetOutput();
  gradientImage->Update();
  gradientImage->DisconnectPipeline();

  // Set up the point set pixel type characteristics
  //    size of pixel type array = intensities + gradientXs + gradientYs + ...
  //                             = ( numberOfNeighborhoodVoxels + 1 ) * Dimension

  SizeValueType numberOfNeighborhoodVoxels = 1;
  for( SizeValueType d = 0; d < Dimension; d++ )
    {
    numberOfNeighborhoodVoxels *= ( 2 * this->m_NeighborhoodRadius[d] + 1 );
    }
  const SizeValueType sizeOfPixelTypeArray = ( numberOfNeighborhoodVoxels + 1 ) * Dimension;

  // Initialize output point set

  typename OutputMeshType::Pointer output = OutputMeshType::New();
  output->Initialize();

  SizeValueType count = 0;

  ConstNeighborhoodIteratorType ItN( this->m_NeighborhoodRadius, gradientImage,
    gradientImage->GetRequestedRegion() );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    typename InputImageType::IndexType index = ItN.GetIndex();

    if( maskImage->GetPixel( index ) && ItN.InBounds() )
      {
      typename InputImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint( index, imagePoint );

      PointType point;
      point.CastFrom( imagePoint );

      PointSetPixelType array( sizeOfPixelTypeArray );

      unsigned int arrayIndex = 0;
      for( SizeValueType n = 0; n < numberOfNeighborhoodVoxels; n++ )
        {
        typename InputImageType::IndexType offsetIndex = ItN.GetIndex( n );

        InputImagePixelType intensity = inputImage->GetPixel( offsetIndex );
        GradientPixelType gradient = gradientImage->GetPixel( offsetIndex );

        array[arrayIndex++] = intensity;
        for( SizeValueType d = 0; d < Dimension; d++ )
          {
          array[arrayIndex++] = gradient[d];
          }

        output->SetPoint( count, point );
        output->SetPointData( count++, array );
        }
      }
    }

  this->GraftOutput( output );
}

template <typename TInputImage, typename TMaskImage, typename TOutputMesh>
void
ImageIntensityAndGradientToPointSetFilter<TInputImage, TMaskImage, TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "Sigma = " << this->m_Sigma << std::endl;
  os << "Neighborhood radius = " << this->m_NeighborhoodRadius << std::endl;
}
} // end of namespace itk

#endif
