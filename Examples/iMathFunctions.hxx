/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "iMathFunctions.h"
#include "ReadWriteData.h"
#include "antsUtilities.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkStatisticsImageFilter.h"

namespace ants
{

template <class ImageType>
typename ImageType::Pointer
iMathCanny( typename ImageType::Pointer image,
            double sigma,
            double lowerThreshold,
            double upperThreshold )
{

  typedef typename ImageType::PixelType            PixelType;
  typedef itk::CannyEdgeDetectionImageFilter< ImageType, ImageType >  FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetVariance( sigma );
  filter->SetUpperThreshold( (PixelType) upperThreshold );
  filter->SetLowerThreshold( (PixelType) lowerThreshold );
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathGetLargestComponent( typename ImageType::Pointer image,
                     unsigned long smallest )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  Iterator;

  // compute the voxel volume
  typename ImageType::SpacingType spacing = image->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= spacing[i];
    }

  typedef itk::Image<unsigned long, ImageDimension>                          LabelImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType, LabelImageType>         ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> FilterType;
  typedef itk::RelabelComponentImageFilter<LabelImageType, ImageType>        RelabelType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer filter = FilterType::New();
  typename RelabelType::Pointer relabel = RelabelType::New();

  threshold->SetInput(image);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(0.25);  //FIXME - why these values?
  threshold->SetUpperThreshold(1.e9);
  threshold->Update();

  filter->SetInput(threshold->GetOutput() );
  filter->SetFullyConnected( 0 );
  filter->Update();
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( smallest );
  //    relabel->SetUseHistograms(true);

  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
    }

  //  WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
  //  return 0;
  typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(relabel->GetOutput(), 0);
  // typename ImageType::Pointer Clusters=relabel->GetOutput();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  float                     maxtstat = 0;
  std::vector<unsigned int> histogram( (int)maximum + 1);
  std::vector<float>        clustersum( (int)maximum + 1);
  for( int i = 0; i <= maximum; i++ )
    {
    histogram[i] = 0;
    clustersum[i] = 0;
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      float vox = image->GetPixel(vfIter.GetIndex() );
      histogram[(unsigned int)vfIter.Get()] = histogram[(unsigned int)vfIter.Get()] + 1;
      clustersum[(unsigned int)vfIter.Get()] += vox;
      if( vox > maxtstat )
        {
        maxtstat = vox;
        }
      }
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      Clusters->SetPixel( vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]  );
      //  if ( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      //    maximgval=Clusters->GetPixel( vfIter.GetIndex());
      }
    else
      {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
      }
    }

  float maximgval = 0;
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      {
      maximgval = Clusters->GetPixel( vfIter.GetIndex() );
      }
    }

  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) >= maximgval )
      {
      image->SetPixel( vfIter.GetIndex(), 1);
      }
    else
      {
      image->SetPixel( vfIter.GetIndex(), 0);
      }
    }

  return image;
}

template <class ImageType>
typename ImageType::Pointer
iMathMD( typename ImageType::Pointer image, unsigned long radius )
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathME( typename ImageType::Pointer image, unsigned long radius )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType >   FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathMO( typename ImageType::Pointer image, unsigned long radius )
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathNormalize( typename ImageType::Pointer image )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  Iterator;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::StatisticsImageFilter<ImageType> StatsFilterType;
  typename StatsFilterType::Pointer statsFilter = StatsFilterType::New();
  statsFilter->SetInput( image );
  statsFilter->Update();
  PixelType max = statsFilter->GetMaximum();
  PixelType min = statsFilter->GetMinimum();

  typedef itk::ShiftScaleImageFilter<ImageType,ImageType> NormFilterType;
  typename NormFilterType::Pointer normFilter = NormFilterType::New();
  normFilter->SetInput( image );
  normFilter->SetShift( -min );
  normFilter->SetScale( 1.0/(max-min) );
  normFilter->Update();

  return normFilter->GetOutput();
}



} // namespace ants
