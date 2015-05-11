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
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"


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
iMathGC(typename ImageType::Pointer image, unsigned long radius)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

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
iMathGD(typename ImageType::Pointer image, unsigned long radius)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

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
iMathGE( typename ImageType::Pointer image, unsigned long radius)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType >   FilterType;

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
iMathGO( typename ImageType::Pointer image, unsigned long radius)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::GrayscaleMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

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
iMathLaplacian(typename ImageType::Pointer image, double sigma, bool normalize )
{

  typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer laplacian = FilterType::New();
  laplacian->SetInput( image );
  laplacian->SetSigma( sigma );
  laplacian->Update();

  typename ImageType::Pointer output = laplacian->GetOutput();
  if ( normalize )
    {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum( 0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->SetInput( laplacian->GetOutput() );
    rescaler->Update();
    output = rescaler->GetOutput();
    }

  return output;
}

template <class ImageType>
typename ImageType::Pointer
iMathMaurerDistance(typename ImageType::Pointer image,
                    typename ImageType::PixelType foreground )
{

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image);
  thresholder->SetLowerThreshold( foreground );
  thresholder->SetUpperThreshold( foreground );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( thresholder->GetOutput() );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false );
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathMC(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType closeValue)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType                         PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension>
    StructuringElementType;

  typedef itk::BinaryMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType >  FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetKernel( structuringElement );
  filter->SetForegroundValue( closeValue );
  //filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType dilateValue)
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
  filter->SetDilateValue( dilateValue );
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();

}

template <class ImageType>
typename ImageType::Pointer
iMathME( typename ImageType::Pointer image, unsigned long radius,
         typename ImageType::PixelType erodeValue )
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
  filter->SetErodeValue( erodeValue );
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathMO( typename ImageType::Pointer image, unsigned long radius,
         typename ImageType::PixelType openValue )
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
  filter->SetForegroundValue( openValue );
  filter->SetBackgroundValue( 0 );
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
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> NormFilterType;
  typename NormFilterType::Pointer normFilter = NormFilterType::New();
  normFilter->SetInput( image );
  normFilter->SetOutputMinimum( itk::NumericTraits<PixelType>::ZeroValue() );
  normFilter->SetOutputMaximum( itk::NumericTraits<PixelType>::OneValue() );
  normFilter->Update();

  return normFilter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathPeronaMalik( typename ImageType::Pointer image, unsigned long nIterations,
  double conductance )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType, ImageType >
    FilterType;
  typedef typename FilterType::TimeStepType             TimeStepType;

  // Select time step size.
  TimeStepType  spacingsize = 0;
  for( unsigned int d = 0; d < ImageType::ImageDimension; d++ )
    {
    TimeStepType sp = image->GetSpacing()[d];
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );

  // FIXME - cite reason for this step
  double dimPlusOne = ImageType::ImageDimension + 1;
  TimeStepType mytimestep = spacingsize / vcl_pow( 2.0 , dimPlusOne );
  TimeStepType reftimestep = 0.4 / vcl_pow( 2.0 , dimPlusOne );
  if ( mytimestep > reftimestep )
    {
    mytimestep = reftimestep;
    }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetConductanceParameter( conductance ); // might need to change this
  filter->SetNumberOfIterations( nIterations );
  filter->SetTimeStep( mytimestep );

  filter->Update();
  return filter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathSharpen( typename ImageType::Pointer image )
{
  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::Pointer                   ImagePointerType;

  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer sharpenFilter = FilterType::New();
  sharpenFilter->SetInput( image );
  sharpenFilter->Update();

  return sharpenFilter->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer
iMathTruncateIntensity( typename ImageType::Pointer image, double lowerQ, double upperQ, int nBins,
                        typename itk::Image<unsigned int, ImageType::ImageDimension>::Pointer mask )
{

  typedef typename ImageType::PixelType                     PixelType;
  typedef typename ImageType::Pointer                       ImagePointerType;
  typedef unsigned int                                      LabelType;
  typedef itk::Image<LabelType, ImageType::ImageDimension>  MaskType;

  if( mask.IsNull() )
    {
    typedef itk::BinaryThresholdImageFilter<ImageType, MaskType> ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresh = ThresholdFilterType::New();
    thresh->SetInput( image );
    thresh->SetLowerThreshold( 1e-6 );
    thresh->SetUpperThreshold( itk::NumericTraits<PixelType>::max() );
    thresh->SetInsideValue(1);
    thresh->SetOutsideValue(0);
    thresh->Update();
    mask = thresh->GetOutput();
    }
  typedef itk::LabelStatisticsImageFilter<ImageType, MaskType> HistogramFilterType;
  typename HistogramFilterType::Pointer stats = HistogramFilterType::New();

  stats->SetInput( image );
  stats->SetLabelInput( mask );
  stats->Update();
  PixelType minValue = stats->GetMinimum( 1 );
  PixelType maxValue = stats->GetMaximum( 1 );

  // Hack increment by delta
  if (minValue == 0)
    {
    minValue = (PixelType) (minValue + 1e-6);
    }
  if (minValue == 0)
    {
    minValue++;
    }

  stats->SetUseHistograms( true );
  stats->SetHistogramParameters( nBins, minValue, maxValue );
  stats->Update();

  typedef typename HistogramFilterType::HistogramPointer HistogramPointer;
  HistogramPointer histogram = stats->GetHistogram( 1 );

  PixelType lowerQuantile = histogram->Quantile( 0, lowerQ );
  PixelType upperQuantile = histogram->Quantile( 0, upperQ );

  typedef itk::IntensityWindowingImageFilter<ImageType,ImageType> WindowFilterType;
  typename WindowFilterType::Pointer windowFilter = WindowFilterType::New();
  windowFilter->SetInput( image );
  windowFilter->SetWindowMinimum( lowerQuantile );
  windowFilter->SetOutputMinimum( lowerQuantile );
  windowFilter->SetWindowMaximum( upperQuantile );
  windowFilter->SetOutputMaximum( upperQuantile );
  windowFilter->Update();

  return windowFilter->GetOutput();

}


} // namespace ants
