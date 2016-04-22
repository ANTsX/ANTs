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

#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"


namespace ants
{

/*
template <class ImageType>
typename ImageType::Pointer
BlobCorrespondence( typename ImageType::Pointer image, unsigned int nBlobs,
              typename ImageType::Pointer itkNotUsed(image2),
              double itkNotUsed(corrThresh), double itkNotUsed(radius), double itkNotUsed(distanceThresh) )
{
  typedef float RealType;

  // sensitive parameters are set here - begin
  //RealType     gradsig = 1.0;      // sigma for gradient filter
  unsigned int stepsperoctave = 10; // number of steps between doubling of scale
  RealType     minscale = std::pow( 1.0, 1.0 );
  RealType     maxscale = std::pow( 2.0, 10.0 );
  //RealType     uniqfeat_thresh = 0.01;
  //RealType     smallval = 1.e-2; // assumes images are normalizes in [ 0, 1 ]
  //bool         dosinkhorn = false;
  //RealType     maxradiusdiffallowed = 0.25; // IMPORTANT feature size difference
  //RealType     kneighborhoodval = 3;        // IMPORTANT - defines how many nhood nodes to use in k-hood definition
  //unsigned int radval = 20;                 // IMPORTANT radius for correlation
  //RealType     dthresh = 0.02;              // IMPORTANT distance preservation threshold
  // sensitive parameters are set here - end
}
*/

template <class ImageType>
typename ImageType::Pointer
iMathBlobDetector( typename ImageType::Pointer image, unsigned int nBlobs )
{
  typedef float RealType;

  unsigned int stepsperoctave = 10; // number of steps between doubling of scale
  RealType     minscale = std::pow( 1.0, 1.0 );
  RealType     maxscale = std::pow( 2.0, 10.0 );

  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typename BlobFilterType::Pointer blobFilter = BlobFilterType::New();
  blobFilter->SetStartT( minscale );
  blobFilter->SetEndT( maxscale );
  blobFilter->SetStepsPerOctave( stepsperoctave );
  blobFilter->SetNumberOfBlobs( nBlobs );
  blobFilter->SetInput( image );
  blobFilter->Update();

  typedef typename BlobFilterType::BlobRadiusImageType BlobRadiusImageType;
  typename BlobRadiusImageType::Pointer labimg = blobFilter->GetBlobRadiusImage();

  return( labimg );
}

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
iMathDistanceMap( typename ImageType::Pointer image, bool useSpacing )
{
  typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;

  typename  FilterType::Pointer filter = FilterType::New();
  filter->InputIsBinaryOff();
  filter->SetUseImageSpacing(useSpacing);
  filter->SetInput(image);
  filter->Update();

  return filter->GetOutput();
}


// algorithm :
// 1. get distance map of object
// 2. threshold
// 3. label connected components
// 4. label surface
// 5. if everywhere on surface is next to object then it's a hole
// 6. make sure it's not the background
template <class ImageType>
typename ImageType::Pointer
iMathFillHoles( typename ImageType::Pointer image, double holeParam )
{

  if ( (holeParam < 0) || (holeParam > 2) )
    {
    //itk::itkExceptionMacro("FillHoles: holeParam value must lie in [0,2]");
    }

  typedef typename ImageType::Pointer                ImagePointerType;
  typedef itk::Image<int, ImageType::ImageDimension> MaskType;
  typedef typename ImageType::PixelType              PixelType;
  typedef typename MaskType::PixelType               LabelType;

  const PixelType imageMax = itk::NumericTraits<PixelType>::max();
  const LabelType labelMax = itk::NumericTraits<LabelType>::max();
  const PixelType objectMin = 0.5;
  const PixelType distanceMin = 0.001;

  typedef itk::CastImageFilter<MaskType,ImageType>            MaskToImage;
  typedef itk::BinaryThresholdImageFilter<ImageType,MaskType> ThresholdFilterType;
  typedef itk::BinaryThresholdImageFilter<MaskType,MaskType>  ThresholdMaskFilterType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  threshold->SetInput( image );
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(objectMin);
  threshold->SetUpperThreshold(imageMax);

  typedef itk::DanielssonDistanceMapImageFilter<MaskType, ImageType> FilterType;
  typename  FilterType::Pointer distance = FilterType::New();
  distance->InputIsBinaryOff();
  distance->SetUseImageSpacing(false);
  distance->SetInput(threshold->GetOutput());

  typename ThresholdFilterType::Pointer dThreshold = ThresholdFilterType::New();
  dThreshold->SetInput( distance->GetOutput() );
  dThreshold->SetInsideValue(1);
  dThreshold->SetOutsideValue(0);
  dThreshold->SetLowerThreshold(distanceMin);
  dThreshold->SetUpperThreshold(imageMax);
  dThreshold->Update();

  typedef itk::ConnectedComponentImageFilter<MaskType,MaskType> ConnectedFilterType;
  typename ConnectedFilterType::Pointer connected = ConnectedFilterType::New();
  connected->SetInput( dThreshold->GetOutput() );
  connected->SetFullyConnected( false );

  typedef itk::RelabelComponentImageFilter<MaskType, MaskType>  RelabelFilterType;
  typename RelabelFilterType::Pointer relabel = RelabelFilterType::New();
  relabel->SetInput( connected->GetOutput() );
  relabel->SetMinimumObjectSize( 0 );
  relabel->Update();

  if( holeParam == 2 )
    {
    typename ThresholdMaskFilterType::Pointer oThreshold = ThresholdMaskFilterType::New();
    oThreshold->SetInput( relabel->GetOutput() );
    oThreshold->SetInsideValue(1);
    oThreshold->SetOutsideValue(0);
    oThreshold->SetLowerThreshold(2);
    oThreshold->SetUpperThreshold(labelMax);

    typedef itk::AddImageFilter<MaskType,MaskType> AddFilterType;
    typename AddFilterType::Pointer add = AddFilterType::New();
    add->SetInput1( threshold->GetOutput() );
    add->SetInput2( oThreshold->GetOutput() );

    typename MaskToImage::Pointer maskToImage = MaskToImage::New();
    maskToImage->SetInput( add->GetOutput() );
    maskToImage->Update();

    return maskToImage->GetOutput();
    }

  // FIXME - add filter for below -- avoid iterators in these functions
  typename MaskToImage::Pointer caster = MaskToImage::New();
  caster->SetInput( threshold->GetOutput() );
  caster->Update();
  ImagePointerType imageout = caster->GetOutput();

  typedef itk::NeighborhoodIterator<MaskType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageType::ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  // now we have the exact number of objects labeled independently
  for( int lab = 2; lab <= maximum; lab++ )
    {
    float erat = 2;
    if( holeParam <= 1 )
      {
      GHood.GoToBegin();

      unsigned long objectedge = 0;
      unsigned long backgroundedge = 0;
      unsigned long totaledge = 0;
      unsigned long volume = 0;

      while( !GHood.IsAtEnd() )
        {
        typename ImageType::PixelType p = GHood.GetCenterPixel();
        typename ImageType::IndexType ind2;
        if( p == lab )
          {
          volume++;
          for( unsigned int i = 0; i < GHood.Size(); i++ )
            {
            ind2 = GHood.GetIndex(i);
            float val2 = threshold->GetOutput()->GetPixel(ind2);
            if( (val2 == 1) && GHood.GetPixel(i) != lab )
              {
              objectedge++;
              totaledge++;
              }
            else if( (val2 == 1) && GHood.GetPixel(i) != lab )
              {
              backgroundedge++;
              totaledge++;
              }
            }
          }
        ++GHood;
        }

      erat = (float)objectedge / (float)totaledge;
      }

    if( erat > holeParam ) // fill the hole
      {
      // std::cout << " Filling " << lab << " of " << maximum <<  std::endl;
      typedef itk::ImageRegionIteratorWithIndex<MaskType> RelabelIterator;
      RelabelIterator vfIter( relabel->GetOutput(),
                              relabel->GetOutput()->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        if( vfIter.Get() == lab )
          {
          imageout->SetPixel(vfIter.GetIndex(), 1);
          }
        }
      }
    }

  return imageout;
}


template <class ImageType>
typename ImageType::Pointer
iMathGC(typename ImageType::Pointer image, unsigned long radius)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType            PixelType;

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
iMathGrad(typename ImageType::Pointer image, double sigma, bool normalize )
{

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer grad = FilterType::New();
  grad->SetInput( image );
  grad->SetSigma( sigma );
  grad->Update();

  typename ImageType::Pointer output = grad->GetOutput();
  if ( normalize )
    {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum( 0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->SetInput( grad->GetOutput() );
    rescaler->Update();
    output = rescaler->GetOutput();
    }

  return output;
}

template <class ImageType>
typename ImageType::Pointer
iMathHistogramEqualization( typename ImageType::Pointer image, double alpha, double beta, unsigned int r )
{

  if ( image->GetNumberOfComponentsPerPixel() != 1 )
    {
    // NOPE
    }

  typedef itk::AdaptiveHistogramEqualizationImageFilter< ImageType > AdaptiveHistogramEqualizationImageFilterType;
  typename AdaptiveHistogramEqualizationImageFilterType::Pointer adaptiveHistogramEqualizationImageFilter = AdaptiveHistogramEqualizationImageFilterType::New();
  adaptiveHistogramEqualizationImageFilter->SetInput( image );
  typename AdaptiveHistogramEqualizationImageFilterType::RadiusType radius;
  radius.Fill( r );
  adaptiveHistogramEqualizationImageFilter->SetRadius(radius);
  adaptiveHistogramEqualizationImageFilter->SetAlpha(alpha);
  adaptiveHistogramEqualizationImageFilter->SetBeta(beta);
  adaptiveHistogramEqualizationImageFilter->Update( );

  return adaptiveHistogramEqualizationImageFilter->GetOutput();
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
iMathPad( typename ImageType::Pointer image1, int padvalue )
{
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  typename ImageType::SizeType size = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::PointType origin2 = image1->GetOrigin();
  typename ImageType::SizeType newsize = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::RegionType newregion;
  // determine new image size
  for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
      {
      float dimsz = (float)size[i];
      newsize[i] = (unsigned int)(dimsz + padvalue * 2);
      }
  newregion.SetSize(newsize);
  newregion.SetIndex(image1->GetLargestPossibleRegion().GetIndex() );
  typename ImageType::Pointer padimage =
      AllocImage<ImageType>(newregion,
                            image1->GetSpacing(),
                            origin2,
                            image1->GetDirection(), 0);

  typename ImageType::IndexType index;
  typename ImageType::IndexType index2;
  if( padvalue > 0 )
      {
      index.Fill(0);
      index2.Fill( (unsigned int)fabs( static_cast<float>( padvalue ) ) );
      }
  else
      {
      index2.Fill(0);
      index.Fill( (unsigned int)fabs( static_cast<float>( padvalue ) ) );
      }

  typename ImageType::PointType point1, pointpad;
  image1->TransformIndexToPhysicalPoint(index, point1);
  padimage->TransformIndexToPhysicalPoint(index2, pointpad);

  for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
      {
      origin2[i] += (point1[i] - pointpad[i]);
      }

  padimage->SetOrigin(origin2);

  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      typename ImageType::IndexType oindex = iter.GetIndex();
      typename ImageType::IndexType padindex = iter.GetIndex();

      bool isinside = true;
      for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
        {
        float shifted = ( (float)oindex[i] + padvalue);
        if( shifted < 0 || shifted > newsize[i] - 1 )
          {
          isinside = false;
          }
        //      if (shifted < 0) shifted=0;
        // padindex[i]=
        }
      if( isinside )
        {
        for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
          {
          float shifted = ( (float)oindex[i] + padvalue);
          padindex[i] = (unsigned int)shifted;
          }
        padimage->SetPixel(padindex, iter.Get() );
        }
      }
return padimage;
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
  TimeStepType mytimestep = spacingsize / std::pow( 2.0 , dimPlusOne );
  TimeStepType reftimestep = 0.4 / std::pow( 2.0 , dimPlusOne );
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
