/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "ReadWriteData.h"
#include "antsUtilities.h"

#include "itkAddImageFilter.h"
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
#include "itkFastMarchingImageFilterBase.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"
#include "itkFlatStructuringElement.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkPadImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "itkImageFileWriter.h"


namespace ants
{

/*
template <typename ImageType>
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


template <typename ImageType>
typename ImageType::Pointer
iMathGE(typename ImageType::Pointer image, unsigned long radius) /*3*/
{
  const unsigned int                    ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;

  typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, StructuringElementType> FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(structuringElement);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathGO(typename ImageType::Pointer image, unsigned long radius) /*3*/
{
  const unsigned int                    ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;

  typedef itk::GrayscaleMorphologicalOpeningImageFilter<ImageType, ImageType, StructuringElementType> FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(structuringElement);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathGetLargestComponent(typename ImageType::Pointer image, /*3*/
                         unsigned long               smallest)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  if (image->GetNumberOfComponentsPerPixel() != 1)
  {
    // NOPE
  }

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef itk::Image<unsigned long, ImageDimension>                          LabelImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType, LabelImageType>         ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> FilterType;
  typedef itk::RelabelComponentImageFilter<LabelImageType, ImageType>        RelabelType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer          filter = FilterType::New();
  typename RelabelType::Pointer         relabel = RelabelType::New();

  threshold->SetInput(image);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(0.25); // FIXME - why these values?
  threshold->SetUpperThreshold(1.e9);
  threshold->Update();

  filter->SetInput(threshold->GetOutput());
  filter->SetFullyConnected(0);
  filter->Update();
  relabel->SetInput(filter->GetOutput());
  relabel->SetMinimumObjectSize(smallest);
  //    relabel->SetUseHistograms(true);

  try
  {
    relabel->Update();
  }
  catch (const itk::ExceptionObject & itkNotUsed(excep))
  {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
  }

  //  ANTs::WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
  //  return 0;
  typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(relabel->GetOutput(), 0);
  // typename ImageType::Pointer Clusters=relabel->GetOutput();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter(relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());

  float                     maximum = relabel->GetNumberOfObjects();
  float                     maxtstat = 0;
  std::vector<unsigned int> histogram((int)maximum + 1);
  std::vector<float>        clustersum((int)maximum + 1);
  for (int i = 0; i <= maximum; i++)
  {
    histogram[i] = 0;
    clustersum[i] = 0;
  }
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (vfIter.Get() > 0)
    {
      float vox = image->GetPixel(vfIter.GetIndex());
      histogram[(unsigned int)vfIter.Get()] = histogram[(unsigned int)vfIter.Get()] + 1;
      clustersum[(unsigned int)vfIter.Get()] += vox;
      if (vox > maxtstat)
      {
        maxtstat = vox;
      }
    }
  }
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (vfIter.Get() > 0)
    {
      Clusters->SetPixel(vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]);
      //  if ( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      //    maximgval=Clusters->GetPixel( vfIter.GetIndex());
    }
    else
    {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
    }
  }

  float maximgval = 0;
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (Clusters->GetPixel(vfIter.GetIndex()) > maximgval)
    {
      maximgval = Clusters->GetPixel(vfIter.GetIndex());
    }
  }

  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (Clusters->GetPixel(vfIter.GetIndex()) >= maximgval)
    {
      image->SetPixel(vfIter.GetIndex(), 1);
    }
    else
    {
      image->SetPixel(vfIter.GetIndex(), 0);
    }
  }

  return image;
}

template <typename ImageType>
typename ImageType::Pointer
iMathGrad(typename ImageType::Pointer image, double sigma, bool normalize) /*3*/
{

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer                                                     grad = FilterType::New();
  grad->SetInput(image);
  grad->SetSigma(sigma);
  grad->Update();

  typename ImageType::Pointer output = grad->GetOutput();
  if (normalize)
  {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer                            rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(1);
    rescaler->SetInput(grad->GetOutput());
    rescaler->Update();
    output = rescaler->GetOutput();
  }

  return output;
}

template <typename ImageType>
typename ImageType::Pointer
iMathHistogramEqualization(typename ImageType::Pointer image, double alpha, double beta, unsigned int r) /*3*/
{

  if (image->GetNumberOfComponentsPerPixel() != 1)
  {
    // NOPE
  }

  typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> AdaptiveHistogramEqualizationImageFilterType;
  typename AdaptiveHistogramEqualizationImageFilterType::Pointer   adaptiveHistogramEqualizationImageFilter =
    AdaptiveHistogramEqualizationImageFilterType::New();
  adaptiveHistogramEqualizationImageFilter->SetInput(image);
  typename AdaptiveHistogramEqualizationImageFilterType::RadiusType radius;
  radius.Fill(r);
  adaptiveHistogramEqualizationImageFilter->SetRadius(radius);
  adaptiveHistogramEqualizationImageFilter->SetAlpha(alpha);
  adaptiveHistogramEqualizationImageFilter->SetBeta(beta);
  adaptiveHistogramEqualizationImageFilter->Update();

  return adaptiveHistogramEqualizationImageFilter->GetOutput();
}


//
// shape (1=ball, 2=box, 3=cross, 4=annulus, 5=polygon)
template <typename ImageType>
typename ImageType::Pointer
iMathGD(typename ImageType::Pointer image, unsigned long radius) /*0*/ /*3*/
{

  const unsigned int                    ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;

  typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType, StructuringElementType> FilterType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(structuringElement);
  filter->Update();

  return filter->GetOutput();
}

} // namespace ants
