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
// #include "itkFlatStructuringElement.h"
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
unsigned int
morph_shape_flag(const char * shape)
{
  std::string shapeStr(shape);
  std::transform(shapeStr.begin(), shapeStr.end(), shapeStr.begin(), ::tolower);

  unsigned int flag = 1;

  if (!shapeStr.compare("ball"))
  {
    flag = 1;
  }
  else if (!shapeStr.compare("box"))
  {
    flag = 2;
  }
  if (!shapeStr.compare("cross"))
  {
    flag = 3;
  }
  if (!shapeStr.compare("annulus"))
  {
    flag = 4;
  }
  if (!shapeStr.compare("polygon"))
  {
    flag = 5;
  }

  return flag;
}

template <typename ImageType>
typename ImageType::Pointer
iMathBlobDetector(typename ImageType::Pointer image, unsigned int nBlobs) /*?????*/
{
  typedef float RealType;

  unsigned int stepsperoctave = 10; // number of steps between doubling of scale
  RealType     minscale = std::pow(1.0, 1.0);
  RealType     maxscale = std::pow(2.0, 10.0);

  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typename BlobFilterType::Pointer                                   blobFilter = BlobFilterType::New();
  blobFilter->SetStartT(minscale);
  blobFilter->SetEndT(maxscale);
  blobFilter->SetStepsPerOctave(stepsperoctave);
  blobFilter->SetNumberOfBlobs(nBlobs);
  blobFilter->SetInput(image);
  blobFilter->Update();

  typedef typename BlobFilterType::BlobRadiusImageType BlobRadiusImageType;
  typename BlobRadiusImageType::Pointer                labimg = blobFilter->GetBlobRadiusImage();

  return (labimg);
}

template <typename ImageType>
typename ImageType::Pointer
iMathCanny(typename ImageType::Pointer image, /*0*/
           double                      sigma,
           double                      lowerThreshold,
           double                      upperThreshold)
{

  typedef typename ImageType::PixelType                            PixelType;
  typedef itk::CannyEdgeDetectionImageFilter<ImageType, ImageType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetVariance(sigma);
  filter->SetUpperThreshold((PixelType)upperThreshold);
  filter->SetLowerThreshold((PixelType)lowerThreshold);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathDistanceMap(typename ImageType::Pointer image, bool useSpacing) /*0*/
{
  typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
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
template <typename ImageType>
typename ImageType::Pointer
iMathFillHoles(typename ImageType::Pointer image, double holeParam) /*0*/
{

  if ((holeParam < 0) || (holeParam > 2))
  {
    // itk::itkExceptionMacro("FillHoles: holeParam value must lie in [0,2]");
  }

  typedef typename ImageType::Pointer                ImagePointerType;
  typedef itk::Image<int, ImageType::ImageDimension> MaskType;
  typedef typename ImageType::PixelType              PixelType;
  typedef typename MaskType::PixelType               LabelType;

  const PixelType imageMax = itk::NumericTraits<PixelType>::max();
  const LabelType labelMax = itk::NumericTraits<LabelType>::max();
  PixelType       objectMin = 0.5;
  PixelType       distanceMin = 0.001;

  typedef itk::CastImageFilter<MaskType, ImageType>            MaskToImage;
  typedef itk::BinaryThresholdImageFilter<ImageType, MaskType> ThresholdFilterType;
  typedef itk::BinaryThresholdImageFilter<MaskType, MaskType>  ThresholdMaskFilterType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  threshold->SetInput(image);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(objectMin);
  threshold->SetUpperThreshold(imageMax);

  typedef itk::DanielssonDistanceMapImageFilter<MaskType, ImageType> FilterType;
  typename FilterType::Pointer                                       distance = FilterType::New();
  distance->InputIsBinaryOff();
  distance->SetUseImageSpacing(false);
  distance->SetInput(threshold->GetOutput());

  typename ThresholdFilterType::Pointer dThreshold = ThresholdFilterType::New();
  dThreshold->SetInput(distance->GetOutput());
  dThreshold->SetInsideValue(1);
  dThreshold->SetOutsideValue(0);
  dThreshold->SetLowerThreshold(distanceMin);
  dThreshold->SetUpperThreshold(imageMax);
  dThreshold->Update();

  typedef itk::ConnectedComponentImageFilter<MaskType, MaskType> ConnectedFilterType;
  typename ConnectedFilterType::Pointer                          connected = ConnectedFilterType::New();
  connected->SetInput(dThreshold->GetOutput());
  connected->SetFullyConnected(false);

  typedef itk::RelabelComponentImageFilter<MaskType, MaskType> RelabelFilterType;
  typename RelabelFilterType::Pointer                          relabel = RelabelFilterType::New();
  relabel->SetInput(connected->GetOutput());
  relabel->SetMinimumObjectSize(0);
  relabel->Update();

  if (itk::Math::FloatAlmostEqual(holeParam, static_cast<double>(2.0)))
  {
    typename ThresholdMaskFilterType::Pointer oThreshold = ThresholdMaskFilterType::New();
    oThreshold->SetInput(relabel->GetOutput());
    oThreshold->SetInsideValue(1);
    oThreshold->SetOutsideValue(0);
    oThreshold->SetLowerThreshold(2);
    oThreshold->SetUpperThreshold(labelMax);

    typedef itk::AddImageFilter<MaskType, MaskType> AddFilterType;
    typename AddFilterType::Pointer                 add = AddFilterType::New();
    add->SetInput1(threshold->GetOutput());
    add->SetInput2(oThreshold->GetOutput());

    typename MaskToImage::Pointer maskToImage = MaskToImage::New();
    maskToImage->SetInput(add->GetOutput());
    maskToImage->Update();

    return maskToImage->GetOutput();
  }

  // FIXME - add filter for below -- avoid iterators in these functions
  typename MaskToImage::Pointer caster = MaskToImage::New();
  caster->SetInput(threshold->GetOutput());
  caster->Update();
  ImagePointerType imageout = caster->GetOutput();

  typedef itk::NeighborhoodIterator<MaskType> iteratorType;
  typename iteratorType::RadiusType           rad;
  for (unsigned int j = 0; j < ImageType::ImageDimension; j++)
  {
    rad[j] = 1;
  }
  iteratorType GHood(rad, relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());

  float maximum = relabel->GetNumberOfObjects();
  // now we have the exact number of objects labeled independently
  for (int lab = 2; lab <= maximum; lab++)
  {
    float erat = 2;
    if (holeParam <= 1)
    {
      GHood.GoToBegin();

      unsigned long objectedge = 0;
      unsigned long totaledge = 0;

      while (!GHood.IsAtEnd())
      {
        typename ImageType::PixelType p = GHood.GetCenterPixel();
        typename ImageType::IndexType ind2;
        if (itk::Math::FloatAlmostEqual(p, static_cast<typename ImageType::PixelType>(lab)))
        {
          for (unsigned int i = 0; i < GHood.Size(); i++)
          {
            ind2 = GHood.GetIndex(i);
            float val2 = threshold->GetOutput()->GetPixel(ind2);
            if (itk::Math::FloatAlmostEqual(val2, itk::NumericTraits<float>::OneValue()) && GHood.GetPixel(i) != lab)
            {
              objectedge++;
              totaledge++;
            }
            else if (itk::Math::FloatAlmostEqual(val2, itk::NumericTraits<float>::OneValue()) &&
                     GHood.GetPixel(i) != lab)
            {
              totaledge++;
            }
          }
        }
        ++GHood;
      }

      erat = (float)objectedge / (float)totaledge;
    }

    if (erat > static_cast<float>(holeParam)) // fill the hole
    {
      // std::cout << " Filling " << lab << " of " << maximum <<  std::endl;
      typedef itk::ImageRegionIteratorWithIndex<MaskType> RelabelIterator;
      RelabelIterator vfIter(relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());
      for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
      {
        if (vfIter.Get() == lab)
        {
          imageout->SetPixel(vfIter.GetIndex(), 1);
        }
      }
    }
  }

  return imageout;
}


template <typename ImageType>
typename ImageType::Pointer
iMathGC(typename ImageType::Pointer image, unsigned long radius) /*0*/
{

  const unsigned int                    ImageDimension = ImageType::ImageDimension;
  typedef typename ImageType::PixelType PixelType;

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;

  typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, ImageType, StructuringElementType> FilterType;

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
