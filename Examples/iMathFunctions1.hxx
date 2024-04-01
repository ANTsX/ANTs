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

//
// shape (1=ball, 2=box, 3=cross, 4=annulus, 5=polygon)

template <unsigned int ImageDimension>
typename itk::FlatStructuringElement<ImageDimension>
iMathGetFlatStructuringElement(unsigned int  shape,
                               unsigned long radius,
                               bool          radiusIsParametric,
                               unsigned int  lines,
                               unsigned int  thickness,
                               bool          includeCenter)
{
  typedef typename itk::FlatStructuringElement<ImageDimension> ElementType;
  ElementType                                                  element;

  typename ElementType::RadiusType elRadius;
  elRadius.Fill(radius);

  switch (shape)
  {
    case 1:
      element = ElementType::Ball(elRadius, radiusIsParametric);
      break;
    case 2:
      element = ElementType::Box(elRadius);
      break;
    case 3:
      element = ElementType::Cross(elRadius);
      break;
    case 4:
      element = ElementType::Annulus(elRadius, thickness, includeCenter, radiusIsParametric);
      break;
    case 5:
      element = ElementType::Polygon(elRadius, lines);
      break;
    default:
      break;
  }

  return element;
}

template <typename ImageType>
typename ImageType::Pointer
iMathLaplacian(typename ImageType::Pointer image, double sigma, bool normalize) /*1*/
{

  typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer                                             laplacian = FilterType::New();
  laplacian->SetInput(image);
  laplacian->SetSigma(sigma);
  laplacian->Update();

  typename ImageType::Pointer output = laplacian->GetOutput();
  if (normalize)
  {
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer                            rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(1);
    rescaler->SetInput(laplacian->GetOutput());
    rescaler->Update();
    output = rescaler->GetOutput();
  }

  return output;
}

template <typename ImageType>
typename ImageType::Pointer
iMathMaurerDistance(typename ImageType::Pointer   image, /*1*/
                    typename ImageType::PixelType foreground)
{

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer                             thresholder = ThresholderType::New();
  thresholder->SetInput(image);
  thresholder->SetLowerThreshold(foreground);
  thresholder->SetUpperThreshold(foreground);
  thresholder->SetInsideValue(1);
  thresholder->SetOutsideValue(0);

  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer                                          filter = FilterType::New();
  filter->SetInput(thresholder->GetOutput());
  filter->SetSquaredDistance(false);
  filter->SetUseImageSpacing(true);
  filter->SetInsideIsPositive(false);
  filter->Update();

  return filter->GetOutput();
}


template <typename ImageType>
typename ImageType::Pointer
iMathMC(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType closeValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned int                  lines,
        unsigned int                  thickness,
        bool                          includeCenter)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef typename itk::FlatStructuringElement<ImageType::ImageDimension> ElementType;
  ElementType                                                             element =
    iMathGetFlatStructuringElement<ImageDimension>(shape, radius, radiusIsParametric, lines, thickness, includeCenter);

  typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, ElementType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(element);
  filter->SetForegroundValue(closeValue);
  // filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType dilateValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned int                  lines,
        unsigned int                  thickness,
        bool                          includeCenter)
{

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef typename itk::FlatStructuringElement<ImageType::ImageDimension> ElementType;
  ElementType                                                             element =
    iMathGetFlatStructuringElement<ImageDimension>(shape, radius, radiusIsParametric, lines, thickness, includeCenter);

  typedef itk::BinaryDilateImageFilter<ImageType, ImageType, ElementType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(element);
  filter->SetDilateValue(dilateValue);
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathME(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType erodeValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned int                  lines,
        unsigned int                  thickness,
        bool                          includeCenter)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef typename itk::FlatStructuringElement<ImageType::ImageDimension> ElementType;
  ElementType                                                             element =
    iMathGetFlatStructuringElement<ImageDimension>(shape, radius, radiusIsParametric, lines, thickness, includeCenter);

  typedef itk::BinaryErodeImageFilter<ImageType, ImageType, ElementType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(element);
  filter->SetErodeValue(erodeValue);
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathMO(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType openValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned int                  lines,
        unsigned int                  thickness,
        bool                          includeCenter)
{
  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef typename itk::FlatStructuringElement<ImageType::ImageDimension> ElementType;
  ElementType                                                             element =
    iMathGetFlatStructuringElement<ImageDimension>(shape, radius, radiusIsParametric, lines, thickness, includeCenter);

  typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType, ElementType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetKernel(element);
  filter->SetForegroundValue(openValue);
  filter->SetBackgroundValue(0);
  filter->Update();

  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathNormalize(typename ImageType::Pointer image) /*1*/
{
  if (image->GetNumberOfComponentsPerPixel() != 1)
  {
    // NOPE
  }

  typedef typename ImageType::PixelType PixelType;

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> NormFilterType;
  typename NormFilterType::Pointer                               normFilter = NormFilterType::New();
  normFilter->SetInput(image);
  normFilter->SetOutputMinimum(itk::NumericTraits<PixelType>::ZeroValue());
  normFilter->SetOutputMaximum(itk::NumericTraits<PixelType>::OneValue());
  normFilter->Update();

  return normFilter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathPad(typename ImageType::Pointer image1, int padvalue) /*1*/
{
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typename ImageType::SizeType                         size = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::PointType                        origin2 = image1->GetOrigin();
  typename ImageType::SizeType                         newsize = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::RegionType                       newregion;
  // determine new image size
  for (unsigned int i = 0; i < ImageType::ImageDimension; i++)
  {
    float dimsz = (float)size[i];
    newsize[i] = (unsigned int)(dimsz + padvalue * 2);
  }
  newregion.SetSize(newsize);
  newregion.SetIndex(image1->GetLargestPossibleRegion().GetIndex());
  typename ImageType::Pointer padimage =
    AllocImage<ImageType>(newregion, image1->GetSpacing(), origin2, image1->GetDirection(), 0);

  typename ImageType::IndexType index;
  typename ImageType::IndexType index2;
  if (padvalue > 0)
  {
    index.Fill(0);
    index2.Fill((unsigned int)fabs(static_cast<float>(padvalue)));
  }
  else
  {
    index2.Fill(0);
    index.Fill((unsigned int)fabs(static_cast<float>(padvalue)));
  }

  typename ImageType::PointType point1, pointpad;
  image1->TransformIndexToPhysicalPoint(index, point1);
  padimage->TransformIndexToPhysicalPoint(index2, pointpad);

  for (unsigned int i = 0; i < ImageType::ImageDimension; i++)
  {
    origin2[i] += (point1[i] - pointpad[i]);
  }

  padimage->SetOrigin(origin2);

  Iterator iter(image1, image1->GetLargestPossibleRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    typename ImageType::IndexType oindex = iter.GetIndex();
    typename ImageType::IndexType padindex = iter.GetIndex();

    bool isinside = true;
    for (unsigned int i = 0; i < ImageType::ImageDimension; i++)
    {
      float shifted = ((float)oindex[i] + padvalue);
      if (shifted < 0 || shifted > newsize[i] - 1)
      {
        isinside = false;
      }
      //      if (shifted < 0) shifted=0;
      // padindex[i]=
    }
    if (isinside)
    {
      for (unsigned int i = 0; i < ImageType::ImageDimension; i++)
      {
        float shifted = ((float)oindex[i] + padvalue);
        padindex[i] = (unsigned int)shifted;
      }
      padimage->SetPixel(padindex, iter.Get());
    }
  }
  return padimage;
}


template <typename ImageType>
typename ImageType::Pointer
iMathPeronaMalik(typename ImageType::Pointer image,
                 unsigned long               nIterations, /*1*/
                 double                      conductance)
{
  if (image->GetNumberOfComponentsPerPixel() != 1)
  {
    // NOPE
  }

  typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> FilterType;
  typedef typename FilterType::TimeStepType                                  TimeStepType;

  // Select time step size.
  TimeStepType spacingsize = 0;
  for (unsigned int d = 0; d < ImageType::ImageDimension; d++)
  {
    TimeStepType sp = image->GetSpacing()[d];
    spacingsize += sp * sp;
  }
  spacingsize = sqrt(spacingsize);

  // FIXME - cite reason for this step
  double       dimPlusOne = ImageType::ImageDimension + 1;
  TimeStepType mytimestep = spacingsize / std::pow(2.0, dimPlusOne);
  TimeStepType reftimestep = 0.4 / std::pow(2.0, dimPlusOne);
  if (mytimestep > reftimestep)
  {
    mytimestep = reftimestep;
  }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetConductanceParameter(conductance); // might need to change this
  filter->SetNumberOfIterations(nIterations);
  filter->SetTimeStep(mytimestep);

  filter->Update();
  return filter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathPropagateLabelsThroughMask(typename ImageType::Pointer speedimage, /*1*/
                                typename ImageType::Pointer labimage,
                                double                      stoppingValue,
                                unsigned int                propagationMethod)
{

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef itk::FastMarchingThresholdStoppingCriterion<ImageType, ImageType> CriterionType;

  typedef typename CriterionType::Pointer                        CriterionPointer;
  typedef itk::FastMarchingImageFilterBase<ImageType, ImageType> FastMarchingFilterType;
  typedef typename FastMarchingFilterType::LabelImageType        LabelImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>  ThresholderType;
  typedef itk::LabelContourImageFilter<ImageType, ImageType>     ContourFilterType;
  typedef typename FastMarchingFilterType::NodePairContainerType NodeContainer;
  typedef typename FastMarchingFilterType::NodePairType          NodePairType;
  typedef itk::CastImageFilter<ImageType, ImageType>             CastFilterType;

  typename ImageType::Pointer fastimage = ImageType::New();
  fastimage->SetRegions(speedimage->GetLargestPossibleRegion());
  fastimage->SetSpacing(speedimage->GetSpacing());
  fastimage->SetOrigin(speedimage->GetOrigin());
  fastimage->SetDirection(speedimage->GetDirection());
  fastimage->Allocate();
  fastimage->FillBuffer(1.e9);

  /*
  typename ImageType::Pointer outlabimage = ImageType::New();
  outlabimage->SetRegions( speedimage->GetLargestPossibleRegion() );
  outlabimage->SetSpacing( speedimage->GetSpacing() );
  outlabimage->SetOrigin( speedimage->GetOrigin() );
  outlabimage->SetDirection( speedimage->GetDirection() );
  outlabimage->AllocateInitialized();
  */

  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput(labimage);
  caster->Update();
  typename ImageType::Pointer outlabimage = caster->GetOutput();


  // FIXME - why is thresh set to 0.5?
  double   maxlabel = 0;
  double   thresh = 0.5;
  Iterator vfIter2(labimage, labimage->GetLargestPossibleRegion());
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    bool   isinside = true;
    double speedval = speedimage->GetPixel(vfIter2.GetIndex());
    double labval = labimage->GetPixel(vfIter2.GetIndex());
    if (speedval < thresh)
    {
      isinside = false;
    }
    if (isinside)
    {
      if (labval > maxlabel)
      {
        maxlabel = labval;
      }
    }
  }

  CriterionPointer criterion = CriterionType::New();
  criterion->SetThreshold(stoppingValue);

  typename FastMarchingFilterType::Pointer fastMarching;
  for (unsigned int lab = 1; lab <= (unsigned int)maxlabel; lab++)
  {

    // Create binary mask for each label
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput(labimage);
    thresholder->SetLowerThreshold(lab);
    thresholder->SetUpperThreshold(lab);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);

    // Get pixels on border of the label mask
    typename ContourFilterType::Pointer contour = ContourFilterType::New();
    contour->SetInput(thresholder->GetOutput());
    contour->FullyConnectedOff();
    contour->SetBackgroundValue(itk::NumericTraits<typename LabelImageType::PixelType>::ZeroValue());
    contour->Update();
    typename ImageType::Pointer contourimage = contour->GetOutput();


    fastMarching = FastMarchingFilterType::New();
    fastMarching->SetInput(speedimage);
    fastMarching->SetStoppingCriterion(criterion);

    if (propagationMethod == 1) // Strict
    {
      // std::cout << " strict " << std::endl;
      fastMarching->SetTopologyCheck(itk::FastMarchingTraitsEnums::TopologyCheck::Strict);
    }
    if (propagationMethod == 2) // No handles
    {
      // std::cout << " no handles " << std::endl;
      fastMarching->SetTopologyCheck(itk::FastMarchingTraitsEnums::TopologyCheck::NoHandles);
    }

    typename NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();

    typename NodeContainer::Pointer alivePoints = NodeContainer::New();
    alivePoints->Initialize();

    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      double labval = labimage->GetPixel(vfIter2.GetIndex());
      double contourval = contourimage->GetPixel(vfIter2.GetIndex());
      if (((unsigned int)contourval == 1) && ((unsigned int)labval == lab))
      {
        seeds->push_back(NodePairType(vfIter2.GetIndex(), 0.0));
      }
      if (((unsigned int)contourval == 0) && ((unsigned int)labval == lab))
      {
        alivePoints->push_back(NodePairType(vfIter2.GetIndex(), 0.0));
      }
    }
    fastMarching->SetTrialPoints(seeds);
    fastMarching->SetAlivePoints(alivePoints);
    fastMarching->Update();

    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      bool   isinside = true;
      double speedval = speedimage->GetPixel(vfIter2.GetIndex());
      double labval = labimage->GetPixel(vfIter2.GetIndex());
      if (speedval < thresh)
      {
        isinside = false;
      }


      if (isinside && itk::Math::FloatAlmostEqual(labval, itk::NumericTraits<double>::ZeroValue()))
      {
        double fmarrivaltime = fastMarching->GetOutput()->GetPixel(vfIter2.GetIndex());
        double mmm = fastimage->GetPixel(vfIter2.GetIndex());
        if (fmarrivaltime < mmm)
        {
          fastimage->SetPixel(vfIter2.GetIndex(), fmarrivaltime);
          outlabimage->SetPixel(vfIter2.GetIndex(), lab);
        }
      }
      else if (!isinside)
      {
        outlabimage->SetPixel(vfIter2.GetIndex(), 0);
      }
    }
  }

  return outlabimage;
}

template <typename ImageType>
typename ImageType::Pointer
iMathSharpen(typename ImageType::Pointer image) /*1*/
{
  if (image->GetNumberOfComponentsPerPixel() != 1)
  {
    // NOPE
  }

  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer                                      sharpenFilter = FilterType::New();
  sharpenFilter->SetInput(image);
  sharpenFilter->Update();

  return sharpenFilter->GetOutput();
}

template <typename ImageType>
typename ImageType::Pointer
iMathTruncateIntensity(typename ImageType::Pointer                                           image,
                       double                                                                lowerQ,
                       double                                                                upperQ,
                       int                                                                   nBins, /*1*/
                       typename itk::Image<unsigned int, ImageType::ImageDimension>::Pointer mask)
{

  typedef typename ImageType::PixelType                    PixelType;
  typedef unsigned int                                     LabelType;
  typedef itk::Image<LabelType, ImageType::ImageDimension> MaskType;

  if (mask.IsNull())
  {
    typedef itk::BinaryThresholdImageFilter<ImageType, MaskType> ThresholdFilterType;
    typename ThresholdFilterType::Pointer                        thresh = ThresholdFilterType::New();
    thresh->SetInput(image);
    thresh->SetLowerThreshold(itk::NumericTraits<PixelType>::NonpositiveMin());
    thresh->SetUpperThreshold(itk::NumericTraits<PixelType>::max());
    thresh->SetInsideValue(1);
    thresh->SetOutsideValue(0);
    thresh->Update();
    mask = thresh->GetOutput();
  }
  typedef itk::LabelStatisticsImageFilter<ImageType, MaskType> HistogramFilterType;
  typename HistogramFilterType::Pointer                        stats = HistogramFilterType::New();

  stats->SetInput(image);
  stats->SetLabelInput(mask);
  stats->Update();
  PixelType minValue = stats->GetMinimum(1);
  PixelType maxValue = stats->GetMaximum(1);

  // Hack increment by delta
  if (itk::Math::FloatAlmostEqual(minValue, itk::NumericTraits<PixelType>::ZeroValue()))
  {
    minValue = minValue + static_cast<PixelType>(1e-6);
  }
  if (itk::Math::FloatAlmostEqual(minValue, itk::NumericTraits<PixelType>::ZeroValue()))
  {
    minValue++;
  }

  typedef typename HistogramFilterType::HistogramPointer HistogramPointer;
  typename HistogramFilterType::Pointer                  stats2 = HistogramFilterType::New();
  stats2->SetInput(image);
  stats2->SetLabelInput(mask);
  stats2->UseHistogramsOn();
  stats2->SetHistogramParameters(nBins, minValue, maxValue);
  stats2->Update();
  HistogramPointer histogram = stats2->GetHistogram(1);

  PixelType lowerQuantile = histogram->Quantile(0, lowerQ);
  PixelType upperQuantile = histogram->Quantile(0, upperQ);

  typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> WindowFilterType;
  typename WindowFilterType::Pointer                               windowFilter = WindowFilterType::New();
  windowFilter->SetInput(image);
  windowFilter->SetWindowMinimum(lowerQuantile);
  windowFilter->SetOutputMinimum(lowerQuantile);
  windowFilter->SetWindowMaximum(upperQuantile);
  windowFilter->SetOutputMaximum(upperQuantile);
  windowFilter->Update();

  return windowFilter->GetOutput();
}

//
// shape (1=ball, 2=box, 3=cross, 4=annulus, 5=polygon)

} // namespace ants
