/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWASPSegmentationImageFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWASPSegmentationImageFilter_h
#define __itkWASPSegmentationImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkArray.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkFixedArray.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include <vector>

namespace itk
{
/** \class ApocritaSegmentationImageFilter
 * \brief Apocrita:  A Priori Classification with Registration Initialized
 *  Template Assistance
 *
 * This filter provides an Expectation-Maximization framework for statistical
 * segmentation where the intensity profile of each class is modeled as a
 * Gaussian (a gaussian mixture model (GMM)) and spatial smoothness is
 * enforced by an MRF prior.
 *
 * Initial labeling can be performed by otsu thresholding, kmeans clustering,
 * a set of user-specified prior probability images, or a prior label image.
 * If specified, the latter two initialization options are also used as
 * priors in the MRF update step.
 *
 * The assumed labeling is such that classes are assigned consecutive
 * indices 1, 2, 3, etc.  Label 0 is reserved for the background when a
 * mask is specified.
 *
 */

template <class TInputImage, class TMaskImage
            = Image<unsigned char, ::itk::GetImageDimension<TInputImage>::ImageDimension>,
          class TClassifiedImage = TMaskImage>
class ITK_EXPORT WASPSegmentationImageFilter :
  public ImageToImageFilter<TInputImage, TClassifiedImage>
{
public:
  /** Standard class typdedefs. */
  typedef WASPSegmentationImageFilter                       Self;
  typedef ImageToImageFilter<TInputImage, TClassifiedImage> Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WASPSegmentationImageFilter, ImageToImageFilter );

  /** Dimension of the images. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );
  itkStaticConstMacro( ClassifiedImageDimension, unsigned int,
                       TClassifiedImage::ImageDimension );
  itkStaticConstMacro( MaskImageDimension, unsigned int,
                       TMaskImage::ImageDimension );

  /** Typedef support of input types. */
  typedef TInputImage                             ImageType;
  typedef typename ImageType::PixelType           PixelType;
  typedef TMaskImage                              MaskImageType;
  typedef typename MaskImageType::PixelType       MaskLabelType;
  typedef TClassifiedImage                        ClassifiedImageType;
  typedef typename ClassifiedImageType::PixelType LabelType;

  /** Some convenient typedefs. */
  typedef float RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( ImageDimension )>         RealImageType;
  typedef Array<double> ParametersType;
  typedef FixedArray<unsigned,
                     itkGetStaticConstMacro( ImageDimension )>         ArrayType;

  /** B-spline fitting typedefs */
  typedef Vector<RealType, 1> ScalarType;
  typedef Image<ScalarType,
                itkGetStaticConstMacro( ImageDimension )>         ScalarImageType;
  typedef PointSet<ScalarType,
                   itkGetStaticConstMacro( ImageDimension )>         PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
  <PointSetType, ScalarImageType>                   BSplineFilterType;
  typedef typename
  BSplineFilterType::PointDataImageType             ControlPointLatticeType;

  enum InitializationStrategyType
        { KMeans, Otsu, PriorProbabilityImages, PriorLabelImage };

  /** ivars Set/Get functionality */

  itkSetClampMacro( NumberOfClasses, unsigned int, 2,
                    NumericTraits<LabelType>::max() );
  itkGetConstMacro( NumberOfClasses, unsigned int );

  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

  itkSetMacro( ConvergenceThreshold, RealType );
  itkGetConstMacro( ConvergenceThreshold, RealType );

  itkGetConstMacro( CurrentConvergenceMeasurement, RealType );

  itkGetConstMacro( ElapsedIterations, unsigned int );

  itkSetMacro( MRFSmoothingFactor, RealType );
  itkGetConstMacro( MRFSmoothingFactor, RealType );

  itkSetMacro( MRFSigmoidAlpha, RealType );
  itkGetConstMacro( MRFSigmoidAlpha, RealType );

  itkSetMacro( MRFSigmoidBeta, RealType );
  itkGetConstMacro( MRFSigmoidBeta, RealType );

  itkSetMacro( MRFRadius, ArrayType );
  itkGetConstMacro( MRFRadius, ArrayType );

  itkSetMacro( InitializationStrategy, InitializationStrategyType );
  itkGetConstMacro( InitializationStrategy, InitializationStrategyType );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( NumberOfLevels, ArrayType );
  itkGetConstMacro( NumberOfLevels, ArrayType );

  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  void SetMaskImage( const MaskImageType * mask );

  const MaskImageType * GetMaskImage() const;

  typename RealImageType::Pointer GetPosteriorImage( unsigned int n )
  {
    if( n > this->m_PosteriorImages.size() )
      {
      return NULL;
      }
    else
      {
      return this->m_PosteriorImages[n];
      }
  }

  itkSetClampMacro( PriorProbabilityWeighting, RealType, 0, 1 );
  itkGetConstMacro( PriorProbabilityWeighting, RealType );

  void AmassDistancePriors();

  void SetPriorLabelSigmas( std::vector<RealType> s )
  {
    this->m_PriorLabelSigmas = s;
    this->Modified();
  }

  void SetPriorProbabilityImage(unsigned int whichClass, const RealImageType * prior );

  const RealImageType * GetPriorProbabilityImage( unsigned int i ) const;

  void SetPriorLabelImage( const ClassifiedImageType * prior );

  const ClassifiedImageType * GetPriorLabelImage() const;

  typename RealImageType::Pointer CalculatePosteriorProbabilityImage( unsigned int, bool calcdist );

  typename RealImageType::Pointer CalculateSmoothIntensityImageFromPriorProbabilityImage( unsigned int );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SameDimensionCheck1,
                   ( Concept::SameDimension<ImageDimension,
                                            ClassifiedImageDimension> ) );
  itkConceptMacro( SameDimensionCheck2,
                   ( Concept::SameDimension<ImageDimension,
                                            MaskImageDimension> ) );
  /** End concept checking */
#endif
protected:
  WASPSegmentationImageFilter();
  ~WASPSegmentationImageFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  WASPSegmentationImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);              // purposely not implemented

  void NormalizePriorProbabilityImages();

  void GenerateInitialClassLabeling();

  void GenerateInitialClassLabelingWithOtsuThresholding();

  void GenerateInitialClassLabelingWithKMeansClustering();

  void GenerateInitialClassLabelingWithPriorProbabilityImages();

//  void SpatiallySmoothClassLabelingWithMRF();

/** this should be used when we cant afford to store the posteriors */
  RealType RecursiveUpdateClassParametersAndLabeling();

/** this should be used for cases when we can store all images */
  RealType StraightUpdateClassParametersAndLabeling();

  unsigned int m_NumberOfClasses;
  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumNumberOfIterations;
  RealType     m_CurrentConvergenceMeasurement;
  RealType     m_ConvergenceThreshold;

  MaskLabelType m_MaskLabel;

  std::vector<ParametersType> m_CurrentClassParameters;
  InitializationStrategyType  m_InitializationStrategy;

  ArrayType m_MRFRadius;
  RealType  m_MRFSmoothingFactor;
  RealType  m_MRFSigmoidAlpha;
  RealType  m_MRFSigmoidBeta;

  unsigned int m_SplineOrder;
  ArrayType    m_NumberOfLevels;
  ArrayType    m_NumberOfControlPoints;
  std::vector<typename
              ControlPointLatticeType::Pointer> m_ControlPointLattices;
  std::vector<RealType>                         m_PriorLabelSigmas;

  RealType m_PriorProbabilityWeighting;
  typename RealImageType::Pointer  m_SumProbabilityImage;
  std::vector<typename RealImageType::Pointer> m_PosteriorImages;
  std::vector<typename RealImageType::Pointer> m_DistanceImages;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWASPSegmentationImageFilter.txx"
#endif

#endif
