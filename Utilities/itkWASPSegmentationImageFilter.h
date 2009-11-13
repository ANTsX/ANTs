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
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

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

  typename RealImageType::Pointer GetLargestComponent(  typename RealImageType::Pointer image1 )
  {
    typedef float                                        PixelType;
    typedef itk::Vector<float, ImageDimension>           VectorType;
    typedef itk::Image<VectorType, ImageDimension>       FieldType;
    typedef itk::Image<PixelType, ImageDimension>        ImageType;
    typedef typename ImageType::IndexType                IndexType;
    typedef typename ImageType::SizeType                 SizeType;
    typedef typename ImageType::SpacingType              SpacingType;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

    unsigned long smallest = 50;
    typename ImageType::SpacingType spacing = image1->GetSpacing();
    float volumeelement = 1.0;
    for( unsigned int i = 0;  i < spacing.Size(); i++ )
      {
      volumeelement *= spacing[i];
      }

    typedef float InternalPixelType;
    //  typedef unsigned long PixelType;
    //  typedef Image<PixelType,ImageDimension>  labelimagetype;
    typedef itk::Image<unsigned long, ImageDimension>                          labelimagetype;
    typedef RealImageType                                                      InternalImageType;
    typedef RealImageType                                                      OutputImageType;
    typedef itk::BinaryThresholdImageFilter<InternalImageType, labelimagetype> ThresholdFilterType;
    typedef itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype> FilterType;
    typedef itk::RelabelComponentImageFilter<labelimagetype, ImageType>        RelabelType;

    typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
    typename FilterType::Pointer filter = FilterType::New();
    typename RelabelType::Pointer relabel = RelabelType::New();

    //  InternalPixelType threshold_low, threshold_hi;

    threshold->SetInput(image1);
    threshold->SetInsideValue(1);
    threshold->SetOutsideValue(0);
    threshold->SetLowerThreshold(0.25);
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
      std::cerr << "Relabel: exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      }

    //  WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
    //  return 0;
    typename ImageType::Pointer Clusters = relabel->GetOutput();
    // typename ImageType::Pointer Clusters=relabel->GetOutput();
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    Iterator vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );

    float maximum = relabel->GetNumberOfObjects();
    std::cout << " #ob " << maximum << std::endl;
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
        float vox = image1->GetPixel(vfIter.GetIndex() );
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

    std::cout << " max float size "
              <<  (maximgval
         * volumeelement) << " long-size: " << (unsigned long) (maximgval * volumeelement)  << std::endl;
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      if( Clusters->GetPixel( vfIter.GetIndex() ) >= maximgval )
        {
        image1->SetPixel( vfIter.GetIndex(), 1);
        }
      else
        {
        image1->SetPixel( vfIter.GetIndex(), 0);
        }
      }

    return image1;
  }

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
