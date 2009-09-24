/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabelOverlapMeasuresImageFilter_h
#define __itkLabelOverlapMeasuresImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkFastMutexLock.h"
#include "itkNumericTraits.h"

#include "itk_hash_map.h"

namespace itk
{
/** \class LabelOverlapMeasuresImageFilter
 * \brief Computes overlap measures between the set same set of labels of
 * pixels of two images.  Background is assumed to be 0.
 *
 * \sa LabelOverlapMeasuresImageFilter
 *
 * \ingroup MultiThreaded
 */
template <class TLabelImage>
class ITK_EXPORT LabelOverlapMeasuresImageFilter :
  public ImageToImageFilter<TLabelImage, TLabelImage>
{
public:
  /** Standard Self typedef */
  typedef LabelOverlapMeasuresImageFilter              Self;
  typedef ImageToImageFilter<TLabelImage, TLabelImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( LabelOverlapMeasuresImageFilter, ImageToImageFilter );

  /** Image related typedefs. */
  typedef TLabelImage                        LabelImageType;
  typedef typename TLabelImage::Pointer      LabelImagePointer;
  typedef typename TLabelImage::ConstPointer LabelImageConstPointer;

  typedef typename TLabelImage::RegionType RegionType;
  typedef typename TLabelImage::SizeType   SizeType;
  typedef typename TLabelImage::IndexType  IndexType;

  typedef typename TLabelImage::PixelType LabelType;

  /** Type to use form computations. */
  typedef typename NumericTraits<LabelType>::RealType RealType;

  /** \class LabelLabelOverlapMeasuress
   * \brief Metrics stored per label */
  class LabelSetMeasures
  {
public:
    // default constructor
    LabelSetMeasures()
    {
      m_VolumeSource = 0;
      m_VolumeTarget = 0;
      m_VolumeUnion = 0;
      m_VolumeIntersection = 0;
      m_VolumeSourceComplement = 0;
      m_VolumeTargetComplement = 0;

      m_SurfaceSource = 0;
      m_SurfaceTarget = 0;
      m_SurfaceUnion = 0;
      m_SurfaceIntersection = 0;
      m_SurfaceSourceComplement = 0;
      m_SurfaceTargetComplement = 0;

      m_ContourMeanDistance = 0.0;
      m_DirectedContourMeanDistance = 0.0;
      m_HausdorffDistance = 0.0;
      m_DirectedHausdorffDistance = 0.0;
    }

    // added for completeness
    LabelSetMeasures & operator=( const LabelSetMeasures& l )
    {
      m_VolumeSource = l.m_VolumeSource;
      m_VolumeTarget = l.m_VolumeTarget;
      m_VolumeUnion = l.m_VolumeUnion;
      m_VolumeIntersection = l.m_VolumeIntersection;
      m_VolumeSourceComplement = l.m_VolumeSourceComplement;
      m_VolumeTargetComplement = l.m_VolumeTargetComplement;

      m_SurfaceSource = l.m_SurfaceSource;
      m_SurfaceTarget = l.m_SurfaceTarget;
      m_SurfaceUnion = l.m_SurfaceUnion;
      m_SurfaceIntersection = l.m_SurfaceIntersection;
      m_SurfaceSourceComplement = l.m_SurfaceSourceComplement;
      m_SurfaceTargetComplement = l.m_SurfaceTargetComplement;

      m_ContourMeanDistance = l.m_ContourMeanDistance;
      m_DirectedContourMeanDistance = l.m_DirectedContourMeanDistance;
      m_HausdorffDistance = l.m_HausdorffDistance;
      m_DirectedHausdorffDistance = l.m_DirectedHausdorffDistance;
    }

    unsigned long m_VolumeSource;
    unsigned long m_VolumeTarget;
    unsigned long m_VolumeUnion;
    unsigned long m_VolumeIntersection;
    unsigned long m_VolumeSourceComplement;
    unsigned long m_VolumeTargetComplement;

    unsigned long m_SurfaceSource;
    unsigned long m_SurfaceTarget;
    unsigned long m_SurfaceUnion;
    unsigned long m_SurfaceIntersection;
    unsigned long m_SurfaceSourceComplement;
    unsigned long m_SurfaceTargetComplement;

    RealType m_ContourMeanDistance;
    RealType m_DirectedContourMeanDistance;
    RealType m_HausdorffDistance;
    RealType m_DirectedHausdorffDistance;
  };

  /** Type of the map used to store data per label */
  typedef hash_map<LabelType, LabelSetMeasures> MapType;
  typedef typename MapType::iterator            MapIterator;
  typedef typename MapType::const_iterator      MapConstIterator;

  /** Image related typedefs. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TLabelImage::ImageDimension );

  /** Set the source image. */
  void SetSourceImage( const LabelImageType * image )
  {
    this->SetNthInput( 0, const_cast<LabelImageType *>( image ) );
  }

  /** Set the target image. */
  void SetTargetImage( const LabelImageType * image )
  {
    this->SetNthInput( 1, const_cast<LabelImageType *>( image ) );
  }

  /** Get the source image. */
  const LabelImageType * GetSourceImage( void )
  {
    return this->GetInput( 0 );
  }

  /** Get the target image. */
  const LabelImageType * GetTargetImage( void )
  {
    return this->GetInput( 1 );
  }

  /** Get the label set measures */
  MapType GetLabelSetMeasures()
  {
    return this->m_LabelSetMeasures;
  }

  /**
   * Volumetric overlap measures
   */
  /** measures over all labels */
  RealType GetVolumeTotalOverlap();

  RealType GetVolumeUnionOverlap();

  RealType GetVolumeMeanOverlap();

  RealType GetVolumeFalseNegativeError();

  RealType GetVolumeFalsePositiveError();

  /** measures over individual labels */
  RealType GetVolumeTotalOverlap( LabelType );
  RealType GetVolumeUnionOverlap( LabelType );
  RealType GetVolumeMeanOverlap( LabelType );
  RealType GetVolumeFalseNegativeError( LabelType );
  RealType GetVolumeFalsePositiveError( LabelType );
  /** alternative names */
  RealType GetVolumeJaccardCoefficient()
  {
    return this->GetVolumeUnionOverlap();
  }

  RealType GetVolumeJaccardCoefficient( LabelType label )
  {
    return this->GetVolumeUnionOverlap( label );
  }

  RealType GetVolumeDiceCoefficient()
  {
    return this->GetVolumeMeanOverlap();
  }

  RealType GetVolumeDiceCoefficient( LabelType label )
  {
    return this->GetVolumeMeanOverlap( label );
  }

  /**
   * Surface overlap measures
   */
  /** measures over all labels */
  RealType GetSurfaceTotalOverlap();

  RealType GetSurfaceUnionOverlap();

  RealType GetSurfaceMeanOverlap();

  RealType GetSurfaceFalseNegativeError();

  RealType GetSurfaceFalsePositiveError();

  /** measures over individual labels */
  RealType GetSurfaceTotalOverlap( LabelType );
  RealType GetSurfaceUnionOverlap( LabelType );
  RealType GetSurfaceMeanOverlap( LabelType );
  RealType GetSurfaceFalseNegativeError( LabelType );
  RealType GetSurfaceFalsePositiveError( LabelType );
  /** alternative names */
  RealType GetSurfaceJaccardCoefficient()
  {
    return this->GetSurfaceUnionOverlap();
  }

  RealType GetSurfaceJaccardCoefficient( LabelType label )
  {
    return this->GetSurfaceUnionOverlap( label );
  }

  RealType GetSurfaceDiceCoefficient()
  {
    return this->GetSurfaceMeanOverlap();
  }

  RealType GetSurfaceDiceCoefficient( LabelType label )
  {
    return this->GetSurfaceMeanOverlap( label );
  }

  /**
   * Other measures
   */
  RealType GetVolumeSimilarity();

  RealType GetVolumeSimilarity( LabelType );
  RealType GetHausdorffDistance( LabelType );
  RealType GetDirectedHausdorffDistance( LabelType );
  RealType GetContourMeanDistance( LabelType );
  RealType GetDirectedContourMeanDistance( LabelType );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( Input1HasNumericTraitsCheck,
                   ( Concept::HasNumericTraits<LabelType> ) );
  /** End concept checking */
#endif
protected:
  LabelOverlapMeasuresImageFilter();
  ~LabelOverlapMeasuresImageFilter()
  {
  };
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /**
   * Pass the input through unmodified. Do this by setting the output to the
   * source this by setting the output to the source image in the
   * AllocateOutputs() method.
   */
  void AllocateOutputs();

  void BeforeThreadedGenerateData();

  void AfterThreadedGenerateData();

  /** Multi-thread version GenerateData. */
  void ThreadedGenerateData( const RegionType &, int );

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion( DataObject *data );

private:
  LabelOverlapMeasuresImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                  // purposely not implemented

  LabelImagePointer m_SourceSurfaceImage;
  LabelImagePointer m_TargetSurfaceImage;

  std::vector<MapType> m_LabelSetMeasuresPerThread;
  MapType              m_LabelSetMeasures;

  SimpleFastMutexLock m_Mutex;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelOverlapMeasuresImageFilter.txx"
#endif

#endif
