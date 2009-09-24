/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLabelOverlapMeasuresImageFilter_txx
#define _itkLabelOverlapMeasuresImageFilter_txx

#include "itkLabelOverlapMeasuresImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkContourDirectedMeanDistanceImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkDirectedHausdorffDistanceImageFilter.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkLabelContourImageFilter.h"

#include "itkImageRegionConstIterator.h"

namespace itk
{
#if defined(__GNUC__) && (__GNUC__ <= 2) // NOTE: This class needs a mutex for gnu 2.95
/** Used for mutex locking */
#define LOCK_HASHMAP this->m_Mutex.Lock()
#define UNLOCK_HASHMAP this->m_Mutex.Unlock()
#else
#define LOCK_HASHMAP
#define UNLOCK_HASHMAP
#endif

template <class TLabelImage>
LabelOverlapMeasuresImageFilter<TLabelImage>
::LabelOverlapMeasuresImageFilter()
{
  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if( this->GetSourceImage() )
    {
    LabelImagePointer source = const_cast
      <LabelImageType *>( this->GetSourceImage() );
    source->SetRequestedRegionToLargestPossibleRegion();
    }
  if( this->GetTargetImage() )
    {
    LabelImagePointer target = const_cast
      <LabelImageType *>( this->GetTargetImage() );
    target->SetRequestedRegionToLargestPossibleRegion();
    }
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::EnlargeOutputRequestedRegion( DataObject *data )
{
  Superclass::EnlargeOutputRequestedRegion( data );
  data->SetRequestedRegionToLargestPossibleRegion();
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::AllocateOutputs()
{
  // Pass the source through as the output
  LabelImagePointer image =
    const_cast<TLabelImage *>( this->GetSourceImage() );

  this->SetNthOutput( 0, image );

  // Nothing that needs to be allocated for the remaining outputs
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::BeforeThreadedGenerateData()
{
  int numberOfThreads = this->GetNumberOfThreads();

  // Resize the thread temporaries
  this->m_LabelSetMeasuresPerThread.resize( numberOfThreads );
  // Initialize the temporaries
  for( int n = 0; n < numberOfThreads; n++ )
    {
    this->m_LabelSetMeasuresPerThread[n].clear();
    }

  // Initialize the final map
  this->m_LabelSetMeasures.clear();

  // Create the source/target surface images
  typedef LabelContourImageFilter
  <LabelImageType, LabelImageType> ContourFilterType;

  typename ContourFilterType::Pointer sourceContours = ContourFilterType::New();
  sourceContours->SetInput( this->GetSourceImage() );
  sourceContours->SetFullyConnected( true );
  sourceContours->SetBackgroundValue( NumericTraits<LabelType>::Zero );
  sourceContours->Update();
  this->m_SourceSurfaceImage = sourceContours->GetOutput();

  typename ContourFilterType::Pointer targetContours = ContourFilterType::New();
  targetContours->SetInput( this->GetTargetImage() );
  targetContours->SetFullyConnected( true );
  targetContours->SetBackgroundValue( NumericTraits<LabelType>::Zero );
  targetContours->Update();
  this->m_TargetSurfaceImage = targetContours->GetOutput();
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::AfterThreadedGenerateData()
{
  // Run through the map for each thread and accumulate the set measures.
  for( int n = 0; n < this->GetNumberOfThreads(); n++ )
    {
    // iterate over the map for this thread
    for( MapConstIterator threadIt = this->m_LabelSetMeasuresPerThread[n].begin();
         threadIt != this->m_LabelSetMeasuresPerThread[n].end();
         ++threadIt )
      {
      // does this label exist in the cumulative stucture yet?
      MapIterator mapIt = this->m_LabelSetMeasures.find( ( *threadIt ).first );
      if( mapIt == this->m_LabelSetMeasures.end() )
        {
        // create a new entry
        typedef typename MapType::value_type MapValueType;
        mapIt = this->m_LabelSetMeasures.insert( MapValueType(
                                                   (*threadIt).first, LabelSetMeasures() ) ).first;
        }

      // accumulate the information from this thread
      (*mapIt).second.m_VolumeSource += (*threadIt).second.m_VolumeSource;
      (*mapIt).second.m_VolumeTarget += (*threadIt).second.m_VolumeTarget;
      (*mapIt).second.m_VolumeUnion += (*threadIt).second.m_VolumeUnion;
      (*mapIt).second.m_VolumeIntersection +=
        (*threadIt).second.m_VolumeIntersection;
      (*mapIt).second.m_VolumeSourceComplement +=
        (*threadIt).second.m_VolumeSourceComplement;
      (*mapIt).second.m_VolumeTargetComplement +=
        (*threadIt).second.m_VolumeTargetComplement;

      (*mapIt).second.m_SurfaceSource += (*threadIt).second.m_SurfaceSource;
      (*mapIt).second.m_SurfaceTarget += (*threadIt).second.m_SurfaceTarget;
      (*mapIt).second.m_SurfaceUnion += (*threadIt).second.m_SurfaceUnion;
      (*mapIt).second.m_SurfaceIntersection +=
        (*threadIt).second.m_SurfaceIntersection;
      (*mapIt).second.m_SurfaceSourceComplement +=
        (*threadIt).second.m_SurfaceSourceComplement;
      (*mapIt).second.m_SurfaceTargetComplement +=
        (*threadIt).second.m_SurfaceTargetComplement;
      } // end of thread map iterator loop
    }   // end of thread loop
  // compute the remainder of the overlap measures
  for( MapIterator mapIt = m_LabelSetMeasures.begin();
       mapIt != m_LabelSetMeasures.end(); ++mapIt )
    {
    LabelType label = (*mapIt).first;

    typedef BinaryThresholdImageFilter<LabelImageType, LabelImageType>
    ThresholderType;
    typename ThresholderType::Pointer sourceThresholder =
      ThresholderType::New();
    sourceThresholder->SetInput( this->GetSourceImage() );
    sourceThresholder->SetLowerThreshold( label );
    sourceThresholder->SetUpperThreshold( label );
    sourceThresholder->SetInsideValue( NumericTraits<LabelType>::One );
    sourceThresholder->SetOutsideValue( NumericTraits<LabelType>::Zero );
    sourceThresholder->Update();

    typename ThresholderType::Pointer targetThresholder =
      ThresholderType::New();
    targetThresholder->SetInput( this->GetTargetImage() );
    targetThresholder->SetLowerThreshold( label );
    targetThresholder->SetUpperThreshold( label );
    targetThresholder->SetInsideValue( NumericTraits<LabelType>::One );
    targetThresholder->SetOutsideValue( NumericTraits<LabelType>::Zero );
    targetThresholder->Update();

    // Contour distances

    typedef ContourDirectedMeanDistanceImageFilter
    <LabelImageType, LabelImageType> DirectedContourFilterType;
    typename DirectedContourFilterType::Pointer directedContourFilter =
      DirectedContourFilterType::New();
    directedContourFilter->SetInput1( sourceThresholder->GetOutput() );
    directedContourFilter->SetInput2( targetThresholder->GetOutput() );
    directedContourFilter->Update();

    typedef ContourMeanDistanceImageFilter
    <LabelImageType, LabelImageType> ContourFilterType;
    typename ContourFilterType::Pointer contourFilter =
      ContourFilterType::New();
    contourFilter->SetInput1( sourceThresholder->GetOutput() );
    contourFilter->SetInput2( targetThresholder->GetOutput() );
    contourFilter->Update();

    (*mapIt).second.m_DirectedContourMeanDistance =
      directedContourFilter->GetContourDirectedMeanDistance();
    (*mapIt).second.m_ContourMeanDistance = contourFilter->GetMeanDistance();

    // Hausdorff distances

    typedef DirectedHausdorffDistanceImageFilter
    <LabelImageType, LabelImageType> DirectedHausdorffFilterType;
    typename DirectedHausdorffFilterType::Pointer directedHausdorffFilter =
      DirectedHausdorffFilterType::New();
    directedHausdorffFilter->SetInput1( sourceThresholder->GetOutput() );
    directedHausdorffFilter->SetInput2( targetThresholder->GetOutput() );
    directedHausdorffFilter->Update();

    typedef HausdorffDistanceImageFilter
    <LabelImageType, LabelImageType> HausdorffFilterType;
    typename HausdorffFilterType::Pointer hausdorffFilter =
      HausdorffFilterType::New();
    hausdorffFilter->SetInput1( sourceThresholder->GetOutput() );
    hausdorffFilter->SetInput2( targetThresholder->GetOutput() );
    hausdorffFilter->Update();

    (*mapIt).second.m_DirectedHausdorffDistance =
      directedHausdorffFilter->GetDirectedHausdorffDistance();
    (*mapIt).second.m_HausdorffDistance =
      hausdorffFilter->GetHausdorffDistance();
    }
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::ThreadedGenerateData( const RegionType& outputRegionForThread,
                        int threadId )
{
  ImageRegionConstIterator<LabelImageType> ItS( this->GetSourceImage(),
                                                outputRegionForThread );
  ImageRegionConstIterator<LabelImageType> ItT( this->GetTargetImage(),
                                                outputRegionForThread );

  // support progress methods/callbacks
  ProgressReporter progress( this, threadId,
                             2 * outputRegionForThread.GetNumberOfPixels() );
  for( ItS.GoToBegin(), ItT.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItT )
    {
    LabelType sourceLabel = ItS.Get();
    LabelType targetLabel = ItT.Get();

    // is the label already in this thread?
    MapIterator mapItS =
      this->m_LabelSetMeasuresPerThread[threadId].find( sourceLabel );
    MapIterator mapItT =
      this->m_LabelSetMeasuresPerThread[threadId].find( targetLabel );

    if( mapItS == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItS = this->m_LabelSetMeasuresPerThread[threadId].insert(
          MapValueType( sourceLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    if( mapItT == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItT = this->m_LabelSetMeasuresPerThread[threadId].insert(
          MapValueType( targetLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    (*mapItS).second.m_VolumeSource++;
    (*mapItT).second.m_VolumeTarget++;

    if( sourceLabel == targetLabel )
      {
      (*mapItS).second.m_VolumeIntersection++;
      (*mapItS).second.m_VolumeUnion++;
      }
    else
      {
      (*mapItS).second.m_VolumeUnion++;
      (*mapItT).second.m_VolumeUnion++;

      (*mapItS).second.m_VolumeSourceComplement++;
      (*mapItT).second.m_VolumeTargetComplement++;
      }

    progress.CompletedPixel();
    }

  ImageRegionConstIterator<LabelImageType> Its( this->m_SourceSurfaceImage,
                                                outputRegionForThread );
  ImageRegionConstIterator<LabelImageType> Itt( this->m_TargetSurfaceImage,
                                                outputRegionForThread );
  for( Its.GoToBegin(), Itt.GoToBegin(); !Its.IsAtEnd(); ++Its, ++Itt )
    {
    LabelType sourceLabel = Its.Get();
    LabelType targetLabel = Itt.Get();

    // is the label already in this thread?
    MapIterator mapIts =
      this->m_LabelSetMeasuresPerThread[threadId].find( sourceLabel );
    MapIterator mapItt =
      this->m_LabelSetMeasuresPerThread[threadId].find( targetLabel );

    if( mapIts == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapIts = this->m_LabelSetMeasuresPerThread[threadId].insert(
          MapValueType( sourceLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    if( mapItt == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItt = this->m_LabelSetMeasuresPerThread[threadId].insert(
          MapValueType( targetLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    (*mapIts).second.m_SurfaceSource++;
    (*mapItt).second.m_SurfaceTarget++;

    if( sourceLabel == targetLabel )
      {
      (*mapIts).second.m_SurfaceIntersection++;
      (*mapIts).second.m_SurfaceUnion++;
      }
    else
      {
      (*mapIts).second.m_SurfaceUnion++;
      (*mapItt).second.m_SurfaceUnion++;

      (*mapIts).second.m_SurfaceSourceComplement++;
      (*mapItt).second.m_SurfaceTargetComplement++;
      }

    progress.CompletedPixel();
    }
}

/**
 * Volume measures
 */
template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeTotalOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_VolumeIntersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_VolumeTarget );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeTotalOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_VolumeIntersection )
    / static_cast<RealType>( (*mapIt).second.m_VolumeTarget );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeUnionOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_VolumeIntersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_VolumeUnion );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeUnionOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_VolumeIntersection )
    / static_cast<RealType>( (*mapIt).second.m_VolumeUnion );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeMeanOverlap()
{
  RealType uo = this->GetVolumeUnionOverlap();

  return 2.0 * uo / ( 1.0 + uo );
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeMeanOverlap( LabelType label )
{
  RealType uo = this->GetVolumeUnionOverlap( label );

  return 2.0 * uo / ( 1.0 + uo );
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeFalseNegativeError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_VolumeTargetComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_VolumeTarget );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeFalseNegativeError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_VolumeTargetComplement )
    / static_cast<RealType>( (*mapIt).second.m_VolumeTarget );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeFalsePositiveError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_VolumeSourceComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_VolumeSource );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeFalsePositiveError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_VolumeSourceComplement )
    / static_cast<RealType>( (*mapIt).second.m_VolumeSource );
  return value;
}

/**
 * Surface measures
 */
template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceTotalOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_SurfaceIntersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_SurfaceTarget );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceTotalOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_SurfaceIntersection )
    / static_cast<RealType>( (*mapIt).second.m_SurfaceTarget );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceUnionOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_SurfaceIntersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_SurfaceUnion );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceUnionOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_SurfaceIntersection )
    / static_cast<RealType>( (*mapIt).second.m_SurfaceUnion );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceMeanOverlap()
{
  RealType uo = this->GetSurfaceUnionOverlap();

  return 2.0 * uo / ( 1.0 + uo );
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceMeanOverlap( LabelType label )
{
  RealType uo = this->GetSurfaceUnionOverlap( label );

  return 2.0 * uo / ( 1.0 + uo );
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceFalseNegativeError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_SurfaceTargetComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_SurfaceTarget );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceFalseNegativeError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_SurfaceTargetComplement )
    / static_cast<RealType>( (*mapIt).second.m_SurfaceTarget );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceFalsePositiveError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_SurfaceSourceComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_SurfaceSource );
    }
  return numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetSurfaceFalsePositiveError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_SurfaceSourceComplement )
    / static_cast<RealType>( (*mapIt).second.m_SurfaceSource );
  return value;
}

/**
 * Other measures
 */
template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeSimilarity()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;

  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
       mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += ( static_cast<RealType>( (*mapIt).second.m_VolumeSource )
                   - static_cast<RealType>( (*mapIt).second.m_VolumeTarget ) );
    denominator += ( static_cast<RealType>( (*mapIt).second.m_VolumeSource )
                     + static_cast<RealType>( (*mapIt).second.m_VolumeTarget ) );
    }
  return 2.0 * numerator / denominator;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeSimilarity( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value = 2.0
    * ( static_cast<RealType>( (*mapIt).second.m_VolumeSource )
        - static_cast<RealType>( (*mapIt).second.m_VolumeTarget ) )
    / ( static_cast<RealType>( (*mapIt).second.m_VolumeSource )
        + static_cast<RealType>( (*mapIt).second.m_VolumeTarget ) );
  return value;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetHausdorffDistance( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  return (*mapIt).second.m_HausdorffDistance;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetDirectedHausdorffDistance( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  return (*mapIt).second.m_DirectedHausdorffDistance;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetContourMeanDistance( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  return (*mapIt).second.m_ContourMeanDistance;
}

template <class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetDirectedContourMeanDistance( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );

  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  return (*mapIt).second.m_DirectedContourMeanDistance;
}

template <class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk
#endif
