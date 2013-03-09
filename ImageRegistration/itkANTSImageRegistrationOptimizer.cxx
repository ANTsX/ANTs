/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkANTSImageRegistrationOptimizer.cxx,v $
  Language:  C++
  Date:      $Date: 2009/04/22 01:00:16 $
  Version:   $Revision: 1.47 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkANTSImageRegistrationOptimizer_hxx_
#define _itkANTSImageRegistrationOptimizer_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif
#include "antsAllocImage.h"
#include "itkVectorParameterizedNeighborhoodOperatorImageFilter.h"
#include "itkANTSImageRegistrationOptimizer.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkVectorGaussianInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "vnl/vnl_math.h"
#include "ANTS_affine_registration2.h"
#include "itkWarpImageMultiTransformFilter.h"
// #include "itkVectorImageFileWriter.h"

namespace itk
{
template <unsigned int TDimension, class TReal>
ANTSImageRegistrationOptimizer<TDimension, TReal>
::ANTSImageRegistrationOptimizer()
{
  this->m_DisplacementField = NULL;
  this->m_InverseDisplacementField = NULL;
  this->m_AffineTransform = NULL;
  itk::TransformFactory<TransformType>::RegisterTransform();
  itk::TransformFactory<itk::ANTSAffine3DTransform<TReal> >::RegisterTransform();
  itk::TransformFactory<itk::ANTSCenteredAffine2DTransform<TReal> >::RegisterTransform();
  this->m_FixedPointSet = NULL;
  this->m_MovingPointSet = NULL;

  this->m_UseMulti = true;
  this->m_UseROI = false;
  this->m_MaskImage = NULL;
  this->m_ReferenceSpaceImage = NULL;
  this->m_Debug = false;

  this->m_ScaleFactor = 1.0;
  this->m_SubsamplingFactors.SetSize( 0 );
  this->m_GaussianSmoothingSigmas.SetSize( 0 );

  this->m_SyNF = NULL;
  this->m_SyNFInv = NULL;
  this->m_SyNM = NULL;
  this->m_SyNMInv = NULL;
  this->m_Parser = NULL;
  this->m_GaussianTruncation = 256;
  this->m_TimeVaryingVelocity = NULL;
  this->m_LastTimeVaryingVelocity = NULL;
  this->m_LastTimeVaryingUpdate = NULL;
  this->m_DeltaTime = 0.1;
  this->m_SyNType = 0;
  this->m_UseNN = false;
  this->m_UseBSplineInterpolation = false;
  this->m_VelocityFieldInterpolator = VelocityFieldInterpolatorType::New();
  this->m_HitImage = NULL;
  this->m_ThickImage = NULL;
  this->m_SyNFullTime = 0;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::ImagePointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SubsampleImage( ImagePointer image, RealType /* scalingFactor */, typename ImageType::PointType outputOrigin,
                  typename ImageType::DirectionType outputDirection,
                  AffineTransformPointer aff )
{
  typename ImageType::SpacingType outputSpacing = this->m_CurrentDomainSpacing;
  typename ImageType::RegionType::SizeType outputSize = this->m_CurrentDomainSize;

//    RealType minimumSpacing = inputSpacing.GetVnlVector().min_value();
//    RealType maximumSpacing = inputSpacing.GetVnlVector().max_value();

  typedef ResampleImageFilter<ImageType, ImageType> ResamplerType;
  typename ResamplerType::Pointer resampler = ResamplerType::New();
  typedef LinearInterpolateImageFunction<ImageType, TComp> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image );
  resampler->SetInterpolator( interpolator );
  typedef itk::IdentityTransform<TComp, TDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer transform = IdentityTransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );
  if( aff )
    {
//      ::ants::antscout << " Setting Aff to " << this->m_AffineTransform << std::endl;
    resampler->SetTransform( aff );
    }
  resampler->SetInput( image );
  resampler->SetOutputSpacing( outputSpacing );
  resampler->SetOutputOrigin( outputOrigin );
  resampler->SetOutputDirection( outputDirection );
  resampler->SetSize( outputSize );
  resampler->Update();

  ImagePointer outimage = resampler->GetOutput();

  if( this->m_UseROI )
    {
    outimage = this->MakeSubImage(outimage);
//       WriteImage<ImageType>(outimage,"temps.hdr");
    // warp with affine & deformable
    }

  return outimage;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::CopyDisplacementField(  DisplacementFieldPointer input  )
{
  DisplacementFieldPointer output =
    AllocImage<DisplacementFieldType>(input);

  typedef ImageRegionIterator<DisplacementFieldType> Iterator;
  Iterator inIter( input, input->GetBufferedRegion() );
  Iterator outIter( output, output->GetBufferedRegion() );
  inIter.GoToBegin();
  outIter.GoToBegin();
  for( ; !inIter.IsAtEnd(); ++inIter, ++outIter )
    {
    outIter.Set( inIter.Get() );
    }

  return output;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SmoothDisplacementFieldGauss(DisplacementFieldPointer field, TReal sig, bool /* useparamimage */, unsigned int lodim)
{
  if( this->m_Debug )
    {
    ::ants::antscout << " enter gauss smooth " <<  sig  << std::endl;
    }
  if( sig <= 0 )
    {
    return;
    }
  if( !field )
    {
    ::ants::antscout << " No Field in gauss Smoother " << std::endl; return;
    }
  DisplacementFieldPointer tempField =
    AllocImage<DisplacementFieldType>(field);

  typedef typename DisplacementFieldType::PixelType    DispVectorType;
  typedef typename DispVectorType::ValueType           ScalarType;
  typedef GaussianOperator<ScalarType, ImageDimension> OperatorType;
  // typedef VectorNeighborhoodOperatorImageFilter< DisplacementFieldType,    DisplacementFieldType> SmootherType;
  typedef VectorParameterizedNeighborhoodOperatorImageFilter<
      DisplacementFieldType,
      DisplacementFieldType, ImageType> SmootherType;

  OperatorType * oper = new OperatorType;
  typename SmootherType::Pointer smoother = SmootherType::New();

  typedef typename DisplacementFieldType::PixelContainerPointer
    PixelContainerPointer;
  PixelContainerPointer swapPtr;

  // graft the output field onto the mini-pipeline
  smoother->GraftOutput( tempField );
  for( unsigned int j = 0; j < lodim; j++ )
    {
    // smooth along this dimension
    oper->SetDirection( j );
    TReal sigt = sig;
    oper->SetVariance( sigt );
    oper->SetMaximumError(0.001 );
    oper->SetMaximumKernelWidth( (unsigned int) this->m_GaussianTruncation );
    oper->CreateDirectional();

    // todo: make sure we only smooth within the buffered region
    smoother->SetOperator( *oper );
    smoother->SetInput( field );
    smoother->Update();

    if( j < lodim - 1 )
      {
      // swap the containers
      swapPtr = smoother->GetOutput()->GetPixelContainer();
      smoother->GraftOutput( field );
      field->SetPixelContainer( swapPtr );
      smoother->Modified();
      }
    }

  // graft the output back to this filter
  tempField->SetPixelContainer( field->GetPixelContainer() );

  // make sure boundary does not move
  TReal weight = 1.0;
  if( sig < 0.5 )
    {
    weight = 1.0 - 1.0 * (sig / 0.5);
    }
  TReal weight2 = 1.0 - weight;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
  typename DisplacementFieldType::SizeType size = field->GetLargestPossibleRegion().GetSize();
  Iterator outIter( field, field->GetLargestPossibleRegion() );
  for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
    {
    bool onboundary = false;
    typename DisplacementFieldType::IndexType index = outIter.GetIndex();
    for( int i = 0; i < ImageDimension; i++ )
      {
      if( index[i] < 1 || index[i] >= static_cast<int>( size[i] ) - 1 )
        {
        onboundary = true;
        }
      }
    if( onboundary )
      {
      DispVectorType vec;
      vec.Fill(0.0);
      outIter.Set(vec);
      }
    else
      {
      // field=this->CopyDisplacementField(
      DispVectorType svec = smoother->GetOutput()->GetPixel(index);
      outIter.Set( svec * weight + outIter.Get() * weight2);
      }
    }

  if( this->m_Debug )
    {
    ::ants::antscout << " done gauss smooth " << std::endl;
    }

  delete oper;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SmoothVelocityGauss(TimeVaryingVelocityFieldPointer field, TReal sig, unsigned int lodim)
{
  if( sig <= 0 )
    {
    return;
    }
  if( !field )
    {
    ::ants::antscout << " No Field in gauss Smoother " << std::endl; return;
    }
  TimeVaryingVelocityFieldPointer tempField =
    AllocImage<TimeVaryingVelocityFieldType>(field);

  typedef typename TimeVaryingVelocityFieldType::PixelType TVVFVectorType;
  typedef typename TVVFVectorType::ValueType               ScalarType;
  typedef GaussianOperator<ScalarType, ImageDimension + 1> OperatorType;
  typedef VectorNeighborhoodOperatorImageFilter<TimeVaryingVelocityFieldType,
                                                TimeVaryingVelocityFieldType> SmootherType;

  OperatorType * oper = new OperatorType;
  typename SmootherType::Pointer smoother = SmootherType::New();

  typedef typename TimeVaryingVelocityFieldType::PixelContainerPointer
    PixelContainerPointer;
  PixelContainerPointer swapPtr;

  // graft the output field onto the mini-pipeline
  smoother->GraftOutput( tempField );
  for( unsigned int j = 0; j < lodim; j++ )
    {
    // smooth along this dimension
    oper->SetDirection( j );
    oper->SetVariance( sig );
    oper->SetMaximumError(0.001 );
    oper->SetMaximumKernelWidth( (unsigned int) this->m_GaussianTruncation );
    oper->CreateDirectional();

    // todo: make sure we only smooth within the buffered region
    smoother->SetOperator( *oper );
    smoother->SetInput( field );
    smoother->Update();

    if( j < lodim - 1 )
      {
      // swap the containers
      swapPtr = smoother->GetOutput()->GetPixelContainer();
      smoother->GraftOutput( field );
      field->SetPixelContainer( swapPtr );
      smoother->Modified();
      }
    }

  // graft the output back to this filter
  tempField->SetPixelContainer( field->GetPixelContainer() );

  // make sure boundary does not move
  TReal weight = 1.0;
  if( sig < 0.5 )
    {
    weight = 1.0 - 1.0 * (sig / 0.5);
    }
  TReal weight2 = 1.0 - weight;
  typedef itk::ImageRegionIteratorWithIndex<TimeVaryingVelocityFieldType> Iterator;
  typename TimeVaryingVelocityFieldType::SizeType size = field->GetLargestPossibleRegion().GetSize();
  Iterator outIter( field, field->GetLargestPossibleRegion() );
  for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
    {
    bool onboundary = false;
    typename TimeVaryingVelocityFieldType::IndexType index = outIter.GetIndex();
    for( int i = 0; i < ImageDimension; i++ )
      {
      if( index[i] < 1 || index[i] >= static_cast<int>( size[i] ) - 1 )
        {
        onboundary = true;
        }
      }
    if( onboundary )
      {
      TVVFVectorType vec;
      vec.Fill(0.0);
      outIter.Set(vec);
      }
    else
      {
      // field=this->CopyDisplacementField(
      TVVFVectorType svec = smoother->GetOutput()->GetPixel(index);
      outIter.Set( svec * weight + outIter.Get() * weight2);
      }
    }

  if( this->m_Debug )
    {
    ::ants::antscout << " done gauss smooth " << std::endl;
    }

  delete oper;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SmoothDisplacementFieldBSpline( DisplacementFieldPointer field, ArrayType meshsize,
                                  unsigned int splineorder, unsigned int numberoflevels )
{
  if( this->m_Debug )
    {
    ::ants::antscout << " enter bspline smooth " << std::endl;
    }
  if( !field )
    {
    ::ants::antscout << " No Field in bspline Smoother " << std::endl; return;
    }

  if( splineorder <= 0 )
    {
    return;
    }

  typename BSplineFilterType::ArrayType numberofcontrolpoints;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if( meshsize[d] <= 0 )
      {
      return;
      }

    numberofcontrolpoints[d] = static_cast<unsigned int>( meshsize[d] ) + splineorder;
    }
  VectorType zeroVector;
  zeroVector.Fill( 0.0 );

//  typedef VectorImageFileWriter<DisplacementFieldType, ImageType>
//    DisplacementFieldWriterType;
//  typename DisplacementFieldWriterType::Pointer writer = DisplacementFieldWriterType::New();
//  writer->SetInput( field );
//  writer->SetFileName( "field.nii.gz" );
//  writer->Update();
//  std::exception();

  typename ImageType::DirectionType originalDirection = field->GetDirection();
  typename ImageType::DirectionType identityDirection;
  identityDirection.SetIdentity();
  field->SetDirection( identityDirection );

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetInput( field );
  bspliner->SetNumberOfLevels( numberoflevels );
  bspliner->SetSplineOrder( splineorder );
  bspliner->SetNumberOfControlPoints( numberofcontrolpoints );
  bspliner->SetIgnorePixelValue( zeroVector );
  bspliner->Update();

  field->SetDirection( originalDirection );

  // make sure boundary does not move
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;

  Iterator bIter( bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion() );
  Iterator outIter( field, field->GetLargestPossibleRegion() );
  for( outIter.GoToBegin(), bIter.GoToBegin();
       !outIter.IsAtEnd(); ++outIter, ++bIter )
    {
//    bool onboundary=false;
//    typename DisplacementFieldType::IndexType index = outIter.GetIndex();
//    for( int i = 0; i < ImageDimension; i++ )
//      {
//      if ( index[i] < 1 || index[i] >= static_cast<int>( size[i] )-1 )
//        onboundary = true;
//      }
//    if (onboundary)
//      {
//      VectorType vec;
//      vec.Fill(0.0);
//      outIter.Set(vec);
//      }
//    else
//      {
    outIter.Set( bIter.Get() );
//      }
    }

  if( this->m_Debug )
    {
    ::ants::antscout << " done bspline smooth " << std::endl;
    }
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::ComposeDiffs(DisplacementFieldPointer fieldtowarpby, DisplacementFieldPointer field,
               DisplacementFieldPointer fieldout,
               TReal timesign)
{
  typedef Point<TReal, itkGetStaticConstMacro(ImageDimension)> VPointType;

//  field->SetSpacing( fieldtowarpby->GetSpacing() );
//  field->SetOrigin( fieldtowarpby->GetOrigin() );
//  field->SetDirection( fieldtowarpby->GetDirection() );

  if( !fieldout )
    {
    VectorType zero;  zero.Fill(0);
    fieldout = AllocImage<DisplacementFieldType>(fieldtowarpby);
    }
  typedef typename DisplacementFieldType::PixelType DispVectorType;

  typedef itk::WarpImageFilter<ImageType, ImageType, DisplacementFieldType> WarperType;
  typedef DisplacementFieldType                                             FieldType;

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef ImageType                                                TRealImageType;

  typedef itk::ImageFileWriter<ImageType> writertype;

  typedef typename DisplacementFieldType::IndexType DispIndexType;

  typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, TReal>   DefaultInterpolatorType;
  typedef itk::VectorGaussianInterpolateImageFunction<DisplacementFieldType, TReal> DefaultInterpolatorType2;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  vinterp->SetInputImage(field);
  //    vinterp->SetParameters(NULL,1);

  VPointType pointIn1;
  VPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  contind;   // married to pointIn2
  VPointType   pointIn3;
  unsigned int ct = 0;
  // iterate through fieldtowarpby finding the points that it maps to via field.
  // then take the difference from the original point and put it in the output field.
  //      ::ants::antscout << " begin iteration " << std::endl;
  FieldIterator m_FieldIter( fieldtowarpby, fieldtowarpby->GetLargestPossibleRegion() );
  for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
    {
    DispIndexType index = m_FieldIter.GetIndex();
    bool          dosample = true;
    //      if (sub && m_TRealImage->GetPixel(index) < 0.5) dosample=false;
    if( dosample )
      {
      fieldtowarpby->TransformIndexToPhysicalPoint( index, pointIn1 );
      DispVectorType disp = m_FieldIter.Get();
      for( int jj = 0; jj < ImageDimension; jj++ )
        {
        pointIn2[jj] = disp[jj] + pointIn1[jj];
        }
      typename DefaultInterpolatorType::OutputType disp2;
      if( vinterp->IsInsideBuffer(pointIn2) )
        {
        disp2 = vinterp->Evaluate( pointIn2 );
        }
      else
        {
        disp2.Fill(0);
        }
      for( int jj = 0; jj < ImageDimension; jj++ )
        {
        pointIn3[jj] = disp2[jj] * timesign + pointIn2[jj];
        }

      DispVectorType out;
      for( int jj = 0; jj < ImageDimension; jj++ )
        {
        out[jj] = pointIn3[jj] - pointIn1[jj];
        }

      fieldout->SetPixel(m_FieldIter.GetIndex(), out);
      ct++;
      } // endif
    }   // end iteration
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::IntegrateConstantVelocity(DisplacementFieldPointer totalField, unsigned int ntimesteps, TReal timestep)
{
  VectorType zero;

  zero.Fill(0);
  DisplacementFieldPointer diffmap =
    AllocImage<DisplacementFieldType>(totalField, zero);
  for( unsigned int nts = 0; nts < ntimesteps; nts++ )
    {
    this->ComposeDiffs(diffmap, totalField, diffmap, timestep);
    }
  return diffmap;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::ComputeUpdateField(DisplacementFieldPointer fixedwarp, DisplacementFieldPointer movingwarp,   PointSetPointer fpoints,
                     PointSetPointer wpoints, DisplacementFieldPointer totalUpdateInvField,
                     bool updateenergy)
{
  ImagePointer mask = NULL;

  if( movingwarp && this->m_MaskImage && !this->m_ComputeThickness )
    {
    mask = this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_MaskImage, NULL, movingwarp, false,
                                     this->m_FixedImageAffineTransform );
    }
  else if( this->m_MaskImage && !this->m_ComputeThickness  )
    {
    mask = this->SubsampleImage( this->m_MaskImage, this->m_ScaleFactor,
                                 this->m_MaskImage->GetOrigin(), this->m_MaskImage->GetDirection(),  NULL);
    }

  if( !fixedwarp )
    {
    ::ants::antscout << " NO F WARP " << std::endl;  fixedwarp = this->m_DisplacementField;
    }
  // if ( !movingwarp) ::ants::antscout<< " NO M WARP " << std::endl;

  // /     ::ants::antscout << " get upd field " << std::endl;
  typename ImageType::SpacingType spacing = fixedwarp->GetSpacing();
  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer updateField;
  DisplacementFieldPointer updateFieldInv;

  DisplacementFieldPointer totalUpdateField =
    AllocImage<DisplacementFieldType>(fixedwarp, zero);

//    bool hadpointsetmetric=false;

  RealType sumWeights = 0.0;
  for( unsigned int n = 0; n < this->m_SimilarityMetrics.size(); n++ )
    {
    sumWeights += this->m_SimilarityMetrics[n]->GetWeightScalar();
    }
  sumWeights = 1;
  for( unsigned int metricCount = 0; metricCount < this->m_SimilarityMetrics.size(); metricCount++ )
    {
    bool ispointsetmetric = false;

    /** build an update field */
    if(  this->m_SimilarityMetrics.size() == 1 )
      {
      updateField = totalUpdateField;
      if( totalUpdateInvField )
        {
        updateFieldInv = totalUpdateInvField;
        }
      }
    else
      {
      updateField = AllocImage<DisplacementFieldType>(fixedwarp, zero);
      if( totalUpdateInvField )
        {
        updateFieldInv = AllocImage<DisplacementFieldType>(fixedwarp, zero);
        }
      }

    /** get the update */
    typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodIteratorType;
    typedef ImageRegionIterator<DisplacementFieldType> UpdateIteratorType;

//        TimeStepType timeStep;
    void *globalData;
//    ::ants::antscout << " B " << std::endl;

    AffineTransformPointer faffinverse = NULL;
    if( this->m_FixedImageAffineTransform )
      {
      faffinverse = AffineTransformType::New();
      this->m_FixedImageAffineTransform->GetInverse(faffinverse);
      }
    AffineTransformPointer affinverse = NULL;
    if( this->m_AffineTransform )
      {
      affinverse = AffineTransformType::New();
      this->m_AffineTransform->GetInverse(affinverse);
      }

// for each metric, warp the assoc. Images
/** We loop Over This To Do MultiVariate */
/** FIXME really should pass an image list and then warp each one in
      turn  then expand the update field to fit size of total
      deformation */
    ImagePointer wmimage = NULL;
    if( fixedwarp )
      {
      wmimage =
        this->WarpMultiTransform(  this->m_ReferenceSpaceImage, this->m_SmoothMovingImages[metricCount],
                                   this->m_AffineTransform,
                                   fixedwarp, false,
                                   NULL );
      }
    else
      {
      wmimage = this->SubsampleImage( this->m_SmoothMovingImages[metricCount], this->m_ScaleFactor,
                                      this->m_SmoothMovingImages[metricCount]->GetOrigin(),
                                      this->m_SmoothMovingImages[metricCount]->GetDirection(),  NULL);
      }

//    ::ants::antscout << " C " << std::endl;
    ImagePointer wfimage = NULL;
    if( movingwarp )
      {
      wfimage =
        this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_SmoothFixedImages[metricCount], NULL, movingwarp,
                                  false,
                                  this->m_FixedImageAffineTransform );
      }
    else
      {
      wfimage = this->SubsampleImage( this->m_SmoothFixedImages[metricCount], this->m_ScaleFactor,
                                      this->m_SmoothFixedImages[metricCount]->GetOrigin(),
                                      this->m_SmoothFixedImages[metricCount]->GetDirection(),  NULL);
      }
    /*
    if (this->m_TimeVaryingVelocity && ! this->m_MaskImage ) {
      std::string outname=this->localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str())+std::string("thick.nii.gz");
      ///      WriteImage<ImageType>(wmimage,outname.c_str());
      outname=this->localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str())+std::string("thick2.nii.gz");
      WriteImage<ImageType>(wfimage,outname.c_str());
    }
    */
    //      std::string
    // outname=this->localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str())+std::string("temp.nii.gz");
    //      WriteImage<ImageType>(wmimage,outname.c_str());
    //      std::string
    // outname2=this->localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str())+std::string("temp2.nii.gz");
    //      WriteImage<ImageType>(wfimage,outname2.c_str());

/** MV Loop END -- Would have to collect update fields then add them
* together somehow -- Would also have to eliminate the similarity
* metric loop within ComputeUpdateField */

    // Get the FiniteDifferenceFunction to use in calculations.
    MetricBaseTypePointer df = this->m_SimilarityMetrics[metricCount]->GetModifiableMetric();
    df->SetFixedImage(wfimage);
    df->SetMovingImage(wmimage);
    if( df->ThisIsAPointSetMetric() )
      {
      ispointsetmetric = true;
      }
    if( fpoints && ispointsetmetric )
      {
      df->SetFixedPointSet(fpoints);
      }
    else if( ispointsetmetric )
      {
      ::ants::antscout << "NO POINTS!! " << std::endl;
      }
    if( wpoints && ispointsetmetric )
      {
      df->SetMovingPointSet(wpoints);
      }
    else if( ispointsetmetric )
      {
      ::ants::antscout << "NO POINTS!! " << std::endl;
      }
    typename ImageType::SizeType  radius = df->GetRadius();
    df->InitializeIteration();
    typename DisplacementFieldType::Pointer output = updateField;
    typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<DisplacementFieldType>
      FaceCalculatorType;
    typedef typename FaceCalculatorType::FaceListType FaceListType;
    FaceCalculatorType faceCalculator;
    FaceListType       faceList = faceCalculator(updateField, updateField->GetLargestPossibleRegion(), radius);
    typename FaceListType::iterator fIt = faceList.begin();
    globalData = df->GetGlobalDataPointer();

    // Process the non-boundary region.
    NeighborhoodIteratorType nD(radius, updateField, *fIt);
    UpdateIteratorType       nU(updateField,  *fIt);
    nD.GoToBegin();
    nU.GoToBegin();
    while( !nD.IsAtEnd() )
      {
      bool  oktosample = true;
      TReal maskprob = 1.0;
      if( mask )
        {
        maskprob = mask->GetPixel( nD.GetIndex() );
        if( maskprob > 1.0 )
          {
          maskprob = 1.0;
          }
        if( maskprob < 0.1 )
          {
          oktosample = false;
          }
        }
      if( oktosample )
        {
        VectorType temp = df->ComputeUpdate(nD, globalData) * maskprob;
        nU.Value() += temp;
        if( totalUpdateInvField )
          {
          typename ImageType::IndexType index = nD.GetIndex();
          temp = df->ComputeUpdateInv(nD, globalData) * maskprob + updateFieldInv->GetPixel(index);
          updateFieldInv->SetPixel(index, temp);
          }       // else nU.Value() -= df->ComputeUpdateInv(nD, globalData)*maskprob;

        ++nD;
        ++nU;
        }
      else
        {
        ++nD;
        ++nU;
        }
      }

    // begin restriction of deformation field
    bool restrict = false;
    for( unsigned int jj = 0; jj < this->m_RestrictDeformation.size();  jj++ )
      {
      TReal temp = this->m_RestrictDeformation[jj];
      if(  fabs( temp - 1 ) > 1.e-5   )
        {
        restrict = true;
        }
      }
    if( restrict && this->m_RestrictDeformation.size() == ImageDimension )
      {
      nU.GoToBegin();
      while( !nU.IsAtEnd() )
        {
        for( unsigned int jj = 0; jj < this->m_RestrictDeformation.size();  jj++ )
          {
          typename ImageType::IndexType index = nU.GetIndex();
          VectorType temp = updateField->GetPixel(index);
          temp[jj] *= this->m_RestrictDeformation[jj];
          updateField->SetPixel(index, temp);
          if( updateFieldInv )
            {
            temp = updateFieldInv->GetPixel(index);
            temp[jj] *= this->m_RestrictDeformation[jj];
            updateFieldInv->SetPixel(index, temp);
            }
          }
        ++nU;
        }
      } // end restrict deformation field

    if( updateenergy )
      {
      this->m_LastEnergy[metricCount] = this->m_Energy[metricCount];
      this->m_Energy[metricCount] = df->GetEnergy();
      // *this->m_SimilarityMetrics[metricCount]->GetWeightScalar()/sumWeights;
      }

    // smooth the fields
    // if (!ispointsetmetric || ImageDimension == 2 ){
    this->SmoothDisplacementField(updateField, true);
    if( updateFieldInv )
      {
      this->SmoothDisplacementField(updateFieldInv, true);
      }
    // /}
    /*
    else // use another strategy -- exact lm? / something like Laplacian
  {
    TReal tmag=0;
    for (unsigned int ff=0; ff<5; ff++)
      {
        tmag=0;
        this->SmoothDisplacementField(updateField,true);
        if (updateFieldInv) this->SmoothDisplacementField(updateFieldInv,true);
        nD.GoToBegin();
        nU.GoToBegin();
        while( !nD.IsAtEnd() )
      {
        typename ImageType::IndexType index=nD.GetIndex();
        bool oktosample=true;
        TReal maskprob=1.0;
        if (mask)
          {
            maskprob=mask->GetPixel( nD.GetIndex() );
            if (maskprob > 1.0) maskprob=1.0;
            if ( maskprob < 0.1) oktosample=false;
          }
        VectorType F1;
        F1.Fill(0);
        VectorType F2;
        F2.Fill(0);
        if ( oktosample )
          {
            F1 = df->ComputeUpdate(nD, globalData)*maskprob;
            if (totalUpdateInvField)
          {
            F2 = df->ComputeUpdateInv(nD, globalData)*maskprob;
          }
            ++nD;
            ++nU;
          }
        else
          {
            ++nD;
            ++nU;
          }

        // compute mags of F1 and F2 -- if large enough, reset them
        TReal f1mag=0,f2mag=0,umag=0;
        for (unsigned int dim=0; dim<ImageDimension; dim++)
          {
            f1mag+=F1[dim]/spacing[dim]*F1[dim]/spacing[dim];
            f2mag+=F2[dim]/spacing[dim]*F2[dim]/spacing[dim];
            umag+=updateField->GetPixel(index)[dim]/spacing[dim]*updateField->GetPixel(index)[dim]/spacing[dim];
          }
        f1mag=sqrt(f1mag); f2mag=sqrt(f2mag); umag=sqrt(umag);
        if ( f1mag > 0.05 ) updateField->SetPixel(index,F1);
        if ( f2mag > 0.05 ) updateFieldInv->SetPixel(index,F2);
        tmag+=umag;
      }
        //           ::ants::antscout << " total mag " << tmag << std::endl;
      }
    //smooth the total field
    this->SmoothDisplacementField(updateField,true);
    if (updateFieldInv) this->SmoothDisplacementField(updateFieldInv,true);
  }

    */
    // normalize update field then add to total field
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    Iterator      dIter(totalUpdateField, totalUpdateField->GetLargestPossibleRegion() );
    TReal         mag = 0.0;
    TReal         max = 0.0;
    unsigned long ct = 0;
    TReal         total = 0;
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType vec = updateField->GetPixel(index);
      mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
        }
      mag = sqrt(mag);
//            if (mag > 0. ) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << std::endl;
      if( mag > max )
        {
        max = mag;
        }
      ct++;
      total += mag;
      //    ::ants::antscout << " mag " << mag << std::endl;
      }
    if( this->m_Debug )
      {
      ::ants::antscout << "PRE MAX " << max << std::endl;
      }
    TReal max2 = 0;
    if( max <= 0 )
      {
      max = 1;
      }
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType vec = updateField->GetPixel(index);
      vec = vec / max;
      mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
        }
      mag = sqrt(mag);
      if( mag >  max2 )
        {
        max2 = mag;
        }
//            if (mag > 0.95) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << " ind " << index
// <<
// std::endl;
      /** FIXME need weights between metrics */

      RealType normalizedWeight
        = this->m_SimilarityMetrics[metricCount]->GetWeightScalar() / sumWeights;
//            RealType weight = this->m_SimilarityMetrics[metricCount]->GetWeightImage()->GetPixel( diter.GetIndex() );
      if( ispointsetmetric )
        {
        VectorType intensityupdate = dIter.Get();
        VectorType lmupdate = vec;
        TReal      lmag = 0;
        for( unsigned int li = 0; li < ImageDimension; li++ )
          {
          lmag += (lmupdate[li] / spacing[li]) * (lmupdate[li] / spacing[li]);
          }
        lmag = sqrt(lmag);
        TReal modi = 1;
        if( lmag > 1 )
          {
          modi = 0;
          }
        else
          {
          modi = 1.0 - lmag;
          }
        TReal      iwt = 1 * modi;
        TReal      lmwt = normalizedWeight;
        VectorType totalv = intensityupdate * iwt + lmupdate * lmwt;
        dIter.Set(totalv);
        }
      else
        {
        dIter.Set(dIter.Get() + vec * normalizedWeight);
        }
      }

    if( totalUpdateInvField )
      {
      Iterator invDIter(totalUpdateInvField, totalUpdateInvField->GetLargestPossibleRegion() );
      mag = 0.0;
      max = 0.0;
      ct = 0;
      total = 0;
      for( invDIter.GoToBegin(); !invDIter.IsAtEnd(); ++invDIter )
        {
        typename ImageType::IndexType index = invDIter.GetIndex();
        VectorType vec = updateFieldInv->GetPixel(index);
        mag = 0;
        for( unsigned int jj = 0; jj < ImageDimension; jj++ )
          {
          mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
          }
        mag = sqrt(mag);
//            if (mag > 0. ) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << std::endl;
        if( mag > max )
          {
          max = mag;
          }
        ct++;
        total += mag;
        //    ::ants::antscout << " mag " << mag << std::endl;
        }
      if( this->m_Debug )
        {
        ::ants::antscout << "PRE MAX " << max << std::endl;
        }
      max2 = 0;
      if( max <= 0 )
        {
        max = 1;
        }
      for( invDIter.GoToBegin(); !invDIter.IsAtEnd(); ++invDIter )
        {
        typename ImageType::IndexType index = invDIter.GetIndex();
        VectorType vec = updateFieldInv->GetPixel(index);
        vec = vec / max;
        mag = 0;
        for( unsigned int jj = 0; jj < ImageDimension; jj++ )
          {
          mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
          }
        mag = sqrt(mag);
        if( mag >  max2 )
          {
          max2 = mag;
          }
//            if (mag > 0.95) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << " ind " << index
// <<
// std::endl;
        /** FIXME need weights between metrics */

        RealType normalizedWeight
          = this->m_SimilarityMetrics[metricCount]->GetWeightScalar() / sumWeights;
//            RealType weight = this->m_SimilarityMetrics[metricCount]->GetWeightImage()->GetPixel( diter.GetIndex() );

        if( ispointsetmetric )
          {
          VectorType intensityupdate = invDIter.Get();
          VectorType lmupdate = vec;
          TReal      lmag = 0;
          for( unsigned int li = 0; li < ImageDimension; li++ )
            {
            lmag += (lmupdate[li] / spacing[li]) * (lmupdate[li] / spacing[li]);
            }
          lmag = sqrt(lmag);
          TReal modi = 1;
          if( lmag > 1 )
            {
            modi = 0;
            }
          else
            {
            modi = 1.0 - lmag;
            }
          TReal      iwt = 1 * modi;
          TReal      lmwt = normalizedWeight;
          VectorType totalv = intensityupdate * iwt + lmupdate * lmwt;
          invDIter.Set(totalv);
          }
        else
          {
          invDIter.Set(invDIter.Get() + vec * normalizedWeight);
          }
        }
      }
    if( this->m_Debug )
      {
      ::ants::antscout << "PO MAX " << max2 << " sz" << totalUpdateField->GetLargestPossibleRegion().GetSize()
                       << std::endl;
      }
    }

//    this->SmoothDisplacementField( totalUpdateField,true);
//    if (totalUpdateInvField) this->SmoothDisplacementField( totalUpdateInvField,true);

  return totalUpdateField;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::ComputeUpdateFieldAlternatingMin(DisplacementFieldPointer fixedwarp, DisplacementFieldPointer movingwarp,
                                   PointSetPointer fpoints, PointSetPointer wpoints,
                                   DisplacementFieldPointer totalUpdateInvField,
                                   bool updateenergy)
{
  ImagePointer mask = NULL;

  if( movingwarp && this->m_MaskImage )
    {
    mask = this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_MaskImage, NULL, movingwarp, false,
                                     this->m_FixedImageAffineTransform );
    }
  else if( this->m_MaskImage )
    {
    mask = this->SubsampleImage( this->m_MaskImage, this->m_ScaleFactor,
                                 this->m_MaskImage->GetOrigin(), this->m_MaskImage->GetDirection(),  NULL);
    }

  if( !fixedwarp )
    {
    ::ants::antscout << " NO F WARP " << std::endl;  fixedwarp = this->m_DisplacementField;
    }
  // if ( !movingwarp) ::ants::antscout<< " NO M WARP " << std::endl;

  // /     ::ants::antscout << " get upd field " << std::endl;
  typename ImageType::SpacingType spacing = fixedwarp->GetSpacing();
  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer updateField;
  DisplacementFieldPointer updateFieldInv;
  DisplacementFieldPointer totalUpdateField =
    AllocImage<DisplacementFieldType>(fixedwarp, zero);

  RealType sumWeights = 0.0;
  for( unsigned int n = 0; n < this->m_SimilarityMetrics.size(); n++ )
    {
    sumWeights += this->m_SimilarityMetrics[n]->GetWeightScalar();
    }
  sumWeights = 1;

  // for ( unsigned int metricCount = 0; metricCount < this->m_SimilarityMetrics.size(); metricCount++ )
  // {
  // MetricBaseTypePointer df = this->m_SimilarityMetrics[metricCount]->GetMetric();
  //  if (df->ThisIsAPointSetMetric()) hadpointsetmetric=true;
  // }

  // for ( unsigned int metricCount = 0; metricCount < this->m_SimilarityMetrics.size(); metricCount++ )
  unsigned int metricCount = this->m_CurrentIteration %  this->m_SimilarityMetrics.size();
    {
    bool ispointsetmetric = false;

    /** build an update field */
    if( true )  // this->m_SimilarityMetrics.size() == 1 )
      {
      updateField = totalUpdateField;
      if( totalUpdateInvField )
        {
        updateFieldInv = totalUpdateInvField;
        }
      }
    else
      {
      updateField = AllocImage<DisplacementFieldType>(fixedwarp, zero);
      if( totalUpdateInvField )
        {
        updateFieldInv = AllocImage<DisplacementFieldType>(fixedwarp, zero);
        }
      }

    /** get the update */
    typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodIteratorType;
    typedef ImageRegionIterator<DisplacementFieldType> UpdateIteratorType;

//        TimeStepType timeStep;
    void *globalData;
//    ::ants::antscout << " B " << std::endl;

// for each metric, warp the assoc. Images
/** We loop Over This To Do MultiVariate */
/** FIXME really should pass an image list and then warp each one in
      turn  then expand the update field to fit size of total
      deformation */
    ImagePointer wmimage = NULL;
    if( fixedwarp )
      {
      wmimage =
        this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_SmoothMovingImages[metricCount],
                                  this->m_AffineTransform,
                                  fixedwarp, false,
                                  this->m_FixedImageAffineTransform );
      }
    else
      {
      wmimage = this->SubsampleImage( this->m_SmoothMovingImages[metricCount], this->m_ScaleFactor,
                                      this->m_SmoothMovingImages[metricCount]->GetOrigin(),
                                      this->m_SmoothMovingImages[metricCount]->GetDirection(),  NULL);
      }

//    ::ants::antscout << " C " << std::endl;
    ImagePointer wfimage = NULL;
    if( movingwarp )
      {
      wfimage =
        this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_SmoothFixedImages[metricCount], NULL, movingwarp,
                                  false,
                                  this->m_FixedImageAffineTransform );
      }
    else
      {
      wfimage = this->SubsampleImage( this->m_SmoothFixedImages[metricCount], this->m_ScaleFactor,
                                      this->m_SmoothFixedImages[metricCount]->GetOrigin(),
                                      this->m_SmoothFixedImages[metricCount]->GetDirection(),  NULL);
      }

//    ::ants::antscout << " D " << std::endl;

/** MV Loop END -- Would have to collect update fields then add them
* together somehow -- Would also have to eliminate the similarity
* metric loop within ComputeUpdateField */

    // Get the FiniteDifferenceFunction to use in calculations.
    MetricBaseTypePointer df = this->m_SimilarityMetrics[metricCount]->GetMetric();
    df->SetFixedImage(wfimage);
    df->SetMovingImage(wmimage);
    if( df->ThisIsAPointSetMetric() )
      {
      ispointsetmetric = true;
      }
    if( fpoints && ispointsetmetric )
      {
      df->SetFixedPointSet(fpoints);
      }
    else if( ispointsetmetric )
      {
      ::ants::antscout << "NO POINTS!! " << std::endl;
      }
    if( wpoints && ispointsetmetric )
      {
      df->SetMovingPointSet(wpoints);
      }
    else if( ispointsetmetric )
      {
      ::ants::antscout << "NO POINTS!! " << std::endl;
      }
    typename ImageType::SizeType  radius = df->GetRadius();
    df->InitializeIteration();
    typename DisplacementFieldType::Pointer output = updateField;
    typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<DisplacementFieldType>
      FaceCalculatorType;
    typedef typename FaceCalculatorType::FaceListType FaceListType;
    FaceCalculatorType faceCalculator;
    FaceListType       faceList = faceCalculator(updateField, updateField->GetLargestPossibleRegion(), radius);
    typename FaceListType::iterator fIt = faceList.begin();
    globalData = df->GetGlobalDataPointer();

    // Process the non-boundary region.
    NeighborhoodIteratorType nD(radius, updateField, *fIt);
    UpdateIteratorType       nU(updateField,  *fIt);
    nD.GoToBegin();
    nU.GoToBegin();
    while( !nD.IsAtEnd() )
      {
      bool  oktosample = true;
      TReal maskprob = 1.0;
      if( mask )
        {
        maskprob = mask->GetPixel( nD.GetIndex() );
        if( maskprob > 1.0 )
          {
          maskprob = 1.0;
          }
        if( maskprob < 0.1 )
          {
          oktosample = false;
          }
        }
      if( oktosample )
        {
        nU.Value() += df->ComputeUpdate(nD, globalData) * maskprob;
        if( totalUpdateInvField )
          {
          typename ImageType::IndexType index = nD.GetIndex();
          VectorType temp = df->ComputeUpdateInv(nD, globalData) * maskprob;
          updateFieldInv->SetPixel(index, temp);
          }        // else nU.Value() -= df->ComputeUpdateInv(nD, globalData)*maskprob;
        ++nD;
        ++nU;
        }
      else
        {
        ++nD;
        ++nU;
        }
      }

    if( updateenergy )
      {
      this->m_LastEnergy[metricCount] = this->m_Energy[metricCount];
      this->m_Energy[metricCount] = df->GetEnergy();
      // *this->m_SimilarityMetrics[metricCount]->GetWeightScalar()/sumWeights;
      }

    // smooth the fields
    //       if (!ispointsetmetric || ImageDimension == 2 ){
    this->SmoothDisplacementField(updateField, true);
    if( updateFieldInv )
      {
      this->SmoothDisplacementField(updateFieldInv, true);
      }
    //  }
    /*
    else // use another strategy -- exact lm? / something like Laplacian
  {
    TReal tmag=0;
    for (unsigned int ff=0; ff<5; ff++)
      {
        tmag=0;
        this->SmoothDisplacementField(updateField,true);
        if (updateFieldInv) this->SmoothDisplacementField(updateFieldInv,true);
        nD.GoToBegin();
        nU.GoToBegin();
        while( !nD.IsAtEnd() )
      {
        typename ImageType::IndexType index=nD.GetIndex();
        bool oktosample=true;
        TReal maskprob=1.0;
        if (mask)
          {
            maskprob=mask->GetPixel( nD.GetIndex() );
            if (maskprob > 1.0) maskprob=1.0;
            if ( maskprob < 0.1) oktosample=false;
          }
        VectorType F1;
        F1.Fill(0);
        VectorType F2;
        F2.Fill(0);
        if ( oktosample )
          {
            F1 = df->ComputeUpdate(nD, globalData)*maskprob;
            if (totalUpdateInvField)
          {
            F2 = df->ComputeUpdateInv(nD, globalData)*maskprob;
          }
            ++nD;
            ++nU;
          }
        else
          {
            ++nD;
            ++nU;
          }

        // compute mags of F1 and F2 -- if large enough, reset them
        TReal f1mag=0,f2mag=0,umag=0;
        for (unsigned int dim=0; dim<ImageDimension; dim++)
          {
            f1mag+=F1[dim]/spacing[dim]*F1[dim]/spacing[dim];
            f2mag+=F2[dim]/spacing[dim]*F2[dim]/spacing[dim];
            umag+=updateField->GetPixel(index)[dim]/spacing[dim]*updateField->GetPixel(index)[dim]/spacing[dim];
          }
        f1mag=sqrt(f1mag); f2mag=sqrt(f2mag); umag=sqrt(umag);
        if ( f1mag > 0.05 ) updateField->SetPixel(index,F1);
        if ( f2mag > 0.05 ) updateFieldInv->SetPixel(index,F2);
        tmag+=umag;
      }
        //           ::ants::antscout << " total mag " << tmag << std::endl;
      }
    //smooth the total field
    this->SmoothDisplacementField(updateField,true);
    if (updateFieldInv) this->SmoothDisplacementField(updateFieldInv,true);
  }

    */
    // normalize update field then add to total field
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    Iterator      dIter(totalUpdateField, totalUpdateField->GetLargestPossibleRegion() );
    TReal         mag = 0.0;
    TReal         max = 0.0;
    unsigned long ct = 0;
    TReal         total = 0;
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType vec = updateField->GetPixel(index);
      mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
        }
      mag = sqrt(mag);
//            if (mag > 0. ) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << std::endl;
      if( mag > max )
        {
        max = mag;
        }
      ct++;
      total += mag;
      //    ::ants::antscout << " mag " << mag << std::endl;
      }
    if( this->m_Debug )
      {
      ::ants::antscout << "PRE MAX " << max << std::endl;
      }
    TReal max2 = 0;
    if( max <= 0 )
      {
      max = 1;
      }
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType vec = updateField->GetPixel(index);
      vec = vec / max;
      mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
        }
      mag = sqrt(mag);
      if( mag >  max2 )
        {
        max2 = mag;
        }
//            if (mag > 0.95) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << " ind " << index
// <<
// std::endl;
      /** FIXME need weights between metrics */

      RealType normalizedWeight
        = this->m_SimilarityMetrics[metricCount]->GetWeightScalar() / sumWeights;
      /*            RealType weight = this->m_SimilarityMetrics[metricCount]->GetWeightImage()->GetPixel( diter.GetIndex() );
      if (ispointsetmetric )
        {
      VectorType intensityupdate=dIter.Get();
      VectorType lmupdate=vec;
      TReal lmag=0;
      for (unsigned int li=0; li<ImageDimension; li++) lmag+=(lmupdate[li]/spacing[li])*(lmupdate[li]/spacing[li]);
      lmag=sqrt(lmag);
      TReal modi=1;
      if (lmag > 1) modi=0;
      else modi=1.0-lmag;
      TReal iwt=1*modi;
      TReal lmwt=normalizedWeight;
      VectorType totalv=intensityupdate*iwt+lmupdate*lmwt;
      dIter.Set(totalv);
        }
        else */
      dIter.Set(dIter.Get() + vec * normalizedWeight);
      }

    if( totalUpdateInvField )
      {
      Iterator invDIter(totalUpdateInvField, totalUpdateInvField->GetLargestPossibleRegion() );
      mag = 0.0;
      max = 0.0;
      ct = 0;
      total = 0;
      for( invDIter.GoToBegin(); !invDIter.IsAtEnd(); ++invDIter )
        {
        typename ImageType::IndexType index = invDIter.GetIndex();
        VectorType vec = updateFieldInv->GetPixel(index);
        mag = 0;
        for( unsigned int jj = 0; jj < ImageDimension; jj++ )
          {
          mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
          }
        mag = sqrt(mag);
//            if (mag > 0. ) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << std::endl;
        if( mag > max )
          {
          max = mag;
          }
        ct++;
        total += mag;
        //    ::ants::antscout << " mag " << mag << std::endl;
        }
      if( this->m_Debug )
        {
        ::ants::antscout << "PRE MAX " << max << std::endl;
        }
      max2 = 0;
      if( max <= 0 )
        {
        max = 1;
        }
      for( invDIter.GoToBegin(); !invDIter.IsAtEnd(); ++invDIter )
        {
        typename ImageType::IndexType index = invDIter.GetIndex();
        VectorType vec = updateFieldInv->GetPixel(index);
        vec = vec / max;
        mag = 0;
        for( unsigned int jj = 0; jj < ImageDimension; jj++ )
          {
          mag += vec[jj] / spacing[jj] * vec[jj] / spacing[jj];
          }
        mag = sqrt(mag);
        if( mag >  max2 )
          {
          max2 = mag;
          }
//            if (mag > 0.95) ::ants::antscout << " mag " << mag << " max " << max << " vec " << vec << " ind " << index
// <<
// std::endl;
        /** FIXME need weights between metrics */

        RealType normalizedWeight
          = this->m_SimilarityMetrics[metricCount]->GetWeightScalar() / sumWeights;
//            RealType weight = this->m_SimilarityMetrics[metricCount]->GetWeightImage()->GetPixel( diter.GetIndex() );
        /*
        if (ispointsetmetric )
      {
    VectorType intensityupdate=invDIter.Get();
    VectorType lmupdate=vec;
    TReal lmag=0;
    for (unsigned int li=0; li<ImageDimension; li++) lmag+=(lmupdate[li]/spacing[li])*(lmupdate[li]/spacing[li]);
    lmag=sqrt(lmag);
    TReal modi=1;
    if (lmag > 1) modi=0;
    else modi=1.0-lmag;
    TReal iwt=1*modi;
    TReal lmwt=normalizedWeight;
    VectorType totalv=intensityupdate*iwt+lmupdate*lmwt;
    invDIter.Set(totalv);
      }
      else */
        invDIter.Set(invDIter.Get() + vec * normalizedWeight);
        }
      }
    if( this->m_Debug )
      {
      ::ants::antscout << "PO MAX " << max2 << " sz" << totalUpdateField->GetLargestPossibleRegion().GetSize()
                       << std::endl;
      }
    }

//    this->SmoothDisplacementField( totalUpdateField,true);
//    if (totalUpdateInvField) this->SmoothDisplacementField( totalUpdateInvField,true);

  return totalUpdateField;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::DiffeomorphicExpRegistrationUpdate(ImagePointer /* fixedImage */, ImagePointer movingImage, PointSetPointer fpoints,
                                     PointSetPointer mpoints)
{
  //  this function computes a velocity field that --- when composed with itself --- optimizes the registration
  // solution.
  // it's very different than the Exp approach used by DiffeomorphicDemons which just composes small deformation over
  // time.
  //  DiffDem cannot reconstruct the path between images without recomputing the registration.
  // DiffDem also cannot create an inverse mapping.
  /** FIXME really should pass an image list and then warp each one in
    turn  then expand the update field to fit size of total
    deformation */

  VectorType zero;

  zero.Fill(0);
  DisplacementFieldPointer totalUpdateField = NULL;
  DisplacementFieldPointer totalField = this->m_DisplacementField;
  /** generate phi and phi gradient */
  //    TReal timestep=1.0/(TReal)this->m_NTimeSteps;
  //    for (unsigned int nts=0; nts<=this->m_NTimeSteps; nts+=this->m_NTimeSteps)
  unsigned int nts = (unsigned int)this->m_NTimeSteps;
    {
    DisplacementFieldPointer diffmap = this->IntegrateConstantVelocity(totalField, nts, 1);
    // DisplacementFieldPointer invdiffmap = this->IntegrateConstantVelocity(totalField,(unsigned int)(
    // this->m_NTimeSteps)-nts, (-1.));

    ImagePointer           wfimage, wmimage;
    PointSetPointer        wfpoints = NULL, wmpoints = NULL;
    AffineTransformPointer aff = this->m_AffineTransform;
    if( mpoints )
      {       // need full inverse map
      DisplacementFieldPointer tinvdiffmap = this->IntegrateConstantVelocity(totalField, nts, (-1.) );
      wmpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  mpoints,  aff, tinvdiffmap, true,
                                          this->m_FixedImageAffineTransform );
      }

    DisplacementFieldPointer updateField = this->ComputeUpdateField( diffmap, NULL, fpoints, wmpoints);
    //    updateField = this->IntegrateConstantVelocity( updateField, nts, timestep);
    TReal maxl = this->MeasureDeformation(updateField);
    if( maxl <= 0 )
      {
      maxl = 1;
      }
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    Iterator dIter(updateField, updateField->GetLargestPossibleRegion() );
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      dIter.Set( dIter.Get() * this->m_GradstepAltered / maxl);
      }
    this->ComposeDiffs(updateField, totalField, totalField, 1);
    /*    TReal maxl= this->MeasureDeformation(updateField);
    if (maxl <= 0) maxl=1;
    typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
    Iterator dIter(updateField,updateField->GetLargestPossibleRegion() );
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter ) {
      dIter.Set( dIter.Get()*this->m_GradstepAltered/this->m_NTimeSteps );
      totalField->SetPixel(dIter.GetIndex(), dIter.Get() +  totalField->GetPixel(dIter.GetIndex()) );
        } */
    }
  this->SmoothDisplacementField(totalField, false);

  return;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::GreedyExpRegistrationUpdate(ImagePointer /* fixedImage */, ImagePointer /* movingImage */, PointSetPointer fpoints,
                              PointSetPointer /* mpoints */)
{
  //  similar approach to christensen 96 and diffeomorphic demons

  VectorType zero;

  zero.Fill(0);
  DisplacementFieldPointer totalUpdateField = NULL;

  // we compose the update with this field.
  DisplacementFieldPointer totalField = this->m_DisplacementField;

  TReal        timestep = 1.0 / (TReal) this->m_NTimeSteps;
  unsigned int nts = (unsigned int)this->m_NTimeSteps;

  ImagePointer             wfimage, wmimage;
  PointSetPointer          wfpoints = NULL, wmpoints = NULL;
  AffineTransformPointer   aff = this->m_AffineTransform;
  DisplacementFieldPointer updateField = this->ComputeUpdateField( totalField, NULL, fpoints, wmpoints);
  updateField = this->IntegrateConstantVelocity( updateField, nts, timestep);
  TReal maxl = this->MeasureDeformation(updateField);
  if( maxl <= 0 )
    {
    maxl = 1;
    }
  typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
  Iterator dIter(updateField, updateField->GetLargestPossibleRegion() );
  for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
    {
    dIter.Set( dIter.Get() * this->m_GradstepAltered / maxl );
    }
  this->ComposeDiffs(updateField, totalField, totalField, 1);
  //    maxl= this->MeasureDeformation(totalField);
  //    ::ants::antscout << " maxl " << maxl << " gsa " << this->m_GradstepAltered   << std::endl;
  //    totalField=this->CopyDisplacementField(totalUpdateField);
  this->SmoothDisplacementField(totalField, false);

  return;
}

// added by songgang
template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::AffineTransformPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>::AffineOptimization(OptAffineType & affine_opt)
{
  //    typedef itk::Image<TReal, 3> TempImageType;
  //    typename TempImageType::Pointer fixedImage = TempImageType::New();
  //    typename TempImageType::Pointer movingImage = TempImageType::New();

  ImagePointer fixedImage;
  ImagePointer movingImage;

  /** FIXME -- here we assume the metrics all have the same image */
  fixedImage = this->m_SimilarityMetrics[0]->GetFixedImage();
  movingImage = this->m_SimilarityMetrics[0]->GetMovingImage();

  // TODO: get mask image pointer / type for mask image
  if( this->m_MaskImage )
    {
    affine_opt.mask_fixed = this->m_MaskImage;
    }

  // AffineTransformPointer &transform_init = affine_opt.transform_initial;
  // ImagePointer &maskImage = affine_opt.mask_fixed;

  AffineTransformPointer transform = AffineTransformType::New();

  // ::ants::antscout << "In AffineOptimization: transform_init.IsNotNull()=" << transform_init.IsNotNull() <<
  // std::endl;
  // compute_single_affine_transform(fixedImage, movingImage, maskImage, transform, transform_init);

  // OptAffine<AffineTransformPointer, ImagePointer> opt;
  ComputeSingleAffineTransform<ImageType, TransformType, OptAffineType>(fixedImage,
                                                                        movingImage,
                                                                        affine_opt,
                                                                        transform);

  return transform;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SyNRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints,
                        PointSetPointer mpoints)
{
  VectorType zero;

  zero.Fill(0);
  DisplacementFieldPointer totalUpdateField, totalUpdateInvField =
    AllocImage<DisplacementFieldType>(this->m_DisplacementField, zero);

  if( !this->m_SyNF )
    {
    ::ants::antscout << " Allocating " << std::endl;
    this->m_SyNF = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNFInv = this->CopyDisplacementField(this->m_SyNF);
    this->m_SyNM = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNMInv = this->CopyDisplacementField(this->m_SyNF);
    // this->m_Debug=true;
    if( this->m_Debug )
      {
      ::ants::antscout << " SyNFInv" << this->m_SyNFInv->GetLargestPossibleRegion().GetSize() << std::endl;
      }
    if( this->m_Debug )
      {
      ::ants::antscout << " t updIf " << totalUpdateInvField->GetLargestPossibleRegion().GetSize() << std::endl;
      }
    if( this->m_Debug )
      {
      ::ants::antscout << " synf " << this->m_SyNF->GetLargestPossibleRegion().GetSize() << std::endl;
      }
    // this->m_Debug=false;
    ::ants::antscout << " Allocating Done " << std::endl;
    }

  if( !this->m_SyNF )
    {
    ::ants::antscout << " F'D UP " << std::endl;
    }

  PointSetPointer        wfpoints = NULL, wmpoints = NULL;
  AffineTransformPointer aff = this->m_AffineTransform;
  AffineTransformPointer affinverse = NULL;
  if( aff )
    {
    affinverse = AffineTransformType::New();
    aff->GetInverse(affinverse);
    }

  if( mpoints )
    {
    wmpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  mpoints,  aff, this->m_SyNM, true,
                                        this->m_FixedImageAffineTransform );
    }

  if( fpoints )
    {  // need full inverse map
    wfpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, fixedImage, fpoints,  NULL, this->m_SyNF, false,
                                        this->m_FixedImageAffineTransform  );
    }
  // syncom
  totalUpdateField =
    this->ComputeUpdateField(this->m_SyNMInv, this->m_SyNFInv, wfpoints, wmpoints, totalUpdateInvField);

  this->ComposeDiffs(this->m_SyNF, totalUpdateField, this->m_SyNF, this->m_GradstepAltered);
  this->ComposeDiffs(this->m_SyNM, totalUpdateInvField, this->m_SyNM, this->m_GradstepAltered);

  if( this->m_TotalSmoothingparam > 0 || this->m_TotalSmoothingMeshSize[0] > 0 )
    {
    this->SmoothDisplacementField( this->m_SyNF, false);
    this->SmoothDisplacementField( this->m_SyNM, false);
    }

  this->InvertField(this->m_SyNF, this->m_SyNFInv);
  this->InvertField(this->m_SyNM, this->m_SyNMInv);
  this->InvertField(this->m_SyNFInv, this->m_SyNF);
  this->InvertField(this->m_SyNMInv, this->m_SyNM);

//      ::ants::antscout <<  " F " << this->MeasureDeformation(this->m_SyNF) << " F1 " <<
// this->MeasureDeformation(this->m_SyNFInv) << std::endl;
//      ::ants::antscout <<  " M " << this->MeasureDeformation(this->m_SyNM) << " M1 " <<
// this->MeasureDeformation(this->m_SyNMInv) << std::endl;
  return;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SyNExpRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints,
                           PointSetPointer mpoints)
{
  ::ants::antscout << " SyNEX";

  typename ImageType::SpacingType spacing = fixedImage->GetSpacing();
  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer totalUpdateField, totalUpdateInvField =
    AllocImage<DisplacementFieldType>(this->m_DisplacementField, zero);

  if( !this->m_SyNF )
    {
    ::ants::antscout << " Allocating " << std::endl;
    this->m_SyNF = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNFInv = this->CopyDisplacementField(this->m_SyNF);
    this->m_SyNM = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNMInv = this->CopyDisplacementField(this->m_SyNF);
    ::ants::antscout << " Allocating Done " << std::endl;
    }

  if( !this->m_SyNF )
    {
    ::ants::antscout << " F'D UP " << std::endl;
    }

  ImagePointer           wfimage, wmimage;
  PointSetPointer        wfpoints = NULL, wmpoints = NULL;
  AffineTransformPointer aff = this->m_AffineTransform;
  AffineTransformPointer affinverse = NULL;

// here, SyNF holds the moving velocity field, SyNM holds the fixed
// velocity field and we integrate both to generate the inv/fwd fields

  TReal                    timestep = 1.0 / (TReal) this->m_NTimeSteps;
  unsigned int             nts = this->m_NTimeSteps;
  DisplacementFieldPointer fdiffmap = this->IntegrateConstantVelocity(this->m_SyNF, nts, 1);
  this->m_SyNFInv = this->IntegrateConstantVelocity(this->m_SyNF, nts, (-1.) );
  DisplacementFieldPointer mdiffmap = this->IntegrateConstantVelocity(this->m_SyNM, nts, 1);
  this->m_SyNMInv = this->IntegrateConstantVelocity(this->m_SyNM, nts, (-1.) );

  if( aff )
    {
    affinverse = AffineTransformType::New();
    aff->GetInverse(affinverse);
    }
  if( mpoints )
    {
    wmpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  mpoints,  aff, this->m_SyNM, true,
                                        this->m_FixedImageAffineTransform );
    }
  if( fpoints )
    {  // need full inverse map
    wfpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, fixedImage, fpoints,  NULL, this->m_SyNF, false,
                                        this->m_FixedImageAffineTransform );
    }

  totalUpdateField = this->ComputeUpdateField( this->m_SyNMInv, this->m_SyNFInv, wfpoints, wmpoints,
                                               totalUpdateInvField);
// then addd
  typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
  Iterator dIter(this->m_SyNF, this->m_SyNF->GetLargestPossibleRegion() );
  TReal    max = 0, max2 = 0;
  for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
    {
    typename ImageType::IndexType index = dIter.GetIndex();
    VectorType vecf = totalUpdateField->GetPixel(index);
    VectorType vecm = totalUpdateInvField->GetPixel(index);
    dIter.Set(dIter.Get() + vecf * this->m_GradstepAltered);
    this->m_SyNM->SetPixel( index, this->m_SyNM->GetPixel( index ) + vecm * this->m_GradstepAltered);
// min field difference => geodesic => DV/dt=0
    TReal      geowt1 = 0.99;
    TReal      geowt2 = 1.0 - geowt1;
    VectorType synmv = this->m_SyNM->GetPixel( index );
    VectorType synfv   = this->m_SyNF->GetPixel( index );
    this->m_SyNM->SetPixel( index, synmv * geowt1 - synfv * geowt2);
    this->m_SyNF->SetPixel( index, synfv * geowt1 - synmv * geowt2);
    }

  if( this->m_TotalSmoothingparam > 0 || this->m_TotalSmoothingMeshSize[0] > 0 )
    {
    this->SmoothDisplacementField( this->m_SyNF, false);
    this->SmoothDisplacementField( this->m_SyNM, false);
    }
//      ::ants::antscout <<  " TUF " << this->MeasureDeformation(this->m_SyNF) << " TUM " <<
// this->MeasureDeformation(this->m_SyNM) << std::endl;

  return;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::UpdateTimeVaryingVelocityFieldWithSyNFandSyNM()
{
  typedef TimeVaryingVelocityFieldType tvt;

  bool generatetvfield = false;
  if( !this->m_TimeVaryingVelocity )
    {
    generatetvfield = true;
    }
  else
    {
    for( int jj = 0; jj < ImageDimension; jj++ )
      {
      if( this->m_CurrentDomainSize[jj] !=  this->m_TimeVaryingVelocity->GetLargestPossibleRegion().GetSize()[jj] )
        {
        generatetvfield = true;
        }
      }
    }
  VectorType zero;
  zero.Fill(0);
  if( generatetvfield )
    {
    typename tvt::RegionType gregion;
    typename tvt::SizeType gsize;
    typename tvt::SpacingType gspace;
    typename tvt::PointType gorigin;
    gorigin.Fill(0);
    for( unsigned int dim = 0; dim < TDimension; dim++ )
      {
      gsize[dim] = this->m_CurrentDomainSize[dim];
      gspace[dim] = this->m_CurrentDomainSpacing[dim];
      gorigin[dim] = this->m_CurrentDomainOrigin[dim];
      }
    gsize[TDimension] = (long unsigned int) this->m_NTimeSteps;
    gspace[TDimension] = 1;
    gregion.SetSize(gsize);

/** The TV Field has the direction of the sub-image -- the time domain
    has identity transform */
    typename tvt::DirectionType iddir;
    iddir.Fill(0);
    iddir[ImageDimension][ImageDimension] = 1;
    for( unsigned int i = 0; i < ImageDimension + 1; i++ )
      {
      for( unsigned int j = 0; j < ImageDimension + 1; j++ )
        {
//    if (i == j) iddir[i][j]=1;
        if( i < ImageDimension && j < ImageDimension )
          {
          iddir[i][j] = this->GetDisplacementField()->GetDirection()[i][j];
          }
        }
      }

    this->m_TimeVaryingVelocity =
      AllocImage<tvt>(gregion,
                      gspace,
                      gorigin,
                      iddir,
                      zero);
    }

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;

  TVFieldIterator m_FieldIter( this->m_TimeVaryingVelocity, this->m_TimeVaryingVelocity->GetLargestPossibleRegion() );
  for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
    {
    typename tvt::IndexType velind = m_FieldIter.GetIndex();
    IndexType ind;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      ind[j] = velind[j];
      }
    if( velind[ImageDimension] == 0 )
      {
      VectorType vel = this->m_SyNF->GetPixel(ind);
      m_FieldIter.Set(vel);
      }
    else if( velind[ImageDimension] == 1 )
      {
      VectorType vel = this->m_SyNM->GetPixel(ind) * (-1.0);
      m_FieldIter.Set(vel);
      }
    }

//  ::ants::antscout <<" ALlocated TV F "<< std::endl;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::CopyOrAddToVelocityField( TimeVaryingVelocityFieldPointer velocity,  DisplacementFieldPointer update1,
                            DisplacementFieldPointer update2,
                            TReal timept)
{
  typedef TimeVaryingVelocityFieldType tvt;
  VectorType zero;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;
  int tpupdate = (unsigned int) ( ( (TReal) this->m_NTimeSteps - 1.0) * timept + 0.5);
  // ::ants::antscout <<"  add to " << tpupdate << std::endl;
  TReal           tmag = 0;
  TVFieldIterator m_FieldIter(velocity, velocity->GetLargestPossibleRegion() );
  for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
    {
    typename tvt::IndexType velind = m_FieldIter.GetIndex();
    typename ImageType::IndexType ind;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      ind[j] = velind[j];
      }
    if( velind[ImageDimension] == tpupdate && update1 )
      {
      VectorType vel = update1->GetPixel(ind);
      TReal      mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += vel[jj] * vel[jj];
        }
      tmag += sqrt(mag);
      m_FieldIter.Set(vel + m_FieldIter.Get() );
      }
    if( velind[ImageDimension] == tpupdate && update2 )
      {
      VectorType vel = update2->GetPixel(ind) * (-1);
      m_FieldIter.Set(vel + m_FieldIter.Get() );
      }
    }
  //  ::ants::antscout << " tmag " << tmag << std::endl;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::SyNTVRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints,
                          PointSetPointer mpoints)
{
  VectorType zero;

  zero.Fill(0);
  DisplacementFieldPointer totalUpdateField;
  DisplacementFieldPointer totalUpdateInvField =
    AllocImage<DisplacementFieldType>(this->m_DisplacementField, zero);
  if( !this->m_SyNF )
    {
    ::ants::antscout << " Allocating " << std::endl;
    this->m_SyNF = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNFInv = this->CopyDisplacementField(this->m_SyNF);
    this->m_SyNM = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNMInv = this->CopyDisplacementField(this->m_SyNF);
    ::ants::antscout << " Allocating Done " << std::endl;
    }

  if( !this->m_SyNF )
    {
    ::ants::antscout << " F'D UP " << std::endl;
    }

  ImagePointer           wfimage, wmimage;
  PointSetPointer        wfpoints = NULL, wmpoints = NULL;
  AffineTransformPointer aff = this->m_AffineTransform;
  AffineTransformPointer affinverse = NULL;

  typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
  Iterator dIter(this->m_SyNF, this->m_SyNF->GetLargestPossibleRegion() );

// here, SyNF holds the moving velocity field, SyNM holds the fixed
// velocity field and we integrate both to generate the inv/fwd fields
  typename JacobianFunctionType::Pointer jfunction = JacobianFunctionType::New();
  TReal lot = 0, hit = 0.5;
  TReal lot2 = 1.0;
  this->UpdateTimeVaryingVelocityFieldWithSyNFandSyNM();   // sets tvt to SyNF and SyNM -- works only for 2 time points!
  this->m_SyNFInv = this->IntegrateVelocity(hit, lot);
  this->m_SyNMInv = this->IntegrateVelocity(hit, lot2);
  if( aff )
    {
    affinverse = AffineTransformType::New();
    aff->GetInverse(affinverse);
    }
  if( mpoints )
    {
/**FIXME -- NEED INTEGRATION FOR POINTS ONLY  -- warp landmarks for
* tv-field */
//      ::ants::antscout <<" aff " << std::endl;
/** NOte, totalUpdateInvField is filled with zeroes! -- we only want
      affine mapping */
    wmpoints =
      this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  mpoints,  aff, totalUpdateInvField, true,
                               NULL );
    DisplacementFieldPointer mdiffmap = this->IntegrateLandmarkSetVelocity(lot2, hit, wmpoints, movingImage);
    wmpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  wmpoints,  NULL, mdiffmap, true,
                                        NULL );
    }
  if( fpoints )
    {  // need full inverse map
    wfpoints =
      this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  fpoints, NULL, totalUpdateInvField, true,
                               this->m_FixedImageAffineTransform );
    DisplacementFieldPointer fdiffmap = this->IntegrateLandmarkSetVelocity(lot, hit, wfpoints, fixedImage);
    wfpoints =
      this->WarpMultiTransform(this->m_ReferenceSpaceImage, fixedImage, wfpoints,  NULL, fdiffmap, false, NULL );
    }
  totalUpdateField =
    this->ComputeUpdateField( this->m_SyNMInv, this->m_SyNFInv, wfpoints, wmpoints, totalUpdateInvField,
                              true);
  for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
    {
    typename ImageType::IndexType index = dIter.GetIndex();
    VectorType vecf = totalUpdateField->GetPixel(index) * 1;
    VectorType vecm = totalUpdateInvField->GetPixel(index);
// update time components 1 & 2
    this->m_SyNF->SetPixel( index, this->m_SyNF->GetPixel( index ) + vecf * this->m_GradstepAltered);
    this->m_SyNM->SetPixel( index, this->m_SyNM->GetPixel( index ) + vecm * this->m_GradstepAltered);
// min field difference => geodesic => DV/dt=0
    TReal      geowt1 = 0.95;
    TReal      geowt2 = 1.0 - geowt1;
    VectorType synmv = this->m_SyNM->GetPixel( index );
    VectorType synfv   = this->m_SyNF->GetPixel( index );
    this->m_SyNM->SetPixel( index, synmv * geowt1 - synfv * geowt2);
    this->m_SyNF->SetPixel( index, synfv * geowt1 - synmv * geowt2);
    }

  if(  this->m_TotalSmoothingparam > 0
       || this->m_TotalSmoothingMeshSize[0] > 0 )
    {
    // smooth time components separately
    this->SmoothDisplacementField( this->m_SyNF, false);
    this->SmoothDisplacementField( this->m_SyNM, false);
    }

  return;
}

template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::DiReCTUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints, PointSetPointer mpoints)
{
  typedef TimeVaryingVelocityFieldType tvt;
  TimeVaryingVelocityFieldPointer velocityUpdate = NULL;

  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer totalUpdateInvField
    = AllocImage<DisplacementFieldType>(this->m_DisplacementField, zero);

  unsigned long numpx = this->m_DisplacementField->GetBufferedRegion().GetNumberOfPixels();

  bool generatetvfield = false;
  bool enlargefield = false;
  if( !this->m_TimeVaryingVelocity )
    {
    generatetvfield = true;
    }
  else
    {
    for( int jj = 0; jj < ImageDimension; jj++ )
      {
      if( this->m_CurrentDomainSize[jj] !=  this->m_TimeVaryingVelocity->GetLargestPossibleRegion().GetSize()[jj] )
        {
        enlargefield = true;
        }
      }
    }

  velocityUpdate = tvt::New();
  typename tvt::RegionType gregion;
  typename tvt::SizeType gsize;
  typename tvt::SpacingType gspace;
  typename tvt::PointType gorigin;
  gorigin.Fill(0);
  for( unsigned int dim = 0; dim < TDimension; dim++ )
    {
    gsize[dim] = this->m_CurrentDomainSize[dim];
    gspace[dim] = this->m_CurrentDomainSpacing[dim];
    gorigin[dim] = this->m_CurrentDomainOrigin[dim];
    }
  if( this->m_NTimeSteps < 2 )
    {
    this->m_NTimeSteps = 2;
    }
  gsize[TDimension] = (unsigned long) this->m_NTimeSteps;
  TReal hitstep = 1.0 / ( (TReal) this->m_NTimeSteps - 1);
  gspace[TDimension] = 1;
  gregion.SetSize(gsize);
  velocityUpdate->SetSpacing( gspace );
  velocityUpdate->SetOrigin( gorigin );

/** The TV Field has the direction of the sub-image -- the time domain
    has identity transform */
  typename tvt::DirectionType iddir;
  iddir.Fill(0);
  iddir[ImageDimension][ImageDimension] = 1;
  for( unsigned int i = 0; i < ImageDimension + 1; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension + 1; j++ )
      {
//    if (i == j) iddir[i][j]=1;
      if( i < ImageDimension && j < ImageDimension )
        {
        iddir[i][j] = this->GetDisplacementField()->GetDirection()[i][j];
        }
      }
    }

  velocityUpdate->SetDirection( iddir );
  velocityUpdate->SetLargestPossibleRegion(gregion);
  velocityUpdate->SetRequestedRegion( gregion);
  velocityUpdate->SetBufferedRegion( gregion  );
  velocityUpdate->Allocate();
  velocityUpdate->FillBuffer(zero);

  if( generatetvfield )
    {
    this->m_TimeVaryingVelocity = tvt::New();
    this->m_TimeVaryingVelocity->SetSpacing( gspace );
    this->m_TimeVaryingVelocity->SetOrigin( gorigin );
    this->m_TimeVaryingVelocity->SetDirection( iddir );
    this->m_TimeVaryingVelocity->SetLargestPossibleRegion(gregion);
    this->m_TimeVaryingVelocity->SetRequestedRegion( gregion);
    this->m_TimeVaryingVelocity->SetBufferedRegion( gregion  );
    this->m_TimeVaryingVelocity->Allocate();
    this->m_TimeVaryingVelocity->FillBuffer(zero);
    /*    this->m_LastTimeVaryingVelocity=tvt::New();
    this->m_LastTimeVaryingVelocity->SetSpacing( gspace );
    this->m_LastTimeVaryingVelocity->SetOrigin( gorigin );
    this->m_LastTimeVaryingVelocity->SetDirection( iddir );
    this->m_LastTimeVaryingVelocity->SetLargestPossibleRegion(gregion);
    this->m_LastTimeVaryingVelocity->SetRequestedRegion( gregion);
    this->m_LastTimeVaryingVelocity->SetBufferedRegion( gregion  );
    this->m_LastTimeVaryingVelocity->Allocate();
    this->m_LastTimeVaryingVelocity->FillBuffer(zero); */
    this->m_LastTimeVaryingUpdate = tvt::New();
    this->m_LastTimeVaryingUpdate->SetSpacing( gspace );
    this->m_LastTimeVaryingUpdate->SetOrigin( gorigin );
    this->m_LastTimeVaryingUpdate->SetDirection( iddir );
    this->m_LastTimeVaryingUpdate->SetLargestPossibleRegion(gregion);
    this->m_LastTimeVaryingUpdate->SetRequestedRegion( gregion);
    this->m_LastTimeVaryingUpdate->SetBufferedRegion( gregion  );
    this->m_LastTimeVaryingUpdate->Allocate();
    this->m_LastTimeVaryingUpdate->FillBuffer(zero);
    }
  else if( enlargefield )
    {
    this->m_TimeVaryingVelocity = this->ExpandVelocity();
    this->m_TimeVaryingVelocity->SetSpacing(gspace);
    this->m_TimeVaryingVelocity->SetOrigin(gorigin);
    /*        this->m_LastTimeVaryingVelocity=tvt::New();
    this->m_LastTimeVaryingVelocity->SetSpacing( gspace );
    this->m_LastTimeVaryingVelocity->SetOrigin( gorigin );
    this->m_LastTimeVaryingVelocity->SetDirection( iddir );
    this->m_LastTimeVaryingVelocity->SetLargestPossibleRegion(gregion);
    this->m_LastTimeVaryingVelocity->SetRequestedRegion( gregion);
    this->m_LastTimeVaryingVelocity->SetBufferedRegion( gregion  );
    this->m_LastTimeVaryingVelocity->Allocate();
    this->m_LastTimeVaryingVelocity->FillBuffer(zero);*/
    this->m_LastTimeVaryingUpdate = tvt::New();
    this->m_LastTimeVaryingUpdate->SetSpacing( gspace );
    this->m_LastTimeVaryingUpdate->SetOrigin( gorigin );
    this->m_LastTimeVaryingUpdate->SetDirection( iddir );
    this->m_LastTimeVaryingUpdate->SetLargestPossibleRegion(gregion);
    this->m_LastTimeVaryingUpdate->SetRequestedRegion( gregion);
    this->m_LastTimeVaryingUpdate->SetBufferedRegion( gregion  );
    this->m_LastTimeVaryingUpdate->Allocate();
    this->m_LastTimeVaryingUpdate->FillBuffer(zero);
    }
  if( !this->m_SyNF )
    {
    ::ants::antscout << " Allocating " << std::endl;
    this->m_SyNF = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNFInv = this->CopyDisplacementField(this->m_SyNF);
    this->m_SyNM = this->CopyDisplacementField(totalUpdateInvField);
    this->m_SyNMInv = this->CopyDisplacementField(this->m_SyNF);
    ::ants::antscout << " Allocating Done " << std::endl;
    }

  if( !this->m_SyNF )
    {
    ::ants::antscout << " F'D UP " << std::endl;
    }

  ImagePointer           wfimage, wmimage;
  PointSetPointer        wfpoints = NULL, wmpoints = NULL;
  AffineTransformPointer aff = this->m_AffineTransform;
  AffineTransformPointer affinverse = NULL;

  typedef ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;
  Iterator dIter(this->m_SyNF, this->m_SyNF->GetLargestPossibleRegion() );

// here, SyNF holds the moving velocity field, SyNM holds the fixed
// velocity field and we integrate both to generate the inv/fwd fields
  typename JacobianFunctionType::Pointer jfunction = JacobianFunctionType::New();
  TReal        lot = 0, lot2 = 1.0;
  unsigned int fct = 100;
  for( TReal hit = 0; hit <= 1; hit = hit + hitstep )
    {
    this->m_SyNFInv = this->IntegrateVelocity(hit, lot);
    this->m_SyNMInv = this->IntegrateVelocity(hit, lot2);

    if( false && this->m_CurrentIteration == 1 &&  this->m_SyNFInv  )
      {
      typedef itk::ImageFileWriter<DisplacementFieldType>
        DisplacementFieldWriterType;
      typename DisplacementFieldWriterType::Pointer writer = DisplacementFieldWriterType::New();
      std::ostringstream osstream;
      osstream << fct;
      fct++;
      std::string fnm = std::string("field1") + osstream.str() + std::string("warp.nii.gz");
      writer->SetInput( this->m_SyNFInv );
      writer->SetFileName( fnm.c_str() );
      ::ants::antscout << " write " << fnm << std::endl;
      writer->Update();
      // writer->SetInput( this->m_SyNMInv );
      // writer->SetFileName( fnm2.c_str() );
      //      writer->Update();
      }

    if( aff )
      {
      affinverse = AffineTransformType::New();
      aff->GetInverse(affinverse);
      }
    if( mpoints )
      {
/**FIXME -- NEED INTEGRATION FOR POINTS ONLY  -- warp landmarks for
* tv-field */
//      ::ants::antscout <<" aff " << std::endl;
/** NOte, totalUpdateInvField is filled with zeroes! -- we only want
      affine mapping */
      wmpoints =
        this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  mpoints,  aff, totalUpdateInvField, true,
                                 NULL );
      DisplacementFieldPointer mdiffmap = this->IntegrateLandmarkSetVelocity(lot2, hit, wmpoints, movingImage);
      wmpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  wmpoints,  NULL, mdiffmap, true,
                                          NULL );
      }
    if( fpoints )
      { // need full inverse map
      wfpoints =
        this->WarpMultiTransform(this->m_ReferenceSpaceImage, movingImage,  fpoints, NULL, totalUpdateInvField, true,
                                 this->m_FixedImageAffineTransform );
      DisplacementFieldPointer fdiffmap = this->IntegrateLandmarkSetVelocity(lot, hit, wfpoints, fixedImage);
      wfpoints = this->WarpMultiTransform(this->m_ReferenceSpaceImage, fixedImage, wfpoints,  NULL, fdiffmap, false,
                                          NULL );
      }
    DisplacementFieldPointer totalUpdateField =
      this->ComputeUpdateField( this->m_SyNMInv, this->m_SyNFInv,
                                wfpoints, wmpoints, totalUpdateInvField,
                                true);
    if( this->m_SyNFullTime == 2 )
      {
      totalUpdateInvField = NULL;
      }
    this->CopyOrAddToVelocityField( velocityUpdate, totalUpdateField,  totalUpdateInvField,  hit );
    }
  // http://en.wikipedia.org/wiki/Conjugate_gradient_method
  // below is r_k+1
  this->SmoothVelocityGauss( velocityUpdate,  this->m_GradSmoothingparam, ImageDimension );
  // update total velocity with v-update
  TReal tmag = 0;
  typedef itk::ImageRegionIteratorWithIndex<tvt> TVFieldIterator;
  TVFieldIterator m_FieldIter( this->m_TimeVaryingVelocity, this->m_TimeVaryingVelocity->GetLargestPossibleRegion() );
  for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
    {
    TReal      A = 1;
    TReal      alpha = 0, alpha1 = 0, alpha2 = 0, beta = 0, beta1 = 0, beta2 = 0;
    VectorType vec1 = velocityUpdate->GetPixel(m_FieldIter.GetIndex() );               // r_k+1
    VectorType vec2 = vec1;
    // this->m_LastTimeVaryingVelocity->GetPixel(m_FieldIter.GetIndex());
    // r_k
    VectorType upd = this->m_LastTimeVaryingUpdate->GetPixel(m_FieldIter.GetIndex() ); // p_k
    for( unsigned int ii = 0; ii < ImageDimension; ii++ )
      {
      alpha1 = vec2[ii] * vec2[ii];
      alpha2 = upd[ii] * upd[ii];
      beta1 = vec1[ii] * vec1[ii];
      beta2 = vec2[ii] * vec2[ii];
      }
    if( alpha2 > 0 )
      {
      alpha = alpha1 / (A * alpha2 + 0.001);
      }
    if( beta2 > 0 )
      {
      beta = beta1 / (beta2 + 0.001);
      }
    if( beta > 1 )
      {
      beta = 1;
      }
    if( alpha > 1 )
      {
      alpha = 1;
      }
    //        ::ants::antscout <<" beta " << beta << " alpha " << alpha << " it " << this->m_CurrentIteration <<
    // std::endl;
    VectorType newupd = (vec1);
    if( this->m_CurrentIteration  > 2 )
      {
      newupd = (vec1 + upd) * 0.5;
      }
    VectorType newsoln = m_FieldIter.Get() + this->m_GradstepAltered * newupd;
    m_FieldIter.Set( newsoln );
    //    VectorType vec2u=vec2 - alpha*A*upd;
    //    this->m_LastTimeVaryingVelocity->SetPixel(m_FieldIter.GetIndex(), vec1 );
    this->m_LastTimeVaryingUpdate->SetPixel(m_FieldIter.GetIndex(), newupd);
    TReal      mag = 0;
    VectorType vv = m_FieldIter.Get();
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      mag += vv[jj] * vv[jj];
      }
    tmag += sqrt(mag);
    }
  tmag /= ( (TReal) this->m_NTimeSteps * (TReal)numpx);
  ::ants::antscout << " DiffLength " << tmag << std::endl;
  if(  this->m_TotalSmoothingparam > 0
       || this->m_TotalSmoothingMeshSize[0] > 0 )
    {
    this->SmoothVelocityGauss( this->m_TimeVaryingVelocity,  this->m_TotalSmoothingparam, ImageDimension );
    //      this->SmoothDisplacementField( this->m_SyNF,false);
    //      this->SmoothDisplacementField( this->m_SyNM,false);
    }

  return;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::IntegrateVelocity(TReal starttimein, TReal finishtimein )
{
  ImagePointer mask = NULL;

  if( this->m_SyNMInv && this->m_MaskImage )
    {
    mask = this->WarpMultiTransform( this->m_ReferenceSpaceImage, this->m_MaskImage, NULL, this->m_SyNMInv, false,
                                     this->m_FixedImageAffineTransform );
    }
  else if( this->m_MaskImage )
    {
    mask = this->SubsampleImage( this->m_MaskImage, this->m_ScaleFactor,
                                 this->m_MaskImage->GetOrigin(), this->m_MaskImage->GetDirection(),  NULL);
    }

//  ::ants::antscout << " st " << starttimein << " ft " << finishtimein << std::endl;
  typedef TReal PixelType;
//   typedef itk::Vector<TReal, TDimension>     VectorType;
//   typedef itk::Image<VectorType, TDimension> DisplacementFieldType;
//   typedef itk::Image<PixelType, TDimension>  ImageType;
//   typedef typename  ImageType::IndexType     IndexType;
  typedef typename  ImageType::SizeType    SizeType;
  typedef typename  ImageType::SpacingType SpacingType;
  typedef TimeVaryingVelocityFieldType     tvt;

  bool dothick = false;
  if(  finishtimein > starttimein  && this->m_ComputeThickness )
    {
    dothick = true;
    }
  if( dothick && this->m_CurrentIteration  > 2 )
    {
    this->m_ThickImage =
      AllocImage<ImageType>(this->m_SyNF, 0);

    this->m_HitImage =
      AllocImage<ImageType>(this->m_SyNF, 0);
    }
  else
    {
    this->m_HitImage = NULL;  this->m_ThickImage = NULL;
    }

  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer intfield =
    AllocImage<DisplacementFieldType>
      (this->m_DisplacementField->GetLargestPossibleRegion(), zero);
  intfield->SetSpacing( this->m_CurrentDomainSpacing );
  intfield->SetOrigin(  this->m_DisplacementField->GetOrigin() );
  intfield->SetDirection(  this->m_DisplacementField->GetDirection() );
  if( starttimein == finishtimein )
    {
    return intfield;
    }
  if( !this->m_TimeVaryingVelocity )
    {
    ::ants::antscout << " No TV Field " << std::endl;  return intfield;
    }
  this->m_VelocityFieldInterpolator->SetInputImage(this->m_TimeVaryingVelocity);

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;

  if( starttimein < 0 )
    {
    starttimein = 0;
    }
  if( starttimein > 1 )
    {
    starttimein = 1;
    }
  if( finishtimein < 0 )
    {
    finishtimein = 0;
    }
  if( finishtimein > 1 )
    {
    finishtimein = 1;
    }

  FieldIterator m_FieldIter(this->GetDisplacementField(), this->GetDisplacementField()->GetLargestPossibleRegion() );
//  ::ants::antscout << " Start Int " << starttimein <<  std::endl;
  if( mask  && !this->m_ComputeThickness )
    {
    for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
      {
      IndexType  velind = m_FieldIter.GetIndex();
      VectorType disp;
      if( mask->GetPixel(velind) > 0.05 )
        {
        disp = this->IntegratePointVelocity(starttimein, finishtimein, velind) * mask->GetPixel(velind);
        }
      else
        {
        disp.Fill(0);
        }
      intfield->SetPixel(velind, disp);
      }
    }
  else
    {
    for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
      {
      IndexType  velind = m_FieldIter.GetIndex();
      VectorType disp = this->IntegratePointVelocity(starttimein, finishtimein, velind);
      intfield->SetPixel(velind, disp);
      }
    }
  if( this->m_ThickImage && this->m_MaskImage )
    {
    std::string outname = this->localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str() ) + std::string(
        "thick.nii.gz");
    ::ants::antscout << " write " << outname << std::endl;
    WriteImage<ImageType>(this->m_ThickImage, outname.c_str() );
    }

  return intfield;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DisplacementFieldPointer
ANTSImageRegistrationOptimizer<TDimension, TReal>
::IntegrateLandmarkSetVelocity(TReal starttimein, TReal finishtimein,
                               typename ANTSImageRegistrationOptimizer<TDimension,
                                                                       TReal>::PointSetPointer mypoints,
                               typename ANTSImageRegistrationOptimizer<TDimension,
                                                                       TReal>
                               ::ImagePointer /* refimage */ )
{
  typedef TimeVaryingVelocityFieldType tvt;

  VectorType zero;
  zero.Fill(0);
  DisplacementFieldPointer intfield = AllocImage<DisplacementFieldType>
      (this->m_DisplacementField->GetLargestPossibleRegion(), zero);
  intfield->SetSpacing( this->m_CurrentDomainSpacing );
  intfield->SetOrigin(  this->m_DisplacementField->GetOrigin() );
  intfield->SetDirection(  this->m_DisplacementField->GetDirection() );

  if( starttimein == finishtimein )
    {
    return intfield;
    }
  if( !this->m_TimeVaryingVelocity )
    {
    ::ants::antscout << " No TV Field " << std::endl;  return intfield;
    }
  this->m_VelocityFieldInterpolator->SetInputImage(this->m_TimeVaryingVelocity);

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;

  if( starttimein < 0 )
    {
    starttimein = 0;
    }
  if( starttimein > 1 )
    {
    starttimein = 1;
    }
  if( finishtimein < 0 )
    {
    finishtimein = 0;
    }
  if( finishtimein > 1 )
    {
    finishtimein = 1;
    }

  unsigned long sz1 = mypoints->GetNumberOfPoints();
  for( unsigned long ii = 0; ii < sz1; ii++ )
    {
    PointType point;
    // ::ants::antscout <<" get point " << std::endl;
    const bool isValidPoint = mypoints->GetPoint(ii, &point);
    if( !isValidPoint )
      {
      itkExceptionMacro( << "Invalid Point found at " << ii );
      }
    else
      {
      // ::ants::antscout <<" get point index " << point << std::endl;

      ImagePointType pt, wpt;
      for( unsigned int jj = 0;  jj < ImageDimension; jj++ )
        {
        pt[jj] = point[jj];
        }
      IndexType  velind;
      const bool bisinside = intfield->TransformPhysicalPointToIndex( pt, velind );
      // ::ants::antscout <<" inside? " << bisinside  << std::endl;
      if( bisinside )
        {
//      ::ants::antscout <<  "integrate " << std::endl;
        VectorType disp = this->IntegratePointVelocity(starttimein, finishtimein, velind);
//      ::ants::antscout <<  "put inside " << std::endl;
        intfield->SetPixel(velind, disp);
        }
      }
    }

  return intfield;
}

template <unsigned int TDimension, class TReal>
typename ANTSImageRegistrationOptimizer<TDimension, TReal>::VectorType
ANTSImageRegistrationOptimizer<TDimension, TReal>
::IntegratePointVelocity(TReal starttimein, TReal finishtimein, IndexType velind)
{
  typedef Point<TReal, itkGetStaticConstMacro(ImageDimension + 1)> xPointType;
  this->m_Debug = false;
//  ::ants::antscout <<"Enter IP "<< std::endl;
  typedef typename VelocityFieldInterpolatorType::OutputType InterpPointType;

  typedef TimeVaryingVelocityFieldType tvt;

  VectorType zero;
  zero.Fill(0);
  if( starttimein == finishtimein )
    {
    return zero;
    }

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;

  this->m_VelocityFieldInterpolator->SetInputImage(this->m_TimeVaryingVelocity);

  TReal        dT = this->m_DeltaTime;
  unsigned int m_NumberOfTimePoints = this->m_TimeVaryingVelocity->GetLargestPossibleRegion().GetSize()[TDimension];
  if( starttimein < 0 )
    {
    starttimein = 0;
    }
  if( starttimein > 1 )
    {
    starttimein = 1;
    }
  if( finishtimein < 0 )
    {
    finishtimein = 0;
    }
  if( finishtimein > 1 )
    {
    finishtimein = 1;
    }

  TReal timesign = 1.0;
  if( starttimein  >  finishtimein )
    {
    timesign = -1.0;
    }

  VectorType velo;
  velo.Fill(0);
  typename VelocityFieldInterpolatorType::ContinuousIndexType  vcontind;

  TReal         itime = starttimein;
  unsigned long ct = 0;
  TReal         thislength = 0, euclideandist = 0;
  bool          timedone = false;
  VectorType    disp;
  TReal         deltaTime = dT, vecsign = 1.0;
  typename ImageType::SpacingType spacing = this->m_DisplacementField->GetSpacing();
  if( starttimein  > finishtimein )
    {
    vecsign = -1.0;
    }
  VIndexType vind;
  vind.Fill(0);
  xPointType pointIn1;
  pointIn1[TDimension] = 0; // Note: Set the highest value to zero;
  xPointType pointIn2;
  pointIn2[TDimension] = 0; // Note: Set the highest value to zero;
  xPointType pointIn3;
  pointIn3[TDimension] = 0; // Note: Set the highest value to zero;
  for( unsigned int jj = 0; jj < TDimension; jj++ )
    {
    vind[jj] = velind[jj];
    pointIn1[jj] = velind[jj] * spacing[jj];
    }
  this->m_TimeVaryingVelocity->TransformIndexToPhysicalPoint( vind, pointIn1);
// time is in [0,1]
  pointIn1[TDimension] = starttimein * (m_NumberOfTimePoints - 1);
  xPointType Y1x;
  xPointType Y2x;
  xPointType Y3x;
  xPointType Y4x;
  // set up parameters for start of integration
  disp.Fill(0.0);
  timedone = false;
  itime = starttimein;
  ct = 0;
  while( !timedone )
    {
    TReal itimetn1 = itime - timesign * deltaTime;
    TReal itimetn1h = itime - timesign * deltaTime * 0.5;
    if( itimetn1h < 0 )
      {
      itimetn1h = 0;
      }
    if( itimetn1h > 1 )
      {
      itimetn1h = 1;
      }
    if( itimetn1 < 0 )
      {
      itimetn1 = 0;
      }
    if( itimetn1 > 1 )
      {
      itimetn1 = 1;
      }

    TReal totalmag = 0;
    // first get current position of particle
    typename VelocityFieldInterpolatorType::OutputType f1;  f1.Fill(0);
    typename VelocityFieldInterpolatorType::OutputType f2;  f2.Fill(0);
    typename VelocityFieldInterpolatorType::OutputType f3;  f3.Fill(0);
    typename VelocityFieldInterpolatorType::OutputType f4;  f4.Fill(0);
    for( unsigned int jj = 0; jj < TDimension; jj++ )
      {
      pointIn2[jj] = disp[jj] + pointIn1[jj];
      Y1x[jj] = pointIn2[jj];
      Y2x[jj] = pointIn2[jj];
      Y3x[jj] = pointIn2[jj];
      Y4x[jj] = pointIn2[jj];
      }
    if( this->m_Debug )
      {
      ::ants::antscout << " p2 " << pointIn2 << std::endl;
      }

    Y1x[TDimension] = itimetn1 * (TReal)(m_NumberOfTimePoints - 1);
    Y2x[TDimension] = itimetn1h * (TReal)(m_NumberOfTimePoints - 1);
    Y3x[TDimension] = itimetn1h * (TReal)(m_NumberOfTimePoints - 1);
    Y4x[TDimension] = itime * (TReal)(m_NumberOfTimePoints - 1);

    if( this->m_Debug )
      {
      ::ants::antscout << " p2 " << pointIn2 << " y1 " <<  Y1x[TDimension] <<  " y4 " <<   Y4x[TDimension]
                       << std::endl;
      }

    if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y1x) )
      {
      f1 = this->m_VelocityFieldInterpolator->Evaluate( Y1x );
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        Y2x[jj] += f1[jj] * deltaTime * 0.5;
        }
      }
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y2x) )
      {
      f2 = this->m_VelocityFieldInterpolator->Evaluate( Y2x );
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        Y3x[jj] += f2[jj] * deltaTime * 0.5;
        }
      }
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y3x) )
      {
      f3 = this->m_VelocityFieldInterpolator->Evaluate( Y3x );
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        Y4x[jj] += f3[jj] * deltaTime;
        }
      }
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y4x) )
      {
      f4 = this->m_VelocityFieldInterpolator->Evaluate( Y4x );
      }
    for( unsigned int jj = 0; jj < TDimension; jj++ )
      {
      pointIn3[jj] = pointIn2[jj] + vecsign * deltaTime / 6.0 * ( f1[jj] + 2.0 * f2[jj] + 2.0 * f3[jj] + f4[jj] );
      }
    pointIn3[TDimension] = itime * (TReal)(m_NumberOfTimePoints - 1);

    VectorType out;
    TReal      mag = 0, dmag = 0;
    for( unsigned int jj = 0; jj < TDimension; jj++ )
      {
      out[jj] = pointIn3[jj] - pointIn1[jj];
      mag += (pointIn3[jj] - pointIn2[jj]) * (pointIn3[jj] - pointIn2[jj]);
      dmag += (pointIn3[jj] - pointIn1[jj]) * (pointIn3[jj] - pointIn1[jj]);
      disp[jj] = out[jj];
      }

//      ::ants::antscout << " p3 " << pointIn3 << std::endl;
    dmag = sqrt(dmag);
    totalmag += sqrt(mag);
    ct++;
    thislength += totalmag;
    euclideandist = dmag;
    itime = itime + deltaTime * timesign;
    if( starttimein > finishtimein )
      {
      if( itime <= finishtimein  )
        {
        timedone = true;
        }
      }
    else if( thislength ==  0 )
      {
      timedone = true;
      }
    else
      {
      if( itime >= finishtimein )
        {
        timedone = true;
        }
      }
    }

  // now we have the thickness value stored in thislength
  if( this->m_ThickImage && this->m_HitImage )
    {
    // set up parameters for start of integration
    velo.Fill(0);
    itime = starttimein;
    timedone = false;
    vind.Fill(0);
    for( unsigned int jj = 0; jj < TDimension; jj++ )
      {
      vind[jj] = velind[jj];
      pointIn1[jj] = velind[jj] * spacing[jj];
      }
    this->m_TimeVaryingVelocity->TransformIndexToPhysicalPoint( vind, pointIn1);
// time is in [0,1]
    pointIn1[TDimension] = starttimein * (m_NumberOfTimePoints - 1);
    // set up parameters for start of integration
    disp.Fill(0.0);
    timedone = false;
    itime = starttimein;
    ct = 0;
    while( !timedone )
      {
      TReal itimetn1 = itime - timesign * deltaTime;
      TReal itimetn1h = itime - timesign * deltaTime * 0.5;
      if( itimetn1h < 0 )
        {
        itimetn1h = 0;
        }
      if( itimetn1h > 1 )
        {
        itimetn1h = 1;
        }
      if( itimetn1 < 0 )
        {
        itimetn1 = 0;
        }
      if( itimetn1 > 1 )
        {
        itimetn1 = 1;
        }

      //      TReal totalmag=0;
      // first get current position of particle
      typename VelocityFieldInterpolatorType::OutputType f1;  f1.Fill(0);
      typename VelocityFieldInterpolatorType::OutputType f2;  f2.Fill(0);
      typename VelocityFieldInterpolatorType::OutputType f3;  f3.Fill(0);
      typename VelocityFieldInterpolatorType::OutputType f4;  f4.Fill(0);
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        pointIn2[jj] = disp[jj] + pointIn1[jj];
        Y1x[jj] = pointIn2[jj];
        Y2x[jj] = pointIn2[jj];
        Y3x[jj] = pointIn2[jj];
        Y4x[jj] = pointIn2[jj];
        }
      if( this->m_Debug )
        {
        ::ants::antscout << " p2 " << pointIn2 << std::endl;
        }

      Y1x[TDimension] = itimetn1 * (TReal)(m_NumberOfTimePoints - 1);
      Y2x[TDimension] = itimetn1h * (TReal)(m_NumberOfTimePoints - 1);
      Y3x[TDimension] = itimetn1h * (TReal)(m_NumberOfTimePoints - 1);
      Y4x[TDimension] = itime * (TReal)(m_NumberOfTimePoints - 1);
      if( this->m_Debug )
        {
        ::ants::antscout << " p2 " << pointIn2 << " y1 " <<  Y1x[TDimension] <<  " y4 " <<   Y4x[TDimension]
                         << std::endl;
        }

      if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y1x) )
        {
        f1 = this->m_VelocityFieldInterpolator->Evaluate( Y1x );
        for( unsigned int jj = 0; jj < TDimension; jj++ )
          {
          Y2x[jj] += f1[jj] * deltaTime * 0.5;
          }
        }
      if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y2x) )
        {
        f2 = this->m_VelocityFieldInterpolator->Evaluate( Y2x );
        for( unsigned int jj = 0; jj < TDimension; jj++ )
          {
          Y3x[jj] += f2[jj] * deltaTime * 0.5;
          }
        }
      if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y3x) )
        {
        f3 = this->m_VelocityFieldInterpolator->Evaluate( Y3x );
        for( unsigned int jj = 0; jj < TDimension; jj++ )
          {
          Y4x[jj] += f3[jj] * deltaTime;
          }
        }
      if( this->m_VelocityFieldInterpolator->IsInsideBuffer(Y4x) )
        {
        f4 = this->m_VelocityFieldInterpolator->Evaluate( Y4x );
        }
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        pointIn3[jj] = pointIn2[jj] + vecsign * deltaTime / 6.0 * ( f1[jj] + 2.0 * f2[jj] + 2.0 * f3[jj] + f4[jj] );
        }
      pointIn3[TDimension] = itime * (TReal)(m_NumberOfTimePoints - 1);

      VectorType out;
      TReal      mag = 0, dmag = 0;
      for( unsigned int jj = 0; jj < TDimension; jj++ )
        {
        out[jj] = pointIn3[jj] - pointIn1[jj];
        mag += (pointIn3[jj] - pointIn2[jj]) * (pointIn3[jj] - pointIn2[jj]);
        dmag += (pointIn3[jj] - pointIn1[jj]) * (pointIn3[jj] - pointIn1[jj]);
        disp[jj] = out[jj];
        }
      itime = itime + deltaTime * timesign;
      if( starttimein > finishtimein )
        {
        if( itime <= finishtimein  )
          {
          timedone = true;
          }
        }
      else
        {
        if( itime >= finishtimein )
          {
          timedone = true;
          }
        }

      //        bool isingray=true;
      if( this->m_MaskImage )
        {
        if( this->m_MaskImage->GetPixel(velind) )
          {
          VIndexType thind2;
          IndexType  thind;
          bool       isin = this->m_TimeVaryingVelocity->TransformPhysicalPointToIndex( pointIn3, thind2 );
          for( unsigned int ij = 0; ij < ImageDimension; ij++ )
            {
            thind[ij] = thind2[ij];
            }
          if( isin )
            {
            unsigned long lastct = (unsigned long) this->m_HitImage->GetPixel(thind);
            unsigned long newct = lastct + 1;
            TReal         oldthick = this->m_ThickImage->GetPixel(thind);
            TReal         newthick = (TReal)lastct / (TReal)newct * oldthick + 1.0 / (TReal)newct * euclideandist;
            this->m_HitImage->SetPixel( thind,  newct );
            this->m_ThickImage->SetPixel(thind, newthick );
            }
          else
            {
            ::ants::antscout << " thind " << thind << " edist " << euclideandist << " p3 " << pointIn3 << " p1 "
                             << pointIn1 << std::endl;
            }
          //        this->m_ThickImage->SetPixel(thind, thislength );
          }
        }
      }
    }

//  if (!isinside) { ::ants::antscout << " velind " << velind << " not inside " << Y1 << std::endl;   }

  if( this->m_Debug )
    {
    ::ants::antscout << " Length " << thislength << std::endl;
    }
  this->m_Debug = false;
  return disp;
}

/**
 * Standard "PrintSelf" method
 */
template <unsigned int TDimension, class TReal>
void
ANTSImageRegistrationOptimizer<TDimension, TReal>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk
#endif
