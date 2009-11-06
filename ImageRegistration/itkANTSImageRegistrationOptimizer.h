/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkANTSImageRegistrationOptimizer.h,v $
 Language:  C++
 Date:      $Date: 2009/04/22 01:00:17 $
 Version:   $Revision: 1.44 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkANTSImageRegistrationOptimizer_h
#define __itkANTSImageRegistrationOptimizer_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkCommandLineParser.h"
#include "itkShiftScaleImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "ReadWriteImage.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkANTSAffine3DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkFiniteDifferenceFunction.h"
#include "itkFixedArray.h"
#include "itkANTSSimilarityMetric.h"
#include "itkVectorExpandImageFilter.h"
// #include "itkNeighbohoodAlgorithm.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkDeformationFieldFromMultiTransformFilter.h"
#include "itkWarpImageWAffineFilter.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkGeneralToBSplineDeformationFieldFilter.h"
#include "ANTS_affine_registration2.h"
#include "itkVectorFieldGradientImageFunction.h"

namespace itk
{
template <unsigned int TDimension = 3, class TReal = float>
class ITK_EXPORT ANTSImageRegistrationOptimizer
  : public Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSImageRegistrationOptimizer Self;
  typedef Object                         Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ANTSImageRegistrationOptimizer, Object );
  itkStaticConstMacro( Dimension, unsigned int, TDimension );
  itkStaticConstMacro( ImageDimension, unsigned int, TDimension );

  typedef TReal RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( Dimension )> ImageType;
  typedef typename ImageType::Pointer                                            ImagePointer;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> TransformType;

  /** Point Types  for landmarks and labeled point-sets */
  typedef itk::ANTSLabeledPointSet<Dimension>        LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer      LabeledPointSetPointer;
  typedef typename LabeledPointSetType::PointSetType PointSetType;
  typedef typename PointSetType::Pointer             PointSetPointer;
  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PixelType           PointDataType;
  typedef typename ImageType::PointType              ImagePointType;

  typedef itk::MatrixOffsetTransformBase<double, TDimension, TDimension> AffineTransformType;
  typedef typename AffineTransformType::Pointer                          AffineTransformPointer;
  typedef OptAffine<AffineTransformPointer, ImagePointer>                OptAffineType;

  typedef itk::Vector<float, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef typename DeformationFieldType::Pointer DeformationFieldPointer;

  typedef itk::Image<VectorType, ImageDimension + 1> TimeVaryingVelocityFieldType;
  typedef typename TimeVaryingVelocityFieldType::Pointer
  TimeVaryingVelocityFieldPointer;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float> VelocityFieldInterpolatorType;
  typedef typename DeformationFieldType::IndexType                                       IndexType;

  typedef CommandLineParser               ParserType;
  typedef typename ParserType::OptionType OptionType;

  typedef GeneralToBSplineDeformationFieldFilter<DeformationFieldType> BSplineFilterType;
  typedef FixedArray<RealType,
                     itkGetStaticConstMacro( ImageDimension )>                          ArrayType;

  /** Typedefs for similarity metrics */
  typedef ANTSSimilarityMetric<itkGetStaticConstMacro( Dimension ), float> SimilarityMetricType;
  typedef typename SimilarityMetricType::Pointer                           SimilarityMetricPointer;
  typedef std::vector<SimilarityMetricPointer>                             SimilarityMetricListType;

  /** FiniteDifferenceFunction type. */
  typedef FiniteDifferenceFunction<DeformationFieldType>      FiniteDifferenceFunctionType;
  typedef typename FiniteDifferenceFunctionType::TimeStepType TimeStepType;
  typedef typename
  FiniteDifferenceFunctionType::Pointer FiniteDifferenceFunctionPointer;
  typedef AvantsPDEDeformableRegistrationFunction<ImageType, ImageType,
                                                  DeformationFieldType> MetricBaseType;
  typedef typename MetricBaseType::Pointer MetricBaseTypePointer;

  /* Jacobian and other calculations */
  typedef itk::VectorFieldGradientImageFunction<DeformationFieldType> JacobianFunctionType;

  /** Set functions */
  void SetAffineTransform(AffineTransformPointer A)
  {
    this->m_AffineTransform = A;
  }

  void SetDeformationField(DeformationFieldPointer A)
  {
    this->m_DeformationField = A;
  }

  void SetInverseDeformationField(DeformationFieldPointer A)
  {
    this->m_InverseDeformationField = A;
  }

  void SetMaskImage( ImagePointer m)
  {
    this->m_MaskImage = m;
  }

  void SetFixedImageAffineTransform(AffineTransformPointer A)
  {
    this->m_FixedImageAffineTransform = A;
  }

  AffineTransformPointer GetFixedImageAffineTransform()
  {
    return this->m_FixedImageAffineTransform;
  }

  /** Get functions */
  AffineTransformPointer GetAffineTransform()
  {
    return this->m_AffineTransform;
  }

  DeformationFieldPointer GetDeformationField()
  {
    return this->m_DeformationField;
  }

  DeformationFieldPointer GetInverseDeformationField()
  {
    return this->m_InverseDeformationField;
  }

  /** Initialize all parameters */

  void SetNumberOfLevels(unsigned int i)
  {
    this->m_NumberOfLevels = i;
  }

  void SetParser( typename ParserType::Pointer P )
  {
    this->m_Parser = P;
  }

  /** Basic operations */
  DeformationFieldPointer CopyDeformationField( DeformationFieldPointer input );

  void SmoothDeformationField(DeformationFieldPointer field,  bool TrueEqualsGradElseTotal )
  {
    typename ParserType::OptionType::Pointer regularizationOption
      = this->m_Parser->GetOption( "regularization" );

    if( ( regularizationOption->GetValue() ).find( "DMFFD" )
        != std::string::npos )
      {
      if( ( !TrueEqualsGradElseTotal && this->m_TotalSmoothingparam == 0.0 ) ||
          ( TrueEqualsGradElseTotal && this->m_GradSmoothingparam == 0.0 ) )
        {
        return;
        }
      ArrayType    meshSize;
      unsigned int splineOrder = this->m_BSplineFieldOrder;
      float        bsplineKernelVariance = static_cast<float>( splineOrder + 1 ) / 12.0;
      unsigned int numberOfLevels = 1;

      if( TrueEqualsGradElseTotal )
        {
        if( this->m_GradSmoothingparam < 0.0 )
          {
          meshSize = this->m_GradSmoothingMeshSize;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            meshSize[d] *= static_cast<unsigned int>(
                vcl_pow( 2.0, static_cast<int>( this->m_CurrentLevel ) ) );
            }
          }
        else
          {
          float spanLength = vcl_sqrt( this->m_GradSmoothingparam
                                       / bsplineKernelVariance );
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            meshSize[d] = static_cast<unsigned int>(
                field->GetLargestPossibleRegion().GetSize()[d]
                / spanLength + 0.5 );
            }
          }
        this->SmoothDeformationFieldBSpline( field, meshSize, splineOrder,
                                             numberOfLevels );
        }
      else
        {
        if( this->m_TotalSmoothingparam < 0.0 )
          {
          meshSize = this->m_TotalSmoothingMeshSize;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            meshSize[d] *= static_cast<unsigned int>(
                vcl_pow( 2.0, static_cast<int>( this->m_CurrentLevel ) ) );
            }
          }
        else
          {
          float spanLength = vcl_sqrt( this->m_TotalSmoothingparam
                                       / bsplineKernelVariance );
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            meshSize[d] = static_cast<unsigned int>(
                field->GetLargestPossibleRegion().GetSize()[d]
                / spanLength + 0.5 );
            }
          }

        RealType maxMagnitude = 0.0;

        ImageRegionIterator<DeformationFieldType> It( field,
                                                      field->GetLargestPossibleRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          RealType magnitude = ( It.Get() ).GetNorm();
          if( magnitude > maxMagnitude )
            {
            maxMagnitude = magnitude;
            }
          }
        this->SmoothDeformationFieldBSpline( field, meshSize, splineOrder,
                                             numberOfLevels );

        if( maxMagnitude > 0.0 )
          {
          for( It.GoToBegin(); !It.IsAtEnd(); ++It )
            {
            It.Set( It.Get() / maxMagnitude );
            }
          }
        }
      }
    else // Gaussian
      {
      float sig = 0;
      if( TrueEqualsGradElseTotal )
        {
        sig = this->m_GradSmoothingparam;
        }
      else
        {
        sig = this->m_TotalSmoothingparam;
        }
      this->SmoothDeformationFieldGauss(field, sig);
      }
  }

  void SmoothDeformationFieldGauss(DeformationFieldPointer field = NULL, float sig = 0.0, bool useparamimage = false,
                                   unsigned int lodim = ImageDimension);

//  float = smoothingparam, int = maxdim to smooth
  void SmoothVelocityGauss(TimeVaryingVelocityFieldPointer field, float, unsigned int);

  void SmoothDeformationFieldBSpline(DeformationFieldPointer field, ArrayType meshSize, unsigned int splineorder,
                                     unsigned int numberoflevels );

  DeformationFieldPointer ComputeUpdateFieldAlternatingMin(DeformationFieldPointer fixedwarp,
                                                           DeformationFieldPointer movingwarp,
                                                           PointSetPointer  fpoints = NULL,  PointSetPointer wpoints =
                                                             NULL,
                                                           DeformationFieldPointer updateFieldInv = NULL,
                                                           bool updateenergy = true);

  DeformationFieldPointer ComputeUpdateField(DeformationFieldPointer fixedwarp, DeformationFieldPointer movingwarp,
                                             PointSetPointer  fpoints = NULL,  PointSetPointer wpoints = NULL,
                                             DeformationFieldPointer updateFieldInv = NULL,
                                             bool updateenergy = true);

  TimeVaryingVelocityFieldPointer ExpandVelocity()
  {
    float expandFactors[ImageDimension + 1];

    expandFactors[ImageDimension] = 1;
    m_Debug = false;
    for( int idim = 0; idim < ImageDimension; idim++ )
      {
      expandFactors[idim] = (float)this->m_CurrentDomainSize[idim]
        / (float) this->m_TimeVaryingVelocity->GetLargestPossibleRegion().GetSize()[idim];
      if( expandFactors[idim] < 1 )
        {
        expandFactors[idim] = 1;
        }
      if( this->m_Debug )
        {
        std::cout << " ExpFac " << expandFactors[idim] << " curdsz " << this->m_CurrentDomainSize[idim] << std::endl;
        }
      }
    VectorType pad;  pad.Fill(0);
    typedef VectorExpandImageFilter<TimeVaryingVelocityFieldType, TimeVaryingVelocityFieldType> ExpanderType;
    typename ExpanderType::Pointer m_FieldExpander = ExpanderType::New();
    m_FieldExpander->SetInput(this->m_TimeVaryingVelocity);
    m_FieldExpander->SetExpandFactors( expandFactors );
    m_FieldExpander->SetEdgePaddingValue( pad );
    m_FieldExpander->UpdateLargestPossibleRegion();
    return m_FieldExpander->GetOutput();
  }

  DeformationFieldPointer ExpandField(DeformationFieldPointer field,  typename ImageType::SpacingType targetSpacing)
  {
//      this->m_Debug=true;
    float expandFactors[ImageDimension];

    for( int idim = 0; idim < ImageDimension; idim++ )
      {
      expandFactors[idim] = (float)this->m_CurrentDomainSize[idim]
        / (float)field->GetLargestPossibleRegion().GetSize()[idim];
      if( expandFactors[idim] < 1 )
        {
        expandFactors[idim] = 1;
        }
      //             if (this->m_Debug)  std::cout << " ExpFac " << expandFactors[idim] << " curdsz " <<
      // this->m_CurrentDomainSize[idim] << std::endl;
      }

    VectorType pad;
    pad.Fill(0);
    typedef VectorExpandImageFilter<DeformationFieldType, DeformationFieldType> ExpanderType;
    typename ExpanderType::Pointer m_FieldExpander = ExpanderType::New();
    m_FieldExpander->SetInput(field);
    m_FieldExpander->SetExpandFactors( expandFactors );
    // use default
    m_FieldExpander->SetEdgePaddingValue( pad );
    m_FieldExpander->UpdateLargestPossibleRegion();

    typename DeformationFieldType::Pointer fieldout = m_FieldExpander->GetOutput();
    fieldout->SetSpacing(targetSpacing);
    fieldout->SetOrigin(field->GetOrigin() );
    if( this->m_Debug )
      {
      std::cout << " Field size " << fieldout->GetLargestPossibleRegion().GetSize() << std::endl;
      }
    // this->m_Debug=false;

    return fieldout;
  }

  ImagePointer GetVectorComponent(DeformationFieldPointer field, unsigned int index)
  {
    // Initialize the Moving to the displacement field
    typedef DeformationFieldType FieldType;

    typename ImageType::Pointer sfield = ImageType::New();
    sfield->SetSpacing( field->GetSpacing() );
    sfield->SetOrigin( field->GetOrigin() );
    sfield->SetDirection( field->GetDirection() );
    sfield->SetLargestPossibleRegion(field->GetLargestPossibleRegion() );
    sfield->SetRequestedRegion(field->GetRequestedRegion() );
    sfield->SetBufferedRegion( field->GetBufferedRegion() );
    sfield->Allocate();

    typedef itk::ImageRegionIteratorWithIndex<FieldType> Iterator;
    Iterator vfIter( field, field->GetLargestPossibleRegion() );
    for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      VectorType v1 = vfIter.Get();
      sfield->SetPixel(vfIter.GetIndex(), v1[index]);
      }

    return sfield;
  }

  ImagePointer SubsampleImage( ImagePointer, RealType, typename ImageType::PointType outputOrigin,
                               typename ImageType::DirectionType outputDirection,
                               AffineTransformPointer aff = NULL);

  DeformationFieldPointer SubsampleField( DeformationFieldPointer field, typename ImageType::SizeType
                                          targetSize, typename ImageType::SpacingType targetSpacing )
  {
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << "FIXME -- NOT DONE CORRECTLY " << std::endl;
    std::cout << " SUBSAM FIELD SUBSAM FIELD SUBSAM FIELD " << std::endl;
    typename DeformationFieldType::Pointer sfield = DeformationFieldType::New();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      typename ImageType::Pointer precomp = this->GetVectorComponent(field, i);
      typename ImageType::Pointer comp = this->SubsampleImage(precomp, targetSize, targetSpacing);
      if( i == 0 )
        {
        sfield->SetSpacing( comp->GetSpacing() );
        sfield->SetOrigin( comp->GetOrigin() );
        sfield->SetDirection( comp->GetDirection() );
        sfield->SetLargestPossibleRegion(comp->GetLargestPossibleRegion() );
        sfield->SetRequestedRegion(comp->GetRequestedRegion() );
        sfield->SetBufferedRegion( comp->GetBufferedRegion() );
        sfield->Allocate();
        }

      typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> Iterator;
      typedef typename DeformationFieldType::PixelType                VectorType;
      VectorType v1;
      VectorType zero;
      zero.Fill(0.0);
      Iterator vfIter( sfield, sfield->GetLargestPossibleRegion() );
      for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        v1 = vfIter.Get();
        v1[i] = comp->GetPixel(vfIter.GetIndex() );
        vfIter.Set(v1);
        }
      }

    return sfield;
  }

  PointSetPointer  WarpMultiTransform(ImagePointer referenceimage, ImagePointer movingImage,
                                      PointSetPointer movingpoints,  AffineTransformPointer aff,
                                      DeformationFieldPointer totalField, bool doinverse,
                                      AffineTransformPointer  fixedaff  )
  {
    if( !movingpoints )
      {
      std::cout << " NULL POINTS " << std::endl;  return NULL;
      }

    AffineTransformPointer affinverse = NULL;
    if( aff )
      {
      affinverse = AffineTransformType::New();
      aff->GetInverse(affinverse);
      }
    AffineTransformPointer fixedaffinverse = NULL;
    if( fixedaff )
      {
      fixedaffinverse = AffineTransformType::New();
      fixedaff->GetInverse(fixedaffinverse);
      }

    typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DeformationFieldType, TransformType> WarperType;
    typename WarperType::Pointer  warper = WarperType::New();
    warper->SetInput(movingImage);
    warper->SetEdgePaddingValue( 0);
    warper->SetSmoothScale(1);
    if( !doinverse )
      {
      if( totalField )
        {
        warper->PushBackDeformationFieldTransform(totalField);
        }
      if( fixedaff )
        {
        warper->PushBackAffineTransform(fixedaff);
        }
      else if( aff )
        {
        warper->PushBackAffineTransform(aff);
        }
      }
    else
      {
      if( aff )
        {
        warper->PushBackAffineTransform( affinverse );
        }
      else if( fixedaff )
        {
        warper->PushBackAffineTransform(fixedaffinverse);
        }
      if( totalField )
        {
        warper->PushBackDeformationFieldTransform(totalField);
        }
      }

    warper->SetOutputOrigin(referenceimage->GetOrigin() );
    typename ImageType::SizeType size = referenceimage->GetLargestPossibleRegion().GetSize();
    if( totalField )
      {
      size = totalField->GetLargestPossibleRegion().GetSize();
      }
    warper->SetOutputSize(size);
    typename ImageType::SpacingType spacing = referenceimage->GetSpacing();
    if( totalField )
      {
      spacing = totalField->GetSpacing();
      }
    warper->SetOutputSpacing(spacing);
    warper->SetOutputDirection(referenceimage->GetDirection() );
    totalField->SetOrigin(referenceimage->GetOrigin() );
    totalField->SetDirection(referenceimage->GetDirection() );

    // warper->Update();
//      std::cout << " updated in point warp " << std::endl;
    PointSetPointer outputMesh = PointSetType::New();
    unsigned long   count = 0;
    unsigned long   sz1 = movingpoints->GetNumberOfPoints();
    if( this->m_Debug )
      {
      std::cout << " BEFORE #points " << sz1 << std::endl;
      }
    for( unsigned long ii = 0; ii < sz1; ii++ )
      {
      PointType     point, wpoint;
      PointDataType label = 0;
      movingpoints->GetPoint(ii, &point);
      movingpoints->GetPointData(ii, &label);
// convert pointtype to imagepointtype
      ImagePointType pt, wpt;
      for( unsigned int jj = 0;  jj < ImageDimension; jj++ )
        {
        pt[jj] = point[jj];
        }
      bool bisinside = warper->MultiTransformSinglePoint(pt, wpt);
      if( bisinside )
        {
        for( unsigned int jj = 0;  jj < ImageDimension; jj++ )
          {
          wpoint[jj] = wpt[jj];
          }
        outputMesh->SetPointData( count, label );
        outputMesh->SetPoint( count, wpoint );
//	  if (ii % 100 == 0) std::cout << " pt " << pt << " wpt " << wpt << std::endl;
        count++;
        }
      }
    if( this->m_Debug )
      {
      std::cout << " AFTER #points " << count << std::endl;
      }
//      if (count != sz1 ) std::cout << " WARNING:  POINTS ARE MAPPING OUT OF IMAGE DOMAIN " << 1.0 - (float)
// count/(float)(sz1+1) << std::endl;
    return outputMesh;
  }

  ImagePointer WarpMultiTransform( ImagePointer referenceimage,  ImagePointer movingImage,  AffineTransformPointer aff,
                                   DeformationFieldPointer totalField, bool doinverse,
                                   AffineTransformPointer  fixedaff  )
  {
    typedef typename ImageType::DirectionType DirectionType;
    DirectionType rdirection = referenceimage->GetDirection();
    DirectionType mdirection = movingImage->GetDirection();

    AffineTransformPointer affinverse = NULL;
    if( aff )
      {
      affinverse = AffineTransformType::New();
      aff->GetInverse(affinverse);
      }
    AffineTransformPointer fixedaffinverse = NULL;
    if( fixedaff )
      {
      fixedaffinverse = AffineTransformType::New();
      fixedaff->GetInverse(fixedaffinverse);
      }

    DirectionType iddir;
    iddir.Fill(0);
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      iddir[i][i] = 1;
      }

    typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
    typename InterpolatorType1::Pointer interp1 = InterpolatorType1::New();
    typename InterpolatorType2::Pointer interpnn = InterpolatorType2::New();

    this->m_UseMulti = true;

    if( !this->m_UseMulti )
      {
      ImagePointer wmimage = this->SubsampleImage(movingImage, this->m_ScaleFactor,
                                                  movingImage->GetOrigin(), movingImage->GetDirection(), aff );
      typedef itk::WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;
      typename WarperType::Pointer  warper;
      warper = WarperType::New();
      warper->SetInput( wmimage);
      warper->SetDeformationField(totalField);
      warper->SetOutputSpacing(totalField->GetSpacing() );
      warper->SetOutputOrigin(totalField->GetOrigin() );
      warper->SetInterpolator(interp1);
      if( this->m_UseNN )
        {
        warper->SetInterpolator(interpnn);
        }
//    warper->SetOutputSize(this->m_CurrentDomainSize);
//    warper->SetEdgePaddingValue( 0 );
      warper->Update();
      return warper->GetOutput();
      }

    typedef itk::WarpImageMultiTransformFilter<ImageType, ImageType, DeformationFieldType, TransformType> WarperType;
    typename WarperType::Pointer  warper = WarperType::New();
    warper->SetInput(movingImage);
    warper->SetEdgePaddingValue( 0);
    warper->SetSmoothScale(1);
    warper->SetInterpolator(interp1);
    if( this->m_UseNN )
      {
      warper->SetInterpolator(interpnn);
      }
    if( !doinverse )
      {
      if( totalField )
        {
        warper->PushBackDeformationFieldTransform(totalField);
        }
      if( fixedaff )
        {
        warper->PushBackAffineTransform(fixedaff);
        }
      else if( aff )
        {
        warper->PushBackAffineTransform(aff);
        }
      }
    else
      {
      if( aff )
        {
        warper->PushBackAffineTransform( affinverse );
        }
      else if( fixedaff )
        {
        warper->PushBackAffineTransform(fixedaffinverse);
        }
      if( totalField )
        {
        warper->PushBackDeformationFieldTransform(totalField);
        }
      }

    warper->SetOutputOrigin(referenceimage->GetOrigin() );
    typename ImageType::SizeType size = referenceimage->GetLargestPossibleRegion().GetSize();
    if( totalField )
      {
      size = totalField->GetLargestPossibleRegion().GetSize();
      }
    warper->SetOutputSize(size);
    typename ImageType::SpacingType spacing = referenceimage->GetSpacing();
    if( totalField )
      {
      spacing = totalField->GetSpacing();
      }
    warper->SetOutputSpacing(spacing);
    warper->SetOutputDirection(referenceimage->GetDirection() );
    totalField->SetOrigin(referenceimage->GetOrigin() );
    totalField->SetDirection(referenceimage->GetDirection() );

    warper->Update();
    if( this->m_Debug )
      {
      std::cout << " updated ok -- warped image output size "
                << warper->GetOutput()->GetLargestPossibleRegion().GetSize() << " requested size "
                <<  totalField->GetLargestPossibleRegion().GetSize() <<  std::endl;
      }

    typename ImageType::Pointer outimg = warper->GetOutput();

    return outimg;
  }

  ImagePointer  SmoothImageToScale(ImagePointer image,  float scalingFactor )
  {
    typename ImageType::SpacingType inputSpacing = image->GetSpacing();
    typename ImageType::RegionType::SizeType inputSize = image->GetRequestedRegion().GetSize();

    typename ImageType::SpacingType outputSpacing;
    typename ImageType::RegionType::SizeType outputSize;

    RealType minimumSpacing = inputSpacing.GetVnlVector().min_value();
//    RealType maximumSpacing = inputSpacing.GetVnlVector().max_value();
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      RealType scaling = vnl_math_min( scalingFactor * minimumSpacing / inputSpacing[d],
                                       static_cast<RealType>( inputSize[d] ) / 32.0 );
      outputSpacing[d] = inputSpacing[d] * scaling;
      outputSize[d] = static_cast<unsigned long>( inputSpacing[d]
                                                  * static_cast<RealType>( inputSize[d] ) / outputSpacing[d] + 0.5 );

      typedef RecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
      typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
      smoother->SetInputImage( image );
      smoother->SetDirection( d );
      smoother->SetNormalizeAcrossScale( false );
      float sig = (outputSpacing[d] / inputSpacing[d] - 1.0) * 0.2; // /(float)ImageDimension;
      smoother->SetSigma(sig );

      if( smoother->GetSigma() > 0.0 )
        {
        smoother->Update();
        image = smoother->GetOutput();
        }
      }

    image = this->NormalizeImage(image);

    return image;
  }

  typename ANTSImageRegistrationOptimizer<TDimension, TReal>::DeformationFieldPointer IntegrateConstantVelocity(
    DeformationFieldPointer totalField, unsigned int ntimesteps, float timeweight);

  /** Base optimization functions */
  // AffineTransformPointer AffineOptimization(AffineTransformPointer &aff_init, OptAffine &affine_opt); // {return
  // NULL;}
  AffineTransformPointer AffineOptimization(OptAffineType & affine_opt);  // {return NULL;}

  std::string GetTransformationModel()
  {
    return this->m_TransformationModel;
  }

  void SetTransformationModel( std::string  s)
  {
    this->m_TransformationModel = s;
    std::cout << " Requested Transformation Model:  " << this->m_TransformationModel << " : Using " << std::endl;
    if( this->m_TransformationModel  == std::string("Elast") )
      {
      std::cout << "Elastic model for transformation. " << std::endl;
      }
    else if( this->m_TransformationModel  == std::string("SyN") )
      {
      std::cout << "SyN diffeomorphic model for transformation. " << std::endl;
      }
    else
      {
      std::cout << "Exp Diff model for transformation. " << std::endl;
      this->m_TransformationModel = std::string("Exp");
      }
  }

  void SetUpParameters()
  {
    /** Univariate Deformable Mapping */
//    this->SetTransformationModel(this->m_Parser->GetOption( "transformation-model" )->GetValue());
//    std::string temp=this->m_Parser->GetOption( "number-of-iterations" )->GetValue();
//    this->m_Iterations = this->m_Parser->template ConvertVector<unsigned int>(temp);
//    this->SetNumberOfLevels(this->m_Iterations.size());
//    this->m_UseROI=false;
//    if ( typename OptionType::Pointer option = this->m_Parser->GetOption( "r" ) )
//      {
//      temp=this->m_Parser->GetOption( "r" )->GetValue();
//      this->m_RoiNumbers = this->m_Parser->template ConvertVector<float>(temp);
//      if (temp.length() > 3) this->m_UseROI=true;
//      }
//
//    temp=this->m_Parser->GetOption( "gradient-step-length" )->GetValue();
//    this->m_Gradstep = this->m_Parser->template Convert<float>( temp );
//    temp=this->m_Parser->GetOption( "def-field-sigma")->GetValue();
//    this->m_TotalSmoothingparam=this->m_Parser->template Convert<float>( temp );
//    temp=this->m_Parser->GetOption( "gradient-field-sigma")->GetValue();
//    this->m_GradSmoothingparam=this->m_Parser->template Convert<float>( temp );
//    std::cout <<"  Grad Step " << this->m_Gradstep << " total-smoothing " << this->m_TotalSmoothingparam << "
// gradient-smoothing " << this->m_GradSmoothingparam << std::endl;

    std::string temp = this->m_Parser->GetOption( "number-of-iterations" )->GetValue();

    this->m_Iterations = this->m_Parser->template ConvertVector<unsigned int>(temp);
    this->SetNumberOfLevels(this->m_Iterations.size() );
    this->m_UseROI = false;
    if( typename OptionType::Pointer option = this->m_Parser->GetOption( "roi" ) )
      {
      temp = this->m_Parser->GetOption( "roi" )->GetValue();
      this->m_RoiNumbers = this->m_Parser->template ConvertVector<float>(temp);
      if( temp.length() > 3 )
        {
        this->m_UseROI = true;
        }
      }

    typename ParserType::OptionType::Pointer oOption
      = this->m_Parser->GetOption( "output-naming" );
    this->m_OutputNamingConvention = oOption->GetValue();

    typename ParserType::OptionType::Pointer thicknessOption
      = this->m_Parser->GetOption( "compute-thickness" );
    if( thicknessOption->GetValue() == "true" ||  thicknessOption->GetValue() == "1" )
      {
      this->m_ComputeThickness = 1; this->m_SyNFullTime = 2;
      }                                                                                                                                      //
                                                                                                                                             // asymm
                                                                                                                                             // forces
    else if(  thicknessOption->GetValue() == "2" )
      {
      this->m_ComputeThickness = 1; this->m_SyNFullTime = 1;
      }                                                                                                    // symmetric
                                                                                                           // forces
    else
      {
      this->m_ComputeThickness = 0;  // not full time varying stuff
      }
    /**
     * Get transformation model and associated parameters
     */
    typename ParserType::OptionType::Pointer transformOption
      = this->m_Parser->GetOption( "transformation-model" );
    this->SetTransformationModel( transformOption->GetValue() );
    if( transformOption->GetNumberOfParameters() >= 1 )
      {
      std::string parameter = transformOption->GetParameter( 0, 0 );
      float       temp = this->m_Parser->template Convert<float>( parameter );
      this->m_Gradstep = temp;
      this->m_GradstepAltered = temp;
      }
    else
      {
      this->m_Gradstep = 0.5;  this->m_GradstepAltered = 0.5;
      }
    if( transformOption->GetNumberOfParameters() >= 2 )
      {
      std::string parameter = transformOption->GetParameter( 0, 1 );
      this->m_NTimeSteps = this->m_Parser->template Convert<unsigned int>( parameter );
      }
    else
      {
      this->m_NTimeSteps = 1;
      }
    if( transformOption->GetNumberOfParameters() >= 3 )
      {
      std::string parameter = transformOption->GetParameter( 0, 2 );
      this->m_DeltaTime
        = this->m_Parser->template Convert<float>( parameter );
      if( this->m_DeltaTime  > 1 )
        {
        this->m_DeltaTime = 1;
        }
      if( this->m_DeltaTime  <= 0 )
        {
        this->m_DeltaTime = 0.001;
        }
      std::cout << " set DT " << this->m_DeltaTime << std::endl;
      this->m_SyNType = 1;
      }
    else
      {
      this->m_DeltaTime = 0.1;
      }
//    if ( transformOption->GetNumberOfParameters() >= 3 )
//      {
//      std::string parameter = transformOption->GetParameter( 0, 2 );
//      this->m_SymmetryType
//        = this->m_Parser->template Convert<unsigned int>( parameter );
//      }

    /**
     * Get regularization and associated parameters
     */
    this->m_GradSmoothingparam = -1;
    this->m_TotalSmoothingparam = -1;
    this->m_GradSmoothingMeshSize.Fill( 0 );
    this->m_TotalSmoothingMeshSize.Fill( 0 );

    typename ParserType::OptionType::Pointer regularizationOption
      = this->m_Parser->GetOption( "regularization" );
    if( regularizationOption->GetValue() == "Gauss" )
      {
      if( regularizationOption->GetNumberOfParameters() >= 1 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 0 );
        this->m_GradSmoothingparam = this->m_Parser->template Convert<float>( parameter );
        }
      else
        {
        this->m_GradSmoothingparam = 3;
        }
      if( regularizationOption->GetNumberOfParameters() >= 2 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 1 );
        this->m_TotalSmoothingparam = this->m_Parser->template Convert<float>( parameter );
        }
      else
        {
        this->m_TotalSmoothingparam = 0.5;
        }
      if( regularizationOption->GetNumberOfParameters() >= 3 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 2 );
        this->m_GaussianTruncation = this->m_Parser->template Convert<float>( parameter );
        }
      else
        {
        this->m_GaussianTruncation = 256;
        }
      std::cout << "  Grad Step " << this->m_Gradstep << " total-smoothing " << this->m_TotalSmoothingparam
                << " gradient-smoothing " << this->m_GradSmoothingparam << std::endl;
      }
    else if( ( regularizationOption->GetValue() ).find( "DMFFD" )
             != std::string::npos )
      {
      if( regularizationOption->GetNumberOfParameters() >= 1 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 0 );
        if( parameter.find( "x" ) != std::string::npos )
          {
          std::vector<unsigned int> gradMeshSize
            = this->m_Parser->template ConvertVector<unsigned int>( parameter );
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            this->m_GradSmoothingMeshSize[d] = gradMeshSize[d];
            }
          }
        else
          {
          this->m_GradSmoothingparam
            = this->m_Parser->template Convert<float>( parameter );
          }
        }
      else
        {
        this->m_GradSmoothingparam = 3.0;
        }
      if( regularizationOption->GetNumberOfParameters() >= 2 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 1 );
        if( parameter.find( "x" ) != std::string::npos )
          {
          std::vector<unsigned int> totalMeshSize
            = this->m_Parser->template ConvertVector<unsigned int>( parameter );
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            this->m_TotalSmoothingMeshSize[d] = totalMeshSize[d];
            }
          }
        else
          {
          this->m_TotalSmoothingparam
            = this->m_Parser->template Convert<float>( parameter );
          }
        }
      else
        {
        this->m_TotalSmoothingparam = 0.5;
        }
      if( regularizationOption->GetNumberOfParameters() >= 3 )
        {
        std::string parameter = regularizationOption->GetParameter( 0, 2 );
        this->m_BSplineFieldOrder
          = this->m_Parser->template Convert<unsigned int>( parameter );
        }
      else
        {
        this->m_BSplineFieldOrder = 3;
        }
      std::cout << "  Grad Step " << this->m_Gradstep
                << " total-smoothing " << this->m_TotalSmoothingparam
                << " gradient-smoothing " << this->m_GradSmoothingparam
                << " bspline-field-order " << this->m_BSplineFieldOrder
                << std::endl;
      }
    else
      {
      this->m_GradSmoothingparam = 3;
      this->m_TotalSmoothingparam = 0.5;
      std::cout << " Default Regularization is Gaussian smoothing with : " << this->m_GradSmoothingparam  << " & "
                << this->m_TotalSmoothingparam << std::endl;
//      itkExceptionMacro( "Invalid regularization: " << regularizationOption->GetValue() );
      }
  }

  void ComputeMultiResolutionParameters(ImagePointer fixedImage )
  {
    VectorType zero;

    zero.Fill(0);
    /** Compute scale factors */
    this->m_FullDomainSpacing = fixedImage->GetSpacing();
    this->m_FullDomainSize = fixedImage->GetRequestedRegion().GetSize();
    this->m_CurrentDomainSpacing = fixedImage->GetSpacing();
    this->m_CurrentDomainSize = fixedImage->GetRequestedRegion().GetSize();
    this->m_CurrentDomainDirection = fixedImage->GetDirection();
    this->m_FullDomainOrigin.Fill(0);
    this->m_CurrentDomainOrigin.Fill(0);
    /** alter the input size based on information gained from the ROI information - if available */
    if( this->m_UseROI )
      {
      for( unsigned int ii = 0; ii < ImageDimension; ii++ )
        {
        this->m_FullDomainSize[ii] =
          (typename ImageType::SizeType::SizeValueType) this->m_RoiNumbers[ii + ImageDimension];
        this->m_FullDomainOrigin[ii] = this->m_RoiNumbers[ii];
        }
      std::cout << " ROI #s : size " << this->m_FullDomainSize << " orig " <<   this->m_FullDomainOrigin  << std::endl;
      }

    RealType minimumSpacing = this->m_FullDomainSpacing.GetVnlVector().min_value();
//            RealType maximumSpacing = this->m_FullDomainSpacing.GetVnlVector().max_value();
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      RealType scaling = vnl_math_min( this->m_ScaleFactor * minimumSpacing / this->m_FullDomainSpacing[d],
                                       static_cast<RealType>( this->m_FullDomainSize[d] ) / 32.0 );
      if( scaling < 1 )
        {
        scaling = 1;
        }
      this->m_CurrentDomainSpacing[d] = this->m_FullDomainSpacing[d] * scaling;
      this->m_CurrentDomainSize[d] =
        static_cast<unsigned long>( this->m_FullDomainSpacing[d] * static_cast<RealType>( this->m_FullDomainSize[d] )
                                    / this->m_CurrentDomainSpacing[d] + 0.5 );
      this->m_CurrentDomainOrigin[d] =
        static_cast<unsigned long>( this->m_FullDomainSpacing[d]
                                    * static_cast<RealType>( this->m_FullDomainOrigin[d] )
                                    / this->m_CurrentDomainSpacing[d] + 0.5 );
      }
//            this->m_Debug=true;
    if( this->m_Debug )
      {
      std::cout << " outsize " << this->m_CurrentDomainSize <<  " curspc " << this->m_CurrentDomainSpacing
                << " fullspc " << this->m_FullDomainSpacing << " fullsz " <<  this->m_FullDomainSize   << std::endl;
      }
//            this->m_Debug=false;

    if( !this->m_DeformationField )
      {      /*FIXME -- need initial deformation strategy */
      this->m_DeformationField = DeformationFieldType::New();
      this->m_DeformationField->SetSpacing( this->m_CurrentDomainSpacing);
      this->m_DeformationField->SetOrigin( fixedImage->GetOrigin() );
      this->m_DeformationField->SetDirection( fixedImage->GetDirection() );
      typename ImageType::RegionType region;
      region.SetSize( this->m_CurrentDomainSize);
      this->m_DeformationField->SetLargestPossibleRegion(region);
      this->m_DeformationField->SetRequestedRegion(region);
      this->m_DeformationField->SetBufferedRegion(region);
      this->m_DeformationField->Allocate();
      this->m_DeformationField->FillBuffer(zero);
      std::cout <<  " allocated def field " << this->m_DeformationField->GetDirection() << std::endl;
      // exit(0);
      }
    else
      {
      this->m_DeformationField = this->ExpandField(this->m_DeformationField, this->m_CurrentDomainSpacing);
      if( this->m_TimeVaryingVelocity )
        {
        this->ExpandVelocity();
        }
      }
  }

  ImagePointer NormalizeImage( ImagePointer image)
  {
    typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput( image );
    minMaxFilter->Update();

    double min = minMaxFilter->GetMinimum();
    double shift = -1.0 * static_cast<double>( min );
    double scale = static_cast<double>( minMaxFilter->GetMaximum() );
    scale += shift;
    scale = 1.0 / scale;

    typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput( image );
    filter->SetShift( shift );
    filter->SetScale( scale );
    filter->Update();

    return filter->GetOutput();
  }

  void DeformableOptimization()
  {
    DeformationFieldPointer updateField = NULL;

    this->SetUpParameters();
    typename ImageType::SpacingType spacing;
    VectorType zero;
    zero.Fill(0);
    std::cout << " setting N-TimeSteps = "
              << this->m_NTimeSteps << " trunc " << this->m_GaussianTruncation << std::endl;

    unsigned int maxits = 0;
    for( unsigned int currentLevel = 0; currentLevel < this->m_NumberOfLevels; currentLevel++ )
      {
      if( this->m_Iterations[currentLevel] > maxits )
        {
        maxits = this->m_Iterations[currentLevel];
        }
      }
    if( maxits == 0 )
      {
      this->m_DeformationField = NULL;
      this->m_InverseDeformationField = NULL;
      //	this->ComputeMultiResolutionParameters(this->m_SimilarityMetrics[0]->GetFixedImage());
      return;
      }

    /* this is a hack to force univariate mappings   in the future,
       we will re-cast this framework  s.t. multivariate images can be used */
    unsigned int numberOfMetrics = this->m_SimilarityMetrics.size();
    for( unsigned int metricCount = 1;  metricCount < numberOfMetrics;  metricCount++ )
      {
      this->m_SimilarityMetrics[metricCount]->GetFixedImage()->SetOrigin(
        this->m_SimilarityMetrics[0]->GetFixedImage()->GetOrigin() );
      this->m_SimilarityMetrics[metricCount]->GetFixedImage()->SetDirection(
        this->m_SimilarityMetrics[0]->GetFixedImage()->GetDirection() );
      this->m_SimilarityMetrics[metricCount]->GetMovingImage()->SetOrigin( this->m_SimilarityMetrics[0]->GetMovingImage(
                                                                             )->GetOrigin() );
      this->m_SimilarityMetrics[metricCount]->GetMovingImage()->SetDirection(
        this->m_SimilarityMetrics[0]->GetMovingImage()->GetDirection() );
      }
    /* here, we assign all point set pointers to any single
       non-null point-set pointer */
    for( unsigned int metricCount = 0;  metricCount <   numberOfMetrics;  metricCount++ )
      {
      for( unsigned int metricCount2 = 0;  metricCount2 <   numberOfMetrics;  metricCount2++ )
        {
        if( this->m_SimilarityMetrics[metricCount]->GetFixedPointSet() )
          {
          this->m_SimilarityMetrics[metricCount2]->SetFixedPointSet(
            this->m_SimilarityMetrics[metricCount]->GetFixedPointSet() );
          }
        if( this->m_SimilarityMetrics[metricCount]->GetMovingPointSet() )
          {
          this->m_SimilarityMetrics[metricCount2]->SetMovingPointSet(
            this->m_SimilarityMetrics[metricCount]->GetMovingPointSet() );
          }
        }
      }
    this->m_SmoothFixedImages.resize(numberOfMetrics, NULL);
    this->m_SmoothMovingImages.resize(numberOfMetrics, NULL);
    for( unsigned int currentLevel = 0; currentLevel < this->m_NumberOfLevels; currentLevel++ )
      {
      this->m_CurrentLevel = currentLevel;
      typedef Vector<float, 1>                      ProfilePointDataType;
      typedef Image<ProfilePointDataType, 1>        CurveType;
      typedef PointSet<ProfilePointDataType, 1>     EnergyProfileType;
      typedef typename EnergyProfileType::PointType ProfilePointType;

      std::vector<EnergyProfileType::Pointer> energyProfiles;
      energyProfiles.resize( numberOfMetrics );
      for( unsigned int qq = 0; qq < numberOfMetrics; qq++ )
        {
        energyProfiles[qq] = EnergyProfileType::New();
        energyProfiles[qq]->Initialize();
        }

      ImagePointer fixedImage;
      ImagePointer movingImage;
      this->m_GradstepAltered = this->m_Gradstep;
      this->m_ScaleFactor = pow( 2.0, static_cast<RealType>( this->m_NumberOfLevels - currentLevel - 1 ) );
      std::cout << " this->m_ScaleFactor " << this->m_ScaleFactor
                << " nlev " << this->m_NumberOfLevels << " curl " << currentLevel << std::endl;
      /** FIXME -- here we assume the metrics all have the same image */
      fixedImage = this->m_SimilarityMetrics[0]->GetFixedImage();
      movingImage = this->m_SimilarityMetrics[0]->GetMovingImage();
      spacing = fixedImage->GetSpacing();
      this->ComputeMultiResolutionParameters(fixedImage);
      std::cout << " Its at this level " << this->m_Iterations[currentLevel] << std::endl;
      /*  generate smoothed images for all metrics */
      for( unsigned int metricCount = 0;  metricCount < numberOfMetrics;  metricCount++ )
        {
        this->m_SmoothFixedImages[metricCount] = this->SmoothImageToScale(
            this->m_SimilarityMetrics[metricCount]->GetFixedImage(), this->m_ScaleFactor);
        this->m_SmoothMovingImages[metricCount] = this->SmoothImageToScale(
            this->m_SimilarityMetrics[metricCount]->GetMovingImage(), this->m_ScaleFactor);
        }
      fixedImage = this->m_SmoothFixedImages[0];
      movingImage = this->m_SmoothMovingImages[0];

      unsigned int nmet = this->m_SimilarityMetrics.size();
      this->m_LastEnergy.resize(nmet, 1.e12);
      this->m_Energy.resize(nmet, 1.e9);
      this->m_EnergyBad.resize(nmet, 0);
      bool converged = false;
      this->m_CurrentIteration = 0;

      if( this->GetTransformationModel() != std::string("SyN") )
        {
        this->m_FixedImageAffineTransform = NULL;
        }
      while( !converged )
        {
        for( unsigned int metricCount = 0;  metricCount <   numberOfMetrics;  metricCount++ )
          {
          this->m_SimilarityMetrics[metricCount]->GetMetric()->SetIterations(this->m_CurrentIteration);
          }

        if( this->GetTransformationModel() == std::string("Elast") )
          {
          if( this->m_Iterations[currentLevel] > 0 )
            {
            this->ElasticRegistrationUpdate(fixedImage, movingImage);
            }
          }
        else if( this->GetTransformationModel() == std::string("SyN") )
          {
          if( currentLevel > 0  )
            {
            this->m_SyNF = this->ExpandField(this->m_SyNF, this->m_CurrentDomainSpacing);
            this->m_SyNFInv = this->ExpandField(this->m_SyNFInv, this->m_CurrentDomainSpacing);
            this->m_SyNM = this->ExpandField(this->m_SyNM, this->m_CurrentDomainSpacing);
            this->m_SyNMInv = this->ExpandField(this->m_SyNMInv, this->m_CurrentDomainSpacing);
            }
          if( this->m_Iterations[currentLevel] > 0 )
            {
            if( this->m_SyNType && this->m_ComputeThickness )
              {
              this->DiReCTUpdate(fixedImage, movingImage,
                                 this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                 this->m_SimilarityMetrics[0]->GetMovingPointSet() );
              }
            else if( this->m_SyNType )
              {
              this->SyNTVRegistrationUpdate(fixedImage, movingImage,
                                            this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                            this->m_SimilarityMetrics[0]->GetMovingPointSet() );
              }
            else
              {
              this->SyNRegistrationUpdate(fixedImage, movingImage,
                                          this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                          this->m_SimilarityMetrics[0]->GetMovingPointSet() );
              }
            }
          else if( this->m_SyNType )
            {
            this->UpdateTimeVaryingVelocityFieldWithSyNFandSyNM();
            }
          //            this->CopyOrAddToVelocityField( this->m_SyNF,  0 , false);
          }
        else if( this->GetTransformationModel() == std::string("Exp") )
          {
          if( this->m_Iterations[currentLevel] > 0 )
            {
            this->DiffeomorphicExpRegistrationUpdate(fixedImage, movingImage,
                                                     this->m_SimilarityMetrics[0]->GetFixedPointSet(),
                                                     this->m_SimilarityMetrics[0]->GetMovingPointSet() );
            }
          }

        this->m_CurrentIteration++;
        /**
         * This is where we track the energy profile to check for convergence.
         */
        for( unsigned int qq = 0; qq < numberOfMetrics; qq++ )
          {
          ProfilePointType point;
          point[0] = this->m_CurrentIteration - 1;

          ProfilePointDataType energy;
          energy[0] = this->m_Energy[qq];

          energyProfiles[qq]->SetPoint( this->m_CurrentIteration - 1, point );
          energyProfiles[qq]->SetPointData( this->m_CurrentIteration - 1, energy );
          }

        /**
         * If there are a sufficent number of iterations, fit a quadratic
         * single B-spline span to the number of energy profile points
         * in the first metric. To test convergence, evaluate the derivative
         * at the end of the profile to determine if >= 0.  To change to a
         * window of the energy profile, simply change the origin (assuming that
         * the desired window will start at the user-specified origin and
         * end at the current iteration).
         */
        if( this->m_CurrentIteration >= 3 )
          {
          typedef BSplineScatteredDataPointSetToImageFilter
          <EnergyProfileType, CurveType> BSplinerType;
          typename BSplinerType::Pointer bspliner
            = BSplinerType::New();

          typename CurveType::PointType origin;
          unsigned int domainorigin = 0;
          unsigned int domainsize = this->m_CurrentIteration - domainorigin;
          unsigned int domtar = 10;
          if( this->m_CurrentIteration > domtar )
            {
            domainsize = domtar;  domainorigin = this->m_CurrentIteration - domainsize;
            }
          origin.Fill( domainorigin );
          typename CurveType::SizeType size;
          size.Fill( domainsize );
          typename CurveType::SpacingType spacing;
          spacing.Fill( 1 );

          typename EnergyProfileType::Pointer energyProfileWindow = EnergyProfileType::New();
          energyProfileWindow->Initialize();

          unsigned int windowBegin = static_cast<unsigned int>( origin[0] );
          float        totale = 0;
          for( unsigned int qq = windowBegin; qq < this->m_CurrentIteration; qq++ )
            {
            ProfilePointType point;
            point[0] = qq;
            ProfilePointDataType energy;
            energy.Fill( 0 );
            energyProfiles[0]->GetPointData( qq, &energy );
            totale += energy[0];
            energyProfileWindow->SetPoint( qq - windowBegin, point );
            energyProfileWindow->SetPointData( qq - windowBegin, energy );
            }
//	  std::cout <<" totale " << totale << std::endl;
          if( totale > 0 )
            {
            totale *= (-1.0);
            }
          for( unsigned int qq = windowBegin; qq < this->m_CurrentIteration; qq++ )
            {
            ProfilePointDataType energy;  energy.Fill(0);
            energyProfiles[0]->GetPointData( qq, &energy );
            energyProfileWindow->SetPointData( qq - windowBegin, energy / totale);
            }

          bspliner->SetInput( energyProfileWindow );
          bspliner->SetOrigin( origin );
          bspliner->SetSpacing( spacing );
          bspliner->SetSize( size );
          bspliner->SetNumberOfLevels( 1 );
          unsigned int order = 1;
          bspliner->SetSplineOrder( order );
          typename BSplinerType::ArrayType ncps;
          ncps.Fill( order + 1);  // single span, order = 2
          bspliner->SetNumberOfControlPoints( ncps );
          bspliner->Update();

          ProfilePointType endPoint;
          //  endPoint[0] = static_cast<float>( this->m_CurrentIteration-1 );
          endPoint[0] = 0.9 * static_cast<float>( this->m_CurrentIteration - 1 );
          typename BSplinerType::GradientType gradient;
          bspliner->EvaluateGradientAtPoint( endPoint, gradient );
          if(  gradient[0][0]  < 0.0001 && this->m_CurrentIteration > 9 )
            {
            converged = true;
            }
          std::cout << " E-Slope " << gradient[0][0]; // << std::endl;
          }
        for( unsigned int qq = 0; qq < this->m_Energy.size(); qq++ )
          {
          if( qq == 0 )
            {
            std::cout << " iteration " << this->m_CurrentIteration;
            }

          std::cout << " energy " << qq << " : " << this->m_Energy[qq]; //  << " Last " << this->m_LastEnergy[qq];
          if( this->m_LastEnergy[qq] < this->m_Energy[qq] )
            {
            this->m_EnergyBad[qq]++;
            }
          }
        unsigned int numbade = 0;
        for( unsigned int qq = 0; qq < this->m_Energy.size(); qq++ )
          {
          if( this->m_CurrentIteration <= 1 )
            {
            this->m_EnergyBad[qq] = 0;
            }
          else if( this->m_EnergyBad[qq] > 1 )
            {
            numbade += this->m_EnergyBad[qq];
            }
          }

        // if ( this->m_EnergyBad[0] > 2)
        //         {
        //          this->m_GradstepAltered*=0.8;
        //	    std::cout <<" reducing gradstep " <<  this->m_GradstepAltered;
        // this->m_EnergyBad[this->m_Energy.size()-1]=0;
        // }
        std::cout << std::endl;

        if( this->m_CurrentIteration >= this->m_Iterations[currentLevel] )
          {
          converged = true;
          }
        //        || this->m_EnergyBad[0] >= 6 )
        //
        if( converged && this->m_CurrentIteration >= this->m_Iterations[currentLevel] )
          {
          std::cout << " tired convergence: reached max iterations " << std::endl;
          }
        else if( converged )
          {
          std::cout << " Converged due to oscillation in optimization ";
          for( unsigned int qq = 0; qq < this->m_Energy.size(); qq++ )
            {
            std::cout << " metric " << qq << " bad " << this->m_EnergyBad[qq] << "  ";
            }
          std::cout << std::endl;
          }
        }
      }

    if( this->GetTransformationModel() == std::string("SyN") )
      {
      //         float timestep=1.0/(float)this->m_NTimeSteps;
      //         unsigned int nts=this->m_NTimeSteps;
      if( this->m_SyNType )
        {
        //         this->m_SyNFInv = this->IntegrateConstantVelocity(this->m_SyNF, nts, timestep*(-1.));
        //         this->m_SyNMInv = this->IntegrateConstantVelocity(this->m_SyNM, nts, timestep*(-1.));
        //         this->m_SyNF= this->IntegrateConstantVelocity(this->m_SyNF, nts, timestep);
        //         this->m_SyNM= this->IntegrateConstantVelocity(this->m_SyNM,
        //         nts, timestep);
        //	 DeformationFieldPointer fdiffmap = this->IntegrateVelocity(0,0.5);
        //	 this->m_SyNFInv = this->IntegrateVelocity(0.5,0);
        //	 DeformationFieldPointer mdiffmap = this->IntegrateVelocity(0.5,1);
        //	 this->m_SyNMInv = this->IntegrateVelocity(1,0.5);
        //	 this->m_SyNM=this->CopyDeformationField(mdiffmap);
        //	 this->m_SyNF=this->CopyDeformationField(fdiffmap);
        this->m_DeformationField = this->IntegrateVelocity(0, 1);
        //	 ImagePointer wmimage= this->WarpMultiTransform(
        //  this->m_SmoothFixedImages[0],this->m_SmoothMovingImages[0], this->m_AffineTransform,
        // this->m_DeformationField, false , this->m_ScaleFactor );
        this->m_InverseDeformationField = this->IntegrateVelocity(1, 0);
        }
      else
        {
        this->m_InverseDeformationField = this->CopyDeformationField( this->m_SyNM);
        this->ComposeDiffs(this->m_SyNF, this->m_SyNMInv, this->m_DeformationField, 1);
        this->ComposeDiffs(this->m_SyNM, this->m_SyNFInv, this->m_InverseDeformationField, 1);
        }
      }
    else if( this->GetTransformationModel() == std::string("Exp") )
      {
      DeformationFieldPointer diffmap =
        this->IntegrateConstantVelocity( this->m_DeformationField, (unsigned int)this->m_NTimeSteps,
                                         1.0);                                                                                             //
                                                                                                                                           // /(float)this->m_NTimeSteps
      DeformationFieldPointer invdiffmap =
        this->IntegrateConstantVelocity(this->m_DeformationField, (unsigned int) this->m_NTimeSteps,
                                        -1.0);                                                                                               //
                                                                                                                                             // /(float)this->m_NTimeSteps*(-1.0)
      this->m_InverseDeformationField = invdiffmap;
      this->m_DeformationField = diffmap;
      AffineTransformPointer invaff = NULL;
      if( this->m_AffineTransform )
        {
        invaff = AffineTransformType::New();
        this->m_AffineTransform->GetInverse(invaff);
        if( this->m_Debug )
          {
          std::cout << " ??????invaff " << this->m_AffineTransform << std::endl << std::endl;
          }
        if( this->m_Debug )
          {
          std::cout << " invaff?????? " << invaff << std::endl << std::endl;
          }
        }
      }

    this->m_DeformationField->SetOrigin( this->m_SimilarityMetrics[0]->GetFixedImage()->GetOrigin() );
    this->m_DeformationField->SetDirection( this->m_SimilarityMetrics[0]->GetFixedImage()->GetDirection() );
    if( this->m_InverseDeformationField )
      {
      this->m_InverseDeformationField->SetOrigin( this->m_SimilarityMetrics[0]->GetFixedImage()->GetOrigin() );
      this->m_InverseDeformationField->SetDirection( this->m_SimilarityMetrics[0]->GetFixedImage()->GetDirection() );
      }

    if( this->m_TimeVaryingVelocity  )
      {
      std::string outname = localANTSGetFilePrefix(this->m_OutputNamingConvention.c_str() ) + std::string(
          "velocity.mhd");
      typename itk::ImageFileWriter<TimeVaryingVelocityFieldType>::Pointer writer =
        itk::ImageFileWriter<TimeVaryingVelocityFieldType>::New();
      writer->SetFileName(outname.c_str() );
      writer->SetInput( this->m_TimeVaryingVelocity);
      writer->UpdateLargestPossibleRegion();
      //	writer->Write();
      std::cout << " write tv field " << outname << std::endl;
      //        WriteImage<TimeVaryingVelocityFieldType>( this->m_TimeVaryingVelocity , outname.c_str());
      }
  }

  void DiffeomorphicExpRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints =
                                            NULL,
                                          PointSetPointer mpoints = NULL);

  void SyNRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints = NULL,
                             PointSetPointer mpoints = NULL);

  void SyNExpRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints = NULL,
                                PointSetPointer mpoints = NULL);

  void SyNTVRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints = NULL,
                               PointSetPointer mpoints = NULL);

  void DiReCTUpdate(ImagePointer fixedImage, ImagePointer movingImage, PointSetPointer fpoints = NULL,
                    PointSetPointer mpoints = NULL);

  /** allows one to copy or add a field to a time index within the velocity
* field
*/
  void UpdateTimeVaryingVelocityFieldWithSyNFandSyNM();

  void CopyOrAddToVelocityField( TimeVaryingVelocityFieldPointer velocity,  DeformationFieldPointer update1,
                                 DeformationFieldPointer update2,
                                 float timept);

// void CopyOrAddToVelocityField( DeformationFieldPointer update,  unsigned int timeindex,  bool
// CopyIsTrueOtherwiseAdd);

  void ElasticRegistrationUpdate(ImagePointer fixedImage, ImagePointer movingImage)
  {
    typename ImageType::SpacingType spacing;
    VectorType zero;
    zero.Fill(0);
    DeformationFieldPointer updateField;

    updateField = this->ComputeUpdateField(this->m_DeformationField, NULL, NULL, NULL, NULL);

    typedef ImageRegionIteratorWithIndex<DeformationFieldType> Iterator;
    Iterator dIter(this->m_DeformationField, this->m_DeformationField->GetLargestPossibleRegion() );
    for( dIter.GoToBegin(); !dIter.IsAtEnd(); ++dIter )
      {
      typename ImageType::IndexType index = dIter.GetIndex();
      VectorType vec = updateField->GetPixel(index);
      dIter.Set(dIter.Get() + vec * this->m_Gradstep);
      }

    if( this->m_Debug )
      {
      std::cout << " updated elast " << " up-sz " << updateField->GetLargestPossibleRegion() <<  std::endl;
      std::cout <<  " t-sz " << this->m_DeformationField->GetLargestPossibleRegion() <<  std::endl;
      }
    this->SmoothDeformationField(this->m_DeformationField, false);

    return;
  }

  ImagePointer WarpImageBackward( ImagePointer image, DeformationFieldPointer field )
  {
    typedef WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;
    typename WarperType::Pointer  warper = WarperType::New();
    typedef NearestNeighborInterpolateImageFunction<ImageType, double>
    InterpolatorType;
    warper->SetInput(image);
    warper->SetDeformationField( field );
    warper->SetEdgePaddingValue( 0);
    warper->SetOutputSpacing(field->GetSpacing() );
    warper->SetOutputOrigin( field->GetOrigin() );
    warper->Update();
    return warper->GetOutput();
  }

  void ComposeDiffs(DeformationFieldPointer fieldtowarpby, DeformationFieldPointer field,
                    DeformationFieldPointer fieldout,
                    float sign);

  void SetSimilarityMetrics( SimilarityMetricListType S )
  {
    this->m_SimilarityMetrics = S;
  }

  void SetFixedPointSet(  PointSetPointer p )
  {
    this->m_FixedPointSet = p;
  }

  void SetMovingPointSet(  PointSetPointer p )
  {
    this->m_MovingPointSet = p;
  }

  void SetDeltaTime( float t)
  {
    this->m_DeltaTime = t;
  }

  float InvertField(DeformationFieldPointer field,
                    DeformationFieldPointer inverseField, float weight = 1.0,
                    float toler = 0.5, int maxiter = 20, bool print = false)
  {
    VectorType zero; zero.Fill(0);

    //  if (this->GetElapsedIterations() < 2 ) maxiter=10;

    ImagePointer floatImage = ImageType::New();
    floatImage->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
    floatImage->SetBufferedRegion( field->GetLargestPossibleRegion().GetSize() );
    floatImage->SetSpacing(field->GetSpacing() );
    floatImage->SetOrigin(field->GetOrigin() );
    floatImage->SetDirection(field->GetDirection() );
    floatImage->Allocate();

    typedef typename DeformationFieldType::PixelType           VectorType;
    typedef typename DeformationFieldType::IndexType           IndexType;
    typedef typename VectorType::ValueType                     ScalarType;
    typedef ImageRegionIteratorWithIndex<DeformationFieldType> Iterator;

    DeformationFieldPointer lagrangianInitCond = DeformationFieldType::New();
    lagrangianInitCond->SetSpacing( field->GetSpacing() );
    lagrangianInitCond->SetOrigin( field->GetOrigin() );
    lagrangianInitCond->SetDirection( field->GetDirection() );
    lagrangianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
    lagrangianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
    lagrangianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
    lagrangianInitCond->Allocate();
    DeformationFieldPointer eulerianInitCond = DeformationFieldType::New();
    eulerianInitCond->SetSpacing( field->GetSpacing() );
    eulerianInitCond->SetOrigin( field->GetOrigin() );
    eulerianInitCond->SetDirection( field->GetDirection() );
    eulerianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
    eulerianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
    eulerianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
    eulerianInitCond->Allocate();

    typedef typename DeformationFieldType::SizeType SizeType;
    SizeType size = field->GetLargestPossibleRegion().GetSize();

    typename ImageType::SpacingType spacing = field->GetSpacing();
    float         subpix = 0.0;
    unsigned long npix = 1;
    for( int j = 0; j < ImageDimension; j++ ) // only use in-plane spacing
      {
      npix *= field->GetLargestPossibleRegion().GetSize()[j];
      }
    subpix = pow( (float)ImageDimension, (float)ImageDimension) * 0.5;

    float    max = 0;
    Iterator iter( field, field->GetLargestPossibleRegion() );
    for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      IndexType  index = iter.GetIndex();
      VectorType vec1 = iter.Get();
      VectorType newvec = vec1 * weight;
      lagrangianInitCond->SetPixel(index, newvec);
      float mag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        mag += newvec[jj] * newvec[jj];
        }
      mag = sqrt(mag);
      if( mag > max )
        {
        max = mag;
        }
      }

    eulerianInitCond->FillBuffer(zero);

    float scale = (1.) / max;
    if( scale > 1. )
      {
      scale = 1.0;
      }
//    float initscale=scale;
    Iterator vfIter( inverseField, inverseField->GetLargestPossibleRegion() );

//  int num=10;
//  for (int its=0; its<num; its++)
    float difmag = 10.0;
    int   ct = 0;
    float denergy = 10;
    float denergy2 = 10;
    float laste = 1.e9;
    float meandif = 1.e8;
//    int badct=0;
//  while (difmag > subpix && meandif > subpix*0.1 && badct < 2 )//&& ct < 20 && denergy > 0)
//    float length=0.0;
    float stepl = 2.;
    float lastdifmag = 0;

    float epsilon = (float)size[0] / 256;
    if( epsilon > 1 )
      {
      epsilon = 1;
      }

    while( difmag > toler && ct<maxiter && meandif> 0.001 )
      {
      denergy = laste - difmag; // meandif;
      denergy2 = laste - meandif;
      laste = difmag; // meandif;
      meandif = 0.0;

      // this field says what position the eulerian field should contain in the E domain
      this->ComposeDiffs(inverseField, lagrangianInitCond,    eulerianInitCond, 1);
      difmag = 0.0;
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        IndexType  index = vfIter.GetIndex();
        VectorType update = eulerianInitCond->GetPixel(index);
        float      mag = 0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          update[j] *= (-1.0);
          mag += (update[j] / spacing[j]) * (update[j] / spacing[j]);
          }
        mag = sqrt(mag);
        meandif += mag;
        if( mag > difmag )
          {
          difmag = mag;
          }
        //	  if (mag < 1.e-2) update.Fill(0);

        eulerianInitCond->SetPixel(index, update);
        floatImage->SetPixel(index, mag);
        }
      meandif /= (float)npix;
      if( ct == 0 )
        {
        epsilon = 0.75;
        }
      else
        {
        epsilon = 0.5;
        }
      stepl = difmag * epsilon;
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        float      val = floatImage->GetPixel(vfIter.GetIndex() );
        VectorType update = eulerianInitCond->GetPixel(vfIter.GetIndex() );
        if( val > stepl )
          {
          update = update * (stepl / val);
          }
        VectorType upd = vfIter.Get() + update * (epsilon);
        vfIter.Set(upd);
        }
      ct++;
      lastdifmag = difmag;
      }

    // std::cout <<" difmag " << difmag << ": its " << ct <<  std::endl;

    return difmag;
  }

  void SetUseNearestNeighborInterpolation( bool useNN)
  {
    this->m_UseNN = useNN;
  }

protected:

  DeformationFieldPointer IntegrateVelocity(float, float);

  DeformationFieldPointer IntegrateLandmarkSetVelocity(float, float, PointSetPointer movingpoints,
                                                       ImagePointer referenceimage );

  VectorType IntegratePointVelocity(float starttimein, float finishtimein, IndexType startPoint);

  ImagePointer  MakeSubImage( ImagePointer bigimage)
  {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    ImagePointer varimage = ImageType::New();

    typename ImageType::RegionType region;
    typename ImageType::SizeType size = bigimage->GetLargestPossibleRegion().GetSize();
    region.SetSize( this->m_CurrentDomainSize);
    typename ImageType::IndexType index;  index.Fill(0);
    region.SetIndex(index);
    varimage->SetRegions( region );
    varimage->SetSpacing(this->m_CurrentDomainSpacing);
    varimage->SetOrigin(bigimage->GetOrigin() );
    varimage->SetDirection(bigimage->GetDirection() );
    varimage->Allocate();
    varimage->FillBuffer(0);

    typename ImageType::IndexType cornerind;
    cornerind.Fill(0);
    for( unsigned int ii = 0; ii < ImageDimension; ii++ )
      {
      float diff = (float)this->m_CurrentDomainOrigin[ii] - (float)this->m_CurrentDomainSize[ii] / 2;
      if( diff < 0 )
        {
        diff = 0;
        }
      cornerind[ii] = (unsigned long) diff;
      }
    //  std::cout << " corner index " << cornerind << std::endl;
    Iterator vfIter2( bigimage,  bigimage->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      typename ImageType::IndexType origindex = vfIter2.GetIndex();
      typename ImageType::IndexType index = vfIter2.GetIndex();
      bool oktosample = true;
      for( unsigned int ii = 0; ii < ImageDimension; ii++ )
        {
        float centerbasedcoord        = (origindex[ii] - this->m_CurrentDomainOrigin[ii]);
//                float diff =
//                    index[ii]=origindex[ii]-cornerind[ii];
        if( fabs(centerbasedcoord) >  (this->m_CurrentDomainSize[ii] / 2 - 1) )
          {
          oktosample = false;
          }
        }
      if( oktosample )
        {
        //      std::cout << " index " << index <<  " origindex " << origindex << " ok? " << oktosample << std::endl;
        varimage->SetPixel(index, bigimage->GetPixel(origindex) );
        }
      }
    // std::cout << " sizes " << varimage->GetLargestPossibleRegion().GetSize() << " bigimage " <<
    // bigimage->GetLargestPossibleRegion().GetSize() << std::endl;
    return varimage;
  }

  float MeasureDeformation(DeformationFieldPointer field, int option = 0)
  {
    typedef typename DeformationFieldType::PixelType           VectorType;
    typedef typename DeformationFieldType::IndexType           IndexType;
    typedef typename DeformationFieldType::SizeType            SizeType;
    typedef typename VectorType::ValueType                     ScalarType;
    typedef ImageRegionIteratorWithIndex<DeformationFieldType> Iterator;
    // all we have to do here is add the local field to the global field.
    Iterator      vfIter( field,  field->GetLargestPossibleRegion() );
    SizeType      size = field->GetLargestPossibleRegion().GetSize();
    unsigned long ct = 1;
    double        totalmag = 0;
    float         maxstep = 0;
//  this->m_EuclideanNorm=0;

    typename ImageType::SpacingType myspacing = field->GetSpacing();
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      IndexType  index = vfIter.GetIndex();
      IndexType  rindex = vfIter.GetIndex();
      IndexType  lindex = vfIter.GetIndex();
      VectorType update = vfIter.Get();
      float      mag = 0.0;
      float      stepl = 0.0;
      for( int i = 0; i < ImageDimension; i++ )
        {
        rindex = index;
        lindex = index;
        if( (int)index[i] < (int)(size[i] - 2) )
          {
          rindex[i] = rindex[i] + 1;
          }
        if( index[i] > 2 )
          {
          lindex[i] = lindex[i] - 1;
          }
        VectorType rupdate = field->GetPixel(rindex);
        VectorType lupdate = field->GetPixel(lindex);
        VectorType dif = rupdate - lupdate;
        for( int tt = 0; tt < ImageDimension; tt++ )
          {
          stepl += update[tt] * update[tt] / (myspacing[tt] * myspacing[tt]);
          mag += dif[tt] * dif[tt] / (myspacing[tt] * myspacing[tt]);
          }
        }
      stepl = sqrt(stepl);
      mag = sqrt(mag);
      if( stepl > maxstep )
        {
        maxstep = stepl;
        }
      ct++;
      totalmag += mag;
//    this->m_EuclideanNorm+=stepl;
      }
    // this->m_EuclideanNorm/=ct;
    // this->m_ElasticPathLength = totalmag/ct;
    // this->m_LinftyNorm = maxstep;
//  std::cout << " Elast path length " << this->m_ElasticPathLength <<  " L inf norm " << this->m_LinftyNorm <<
// std::endl;
// if (this->m_ElasticPathLength >= this->m_ArcLengthGoal)
//  if (maxstep >= this->m_ArcLengthGoal)
      {
//  this->StopRegistration();
// scale the field to the right length
//  float scale=this->m_ArcLengthGoal/this->m_ElasticPathLength;
//  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )vfIter.Set(vfIter.Get()*scale);
      }
    // if (this->m_LinftyNorm <= 0) this->m_LinftyNorm=1;
    // if (this->m_ElasticPathLength <= 0) this->m_ElasticPathLength=0;
    // if (this->m_EuclideanNorm <= 0) this->m_EuclideanNorm=0;

    // if (option==0) return this->m_ElasticPathLength;
    // else if (option==2) return this->m_EuclideanNorm;
    // else
    return maxstep;
  }

  ANTSImageRegistrationOptimizer();
  virtual ~ANTSImageRegistrationOptimizer()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  ANTSImageRegistrationOptimizer( const Self & ); // purposely not implemented
  void operator=( const Self & );                 // purposely not implemented

  typename VelocityFieldInterpolatorType::Pointer m_VelocityFieldInterpolator;

  typename ImageType::SizeType   m_CurrentDomainSize;
  typename ImageType::PointType   m_CurrentDomainOrigin;
  typename ImageType::SpacingType   m_CurrentDomainSpacing;
  typename ImageType::DirectionType   m_CurrentDomainDirection;
  typename ImageType::SizeType   m_FullDomainSize;
  typename ImageType::PointType   m_FullDomainOrigin;
  typename ImageType::SpacingType   m_FullDomainSpacing;

  AffineTransformPointer  m_AffineTransform;
  AffineTransformPointer  m_FixedImageAffineTransform;
  DeformationFieldPointer m_DeformationField;
  DeformationFieldPointer m_InverseDeformationField;
  DeformationFieldPointer m_StaticVelocityField;

  std::vector<float>        m_GradientDescentParameters;
  std::vector<float>        m_MetricScalarWeights;
  std::vector<ImagePointer> m_SmoothFixedImages;
  std::vector<ImagePointer> m_SmoothMovingImages;

  bool         m_Debug;
  unsigned int m_NumberOfLevels;
  typename ParserType::Pointer m_Parser;
  SimilarityMetricListType  m_SimilarityMetrics;
  ImagePointer              m_MaskImage;
  float                     m_ScaleFactor;
  bool                      m_UseMulti;
  bool                      m_UseROI;
  bool                      m_UseNN;
  unsigned int              m_CurrentIteration;
  unsigned int              m_CurrentLevel;
  std::string               m_TransformationModel;
  std::string               m_OutputNamingConvention;
  PointSetPointer           m_FixedPointSet;
  PointSetPointer           m_MovingPointSet;
  std::vector<unsigned int> m_Iterations;
  std::vector<float>        m_RoiNumbers;
  float                     m_GradSmoothingparam;
  float                     m_TotalSmoothingparam;
  float                     m_Gradstep;
  float                     m_GradstepAltered;
  float                     m_NTimeSteps;
  float                     m_GaussianTruncation;
  float                     m_DeltaTime;

/** energy stuff */
  std::vector<float>        m_Energy;
  std::vector<float>        m_LastEnergy;
  std::vector<unsigned int> m_EnergyBad;

/** for SyN only */
  DeformationFieldPointer         m_SyNF;
  DeformationFieldPointer         m_SyNFInv;
  DeformationFieldPointer         m_SyNM;
  DeformationFieldPointer         m_SyNMInv;
  TimeVaryingVelocityFieldPointer m_TimeVaryingVelocity;
  unsigned int                    m_SyNType;

/** for BSpline stuff */
  unsigned int m_BSplineFieldOrder;
  ArrayType    m_GradSmoothingMeshSize;
  ArrayType    m_TotalSmoothingMeshSize;

/** For thickness calculation */
  ImagePointer m_HitImage;
  ImagePointer m_ThickImage;
  unsigned int m_ComputeThickness;
  unsigned int m_SyNFullTime;
};
}
// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkANTSImageRegistrationOptimizer.cxx"
#endif

#endif
