/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCrossCorrelationRegistrationFunction.h,v $
  Language:  C++
  Date:      $Date: 2009/01/05 20:09:47 $
  Version:   $Revision: 1.19 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCrossCorrelationRegistrationFunction_h_
#define _itkCrossCorrelationRegistrationFunction_h_

#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "itkAvantsMutualInformationRegistrationFunction.h"

namespace itk
{
/**
 * \class CrossCorrelationRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by CrossCorrelationRegistrationFilter to compute the
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa CrossCorrelationRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT CrossCorrelationRegistrationFunction :
  public         AvantsPDEDeformableRegistrationFunction<TFixedImage,
                                                         TMovingImage, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef CrossCorrelationRegistrationFunction Self;
  typedef AvantsPDEDeformableRegistrationFunction<TFixedImage,
                                                  TMovingImage, TDeformationField>    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( CrossCorrelationRegistrationFunction,
                AvantsPDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::MetricImageType          MetricImageType;
  typedef typename Superclass::MetricImageType::Pointer MetricImagePointer;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;
  typedef typename FixedImageType::IndexType            IndexType;
  typedef typename FixedImageType::SizeType             SizeType;

  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer
    DeformationFieldTypePointer;
  typedef typename TDeformationField::PixelType VectorType;

  typedef CovariantVector<float,
                          itkGetStaticConstMacro(ImageDimension)> GradientPixelType;
  typedef Image<GradientPixelType,
                itkGetStaticConstMacro(ImageDimension)> GradientImageType;
  typedef SmartPointer<GradientImageType> GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter<MetricImageType, GradientImageType>
    GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer            GradientImageFilterPointer;
  typedef Image<float, itkGetStaticConstMacro(ImageDimension)> BinaryImageType;
  typedef typename BinaryImageType::Pointer                    BinaryImagePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
//  typedef typename Superclass::NeighborhoodType    BoundaryNeighborhoodType;
  typedef typename Superclass::FloatOffsetType FloatOffsetType;
  typedef typename Superclass::TimeStepType    TimeStepType;

  /** Interpolator type. */
  typedef double                                                  CoordRepType;
  typedef InterpolateImageFunction<MovingImageType, CoordRepType> InterpolatorType;
  typedef typename InterpolatorType::Pointer                      InterpolatorPointer;
  typedef typename InterpolatorType::PointType                    PointType;
  typedef LinearInterpolateImageFunction<MovingImageType, CoordRepType>
    DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double, itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
  {
    m_MovingImageInterpolator = ptr;
  }

  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator(void)
  {
    return m_MovingImageInterpolator;
  }

  typename TDeformationField::PixelType ComputeMetricAtPairB(IndexType fixedindex,
                                                             typename TDeformationField::PixelType vec );
  typename TDeformationField::PixelType ComputeMetricAtPairC(IndexType fixedindex,
                                                             typename TDeformationField::PixelType vec );

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void *GlobalData) const
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void * GetGlobalDataPointer() const
  {
    GlobalDataStruct *global = new GlobalDataStruct();

    return global;
  }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
  {
    delete (GlobalDataStruct *) GlobalData;
  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  virtual VectorType  OpticalFlowUpdate(const NeighborhoodType & neighborhood, bool returninv = false)
  {
    // Get fixed image related information
    IndexType index = neighborhood.GetIndex();

    typename TDeformationField::PixelType vec;
    vec.Fill(0);
    VectorType update;
    update.Fill(0.0);
    double              fixedValue;
    CovariantVectorType fixedGradient, movingGradient;
    double              fixedGradientSquaredMagnitude = 0;
    double              movingGradientSquaredMagnitude = 0;
    fixedValue = (double) this->finitediffimages[0]->GetPixel( index );
    double movingValue = (double) this->finitediffimages[1]->GetPixel( index );

    if( !returninv )
      {
      fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] );
        }

      double speedValue = fixedValue - movingValue;
      if( fabs(speedValue) < this->m_RobustnessParameter )
        {
        speedValue = 0;
        }
      double denominator = vnl_math_sqr( speedValue ) / m_Normalizer
        + fixedGradientSquaredMagnitude;
      double m_DenominatorThreshold = 1e-9;
      double m_IntensityDifferenceThreshold = 0.001;
      if( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold ||
          denominator < m_DenominatorThreshold )
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          update[j] = 0.0;
          }
        return update;
        }
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        update[j] = speedValue * fixedGradient[j] / denominator;
        }
      this->m_Energy += speedValue * speedValue;

      return update;
      }
    else
      {
      // Get fixed image related information
      movingGradient = m_MovingImageGradientCalculator->EvaluateAtIndex( index );
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        movingGradientSquaredMagnitude += vnl_math_sqr( movingGradient[j] );
        }

      double speedValue = movingValue - fixedValue;
      double denominator = vnl_math_sqr( speedValue ) / m_Normalizer
        + movingGradientSquaredMagnitude;
      double m_DenominatorThreshold = 1e-9;
      double m_IntensityDifferenceThreshold = 0.001;
      if( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold ||
          denominator < m_DenominatorThreshold )
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          update[j] = 0.0;
          }
        return update;
        }
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        update[j] = speedValue * movingGradient[j] / denominator;
        }
      return update;
      }
  }

  void SetFixedImageMask( MetricImageType* img)
  {
    m_FixedImageMask = img;
  }

  virtual VectorType ComputeUpdate(const NeighborhoodType & neighborhood,
                                   void *globalData,
                                   const FloatOffsetType & offset = FloatOffsetType(0.0) )
  {
    VectorType update;

    update.Fill(0.0);
    //  std::cout << " updating " << std::endl;
    IndexType oindex = neighborhood.GetIndex();
//    std::cout <<" Index " << oindex << std::endl;
//    if (!this->GetFixedImage() ) { std::cout << " no image " << std::endl; return update;}
    update = this->OpticalFlowUpdate(neighborhood, false);
    // update=this->ComputeMetricAtPairB(oindex,update);

    return update;
  }

  virtual VectorType ComputeUpdateInv(const NeighborhoodType & neighborhood,
                                      void *globalData,
                                      const FloatOffsetType & offset = FloatOffsetType(0.0) )
  {
    VectorType update;

    update.Fill(0.0);
    IndexType       oindex = neighborhood.GetIndex();
    FixedImageType* img = const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer() );
    if( !img )
      {
      return update;
      }

    update = this->OpticalFlowUpdate(neighborhood, true);
    // update=this->ComputeMetricAtPairC(oindex,update);

    return update;
  }

  double ComputeCrossCorrelation()
  {
    if( finitediffimages[0] )
      {
      double        totalcc = 0;
      unsigned long ct = 0;
      typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
      ittype it(this->finitediffimages[0], this->finitediffimages[0]->GetLargestPossibleRegion().GetSize() );
      for( it.GoToBegin(); !it.IsAtEnd(); ++it )
        {
        IndexType oindex = it.GetIndex();
        double    sfm = finitediffimages[2]->GetPixel(oindex);
        double    sff = finitediffimages[3]->GetPixel(oindex);
        double    smm = finitediffimages[4]->GetPixel(oindex);
        double    cc = 0;
        if( fabs(sff * smm) > 0 )
          {
          cc += sfm * sfm / (sff * smm); ct++;
          }
        totalcc += cc;
        }
      this->m_Energy = totalcc / (float)ct * (-1.0);
      return this->m_Energy * (-1.0);
      }
    else
      {
      return 0;
      }
  }

  void SetFullyRobust(bool b)
  {
    m_FullyRobust = b;
  }

  void GetProbabilities();

  double localCrossCorrelation;
  float  m_TEMP;
protected:
  CrossCorrelationRegistrationFunction();
  ~CrossCorrelationRegistrationFunction()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
    };

  MetricImagePointer MakeImage()
  {
    typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
    typedef ImageRegionIteratorWithIndex<BinaryImageType> ittype2;
    FixedImageType* img = const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer() );
    typename FixedImageType::SizeType imagesize = img->GetLargestPossibleRegion().GetSize();

      {
      Superclass::m_MetricImage = MetricImageType::New();
      Superclass::m_MetricImage->SetLargestPossibleRegion(img->GetLargestPossibleRegion()  );
      Superclass::m_MetricImage->SetBufferedRegion(img->GetLargestPossibleRegion() );
      Superclass::m_MetricImage->SetSpacing(img->GetSpacing() );
      Superclass::m_MetricImage->SetOrigin(img->GetOrigin() );
      Superclass::m_MetricImage->Allocate();
      //    ittype it(m_MetricImage,Superclass::m_MetricImage->GetLargestPossibleRegion().GetSize());
      //    for( it.GoToBegin(); !it.IsAtEnd(); ++it ) it.Set(0);
      }
    bool makebinimg = false;
    //    if (!binaryimage) makebinimg=true;
    // else if ( binaryimage->GetLargestPossibleRegion().GetSize()[0] !=
    //                img->GetLargestPossibleRegion().GetSize()[0] )makebinimg=true;

    if( makebinimg )
      {
      m_Iteration = 0;
      binaryimage = BinaryImageType::New();
      binaryimage->SetLargestPossibleRegion(img->GetLargestPossibleRegion()  );
      binaryimage->SetBufferedRegion(img->GetLargestPossibleRegion() );
      binaryimage->SetSpacing(img->GetSpacing() );
      binaryimage->SetOrigin(img->GetOrigin() );
      binaryimage->Allocate();
      ittype2 it(binaryimage, binaryimage->GetLargestPossibleRegion().GetSize() );
      for( it.GoToBegin(); !it.IsAtEnd(); ++it )
        {
        it.Set(1);
        }
      }
    return Superclass::m_MetricImage;
  }

private:
  CrossCorrelationRegistrationFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                       // purposely not implemented

  /** Cache fixed image information. */
  typename TFixedImage::SpacingType                  m_FixedImageSpacing;
  typename TFixedImage::PointType                  m_FixedImageOrigin;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer m_FixedImageGradientCalculator;

  GradientCalculatorPointer m_MovingImageGradientCalculator;

  /** Function to interpolate the moving image. */
  InterpolatorPointer m_MovingImageInterpolator;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_IntensityDifferenceThreshold;

  mutable double m_MetricTotal;
  mutable float  m_MinMag;
  mutable float  m_MaxMag;
  mutable float  m_AvgMag;
  mutable float  m_Thresh;

  GradientImagePointer m_MetricGradientImage;

  MetricImagePointer finitediffimages[5];
  BinaryImagePointer binaryimage;

  MetricImagePointer m_FixedImageMask;
  MetricImagePointer m_MovingImageMask;

  typedef itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType,
                                                           DeformationFieldType> MetricType2;
  typename MetricType2::Pointer m_MIFunct;
  unsigned int m_NumberOfHistogramBins;
  bool         m_FullyRobust;
  unsigned int m_Iteration;
  float        m_Normalizer;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCrossCorrelationRegistrationFunction.cxx"
#endif

#endif
