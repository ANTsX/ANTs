/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAvantsMutualInformationRegistrationFunction.h,v $
  Language:  C++
  Date:      $Date: 2009/01/08 15:14:48 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAvantsMutualInformationRegistrationFunction_h
#define __itkAvantsMutualInformationRegistrationFunction_h

#include "itkImageFileWriter.h"
#include "itkImageToImageMetric.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkBSplineKernelFunction.h"
#include "itkBSplineDerivativeKernelFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineDeformableTransform.h"
#include "itkTranslationTransform.h"
#include "itkArray2D.h"
#include "itkImageBase.h"
#include "itkTransform.h"
#include "itkInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkExceptionObject.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSpatialObject.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
/** \class AvantsMutualInformationRegistrationFunction
 * \brief Computes the mutual information between two images to be
 * registered using the method of Avants et al.
 *
 * AvantsMutualInformationRegistrationFunction computes the mutual
 * information between a fixed and moving image to be registered.
 *
 * This class is templated over the FixedImage type and the MovingImage
 * type.
 *
 * The fixed and moving images are set via methods SetFixedImage() and
 * SetMovingImage(). This metric makes use of user specified Transform and
 * Interpolator. The Transform is used to map points from the fixed image to
 * the moving image domain. The Interpolator is used to evaluate the image
 * intensity at user specified geometric points in the moving image.
 * The Transform and Interpolator are set via methods SetTransform() and
 * SetInterpolator().
 *
 * If a BSplineInterpolationFunction is used, this class obtain
 * image derivatives from the BSpline interpolator. Otherwise,
 * image derivatives are computed using central differencing.
 *
 * \warning This metric assumes that the moving image has already been
 * connected to the interpolator outside of this class.
 *
 * The method GetValue() computes of the mutual information
 * while method GetValueAndDerivative() computes
 * both the mutual information and its derivatives with respect to the
 * transform parameters.
 *
 * The calculations are based on the method of Avants et al [1,2]
 * where the probability density distribution are estimated using
 * Parzen histograms. Since the fixed image PDF does not contribute
 * to the derivatives, it does not need to be smooth. Hence,
 * a zero order (box car) BSpline kernel is used
 * for the fixed image intensity PDF. On the other hand, to ensure
 * smoothness a third order BSpline kernel is used for the
 * moving image intensity PDF.
 *
 * On Initialize(), the FixedImage is uniformly sampled within
 * the FixedImageRegion. The number of samples used can be set
 * via SetNumberOfSpatialSamples(). Typically, the number of
 * spatial samples used should increase with the image size.
 *
 * During each call of GetValue(), GetDerivatives(),
 * GetValueAndDerivatives(), marginal and joint intensity PDF's
 * values are estimated at discrete position or bins.
 * The number of bins used can be set via SetNumberOfHistogramBins().
 * To handle data with arbitray magnitude and dynamic range,
 * the image intensity is scale such that any contribution to the
 * histogram will fall into a valid bin.
 *
 * One the PDF's have been contructed, the mutual information
 * is obtained by doubling summing over the discrete PDF values.
 *
 *
 * Notes:
 * 1. This class returns the negative mutual information value.
 * 2. This class in not thread safe due the private data structures
 *     used to the store the sampled points and the marginal and joint pdfs.
 *
 * References:
 * [1] "Nonrigid multimodality image registration"
 *      D. Avants, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      Medical Imaging 2001: Image Processing, 2001, pp. 1609-1620.
 * [2] "PET-CT Image Registration in the Chest Using Free-form Deformations"
 *      D. Avants, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      IEEE Transactions in Medical Imaging. Vol.22, No.1,
        January 2003. pp.120-128.
 * [3] "Optimization of Mutual Information for MultiResolution Image
 *      Registration"
 *      P. Thevenaz and M. Unser
 *      IEEE Transactions in Image Processing, 9(12) December 2000.
 *
 * \ingroup RegistrationMetrics
 * \ingroup ThreadUnSafe
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT AvantsMutualInformationRegistrationFunction :
  public         AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef AvantsMutualInformationRegistrationFunction Self;
  typedef AvantsPDEDeformableRegistrationFunction<TFixedImage,
                                                  TMovingImage, TDeformationField>    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( AvantsMutualInformationRegistrationFunction,
                AvantsPDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;
  typedef typename FixedImageType::IndexType     IndexType;
  typedef typename FixedImageType::SizeType      SizeType;
  typedef typename FixedImageType::SpacingType   SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::VectorType           VectorType;
  typedef typename Superclass::DeformationFieldType DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer
    DeformationFieldTypePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Interpolator type. */
  typedef double CoordRepType;
  typedef    //       //    LinearInterpolateImageFunction<MovingImageType,CoordRepType>
    BSplineInterpolateImageFunction<MovingImageType, CoordRepType>
    InterpolatorType;
  typedef typename InterpolatorType::Pointer   InterpolatorPointer;
  typedef typename InterpolatorType::PointType PointType;
  typedef InterpolatorType                     DefaultInterpolatorType;
  //  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
  // DefaultInterpolatorType;

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

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void *itkNotUsed(GlobalData) ) const
  {
    return 1;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void * GetGlobalDataPointer() const
  {
    GlobalDataStruct *global = new GlobalDataStruct();

    //    global->m_SumOfSquaredDifference  = 0.0;
    /// global->m_NumberOfPixelsProcessed = 0L;
    // global->m_SumOfSquaredChange      = 0;
    return global;
  }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
  {
    delete (GlobalDataStruct *) GlobalData;
  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  typedef double CoordinateRepresentationType;

  /** Types inherited from Superclass. */
  typedef TranslationTransform<CoordinateRepresentationType,
                               //                    itkGetStaticConstMacro(ImageDimension),
                               itkGetStaticConstMacro(ImageDimension)> TransformType;

  typedef ImageToImageMetric<TFixedImage, TMovingImage> Metricclass;

  typedef typename TransformType::Pointer             TransformPointer;
  typedef typename Metricclass::TransformJacobianType TransformJacobianType;
  //  typedef typename Metricclass::InterpolatorType         InterpolatorType;
  typedef typename Metricclass::MeasureType             MeasureType;
  typedef typename Metricclass::DerivativeType          DerivativeType;
  typedef typename TransformType::ParametersType        ParametersType;
  typedef typename Metricclass::FixedImageConstPointer  FixedImageConstPointer;
  typedef typename Metricclass::MovingImageConstPointer MovingImageCosntPointer;
  // typedef typename TransformType::CoordinateRepresentationType  CoordinateRepresentationType;

  /** Index and Point typedef support. */
  typedef typename FixedImageType::IndexType           FixedImageIndexType;
  typedef typename FixedImageIndexType::IndexValueType FixedImageIndexValueType;
  typedef typename MovingImageType::IndexType          MovingImageIndexType;
  typedef typename TransformType::InputPointType       FixedImagePointType;
  typedef typename TransformType::OutputPointType      MovingImagePointType;

  /** The marginal PDFs are stored as std::vector. */
  typedef float PDFValueType;
  //  typedef std::vector<PDFValueType> MarginalPDFType;
  typedef Image<PDFValueType, 1>              MarginalPDFType;
  typedef typename MarginalPDFType::IndexType MarginalPDFIndexType;
  /** Typedef for the joint PDF and PDF derivatives are stored as ITK Images. */
  typedef Image<PDFValueType, 2>              JointPDFType;
  typedef Image<PDFValueType, 3>              JointPDFDerivativesType;
  typedef JointPDFType::IndexType             JointPDFIndexType;
  typedef JointPDFType::PixelType             JointPDFValueType;
  typedef JointPDFType::RegionType            JointPDFRegionType;
  typedef JointPDFType::SizeType              JointPDFSizeType;
  typedef JointPDFDerivativesType::IndexType  JointPDFDerivativesIndexType;
  typedef JointPDFDerivativesType::PixelType  JointPDFDerivativesValueType;
  typedef JointPDFDerivativesType::RegionType JointPDFDerivativesRegionType;
  typedef JointPDFDerivativesType::SizeType   JointPDFDerivativesSizeType;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const ParametersType& parameters, DerivativeType & Derivative ) const;

  /**  Get the value and derivatives for single valued optimizers. */
  double GetValueAndDerivative( IndexType index, MeasureType& Value, DerivativeType& Derivative1,
                                DerivativeType& Derivative2 );

  /**  Get the value and derivatives for single valued optimizers. */
  double GetValueAndDerivativeInv( IndexType index, MeasureType& Value, DerivativeType& Derivative1,
                                   DerivativeType& Derivative2 );

  /**  Get the value and derivatives for single valued optimizers. */
  double GetValueAndDerivative2( IndexType index, MeasureType& Value, DerivativeType& Derivative );

  /** Number of spatial samples to used to compute metric */
  //  itkSetClampMacro( NumberOfSpatialSamples, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  // itkGetConstReferenceMacro( NumberOfSpatialSamples, unsigned long);

  /** Number of bins to used in the histogram. Typical value is 50. */
  //  itkSetClampMacro( NumberOfHistogramBins, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  // itkGetConstReferenceMacro( NumberOfHistogramBins, unsigned long);
  void SetNumberOfHistogramBins(unsigned long nhb)
  {
    m_NumberOfHistogramBins = nhb + 2 * this->m_Padding;
  }

  unsigned long GetNumberOfHistogramBins()
  {
    return m_NumberOfHistogramBins;
  }

  void SetNumberOfSpatialSamples(unsigned long nhb)
  {
    m_NumberOfSpatialSamples = nhb;
  }

  unsigned long GetNumberOfSpatialSamples()
  {
    return m_NumberOfSpatialSamples;
  }

  /** Provide API to reinitialize the seed of the random number generator */
  static void ReinitializeSeed();

  static void ReinitializeSeed(int);

  void SetTransform(TransformPointer t)
  {
    m_Transform = t;
  }

  TransformPointer GetTransform()
  {
    return m_Transform;
  }

  void SetInterpolator(InterpolatorPointer t)
  {
    m_Interpolator = t;
  }

  InterpolatorPointer GetInterpolator()
  {
    return m_Interpolator;
  }

  void GetProbabilities();

  void GetProbabilitiesLocal(typename TFixedImage::IndexType centerindex, float radius);

  virtual CovariantVectorType  OpticalFlowUpdate(const NeighborhoodType & neighborhood, bool mov = false)
  {
    // Get fixed image related information
    IndexType index = neighborhood.GetIndex();

    typename TDeformationField::PixelType vec = Superclass::m_DeformationField->GetPixel(index);
    CovariantVectorType update;
    double              fixedValue;
    CovariantVectorType fixedGradient;
    double              fixedGradientSquaredMagnitude = 0;
    fixedValue = (double) Superclass::m_FixedImage->GetPixel( index );
    if( mov )
      {
      fixedGradient = m_MovingImageGradientCalculator->EvaluateAtIndex( index ) * (-1.0);
      }
    fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] );
      }
    double    movingValue;
    int       j;
    PointType mappedPoint;
    for( j = 0; j < ImageDimension; j++ )
      {
      mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j]
        + m_FixedImageOrigin[j];
      mappedPoint[j] += vec[j];
      }
    if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
      movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
      }
    else
      {
      for( j = 0; j < ImageDimension; j++ )
        {
        update[j] = 0.0;
        }
      return update;
      }
    double speedValue;
    if( mov )
      {
      speedValue = movingValue - fixedValue;
      }
    else
      {
      speedValue = fixedValue - movingValue;
      }
    double denominator = vnl_math_sqr( speedValue ) / m_Normalizer
      + fixedGradientSquaredMagnitude;
    double m_DenominatorThreshold = 1e-9;
    double m_IntensityDifferenceThreshold = 0.001;
    if( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold ||
        denominator < m_DenominatorThreshold )
      {
      for( j = 0; j < ImageDimension; j++ )
        {
        update[j] = 0.0;
        }
      return update;
      }
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = speedValue * fixedGradient[j] / denominator;
      }
    return update;
  }

  double ComputeMutualInformation()
  {
    //      typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
    //      JointPDFIteratorType jointPDFIterator ( m_JointPDF, m_JointPDF->GetBufferedRegion() );

    float         px, py, pxy;
    double        mival = 0;
    double        mi;
    unsigned long ct = 0;

    typename JointPDFType::IndexType index;
    for( unsigned int ii = this->m_Padding + 1; ii < m_NumberOfHistogramBins - this->m_Padding - 2; ii++ )
      {
      MarginalPDFIndexType mind;
      mind[0] = ii;
      px = m_FixedImageMarginalPDF->GetPixel(mind);
      for( unsigned int jj = this->m_Padding + 1; jj < m_NumberOfHistogramBins - this->m_Padding - 2; jj++ )
        {
        mind[0] = jj;
        py = m_MovingImageMarginalPDF->GetPixel(mind);
        float denom = px * py;
        index[0] = ii;
        index[1] = jj;
        // pxy=m_JointPDF->GetPixel(index);

        JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer()
          + ( ii * m_NumberOfHistogramBins );
        // Move the pointer to the first affected bin
        int pdfMovingIndex = static_cast<int>( jj );
        pdfPtr += pdfMovingIndex;
        pxy = *(pdfPtr);

        mi = 0;
        if( fabs(denom) > 0 )
          {
          if( pxy / denom > 0 )
            {
            // true mi
            mi = pxy * log(pxy / denom);
            // test mi
            // mi = 1.0 + log(pxy/denom);
            ct++;
            }
          //  std::cout << " II " << ii << " JJ " << jj << " mi " << mi << " mival " << mival << std::endl;
          }

        mival += mi;
        }
      }
    this->m_Energy = mival / ( (float)ct) * (-1.0);
    return this->m_Energy;
  }

  double ComputeLocalMutualInformation(IndexType imageindex, double fixedImageValue, double movingImageValue,
                                       bool deriv = false)
  {
    // first find corresponding parzen indices
    double windowTerm =
      static_cast<double>( fixedImageValue ) / m_FixedImageBinSize
      - m_FixedImageNormalizedMin;
    unsigned int pindex = static_cast<unsigned int>( floor( windowTerm ) );

    // Make sure the extreme values are in valid bins
    if( pindex < this->m_Padding )
      {
      pindex = this->m_Padding;
      }
    else if( pindex > (m_NumberOfHistogramBins -  this->m_Padding - 1) )
      {
      pindex = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).
    double movingImageParzenWindowTerm =
      movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
    unsigned int movingImageParzenWindowIndex =
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

    // Make sure the extreme values are in valid bins
    if( movingImageParzenWindowIndex < this->m_Padding )
      {
      movingImageParzenWindowIndex = this->m_Padding;
      }
    else if( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      movingImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
      }
    unsigned int  ii = pindex;
    unsigned int  jj = movingImageParzenWindowIndex;
    float         px, py, pxy;
    double        mi;
    unsigned long ct = 0;
    typename JointPDFType::IndexType index;
    typename MarginalPDFType::IndexType mind; mind[0] = ii;
    px = m_FixedImageMarginalPDF->GetPixel(mind); mind[0] = jj;
    py = m_MovingImageMarginalPDF->GetPixel(mind);
    float denom = px * py;
    index[0] = ii;
    index[1] = jj;
    pxy = m_JointPDF->GetPixel(index);

    mi = 0;
    if( fabs(denom) > 0 )
      {
      if( pxy / denom > 0 )
        {
        mi = pxy * log(pxy / denom);
        if( deriv )
          {
          mi = 1.0 + log(pxy / denom);
          }
        ct++;
        }
      }

    return mi;
  }

  double ComputeProbabilities(IndexType oindex, unsigned int option = 0)
  {
    double         fixedImageValue = (double)this->Superclass::m_FixedImage->GetPixel(oindex);
    DerivativeType zero(ImageDimension);

    zero.Fill(0);
    double       movingImageValue = this->GetMovingImageValue(oindex, zero);
    unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue);
    unsigned int movingIndex = this->GetMovingValueIndex(movingImageValue);

    double dJPDF = 0, dFmPDF = 0, dMmPDF = 0, jointPDFValue = 0, fixedImagePDFValue = 0, movingImagePDFValue = 0;

    // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
    double movingImageParzenWindowTerm =
      movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
    unsigned int movingImageParzenWindowIndex =
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

    // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).
    double fixedImageParzenWindowTerm =
      fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;
    unsigned int fixedImageParzenWindowIndex =
      static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

    if( movingImageParzenWindowTerm < this->m_Padding )
      {
      movingImageParzenWindowTerm = this->m_Padding;
      }
    else if( movingImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      movingImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    // Make sure the extreme values are in valid bins
    if( fixedImageParzenWindowTerm < this->m_Padding )
      {
      fixedImageParzenWindowTerm = this->m_Padding;
      }
    else if( fixedImageParzenWindowTerm > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      fixedImageParzenWindowTerm = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

      {
      typename JointPDFType::IndexType pdfind2;
      pdfind2[1] = fixedIndex;
      pdfind2[0] = movingIndex;

      typename JointPDFType::PointType pdfind;
      pdfind[1] = fixedImageParzenWindowTerm;
      pdfind[0] = movingImageParzenWindowTerm;
      jointPDFValue = pdfinterpolator->Evaluate(pdfind);
      //	dJPDF = (1.0)*(pdfinterpolator->EvaluateDerivative( pdfind ))[1];

      //    dJPDF += (-0.25)*(pdfinterpolator->EvaluateDerivative( pdfind ))[1];
      // dJPDF = dpdfinterpolator->Evaluate(dpdfind);
      }
      {
      typename MarginalPDFType::PointType mind;
      mind[0] = fixedImageParzenWindowTerm;
      // dFmPDF =(1.0)*(pdfinterpolator2->EvaluateDerivative( mind ))[0];
      fixedImagePDFValue = pdfinterpolator2->Evaluate(mind);
      typename MarginalPDFType::IndexType mind2;
      mind2[0] = fixedIndex;
      //    fixedImagePDFValue = m_FixedImageMarginalPDF->GetPixel(mind2);
      }

      {
      typename MarginalPDFType::PointType mind;
      mind[0] = movingImageParzenWindowTerm;
      movingImagePDFValue = pdfinterpolator3->Evaluate(mind);
      // dMmPDF = (pdfinterpolator3->EvaluateDerivative(mind))[0];

      //    mind[0]=movingImageParzenWindowTerm2;
      // xsdMmPDF = pdfinterpolator3->Evaluate(mind);
      //  mind[0]=movingImageParzenWindowTerm1;
      }

    if( option == 0 )
      {
      return jointPDFValue;
      }
    else if( option == 1 )
      {
      return fixedImagePDFValue;
      }
    else if( option == 2 )
      {
      return movingImagePDFValue;
      }
  }

  virtual VectorType ComputeUpdate(const NeighborhoodType & neighborhood,
                                   void *globalData,
                                   const FloatOffsetType & offset = FloatOffsetType(0.0) )
  {
    VectorType update;

    update.Fill(0.0);
    IndexType oindex = neighborhood.GetIndex();

    FixedImageType* img = const_cast<FixedImageType *>(this->Superclass::m_FixedImage.GetPointer() );
    if( !img )
      {
      return update;
      }
    typename FixedImageType::SpacingType spacing = img->GetSpacing();
    typename FixedImageType::SizeType imagesize = img->GetLargestPossibleRegion().GetSize();
    // bool inimage=true;
    for( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      if( oindex[dd] < 1 ||
          oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1) )
        {
        return update;
        }
      }

    CovariantVectorType fixedGradient;
    //    ImageDerivativesType fixedGradient;
    CovariantVectorType fixedGradientNorm;
    // std::cout << " grad " << std::endl;

    double loce = 0.0;
//    double nccp1=0,nccm1=0;
    ParametersType fdvec1(ImageDimension);
    ParametersType fdvec2(ImageDimension);
    fdvec1.Fill(0);
    fdvec2.Fill(0);

    fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( oindex );
    float mag1 = 0;
    for( int imd = 0; imd < ImageDimension; imd++ )
      {
      mag1 += fixedGradient[imd] * fixedGradient[imd];
      fixedGradientNorm[imd] = fixedGradient[imd];
      }
    mag1 = sqrt(mag1);
    if( mag1 > 1.e-5 )
      {
      fixedGradientNorm /= mag1;
      }

    double nccm1 = 0;
    loce = this->GetValueAndDerivative(oindex, nccm1, fdvec1, fdvec2);
    float sign = (1.) * loce; // nccp1-nccm1;

    double denominator = mag1;
    if( denominator < 1.e-12 )
      {
      denominator = 1.0;
      }
    denominator = 1.;
    for( int imd = 0; imd < ImageDimension; imd++ )
      {
      update[imd] = sign * fixedGradient[imd] / denominator * spacing[imd];
      }

//      if (this->m_MetricImage) this->m_MetricImage->SetPixel(oindex,(loce));//+this->m_MetricImage->GetPixel(oindex));
    if( this->m_MetricImage )
      {
      this->m_MetricImage->SetPixel(oindex, loce);
      }
//	else std::cout << " no " << std::endl;
//  if ( loce < this->m_RobustnessParameter) update.Fill(0);

    if( ImageDimension == 2 )
      {
      if( this->m_MetricImage &&
          (unsigned int)oindex[0] ==
          (unsigned int)this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0] - 2 &&
          (unsigned int)oindex[1] ==
          (unsigned int)this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1] - 2 )
        {
        this->WriteImages();
        }
      }
    /*
    else if (ImageDimension == 3)
      {

  if (//this->m_MetricImage &&
      oindex[0] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0]-2 &&
      oindex[1] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1]-2 &&
      oindex[2] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[2]-2 )
    {
      this->ComputeMutualInformation();
      this->WriteImages();
    }
      }
    */
    mag1 = 0;
    // for (int imd=0; imd<ImageDimension; imd++) mag1+=update[imd]*update[imd];
    //     mag1=sqrt(mag1);
//      if (mag1 > 1.e-5) std::cout << " mag1 " << mag1 << std::endl;//update=update*(1.0/mag1);

    return update;
  }

  virtual VectorType ComputeUpdateInv(const NeighborhoodType & neighborhood,
                                      void *globalData,
                                      const FloatOffsetType & offset = FloatOffsetType(0.0) )
  {
    VectorType update;

    update.Fill(0.0);
    ///   return update;
    IndexType oindex = neighborhood.GetIndex();

    FixedImageType* img = const_cast<FixedImageType *>(this->Superclass::m_MovingImage.GetPointer() );
    if( !img )
      {
      return update;
      }
    typename FixedImageType::SpacingType spacing = img->GetSpacing();
    typename FixedImageType::SizeType imagesize = img->GetLargestPossibleRegion().GetSize();
//    bool inimage=true;
    for( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      if( oindex[dd] < 1 ||
          oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1) )
        {
        return update;
        }
      }

    CovariantVectorType movingGradient;
    //    ImageDerivativesType movingGradient;
    CovariantVectorType movingGradientNorm;
    // std::cout << " grad " << std::endl;

    double loce = 0.0;
//    double nccp1=0,nccm1=0;
    ParametersType fdvec1(ImageDimension);
    ParametersType fdvec2(ImageDimension);
    fdvec1.Fill(0);
    fdvec2.Fill(0);

    movingGradient = m_MovingImageGradientCalculator->EvaluateAtIndex( oindex );
    float mag1 = 0;
    for( int imd = 0; imd < ImageDimension; imd++ )
      {
      mag1 += movingGradient[imd] * movingGradient[imd];
      movingGradientNorm[imd] = movingGradient[imd];
      }
    mag1 = sqrt(mag1);
    if( mag1 > 1.e-5 )
      {
      movingGradientNorm /= mag1;
      }

    //	for (int imd=0; imd<ImageDimension; imd++)
      {
      //  fdvec1[imd]=movingGradientNorm[imd]*spacing[imd]*0.5;
      // fdvec2[imd]=movingGradientNorm[imd]*(-1.)*spacing[imd]*0.5;
      }
    double nccm1 = 0;
    loce = this->GetValueAndDerivativeInv(oindex, nccm1, fdvec1, fdvec2);
    float  sign = (1.0) * loce; // nccp1-nccm1;
    double denominator = mag1;
    if( denominator < 1.e-12 )
      {
      denominator = 1.0;
      }
    denominator = 1.;
    for( int imd = 0; imd < ImageDimension; imd++ )
      {
      update[imd] = sign * movingGradient[imd] / denominator * spacing[imd];
      }

//	if (this->m_MetricImage) this->m_MetricImage->SetPixel(oindex,loce); //this->m_MetricImage->GetPixel(oindex));
//    if ( loce < this->m_RobustnessParameter) update.Fill(0);

/*
    if (ImageDimension == 2)
      {

  if (//this->m_MetricImage &&
      oindex[0] == static_cast<int>(
        this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0]-2 ) &&
      oindex[1] == static_cast<int>(
        this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1]-2 ) )
    {
      //	    this->GetProbabilities();
      this->ComputeMutualInformation();
      this->WriteImages();
    }
      }
    else if (ImageDimension == 3)
      {

  if (//this->m_MetricImage &&
      oindex[0] == static_cast<int>(
        this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0]-2 ) &&
      oindex[1] == static_cast<int>(
        this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1]-2 ) &&
      oindex[2] == static_cast<int>(
        this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[2]-2 ) )
    {
      this->ComputeMutualInformation();
      this->WriteImages();
    }
      }
*/
    //    for (int imd=0; imd<ImageDimension; imd++) mag1+=update[imd]*update[imd];
    //  mag1=sqrt(mag1);
    // if (mag1 > 1.e-5) update=update*(1.0/mag1);

    return update;
  }

  void WriteImages()
  {
    if( this->m_MetricImage )
      {
      typedef ImageFileWriter<FixedImageType> writertype;
      typename writertype::Pointer w = writertype::New();
      w->SetInput(  this->m_MetricImage);
      w->SetFileName("ZZmetric.nii");
      w->Write();
      }
  }

  void SetOpticalFlow(bool b)
  {
    m_OpticalFlow = b;
  }

  typename JointPDFType::Pointer GetJointPDF()
  {
    return m_JointPDF;
  }

  typename JointPDFType::Pointer GetJointHist()
  {
    return m_JointHist;
  }

  void SetFixedImageMask( FixedImageType* img)
  {
    m_FixedImageMask = img;
  }

  void NormalizeJointHist();

  void ChangeJointHist( float oldfixval, float oldmovval, float newfixval, float newmovval)
  {
    // first find corresponding parzen indices
    double windowTerm1 =
      static_cast<double>( oldfixval ) / m_FixedImageBinSize - m_FixedImageNormalizedMin;
    unsigned int pindex1 = static_cast<unsigned int>( floor( windowTerm1 ) );

    // first find corresponding parzen indices
    double windowTerm2 =
      static_cast<double>( newfixval ) / m_FixedImageBinSize - m_FixedImageNormalizedMin;
    unsigned int pindex2 = static_cast<unsigned int>( floor( windowTerm2 ) );

    // Make sure the extreme values are in valid bins
    if( pindex1 < this->m_Padding )
      {
      pindex1 = this->m_Padding;
      }
    else if( pindex1 > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      pindex1 = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    // Make sure the extreme values are in valid bins
    if( pindex2 < this->m_Padding )
      {
      pindex2 = this->m_Padding;
      }
    else if( pindex2 > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      pindex2 = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).
    double movingImageParzenWindowTerm1 =
      oldmovval / m_MovingImageBinSize - m_MovingImageNormalizedMin;
    unsigned int movingImageParzenWindowIndex1 =
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm1 ) );
    // Make sure the extreme values are in valid bins
    if( movingImageParzenWindowIndex1 < this->m_Padding )
      {
      movingImageParzenWindowIndex1 = this->m_Padding;
      }
    else if( movingImageParzenWindowIndex1 > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      movingImageParzenWindowIndex1 = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).
    double movingImageParzenWindowTerm2 =
      newmovval / m_MovingImageBinSize - m_MovingImageNormalizedMin;
    unsigned int movingImageParzenWindowIndex2 =
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm2 ) );
    // Make sure the extreme values are in valid bins
    if( movingImageParzenWindowIndex2 < this->m_Padding )
      {
      movingImageParzenWindowIndex2 = this->m_Padding;
      }
    else if( movingImageParzenWindowIndex2 > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      movingImageParzenWindowIndex2 = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    unsigned int ii1 = pindex1;
    unsigned int jj1 = movingImageParzenWindowIndex1;
    unsigned int ii2 = pindex2;
    unsigned int jj2 = movingImageParzenWindowIndex2;

    typename JointPDFType::IndexType index1;
    index1[0] = ii1;
    index1[1] = jj1;
    float temp1 = m_JointHist->GetPixel(index1);
    m_JointHist->SetPixel(index1, temp1 - 1);

    typename JointPDFType::IndexType index2;
    index2[0] = ii2;
    index2[1] = jj2;
    float temp2 = m_JointHist->GetPixel(index2);
    m_JointHist->SetPixel(index2, temp2 + 1);

    this->NormalizeJointHist();
    this->ComputeMutualInformation();

    // then change back
    temp1 = m_JointHist->GetPixel(index1);
    m_JointHist->SetPixel(index1, temp1 + 1);

    temp2 = m_JointHist->GetPixel(index2);
    m_JointHist->SetPixel(index2, temp2 - 1);
  }

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
    };

  /**
   * A fixed image spatial sample consists of the fixed domain point
   * and the fixed image value at that point. */
  class FixedImageSpatialSample
  {
public:
    FixedImageSpatialSample() : FixedImageValue(0.0)
    {
      FixedImagePointValue.Fill(0.0);
    }

    ~FixedImageSpatialSample()
    {
    };

    IndexType           FixedImageIndex;
    FixedImagePointType FixedImagePointValue;
    double              FixedImageValue;
    unsigned int        FixedImageParzenWindowIndex;
  };

  /** FixedImageSpatialSample typedef support. */
  typedef std::vector<FixedImageSpatialSample>
    FixedImageSpatialSampleContainer;

  /** Container to store a set of points and fixed image values. */
  FixedImageSpatialSampleContainer m_FixedImageSamples;

  /** Uniformly select a sample set from the fixed image domain. */
  virtual void SampleFixedImageDomain(FixedImageSpatialSampleContainer& samples);

  void SampleFixedImageDomainLocal(
    FixedImageSpatialSampleContainer & samples, typename TFixedImage::IndexType index);

  /** Transform a point from FixedImage domain to MovingImage domain.
   * This function also checks if mapped point is within support region. */
  virtual void TransformPoint( unsigned int sampleNumber, MovingImagePointType& mappedPoint,
                               bool& sampleWithinSupportRegion, double& movingImageValue ) const;

  unsigned int GetFixedValueIndex(double fixedImageValue)
  {
    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).
    double fixedImageParzenWindowTerm =
      fixedImageValue / m_FixedImageBinSize - m_FixedImageNormalizedMin;
    unsigned int fixedImageParzenWindowIndex =
      static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

    // Make sure the extreme values are in valid bins
    if( fixedImageParzenWindowIndex < this->m_Padding )
      {
      fixedImageParzenWindowIndex = this->m_Padding;
      }
    else if( fixedImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      fixedImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
      }

    return fixedImageParzenWindowIndex;
  }

  double GetFixedImageValue(IndexType oindex, DerivativeType derivative)
  {
    double    fixedImageValue = 0;
    PointType mappedPoint;

    for( int j = 0; j < ImageDimension; j++ )
      {
      mappedPoint[j] = double( oindex[j] ) * this->Superclass::m_FixedImage->GetSpacing()[j]
        + this->Superclass::m_FixedImage->GetOrigin()[j];
      mappedPoint[j] += derivative[j];
      }
    if( m_FixedImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
      fixedImageValue = m_FixedImageInterpolator->Evaluate( mappedPoint );
      }

    return fixedImageValue;
  }

  double GetMovingImageValue(IndexType oindex, DerivativeType derivative)
  {
    double    movingImageValue = 0;
    PointType mappedPoint;

    for( int j = 0; j < ImageDimension; j++ )
      {
      mappedPoint[j] = double( oindex[j] ) * this->Superclass::m_FixedImage->GetSpacing()[j]
        + this->Superclass::m_FixedImage->GetOrigin()[j];
      mappedPoint[j] += derivative[j];
      }
    if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
      movingImageValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
      }

    return movingImageValue;
  }

  unsigned int GetMovingValueIndex(double movingImageValue)
  {
    double movingImageParzenWindowTerm =
      movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
    unsigned int movingImageParzenWindowIndex =
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

    // Make sure the extreme values are in valid bins
    if( movingImageParzenWindowIndex < this->m_Padding )
      {
      movingImageParzenWindowIndex = this->m_Padding;
      }
    else if( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding - 1) )
      {
      movingImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
      }
    return movingImageParzenWindowIndex;
  }

  void ComputePDFDerivatives( int pdfFixedIndex, int pdfMovingIndex, CovariantVectorType& fixedImageGradientValue,
                              double cubicBSplineDerivativeValue ) const;

  void GetValueAndDerivative3(IndexType oindex, MeasureType& value, DerivativeType& derivative);

  /** The fixed image marginal PDF. */
  mutable MarginalPDFType::Pointer m_FixedImageMarginalPDF;

  /** The moving image marginal PDF. */
  mutable MarginalPDFType::Pointer m_MovingImageMarginalPDF;

  /** The joint PDF and PDF derivatives. */
  typename JointPDFType::Pointer m_JointPDF;

  unsigned long m_NumberOfSpatialSamples;
  unsigned long m_NumberOfParameters;

  /** Variables to define the marginal and joint histograms. */
  unsigned long m_NumberOfHistogramBins;
  double        m_MovingImageNormalizedMin;
  double        m_FixedImageNormalizedMin;
  double        m_MovingImageTrueMin;
  double        m_MovingImageTrueMax;
  double        m_FixedImageBinSize;
  double        m_MovingImageBinSize;
protected:

  AvantsMutualInformationRegistrationFunction();
  virtual ~AvantsMutualInformationRegistrationFunction()
  {
  };
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  AvantsMutualInformationRegistrationFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                              // purposely not implemented

  typename JointPDFType::Pointer m_JointHist;
  typename JointPDFDerivativesType::Pointer m_JointPDFDerivatives;

  /** Typedefs for BSpline kernel and derivative functions. */
  typedef BSplineKernelFunction<3> CubicBSplineFunctionType;
  typedef BSplineDerivativeKernelFunction<3>
    CubicBSplineDerivativeFunctionType;

  /** Cubic BSpline kernel for computing Parzen histograms. */
  typename CubicBSplineFunctionType::Pointer m_CubicBSplineKernel;
  typename CubicBSplineDerivativeFunctionType::Pointer
  m_CubicBSplineDerivativeKernel;

  /** Precompute fixed image parzen window indices. */
  virtual void ComputeFixedImageParzenWindowIndices( FixedImageSpatialSampleContainer& samples );

  /**
   * Types and variables related to image derivative calculations.
   * If a BSplineInterpolationFunction is used, this class obtain
   * image derivatives from the BSpline interpolator. Otherwise,
   * image derivatives are computed using central differencing.
   */
  typedef CovariantVector<double,
                          itkGetStaticConstMacro(ImageDimension)> ImageDerivativesType;

  /** Compute image derivatives at a point. */
  virtual void ComputeImageDerivatives( const MovingImagePointType& mappedPoint, ImageDerivativesType& gradient ) const;

  /** Boolean to indicate if the interpolator BSpline. */
  bool m_InterpolatorIsBSpline;

  // boolean to determine if we use mono-modality assumption
  bool m_OpticalFlow;

  /** Typedefs for using BSpline interpolator. */
  typedef
    BSplineInterpolateImageFunction<MovingImageType,
                                    CoordinateRepresentationType> BSplineInterpolatorType;

  /** Pointer to BSplineInterpolator. */
  typename BSplineInterpolatorType::Pointer m_BSplineInterpolator;

  /** Typedefs for using central difference calculator. */
  typedef CentralDifferenceImageFunction<MovingImageType,
                                         CoordinateRepresentationType> DerivativeFunctionType;

  /** Pointer to central difference calculator. */
  typename DerivativeFunctionType::Pointer m_DerivativeCalculator;

  /**
   * Types and variables related to BSpline deformable transforms.
   * If the transform is of type third order BSplineDeformableTransform,
   * then we can speed up the metric derivative calculation by
   * only inspecting the parameters within the support region
   * of a mapped point.
   */

  /** Boolean to indicate if the transform is BSpline deformable. */
  bool m_TransformIsBSpline;

  /** The number of BSpline parameters per image dimension. */
  long m_NumParametersPerDim;

  /**
   * The number of BSpline transform weights is the number of
   * of parameter in the support region (per dimension ). */
  unsigned long m_NumBSplineWeights;

  /**
   * Enum of the deformabtion field spline order.
   */
  enum { DeformationSplineOrder = 3 };

  /**
   * Typedefs for the BSplineDeformableTransform.
   */
  typedef BSplineDeformableTransform<
      CoordinateRepresentationType,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension,
      DeformationSplineOrder> BSplineTransformType;
  typedef typename BSplineTransformType::WeightsType
    BSplineTransformWeightsType;
  typedef typename BSplineTransformType::ParameterIndexArrayType
    BSplineTransformIndexArrayType;

  /**
   * Variables used when transform is of type BSpline deformable.
   */
  typename BSplineTransformType::Pointer m_BSplineTransform;

  /**
   * Cache pre-transformed points, weights, indices and
   * within support region flag.
   */
  typedef typename BSplineTransformWeightsType::ValueType    WeightsValueType;
  typedef          Array2D<WeightsValueType>                 BSplineTransformWeightsArrayType;
  typedef typename BSplineTransformIndexArrayType::ValueType IndexValueType;
  typedef          Array2D<IndexValueType>                   BSplineTransformIndicesArrayType;
  typedef          std::vector<MovingImagePointType>         MovingImagePointArrayType;
  typedef          std::vector<bool>                         BooleanArrayType;

  BSplineTransformWeightsArrayType m_BSplineTransformWeightsArray;
  BSplineTransformIndicesArrayType m_BSplineTransformIndicesArray;
  MovingImagePointArrayType        m_PreTransformPointsArray;
  BooleanArrayType                 m_WithinSupportRegionArray;

  typename TFixedImage::SpacingType                  m_FixedImageSpacing;
  typename TFixedImage::PointType                  m_FixedImageOrigin;

  typedef FixedArray<unsigned long,
                     ::itk::GetImageDimension<FixedImageType>::ImageDimension> ParametersOffsetType;
  ParametersOffsetType m_ParametersOffset;

  mutable TransformPointer m_Transform;
  InterpolatorPointer      m_Interpolator;
  InterpolatorPointer      m_FixedImageInterpolator;
  InterpolatorPointer      m_MovingImageInterpolator;

  GradientCalculatorPointer m_FixedImageGradientCalculator;
  GradientCalculatorPointer m_MovingImageGradientCalculator;

  FixedImagePointer  m_FixedImageMask;
  MovingImagePointer m_MovingImageMask;

  double m_NormalizeMetric;
  float  m_Normalizer;

  typedef BSplineInterpolateImageFunction<JointPDFType, double> pdfintType;
  typename pdfintType::Pointer pdfinterpolator;

  typedef BSplineInterpolateImageFunction<JointPDFDerivativesType, double> dpdfintType;
  typename dpdfintType::Pointer dpdfinterpolator;

  typedef BSplineInterpolateImageFunction<MarginalPDFType, double> pdfintType2;
  typename pdfintType2::Pointer pdfinterpolator2;
  typename pdfintType2::Pointer pdfinterpolator3;

  unsigned int m_Padding;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAvantsMutualInformationRegistrationFunction.cxx"
#endif

#endif
