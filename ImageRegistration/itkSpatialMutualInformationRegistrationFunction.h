/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSpatialMutualInformationRegistrationFunction_h
#define __itkSpatialMutualInformationRegistrationFunction_h

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
#include "itkTranslationTransform.h"
#include "itkArray2D.h"
#include "itkImageBase.h"
#include "itkTransform.h"
#include "itkInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkMacro.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSpatialObject.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
/** \class SpatialMutualInformationRegistrationFunction
 * \brief Computes the mutual information between two images to be
 * registered using the method of Spatial et al.
 *
 * SpatialMutualInformationRegistrationFunction computes the mutual
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
 * The calculations are based on the method of Spatial et al [1,2]
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
 *      D. Spatial, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      Medical Imaging 2001: Image Processing, 2001, pp. 1609-1620.
 * [2] "PET-CT Image Registration in the Chest Using Free-form Deformations"
 *      D. Spatial, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
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
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
class SpatialMutualInformationRegistrationFunction final
  : public AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef SpatialMutualInformationRegistrationFunction                                           Self;
  typedef AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField> Superclass;
  typedef SmartPointer<Self>                                                                     Pointer;
  typedef SmartPointer<const Self>                                                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SpatialMutualInformationRegistrationFunction);

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
  typedef typename Superclass::VectorType                   VectorType;
  typedef typename Superclass::DisplacementFieldType        DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldTypePointer DisplacementFieldTypePointer;

  /** Inherit some enums from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Interpolator type. */
  typedef double CoordRepType;
  typedef //       //    LinearInterpolateImageFunction<MovingImageType,CoordRepType>
    BSplineInterpolateImageFunction<MovingImageType, CoordRepType>
                                               InterpolatorType;
  typedef typename InterpolatorType::Pointer   InterpolatorPointer;
  typedef typename InterpolatorType::PointType PointType;
  typedef InterpolatorType                     DefaultInterpolatorType;
  //  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
  // DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double, Self::ImageDimension> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void
  SetMovingImageInterpolator(InterpolatorType * ptr)
  {
    m_MovingImageInterpolator = ptr;
  }

  /** Get the moving image interpolator. */
  InterpolatorType *
  GetMovingImageInterpolator(void)
  {
    return m_MovingImageInterpolator;
  }

  /** This class uses a constant timestep of 1. */
  TimeStepType
  ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const override
  {
    return 1;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  void *
  GetGlobalDataPointer() const override
  {
    GlobalDataStruct * global = new GlobalDataStruct();

    //    global->m_SumOfSquaredDifference  = 0.0;
    // / global->m_NumberOfPixelsProcessed = 0L;
    // global->m_SumOfSquaredChange      = 0;
    return global;
  }

  /** Release memory for global data structure. */
  void
  ReleaseGlobalDataPointer(void * GlobalData) const override
  {
    delete (GlobalDataStruct *)GlobalData;
  }

  /** Set the object's state before each iteration. */
  void
  InitializeIteration() override;

  typedef double CoordinateRepresentationType;

  /** Types inherited from Superclass. */
  typedef TranslationTransform<CoordinateRepresentationType,
                               //                    Self::ImageDimension,
                               Self::ImageDimension>
    TransformType;

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

  /**  Get the value and derivatives for single valued optimizers. */
  double
  GetValueAndDerivative(IndexType        index,
                        MeasureType &    Value,
                        DerivativeType & Derivative1,
                        DerivativeType & Derivative2);

  /**  Get the value and derivatives for single valued optimizers. */
  double
  GetValueAndDerivativeInv(IndexType        index,
                           MeasureType &    Value,
                           DerivativeType & Derivative1,
                           DerivativeType & Derivative2);

  /** Number of spatial samples to used to compute metric */
  //  itkSetClampMacro( NumberOfSpatialSamples, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  // itkGetConstReferenceMacro( NumberOfSpatialSamples, unsigned long);

  /** Number of bins to used in the histogram. Typical value is 50. */
  //  itkSetClampMacro( NumberOfHistogramBins, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  // itkGetConstReferenceMacro( NumberOfHistogramBins, unsigned long);
  void
  SetNumberOfHistogramBins(unsigned long nhb)
  {
    m_NumberOfHistogramBins = nhb + 2 * this->m_Padding;
  }

  unsigned long
  GetNumberOfHistogramBins()
  {
    return m_NumberOfHistogramBins;
  }

  void
  SetTransform(TransformPointer t)
  {
    m_Transform = t;
  }

  TransformPointer
  GetTransform()
  {
    return m_Transform;
  }

  void
  SetInterpolator(InterpolatorPointer t)
  {
    m_Interpolator = t;
  }

  InterpolatorPointer
  GetInterpolator()
  {
    return m_Interpolator;
  }

  void
  GetProbabilities();

  double
  ComputeMutualInformation()
  {
    //      typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
    //      JointPDFIteratorType jointPDFIterator ( m_JointPDF, m_JointPDF->GetBufferedRegion() );

    float         px, py, pxy;
    double        mival = 0;
    double        mi;
    // unsigned long ct = 0;

    typename JointPDFType::IndexType index;
    //    for (unsigned int ii=this->m_Padding+1; ii<m_NumberOfHistogramBins-this->m_Padding-2; ii++)
    for (unsigned int ii = 0; ii < m_NumberOfHistogramBins; ii++)
    {
      MarginalPDFIndexType mind;
      mind[0] = ii;
      px = m_FixedImageMarginalPDF->GetPixel(mind);
      //  for (unsigned int jj=this->m_Padding+1; jj<m_NumberOfHistogramBins-this->m_Padding-2; jj++)
      for (unsigned int jj = 0; jj < m_NumberOfHistogramBins; jj++)
      {
        mind[0] = jj;
        py = m_MovingImageMarginalPDF->GetPixel(mind);
        float denom = px * py;
        index[0] = ii;
        index[1] = jj;
        // pxy=m_JointPDF->GetPixel(index);

        JointPDFValueType * pdfPtr = m_JointPDF->GetBufferPointer() + (ii * m_NumberOfHistogramBins);
        // Move the pointer to the first affected bin
        int pdfMovingIndex = static_cast<int>(jj);
        pdfPtr += pdfMovingIndex;
        pxy = *(pdfPtr);

        mi = 0;
        if (fabs(denom) > 0)
        {
          if (pxy / denom > 0)
          {
            // true mi
            mi = static_cast<double>(pxy * std::log(pxy / denom));
            // test mi
            // mi = 1.0 + log(pxy/denom);
            // ct++;
          }
        }

        mival += mi;
      }
      //      std::cout << " II " << ii << " JJ " << ii << " pxy " << pxy << " px " << px << std::endl;
    }
    this->m_Energy = (-1.0) * mival / log(2);
    return this->m_Energy;
  }

  double
  ComputeSpatialMutualInformation()
  {
    float  pxy;
    double SMI = 0;
    double JointEntropy = 0, JointEntropyXuY = 0, JointEntropyXYu = 0, JointEntropyXlY = 0, JointEntropyXYl = 0;
    double JointEntropyXuYl = 0, JointEntropyXlYu = 0, JointEntropyXrYu = 0, JointEntropyXuYr = 0;

    for (unsigned int ii = 0; ii < m_NumberOfHistogramBins; ii++)
    {
      for (unsigned int jj = 0; jj < m_NumberOfHistogramBins; jj++)
      {
        int                 pdfMovingIndex = static_cast<int>(jj);
        JointPDFValueType * pdfPtr = m_JointPDF->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropy -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXuY->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXuY -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXYu->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXYu -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXlY->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXlY -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXYl->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXYl -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXuYl->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXuYl -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXlYu->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXlYu -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXrYu->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXrYu -= static_cast<double>(pxy * std::log(pxy));
        }

        pdfPtr = m_JointPDFXuYr->GetBufferPointer() + (ii * m_NumberOfHistogramBins) + pdfMovingIndex;
        pxy = *(pdfPtr);
        if (std::fabs(pxy) > itk::NumericTraits<float>::ZeroValue())
        {
          JointEntropyXuYr -= static_cast<double>(pxy * std::log(pxy));
        }
      }
    }
    SMI = (0.5) * (JointEntropyXuY + JointEntropyXYu + JointEntropyXlY + JointEntropyXYl) -
          (0.25) * (4 * JointEntropy + JointEntropyXuYr + JointEntropyXrYu + JointEntropyXlYu + JointEntropyXuYl);

    //    std::cout << " JE " << JointEntropy << " JEXuY " << JointEntropyXuY << " JEXYu " << JointEntropyXYu <<
    // "
    // JEXuYr " << JointEntropyXuYr << std::endl;

    this->m_Energy = (-1.0) * fabs(SMI);
    return this->m_Energy;
  }

  VectorType
  ComputeUpdateInv(const NeighborhoodType & neighborhood,
                   void * /* globalData */,
                   const FloatOffsetType & /* offset */ = FloatOffsetType(0.0)) override
  {
    VectorType update;

    update.Fill(0.0);
    IndexType oindex = neighborhood.GetIndex();

    FixedImageType * img = const_cast<FixedImageType *>(this->Superclass::m_FixedImage.GetPointer());
    if (!img)
    {
      return update;
    }
    typename FixedImageType::SpacingType spacing = img->GetSpacing();
    typename FixedImageType::SizeType    imagesize = img->GetLargestPossibleRegion().GetSize();
    for (unsigned int dd = 0; dd < ImageDimension; dd++)
    {
      if (oindex[dd] < 1 || oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1))
      {
        return update;
      }
    }

    CovariantVectorType fixedGradient;
    ParametersType      fdvec1(ImageDimension);
    ParametersType      fdvec2(ImageDimension);
    fdvec1.Fill(0);
    fdvec2.Fill(0);
    fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex(oindex);
    double nccm1 = 0;
    double loce = this->GetValueAndDerivativeInv(oindex, nccm1, fdvec1, fdvec2);
    double eps = 10.0;
    if (loce > eps)
    {
      loce = eps;
    }
    if (loce < -eps)
    {
      loce = -eps;
    }
    for (unsigned int imd = 0; imd < ImageDimension; imd++)
    {
      update[imd] = loce * fixedGradient[imd] * spacing[imd];
    }
    if (this->m_MetricImage)
    {
      this->m_MetricImage->SetPixel(oindex, loce);
    }
    return update;
  }

  /*   Normalizing the image to the range of [0 1] */
  double
  GetMovingParzenTerm(double intensity)
  {
    double windowTerm = static_cast<double>(intensity) - this->m_MovingImageTrueMin;

    windowTerm = windowTerm / (this->m_MovingImageTrueMax - this->m_MovingImageTrueMin);
    return windowTerm;
  }

  double
  GetFixedParzenTerm(double intensity)
  {
    double windowTerm = static_cast<double>(intensity) - this->m_FixedImageTrueMin;

    windowTerm = windowTerm / (this->m_FixedImageTrueMax - this->m_FixedImageTrueMin);
    return windowTerm;
  }

  /* find the image index in the number of histogram bins */
  unsigned int
  FitIndexInBins(double windowTerm)
  {
    unsigned int movingImageParzenWindowIndex = static_cast<unsigned int>(
      this->m_Padding +
      static_cast<unsigned int>(static_cast<float>(windowTerm) *
                                  static_cast<float>(this->m_NumberOfHistogramBins - 1 - this->m_Padding) +
                                0.5f));

    // Make sure the extreme values are in valid bins
    if (movingImageParzenWindowIndex < this->m_Padding)
    {
      movingImageParzenWindowIndex = this->m_Padding - 1;
    }
    else if (movingImageParzenWindowIndex > (m_NumberOfHistogramBins - this->m_Padding))
    {
      movingImageParzenWindowIndex = m_NumberOfHistogramBins - this->m_Padding - 1;
    }

    return movingImageParzenWindowIndex;
  }

  double
  FitContIndexInBins(double windowTerm)
  {
    return (static_cast<double>(this->m_Padding) +
            windowTerm * static_cast<double>(this->m_NumberOfHistogramBins - this->m_Padding));
  }

  VectorType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void * /* globalData */,
                const FloatOffsetType & /* offset */ = FloatOffsetType(0.0)) override
  {
    VectorType update;

    update.Fill(0.0);
    IndexType oindex = neighborhood.GetIndex();

    FixedImageType * img = const_cast<FixedImageType *>(this->Superclass::m_MovingImage.GetPointer());
    if (!img)
    {
      return update;
    }
    typename FixedImageType::SpacingType spacing = img->GetSpacing();
    typename FixedImageType::SizeType    imagesize = img->GetLargestPossibleRegion().GetSize();
    for (unsigned int dd = 0; dd < ImageDimension; dd++)
    {
      if (oindex[dd] < 1 || oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1))
      {
        return update;
      }
    }

    CovariantVectorType movingGradient;
    ParametersType      fdvec1(ImageDimension);
    ParametersType      fdvec2(ImageDimension);
    fdvec1.Fill(0);
    fdvec2.Fill(0);
    movingGradient = m_MovingImageGradientCalculator->EvaluateAtIndex(oindex);

    double nccm1 = 0;
    double loce = this->GetValueAndDerivative(oindex, nccm1, fdvec1, fdvec2);

    double eps = 10;
    if (loce > eps)
    {
      loce = eps;
    }
    if (loce < -eps)
    {
      loce *= -1.0;
    }
    for (unsigned int imd = 0; imd < ImageDimension; imd++)
    {
      update[imd] = static_cast<typename VectorType::ComponentType>(loce * movingGradient[imd] * spacing[imd]);
    }
    return update;
  }

  void
  WriteImages()
  {
    if (this->m_MetricImage)
    {
      typedef ImageFileWriter<FixedImageType> writertype;
      typename writertype::Pointer            w = writertype::New();
      w->SetInput(this->m_MetricImage);
      w->SetFileName("ZZmetric.nii.gz");
      w->Write();
    }
  }

  void
  SetOpticalFlow(bool b)
  {
    m_OpticalFlow = b;
  }

  typename JointPDFType::Pointer
  GetJointPDF()
  {
    return m_JointPDF;
  }

  typename JointPDFType::Pointer
  GetJointHist()
  {
    return m_JointHist;
  }

  void
  SetFixedImageMask(FixedImageType * img)
  {
    m_FixedImageMask = img;
  }

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
  {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
  };

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
  double        m_FixedImageTrueMin;
  double        m_FixedImageTrueMax;
  double        m_FixedImageBinSize;
  double        m_MovingImageBinSize;

protected:
  SpatialMutualInformationRegistrationFunction();
  virtual ~SpatialMutualInformationRegistrationFunction() override{};
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  SpatialMutualInformationRegistrationFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  typename JointPDFType::Pointer            m_JointHist;
  typename JointPDFDerivativesType::Pointer m_JointPDFDerivatives;

  typedef BSplineInterpolateImageFunction<JointPDFType, double> pdfintType;
  typename pdfintType::Pointer                                  pdfinterpolator;

  //  typename JointPDFType::Pointer m_JointEntropy;
  typename JointPDFType::Pointer m_JointPDFXuY;
  typename JointPDFType::Pointer m_JointPDFXYu;
  typename JointPDFType::Pointer m_JointPDFXlY;
  typename JointPDFType::Pointer m_JointPDFXYl;
  typename JointPDFType::Pointer m_JointPDFXuYl;
  typename JointPDFType::Pointer m_JointPDFXlYu;
  typename JointPDFType::Pointer m_JointPDFXrYu;
  typename JointPDFType::Pointer m_JointPDFXuYr;

  typename pdfintType::Pointer pdfinterpolatorXuY;
  typename pdfintType::Pointer pdfinterpolatorXYu;
  typename pdfintType::Pointer pdfinterpolatorXlY;
  typename pdfintType::Pointer pdfinterpolatorXYl;
  typename pdfintType::Pointer pdfinterpolatorXuYl;
  typename pdfintType::Pointer pdfinterpolatorXlYu;
  typename pdfintType::Pointer pdfinterpolatorXuYr;
  typename pdfintType::Pointer pdfinterpolatorXrYu;

  /** Typedefs for BSpline kernel and derivative functions. */
  typedef BSplineKernelFunction<3>           CubicBSplineFunctionType;
  typedef BSplineDerivativeKernelFunction<3> CubicBSplineDerivativeFunctionType;

  /** Cubic BSpline kernel for computing Parzen histograms. */
  typename CubicBSplineFunctionType::Pointer           m_CubicBSplineKernel;
  typename CubicBSplineDerivativeFunctionType::Pointer m_CubicBSplineDerivativeKernel;

  /**
   * Types and variables related to image derivative calculations.
   * If a BSplineInterpolationFunction is used, this class obtain
   * image derivatives from the BSpline interpolator. Otherwise,
   * image derivatives are computed using central differencing.
   */
  typedef CovariantVector<double, Self::ImageDimension> ImageDerivativesType;

  /** Boolean to indicate if the interpolator BSpline. */
  bool m_InterpolatorIsBSpline;

  // boolean to determine if we use mono-modality assumption
  bool m_OpticalFlow;

  /** Typedefs for using BSpline interpolator. */
  typedef BSplineInterpolateImageFunction<MovingImageType, CoordinateRepresentationType> BSplineInterpolatorType;

  /** Pointer to BSplineInterpolator. */
  typename BSplineInterpolatorType::Pointer m_BSplineInterpolator;

  /** Typedefs for using central difference calculator. */
  typedef CentralDifferenceImageFunction<MovingImageType, CoordinateRepresentationType> DerivativeFunctionType;

  /** Pointer to central difference calculator. */
  typename DerivativeFunctionType::Pointer m_DerivativeCalculator;

  /**
   * Types and variables related to BSpline deformable transforms.
   * If the transform is of type third order BSplineTransform,
   * then we can speed up the metric derivative calculation by
   * only inspecting the parameters within the support region
   * of a mapped point.
   */

  /** Boolean to indicate if the transform is BSpline deformable. */

  /** The number of BSpline parameters per image dimension. */
  long m_NumParametersPerDim;

  /**
   * The number of BSpline transform weights is the number of
   * of parameter in the support region (per dimension ). */
  unsigned long m_NumBSplineWeights;

  /**
   * Enum of the deformabtion field spline order.
   */
  enum
  {
    DeformationSplineOrder = 3
  };

  typename TFixedImage::SpacingType m_FixedImageSpacing;
  typename TFixedImage::PointType   m_FixedImageOrigin;

  typedef FixedArray<unsigned long, FixedImageType::ImageDimension> ParametersOffsetType;
  ParametersOffsetType                                              m_ParametersOffset;

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

  //  typedef BSplineInterpolateImageFunction<JointPDFDerivativesType,double> dpdfintType;
  // typename dpdfintType::Pointer dpdfinterpolator;

  typedef BSplineInterpolateImageFunction<MarginalPDFType, double> pdfintType2;
  typename pdfintType2::Pointer                                    pdfinterpolator2;
  typename pdfintType2::Pointer                                    pdfinterpolator3;

  unsigned int m_Padding;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpatialMutualInformationRegistrationFunction.cxx"
#endif

#endif
