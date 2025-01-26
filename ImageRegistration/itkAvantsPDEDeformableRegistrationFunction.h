/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAvantsPDEDeformableRegistrationFunction_h_
#define _itkAvantsPDEDeformableRegistrationFunction_h_
#include "antsAllocImage.h"
#include "itkPDEDeformableRegistrationFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
namespace itk
{
/** \class AvantsPDEDeformableRegistrationFunction
 *
 * This is an abstract base class for all PDE functions which drives a
 * deformable registration algorithm. It is used by
 * PDEDeformationRegistrationFilter subclasses to compute the
 * output deformation field which will map a moving image onto
 * a fixed image.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation field type.
 *
 * \sa AvantsPDEDeformableRegistrationFilter
 * \ingroup PDEDeformableRegistrationFunctions
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
class AvantsPDEDeformableRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef AvantsPDEDeformableRegistrationFunction                                          Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField> Superclass;
  typedef SmartPointer<Self>                                                               Pointer;
  typedef SmartPointer<const Self>                                                         ConstPointer;

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(AvantsPDEDeformableRegistrationFunction);

  /** MovingImage image type. */
  typedef TMovingImage                           MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImagePointer;
  typedef typename TMovingImage::IndexType       IndexType;
  enum
  {
    ImageDimension = MovingImageType::ImageDimension
  };

  /** FixedImage image type. */
  typedef TFixedImage                           FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImagePointer;

  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;

  /** Deformation field type. */
  typedef TDisplacementField                      DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer DisplacementFieldTypePointer;

  typedef typename TDisplacementField::PixelType VectorType;

  /**  PointSet Types */
  typedef itk::PointSet<long, ImageDimension> PointSetType;
  typedef typename PointSetType::Pointer      PointSetPointer;
  typedef typename PointSetType::PointType    PointType;
  typedef typename PointSetType::PixelType    PointDataType;

  typedef Image<float, ImageDimension>                   MetricImageType;
  typedef typename Image<float, ImageDimension>::Pointer MetricImagePointer;

  /** Set the moving image.  */
  void
  SetMovingImage(const MovingImageType * ptr)
  {
    Superclass::m_MovingImage = ptr;
  }

  /** Get the moving image. */
  MovingImageType *
  GetMovingImage(void)
  {
    return const_cast<MovingImageType *>(Superclass::m_MovingImage.GetPointer());
  }

  /** Set the fixed image. */
  void
  SetFixedImage(const FixedImageType * ptr)
  {
    Superclass::m_FixedImage = ptr;
  }

  /** Get the fixed image. */
  FixedImageType *
  GetFixedImage(void)
  {
    return const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer());
  }

  /** Set the fixed image. */
  void
  SetDisplacementField(DisplacementFieldTypePointer ptr)
  {
    Superclass::m_DisplacementField = ptr;
  }

  /** Get the fixed image. */
  DisplacementFieldTypePointer
  GetDisplacementField(void)
  {
    return Superclass::m_DisplacementField;
  }

  void
  SetEnergy(double /* e */)
  {
    this->m_Energy = 0.0;
  }

  double
  GetEnergy()
  {
    return this->m_Energy;
  }

  void
  SetGradientStep(double e)
  {
    this->m_GradientStep = e;
  }

  void
  SetMaxAllowedStep(double e)
  {
    this->m_MaxAllowedStep = e;
  }

  double
  GetGradientStep()
  {
    return this->m_GradientStep;
  }

  void
  SetNormalizeGradient(bool e)
  {
    this->m_NormalizeGradient = e;
  }

  bool
  GetNormalizeGradient()
  {
    return this->m_NormalizeGradient;
  }

  /*
   * allows one to compute the metric everywhere
   */
  void
  ComputeMetricImage()
  {
    bool makenewimage = false;

    FixedImageType *                  img = const_cast<FixedImageType *>(this->m_FixedImage.GetPointer());
    typename FixedImageType::SizeType imagesize = img->GetLargestPossibleRegion().GetSize();

    if (!m_MetricImage)
    {
      makenewimage = true;
    }
    else if (imagesize[0] != m_MetricImage->GetLargestPossibleRegion().GetSize()[0])
    {
      makenewimage = true;
    }

    if (makenewimage)
    {
      m_MetricImage = AllocImage<MetricImageType>(img->GetLargestPossibleRegion(), 0);
      m_MetricImage->SetSpacing(img->GetSpacing());
      m_MetricImage->SetOrigin(img->GetOrigin());
    }
    return;
  }

  virtual double ComputeMetricAtPair(IndexType /* fixedindex */, typename TDisplacementField::PixelType /* vec */)
  {
    return 0.0;
  }

  virtual VectorType
  ComputeUpdateInv(const NeighborhoodType & neighborhood,
                   void * /* globalData */,
                   const FloatOffsetType & /* offset */ = FloatOffsetType(0.0))
  {
    bool       m_Use1SidedDiff = false;
    VectorType update;

    update.Fill(0.0);
    typename FixedImageType::SpacingType   spacing = this->GetFixedImage()->GetSpacing();
    IndexType                              oindex = neighborhood.GetIndex();
    typename TDisplacementField::PixelType vec = this->m_DisplacementField->GetPixel(oindex);
    float                                  loce = 0.0;
    double                                 nccp1, nccm1;
    if (m_Use1SidedDiff)
    {
      nccp1 = this->ComputeMetricAtPair(oindex, vec);
    }
    for (int imd = 0; imd < ImageDimension; imd++)
    {
      typename TDisplacementField::PixelType fdvec1 = this->m_DisplacementField->GetPixel(oindex);
      typename TDisplacementField::PixelType fdvec2 = this->m_DisplacementField->GetPixel(oindex);
      float                                  step = 0.1 * spacing[imd];
      fdvec1[imd] = vec[imd] + step;
      if (!m_Use1SidedDiff)
      {
        nccp1 = this->ComputeMetricAtPair(oindex, fdvec1);
      }
      fdvec2[imd] = vec[imd] - step;
      nccm1 = this->ComputeMetricAtPair(oindex, fdvec2);
      update[imd] = nccp1 - nccm1;
      loce += static_cast<float>(nccp1 + nccm1);
    }
    loce /= (2.0f * static_cast<float>(ImageDimension)); // this->ComputeMetricAtPair(oindex,vec);
    this->m_Energy += static_cast<double>(loce);
    if (m_MetricImage)
    {
      m_MetricImage->SetPixel(oindex, loce);
    }
    return update * this->m_GradientStep;
  }

  VectorType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void * /* globalData */,
                const FloatOffsetType & /* offset */ = FloatOffsetType(0.0)) override
  {
    bool       m_Use1SidedDiff = false;
    VectorType update;

    update.Fill(0.0);
    typename FixedImageType::SpacingType   spacing = this->GetFixedImage()->GetSpacing();
    IndexType                              oindex = neighborhood.GetIndex();
    typename TDisplacementField::PixelType vec = this->m_DisplacementField->GetPixel(oindex);
    float                                  loce = 0.0;
    double                                 nccp1, nccm1;
    if (m_Use1SidedDiff)
    {
      nccp1 = this->ComputeMetricAtPair(oindex, vec);
    }
    for (int imd = 0; imd < ImageDimension; imd++)
    {
      typename TDisplacementField::PixelType fdvec1 = this->m_DisplacementField->GetPixel(oindex);
      typename TDisplacementField::PixelType fdvec2 = this->m_DisplacementField->GetPixel(oindex);
      float                                  step = 0.1 * spacing[imd];
      fdvec1[imd] = vec[imd] + step;
      if (!m_Use1SidedDiff)
      {
        nccp1 = this->ComputeMetricAtPair(oindex, fdvec1);
      }
      fdvec2[imd] = vec[imd] - step;
      nccm1 = this->ComputeMetricAtPair(oindex, fdvec2);
      update[imd] = nccp1 - nccm1;
      loce += static_cast<float>(nccp1 + nccm1);
    }
    loce /= (2.0f * static_cast<float>(ImageDimension)); // this->ComputeMetricAtPair(oindex,vec);
    this->m_Energy += static_cast<double>(loce);
    if (m_MetricImage)
    {
      m_MetricImage->SetPixel(oindex, loce);
    }
    //    float mag=0;
    //    for (int imd=0; imd<ImageDimension; imd++) mag+=update[imd]*update[imd];
    //    if (mag > 1) update.Fill(0.0);
    //      //std::cout << " update " << update << " ind " << oindex << std::endl;
    return update * this->m_GradientStep;
  }

  void
  SetIterations(unsigned int i)
  {
    this->m_Iterations = i;
  }

  /** this parameter is used to set a minimum value for the metric -- a
     minimum value that is treated as " close enough" --  */
  void
  SetRobustnessParameter(float i)
  {
    this->m_RobustnessParameter = i;
  }

  void
  SetFixedPointSet(PointSetPointer p)
  {
    this->m_FixedPointSet = p;
  }

  void
  SetMovingPointSet(PointSetPointer p)
  {
    this->m_MovingPointSet = p;
  }

  bool
  ThisIsAPointSetMetric()
  {
    return this->m_IsPointSetMetric;
  }

protected:
  AvantsPDEDeformableRegistrationFunction()
  {
    this->m_MovingImage = nullptr;
    m_MetricImage = nullptr;
    this->m_FixedImage = nullptr;
    this->m_DisplacementField = nullptr;
    this->m_Energy = 0.0;
    m_BestEnergy = 0.0;
    this->m_NormalizeGradient = true;
    this->m_GradientStep = 1.0;
    m_AverageStepMag = 0;
    m_MaxStepMag = 0;
    m_AvgCt = 0;
    m_Iterations = 0;

    this->m_FixedPointSet = nullptr;
    this->m_MovingPointSet = nullptr;
    this->m_IsPointSetMetric = false;
    this->m_RobustnessParameter = -1.e12;
  }

  ~AvantsPDEDeformableRegistrationFunction() override = default;

  mutable double m_BestEnergy;
  mutable double m_LastLastEnergy;
  mutable double m_MaxAllowedStep;

  mutable double        m_AverageStepMag;
  mutable unsigned long m_AvgCt;
  mutable unsigned long m_Iterations;
  mutable double        m_MaxStepMag;

  PointSetPointer m_FixedPointSet;
  PointSetPointer m_MovingPointSet;
  bool            m_IsPointSetMetric;

  MetricImagePointer m_MetricImage;

  float m_RobustnessParameter;

private:
  AvantsPDEDeformableRegistrationFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};
} // end namespace itk

#endif
