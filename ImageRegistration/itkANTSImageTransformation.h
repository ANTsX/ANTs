/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkANTSImageTransformation_h
#define __itkANTSImageTransformation_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "antsCommandLineParser.h"
#include "itkImage.h"
#include "itkMacro.h"

#include "itkCenteredEuler3DTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkANTSAffine3DTransform.h"
// #include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
namespace itk
{
template <unsigned int TDimension = 3, typename TReal = float>
class ANTSImageTransformation final : public Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSImageTransformation  Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef double TComp;
  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ANTSImageTransformation);
  static constexpr unsigned int Dimension = TDimension;
  static constexpr unsigned int ImageDimension = TDimension;

  typedef TReal                                              RealType;
  typedef Image<RealType, Self::Dimension> ImageType;
  /** declare transformation types */

  typedef itk::MatrixOffsetTransformBase<TComp, TDimension, TDimension> AffineTransformType;

  typedef typename AffineTransformType::Pointer      AffineTransformPointer;
  typedef itk::Vector<TReal, ImageDimension>         VectorType;
  typedef itk::Image<VectorType, ImageDimension>     DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer    DisplacementFieldPointer;
  typedef typename DisplacementFieldType::RegionType DeformationRegionOfInterestType;
  typedef typename DisplacementFieldType::SizeType   DeformationRegionOfInterestSizeType;
  typedef typename DisplacementFieldType::PointType  DeformationRegionOfInterestCenterType;

  typedef typename ants::CommandLineParser ParserType;
  typedef typename ParserType::OptionType  OptionType;

  /** Set functions */
  void
  SetDeformationRegionOfInterest(DeformationRegionOfInterestCenterType       DRC,
                                 DeformationRegionOfInterestSizeType         DRS,
                                 typename DisplacementFieldType::SpacingType DRSp)
  {
    m_DeformationRegionOfInterestCenter = DRC;
    m_DeformationRegionOfInterestSize = DRS;
    m_DeformationRegionSpacing = DRSp;
  }

  void
  SetAffineTransform(AffineTransformPointer A)
  {
    this->m_AffineTransform = A;
  }

  void
  SetDisplacementField(DisplacementFieldPointer A)
  {
    this->m_DisplacementField = A;
  }

  void
  SetInverseDisplacementField(DisplacementFieldPointer A)
  {
    this->m_InverseDisplacementField = A;
  }

  void
  SetNamingConvention(std::string name)
  {
    this->m_NamingConvention = name;
  }

  /** Get functions */
  AffineTransformPointer
  GetAffineTransform()
  {
    return this->m_AffineTransform;
  }

  DisplacementFieldPointer
  GetDisplacementField()
  {
    return this->m_DisplacementField;
  }

  DisplacementFieldPointer
  GetInverseDisplacementField()
  {
    return this->m_InverseDisplacementField;
  }

  void
  SetFixedImageAffineTransform(AffineTransformPointer A)
  {
    this->m_FixedImageAffineTransform = A;
  }

  AffineTransformPointer
  GetFixedImageAffineTransform()
  {
    return this->m_FixedImageAffineTransform;
  }

  /** Initialize the mapping */
  void
  InitializeTransform()
  {
    if (!this->m_AffineTransform)
    {
      this->m_AffineTransform = AffineTransformType::New();
      this->m_AffineTransform->SetIdentity();
    }
    // deformation fields too
    if (!this->m_DisplacementField)
    {
      VectorType zero;
      zero.Fill(0);
      this->m_DisplacementField = DisplacementFieldType::New();
      /*      m_DeformationRegionOfInterest.SetSize( m_DeformationRegionOfInterestSize );
            this->m_DisplacementField->SetSpacing( m_DeformationRegionSpacing );
            this->m_DisplacementField->SetOrigin( m_DeformationRegionOfInterestCenter );
            this->m_DisplacementField->SetLargestPossibleRegion( m_DeformationRegionOfInterest );
            this->m_DisplacementField->SetRequestedRegion( m_DeformationRegionOfInterest );
            this->m_DisplacementField->SetBufferedRegion( m_DeformationRegionOfInterest );
            this->m_DisplacementField->AllocateInitialized();
      */
    }
  }

  /** Write the transformations out */
  void
  Write();

  /** Concatenate all transformations  */
  void
  Compose();

  /** Concatenate all transformations in inverse direction */
  void
  ComposeInverse();

  itkSetMacro(WriteComponentImages, bool);
  itkGetMacro(WriteComponentImages, bool);
  itkBooleanMacro(WriteComponentImages);

protected:
  ANTSImageTransformation();
  ~ANTSImageTransformation() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  ANTSImageTransformation(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  AffineTransformPointer                      m_AffineTransform;
  AffineTransformPointer                      m_FixedImageAffineTransform;
  DisplacementFieldPointer                    m_DisplacementField;
  DeformationRegionOfInterestType             m_DeformationRegionOfInterest;
  DeformationRegionOfInterestCenterType       m_DeformationRegionOfInterestCenter;
  DeformationRegionOfInterestSizeType         m_DeformationRegionOfInterestSize;
  typename DisplacementFieldType::SpacingType m_DeformationRegionSpacing;
  DisplacementFieldPointer                    m_InverseDisplacementField;
  std::string                                 m_NamingConvention;
  bool                                        m_WriteComponentImages;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSImageTransformation.cxx"
#endif

#endif
