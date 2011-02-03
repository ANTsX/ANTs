/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkANTSImageTransformation.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

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
#include  "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
namespace itk
{
template <unsigned int TDimension = 3, class TReal = float>
class ITK_EXPORT ANTSImageTransformation
  : public       Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSImageTransformation  Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef double TComp;
  /** Run-time type information (and related methods). */
  itkTypeMacro( ANTSImageTransformation, Object );
  itkStaticConstMacro( Dimension, unsigned int, TDimension );
  itkStaticConstMacro( ImageDimension, unsigned int, TDimension );

  typedef TReal RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( Dimension )>                   ImageType;
  /** declare transformation types */

  typedef itk::MatrixOffsetTransformBase<TComp, TDimension, TDimension> AffineTransformType;

  typedef typename AffineTransformType::Pointer     AffineTransformPointer;
  typedef itk::Vector<TReal, ImageDimension>        VectorType;
  typedef itk::Image<VectorType, ImageDimension>    DeformationFieldType;
  typedef typename DeformationFieldType::Pointer    DeformationFieldPointer;
  typedef typename DeformationFieldType::RegionType DeformationRegionOfInterestType;
  typedef typename DeformationFieldType::SizeType   DeformationRegionOfInterestSizeType;
  typedef typename DeformationFieldType::PointType  DeformationRegionOfInterestCenterType;

  typedef typename ants::CommandLineParser ParserType;
  typedef typename ParserType::OptionType  OptionType;

  /** Set functions */
  void SetDeformationRegionOfInterest( DeformationRegionOfInterestCenterType DRC,
                                       DeformationRegionOfInterestSizeType DRS,
                                       typename DeformationFieldType::SpacingType DRSp)
  {
    m_DeformationRegionOfInterestCenter = DRC;
    m_DeformationRegionOfInterestSize = DRS;
    m_DeformationRegionSpacing = DRSp;
  }

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

  void SetNamingConvention(std::string name)
  {
    this->m_NamingConvention = name;
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

  void SetFixedImageAffineTransform(AffineTransformPointer A)
  {
    this->m_FixedImageAffineTransform = A;
  }

  AffineTransformPointer GetFixedImageAffineTransform()
  {
    return this->m_FixedImageAffineTransform;
  }

  /** Initialize the mapping */
  void InitializeTransform()
  {
    if( !this->m_AffineTransform )
      {
      this->m_AffineTransform = AffineTransformType::New();
      this->m_AffineTransform->SetIdentity();
      }
    // deformation fields too
    if( !this->m_DeformationField )
      {
      VectorType zero;
      zero.Fill(0);
      this->m_DeformationField = DeformationFieldType::New();
/*      m_DeformationRegionOfInterest.SetSize( m_DeformationRegionOfInterestSize );
      this->m_DeformationField->SetSpacing( m_DeformationRegionSpacing );
      this->m_DeformationField->SetOrigin( m_DeformationRegionOfInterestCenter );
      this->m_DeformationField->SetLargestPossibleRegion( m_DeformationRegionOfInterest );
      this->m_DeformationField->SetRequestedRegion( m_DeformationRegionOfInterest );
      this->m_DeformationField->SetBufferedRegion( m_DeformationRegionOfInterest );
      this->m_DeformationField->Allocate();
      this->m_DeformationField->FillBuffer(zero);
*/
      }
  }

  /** Write the transformations out */
  void Write();

  /** Concatenate all transformations  */
  void Compose();

  /** Concatenate all transformations in inverse direction */
  void ComposeInverse();

  itkSetMacro( WriteComponentImages, bool );
  itkGetMacro( WriteComponentImages, bool );
  itkBooleanMacro( WriteComponentImages );
protected:
  ANTSImageTransformation();
  virtual ~ANTSImageTransformation()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  ANTSImageTransformation( const Self & ); // purposely not implemented
  void operator=( const Self & );          // purposely not implemented

  AffineTransformPointer                m_AffineTransform;
  AffineTransformPointer                m_FixedImageAffineTransform;
  DeformationFieldPointer               m_DeformationField;
  DeformationRegionOfInterestType       m_DeformationRegionOfInterest;
  DeformationRegionOfInterestCenterType m_DeformationRegionOfInterestCenter;
  DeformationRegionOfInterestSizeType   m_DeformationRegionOfInterestSize;
  typename DeformationFieldType::SpacingType m_DeformationRegionSpacing;
  DeformationFieldPointer m_InverseDeformationField;
  std::string             m_NamingConvention;
  bool                    m_WriteComponentImages;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkANTSImageTransformation.cxx"
#endif

#endif
