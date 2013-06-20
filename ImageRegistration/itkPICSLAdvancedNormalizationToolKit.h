/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkPICSLAdvancedNormalizationToolKit.h,v $
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
#ifndef __itkPICSLAdvancedNormalizationToolKit_h
#define __itkPICSLAdvancedNormalizationToolKit_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkANTSImageTransformation.h"
#include "itkANTSImageRegistrationOptimizer.h"
#include "itkANTSSimilarityMetric.h"
#include "antsCommandLineParser.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "itkANTSLabeledPointSet.h"

namespace itk
{
template <unsigned int TDimension = 3, class TReal = float>
class PICSLAdvancedNormalizationToolKit
  : public       Object
{
public:
  /** Standard class typedefs. */
  typedef PICSLAdvancedNormalizationToolKit Self;
  typedef Object                            Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( PICSLAdvancedNormalizationToolKit, Object );
  itkStaticConstMacro( Dimension, unsigned int, TDimension );
  typedef double TComp;
  typedef TReal  RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( Dimension )>                   ImageType;
  typedef typename ImageType::Pointer   ImagePointer;
  typedef typename ImageType::PixelType PixelType;

  typedef itk::ANTSImageTransformation<Dimension, TReal>          TransformationModelType;
  typedef typename TransformationModelType::Pointer               TransformationModelPointer;
  typedef itk::ANTSImageRegistrationOptimizer<Dimension, TReal>   RegistrationOptimizerType;
  typedef typename RegistrationOptimizerType::Pointer             RegistrationOptimizerPointer;
  typedef typename TransformationModelType::DisplacementFieldType DisplacementFieldType;
  typedef typename TransformationModelType::AffineTransformType   AffineTransformType;
  typedef typename RegistrationOptimizerType::OptAffineType       OptAffineType;

  /** Point Set Type */
  typedef itk::ANTSLabeledPointSet<Dimension>        LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer      LabeledPointSetPointer;
  typedef typename LabeledPointSetType::PointSetType PointSetType;

  /** Typedefs for similarity metrics */
  typedef ANTSSimilarityMetric<itkGetStaticConstMacro( Dimension ), TReal> SimilarityMetricType;
  typedef typename SimilarityMetricType::Pointer                           SimilarityMetricPointer;
  typedef std::vector<SimilarityMetricPointer>                             SimilarityMetricListType;
  typedef typename SimilarityMetricType::MetricType                        MetricBaseType;

  typedef ants::CommandLineParser         ParserType;
  typedef typename ParserType::OptionType OptionType;

  void ParseCommandLine( int argc, char * *argv );

  TransformationModelPointer GetTransformationModel()
  {
    return this->m_TransformationModel;
  }

  RegistrationOptimizerPointer  SetRegistrationOptimizer()
  {
    return this->m_RegistrationOptimizer;
  }

  void SetTransformationModel( TransformationModelPointer T )
  {
    this->m_TransformationModel = T;
  }

  void SetRegistrationOptimizer( RegistrationOptimizerPointer T )
  {
    this->m_RegistrationOptimizer = T;
  }

  void InitializeTransformAndOptimizer();

  void RunRegistration();

protected:
  PICSLAdvancedNormalizationToolKit();
  virtual ~PICSLAdvancedNormalizationToolKit()
  {
  }

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  PICSLAdvancedNormalizationToolKit( const Self & ); // purposely not implemented
  void operator=( const Self & );                    // purposely not implemented

  ImagePointer SubsampleImage( ImagePointer, RealType );
  ImagePointer PreprocessImage( ImagePointer );
  ImagePointer ReplaceProblematicPixelValues( ImagePointer, PixelType );

  void InitializeCommandLineOptions();

  void ReadImagesAndMetrics();

  typename ParserType::Pointer                             m_Parser;

  TransformationModelPointer   m_TransformationModel;
  RegistrationOptimizerPointer m_RegistrationOptimizer;

  SimilarityMetricListType m_SimilarityMetrics;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPICSLAdvancedNormalizationToolKit.hxx"
#endif

#endif
