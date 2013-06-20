/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsManifoldParzenWindowsListSampleFunction.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsManifoldParzenWindowsListSampleFunction_h
#define __antsManifoldParzenWindowsListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkGaussianMembershipFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include <vector>

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class ManifoldParzenWindowsListSampleFunction.h
 * \brief point set filter.
 */

template <class TListSample, class TOutput = double, class TCoordRep = double>
class ManifoldParzenWindowsListSampleFunction
  : public       ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef ManifoldParzenWindowsListSampleFunction Self;
  typedef ListSampleFunction
    <TListSample, TOutput, TCoordRep>                      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ManifoldParzenWindowsListSampleFunction, ListSampleFunction );

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Kd tree typedefs */
  typedef typename itk::Statistics::
    WeightedCentroidKdTreeGenerator<InputListSampleType>    TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType KdTreeType;
  typedef typename KdTreeType
    ::InstanceIdentifierVectorType                          NeighborhoodIdentifierType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  typedef typename itk::Statistics::GaussianMembershipFunction
    <InputMeasurementVectorType>                         GaussianType;
  typedef std::vector<typename GaussianType::Pointer> GaussianContainerType;
  typedef typename GaussianType::CovarianceMatrixType CovarianceMatrixType;

  /** Helper functions */

  itkSetMacro( EvaluationKNeighborhood, unsigned int );
  itkGetConstMacro( EvaluationKNeighborhood, unsigned int );

  itkSetMacro( RegularizationSigma, RealType );
  itkGetConstMacro( RegularizationSigma, RealType );

  itkSetMacro( CovarianceKNeighborhood, unsigned int );
  itkGetConstMacro( CovarianceKNeighborhood, unsigned int );

  itkSetMacro( KernelSigma, RealType );
  itkGetConstMacro( KernelSigma, RealType );

  virtual void SetInputListSample( const InputListSampleType * ptr );

  virtual TOutput Evaluate( const InputMeasurementVectorType& measurement ) const;

protected:
  ManifoldParzenWindowsListSampleFunction();
  virtual ~ManifoldParzenWindowsListSampleFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  // purposely not implemented
  ManifoldParzenWindowsListSampleFunction( const Self & );
  void operator=( const Self & );

  unsigned int m_CovarianceKNeighborhood;
  unsigned int m_EvaluationKNeighborhood;
  RealType     m_RegularizationSigma;
  RealType     m_KernelSigma;
  RealType     m_NormalizationFactor;

  typename TreeGeneratorType::Pointer           m_KdTreeGenerator;

  GaussianContainerType m_Gaussians;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsManifoldParzenWindowsListSampleFunction.hxx"
#endif

#endif
