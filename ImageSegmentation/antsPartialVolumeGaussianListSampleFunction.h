/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsPartialVolumeGaussianListSampleFunction.h,v $
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
#ifndef __antsPartialVolumeGaussianListSampleFunction_h
#define __antsPartialVolumeGaussianListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkGaussianMembershipFunction.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class PartialVolumeGaussianListSampleFunction.h
 * \brief point set filter.
 */

template <class TListSample, class TOutput = double, class TCoordRep = double>
class PartialVolumeGaussianListSampleFunction
  : public       ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef PartialVolumeGaussianListSampleFunction Self;
  typedef ListSampleFunction
    <TListSample, TOutput, TCoordRep>                      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( PartialVolumeGaussianListSampleFunction,
                ListSampleFunction );

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  typedef typename Superclass::ListSampleWeightArrayType ListSampleWeightArrayType;

  /** Gaussian typedefs */
  typedef typename itk::Statistics::GaussianMembershipFunction
    <InputMeasurementVectorType>                            GaussianType;
  typedef typename GaussianType::MeanVectorType       MeanType;
  typedef typename GaussianType::CovarianceMatrixType CovarianceType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  virtual void SetIndexedInputListSample(unsigned int d, const InputListSampleType * ptr );

  virtual TOutput Evaluate( const InputMeasurementVectorType& measurement ) const;

protected:
  PartialVolumeGaussianListSampleFunction();
  virtual ~PartialVolumeGaussianListSampleFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  // purposely not implemented
  PartialVolumeGaussianListSampleFunction( const Self & );
  void operator=( const Self & );

  void CalculateGaussianParametersFromListSample( const InputListSampleType *, const ListSampleWeightArrayType *,
                                                  MeanType & );

  void CalculateGaussianParameters();

  MeanType m_Mean[2];
  bool     m_IsCalculated[2];

  typename GaussianType::Pointer                            m_Gaussian;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsPartialVolumeGaussianListSampleFunction.hxx"
#endif

#endif
