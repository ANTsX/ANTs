/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsLogEuclideanGaussianListSampleFunction.h,v $
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
#ifndef __antsLogEuclideanGaussianListSampleFunction_h
#define __antsLogEuclideanGaussianListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkVariableSizeMatrix.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class LogEuclideanGaussianListSampleFunction.h
 * \brief
 */

template <class TListSample, class TOutput = double, class TCoordRep = double>
class LogEuclideanGaussianListSampleFunction
  : public       ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef LogEuclideanGaussianListSampleFunction Self;
  typedef ListSampleFunction
    <TListSample, TOutput, TCoordRep>                      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( LogEuclideanGaussianListSampleFunction, ListSampleFunction );

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  /** Other typedef */
  typedef TOutput                      RealType;
  typedef TOutput                      OutputType;
  typedef VariableSizeMatrix<RealType> TensorType;

  /** Helper functions */

  virtual void SetInputListSample( const InputListSampleType * ptr );

  virtual TOutput Evaluate( const InputMeasurementVectorType& measurement ) const;

protected:
  LogEuclideanGaussianListSampleFunction();
  virtual ~LogEuclideanGaussianListSampleFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

  TensorType LogTensorTransform( const TensorType & ) const;

  TensorType ExpTensorTransform( const TensorType & ) const;

  RealType CalculateTensorDistance( const TensorType &, const TensorType & ) const;

  TensorType m_MeanTensor;
  RealType   m_Dispersion;
private:
  // purposely not implemented
  LogEuclideanGaussianListSampleFunction( const Self & );
  void operator=( const Self & );
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsLogEuclideanGaussianListSampleFunction.hxx"
#endif

#endif
