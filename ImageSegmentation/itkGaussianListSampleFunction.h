/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianListSampleFunction.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianListSampleFunction_h
#define __itkGaussianListSampleFunction_h

#include "itkListSampleFunction.h"

#include "itkGaussianMembershipFunction.h"

namespace itk
{
namespace Statistics
{
/** \class GaussianListSampleFunction.h
 * \brief point set filter.
 */

template <class TListSample, class TOutput = double, class TCoordRep = double>
class ITK_EXPORT GaussianListSampleFunction
  : public       ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef GaussianListSampleFunction Self;
  typedef ListSampleFunction
    <TListSample, TOutput, TCoordRep>                      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GaussianListSampleFunction, ListSampleFunction );

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Gaussian typedefs */
  typedef typename Statistics::GaussianMembershipFunction
    <InputMeasurementVectorType>                          GaussianType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  /** Helper functions */

  virtual void SetInputListSample( const InputListSampleType * ptr );

  virtual TOutput Evaluate( const InputMeasurementVectorType& measurement ) const;

protected:
  GaussianListSampleFunction();
  virtual ~GaussianListSampleFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  // purposely not implemented
  GaussianListSampleFunction( const Self & );
  void operator=( const Self & );

  typename GaussianType::Pointer                                  m_Gaussian;
};
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianListSampleFunction.txx"
#endif

#endif
