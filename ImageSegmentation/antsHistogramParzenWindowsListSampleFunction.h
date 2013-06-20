/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsHistogramParzenWindowsListSampleFunction.h,v $
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
#ifndef __antsHistogramParzenWindowsListSampleFunction_h
#define __antsHistogramParzenWindowsListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkImage.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class HistogramParzenWindowsListSampleFunction.h
 * \brief point set filter.
 */

template <class TListSample, class TOutput = double, class TCoordRep = double>
class HistogramParzenWindowsListSampleFunction
  : public       ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef HistogramParzenWindowsListSampleFunction Self;
  typedef ListSampleFunction
    <TListSample, TOutput, TCoordRep>                      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( HistogramParzenWindowsListSampleFunction, ListSampleFunction );

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  typedef Image<RealType, 1> HistogramImageType;
  typedef BSplineInterpolateImageFunction<HistogramImageType>
    InterpolatorType;
  typedef LinearInterpolateImageFunction<HistogramImageType>
    LInterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  /** Helper functions */

  itkSetMacro( Sigma, RealType );
  itkGetConstMacro( Sigma, RealType );

  itkSetMacro( NumberOfHistogramBins, unsigned int );
  itkGetConstMacro( NumberOfHistogramBins, unsigned int );

  virtual void SetInputListSample( const InputListSampleType * ptr );

  virtual TOutput Evaluate( const InputMeasurementVectorType& measurement ) const;

protected:
  HistogramParzenWindowsListSampleFunction();
  virtual ~HistogramParzenWindowsListSampleFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  // purposely not implemented
  HistogramParzenWindowsListSampleFunction( const Self & );
  void operator=( const Self & );

  unsigned int                                      m_NumberOfHistogramBins;
  RealType                                          m_Sigma;
  InterpolatorPointer                               m_Interpolator;
  std::vector<typename HistogramImageType::Pointer> m_HistogramImages;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsHistogramParzenWindowsListSampleFunction.hxx"
#endif

#endif
