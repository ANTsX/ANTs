/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsBoxPlotQuantileListSampleFilter.h,v $
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
#ifndef __antsBoxPlotQuantileListSampleFilter_h
#define __antsBoxPlotQuantileListSampleFilter_h

#include "antsListSampleToListSampleFilter.h"

#include <vector>

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class BoxPlotQuantileListSampleFilter
 * \brief Base class of filters intended to generate scalar samples from
 * intensity samples.
 *
 */

template <class TScalarListSample>
class BoxPlotQuantileListSampleFilter
  : public       ListSampleToListSampleFilter<TScalarListSample, TScalarListSample>
{
public:
  /**
   * Standard class typedefs.
   */
  typedef BoxPlotQuantileListSampleFilter Self;
  typedef ListSampleToListSampleFilter
    <TScalarListSample, TScalarListSample>                    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /**
   * Standard macros
   */
  itkTypeMacro( BoxPlotQuantileListSampleFilter,
                ListSampleToScalarListSampleFilter );

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro( Self );

  /**
   * Conveneient typedefs
   */
  typedef double            RealType;
  typedef TScalarListSample ScalarListSampleType;
  typedef typename ScalarListSampleType
    ::MeasurementVectorType                           MeasurementVectorType;
  typedef typename ScalarListSampleType
    ::InstanceIdentifier                              InstanceIdentifierType;
  typedef std::vector<InstanceIdentifierType> InstanceIdentifierContainerType;

  enum OutlierHandlingType { None, Trim, Winsorize };

  itkSetMacro( OutlierHandling, OutlierHandlingType );
  itkGetConstMacro( OutlierHandling, OutlierHandlingType );

  itkSetMacro( WhiskerScalingFactor, RealType );
  itkGetConstMacro( WhiskerScalingFactor, RealType );

  itkSetClampMacro( UpperPercentile, RealType, 0, 1 );
  itkGetConstMacro( UpperPercentile, RealType );

  itkSetClampMacro( LowerPercentile, RealType, 0, 1 );
  itkGetConstMacro( LowerPercentile, RealType );

  InstanceIdentifierContainerType GetOutlierInstanceIdentifiers()
  {
    return this->m_OutlierInstanceIdentifiers;
  }

//   itkGetConstMacro( Outliers, InstanceIdentifierContainerType );
protected:
  BoxPlotQuantileListSampleFilter();
  virtual ~BoxPlotQuantileListSampleFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  virtual void GenerateData();

private:
  BoxPlotQuantileListSampleFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                  // purposely not implemented

  InstanceIdentifierType FindMaximumNonOutlierDeviationValue( RealType, RealType );
  bool                   IsMeasurementAnOutlier( RealType, RealType, RealType, unsigned long );

  OutlierHandlingType             m_OutlierHandling;
  InstanceIdentifierContainerType m_OutlierInstanceIdentifiers;
  RealType                        m_WhiskerScalingFactor;
  RealType                        m_LowerPercentile;
  RealType                        m_UpperPercentile;
};    // end of class
} // end of namespace Statistics
} // end of namespace ants4443
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsBoxPlotQuantileListSampleFilter.hxx"
#endif

#endif
