/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsGrubbsRosnerListSampleFilter_h
#define __antsGrubbsRosnerListSampleFilter_h

#include "antsListSampleToListSampleFilter.h"

#include <vector>

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class GrubbsRosnerListSampleFilter
 * \brief Base class of filters intended to generate scalar samples from
 * intensity samples.
 *
 */

template <typename TScalarListSample>
class GrubbsRosnerListSampleFilter final : public ListSampleToListSampleFilter<TScalarListSample, TScalarListSample>
{
public:
  /**
   * Standard class typedefs.
   */
  typedef GrubbsRosnerListSampleFilter                                       Self;
  typedef ListSampleToListSampleFilter<TScalarListSample, TScalarListSample> Superclass;
  typedef SmartPointer<Self>                                                 Pointer;
  typedef SmartPointer<const Self>                                           ConstPointer;

  /**
   * Standard macros
   */
  itkOverrideGetNameOfClassMacro(GrubbsRosnerListSampleFilter);

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /**
   * Conveneient typedefs
   */
  typedef double                                                RealType;
  typedef TScalarListSample                                     ScalarListSampleType;
  typedef typename ScalarListSampleType ::MeasurementVectorType MeasurementVectorType;
  typedef typename ScalarListSampleType ::InstanceIdentifier    InstanceIdentifierType;
  typedef std::vector<InstanceIdentifierType>                   InstanceIdentifierContainerType;

  enum OutlierHandlingType
  {
    None,
    Trim,
    Winsorize
  };

  itkSetMacro(OutlierHandling, OutlierHandlingType);
  itkGetConstMacro(OutlierHandling, OutlierHandlingType);

  itkSetMacro(WinsorizingLevel, RealType);
  itkGetConstMacro(WinsorizingLevel, RealType);

  itkSetMacro(SignificanceLevel, RealType);
  itkGetConstMacro(SignificanceLevel, RealType);

  InstanceIdentifierContainerType
  GetOutlierInstanceIdentifiers()
  {
    return this->m_OutlierInstanceIdentifiers;
  }

  //   itkGetConstMacro( Outliers, InstanceIdentifierContainerType );
protected:
  GrubbsRosnerListSampleFilter();
  ~GrubbsRosnerListSampleFilter() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  GrubbsRosnerListSampleFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  InstanceIdentifierType FindMaximumNonOutlierDeviationValue(RealType, RealType);
  bool
  IsMeasurementAnOutlier(RealType, RealType, RealType, unsigned long);

  OutlierHandlingType             m_OutlierHandling;
  RealType                        m_WinsorizingLevel;
  InstanceIdentifierContainerType m_OutlierInstanceIdentifiers;
  RealType                        m_SignificanceLevel;
}; // end of class
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsGrubbsRosnerListSampleFilter.hxx"
#endif

#endif
