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
#ifndef __antsJointHistogramParzenWindowsListSampleFunction_h
#define __antsJointHistogramParzenWindowsListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkImage.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class JointHistogramParzenWindowsListSampleFunction.h
 * \brief
 */

template <typename TListSample, typename TOutput = double, typename TCoordRep = double>
class JointHistogramParzenWindowsListSampleFunction : public ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef JointHistogramParzenWindowsListSampleFunction       Self;
  typedef ListSampleFunction<TListSample, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(JointHistogramParzenWindowsListSampleFunction);

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;
  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  typedef Image<RealType, 2>  JointHistogramImageType;
  typedef Vector<RealType, 2> ThetaPsiType;

  /** Helper functions */

  itkSetMacro(Sigma, RealType);
  itkGetConstMacro(Sigma, RealType);

  itkSetMacro(NumberOfJointHistogramBins, unsigned int);
  itkGetConstMacro(NumberOfJointHistogramBins, unsigned int);

  virtual void
  SetInputListSample(const InputListSampleType * ptr);

  virtual TOutput
  Evaluate(const InputMeasurementVectorType & measurement) const;

protected:
  JointHistogramParzenWindowsListSampleFunction();
  virtual ~JointHistogramParzenWindowsListSampleFunction();
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  void
  GenerateData();

  void
  IncrementJointHistogram(RealType e1, RealType e2, unsigned int which);

private:
  // purposely not implemented
  JointHistogramParzenWindowsListSampleFunction(const Self &);
  void
  operator=(const Self &);

  unsigned int                                           m_NumberOfJointHistogramBins;
  RealType                                               m_Sigma;
  bool                                                   m_UseNNforJointHistIncrements;
  std::vector<typename JointHistogramImageType::Pointer> m_JointHistogramImages;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsJointHistogramParzenWindowsListSampleFunction.hxx"
#endif

#endif
