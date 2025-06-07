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
#ifndef __antsGaussianListSampleFunction_h
#define __antsGaussianListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkGaussianMembershipFunction.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class GaussianListSampleFunction.h
 * \brief point set filter.
 */

template <typename TListSample, typename TOutput = double, typename TCoordRep = double>
class GaussianListSampleFunction final : public ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef GaussianListSampleFunction                          Self;
  typedef ListSampleFunction<TListSample, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(GaussianListSampleFunction);

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  typedef typename Superclass::ListSampleWeightArrayType ListSampleWeightArrayType;

  /** Gaussian typedefs */
  typedef typename itk::Statistics::GaussianMembershipFunction<InputMeasurementVectorType> GaussianType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  /** Helper functions */

  void
  SetInputListSample(const InputListSampleType * ptr) override;

  TOutput
  Evaluate(const InputMeasurementVectorType & measurement) const override;

protected:
  GaussianListSampleFunction();
  ~GaussianListSampleFunction() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData();

private:
  // purposely not implemented
  GaussianListSampleFunction(const Self &);
  void
  operator=(const Self &);

  typename GaussianType::Pointer m_Gaussian;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsGaussianListSampleFunction.hxx"
#endif

#endif
