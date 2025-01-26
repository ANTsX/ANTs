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

template <typename TListSample, typename TOutput = double, typename TCoordRep = double>
class PartialVolumeGaussianListSampleFunction final : public ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef PartialVolumeGaussianListSampleFunction             Self;
  typedef ListSampleFunction<TListSample, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PartialVolumeGaussianListSampleFunction);

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;

  typedef typename Superclass::ListSampleWeightArrayType ListSampleWeightArrayType;

  /** Gaussian typedefs */
  typedef typename itk::Statistics::GaussianMembershipFunction<InputMeasurementVectorType> GaussianType;
  typedef typename GaussianType::MeanVectorType                                            MeanType;
  typedef typename GaussianType::CovarianceMatrixType                                      CovarianceType;

  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  void
  SetIndexedInputListSample(unsigned int d, const InputListSampleType * ptr) override;

  TOutput
  Evaluate(const InputMeasurementVectorType & measurement) const override;

protected:
  PartialVolumeGaussianListSampleFunction();
  ~PartialVolumeGaussianListSampleFunction() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData();

private:
  // purposely not implemented
  PartialVolumeGaussianListSampleFunction(const Self &);
  void
  operator=(const Self &);

  void
  CalculateGaussianParametersFromListSample(const InputListSampleType *, const ListSampleWeightArrayType *, MeanType &);

  void
  CalculateGaussianParameters();

  MeanType m_Mean[2];
  bool     m_IsCalculated[2];

  typename GaussianType::Pointer m_Gaussian;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsPartialVolumeGaussianListSampleFunction.hxx"
#endif

#endif
