/*=========================================================================

 Program:   Advanced Normalization Tools

 Copyright (c) ConsortiumOfANTS. All rights reserved.
 See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkAverageAffineTransformFunction_hxx
#define __itkAverageAffineTransformFunction_hxx

#include "itkNumericTraits.h"
#include <limits>

namespace itk
{
/**
 * Default constructor.
 */
template <typename TTransform>
AverageAffineTransformFunction<TTransform>::AverageAffineTransformFunction() = default;

template <typename TTransform>
void
AverageAffineTransformFunction<TTransform>::PrintTransformList()
{
  std::cout << "transform list: " << std::endl;

  typename TransformListType::iterator it = (m_TransformList.begin());
  for (int ii = 0; it != m_TransformList.end(); it++, ii++)
  {
    std::cout << '[' << ii << ":" << it->weight << "]:" << it->aff << std::endl;
  }
}

// /**
// * Standard PrintSelf method.
// */
// template <typename TInputImage,typename TOutputImage,typename TDisplacementField, typename TTransform>
// void
// WarpImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::PrintSelf(std::ostream& os, Indent indent) const
// {
//
//    Superclass::PrintSelf(os, indent);
// }

template <typename TTransform>
void
AverageAffineTransformFunction<TTransform>::PushBackAffineTransform(const GenericAffineTransformType * t, double weight)
{
  if (t)
  {
    SingleTransformItemType item;
    item.aff = const_cast<GenericAffineTransformType *>(t);
    item.weight = weight;
    m_TransformList.push_back(SingleTransformItemType(item));
  }
}

template <typename TTransform>
void
AverageAffineTransformFunction<TTransform>::AverageMultipleAffineTransform(
  const PointType &                   reference_center,
  GenericAffineTransformPointerType & affine_output)
{
  // TransformTypePointer affine_output = TransformType::New();

  affine_output->SetIdentity();
  affine_output->SetCenter(reference_center);

  typename TransformListType::iterator it = m_TransformList.begin();

  typename InternalAffineTransformType::Pointer average_iaff = InternalAffineTransformType::New();
  if (verbose)
  {
    average_iaff->DebugOn();
  }

  typename InternalAffineTransformType::ParametersType average_parameters = average_iaff->GetParameters();
  for (; it != m_TransformList.end(); it++)
  {
    SingleInternalTransformItemType internal_item;
    internal_item.aff = InternalAffineTransformType::New();
    if (verbose)
    {
      internal_item.aff->DebugOn();
    }
    ConvertGenericAffineToInternalAffineByFixingCenter(it->aff, internal_item.aff, reference_center);
    internal_item.weight = it->weight;
    m_InternalTransformList.push_back(internal_item);

    if (verbose)
    {
      std::cout << "internal_transform: " << internal_item.aff << std::endl;
    }
  }

  HelperType::ComputeAverageScaleParameters(m_InternalTransformList, average_parameters, verbose);
  HelperType::ComputeAverageShearingParameters(m_InternalTransformList, average_parameters, verbose);
  if (useRigid)
  {
    HelperType::ComputeAverageRotationParameters(m_InternalTransformList, average_parameters, verbose);
    HelperType::ComputeAverageTranslationParameters(m_InternalTransformList, average_parameters, verbose);
  }

  average_iaff->SetParameters(average_parameters);
  average_iaff->SetCenter(reference_center);

  if (verbose)
  {
    std::cout << "average_iaff" << average_iaff << std::endl;
  }

  ConvertInternalAffineToGenericAffine(average_iaff, affine_output);

  if (verbose)
  {
    std::cout << "affine_output" << affine_output << std::endl;
  }
  return;
}

template <typename TTransform>
void
AverageAffineTransformFunction<TTransform>::ConvertGenericAffineToInternalAffineByFixingCenter(
  GenericAffineTransformPointerType &  aff,
  InternalAffineTransformPointerType & iaff,
  const PointType &                    center)
{
  iaff->SetCenter(center);
  iaff->SetMatrix(aff->GetMatrix());
  iaff->SetTranslation(aff->GetTranslation());

  return;
}

template <typename TTransform>
void
AverageAffineTransformFunction<TTransform>::ConvertInternalAffineToGenericAffine(
  InternalAffineTransformPointerType & iaff,
  GenericAffineTransformPointerType &  aff)
{
  aff->SetCenter(iaff->GetCenter());
  aff->SetTranslation(iaff->GetTranslation());
  aff->SetMatrix(iaff->GetMatrix());

  return;
}

namespace AverageAffineTransformFunctionHelperNameSpace
{
template <typename TAffine>
void
HelperCommonType<TAffine>::ComputeAveragePartialParameters(InternalTransformListType & transform_list,
                                                           ParametersType &            average_parameters,
                                                           unsigned int                istart,
                                                           unsigned int                iend, bool verbose)
{
  double w = 0.0;

  // initialize partial parameters to zero
  for (unsigned int k = istart; k <= iend; k++)
  {
    average_parameters[k] = 0.0;
  }

  typename InternalTransformListType::iterator it = transform_list.begin();
  unsigned int                                 cnt = 0;
  for (; it != transform_list.end(); it++)
  {
    ParametersType current_parameters = it->aff->GetParameters();
    w += it->weight;
    ++cnt;
    for (unsigned int k = istart; k <= iend; k++)
    {
      average_parameters[k] += it->weight * current_parameters[k];
    }

    if (verbose)
    {
      std::cout << "[" << cnt << "]:" << it->weight << "\t";
      for (unsigned int k = istart; k <= iend; k++)
      {
        std::cout << current_parameters[k] << " ";
      }
      std::cout << std::endl;
    }
  }

  if (w <= 0.0)
  {
    if (verbose)
    {
      std::cout << "Total weight smaller than 0!!!" << std::endl;
    }
    // throw std::runtime_error("Total weight smaller than 0");
  }

  if (verbose)
  {
    std::cout << "sum:w=" << w << "\t";
    for (unsigned int k = istart; k <= iend; k++)
    {
      std::cout << average_parameters[k] << " ";
    }
    std::cout << std::endl;
  }

  // normalize by weight
  for (unsigned int k = istart; k <= iend; k++)
  {
    average_parameters[k] /= w;
  }

  if (verbose)
  {
    std::cout << "average\t";
    for (unsigned int k = istart; k <= iend; k++)
    {
      std::cout << average_parameters[k] << " ";
    }
    std::cout << std::endl;
  }

  return;
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<2>, TParametersValueType>::ComputeAverageScaleParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 1;
  unsigned int iend = 2;

  if (verbose)
  {
    std::cout << "average 2D scale parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<2>, TParametersValueType>::ComputeAverageShearingParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 3;
  unsigned int iend = 3;

  if (verbose)
  {
    std::cout << "average 2D shearing parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<2>, TParametersValueType>::ComputeAverageRotationParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 0;
  unsigned int iend = 0;

  if (verbose)
  {
    std::cout << "average 2D rotation parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<2>, TParametersValueType>::ComputeAverageTranslationParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 6;
  unsigned int iend = 7;

  if (verbose)
  {
    std::cout << "average 2D translation parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<3>, TParametersValueType>::ComputeAverageScaleParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 4;
  unsigned int iend = 6;

  if (verbose)
  {
    std::cout << "average 3D scale parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<3>, TParametersValueType>::ComputeAverageShearingParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters, bool verbose)
{
  unsigned int istart = 7;
  unsigned int iend = 9;

  if (verbose)
  {
    std::cout << "average 3D shearing parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<3>, TParametersValueType>::ComputeAverageRotationParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters,
        bool verbose)
{
  unsigned int istart = 0;
  unsigned int iend = 3;

  if (verbose)
  {
    std::cout << "average 3D rotation parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);

  // extra normalization for quaternion

  double quat_mag = 0.0;
  for (unsigned int j = istart; j <= iend; j++)
  {
    quat_mag += average_parameters[j] * average_parameters[j];
  }
  quat_mag = sqrt(quat_mag);
  for (unsigned int j = 0; j < 4; j++)
  {
    average_parameters[j] /= quat_mag;
  }
}

template <typename TParametersValueType>
void
HelperType<Dispatcher<3>, TParametersValueType>::ComputeAverageTranslationParameters(
  InternalTransformListType & transform_list,
  ParametersType &            average_parameters,
        bool verbose)
{
  unsigned int istart = 10;
  unsigned int iend = 12;

  if (verbose)
  {
    std::cout << "average 3D translation parameter " << std::endl;
  }

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend, verbose);
}
} // end namespace AverageAffineTransformFunctionHelperNameSpace
} // end namespace itk

#endif // __itkAverageAffineTransformFunction_hxx
