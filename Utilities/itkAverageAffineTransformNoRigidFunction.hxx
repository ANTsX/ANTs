/*=========================================================================

 Program:   Advanced Normalization Tools
 Module:    $RCSfile: itkWarpImageMultiTransformFilter.hxx,v $
 Language:  C++
 Date:      $Date: 2009/01/08 21:36:48 $
 Version:   $Revision: 1.18 $

 Copyright (c) ConsortiumOfANTS. All rights reserved.
 See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkAverageAffineTransformNoRigidFunction_hxx
#define __itkAverageAffineTransformNoRigidFunction_hxx
#include "itkAverageAffineTransformNoRigidFunction.h"

#include "itkNumericTraits.h"
#include <limits>

namespace itk
{
/**
 * Default constructor.
 */
template <class TTransform>
AverageAffineTransformNoRigidFunction<TTransform>::AverageAffineTransformNoRigidFunction()
{
}

template <class TTransform>
void AverageAffineTransformNoRigidFunction<TTransform>::PrintTransformList()
{
  std::cout << "transform list: " << std::endl;

  typename TransformListType::iterator it = (m_TransformList.begin() );
  for( int ii = 0; it != m_TransformList.end(); it++, ii++ )
    {
    std::cout << '[' << ii << ":" << it->weight << "]:" << it->aff
                     << std::endl;
    }
}

// /**
// * Standard PrintSelf method.
// */
// template <class TInputImage,class TOutputImage,class TDisplacementField, class TTransform>
// void
// WarpImageMultiTransformFilter<TInputImage,TOutputImage,TDisplacementField, TTransform>
// ::PrintSelf(std::ostream& os, Indent indent) const
// {
//
//    Superclass::PrintSelf(os, indent);
// }

template <class TTransform>
void AverageAffineTransformNoRigidFunction<TTransform>::PushBackAffineTransform(
  const GenericAffineTransformType* t, double weight)
{
  if( t )
    {
    SingleTransformItemType item;
    item.aff = const_cast<GenericAffineTransformType *>(t);
    item.weight = weight;
    m_TransformList.push_back(SingleTransformItemType(item) );
    }
}

template <class TTransform>
void AverageAffineTransformNoRigidFunction<TTransform>::AverageMultipleAffineTransform(
  const PointType & reference_center,
  GenericAffineTransformPointerType & affine_output)
{
//    std::cout << "test " ;
//    TransformTypePointer affine_output = TransformType::New();

  affine_output->SetIdentity();
  affine_output->SetCenter(reference_center);

  unsigned int number_of_affine = m_TransformList.size();

  number_of_affine--;

//    std::cout << affine_output;

  typename TransformListType::iterator it = m_TransformList.begin();

//    typename InternalAffineTransformType::InputPointType center_ref = m_ReferenceCenter;
  typename InternalAffineTransformType::Pointer average_iaff =
    InternalAffineTransformType::New();

  typename InternalAffineTransformType::ParametersType average_parameters =
    average_iaff->GetParameters();
  for( ; it != m_TransformList.end(); it++ )
    {
    SingleInternalTransformItemType internal_item;
    internal_item.aff = InternalAffineTransformType::New();
    ConvertGenericAffineToInternalAffineByFixingCenter(it->aff,
                                                       internal_item.aff, reference_center);
    internal_item.weight = it->weight;
    m_InternalTransformList.push_back(internal_item);

    std::cout << "internal_transform: " << internal_item.aff << std::endl;
    }

  HelperType::ComputeAverageScaleParameters(m_InternalTransformList,
                                            average_parameters);
  HelperType::ComputeAverageShearingParameters(m_InternalTransformList,
                                               average_parameters);
  // HelperType::ComputeAverageRotationParameters(m_InternalTransformList,
  //                                             average_parameters);
  // HelperType::ComputeAverageTranslationParameters(m_InternalTransformList,
  //                                                average_parameters);

  average_iaff->SetParameters(average_parameters);
  average_iaff->SetCenter(reference_center);

  std::cout << "average_iaff" << average_iaff << std::endl;

  ConvertInternalAffineToGenericAffine(average_iaff, affine_output);

  std::cout << "affine_output" << affine_output << std::endl;
  return;
}

template <class TTransform>
void AverageAffineTransformNoRigidFunction<TTransform>::ConvertGenericAffineToInternalAffineByFixingCenter(
  GenericAffineTransformPointerType & aff,
  InternalAffineTransformPointerType & iaff, const PointType & center)
{
  iaff->SetCenter(center);
  iaff->SetMatrix(aff->GetMatrix() );
  iaff->SetTranslation(aff->GetTranslation() );

  return;
}

template <class TTransform>
void AverageAffineTransformNoRigidFunction<TTransform>::ConvertInternalAffineToGenericAffine(
  InternalAffineTransformPointerType & iaff,
  GenericAffineTransformPointerType & aff)
{
  aff->SetCenter(iaff->GetCenter() );
  aff->SetTranslation(iaff->GetTranslation() );
  aff->SetMatrix(iaff->GetMatrix() );

  return;
}

namespace AverageAffineTransformNoRigidFunctionHelperNameSpace
{
template <class TAffine>
void HelperCommonType<TAffine>::ComputeAveragePartialParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters, unsigned int istart,
  unsigned int iend)
{
  double w = 0.0;

  // initialize partial parameters to zero
  for( unsigned int k = istart; k <= iend; k++ )
    {
    average_parameters[k] = 0.0;
    }

  typename InternalTransformListType::iterator it = transform_list.begin();
  unsigned int cnt = 0;
  for( ; it != transform_list.end(); it++ )
    {
    ParametersType current_parameters = it->aff->GetParameters();
    w += it->weight;

    std::cout << "[" << cnt++ << "]:" << it->weight << "\t";
    for( unsigned int k = istart; k <= iend; k++ )
      {
      average_parameters[k] += it->weight * current_parameters[k];

      std::cout << current_parameters[k] << " ";
      }

    std::cout << std::endl;
    }

  if( w <= 0.0 )
    {
    std::cout << "Total weight smaller than 0!!!" << std::endl;
    std::exception();
    }

  // normalize by weight
  std::cout << "sum:w=" << w <<  "\t";
  for( unsigned int k = istart; k <= iend; k++ )
    {
    std::cout << average_parameters[k] << " ";
    }
  std::cout << std::endl;

  // normalize by weight
  std::cout << "average" << "\t";
  for( unsigned int k = istart; k <= iend; k++ )
    {
    average_parameters[k] /= w;
    std::cout << average_parameters[k] << " ";
    }

  std::cout << std::endl;
  return;
}

void HelperType<Dispatcher<2> >::ComputeAverageScaleParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 1;
  unsigned int iend = 2;

  std::cout << "average 2D scale parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<2> >::ComputeAverageShearingParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 3;
  unsigned int iend = 3;

  std::cout << "average 2D shearing parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<2> >::ComputeAverageRotationParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 0;
  unsigned int iend = 0;

  std::cout << "average 2D rotation parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<2> >::ComputeAverageTranslationParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 6;
  unsigned int iend = 7;

  std::cout << "average 2D translation parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<3> >::ComputeAverageScaleParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 4;
  unsigned int iend = 6;

  std::cout << "average 3D scale parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<3> >::ComputeAverageShearingParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 7;
  unsigned int iend = 9;

  std::cout << "average 3D shearing parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}

void HelperType<Dispatcher<3> >::ComputeAverageRotationParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 0;
  unsigned int iend = 3;

  std::cout << "average 3D rotation parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);

  // extra normalization for quaternion

  double quat_mag = 0.0;
  for( unsigned int j = istart; j <= iend; j++ )
    {
    quat_mag += average_parameters[j] * average_parameters[j];
    }
  quat_mag = sqrt(quat_mag);
  for( unsigned int j = 0; j < 4; j++ )
    {
    average_parameters[j] /= quat_mag;
    }
}

void HelperType<Dispatcher<3> >::ComputeAverageTranslationParameters(
  InternalTransformListType & transform_list,
  ParametersType & average_parameters)
{
  unsigned int istart = 10;
  unsigned int iend = 12;

  std::cout << "average 3D translation parameter " << std::endl;

  HelperCommonType<InternalAffineTransformType>::ComputeAveragePartialParameters(
    transform_list, average_parameters, istart, iend);
}
} // end namespace AverageAffineTransformNoRigidFunctionHelperNameSpace
} // end namespace itk

#endif  // __itkAverageAffineTransformNoRigidFunction_hxx
