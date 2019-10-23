/*=========================================================================

 Program:   Advanced Normalization Tools

 Copyright (c) ConsortiumOfANTS. All rights reserved.
 See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkAverageAffineTransformNoRigidFunction_h
#define __itkAverageAffineTransformNoRigidFunction_h

#include "itkMacro.h"
#include "itkConceptChecking.h"

#include "itkANTSCenteredAffine2DTransform.h"
#include "itkANTSAffine3DTransform.h"

#include <list>

namespace itk
{
/*
 */

namespace AverageAffineTransformNoRigidFunctionHelperNameSpace
{
template <int D>
struct Dispatcher
  {
  };

template <typename TAffine>
struct HelperCommonType
  {
  typedef TAffine InternalAffineTransformType;

  typedef typename InternalAffineTransformType::Pointer InternalAffineTransformPointerType;
  typedef struct
    {
    InternalAffineTransformPointerType aff;
    double weight;
    } SingleInternalTransformItemType;

  typedef std::list<SingleInternalTransformItemType>           InternalTransformListType;
  typedef typename InternalAffineTransformType::ParametersType ParametersType;

  static void ComputeAveragePartialParameters(InternalTransformListType & transform_list,
                                              ParametersType & average_parameters, unsigned int iStart,
                                              unsigned int iEnd);
  };

template <typename T>
class HelperType;

// {
// // purposely not include any types
// };

// explicit specialization for 2D affine transform
template <>
class HelperType<Dispatcher<2> >
{
public:
  typedef ANTSCenteredAffine2DTransform<double> InternalAffineTransformType;

  typedef HelperCommonType<InternalAffineTransformType>::InternalAffineTransformPointerType
    InternalAffineTransformPointerType;
  typedef HelperCommonType<InternalAffineTransformType>::SingleInternalTransformItemType
    SingleInternalTransformItemType;
  typedef HelperCommonType<InternalAffineTransformType>::InternalTransformListType InternalTransformListType;
  typedef HelperCommonType<InternalAffineTransformType>::ParametersType            ParametersType;

  static void ComputeAverageScaleParameters(InternalTransformListType & transform_list,
                                            ParametersType & average_parameters);

  static void ComputeAverageShearingParameters(InternalTransformListType & transform_list,
                                               ParametersType & average_parameters);

  static void ComputeAverageRotationParameters(InternalTransformListType & transform_list,
                                               ParametersType & average_parameters);

  static void ComputeAverageTranslationParameters(InternalTransformListType & transform_list,
                                                  ParametersType & average_parameters);
};

// explicit specialization for 3D affine transform
template <>
class HelperType<Dispatcher<3> >
{
public:
  typedef ANTSAffine3DTransform<double> InternalAffineTransformType;

  typedef HelperCommonType<InternalAffineTransformType>::InternalAffineTransformPointerType
    InternalAffineTransformPointerType;
  typedef HelperCommonType<InternalAffineTransformType>::SingleInternalTransformItemType
    SingleInternalTransformItemType;
  typedef HelperCommonType<InternalAffineTransformType>::InternalTransformListType InternalTransformListType;
  typedef HelperCommonType<InternalAffineTransformType>::ParametersType            ParametersType;

  static void ComputeAverageScaleParameters(InternalTransformListType & transform_list,
                                            ParametersType & average_parameters);

  static void ComputeAverageShearingParameters(InternalTransformListType & transform_list,
                                               ParametersType & average_parameters);

  static void ComputeAverageRotationParameters(InternalTransformListType & transform_list,
                                               ParametersType & average_parameters);

  static void ComputeAverageTranslationParameters(InternalTransformListType & transform_list,
                                                  ParametersType & average_parameters);
};
}

template <typename TTransform>
class AverageAffineTransformNoRigidFunction
{
public:

  typedef TTransform GenericAffineTransformType;

  static constexpr unsigned int InputSpaceDimension = GenericAffineTransformType::InputSpaceDimension;
  static constexpr unsigned int OutputSpaceDimension = GenericAffineTransformType::OutputSpaceDimension;
  static constexpr unsigned int SpaceDimension = InputSpaceDimension;

  typedef double                                    InternalScalarType;
  typedef Point<InternalScalarType, SpaceDimension> PointType;

  AverageAffineTransformNoRigidFunction();
  ~AverageAffineTransformNoRigidFunction() = default;

  // void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename GenericAffineTransformType::Pointer GenericAffineTransformPointerType;

  typedef struct
    {
    GenericAffineTransformPointerType aff;
    double weight;
    } SingleTransformItemType;

  typedef std::list<SingleTransformItemType> TransformListType;

  TransformListType & GetTransformList()
  {
    return m_TransformList;
  }

  void PrintTransformList();

  void PushBackAffineTransform(const GenericAffineTransformType* t, double weight);

  void AverageMultipleAffineTransform(const PointType & center_output,
                                      GenericAffineTransformPointerType & affine_output);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1,
                  (Concept::SameDimension<InputSpaceDimension, OutputSpaceDimension> ) );
  /** End concept checking */
#endif
protected:
  TransformListType m_TransformList;

// type declaration to include support both 2D and 3D affine transform
protected:

  typedef typename::itk::AverageAffineTransformNoRigidFunctionHelperNameSpace::HelperType<
      ::itk::AverageAffineTransformNoRigidFunctionHelperNameSpace::Dispatcher<
        SpaceDimension> > HelperType;

  typedef typename HelperType::InternalAffineTransformType        InternalAffineTransformType;
  typedef typename HelperType::InternalAffineTransformPointerType InternalAffineTransformPointerType;
  typedef typename HelperType::SingleInternalTransformItemType    SingleInternalTransformItemType;
  typedef typename HelperType::InternalTransformListType          InternalTransformListType;

  InternalTransformListType m_InternalTransformList;

  void ConvertGenericAffineToInternalAffineByFixingCenter(GenericAffineTransformPointerType & aff,
                                                          InternalAffineTransformPointerType & iaff,
                                                          const PointType & center);

  void ConvertInternalAffineToGenericAffine(InternalAffineTransformPointerType & iaff,
                                            GenericAffineTransformPointerType & aff);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAverageAffineTransformNoRigidFunction.hxx"
#endif

#endif /*__itkAverageAffineTransformNoRigidFunction_h*/
