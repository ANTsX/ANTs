/*=========================================================================

 Program:   Advanced Normalization Tools

 Copyright (c) ConsortiumOfANTS. All rights reserved.
 See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkAverageAffineTransformFunction_h
#define __itkAverageAffineTransformFunction_h

#include "itkMacro.h"
#include "itkConceptChecking.h"

#include "itkANTSCenteredAffine2DTransform.h"
#include "itkANTSAffine3DTransform.h"

#include <list>

namespace itk
{
/*
 */

namespace AverageAffineTransformFunctionHelperNameSpace
{
template <int D>
struct Dispatcher
{};

template <typename TAffine>
struct HelperCommonType
{
  typedef TAffine InternalAffineTransformType;

  typedef typename InternalAffineTransformType::Pointer InternalAffineTransformPointerType;
  typedef struct
  {
    InternalAffineTransformPointerType aff;
    double                             weight;
  } SingleInternalTransformItemType;

  typedef std::list<SingleInternalTransformItemType>           InternalTransformListType;
  typedef typename InternalAffineTransformType::ParametersType ParametersType;

  static void
  ComputeAveragePartialParameters(InternalTransformListType & transform_list,
                                  ParametersType &            average_parameters,
                                  unsigned int                iStart,
                                  unsigned int                iEnd, bool verbose);
};

template <typename T, typename TParametersValueType>
class HelperType;

// {
// // purposely not include any types
// };

// explicit specialization for 2D affine transform
template <typename TParametersValueType>
class HelperType<Dispatcher<2>, TParametersValueType>
{
public:
  typedef ANTSCenteredAffine2DTransform<TParametersValueType> InternalAffineTransformType;

  typedef typename HelperCommonType<InternalAffineTransformType>::InternalAffineTransformPointerType
    InternalAffineTransformPointerType;
  typedef typename HelperCommonType<InternalAffineTransformType>::SingleInternalTransformItemType
    SingleInternalTransformItemType;
  typedef typename HelperCommonType<InternalAffineTransformType>::InternalTransformListType InternalTransformListType;
  typedef typename HelperCommonType<InternalAffineTransformType>::ParametersType            ParametersType;

  static void
  ComputeAverageScaleParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageShearingParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageRotationParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageTranslationParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);
};

// explicit specialization for 3D affine transform
template <typename TParametersValueType>
class HelperType<Dispatcher<3>, TParametersValueType>
{
public:
  typedef ANTSAffine3DTransform<TParametersValueType> InternalAffineTransformType;

  typedef typename HelperCommonType<InternalAffineTransformType>::InternalAffineTransformPointerType
    InternalAffineTransformPointerType;
  typedef typename HelperCommonType<InternalAffineTransformType>::SingleInternalTransformItemType
    SingleInternalTransformItemType;
  typedef typename HelperCommonType<InternalAffineTransformType>::InternalTransformListType InternalTransformListType;
  typedef typename HelperCommonType<InternalAffineTransformType>::ParametersType            ParametersType;

  static void
  ComputeAverageScaleParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageShearingParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageRotationParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);

  static void
  ComputeAverageTranslationParameters(InternalTransformListType & transform_list, ParametersType & average_parameters, bool verbose);
};
} // namespace AverageAffineTransformFunctionHelperNameSpace

template <typename TTransform>
class AverageAffineTransformFunction
{
public:
  typedef TTransform GenericAffineTransformType;

  static constexpr unsigned int InputSpaceDimension = GenericAffineTransformType::InputSpaceDimension;
  static constexpr unsigned int OutputSpaceDimension = GenericAffineTransformType::OutputSpaceDimension;
  static constexpr unsigned int SpaceDimension = InputSpaceDimension;

  typedef typename TTransform::ParametersValueType  InternalScalarType;
  typedef Point<InternalScalarType, SpaceDimension> PointType;

  AverageAffineTransformFunction();
  ~AverageAffineTransformFunction() = default;

  // void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename GenericAffineTransformType::Pointer GenericAffineTransformPointerType;

  typedef struct
  {
    GenericAffineTransformPointerType aff;
    double                            weight;
  } SingleTransformItemType;

  typedef std::list<SingleTransformItemType> TransformListType;

  TransformListType &
  GetTransformList()
  {
    return m_TransformList;
  }

  void
  PrintTransformList();

  void
  PushBackAffineTransform(const GenericAffineTransformType * t, double weight);

  void
  AverageMultipleAffineTransform(const PointType & center_output, GenericAffineTransformPointerType & affine_output);

  /** Whether progress and debugging information are printed to standard output (cout). */
  bool verbose = false;

  /** Set to false to ignore rotation and translation. */
  bool useRigid = true;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1, (Concept::SameDimension<InputSpaceDimension, OutputSpaceDimension>));
  /** End concept checking */
#endif
protected:
  TransformListType m_TransformList;

  // type declaration to include support both 2D and 3D affine transform
protected:
  typedef typename ::itk::AverageAffineTransformFunctionHelperNameSpace::HelperType<
    ::itk::AverageAffineTransformFunctionHelperNameSpace::Dispatcher<SpaceDimension>,
    InternalScalarType>
    HelperType;

  typedef typename HelperType::InternalAffineTransformType        InternalAffineTransformType;
  typedef typename HelperType::InternalAffineTransformPointerType InternalAffineTransformPointerType;
  typedef typename HelperType::SingleInternalTransformItemType    SingleInternalTransformItemType;
  typedef typename HelperType::InternalTransformListType          InternalTransformListType;

  InternalTransformListType m_InternalTransformList;

  void
  ConvertGenericAffineToInternalAffineByFixingCenter(GenericAffineTransformPointerType &  aff,
                                                     InternalAffineTransformPointerType & iaff,
                                                     const PointType &                    center);

  void
  ConvertInternalAffineToGenericAffine(InternalAffineTransformPointerType & iaff,
                                       GenericAffineTransformPointerType &  aff);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkAverageAffineTransformFunction.hxx"
#endif

#endif /*__itkAverageAffineTransformFunction_h*/
