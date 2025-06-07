/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPointSetFunction_h
#define __itkPointSetFunction_h

#include "itkFunctionBase.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkContinuousIndex.h"
#include "itkImageBase.h"

namespace itk
{
/** \class PointSetFunction
 * \brief Evaluates a function of an image at specified position.
 *
 * PointSetFunction is a baseclass for all objects that evaluates
 * a function of an image at index, continuous index or point.
 * This class is templated over the input image type, the type
 * of the function output and the coordinate representation type
 * (e.g. float or double).
 *
 * The input image is set via method SetInputPointSet().
 * Methods Evaluate, EvaluateAtIndex and EvaluateAtContinuousIndex
 * respectively evaluates the function at an geometric point,
 * image index and continuous image index.
 *
 * \warning Image BufferedRegion information is cached during
 * in SetInputPointSet( image ). If the image BufferedRegion has changed
 * one must call SetInputPointSet( image ) again to update the cache
 * to the current values.
 *
 * \sa Point
 * \sa Index
 * \sa ContinuousIndex
 *
 * \ingroup PointSetFunctions
 */
template <typename TInputPointSet, typename TOutput, typename TCoordRep = float>
class PointSetFunction : public FunctionBase<typename TInputPointSet::PointType, TOutput>
{
public:
  /** Dimension underlying input point set. */
  static constexpr unsigned int Dimension = TInputPointSet::PointDimension;

  /** Standard class typedefs. */
  typedef PointSetFunction                                          Self;
  typedef FunctionBase<typename TInputPointSet::PointType, TOutput> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PointSetFunction);

  /** InputPointSetType typedef support. */
  typedef TInputPointSet InputPointSetType;

  /** InputPixel typedef support */
  typedef typename InputPointSetType::PointType InputPointType;
  typedef typename InputPointSetType::PixelType InputPixelType;

  /** InputPointSetPointer typedef support */
  typedef typename InputPointSetType::ConstPointer InputPointSetConstPointer;

  /** OutputType typedef support. */
  typedef TOutput OutputType;

  /** CoordRepType typedef support. */
  typedef TCoordRep CoordRepType;

  /** Set the input point set.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputPointSet again to update cached values. */
  virtual void
  SetInputPointSet(const InputPointSetType * ptr);

  /** Get the input image. */
  const InputPointSetType *
  GetInputPointSet() const
  {
    return m_PointSet.GetPointer();
  }

  /** Evaluate the function at specified Point position.
   * Subclasses must provide this method. */
  virtual TOutput
  Evaluate(const InputPointType & point) const = 0;

protected:
  PointSetFunction();
  ~PointSetFunction() {}

  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** Const pointer to the input image. */
  InputPointSetConstPointer m_PointSet;

private:
  PointSetFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_PointSetFunction(_, EXPORT, x, y)            \
  namespace itk                                                   \
  {                                                               \
  _(3(class EXPORT PointSetFunction<ITK_TEMPLATE_3 x>))           \
  namespace Templates                                             \
  {                                                               \
  typedef PointSetFunction<ITK_TEMPLATE_3 x> PointSetFunction##y; \
  }                                                               \
  }

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPointSetFunction.hxx"
#endif

#endif
