/*=========================================================================

  Program:   Advanced Normalization Tools
  Date:      $$

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsListSampleToListSampleFilter_h
#define __antsListSampleToListSampleFilter_h

#include "itkProcessObject.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class ListSampleToListSampleFilter
 * \brief Base class for filters that take a list sample as an input and output
 * another list sample.
 *
 * ListSampleToListSampleFilter is the base class for all process objects that output
 * list sample data, and require list sample data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ListSampleFilters
 *
 */
template <typename TInputListSample, typename TOutputListSample = TInputListSample>
class ListSampleToListSampleFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ListSampleToListSampleFilter Self;
  typedef ProcessObject                Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ListSampleToListSampleFilter);

  /** Some convenient typedefs. */
  typedef TInputListSample  InputListSampleType;
  typedef TOutputListSample OutputListSampleType;

  /** Set the list sample input of this object.  */
  virtual void
  SetInputListSample(const InputListSampleType * input);

  /** Get the list sample input of this object.  */
  InputListSampleType *
  GetInput();

  /** Get the list sample output of this object.  */
  OutputListSampleType *
  GetOutput();

  void
  Update() override
  {
    this->GenerateData();
  }

protected:
  ListSampleToListSampleFilter();
  ~ListSampleToListSampleFilter() override = default;

  void
  GenerateData() override = 0;

  void
  AllocateOutput();

private:
  ListSampleToListSampleFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  //   typename InputListSampleType::ConstPointer              m_InputListSample;
  //   typename OutputListSampleType::Pointer                  m_OutputListSample;
};
} // end namespace Statistics
} // end namespace ants
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsListSampleToListSampleFilter.hxx"
#endif

#endif
