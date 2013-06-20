/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsListSampleToListSampleFilter.h,v $
  Language:  C++
  Date:      $$
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
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
template <class TInputListSample, class TOutputListSample = TInputListSample>
class ListSampleToListSampleFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ListSampleToListSampleFilter Self;
  typedef ProcessObject                Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ListSampleToListSampleFilter, ProcessObject );

  /** Some convenient typedefs. */
  typedef TInputListSample  InputListSampleType;
  typedef TOutputListSample OutputListSampleType;

  /** Set the list sample input of this object.  */
  virtual void SetInputListSample( const InputListSampleType *input );

  /** Get the list sample input of this object.  */
  InputListSampleType * GetInput();

  /** Get the list sample output of this object.  */
  OutputListSampleType * GetOutput();

  virtual void Update()
  {
    this->GenerateData();
  }

protected:
  ListSampleToListSampleFilter();
  ~ListSampleToListSampleFilter()
  {
  };

  virtual void GenerateData() = 0;

  void AllocateOutput();

private:
  ListSampleToListSampleFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );               // purposely not implemented

//   typename InputListSampleType::ConstPointer              m_InputListSample;
//   typename OutputListSampleType::Pointer                  m_OutputListSample;
};
} // end namespace Statistics
} // end namespace ants
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsListSampleToListSampleFilter.hxx"
#endif

#endif
