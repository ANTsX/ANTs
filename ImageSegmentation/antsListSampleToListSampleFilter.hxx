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
#ifndef _itkListSampleToListSampleFilter_hxx
#define _itkListSampleToListSampleFilter_hxx


namespace itk
{
namespace ants
{
namespace Statistics
{
/**
 *
 */
template <typename TInputListSample, typename TOutputListSample>
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::ListSampleToListSampleFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);
}

template <typename TInputListSample, typename TOutputListSample>
void
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::SetInputListSample(const TInputListSample * input)
{
  //   this->m_InputListSample = const_cast<InputListSampleType *>( input );
  this->ProcessObject::SetNthInput(0, reinterpret_cast<DataObject *>(const_cast<InputListSampleType *>(input)));
}

template <typename TInputListSample, typename TOutputListSample>
void
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::AllocateOutput()
{
  typename DataObject::Pointer obj = reinterpret_cast<DataObject *>(TOutputListSample::New().GetPointer());

  //   typename TOutputListSample::Pointer output
  //     = reinterpret_cast<TOutputListSample*>(obj.GetPointer());
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, obj.GetPointer());

  //   this->m_OutputListSample = OutputListSampleType::New();
}

/**
 *
 */
template <typename TInputListSample, typename TOutputListSample>
typename ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::InputListSampleType *
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::GetInput()
{
  return reinterpret_cast<TInputListSample *>(this->ProcessObject::GetInput(0));
}

template <typename TInputListSample, typename TOutputListSample>
typename ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::OutputListSampleType *
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::GetOutput()
{
  if (this->GetNumberOfOutputs() < 1)
  {
    return nullptr;
  }

  // we assume that the first output is of the templated type
  return reinterpret_cast<TOutputListSample *>(this->ProcessObject::GetOutput(0));
  //   return this->m_OutputListSample.GetPointer();
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
