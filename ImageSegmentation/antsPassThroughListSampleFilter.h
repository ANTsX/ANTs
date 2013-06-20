/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsPassThroughListSampleFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsPassThroughListSampleFilter_h
#define __antsPassThroughListSampleFilter_h

#include "antsListSampleToListSampleFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class PassThroughListSampleFilter
 * \brief Simple class which pass the input to the output.
 *
 */

template <class TListSample>
class PassThroughListSampleFilter
  : public       ListSampleToListSampleFilter<TListSample, TListSample>
{
public:
  /**
   * Standard class typedefs.
   */
  typedef PassThroughListSampleFilter                            Self;
  typedef ListSampleToListSampleFilter<TListSample, TListSample> Superclass;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;

  /**
   * Standard macros
   */
  itkTypeMacro( PassThroughListSampleFilter,
                ListSampleToScalarListSampleFilter );

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro( Self );

  /**
   * Conveneient typedefs
   */
  typedef TListSample                                    ListSampleType;
  typedef ListSampleType                                 InputType;
  typedef typename ListSampleType::MeasurementVectorType MeasurementVectorType;
  typedef typename ListSampleType::MeasurementType       MeasurementType;
protected:
  PassThroughListSampleFilter();
  virtual ~PassThroughListSampleFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  virtual void GenerateData();

private:
  PassThroughListSampleFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );              // purposely not implemented
};                                             // end of class
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsPassThroughListSampleFilter.hxx"
#endif

#endif
