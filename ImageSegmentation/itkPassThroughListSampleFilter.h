/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPassThroughListSampleFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPassThroughListSampleFilter_h
#define __itkPassThroughListSampleFilter_h

#include "itkListSampleToListSampleFilter.h"

namespace itk
{
namespace Statistics
{
/** \class PassThroughListSampleFilter
 * \brief Simple class which pass the input to the output.
 *
 */

template <class TListSample>
class ITK_EXPORT PassThroughListSampleFilter
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
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPassThroughListSampleFilter.txx"
#endif

#endif
