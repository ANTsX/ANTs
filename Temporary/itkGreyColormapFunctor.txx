/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkGreyColormapFunctor.txx,v $
  Language:  C++
  Date:      $Date: 2009-05-15 02:47:59 $
  Version:   $Revision: 1.1 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGreyColormapFunctor_txx
#define __itkGreyColormapFunctor_txx

#include "itkGreyColormapFunctor.h"

namespace itk
{
namespace Functor
{
template <class TScalar, class TRGBPixel>
typename GreyColormapFunctor<TScalar, TRGBPixel>::RGBPixelType
GreyColormapFunctor<TScalar, TRGBPixel>
::operator()( const TScalar & v ) const
{
  // Map the input scalar between [0, 1].
  RealType value = this->RescaleInputValue( v );

  // Set the rgb components after rescaling the values.
  RGBPixelType pixel;

  pixel[0] = this->RescaleRGBComponentValue( value );
  pixel[1] = pixel[0];
  pixel[2] = pixel[0];

  return pixel;
}
} // end namespace Functor
} // end namespace itk

#endif
