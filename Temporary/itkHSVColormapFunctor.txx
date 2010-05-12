/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkHSVColormapFunctor.txx,v $
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
#ifndef __itkHSVColormapFunctor_txx
#define __itkHSVColormapFunctor_txx

#include "itkHSVColormapFunctor.h"

namespace itk
{
namespace Functor
{
template <class TScalar, class TRGBPixel>
typename HSVColormapFunctor<TScalar, TRGBPixel>::RGBPixelType
HSVColormapFunctor<TScalar, TRGBPixel>
::operator()( const TScalar & v ) const
{
  // Map the input scalar between [0, 1].
  RealType value = this->RescaleInputValue( v );

  // Apply the color mapping.
  // Apply the color mapping.
  RealType red = vnl_math_abs( 5.0 * ( value - 0.5 ) ) - 5.0 / 6.0;

  red = vnl_math_min( red, 1.0 );
  red = vnl_math_max( 0.0, red );

  RealType green = -vnl_math_abs( 5.0 * ( value - 11.0 / 30.0 ) ) + 11.0 / 6.0;
  green = vnl_math_min( green, 1.0 );
  green = vnl_math_max( 0.0, green );

  RealType blue = -vnl_math_abs( 5.0 * ( value - 19.0 / 30.0 ) ) + 11.0 / 6.0;
  blue = vnl_math_min( blue, 1.0 );
  blue = vnl_math_max( 0.0, blue );

  // Set the rgb components after rescaling the values.
  RGBPixelType pixel;

  pixel[0] = this->RescaleRGBComponentValue( red );
  pixel[1] = this->RescaleRGBComponentValue( green );
  pixel[2] = this->RescaleRGBComponentValue( blue );

  return pixel;
}
} // end namespace Functor
} // end namespace itk

#endif
