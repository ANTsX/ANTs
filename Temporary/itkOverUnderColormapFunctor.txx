/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkOverUnderColormapFunctor.txx,v $
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
#ifndef __itkOverUnderColormapFunctor_txx
#define __itkOverUnderColormapFunctor_txx

#include "itkOverUnderColormapFunctor.h"

namespace itk
{
namespace Functor
{
template <class TScalar, class TRGBPixel>
typename OverUnderColormapFunctor<TScalar, TRGBPixel>::RGBPixelType
OverUnderColormapFunctor<TScalar, TRGBPixel>
::operator()( const TScalar & v ) const
{
  // Map the input scalar between [0, 1].
  RealType value = this->RescaleInputValue( v );

  // Apply the color mapping.
  RealType red = value;
  RealType green = value;
  RealType blue = value;

  if( value == 0.0 )
    {
    // pixel is saturated in the dark
    red = 0.0;
    green = 0.0;
    blue = 1.0;
    }
  else if( value == 1.0 )
    {
    // pixel is saturated in the white
    red = 1.0;
    green = 0.0;
    blue = 0.0;
    }

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
