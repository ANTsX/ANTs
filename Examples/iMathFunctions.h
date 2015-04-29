/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __iMathFunctions_h
#define __iMathFunctions_h

#include "antsUtilities.h"

namespace ants
{
// Templated functions that perform the work for
// iMath.cxx and iMath.cpp (in ANTSR)
template <class ImageType>
typename ImageType::Pointer
iMathCanny(typename ImageType::Pointer image );

template <class ImageType>
typename ImageType::Pointer
iMathGetLargestComponent(typename ImageType::Pointer image,
                         unsigned long minSize );
#define iMathGetLargestComponentMinSize 50;

template <class ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer image, unsigned long radius );
#define iMathMDRadius 1;

template <class ImageType>
typename ImageType::Pointer
iMathME(typename ImageType::Pointer image, unsigned long radius );
#define iMathMERadius 1;

template <class ImageType>
typename ImageType::Pointer
iMathMO(typename ImageType::Pointer image, unsigned long radius );
#define iMathMORadius 1;

template <class ImageType>
typename ImageType::Pointer
iMathNormalize( typename ImageType::Pointer image );
}

#include "iMathFunctions.hxx"

#endif
