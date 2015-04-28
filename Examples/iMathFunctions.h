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
iMathGetLargestComponent(typename ImageType::Pointer image,
                    unsigned long smallest );

}

#include "iMathFunctions.hxx"

#endif
