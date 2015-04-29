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
iMathMC(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType closeValue );
#define iMathMCRadius 1;
#define iMathMCValue 1;

template <class ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType dilateValue );
#define iMathMDRadius 1;
#define iMathMDValue 1;

template <class ImageType>
typename ImageType::Pointer
iMathME(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::Pointer erodeValue );
#define iMathMERadius 1;
#define iMathMEValue 1;

template <class ImageType>
typename ImageType::Pointer
iMathMO(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType openValue );
#define iMathMORadius 1;
#define iMathMOValue 1;

template <class ImageType>
typename ImageType::Pointer
iMathNormalize( typename ImageType::Pointer image );
}

#include "iMathFunctions.hxx"

#endif
