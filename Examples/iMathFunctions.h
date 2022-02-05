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
// after each function, suggested default parameters are defined

unsigned int
morph_shape_flag(const char * shape);

template <typename ImageType>
typename ImageType::Pointer
iMathBlobDetector(typename ImageType::Pointer image, unsigned int nBlobs); /*???*/

// Canny Edge Filter
template <typename ImageType>
typename ImageType::Pointer
iMathCanny(typename ImageType::Pointer image, /*0*/
           double                      sigma,
           double                      lowerThreshold,
           double                      upperThreshold);

// Distance Map
template <typename ImageType>
typename ImageType::Pointer
iMathDistanceMap(typename ImageType::Pointer image, bool useSpacing); /*0*/
#define iMathDistanceMapUseSpacing true;

// Fill Holes in objects
template <typename ImageType>
typename ImageType::Pointer
iMathFillHoles(typename ImageType::Pointer image, double holeParam); /*0*/
#define iMathFillHolesHoleParam 2;


#define iMathGetLargestComponentMinSize 50;

// Grayscale morphological closing
template <typename ImageType>
typename ImageType::Pointer
iMathGC(typename ImageType::Pointer image, unsigned long radius); /*0*/
#define iMathGCRadius 1;
#define iMathGCValue 1;


} // namespace ants
#include "iMathFunctions.hxx"

#endif
