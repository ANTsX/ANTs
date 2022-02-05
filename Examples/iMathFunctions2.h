/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __iMathFunctions2_h
#define __iMathFunctions2_h

#include "antsUtilities.h"
#include "itkFlatStructuringElement.h"

namespace ants
{


#define iMathGetLargestComponentMinSize 50;


// Templated functions that perform the work for
// iMath.cxx and iMath.cpp (in ANTSR)
// after each function, suggested default parameters are defined

// iMathFillHolesHoleParam 2;

// Return the largest connected component in a mask
template <typename ImageType>
typename ImageType::Pointer
iMathGetLargestComponent(typename ImageType::Pointer image, /*3*/
                         unsigned long               minSize);
#define iMathGetLargestComponentMinSize 50;


// Grayscale morphological erosion
template <typename ImageType>
typename ImageType::Pointer
iMathGE(typename ImageType::Pointer image, unsigned long radius); ///*3*/
#define iMathGERadius 1;
#define iMathGEValue 1;

// Grayscale morphological opening
template <typename ImageType>
typename ImageType::Pointer
iMathGO(typename ImageType::Pointer image, unsigned long radius); /*3*/
#define iMathGORadius 1;
#define iMathGOValue 1;

template <typename ImageType>
typename ImageType::Pointer
iMathGrad(typename ImageType::Pointer image, double sigma, bool normalize); /*3*/
#define iMathGradSigma 0.5;
#define iMathGradNormalize false;

template <typename ImageType>
typename ImageType::Pointer
iMathHistogramEqualization(typename ImageType::Pointer image, double, double, unsigned int); /*3*/

// Grayscale morphological dilation
template <typename ImageType>
typename ImageType::Pointer
iMathGD(typename ImageType::Pointer image, unsigned long radius);
  /*0*/ /*3*/
#define iMathGDRadius 1;
#define iMathGDValue 1;


} // namespace ants
#include "iMathFunctions2.hxx"

#endif
