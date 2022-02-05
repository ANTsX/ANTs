/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __iMathFunctions1_h
#define __iMathFunctions1_h

#include "antsUtilities.h"
#include "itkFlatStructuringElement.h"

namespace ants
{
// Templated functions that perform the work for
// iMath.cxx and iMath.cpp (in ANTSR)
// after each function, suggested default parameters are defined

#define iMathGetLargestComponentMinSize 50;

template <unsigned int ImageDimension>
typename itk::FlatStructuringElement<ImageDimension>
iMathGetFlatStructuringElement(unsigned int  shape,
                               unsigned long radius,
                               bool          radiusIsParametric,
                               unsigned int  lines,
                               unsigned int  thickness,
                               bool          includeCenter);
#define iMathGetFlatStructuringElementShape 1;
#define iMathGetFlatStructuringElementRadius 1;
#define iMathGetFlatStructuringElementLines 3;
#define iMathGetFlatStructuringElementThickness 1;
#define iMathGetFlatStructuringElementIncludeCenter false;
#define iMathGetFlatStructuringElementRadiusIsParametric false;


// Morphological Closing
template <typename ImageType>
typename ImageType::Pointer
iMathMC(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType closeValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned                      lines,
        unsigned int                  thickness,
        bool                          includeCenter);
#define iMathMCRadius 1;
#define iMathMCValue 1;

// Morphological dilation
template <typename ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType dilateValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned                      lines,
        unsigned int                  thickness,
        bool                          includeCenter);
#define iMathMDRadius 1;
#define iMathMDValue 1;

// Morphological erosion
template <typename ImageType>
typename ImageType::Pointer
iMathME(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType erodeValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned                      lines,
        unsigned int                  thickness,
        bool                          includeCenter);
#define iMathMERadius 1;
#define iMathMEValue 1;

// Morphological opening
template <typename ImageType>
typename ImageType::Pointer
iMathMO(typename ImageType::Pointer   image,
        unsigned long                 radius, /*1*/
        typename ImageType::PixelType openValue,
        unsigned int                  shape,
        bool                          radiusIsParametric,
        unsigned                      lines,
        unsigned int                  thickness,
        bool                          includeCenter);
#define iMathMORadius 1;
#define iMathMOValue 1;

// Maurer distance - returns Euclidean distance to binary object
template <typename ImageType>
typename ImageType::Pointer
iMathMaurerDistance(typename ImageType::Pointer   image, /*1*/
                    typename ImageType::PixelType foreground);
#define iMathMaurerDistanceForeground 1;


template <typename ImageType>
typename ImageType::Pointer
iMathLaplacian(typename ImageType::Pointer image, double sigma, bool normalize); /*1*/
#define iMathLaplacianSigma 0.5;
#define iMathLaplacianNormalize false;

// Normalize intensity values to lie in [0,1]
template <typename ImageType>
typename ImageType::Pointer
iMathNormalize(typename ImageType::Pointer image); /*1*/

template <typename ImageType>
typename ImageType::Pointer
iMathPad(typename ImageType::Pointer image, int padding); /*1*/

template <typename ImageType>
typename ImageType::Pointer
iMathPeronaMalik(typename ImageType::Pointer image,
                 unsigned long               nIterations, /*1*/
                 double                      conductance);
#define iMathPeronaMalikConductance 0.25;
#define iMathPeronaMalikNIterations 1;

template <typename ImageType>
typename ImageType::Pointer
iMathPropagateLabelsThroughMask(typename ImageType::Pointer mask, /*1*/
                                typename ImageType::Pointer lables,
                                double                      stoppingValue,
                                unsigned int                propagationMethod);
#define iMathPropagateLabelsThroughMaskStoppingValue 100.0;
#define iMathPropagateLabelsThroughMaskMethod 0;

template <typename ImageType>
typename ImageType::Pointer
iMathSharpen(typename ImageType::Pointer image); /*1*/

template <typename ImageType>
typename ImageType::Pointer
iMathTruncateIntensity(typename ImageType::Pointer                                  image,
                       double                                                       lowerQ, /*1*/
                       double                                                       upperQ,
                       int                                                          nBins,
                       typename itk::Image<int, ImageType::ImageDimension>::Pointer maskImage);
#define iMathTruncateIntensityNBins 64;


} // namespace ants
#include "iMathFunctions1.hxx"

#endif
