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
#include "itkFlatStructuringElement.h"

namespace ants
{
// Templated functions that perform the work for
// iMath.cxx and iMath.cpp (in ANTSR)
// after each function, suggested default parameters are defined

unsigned int morph_shape_flag( const char * shape );

template<class ImageType>
typename ImageType::Pointer
iMathBlobDetector( typename ImageType::Pointer image, unsigned int nBlobs);

// Canny Edge Filter
template <class ImageType>
typename ImageType::Pointer
iMathCanny(typename ImageType::Pointer image,
           double sigma,
           double lowerThreshold,
           double upperThreshold );

// Distance Map
template <class ImageType>
typename ImageType::Pointer
iMathDistanceMap(typename ImageType::Pointer image, bool useSpacing );
#define iMathDistanceMapUseSpacing true;

// Fill Holes in objects
template <class ImageType>
typename ImageType::Pointer
iMathFillHoles(typename ImageType::Pointer image, double holeParam );
#define iMathFillHolesHoleParam 2;

// Return the largest connected component in a mask
template <class ImageType>
typename ImageType::Pointer
iMathGetLargestComponent(typename ImageType::Pointer image,
                         unsigned long minSize );
#define iMathGetLargestComponentMinSize 50;

template <unsigned int ImageDimension>
typename itk::FlatStructuringElement<ImageDimension>
iMathGetFlatStructuringElement( unsigned int shape,
                                unsigned long radius,
                                bool radiusIsParametric,
                                unsigned int lines,
                                unsigned int thickness,
                                bool includeCenter );
#define iMathGetFlatStructuringElementShape 1;
#define iMathGetFlatStructuringElementRadius 1;
#define iMathGetFlatStructuringElementLines 3;
#define iMathGetFlatStructuringElementThickness 1;
#define iMathGetFlatStructuringElementIncludeCenter false;
#define iMathGetFlatStructuringElementRadiusIsParametric false;

// Morphological Closing
template <class ImageType>
typename ImageType::Pointer
iMathMC(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType closeValue, unsigned int shape,
        bool radiusIsParametric, unsigned lines, unsigned int thickness,
        bool includeCenter );
#define iMathMCRadius 1;
#define iMathMCValue 1;

// Morphological dilation
template <class ImageType>
typename ImageType::Pointer
iMathMD(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType dilateValue, unsigned int shape,
        bool radiusIsParametric, unsigned lines, unsigned int thickness,
        bool includeCenter );
#define iMathMDRadius 1;
#define iMathMDValue 1;

// Morphological erosion
template <class ImageType>
typename ImageType::Pointer
iMathME(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType erodeValue, unsigned int shape,
        bool radiusIsParametric, unsigned lines, unsigned int thickness,
        bool includeCenter );
#define iMathMERadius 1;
#define iMathMEValue 1;

// Morphological opening
template <class ImageType>
typename ImageType::Pointer
iMathMO(typename ImageType::Pointer image, unsigned long radius,
        typename ImageType::PixelType openValue, unsigned int shape,
        bool radiusIsParametric, unsigned lines, unsigned int thickness,
        bool includeCenter );
#define iMathMORadius 1;
#define iMathMOValue 1;

// Maurer distance - returns Euclidean distance to binary object
template <class ImageType>
typename ImageType::Pointer
iMathMaurerDistance(typename ImageType::Pointer image,
                    typename ImageType::PixelType foreground );
#define iMathMaurerDistanceForeground 1;

// Grayscale morphological closing
template <class ImageType>
typename ImageType::Pointer
iMathGC(typename ImageType::Pointer image, unsigned long radius);
#define iMathGCRadius 1;
#define iMathGCValue 1;

// Grayscale morphological dilation
template <class ImageType>
typename ImageType::Pointer
iMathGD(typename ImageType::Pointer image, unsigned long radius);
#define iMathGDRadius 1;
#define iMathGDValue 1;

// Grayscale morphological erosion
template <class ImageType>
typename ImageType::Pointer
iMathGE(typename ImageType::Pointer image, unsigned long radius);
#define iMathGERadius 1;
#define iMathGEValue 1;

// Grayscale morphological opening
template <class ImageType>
typename ImageType::Pointer
iMathGO(typename ImageType::Pointer image, unsigned long radius);
#define iMathGORadius 1;
#define iMathGOValue 1;

template <class ImageType>
typename ImageType::Pointer
iMathGrad( typename ImageType::Pointer image, double sigma, bool normalize );
#define iMathGradSigma 0.5;
#define iMathGradNormalize false;

template <class ImageType>
typename ImageType::Pointer
iMathHistogramEqualization( typename ImageType::Pointer image, double, double, unsigned int );

template <class ImageType>
typename ImageType::Pointer
iMathLaplacian( typename ImageType::Pointer image, double sigma, bool normalize );
#define iMathLaplacianSigma 0.5;
#define iMathLaplacianNormalize false;

// Normalize intensity values to lie in [0,1]
template <class ImageType>
typename ImageType::Pointer
iMathNormalize( typename ImageType::Pointer image );

template <class ImageType>
typename ImageType::Pointer
iMathPad( typename ImageType::Pointer image, int padding );

template <class ImageType>
typename ImageType::Pointer
iMathPeronaMalik( typename ImageType::Pointer image, unsigned long nIterations,
                  double conductance );
#define iMathPeronaMalikConductance 0.25;
#define iMathPeronaMalikNIterations 1;

template <class ImageType>
typename ImageType::Pointer
iMathSharpen( typename ImageType::Pointer image );

template <class ImageType>
typename ImageType::Pointer
iMathTruncateIntensity( typename ImageType::Pointer image, double lowerQ,
                        double upperQ, int nBins,
                        typename itk::Image<int, ImageType::ImageDimension>::Pointer maskImage );
#define iMathTruncateIntensityNBins 64;


}
#include "iMathFunctions.hxx"

#endif
