/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ThresholdImage.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/27 23:45:44 $
  Version:   $Revision: 1.20 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "ReadWriteData.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMultiplyImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkResampleImageFilter.h"

namespace ants
{
template <class TImage>
typename TImage::Pointer
MultiplyImage(typename TImage::Pointer image1, typename TImage::Pointer image2)
{
  std::cout << " Multiply " << std::endl;
  // Begin Multiply Images
  typedef TImage tImageType;
  //  output will be the speed image for FMM
  typedef itk::MultiplyImageFilter<tImageType,
                                   tImageType, tImageType>  MultFilterType;
  typename MultFilterType::Pointer filter = MultFilterType::New();
  filter->SetInput1( image1 );
  filter->SetInput2( image2 );
  filter->Update();
  return filter->GetOutput();   // this is the speed image

  // write a function to threshold the speedimage so
  // if the dist is g.t. D then speed = 1
}

template <class TImage>
typename TImage::Pointer BinaryThreshold_AltInsideOutside_threashold(
  typename TImage::PixelType low,
  typename TImage::PixelType high,
  typename TImage::PixelType insideval, typename TImage::PixelType outsideval,
  typename TImage::Pointer input )
{
  std::cout << " Binary Thresh " << std::endl;

  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef itk::BinaryThresholdImageFilter<TImage, TImage> InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  inputThresholder->SetInsideValue(  insideval );
  inputThresholder->SetOutsideValue( outsideval );

  if( high < low )
    {
    high = 255;
    }
  float eps = 1.e-6 * low;
  inputThresholder->SetLowerThreshold( (PixelType) low - eps );
  inputThresholder->SetUpperThreshold( (PixelType) high + eps);
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <class TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground,
             typename TImage::PixelType newval, typename TImage::Pointer input)
{
  std::cout << " Label Surf " << std::endl;
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  // ORIENTATION ALERT: Original code set spacing & origin without
  // setting directions.
  typename   ImageType::Pointer Image = AllocImage<ImageType>(input);

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;

  typename iteratorType::RadiusType rad;
  for( int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

//  std::cout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    typename TImage::IndexType ind2;
    if( p == foreground )
      {
      bool atedge = false;
      for( int i = 0; i < GHood.Size(); i++ )
        {
        ind2 = GHood.GetIndex(i);
        float dist = 0.0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          dist += (float)(ind[j] - ind2[j]) * (float)(ind[j] - ind2[j]);
          }
        dist = sqrt(dist);
        if( GHood.GetPixel(i) != foreground && dist < 1.1 )
          {
          atedge = true;
          }
        }
      if( atedge && p == foreground )
        {
        Image->SetPixel(ind, newval);
        }
      else if( p == foreground )
        {
        Image->SetPixel(ind, 0);
        }
      }
    ++GHood;
    }

  return Image;
}

template <class TImage>
typename TImage::Pointer OtsuThreshold(
  int NumberOfThresholds, typename TImage::Pointer input)
{
  std::cout << " Otsu Thresh with " << NumberOfThresholds << " thresholds" << std::endl;

  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef itk::OtsuMultipleThresholdsImageFilter<TImage, TImage> InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  /*
  inputThresholder->SetInsideValue(  replaceval );
  int outval=0;
  if ((float) replaceval == (float) -1) outval=1;
  inputThresholder->SetOutsideValue( outval );
  */
  inputThresholder->SetNumberOfThresholds( NumberOfThresholds );

  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <unsigned int InImageDimension>
int ThresholdImage( int argc, char * argv[] )
{
  //  const     unsigned int   InImageDimension = AvantsImageDimension;
  typedef   float                                   PixelType;
  typedef   itk::Image<PixelType, InImageDimension> FixedImageType;
  typename FixedImageType::Pointer fixed;
  ReadImage<FixedImageType>( fixed, argv[2] );
  // Label the surface of the image
  typename FixedImageType::Pointer thresh;
  std::string threshtype = std::string(argv[4]);
  if( strcmp(threshtype.c_str(), "Otsu") == 0 )
    {
    thresh = OtsuThreshold<FixedImageType>(atoi(argv[5]), fixed );
    }
  else
    {
    PixelType insideValue = 1;
    PixelType outsideValue = 0;
    if( argc > 6 )
      {
      insideValue = static_cast<PixelType>( atof( argv[6] ) );
      }
    if( argc > 7 )
      {
      outsideValue = static_cast<PixelType>( atof( argv[7] ) );
      }
    thresh = BinaryThreshold_AltInsideOutside_threashold<FixedImageType>(atof(argv[4]), atof(
                                                                           argv[5]),
                                                                         insideValue, outsideValue,
                                                                         fixed );
    }

  WriteImage<FixedImageType>( thresh, argv[3] );
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ThresholdImage( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ThresholdImage" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );

  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0];
    std::cout << "   ImageDimension ImageIn.ext outImage.ext  threshlo threshhi <insideValue> <outsideValue>"
             << std::endl;
    std::cout << "   ImageDimension ImageIn.ext outImage.ext  Otsu NumberofThresholds " << std::endl;

    std::cout << " Inclusive thresholds " << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      ThresholdImage<2>(argc, argv);
      }
      break;
    case 3:
      {
      ThresholdImage<3>(argc, argv);
      }
      break;
    case 4:
      {
      ThresholdImage<4>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return 0;
}
} // namespace ants
