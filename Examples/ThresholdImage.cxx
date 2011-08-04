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

#include <cstdlib>
#include <ctime>
#include <iostream>

// Software Guide : BeginCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkExtractImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"

#include <fstream>

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
typename TImage::Pointer BinaryThreshold(
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
  typename   ImageType::Pointer     Image = ImageType::New();
  Image->SetLargestPossibleRegion(input->GetLargestPossibleRegion()  );
  Image->SetBufferedRegion(input->GetLargestPossibleRegion() );
  Image->Allocate();
  Image->SetSpacing(input->GetSpacing() );
  Image->SetOrigin(input->GetOrigin() );
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
typename TImage::Pointer
DanielssonDistanceMap(
  typename TImage::PixelType pixlo,
  typename TImage::PixelType pixhi,
  typename TImage::Pointer input)
{
  std::cout << " DDMap " << std::endl;

  typedef TImage ImageType;

  typedef itk::DanielssonDistanceMapImageFilter<
      ImageType, ImageType>  FilterType;

  typename  FilterType::Pointer filter = FilterType::New();
  filter->InputIsBinaryOn();
  filter->SetUseImageSpacing(true);
  filter->SetInput(BinaryThreshold<TImage>(pixlo, pixhi, pixhi, input) );
  filter->Update();

//  std::string fn="C:\\Data\\temp.img";
//  WriteImage(filter->GetOutput(),fn.c_str());
//  fn="C:\\Data\\temp2.img";
//  WriteImage(filter->GetVoronoiMap(),fn.c_str());

  return filter->GetOutput();
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
  typedef   itk::ImageFileReader<FixedImageType>    FixedReaderType;
  typename FixedReaderType::Pointer fixedReader = FixedReaderType::New();
  fixedReader->SetFileName( argv[2] );

  typedef   itk::ImageFileWriter<FixedImageType> MovingWriterType;
  typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();
  typename MovingWriterType::Pointer movingWriter2 = MovingWriterType::New();
  movingWriter->SetFileName( argv[3] );

  try
    {
    fixedReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  // Label the surface of the image
  typename FixedImageType::Pointer thresh;
  std::string threshtype = std::string(argv[4]);
  if( strcmp(threshtype.c_str(), "Otsu") == 0 )
    {
    thresh = OtsuThreshold<FixedImageType>(atoi(argv[5]), fixedReader->GetOutput() );
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
    thresh = BinaryThreshold<FixedImageType>(atof(argv[4]), atof(argv[5]),
                                             insideValue, outsideValue, fixedReader->GetOutput() );
    }

  movingWriter->SetInput(thresh);
  movingWriter->Write();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << "   ImageDimension ImageIn.ext outImage.ext  threshlo threshhi <insideValue> <outsideValue>"
              << std::endl;
    std::cerr << "   ImageDimension ImageIn.ext outImage.ext  Otsu NumberofThresholds " << std::endl;

    std::cout << " Inclusive thresholds " << std::endl;
    return 1;
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
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
