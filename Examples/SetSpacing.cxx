/*
 * SetSpacing.cxx
 *
 *  Created on: Mar 5, 2009
 *      Author: songgang
 */

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SetSpacing.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/06 03:57:48 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

template <unsigned int ImageDimension>
int SetSpacing(int argc, char *argv[])
{
  typedef  float                                     outPixelType;
  typedef  float                                     floatPixelType;
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::Image<floatPixelType, ImageDimension> IntermediateType;
  typedef itk::Image<outPixelType, ImageDimension>   OutImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;
  typedef itk::ImageFileWriter<OutImageType>         writertype;

  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  typename OutImageType::Pointer outim = reader->GetOutput();
  typename OutImageType::SpacingType spacing = outim->GetSpacing();

  std::cout << " Old Spacing " <<  outim->GetSpacing();
  if( argc > 3 )
    {
    spacing[0] = atof(argv[3]);
    }
  if( argc > 4 )
    {
    spacing[1] = atof(argv[4]);
    }
  if( argc > 5 )
    {
    spacing[2] = atof(argv[5]);
    }
  std::cout << "  New Spacing " << spacing << std::endl;

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typename ImageType::Pointer varimage = ImageType::New();
  varimage->SetLargestPossibleRegion( outim->GetLargestPossibleRegion() );
  varimage->SetBufferedRegion( outim->GetLargestPossibleRegion() );
  varimage->SetLargestPossibleRegion( outim->GetLargestPossibleRegion() );
  varimage->Allocate();
  varimage->SetSpacing(spacing);
  varimage->SetOrigin(outim->GetOrigin() );
  varimage->SetDirection( outim->GetDirection() );
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    vfIter2.Set(outim->GetPixel(vfIter2.GetIndex() ) );
    }

  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(  varimage );
  writer->Update();
  writer->Write();

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Usage:   " << argv[0] << "  Dimension infile.hdr outfile.nii  SpacingX SpacingY {SpacingZ} "
              << std::endl;
    return 1;
    }

  // Get the image dimension
  switch( atoi(argv[1]) )
    {
    case 2:
      {
      SetSpacing<2>(argc - 1, argv + 1);
      }
      break;
    case 3:
      {
      SetSpacing<3>(argc - 1, argv + 1);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
