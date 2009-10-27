/*
 * ResetDirection2.cxx
 *
 *  Created on: Nov 14, 2008
 *      Author: songgang
 */

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ResetDirection2.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/14 20:47:46 $
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

// #include "../itkAvants.DefineDimension"

template <unsigned int ImageDimension>
int ResetDirection(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Useage ex:   " << argv[0] << "  infile.nii outfile.nii direction matrix in a row " << std::endl;
    return 1;
    }

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
  typename OutImageType::DirectionType direction = outim->GetDirection();
  // direction->SetIdentity();
  // direction.Fill(0);
  // for (unsigned int i=0;i<ImageDimension;i++) direction[i][i]=1;

  char* * dtext = argv + 3;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      direction[i][j] = atof(dtext[i * ImageDimension + j]);
      std::cout << "direction[" << i << "][" << j << "]=" << direction[i][j] << std::endl;
      }
    }

  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typename ImageType::Pointer varimage = ImageType::New();
  varimage->SetLargestPossibleRegion( outim->GetLargestPossibleRegion() );
  varimage->SetBufferedRegion( outim->GetLargestPossibleRegion() );
  varimage->SetLargestPossibleRegion( outim->GetLargestPossibleRegion() );
  varimage->Allocate();
  varimage->SetSpacing(outim->GetSpacing() );
  varimage->SetOrigin(outim->GetOrigin() );
  varimage->SetDirection( direction );
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
    std::cout << "Useage ex:   " << argv[0] << "  infile.nii outfile.nii  d01 d02 d03 d10 .... " << std::endl;
    return 1;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();

  switch( imageIO->GetNumberOfDimensions() )
    {
    case 2:
      {
      ResetDirection<2>(argc, argv);
      }
      break;
    case 3:
      {
      ResetDirection<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
