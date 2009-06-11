/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: PrintHeader.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/05 20:09:47 $
  Version:   $Revision: 1.4 $

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
int PrintHeader(int argc, char *argv[])
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
  std::cout << " Spacing " << reader->GetOutput()->GetSpacing() << std::endl;
  std::cout << " Origin " << reader->GetOutput()->GetOrigin() << std::endl;
  std::cout << " Direction " << std::endl << reader->GetOutput()->GetDirection() << std::endl;
  if( ImageDimension == 1 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " " <<  << std::endl;
    }
  else if( ImageDimension == 2 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << " " << std::endl;
    }
  else if( ImageDimension == 2 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << " " <<  " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    }
  else
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
    }
//  std::cout << " Orientation " << reader->GetOutput()->GetOrientation() << std::endl;

  return 1;
}

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cout << "Useage ex:  " << argv[0] << " image.ext " << std::endl;
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
      PrintHeader<2>(argc, argv);
      }
      break;
    case 3:
      {
      PrintHeader<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
