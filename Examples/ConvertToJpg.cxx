/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: ConvertToJpg.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.19 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

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
int ConvertType(int argc, char *argv[])
{
  typedef  unsigned char                             outPixelType;
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
  std::cout << " Updated reader " << std::endl;

  typedef itk::CastImageFilter<ImageType, IntermediateType> castertype;
  typename   castertype::Pointer caster = castertype::New();
  caster->SetInput(reader->GetOutput() );
  caster->Update();

  // Rescale the image intensities so that they fall between 0 and 255
  typedef itk::RescaleIntensityImageFilter<IntermediateType, IntermediateType> FilterType;
  typename   FilterType::Pointer fixedrescalefilter = FilterType::New();
  fixedrescalefilter->SetInput(caster->GetOutput() );
  const double desiredMinimum =  0.0;
  const double desiredMaximum =  255.0;
  fixedrescalefilter->SetOutputMinimum( desiredMinimum );
  fixedrescalefilter->SetOutputMaximum( desiredMaximum );
  fixedrescalefilter->UpdateLargestPossibleRegion();

  typedef itk::CastImageFilter<IntermediateType, OutImageType> castertype2;
  typename castertype2::Pointer caster2 = castertype2::New();
  caster2->SetInput(fixedrescalefilter->GetOutput() );
  caster2->Update();

  typename   OutImageType::Pointer outim = caster2->GetOutput();
  typename   OutImageType::SpacingType spc = outim->GetSpacing();
  outim->SetSpacing(spc);
  std::cout << " Dire in " << reader->GetOutput()->GetDirection() << std::endl;
  std::cout << " Dire out " << outim->GetDirection() << std::endl;
  typename   writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(outim);
  writer->Update();
  writer->Write();

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Usage:   ConvertToJpg infile.nii out.jpg " << std::endl;
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
      ConvertType<2>(argc, argv);
      }
      break;
    case 3:
      {
      ConvertType<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
