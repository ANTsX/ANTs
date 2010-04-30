/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ConvertToJpg.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.19 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <limits.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

template <unsigned int ImageDimension, class TPIXELTYPE>
int ConvertType(int argc, char *argv[], double MINVAL, double MAXVAL)
{
  typedef  TPIXELTYPE                                outPixelType;
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
  const double desiredMinimum =  MINVAL;
  double       desiredMaximum =  MAXVAL;
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
    std::cout << "Usage:   " << argv[0] << " infile.nii out.ext TYPE-OPTION " << std::endl;
    std::cout << " ext is the extension you want, e.g. tif.  " << std::endl;
    std::cout << " TYPE-OPTION  :  TYPE " << std::endl;
    std::cout << "  0  :  char   " << std::endl;
    std::cout << "  1  :  unsigned char   " << std::endl;
    std::cout << "  2  :  short   " << std::endl;
    std::cout << "  3  :  unsigned short   " << std::endl;
    std::cout << "  4  :  int   " << std::endl;
    std::cout << "  5  :  unsigned int   " << std::endl;
    std::cout
      << " Note that some pixel types are not supported by some image formats. e.g.  int is not supported by jpg. "
      << std::endl;
    std::cout << " You can easily extend this for other pixel types with a few lines of code and adding usage info. "
              << std::endl;
    std::cout
      <<
    " The image intensity will be scaled to the dynamic range of the pixel type.  E.g. uchar => 0  (min), 255 (max). "
      << std::endl;
    return 1;
    }
  unsigned int typeoption = 0;
  if( argc > 3 )
    {
    typeoption = atoi(argv[3]);
    }
  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();

  if( typeoption == 0 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, char>(argc, argv, SCHAR_MIN, SCHAR_MAX );
        }
        break;
      case 3:
        {
        ConvertType<3, char>(argc, argv,  SCHAR_MIN, SCHAR_MAX );
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }
  else if( typeoption == 1 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, unsigned char>(argc, argv,  0, UCHAR_MAX );
        }
        break;
      case 3:
        {
        ConvertType<3, unsigned char>(argc, argv,  0, UCHAR_MAX );
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }
  else if( typeoption == 2 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, short>(argc, argv,  SHRT_MIN, SHRT_MAX);
        }
        break;
      case 3:
        {
        ConvertType<3, short>(argc, argv,  SHRT_MIN, SHRT_MAX);
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }
  else if( typeoption == 3 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, unsigned short>(argc, argv,  0, USHRT_MAX );
        }
        break;
      case 3:
        {
        ConvertType<3, unsigned short>(argc, argv,  0, USHRT_MAX );
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }
  else if( typeoption == 4 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, int>(argc, argv,  INT_MIN, INT_MAX);
        }
        break;
      case 3:
        {
        ConvertType<3, int>(argc, argv,  INT_MIN, INT_MAX);
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }
  else if( typeoption == 5 )
    {
    switch( imageIO->GetNumberOfDimensions() )
      {
      case 2:
        {
        ConvertType<2, unsigned int>(argc, argv,  0, UINT_MAX );
        }
        break;
      case 3:
        {
        ConvertType<3, unsigned int>(argc, argv,  0, UINT_MAX );
        }
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }
    }

  return 0;
}
