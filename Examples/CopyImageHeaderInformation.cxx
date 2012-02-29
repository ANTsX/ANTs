/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: CopyImageHeaderInformation.cxx,v $
  Language:  C++
  Date:      $Date: 2009/04/30 18:32:36 $
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
#include "ReadWriteImage.h"
#include "TensorFunctions.h"

template <unsigned int ImageDimension>
int CopyImageHeaderInformation(int argc, char *argv[])
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
  //  std::cout << " Spacing " << reader->GetOutput()->GetSpacing() << std::endl;
  // std::cout << " Origin " << reader->GetOutput()->GetOrigin() << std::endl;
  // std::cout << " Direction " << std::endl << reader->GetOutput()->GetDirection() << std::endl;
  // std::cout << " Size " << std::endl << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

  bool istensor = false;
  if( argc > 7 )
    {
    if( atoi(argv[7]) )
      {
      istensor = true;
      }
    }
  if( istensor )
    {
    typedef itk::Vector<float, 6>                  TensorType;
    typedef itk::Image<TensorType, ImageDimension> TensorFieldType;
    typename TensorFieldType::Pointer timage;
    ReadTensorImage<TensorFieldType>(timage, argv[2], false);
    //      std::cout<< " tim dir " << timage->GetDirection() << std::endl;
    if( argc > 6 )
      {
      if( atoi(argv[6]) )
        {
        timage->SetSpacing(  reader->GetOutput()->GetSpacing()  );
        }
      }
    if( argc > 5 )
      {
      if( atoi(argv[5]) )
        {
        timage->SetOrigin(  reader->GetOutput()->GetOrigin()  );
        }
      }
    if( argc > 4 )
      {
      if( atoi(argv[4]) )
        {
        timage->SetDirection(  reader->GetOutput()->GetDirection()  );
        }
      }

    //      std::cout<< " tim dir " << timage->GetDirection() << std::endl;
    WriteTensorImage<TensorFieldType>( timage, argv[3], false);

    return 0;
    }

  typename readertype::Pointer reader2 = readertype::New();
  reader2->SetFileName(argv[2]);
  reader2->Update();

  // MakeNewImage(typename TImage::Pointer image1, typename TImage::PixelType initval)
  typename ImageType::Pointer newimage = MakeNewImage<ImageType>(reader2->GetOutput(), -1);

  if( argc > 6 )
    {
    if( atoi(argv[6]) )
      {
      newimage->SetSpacing(  reader->GetOutput()->GetSpacing()  );
      }
    }
  if( argc > 5 )
    {
    if( atoi(argv[5]) )
      {
      newimage->SetOrigin(  reader->GetOutput()->GetOrigin()  );
      }
    }
  if( argc > 4 )
    {
    if( atoi(argv[4]) )
      {
      newimage->SetDirection(  reader->GetOutput()->GetDirection()  );
      }
    }

  WriteImage<ImageType>(newimage, argv[3]);

  return 1;
}

int main(int argc, char *argv[])
{
  if( argc < 4  )
    {
    std::cout << "Usage:  " << argv[0]
              <<
      " refimage.ext imagetocopyrefimageinfoto.ext imageout.ext   boolcopydirection  boolcopyorigin boolcopyspacing  {bool-Image2-IsTensor}"
              << std::endl;
    return 1;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim = imageIO->GetNumberOfDimensions();

  switch( dim  )
    {
    case 2:
      {
      CopyImageHeaderInformation<2>(argc, argv);
      }
      break;
    case 3:
      {
      CopyImageHeaderInformation<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension : " << dim << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
