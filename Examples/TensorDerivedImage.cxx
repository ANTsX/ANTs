/*=========================================================================

  Program:   ITK General
  Module:    $RCSfile: TensorDerivedImage.cxx,v $
  Language:  C++
  Date:      $Date: 2009/03/17 18:58:59 $
  Version:   $Revision: 1.2 $

  Author: Jeffrey T. Duda (jtduda@seas.upenn.edu)
  Institution: PICSL

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/
#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// #include "itkDiffusionTensor3D.h"
#include "itkVector.h"
// #include "itkTensorFractionalAnisotropyImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <string>
#include "TensorFunctions.h"
#include "ReadWriteImage.h"
#include "itkRGBPixel.h"

int main( int argc, char *argv[] )
{
  // Pixel and Image typedefs
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, 3>              ScalarImageType;
  typedef itk::ImageFileWriter<ScalarImageType> WriterType;

  // typedef itk::Vector<PixelType, 6>             TensorType;
  typedef itk::SymmetricSecondRankTensor<double, 3> TensorType;
  typedef itk::Image<TensorType, 3>                 TensorImageType;
  typedef itk::ImageFileReader<TensorImageType>     ReaderType;
  typedef itk::RGBPixel<float>                      ColorPixelType;
  typedef itk::Image<ColorPixelType, 3>             ColorImageType;
  typedef itk::ImageFileWriter<ColorImageType>      ColorWriterType;

  // Check for valid input paramters
  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " tensorvolume outputvolume outputtype" << std::endl;
    return 1;
    }

  // Input parameters
  char *      inputName = argv[1];
  char *      outputName = argv[2];
  std::string outType = argv[3];

  TensorImageType::Pointer dtimg;
  ReadTensorImage<TensorImageType>(dtimg, inputName, false);

  std::cout << "tensor_image: " << inputName << std::endl;
  std::cout << "output_image: " << outputName << std::endl;

  ScalarImageType::Pointer outImage;
  ColorImageType::Pointer  colorImage;

  if( outType == "DEC" )
    {
    colorImage = ColorImageType::New();
    colorImage->SetRegions(dtimg->GetLargestPossibleRegion() );
    colorImage->SetSpacing(dtimg->GetSpacing() );
    colorImage->SetOrigin(dtimg->GetOrigin() );
    colorImage->SetDirection(dtimg->GetDirection() );
    colorImage->Allocate();
    }
  else
    {
    outImage = ScalarImageType::New();
    outImage->SetRegions(dtimg->GetLargestPossibleRegion() );
    outImage->SetSpacing(dtimg->GetSpacing() );
    outImage->SetOrigin(dtimg->GetOrigin() );
    outImage->SetDirection(dtimg->GetDirection() );
    outImage->Allocate();
    }

  itk::ImageRegionIteratorWithIndex<TensorImageType> inputIt(dtimg,
                                                             dtimg->GetLargestPossibleRegion() );

  if( (outType == "XX") || (outType == "xx") )
    {
    outType = "0";
    }
  if( (outType == "XY") || (outType == "YX") || (outType == "xy") || (outType == "yx") )
    {
    outType = "1";
    }
  if( (outType == "XZ") || (outType == "ZX") || (outType == "xz") || (outType == "zx") )
    {
    outType = "2";
    }
  if( (outType == "YY") || (outType == "yy") )
    {
    outType = "3";
    }
  if( (outType == "YZ") || (outType == "ZY") || (outType == "yz") || (outType == "zy") )
    {
    outType = "4";
    }
  if( (outType == "ZZ") || (outType == "zz") )
    {
    outType = "5";
    }

  std::cout << "Calculating output..." << std::flush;

  while( !inputIt.IsAtEnd() )
    {
    if( (outType == "0") || (outType == "1") || (outType == "2") || (outType == "3") ||
        (outType == "4") || (outType == "5") )
      {
      int idx = atoi(outType.c_str() );
      outImage->SetPixel(inputIt.GetIndex(), inputIt.Value()[idx]);
      }
    else if( (outType == "TR") || (outType == "MD") )
      {
      ScalarImageType::PixelType tr;
      tr = inputIt.Value()[0] + inputIt.Value()[2] + inputIt.Value()[5];
      if( tr < 0 )
        {
        tr = 0;
        }
      if( tr != tr )
        {
        tr = 0;
        }
      if( outType == "TR" )
        {
        outImage->SetPixel(inputIt.GetIndex(), tr);
        }
      else
        {
        outImage->SetPixel(inputIt.GetIndex(), tr / 3.0);
        }
      }
    else if( outType == "V" )
      {
      // unsigned int invalids = 0;
      bool       success = true;
      TensorType t = TensorLogAndExp<TensorType>(inputIt.Value(), true, success);
      int        current = 1;
      if( !success )
        {
        current = 0;
        }
      outImage->SetPixel(inputIt.GetIndex(), current);
      // std::cout << "Found " << invalids << " invalid tensors" << std::endl;
      }
    else if( outType == "FA" )
      {
      float fa = GetTensorFA<TensorType>(inputIt.Value() );
      outImage->SetPixel(inputIt.GetIndex(), fa);
      }
    else if( outType == "DEC" )
      {
      colorImage->SetPixel(inputIt.GetIndex(),
                           GetTensorRGB<TensorType>(inputIt.Value() ) );
      ++inputIt;
      }
    }

  std::cout << "Done. " << std::endl;

  if( outType == "DEC" )
    {
    std::cout << "output origin: " << colorImage->GetOrigin() << std::endl;
    std::cout << "output size: " << colorImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "output spacing: " << colorImage->GetSpacing() << std::endl;
    std::cout << "output direction: " << colorImage->GetDirection() << std::endl;
    }
  else
    {
    std::cout << "output origin: " << outImage->GetOrigin() << std::endl;
    std::cout << "output size: " << outImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "output spacing: " << outImage->GetSpacing() << std::endl;
    std::cout << "output direction: " << outImage->GetDirection() << std::endl;
    }

  if( outType == "DEC" )
    {
    ColorWriterType::Pointer writer = ColorWriterType::New();
    writer->SetInput( colorImage );
    writer->SetFileName( outputName );
    writer->Update();
    }
  else
    {
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( outImage );
    writer->SetFileName( outputName );
    writer->Update();
    }

  return 0;
}
