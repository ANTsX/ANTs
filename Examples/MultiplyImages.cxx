/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: MultiplyImages.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include "ReadWriteImage.h"

template <unsigned int ImageDimension>
int MultiplyImages(int argc, char *argv[])
{
  typedef float PixelType;
  //  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  std::string fn1 = std::string(argv[2]);
  std::string fn2 = std::string(argv[3]);
  std::string outname = std::string(argv[4]);

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::Pointer varimage = NULL;

  typename readertype::Pointer reader2 = readertype::New();
  typename readertype::Pointer reader1 = readertype::New();
  reader2->SetFileName(fn2.c_str() );

  bool isfloat = false;
  try
    {
    reader2->UpdateLargestPossibleRegion();
    }
  catch( ... )
    {
    std::cout << " Rather than opening " << fn2
              <<
    " as an image file, this program has decided, in its great wisdom, to consider it to be a floating point numerical value, and has acted accordingly -- i.e. read this as a number. "
              << std::endl;
    isfloat = true;
    }

  float floatval = 1.0;
  if( isfloat )
    {
    floatval = atof(argv[3]);
    }
  else
    {
    image2 = reader2->GetOutput();
    }

  reader1->SetFileName(fn1.c_str() );
  try
    {
    reader1->UpdateLargestPossibleRegion();
    image1 = reader1->GetOutput();
    }
  catch( ... )
    {
    std::cout << " read 1 error ";
    }

  varimage = ImageType::New();
  varimage->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
  varimage->SetBufferedRegion( image1->GetLargestPossibleRegion() );
  varimage->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
  varimage->Allocate();
  varimage->SetSpacing(image1->GetSpacing() );
  varimage->SetOrigin(image1->GetOrigin() );
  varimage->SetDirection(image1->GetDirection() );
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    IndexType ind = vfIter2.GetIndex();
    float     pix2;
    if( isfloat )
      {
      pix2 = floatval;
      }
    else
      {
      pix2 = image2->GetPixel(ind);
      }
    float pix1 = image1->GetPixel(ind);
    // std::cout << " p1 " << pix1 << " p2 " << pix2 << std::endl;
    float div = pix1 * pix2;
    vfIter2.Set(div);
    }

  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[4]);
  writer->SetInput( varimage );
  writer->Write();

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << "Usage:  " << std::endl;
    std::cout << argv[0] << " ImageDimension img1.nii img2.nii product.nii {smoothing}" << std::endl;
    return 1;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      MultiplyImages<2>(argc, argv);
      }
      break;
    case 3:
      {
      MultiplyImages<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
