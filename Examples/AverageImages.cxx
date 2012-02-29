/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: AverageImages.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/27 23:25:24 $
  Version:   $Revision: 1.21 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

// Note: could easily add variance computation
// http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Variance.htm

#include "itkArray.h"
#include "itkVariableLengthVector.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkOptimalSharpeningImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"

template <unsigned int ImageDimension, unsigned int NVectorComponents>
int AverageImages1(unsigned int argc, char *argv[])
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageFileReader<ImageType>              ImageFileReader;
  typedef itk::ImageFileWriter<ImageType>              writertype;

  bool  normalizei = atoi(argv[3]);
  float numberofimages = (float)argc - 4.;
  typename ImageType::Pointer averageimage = NULL;
  typename ImageType::Pointer image2 = NULL;

  typename ImageType::SizeType size;
  size.Fill(0);
  unsigned int bigimage = 0;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    std::cout << " fn " << fn << std::endl;
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        bigimage = j;
        std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  std::cout << " largest image " << size << std::endl;
  typename ImageFileReader::Pointer reader = ImageFileReader::New();
  reader->SetFileName(argv[bigimage]);
  reader->Update();
  averageimage = reader->GetOutput();
  unsigned int vectorlength = reader->GetImageIO()->GetNumberOfComponents();
  std::cout << " Averaging " << numberofimages << " images with dim = " << ImageDimension << " vector components "
            << vectorlength << std::endl;
  PixelType meanval = 0;
  averageimage->FillBuffer(meanval);
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::cout << " reading " << std::string(argv[j]) << std::endl;
    typename ImageFileReader::Pointer rdr = ImageFileReader::New();
    rdr->SetFileName(argv[j]);
    rdr->Update();
    image2 = rdr->GetOutput();
    Iterator      vfIter2( image2,  image2->GetLargestPossibleRegion() );
    unsigned long ct = 0;
    if( normalizei )
      {
      meanval = 0;
      for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
        {
        PixelType localp = image2->GetPixel( vfIter2.GetIndex() );
        meanval = meanval + localp;
        ct++;
        }
      if( ct > 0 )
        {
        meanval = meanval / (float)ct;
        }
      if( meanval <= 0 )
        {
        meanval = (1);
        }
      }
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      PixelType val = vfIter2.Get();
      if( normalizei )
        {
        val /= meanval;
        }
      val = val / (float)numberofimages;
      PixelType oldval = averageimage->GetPixel(vfIter2.GetIndex() );
      averageimage->SetPixel(vfIter2.GetIndex(), val + oldval );
      }
    }

  //  typedef itk::OptimalSharpeningImageFilter<ImageType,ImageType > sharpeningFilter;
  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> sharpeningFilter;
  typename sharpeningFilter::Pointer shFilter = sharpeningFilter::New();
  if( normalizei && argc > 3 && vectorlength == 1 )
    {
    shFilter->SetInput( averageimage );
    //    shFilter->SetSValue(0.5);
    averageimage =  shFilter->GetOutput();
    }

  std::cout << " writing output ";
    {
    typename writertype::Pointer writer = writertype::New();
    writer->SetFileName(argv[2]);
    writer->SetInput( averageimage );
    writer->Update();
    }

  return 0;
}

template <unsigned int ImageDimension, unsigned int NVectorComponents>
int AverageImages(unsigned int argc, char *argv[])
{
  typedef itk::Vector<float, NVectorComponents>        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageFileReader<ImageType>              ImageFileReader;
  typedef itk::ImageFileWriter<ImageType>              writertype;

  bool  normalizei = atoi(argv[3]);
  float numberofimages = (float)argc - 4.;
  typename ImageType::Pointer averageimage = NULL;
  typename ImageType::Pointer image2 = NULL;

  typename ImageType::SizeType size;
  size.Fill(0);
  unsigned int bigimage = 0;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    std::cout << " fn " << fn << std::endl;
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        bigimage = j;
        std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  std::cout << " largest image " << size << std::endl;
  typename ImageFileReader::Pointer reader = ImageFileReader::New();
  reader->SetFileName(argv[bigimage]);
  reader->Update();
  averageimage = reader->GetOutput();
  unsigned int vectorlength = reader->GetImageIO()->GetNumberOfComponents();
  std::cout << " Averaging " << numberofimages << " images with dim = " << ImageDimension << " vector components "
            << vectorlength << std::endl;
  typename ImageType::IndexType zindex; zindex.Fill(0);
  PixelType meanval = reader->GetOutput()->GetPixel(zindex);
  meanval.Fill(0);
  averageimage->FillBuffer(meanval);
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::cout << " reading " << std::string(argv[j]) << std::endl;
    typename ImageFileReader::Pointer rdr = ImageFileReader::New();
    rdr->SetFileName(argv[j]);
    rdr->Update();
    image2 = rdr->GetOutput();
    Iterator      vfIter2( image2,  image2->GetLargestPossibleRegion() );
    unsigned long ct = 0;
    if( normalizei )
      {
      meanval.Fill(0);
      for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
        {
        PixelType localp = image2->GetPixel( vfIter2.GetIndex() );
        meanval = meanval + localp;
        ct++;
        }
      if( ct > 0 )
        {
        meanval = meanval / (float)ct;
        }
      if( meanval.GetNorm() <= 0 )
        {
        meanval.Fill(1);
        }
      }
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      PixelType val = vfIter2.Get();
      if( normalizei )
        {
        for( unsigned int k = 0; k < vectorlength; k++ )
          {
          val[k] /= meanval[k];
          }
        }
      val = val / (float)numberofimages;
      PixelType oldval = averageimage->GetPixel(vfIter2.GetIndex() );
      averageimage->SetPixel(vfIter2.GetIndex(), val + oldval );
      }
    }

    {
    typename writertype::Pointer writer = writertype::New();
    writer->SetFileName(argv[2]);
    writer->SetInput( averageimage );
    writer->Update();
    }

  return 0;
}

int main(int argc, char * argv[])
{
  if( argc < 3 )
    {
    std::cout << "\n" << std::endl;
    std::cout << "Usage: \n" << std::endl;
    std::cout << argv[0] << " ImageDimension Outputfname.nii.gz Normalize <images> \n" << std::endl;
    std::cout << " Compulsory arguments: \n" << std::endl;
    std::cout << " ImageDimension: 2 or 3 (for 2 or 3 dimensional input).\n " << std::endl;
    std::cout << " Outputfname.nii.gz: the name of the resulting image.\n" << std::endl;
    std::cout
      <<
    " Normalize: 0 (false) or 1 (true); if true, the 2nd image is divided by its mean. This will select the largest image to average into.\n"
      << std::endl;
    std::cout << " Example Usage:\n" << std::endl;
    std::cout << argv[0] << " 3 average.nii.gz  1  *.nii.gz \n" << std::endl;
    std::cout << " \n" << std::endl;
    return 1;
    }

  int                       dim = atoi( argv[1] );
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(argv[4], itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(argv[4]);
  imageIO->ReadImageInformation();
  unsigned int ncomponents = imageIO->GetNumberOfComponents();

  // Get the image dimension
  switch( atoi(argv[1]) )
    {
    case 2:
      {
      switch( ncomponents )
        {
        case 2:
          {
          AverageImages<2, 2>(argc, argv);
          }
          break;
        default:
          {
          AverageImages1<2, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 3:
      {
      switch( ncomponents )
        {
        case 7:
          {
          AverageImages<3, 7>(argc, argv);
          }
          break;
        case 6:
          {
          AverageImages<3, 6>(argc, argv);
          }
          break;
        case 3:
          {
          AverageImages<3, 3>(argc, argv);
          }
          break;
        case 2:
          {
          AverageImages<3, 2>(argc, argv);
          }
          break;
        default:
          {
          AverageImages1<3, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 4:
      {
      switch( ncomponents )
        {
        case 7:
          {
          AverageImages<4, 7>(argc, argv);
          }
          break;
        case 6:
          {
          AverageImages<4, 6>(argc, argv);
          }
          break;
        case 4:
          {
          AverageImages<4, 4>(argc, argv);
          }
          break;
        case 3:
          {
          AverageImages<4, 3>(argc, argv);
          }
          break;
        default:
          {
          AverageImages1<4, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    default:
      std::cerr << " You passed ImageDimension: " << dim << " . Please use only image domains of 2, 3 or 4  "
                << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
