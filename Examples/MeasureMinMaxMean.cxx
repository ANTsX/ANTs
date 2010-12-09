/*=========================================================================
  Program:   Advanced Normalization Tools
  Module:    $RCSfile: MeasureMinMaxMean.cxx,v $
  Language:  C++
  Date:      $Date: 2008/12/16 17:56:34 $
  Version:   $Revision: 1.20 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// #include "DoSomethingToImage.cxx"

//  RecursiveAverageImages img1  img2  weight

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

// Note: could easily add variance computation
// http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Variance.htm

#include "ReadWriteImage.h"

template <unsigned int ImageDimension, unsigned int NVectorComponents>
int MeasureMinMaxMean(int argc, char *argv[])
{
  typedef itk::Vector<float, NVectorComponents>                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  typename ImageType::Pointer image = NULL;
  typename ImageType::Pointer mask = NULL;
  PixelType     mean;  mean.Fill(0);
  PixelType     max;   max.Fill(9.e9);
  PixelType     min;   min.Fill(-9.e9);
  unsigned long ct = 0;
  ReadImage<ImageType>(image, argv[2]);
  bool takeabsval = false;
  if( argc > 4 )
    {
    takeabsval = atoi(argv[4]);
    }
  if( argc > 5 )
    {
    ReadImage<ImageType>(mask, argv[5]);
    }
  Iterator vfIter2( image,  image->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool isinside = true;
    if( mask )
      {
      if( mask->GetPixel(vfIter2.GetIndex() ).GetNorm() < 0.5 )
        {
        isinside = false;
        }
      }
    if( isinside )
      {
      PixelType val = vfIter2.Get();
      if( takeabsval )
        {
        for( unsigned int k = 0; k < NVectorComponents; k++ )
          {
          val[k] = fabs(val[k]);
          }
        }
      mean += val;
      if( val.GetNorm() > max.GetNorm() )
        {
        max = val;
        }
      else if( val.GetNorm() < min.GetNorm() )
        {
        min = val;
        }
      ct++;
      }
    }
  if( ct > 0 )
    {
    mean /= (float)ct;
    }

  float variance = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool isinside = true;
    if( mask )
      {
      if( mask->GetPixel(vfIter2.GetIndex() ).GetNorm() < 0.5 )
        {
        isinside = false;
        }
      }
    if( isinside )
      {
      PixelType val = vfIter2.Get();
      if( takeabsval )
        {
        for( unsigned int k = 0; k < NVectorComponents; k++ )
          {
          val[k] = fabs(val[k]);
          }
        }
      variance += (val - mean).GetSquaredNorm();
      }
    }
  float temp = (1.0 / (float)ct) * variance;
  std::cout <<  argv[2] << " Max : " << max << " Min : " << min << " Mean : " << mean << " Var : " <<  temp
            << " SD : " << sqrt(1.0 / (float)(ct - 1) * variance) << std::endl;

  if( argc > 3 )
    {
    std::ofstream logfile;
    logfile.open(argv[3], std::ofstream::app);
    if( logfile.good() )
      {
      logfile << argv[2] <<  " Max : " << max << " Min : " << min << " Mean : " << mean << std::endl;
      }
    else
      {
      std::cout << " cant open file " << argv[3] << std::endl;
      }
    logfile.close();
    }

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Basic useage ex: " << std::endl;
    std::cout << argv[0] << " ImageDimension  image.nii {log.txt} {take-absolute-value}  {mask-name} " << std::endl;
    std::cout << "  log.txt is optional  - take-abs-val reports min-max-mean of abs val image " << std::endl;
    return 1;
    }

  int                       dim = atoi( argv[1] );
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(argv[2], itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(argv[2]);
  imageIO->ReadImageInformation();
  unsigned int ncomponents = imageIO->GetNumberOfComponents();
  std::cout << " ncomponents " << ncomponents << " dim " << imageIO->GetNumberOfDimensions() <<  std::endl;

  // Get the image dimension
  switch( atoi(argv[1]) )
    {
    case 2:
      {
      switch( ncomponents )
        {
        case 3:
          {
          MeasureMinMaxMean<2, 3>(argc, argv);
          }
          break;
        case 2:
          {
          MeasureMinMaxMean<2, 2>(argc, argv);
          }
          break;
        default:
          {
          MeasureMinMaxMean<2, 1>(argc, argv);
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
          MeasureMinMaxMean<3, 7>(argc, argv);
          }
          break;
        case 6:
          {
          MeasureMinMaxMean<3, 6>(argc, argv);
          }
          break;
        case 3:
          {
          MeasureMinMaxMean<3, 3>(argc, argv);
          }
          break;
        default:
          {
          MeasureMinMaxMean<3, 1>(argc, argv);
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
          MeasureMinMaxMean<4, 7>(argc, argv);
          }
          break;
        case 6:
          {
          MeasureMinMaxMean<4, 6>(argc, argv);
          }
          break;
        case 4:
          {
          MeasureMinMaxMean<4, 4>(argc, argv);
          }
          break;
        case 3:
          {
          MeasureMinMaxMean<4, 3>(argc, argv);
          }
          break;
        case 2:
          {
          MeasureMinMaxMean<4, 2>(argc, argv);
          }
          break;
        default:
          {
          MeasureMinMaxMean<4, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    default:
      std::cerr << " not supported " << dim  << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
