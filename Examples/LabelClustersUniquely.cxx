/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: LabelClustersUniquely.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include <list>
#include <vector>
#include <fstream>
#include "vnl/vnl_vector.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkCastImageFilter.h"
#include  "ReadWriteImage.h"

template <unsigned int ImageDimension>
int  LabelUniquely(int argc, char *argv[])
{
  typedef float PixelType;
//  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  // typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef float                                           InternalPixelType;
  typedef int                                             ULPixelType;
  typedef itk::Image<ULPixelType, ImageDimension>         labelimagetype;
  typedef itk::CastImageFilter<ImageType, labelimagetype> CastFilterType;

  typedef ImageType                                                          InternalImageType;
  typedef ImageType                                                          OutputImageType;
  typedef itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype> FilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, labelimagetype>   RelabelType;

  // want the average value in each cluster as defined by the mask and the value thresh and the clust thresh

  std::string fn1 = std::string(argv[1]);
  std::string fn2 = std::string(argv[2]);
  float       clusterthresh = atof(argv[3]);

  typename ImageType::Pointer image1 = NULL;

  ReadImage<ImageType>(image1, fn1.c_str() );

  //  typename
  typename FilterType::Pointer filter = FilterType::New();
// typename
  typename RelabelType::Pointer relabel = RelabelType::New();

  typename CastFilterType::Pointer castInput = CastFilterType::New();
  castInput->SetInput(image1);

  filter->SetInput( castInput->GetOutput() );
  int fullyConnected = 0; // atoi( argv[5] );
  filter->SetFullyConnected( fullyConnected );
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( (unsigned int) clusterthresh );

  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

//  float maximum=relabel->GetNumberOfObjects();
  WriteImage<labelimagetype>( relabel->GetOutput(), argv[2]);

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Useage ex:  " << std::endl;
    std::cout << argv[0] << " ImageDimension clustersin.hdr labeledclustersout.hdr   sizethresh " << std::endl;
    return 1;
    }

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      LabelUniquely<2>(argc, argv + 1);
      }
      break;
    case 3:
      {
      LabelUniquely<3>(argc, argv + 1);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
