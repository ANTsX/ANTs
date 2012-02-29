/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: MeasureImageSimilarity.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/05 20:09:47 $
  Version:   $Revision: 1.19 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "ReadWriteImage.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkSpatialMutualInformationRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkCrossCorrelationRegistrationFunction.h"
// #include "itkLandmarkCrossCorrelationRegistrationFunction.h"

template <unsigned int ImageDimension>
int MeasureImageSimilarity(unsigned int argc, char *argv[])
{
  typedef float                                                  PixelType;
  typedef itk::Vector<float, ImageDimension>                     VectorType;
  typedef itk::Image<VectorType, ImageDimension>                 FieldType;
  typedef itk::Image<PixelType, ImageDimension>                  ImageType;
  typedef itk::ImageFileWriter<ImageType>                        writertype;
  typedef typename ImageType::IndexType                          IndexType;
  typedef typename ImageType::SizeType                           SizeType;
  typedef typename ImageType::SpacingType                        SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>           AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType1;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>           Iterator;

  typedef itk::Image<float, 2>                JointHistType;
  typedef itk::ImageFileWriter<JointHistType> jhwritertype;

// get command line params
  unsigned int argct = 2;
  unsigned int whichmetric = atoi(argv[argct]); argct++;
  std::string  fn1 = std::string(argv[argct]); argct++;
  std::string  fn2 = std::string(argv[argct]); argct++;
  std::string  logfilename = "";
  if( argc > argct )
    {
    logfilename = std::string(argv[argct]);
    }
  argct++;
  std::string imgfilename = "";
  if( argc > argct )
    {
    imgfilename = std::string(argv[argct]);
    }
  argct++;
  double targetvalue = 0;
  if( argc > argct )
    {
    targetvalue = atof(argv[argct]);
    }
  argct++;
  double epsilontolerance = 1.e20;
  if( argc > argct )
    {
    epsilontolerance = atof(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::Pointer image2 = NULL;
  ReadImage<ImageType>(image2, fn2.c_str() );

/*
  typedef itk::ImageRegionIteratorWithIndex<FieldType> VIterator;
  typename FieldType::Pointer field=FieldType::New();
  field->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
  field->SetBufferedRegion( image1->GetLargestPossibleRegion() );
  field->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
  field->Allocate();
  field->SetSpacing(image1->GetSpacing());
  field->SetOrigin(image1->GetOrigin());
  VectorType zero;
  zero.Fill(0);
  VIterator vfIter2( field,  field->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
      IndexType index=vfIter2.GetIndex();
      vfIter2.Set(zero);
    }
*/
  typename ImageType::Pointer metricimg = ImageType::New();
  metricimg->SetRegions( image1->GetLargestPossibleRegion() );
  metricimg->Allocate();
  metricimg->SetSpacing(image1->GetSpacing() );
  metricimg->SetOrigin(image1->GetOrigin() );
  metricimg->SetDirection(image1->GetDirection() );
  Iterator iter( metricimg,  metricimg->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    iter.Set(0);
    }

  typedef ImageType FixedImageType;
  typedef ImageType MovingImageType;
  typedef FieldType DisplacementFieldType;

  // Choose the similarity metric
  typedef itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType,
                                                           DisplacementFieldType>  MIMetricType;
  typedef itk::SpatialMutualInformationRegistrationFunction<FixedImageType, MovingImageType,
                                                            DisplacementFieldType> SMIMetricType;
  typedef itk::CrossCorrelationRegistrationFunction<FixedImageType, MovingImageType,
                                                    DisplacementFieldType>         CCMetricType;
  // typedef itk::LandmarkCrossCorrelationRegistrationFunction<FixedImageType,MovingImageType,DisplacementFieldType>
  // MetricType;
  // typename
  typename MIMetricType::Pointer mimet = MIMetricType::New();
  typename SMIMetricType::Pointer smimet = SMIMetricType::New();
  typename CCMetricType::Pointer ccmet = CCMetricType::New();

//  int nbins=32;

  typename CCMetricType::RadiusType hradius;
  typename CCMetricType::RadiusType ccradius;
  ccradius.Fill(4);
  typename MIMetricType::RadiusType miradius;
  miradius.Fill(0);

//  mimet->SetDisplacementField(field);
  mimet->SetFixedImage(image1);
  mimet->SetMovingImage(image2);
  mimet->SetRadius(miradius);
  mimet->SetGradientStep(1.e2);
  mimet->SetNormalizeGradient(false);

//  ccmet->SetDisplacementField(field);
  ccmet->SetFixedImage(image1);
  ccmet->SetMovingImage(image2);
  ccmet->SetRadius(ccradius);
  ccmet->SetGradientStep(1.e2);
  ccmet->SetNormalizeGradient(false);

  double      metricvalue = 0;
  std::string metricname = "";
  if( whichmetric  == 0 )
    {
    hradius = miradius;
    unsigned long ct = 0;
    for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      IndexType index = iter.GetIndex();
      double    fval = image1->GetPixel(index);
      double    mval = image2->GetPixel(index);
      metricvalue += fabs(fval - mval);
      metricimg->SetPixel(index, fabs(fval - mval) );
      ct++;
      }
    metricvalue /= (float)ct;
    metricname = "MSQ ";
    }
  else if( whichmetric == 1 ) // imagedifference
    {
    double ccval = 0;
    hradius = ccradius;
    metricname = "CC ";
    typedef itk::ConstNeighborhoodIterator<FixedImageType> ScanIteratorType;
    typename FixedImageType::RegionType region = image1->GetLargestPossibleRegion();
    ScanIteratorType asamIt( hradius, image1, region);
    unsigned long    ct = 0;
    for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      IndexType index = iter.GetIndex();
      double    val = 0;
      double    fmip = 0, mmip = 0, ffip = 0;
      asamIt.SetLocation(index);
      for( unsigned int i = 0; i < asamIt.Size(); i++ )
        {
        IndexType locind = asamIt.GetIndex(i);
        if( region.IsInside( locind ) )
          {
          double f = image1->GetPixel(locind);
          double m = image2->GetPixel(locind);
          fmip += (f * m);  mmip += (m * m);  ffip += (f * f);
          }
        }
      double denom = mmip * ffip;
      if( denom == 0 )
        {
        val = 1;
        }
      else
        {
        val = fmip / sqrt(denom);
        }
      metricimg->SetPixel(index, val);
      // if (ct % 10000 == 0)
      //        std::cout << val << " index " << index << std::endl;
      //      asamIt.SetLocation(index);
      //      totval+=met->localProbabilistic;
      ccval += val;
      ct++;
      }
    if( ct >  0 )
      {
      metricvalue = ccval / ct;
      }
    else
      {
      metricvalue = 0;
      }
    /*
    std::cout << metricname << std::endl;
    ccmet->InitializeIteration();
    metricimg=ccmet->MakeImage();
    metricvalue=ccmet->ComputeCrossCorrelation()*(1.0);
    */
    }
  else if( whichmetric == 2 )
    {
    hradius = miradius;
    mimet->InitializeIteration();
    metricvalue = mimet->ComputeMutualInformation();
    metricname = "MI ";
    }
  else if( whichmetric == 3 )
    {
    hradius = miradius;
    smimet->InitializeIteration();
    metricvalue = smimet->ComputeSpatialMutualInformation();
    metricname = "SMI ";
    }
  std::cout << fn1 << " : " << fn2 << " => " <<  metricname << metricvalue << std::endl;
  if( logfilename.length() > 3 )
    {
    std::ofstream logfile;
    logfile.open(logfilename.c_str(), std::ofstream::app);
    if( logfile.good() )
      {
      logfile <<  fn1 << " : " << fn2 << " => " << metricname << metricvalue << std::endl;
      }
    else
      {
      std::cout << " cant open file ";
      }
    logfile.close();
    }

  if( imgfilename.length() > 3 )
    {
    std::cout << "Only Implemented for MSQ and CC " << std::endl;
    typedef itk::ImageFileWriter<ImageType> writertype;
    typename writertype::Pointer w = writertype::New();
    w->SetInput(metricimg);
    w->SetFileName(imgfilename.c_str() );
    w->Write(); //  met->WriteImages();

/*
  typename MetricType::NeighborhoodType   asamIt( hradius, field,field->GetLargestPossibleRegion());
  unsigned long ct = 0;
  double totval = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
      IndexType index=vfIter2.GetIndex();
      double val=0;
      asamIt.SetLocation(index);
      //      met->ComputeUpdate( asamIt,  gd);
      met->ComputeMetricAtPairB(index,  zero);
      metricimg->SetPixel(index, val);
      //if (ct % 10000 == 0)
      //        std::cout << val << " index " << index << std::endl;
      //      asamIt.SetLocation(index);
      //      totval+=met->localProbabilistic;
      ct++;
    }

  std::cout << " AvantsMI : " << totval/(double)ct << " E " <<  met->GetEnergy() <<  std::endl;
  std::cout << " write begin " << std::endl;

  std::cout << " write end " << std::endl;
*/
    }

  double diff = ( (double)metricvalue - (double) targetvalue);
  std::cout << " targetvalue " << targetvalue << " metricvalue " << metricvalue << " diff " << diff << " toler "
            << epsilontolerance << std::endl;

  if( diff < epsilontolerance )
    {
    return EXIT_SUCCESS;
    }
  else
    {
    return EXIT_FAILURE;
    }
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << "Basic useage ex: " << std::endl;
    std::cout << argv[0]
              <<
    " ImageDimension whichmetric image1.ext image2.ext {logfile} {outimage.ext}  {target-value}   {epsilon-tolerance}"
              << std::endl;
    std::cout << "  outimage (Not Implemented for MI yet)  and logfile are optional  " << std::endl;
    std::cout
      <<
    " target-value and epsilon-tolerance set goals for the metric value -- if the metric value is within epsilon-tolerance of the target-value, then the test succeeds "
      << std::endl;
    std::cout << "  Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information , 3-SMI " << std::endl;
    return 1;
    }

  int metricsuccess = EXIT_FAILURE;

  // Get the image dimension
  switch( atoi(argv[1]) )
    {
    case 2:
      {
      metricsuccess = MeasureImageSimilarity<2>(argc, argv);
      }
      break;
    case 3:
      {
      metricsuccess = MeasureImageSimilarity<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  std::cout << " Failure? " << metricsuccess << std::endl;

  return metricsuccess;
}
