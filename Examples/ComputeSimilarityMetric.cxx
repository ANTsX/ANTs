/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>
#include <algorithm>
#include "ReadWriteData.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkCrossCorrelationRegistrationFunction.h"

namespace ants
{
template <unsigned int ImageDimension>
int
ComputeSimilarityMetric(int argc, char * argv[])
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
  unsigned int whichmetric = std::stoi(argv[argct]);
  argct++;
  std::string fn1 = std::string(argv[argct]);
  argct++;
  std::string fn2 = std::string(argv[argct]);
  argct++;
  std::string logfilename = "";
  if (argc > argct)
  {
    logfilename = std::string(argv[argct]);
  }
  argct++;
  std::string imgfilename = "";
  if (argc > argct)
  {
    imgfilename = std::string(argv[argct]);
  }
  argct++;

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str());
  typename ImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image2, fn2.c_str());

  /*
    typedef itk::ImageRegionIteratorWithIndex<FieldType> VIterator;
    typename FieldType::Pointer field=FieldType::New();
    field->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
    field->SetBufferedRegion( image1->GetLargestPossibleRegion() );
    field->SetLargestPossibleRegion( image1->GetLargestPossibleRegion() );
    field->AllocateInitialized();
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
  typedef ImageType FixedImageType;
  typedef ImageType MovingImageType;
  typedef FieldType DisplacementFieldType;

  // Choose the similarity metric
  typedef itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType>
    MIMetricType;
  typedef itk::CrossCorrelationRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType>
    CCMetricType;
  // typedef itk::LandmarkCrossCorrelationRegistrationFunction<FixedImageType,MovingImageType,DisplacementFieldType>
  // MetricType;
  // typename
  typename MIMetricType::Pointer mimet = MIMetricType::New();
  typename CCMetricType::Pointer ccmet = CCMetricType::New();

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

  // ORIENTATION ALERT  -- the original code here
  // set the region, spacing, and origin without setting directions.
  typename ImageType::Pointer metricimg = AllocImage<ImageType>(image1, 0);

  if (whichmetric == 0)
  {
    hradius = miradius;
    unsigned long ct = 0;
    Iterator      iter(metricimg, metricimg->GetLargestPossibleRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
    {
      IndexType index = iter.GetIndex();
      double    fval = image1->GetPixel(index);
      double    mval = image2->GetPixel(index);
      metricvalue += fabs(fval - mval);
      ct++;
    }
    metricvalue /= (float)ct;
    metricname = "MSQ ";
  }
  else if (whichmetric == 1) // imagedifference
  {
    hradius = ccradius;
    ccmet->InitializeIteration();
    metricvalue = ccmet->ComputeCrossCorrelation();
    metricname = "CC ";
  }
  else
  {
    hradius = miradius;
    mimet->InitializeIteration();
    metricvalue = mimet->ComputeMutualInformation();
    metricname = "MI ";
  }
  std::cout << fn1 << " : " << fn2 << " => " << metricname << metricvalue << std::endl;
  if (logfilename.length() > 3)
  {
    std::ofstream logfile;
    logfile.open(logfilename.c_str(), std::ofstream::app);
    if (logfile.good())
    {
      logfile << fn1 << " : " << fn2 << " => " << metricname << metricvalue << std::endl;
    }
    else
    {
      std::cout << " cant open file ";
    }
    logfile.close();
  }

  if (imgfilename.length() > 3)
  {
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
      typedef itk::ImageFileWriter<ImageType> writertype;
      writertype::Pointer w= writertype::New();
      w->SetInput(metricimg);
      w->SetFileName(outname.c_str());
      w->Write();  //  met->WriteImages();

      std::cout << " write end " << std::endl;
    */
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ComputeSimilarityMetric(std::vector<std::string> args, std::ostream * out_stream = nullptr)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ComputeSimilarityMetric");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = 0;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  if (argc < 3)
  {
    std::cout << "Basic useage ex: " << std::endl;
    std::cout << argv[0] << " ImageDimension whichmetric image1.ext image2.ext {logfile} {outimage.ext}  " << std::endl;
    std::cout << "  outimage and logfile are optional  " << std::endl;
    return EXIT_FAILURE;
  }

  // Get the image dimension
  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      ComputeSimilarityMetric<2>(argc, argv);
    }
    break;
    case 3:
    {
      ComputeSimilarityMetric<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
} // namespace ants
