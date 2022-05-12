/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include <algorithm>
#include <list>
#include <vector>
#include <fstream>
#include "vnl/vnl_vector.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int
ClusterStatistics(unsigned int argc, char * argv[])
{
  using PixelType = float;
  //  const unsigned int ImageDimension = AvantsImageDimension;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  // typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  using ULPixelType = unsigned long;
  using labelimagetype = itk::Image<ULPixelType, ImageDimension>;
  using FilterType = itk::ConnectedComponentImageFilter<ImageType, labelimagetype>;
  using RelabelType = itk::RelabelComponentImageFilter<labelimagetype, labelimagetype>;

  // want the average value in each cluster as defined by the mask and the value thresh and the clust thresh

  std::string roimaskfn = std::string(argv[2]);
  std::string labelimagefn = std::string(argv[3]);
  std::string outname = std::string(argv[4]);
  float       clusterthresh = atof(argv[5]);
  float       minSize = clusterthresh;
  float       valuethresh = atof(argv[6]);
  //  std::cout << " Cth " << clusterthresh << " Vth " << valuethresh << std::endl;
  typename ImageType::Pointer valimage = nullptr;
  typename ImageType::Pointer roiimage = nullptr;
  typename ImageType::Pointer labelimage = nullptr;

  ReadImage<ImageType>(roiimage, roimaskfn.c_str());
  ReadImage<ImageType>(labelimage, labelimagefn.c_str());

  using MinMaxFilterType = itk::MinimumMaximumImageFilter<ImageType>;
  typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
  minMaxFilter->SetInput(labelimage);
  minMaxFilter->Update();
  double min = minMaxFilter->GetMinimum();
  double max = minMaxFilter->GetMaximum();
  double range = max - min;
  for (unsigned int filecount = 7; filecount < argc; filecount++)
  {
    //    std::cout <<" doing " << std::string(argv[filecount]) << std::endl;

    ReadImage<ImageType>(valimage, argv[filecount]);

    //  first, threshold the value image then get the clusters of min size
    using ThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
    typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
    threshold->SetInput(valimage);
    threshold->SetInsideValue(1);
    threshold->SetOutsideValue(0);
    threshold->SetLowerThreshold(valuethresh);
    threshold->SetUpperThreshold(1.e9);
    threshold->Update();
    typename ImageType::Pointer thresh = threshold->GetOutput();
    using fIterator = itk::ImageRegionIteratorWithIndex<ImageType>;
    using Iterator = itk::ImageRegionIteratorWithIndex<labelimagetype>;
    fIterator tIter(thresh, thresh->GetLargestPossibleRegion());
    for (tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter)
    {
      if (roiimage->GetPixel(tIter.GetIndex()) < static_cast<PixelType>(0.5))
      {
        tIter.Set(0);
      }
    }

    //  typename
    typename FilterType::Pointer filter = FilterType::New();
    // typename
    typename RelabelType::Pointer relabel = RelabelType::New();

    filter->SetInput(thresh);
    int fullyConnected = 0; // std::stoi( argv[5] );
    filter->SetFullyConnected(fullyConnected);
    relabel->SetInput(filter->GetOutput());
    relabel->SetMinimumObjectSize((unsigned int)minSize);

    try
    {
      relabel->Update();
    }
    catch (const itk::ExceptionObject & excep)
    {
      std::cerr << "Relabel: exception caught !" << std::endl;
      std::cerr << excep << std::endl;
    }

    typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(valimage, 0);
    typename ImageType::Pointer Values = MakeNewImage<ImageType>(valimage, 0);
    typename ImageType::Pointer Labels = MakeNewImage<ImageType>(valimage, 0);
    Iterator                    vfIter(relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());

    float maximum = relabel->GetNumberOfObjects();
    //    std::cout << " #object " << maximum << std::endl;
    //    float maxtstat=0;
    std::vector<unsigned long> histogram((int)maximum + 1);
    std::vector<long>          maxlabel((int)maximum + 1);
    std::vector<float>         suminlabel((unsigned long)range + 1);
    std::vector<unsigned long> countinlabel((unsigned long)range + 1);
    std::vector<float>         sumofvalues((int)maximum + 1);
    std::vector<float>         maxvalue((int)maximum + 1);
    for (int i = 0; i <= maximum; i++)
    {
      histogram[i] = 0;
      sumofvalues[i] = 0;
      maxvalue[i] = 0;
      maxlabel[i] = 0;
    }
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      if (vfIter.Get() > 0)
      {
        float vox = valimage->GetPixel(vfIter.GetIndex());
        if (vox >= valuethresh)
        {
          histogram[(unsigned long)vfIter.Get()] = histogram[(unsigned long)vfIter.Get()] + 1;
          sumofvalues[(unsigned long)vfIter.Get()] = sumofvalues[(unsigned long)vfIter.Get()] + vox;
          if (maxvalue[(unsigned long)vfIter.Get()] < vox)
          {
            maxvalue[(unsigned long)vfIter.Get()] = vox;
            maxlabel[(unsigned long)vfIter.Get()] = (long int)labelimage->GetPixel(vfIter.GetIndex());
          }

          suminlabel[(unsigned long)(labelimage->GetPixel(vfIter.GetIndex()) - static_cast<float>(min))] += vox;
          countinlabel[(unsigned long)(labelimage->GetPixel(vfIter.GetIndex()) - static_cast<float>(min))] += 1;
        }
      }
    }
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      if (vfIter.Get() > 0)
      {
        Clusters->SetPixel(vfIter.GetIndex(), histogram[(unsigned long)vfIter.Get()]);
        Values->SetPixel(vfIter.GetIndex(),
                         sumofvalues[(unsigned long)vfIter.Get()] / (float)histogram[(unsigned int)vfIter.Get()]);
        Labels->SetPixel(vfIter.GetIndex(), labelimage->GetPixel(vfIter.GetIndex()));
      }
      else
      {
        Clusters->SetPixel(vfIter.GetIndex(), 0);
        Labels->SetPixel(vfIter.GetIndex(), 0);
        Values->SetPixel(vfIter.GetIndex(), 0);
      }
    }

    //  ANTs::WriteImage<ImageType>(Values,std::string("temp.nii.gz").c_str());
    // ANTs::WriteImage<ImageType>(Clusters,std::string("temp2.nii.gz").c_str());

    float maximgval = 0;
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      if (Clusters->GetPixel(vfIter.GetIndex()) > maximgval)
      {
        maximgval = Clusters->GetPixel(vfIter.GetIndex());
      }
    }
    //  std::cout << " max size " << maximgval << std::endl;
    for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
    {
      if (Clusters->GetPixel(vfIter.GetIndex()) < minSize)
      {
        Clusters->SetPixel(vfIter.GetIndex(), 0);
        Values->SetPixel(vfIter.GetIndex(), 0);
        Labels->SetPixel(vfIter.GetIndex(), 0);
      }
    }

    //  ANTs::WriteImage<ImageType>(Values,(outname+"values.nii.gz").c_str());
    //  ANTs::WriteImage<ImageType>(Labels,(outname+"labels.nii.gz").c_str());
    ANTs::WriteImage<ImageType>(Clusters, (outname + "sizes.nii.gz").c_str());

    // now begin output
    //  std::cout << " Writing Text File " << outname << std::endl;
    std::string   outname2 = outname + std::string("average.csv");
    std::string   outname3 = outname + std::string("volume.csv");
    std::ofstream outf((outname2).c_str(), std::ofstream::out);
    std::ofstream outf2((outname3).c_str(), std::ofstream::out);
    if (outf.good())
    {
      //    outf << std::string(argv[filecount]) << std::endl;
      for (int i = 0; i < maximum + 1; i++)
      {
        if (histogram[i] >= minSize)
        {
          //          outf << " Cluster " << i << " size  " << histogram[i] <<  " average " <<
          // sumofvalues[i]/(float)histogram[i] << " max " << maxvalue[i] << " label " <<  maxlabel[i] <<  std::endl;
          std::cout << " Cluster " << i << " size  " << histogram[i] << " average "
                    << sumofvalues[i] / (float)histogram[i] << " max " << maxvalue[i] << " label " << maxlabel[i]
                    << std::endl;
        }
      }
      for (unsigned int i = 0; i <= range; i++)
      {
        //      if ( countinlabel[i] > 0)
        {
          if (countinlabel[i] == 0)
          {
            countinlabel[i] = 1;
          }
          //          outf << " Label " << i+min <<   " average " << suminlabel[i]/(float)countinlabel[i] <<  std::endl;
          std::cout << " Label " << i + min << " average " << suminlabel[i] / (float)countinlabel[i] << std::endl;
          if (i < range)
          {
            outf << suminlabel[i] / (float)countinlabel[i] << ",";
          }
          else
          {
            outf << suminlabel[i] / (float)countinlabel[i] << std::endl;
          }
        }
      }
    }
    else
    {
      std::cout << " File No Good! " << outname << std::endl;
    }
    outf.close();

    if (outf2.good())
    {
      for (unsigned int i = 0; i <= range; i++)
      {
        if (countinlabel[i] == 0)
        {
          countinlabel[i] = 1;
        }
        if (i < range)
        {
          outf2 << (float)countinlabel[i] << ",";
        }
        else
        {
          outf2 << (float)countinlabel[i] << std::endl;
        }
      }
    }
    else
    {
      std::cout << " File No Good! " << outname << std::endl;
    }
    outf2.close();
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ClusterImageStatistics(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ClusterImageStatistics");
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
  argv[argc] = nullptr;
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

  if (argc < 4)
  {
    std::cout << " Given an ROI and Label Image, find the max and average value   \n in a value image  where the value "
                 "> some user-defined threshold \n and the cluster size  is larger than some min size. \n "
              << std::endl;
    std::cout << "Usage: \n  " << std::endl;
    std::cout << argv[0]
              << "  ImageDimension ROIMask.ext LabelImage.ext  OutPrefix   MinimumClusterSize  ValueImageThreshold  "
                 "Image1WithValuesOfInterest.ext ...  ImageNWithValuesOfInterest.ext  \n \n "
              << std::endl;
    std::cout << " ROIMask.ext -- overall region of interest \n  \n LabelImage.ext -- labels for the sub-regions, e.g. "
                 "Brodmann or just unique labels (see  LabelClustersUniquely ) \n \n  OutputPrefix -- all output  has "
                 "this prefix  \n \n  MinimumClusterSize -- the minimum size of clusters of interest  \n  \n "
                 "ValueImageThreshold -- minimum value of interest \n \n   Image*WithValuesOfInterest.ext  ---  "
                 "image(s) that define the values you want to measure \n ";
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      ClusterStatistics<2>(argc, argv);
    }
    break;
    case 3:
    {
      ClusterStatistics<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
