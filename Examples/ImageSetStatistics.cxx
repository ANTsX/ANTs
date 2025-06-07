/*=========================================================================

  Program:   Advanced Normalization Tools

  Copriyght (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "ReadWriteData.h"

#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkNeighborhoodIterator.h"

//  RecursiveAverageImages img1  img2  weight

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

// Note: could easily add variance computation
// http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Variance.htm
#include "itkDiscreteGaussianImageFilter.h"

namespace ants
{
template <typename TImageType>
void
ReadImage(itk::SmartPointer<TImageType> & target, const char * file, bool copy)
{
  //  std::cout << " reading b " << std::string(file) << std::endl;
  using readertype = itk::ImageFileReader<TImageType>;
  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(file);
  reader->Update();
  if (!copy)
  {
    target = (reader->GetOutput());
  }
  else
  {
    using Iterator = itk::ImageRegionIteratorWithIndex<TImageType>;
    Iterator vfIter2(target, target->GetLargestPossibleRegion());
    for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
    {
      vfIter2.Set(reader->GetOutput()->GetPixel(vfIter2.GetIndex()));
    }
  }
}

double
TProb(double t, int df)
{
  if (itk::Math::FloatAlmostEqual(t, 0.0))
  {
    return 0;
  }
  double a = 0.36338023;
  double w = atan(t / sqrt((double)df));
  double s = sin(w);
  double c = cos(w);

  double t1, t2;
  int    j1, j2, k2;

  if (df % 2 == 0) // even
  {
    t1 = s;
    if (df == 2) // special case df=2
    {
      return 0.5 * (1 + t1);
    }
    t2 = s;
    j1 = -1;
    j2 = 0;
    k2 = (df - 2) / 2;
  }
  else
  {
    t1 = w;
    if (df == 1) // special case df=1
    {
      return 1 - (0.5 * (1 + (t1 * (1 - a))));
    }
    t2 = s * c;
    t1 = t1 + t2;
    if (df == 3) // special case df=3
    {
      return 1 - (0.5 * (1 + (t1 * (1 - a))));
    }
    j1 = 0;
    j2 = 1;
    k2 = (df - 3) / 2;
  }
  for (int i = 1; i >= k2; i++)
  {
    j1 = j1 + 2;
    j2 = j2 + 2;
    t2 = t2 * c * c * j1 / j2;
    t1 = t1 + t2;
  }
  return 1 - (0.5 * (1 + (t1 * (1 - a * (df % 2)))));
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, float sig)
{
  using dgf = itk::DiscreteGaussianImageFilter<TImage, TImage>;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(sig);
  filter->SetUseImageSpacing(true);
  filter->SetMaximumError(.01f);
  filter->SetInput(image);
  filter->Update();
  return filter->GetOutput();
}

template <typename TInputImage>
// typename TInputImage::Pointer
void
HistogramMatch(typename TInputImage::Pointer m_InputFixedImage, typename TInputImage::Pointer m_InputMovingImage)
{
  std::cout << " MATCHING INTENSITIES " << std::endl;

  using FilterType = itk::HistogramMatchingImageFilter<TInputImage, TInputImage>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(m_InputMovingImage);
  filter->SetReferenceImage(m_InputFixedImage);
  filter->SetNumberOfHistogramLevels(256);
  filter->SetNumberOfMatchPoints(10);
  filter->ThresholdAtMeanIntensityOn();
  filter->ThresholdAtMeanIntensityOff();
  filter->Update();
  typename TInputImage::Pointer img = filter->GetOutput();

  using Iterator = itk::ImageRegionIteratorWithIndex<TInputImage>;
  Iterator vfIter(img, img->GetLargestPossibleRegion());
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    m_InputMovingImage->SetPixel(vfIter.GetIndex(), vfIter.Get());
  }
}

template <typename TImage>
void
LocalMean(typename TImage::Pointer image, unsigned int nhood, typename TImage::Pointer meanimage)
{
  typename TImage::Pointer localmean = MakeNewImage<TImage>(image, 0);

  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator                  outIter(image, image->GetLargestPossibleRegion());
  typename TImage::SizeType imagesize = image->GetLargestPossibleRegion().GetSize();

  using iteratorType = itk::NeighborhoodIterator<TImage>;
  typename iteratorType::RadiusType rad;
  rad.Fill(static_cast<itk::SizeValueType>(nhood));

  using IndexType = typename TImage::IndexType;

  for (outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter)
  {
    itk::NeighborhoodIterator<TImage> hoodIt(rad, image, image->GetLargestPossibleRegion());
    typename TImage::IndexType        oindex = outIter.GetIndex();
    hoodIt.SetLocation(oindex);

    double fixedMean = 0;
    // double movingMean=0;

    unsigned int hoodlen = hoodIt.Size();

    // unsigned int inct=0;

    bool takesample = true;

    // double sumj=0;
    // double sumi=0;
    if (takesample)
    {
      // double sumj=0;
      double       sumi = 0;
      unsigned int cter = 0;
      for (unsigned int indct = 0; indct < hoodlen; indct++)
      {
        typename TImage::IndexType index = hoodIt.GetIndex(indct);
        bool                       inimage = true;
        for (itk::SizeValueType dd = 0; dd < TImage::ImageDimension; dd++)
        {
          if (index[dd] < itk::NumericTraits<typename IndexType::IndexValueType>::ZeroValue() ||
              index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd] - 1))
          {
            inimage = false;
          }
        }

        if (inimage)
        {
          sumi += static_cast<double>(image->GetPixel(index));
          cter++;
        }
      }

      if (cter > 0)
      {
        fixedMean = sumi / static_cast<double>(cter);
      }
    }

    float val = image->GetPixel(oindex) - static_cast<float>(fixedMean);
    meanimage->SetPixel(oindex, meanimage->GetPixel(oindex) + static_cast<float>(fixedMean));
    localmean->SetPixel(oindex, val);
  }

  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator vfIter(image, image->GetLargestPossibleRegion());
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    vfIter.Set(localmean->GetPixel(vfIter.GetIndex()));
  }

  // localmean;
}

template <typename TImage>
// std::vector<unsigned int>
float
GetClusterStat(typename TImage::Pointer image,
               float                    Tthreshold,
               unsigned int             minSize,
               unsigned int             whichstat,
               std::string              outfn,
               bool                     TRUTH)
{
  using InternalPixelType = float;
  using InternalImageType = TImage;
  using OutputImageType = TImage;
  using ThresholdFilterType = itk::BinaryThresholdImageFilter<InternalImageType, InternalImageType>;
  using FilterType = itk::ConnectedComponentImageFilter<InternalImageType, OutputImageType>;
  using RelabelType = itk::RelabelComponentImageFilter<OutputImageType, OutputImageType>;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer          filter = FilterType::New();
  typename RelabelType::Pointer         relabel = RelabelType::New();

  InternalPixelType threshold_low, threshold_hi;
  threshold_low = Tthreshold;
  threshold_hi = 1.e9;

  threshold->SetInput(image);
  threshold->SetInsideValue(itk::NumericTraits<InternalPixelType>::OneValue());
  threshold->SetOutsideValue(itk::NumericTraits<InternalPixelType>::ZeroValue());
  threshold->SetLowerThreshold(threshold_low);
  threshold->SetUpperThreshold(threshold_hi);
  threshold->Update();

  filter->SetInput(threshold->GetOutput());
  // if (argc > 5)
  {
    int fullyConnected = 1; // std::stoi( argv[5] );
    filter->SetFullyConnected(fullyConnected);
  }
  relabel->SetInput(filter->GetOutput());
  relabel->SetMinimumObjectSize(minSize);
  //    relabel->SetUseHistograms(true);

  try
  {
    relabel->Update();
  }
  catch (const itk::ExceptionObject & excep)
  {
    std::cout << "Relabel: exception caught !" << std::endl;
    std::cout << excep << std::endl;
  }

  typename TImage::Pointer Clusters = MakeNewImage<TImage>(relabel->GetOutput(), 0);
  // typename TImage::Pointer Clusters=relabel->GetOutput();
  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator vfIter(relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());

  /*

  typename itk::LabelStatisticsImageFilter<TImage,TImage>::Pointer labstat=
    itk::LabelStatisticsImageFilter<TImage,TImage>::New();
  labstat->SetInput(image);
  labstat->SetLabelImage(Clusters);
  labstat->SetUseHistograms(true);
  labstat->Update();

  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator vfIter( Clusters,   Clusters->GetLargestPossibleRegion() );

  float maximum=0;
  // Relabel the Clusters image with the right statistic
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
      if (relabel->GetOutput()->GetPixel(vfIter.GetIndex()) > 0 )
    {
      float pix = relabel->GetOutput()->GetPixel(vfIter.GetIndex());
      float val;
      if (whichstat == 0) val=pix;
      else if (whichstat == 1) val = labstat->GetSum(pix);
      else if (whichstat == 2) val = labstat->GetMean(pix);
      else if (whichstat == 3) val = labstat->GetMaximum(pix);
      if (val > maximum) maximum=val;
      vfIter.Set(val);
    }
    }
  */
  float                     maximum = relabel->GetNumberOfObjects();
  float                     maxtstat = 0;
  std::vector<unsigned int> histogram((int)maximum + 1);
  std::vector<float>        clustersum((int)maximum + 1);
  for (int i = 0; i <= maximum; i++)
  {
    histogram[i] = 0;
    clustersum[i] = 0;
  }
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (vfIter.Get() > 0)
    {
      float vox = image->GetPixel(vfIter.GetIndex());
      histogram[(unsigned int)vfIter.Get()] = histogram[(unsigned int)vfIter.Get()] + 1;
      clustersum[(unsigned int)vfIter.Get()] += vox;
      if (vox > maxtstat)
      {
        maxtstat = vox;
      }
    }
  }
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (vfIter.Get() > 0)
    {
      if (whichstat == 0) // size
      {
        Clusters->SetPixel(vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]);
      }
      if (whichstat == 1) // sum
      {
        Clusters->SetPixel(vfIter.GetIndex(), clustersum[(unsigned int)vfIter.Get()]);
      }
      if (whichstat == 2) // mean
      {
        Clusters->SetPixel(vfIter.GetIndex(),
                           clustersum[(unsigned int)vfIter.Get()] / (float)histogram[(unsigned int)vfIter.Get()]);
      }
      if (whichstat == 3) // max
      {
        Clusters->SetPixel(vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]);
      }
    }
    else
    {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
    }
  }

  //  for (int i=0; i<=maximum; i++)
  //  std::cout << " label " << i << " ct is: " << histogram[i] << std::endl;

  if (TRUTH)
  {
    using writertype = itk::ImageFileWriter<InternalImageType>;
    typename writertype::Pointer writer = writertype::New();
    writer->SetFileName((outfn + std::string("Clusters.nii")).c_str());
    writer->SetInput(Clusters);
    writer->Write();
  }

  if (whichstat == 0)
  {
    return histogram[1];
  }
  else if (whichstat == 1)
  {
    float mx = 0;
    for (int i = 1; i <= maximum; i++)
    {
      if (clustersum[i] > mx)
      {
        mx = clustersum[i];
      }
    }
    return mx;
  }
  else if (whichstat == 2)
  {
    float mx = 0;
    for (int i = 1; i <= maximum; i++)
    {
      if (clustersum[i] / (float)histogram[i] > mx)
      {
        mx = clustersum[i] / (float)histogram[i] * 1000.0f;
      }
    }
    return mx;
  }
  else if (whichstat == 3)
  {
    return maxtstat * 1000.0f;
  }
  else
  {
    return histogram[1];
  }
}

float
median(std::vector<float> vec)
{
  using vec_sz = std::vector<float>::size_type;
  vec_sz size = vec.size();

  if (size == 0)
  {
    return 0;
  }
  //            throw domain_error("median of an empty vector");

  sort(vec.begin(), vec.end());

  vec_sz mid = size / 2;

  return size % 2 == 0 ? (vec[mid] + vec[mid - 1]) / 2 : vec[mid];
}

float
npdf(std::vector<float> vec, bool opt, float www)
{
  using vec_sz = std::vector<float>::size_type;
  vec_sz size = vec.size();

  if (size == 0)
  {
    return 0;
  }
  //            throw domain_error("median of an empty vector");

  float mean = 0, var = 0;
  float max = -1.e9, min = 1.e9;
  for (unsigned int i = 0; i < size; i++)
  {
    float val = vec[i];
    if (val > max)
    {
      max = val;
    }
    else if (val < min)
    {
      min = val;
    }
    auto  n = (float)(i + 1);
    float wt1 = 1.0f / (float)n;
    float wt2 = 1.0f - wt1;
    mean = mean * wt2 + val * wt1;
    if (i > 0)
    {
      float wt3 = 1.0f / ((float)n - 1.0f);
      var = var * wt2 + (val - mean) * (val - mean) * wt3;
    }
  }

  if (itk::Math::FloatAlmostEqual(var, itk::NumericTraits<float>::ZeroValue()))
  {
    return mean;
  }
  //    else std::cout << " Mean " << mean << " var " << var << std::endl;

  // eval parzen probability
  std::vector<float> prob(size);
  float              maxprob = 0;
  //        float maxprobval=0;
  float        weightedmean = 0;
  float        weighttotal = 0;
  unsigned int maxprobind = 0;
  //        float sample=0.0;
  float width;
  if (www > 0)
  {
    width = www;
  }
  else
  {
    width = sqrt(var) / 2.0;
  }
  //        std::cout << " using width " << width << std::endl;
  //        float N=(float)size;
  for (unsigned int j = 0; j < size; j++)
  {
    float sample = vec[j];
    float total = 0.0;
    for (unsigned int i = 0; i < size; i++)
    {
      float delt = vec[i] - sample;
      delt *= delt;
      prob[i] = 1.0f / (2.0f * 3.1214f * width) * static_cast<float>(exp(-0.5f * delt / (width * width)));
      total += prob[i];
      //            maxprobval+=prob[i]
    }
    if (total > maxprob)
    {
      maxprob = total;
      maxprobind = j;
    }

    weightedmean += sample * total;
    weighttotal += total;
    //        for (unsigned int i=0; i<size; i++) prob[i]=prob[i]/total;
  }

  weightedmean /= weighttotal;
  // pxa = 1./N * total ( gaussian )
  maxprob = vec[maxprobind];
  if (opt)
  {
    return maxprob;
  }
  else
  {
    return weightedmean; // vec[maxprobind];
  }
}

float
trimmean(std::vector<float> vec)
{
  using vec_sz = std::vector<float>::size_type;
  vec_sz size = vec.size();

  if (size == 0)
  {
    return 0;
  }
  //            throw domain_error("median of an empty vector");

  sort(vec.begin(), vec.end());

  constexpr unsigned int lo = 0;
  const unsigned int     hi = size;
  const unsigned int     ct = hi - lo;
  float                  total = 0;
  for (unsigned int i = lo; i < hi; i++)
  {
    total += vec[i];
  }
  return total / (float)ct;
}

float
myantsmax(std::vector<float> vec)
{
  using vec_sz = std::vector<float>::size_type;
  vec_sz size = vec.size();
  if (size == 0)
  {
    return 0;
  }

  float max = -1.e9;
  for (unsigned int i = 0; i < size; i++)
  {
    float val = vec[i];
    if (val > max)
    {
      max = val;
    }
  }
  return max;
}

float
myantssimilaritymaxlabel(std::vector<float> labelvec, std::vector<float> similarityvec, bool opt)
{
  using vec_sz = std::vector<float>::size_type;
  vec_sz size = labelvec.size();
  if (size == 0)
  {
    return 0;
  }

  unsigned int max = 0;
  float        maxsim = -1.e9;
  float        totalsim = 0;
  for (unsigned int i = 0; i < size; i++)
  {
    totalsim += similarityvec[i];
  }
  if (fabs(totalsim) <= 0)
  {
    return 0;
  }
  for (unsigned int i = 0; i < size; i++)
  {
    float simval = similarityvec[i];
    if (simval > maxsim)
    {
      maxsim = simval;
      max = i;
    }
  }

  if (opt == true)
  {
    return labelvec[max];
  }
  else
  {
    return max;
  }
}

template <unsigned int ImageDimension>
int
ImageSetStatistics(int argc, char * argv[])
{
  using PixelType = float;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using readertype = itk::ImageFileReader<ImageType>;
  using IndexType = typename ImageType::IndexType;
  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;
  unsigned int mch = 0;
  int          argct = 2;
  std::string  fn1 = std::string(argv[argct]);
  argct++;
  std::string outfn = std::string(argv[argct]);
  argct++;
  unsigned int whichstat = std::stoi(argv[argct]);
  argct++;
  std::string roifn = "";
  if (argc > argct)
  {
    roifn = std::string(argv[argct]);
    argct++;
  }
  std::string simimagelist = std::string("");
  if (argc > argct)
  {
    simimagelist = std::string(argv[argct]);
    argct++;
  }
  float www = 0;
  // if (argc > argct) { www=atof(argv[argct]);argct++;}
  //  unsigned int mchmax= 0;
  // if (argc > argct) { mchmax=atoi(argv[argct]); argct++;}
  unsigned int localmeanrad = 0;
  // if (argc > argct) { localmeanrad=atoi(argv[argct]);argct++;}

  //  std::cout <<" roifn " << roifn << " fn1 " << fn1 << " whichstat " << whichstat << std::endl;

  typename ImageType::Pointer outimage = nullptr;
  typename ImageType::Pointer ROIimg = nullptr;

  if (roifn.length() > 4)
  {
    std::cout << " reading roi image " << roifn << std::endl;
    typename readertype::Pointer reader2 = readertype::New();
    reader2->SetFileName(roifn.c_str());
    reader2->UpdateLargestPossibleRegion();
    try
    {
      ROIimg = reader2->GetOutput();
    }
    catch (...)
    {
      ROIimg = nullptr;
      std::cout << " Error reading ROI image " << std::endl;
      //  return 0;
    }
  }

  // now do the recursive average
  constexpr unsigned int maxChar = 512;
  char                   lineBuffer[maxChar];
  char                   filenm[maxChar];
  unsigned int           filecount1 = 0;
  {
    std::ifstream inputStreamA(fn1.c_str(), std::ios::in);
    if (!inputStreamA.is_open())
    {
      std::cout << "Can't open parameter file: " << fn1 << std::endl;
      return EXIT_FAILURE;
    }
    while (!inputStreamA.eof())
    {
      inputStreamA.getline(lineBuffer, maxChar, '\n');

      if (sscanf(lineBuffer, "%s ", filenm) != 1)
      {
        //      std::cout << "Done.  read " << lineBuffer << " n " << ct1 << " files " << std::endl;
        // std::cout << std::endl;
        continue;
      }
      else
      {
        filecount1++;
      }
    }

    inputStreamA.close();
  }
  std::cout << " NFiles1 " << filecount1 << std::endl;

  unsigned int filecount2 = 0;
  if (simimagelist.length() > 2 && (whichstat == 5 || whichstat == 6))
  {
    std::ifstream inputStreamA(simimagelist.c_str(), std::ios::in);
    if (!inputStreamA.is_open())
    {
      std::cout << "Can't open parameter file: " << fn1 << std::endl;
      return EXIT_FAILURE;
    }
    while (!inputStreamA.eof())
    {
      inputStreamA.getline(lineBuffer, maxChar, '\n');

      if (sscanf(lineBuffer, "%s ", filenm) != 1)
      {
        //      std::cout << "Done.  read " << lineBuffer << " n " << ct1 << " files " << std::endl;
        // std::cout << std::endl;
        continue;
      }
      else
      {
        filecount2++;
      }
    }

    inputStreamA.close();
    if (filecount1 != filecount2)
    {
      std::cout << " the number of similarity images does not match the number of label images --- thus, we have to "
                   "get out of here !! i.e. something's wrong. "
                << std::endl;
      return EXIT_FAILURE;
    }
  } // fi simimagelist
  std::cout << " NFiles2 " << filecount2 << std::endl;

  typename ImageType::Pointer              meanimage;
  std::vector<typename ImageType::Pointer> imagestack;
  imagestack.resize(filecount1);
  //  imagestack.fill(nullptr);
  std::vector<std::string>    filenames(filecount1);
  typename ImageType::Pointer StatImage;
  unsigned int                ct = 0;
  std::ifstream               inputStreamA(fn1.c_str(), std::ios::in);
  if (!inputStreamA.is_open())
  {
    std::cout << "Can't open parameter file: " << fn1 << std::endl;
    return EXIT_FAILURE;
  }
  while (!inputStreamA.eof())
  {
    inputStreamA.getline(lineBuffer, maxChar, '\n');

    if (sscanf(lineBuffer, "%s ", filenm) != 1)
    {
      //      std::cout << "Done.  read " << lineBuffer << " n " << ct1 << " files " << std::endl;
      // std::cout << std::endl;
      continue;
    }
    else
    {
      filenames[ct] = std::string(filenm);
      ReadImage<ImageType>(imagestack[ct], filenm, false);
      if (ct == 0)
      {
        meanimage = MakeNewImage<ImageType>(imagestack[ct], 0);
      }
      if (localmeanrad > 0)
      {
        LocalMean<ImageType>(imagestack[ct], localmeanrad, meanimage);
      }
      std::cout << " done reading " << (float)ct / (float)filecount1 << std::endl;
      ct++;
    }
  }

  inputStreamA.close();

  // read similarity images, if needed
  std::vector<typename ImageType::Pointer> simimagestack;
  simimagestack.resize(filecount2);
  ct = 0;
  if (simimagelist.length() > 2 && (whichstat == 5 || whichstat == 6))
  {
    inputStreamA.open(simimagelist.c_str());
    if (!inputStreamA.is_open())
    {
      std::cout << "Can't open parameter file: " << fn1 << std::endl;
      return -1;
    }
    while (!inputStreamA.eof())
    {
      inputStreamA.getline(lineBuffer, maxChar, '\n');
      if (sscanf(lineBuffer, "%s ", filenm) != 1)
      {
        continue;
      }
      else
      {
        ReadImage<ImageType>(simimagestack[ct], filenm, false);
        ct++;
      }
    }

    inputStreamA.close();
  } // fi read similarity images

  ReadImage<ImageType>(StatImage, filenames[0].c_str(), false);
  Iterator           vfIter(StatImage, StatImage->GetLargestPossibleRegion());
  std::vector<float> voxels(filecount1);
  std::vector<float> similarities(filecount2);
  unsigned long      nvox = 1;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    nvox *= StatImage->GetLargestPossibleRegion().GetSize()[i];
  }

  ct = 0;
  unsigned long prog = nvox / 15;
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (ct % prog == 0)
    {
      std::cout << " % " << (float)ct / (float)nvox << std::endl;
    }
    ct++;
    IndexType    ind = vfIter.GetIndex();
    unsigned int maxval = 0;
    bool         takesample = true;
    if (ROIimg)
    {
      if (ROIimg->GetPixel(ind) < 0.5f)
      {
        takesample = false;
      }
      else
      {
        maxval = (unsigned int)(ROIimg->GetPixel(ind) - 1);
      }
    }
    if (takesample)
    {
      if (mch == 0)
      {
        meanimage->SetPixel(ind, meanimage->GetPixel(ind) / filecount1);
      }
      for (unsigned int j = 0; j < filecount1; j++)
      {
        voxels[j] = imagestack[j]->GetPixel(ind);
      }
      for (unsigned int j = 0; j < filecount2; j++)
      {
        similarities[j] = simimagestack[j]->GetPixel(ind);
      }
      float stat = 0;

      switch (whichstat)
      {
        case 1:
        {
          stat = npdf(voxels, true, www);
          if (ct == 1)
          {
            std::cout << "the max prob appearance \n";
          }
        }
        break;
        case 2:
        {
          stat = npdf(voxels, false, www);
          if (ct == 1)
          {
            std::cout << "the probabilistically weighted appearance " << www << " \n";
          }
        }
        break;

        case 3:
        {
          stat = trimmean(voxels);
          if (ct == 1)
          {
            std::cout << "the trimmed mean appearance \n";
          }
        }
        break;

        case 4:
        {
          stat = myantsmax(voxels);
          if (ct == 1)
          {
            std::cout << "the maximum appearance \n";
          }
        }
        break;
        case 5:
        {
          stat = myantssimilaritymaxlabel(voxels, similarities, true);
          if (ct == 1)
          {
            std::cout << "the maximum similarity-based label \n";
          }
        }
        break;
        case 6:
        {
          stat = myantssimilaritymaxlabel(voxels, similarities, false);
          if (ct == 1)
          {
            std::cout << "which image provides the maximum similarity-based label \n";
          }
        }
        break;
        case 7:
        {
          stat = voxels[maxval];
          if (ct == 1)
          {
            std::cout << "which image provides the maximum similarity-based label \n";
          }
        }
        break;

        default:
        {
          stat = median(voxels);
          if (ct == 1)
          {
            std::cout << "the median appearance \n";
          }
        }
        break;
      }
      float sval = stat;
      if (localmeanrad > 0)
      {
        sval += meanimage->GetPixel(ind);
      }
      StatImage->SetPixel(ind, sval);
    }
    else
    {
      StatImage->SetPixel(ind, 0);
    }
  }
  ANTs::WriteImage<ImageType>(StatImage, outfn.c_str());

  std::cout << " Done " << std::endl;
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ImageSetStatistics(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ImageSetStatistics");

  const int argc = args.size();
  char **   argv = new char *[args.size() + 1];
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
    std::cout << "Usage:  " << std::endl;
    std::cout
      << argv[0]
      << " ImageDimension controlslist.txt outimage.nii whichstat {roi.nii} {imagelist2forsimilarityweightedstats.txt}"
      << std::endl;
    std::cout << " whichstat = 0:  median,  1:  max prob appearance  , 2: weighted mean appearance ,  3: trimmed mean "
                 ", 4 : max value , option 5 : similarity-weighted (must pass imagelist2 as well) else median , option "
                 "6 : same as similarity-weighted option 5 but the label corresponds to the image that provides the "
                 "best local match ... useful if you want to MRF smooth these indices  , option 7 : similar to 5 but "
                 "expects the max-value to be stored in the ROI image and uses it to get the intensity ... "
              << std::endl;
    std::cout << " example:   ImageSetStatistics  3   imagelist.txt  maxvalueimage.nii.gz 4 " << std::endl;
    std::cout << " similarity weighted --- pass in a list of similarity images here which will be used to select the "
                 "best label --- thus, number of similarity images must match the number of label images . "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  // Get the image dimension

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return ImageSetStatistics<2>(argc, argv);
    }
    break;
    case 3:
    {
      return ImageSetStatistics<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
