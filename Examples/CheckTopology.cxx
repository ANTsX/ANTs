/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>

#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "ReadWriteData.h"
#include "itkRandomImageSource.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "BinaryImageToMeshFilter.h"

#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkExtractEdges.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

namespace ants
{
float
random_range(float lowest_number, float highest_number)
{
  float range = highest_number - lowest_number;

  return lowest_number + (float)(range * (float)rand() / (float)(RAND_MAX));
}

float
ComputeGenus(vtkPolyData * pd1)
{
  vtkExtractEdges * edgeex = vtkExtractEdges::New();

  edgeex->SetInputData(pd1);
  edgeex->Update();
  vtkPolyData * edg1 = edgeex->GetOutput();
  vtkIdType     nedg = edg1->GetNumberOfCells();
  vtkIdType     vers = pd1->GetNumberOfPoints();
  int           nfac = pd1->GetNumberOfPolys();
  float         g = 0.5 * (2.0 - vers + nedg - nfac);
  std::cout << " Genus " << g << std::endl;
  std::cout << " face " << nfac << " edg " << nedg << " vert " << vers << std::endl;

  // edg1->Delete(); //caused malloc err
  edgeex->Delete(); // should be deleted b/c of New() above !!
  return g;
}

float
vtkComputeTopology(vtkPolyData * pd)
{
  vtkPolyDataConnectivityFilter * con = vtkPolyDataConnectivityFilter::New();

  con->SetExtractionModeToLargestRegion();
  con->SetInputData(pd);
  float g = ComputeGenus(con->GetOutput());
  con->Delete(); // should be deleted b/c of New() above !!
  return g;
}

template <typename TImage>
float
GetImageTopology(typename TImage::Pointer image)
{
  using ImageType = TImage;
  double aaParm = 0.024;
  using FilterType = BinaryImageToMeshFilter<ImageType>;
  typename FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput(image);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->SetAntiAliasMaxRMSError(-1000.0); // to do nothing
  fltMesh->SetSmoothingIterations(0);
  fltMesh->Update();
  vtkPolyData * vtkmesh = fltMesh->GetMesh();
  std::cout << " start topo " << std::endl;
  float genus = vtkComputeTopology(vtkmesh);
  std::cout << " Genus " << genus << std::endl;
  //  vtkmesh->Delete();
  return genus;
}

template <typename TImage>
void
NormalizeImage(typename TImage::Pointer image)
{
  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  float    max = 0;
  Iterator vfIter2(image, image->GetLargestPossibleRegion());
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    if (vfIter2.Get() > max)
    {
      max = vfIter2.Get();
    }
  }
  if (itk::Math::FloatAlmostEqual(max, itk::NumericTraits<float>::ZeroValue()))
  {
    max = itk::NumericTraits<float>::OneValue();
  }
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    vfIter2.Set(vfIter2.Get() / max);
  }
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, float sig)
{
  using ImageType = TImage;
  enum
  {
    ImageDimension = ImageType::ImageDimension
  };

  using dgf = itk::DiscreteGaussianImageFilter<ImageType, ImageType>;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(sig);
  filter->SetUseImageSpacing(false);
  filter->SetMaximumError(.01f);
  filter->SetInput(image);
  filter->Update();
  return filter->GetOutput();
}

template <typename TImage>
// std::vector<unsigned int>
typename TImage::Pointer
GetLargestComponent(typename TImage::Pointer image)
{
  enum
  {
    ImageDimension = TImage::ImageDimension
  };

  using InternalPixelType = int;
  using InternalImageType = itk::Image<InternalPixelType, ImageDimension>;
  using ThresholdFilterType = itk::BinaryThresholdImageFilter<TImage, InternalImageType>;
  using FilterType = itk::ConnectedComponentImageFilter<InternalImageType, InternalImageType>;
  using RelabelType = itk::RelabelComponentImageFilter<InternalImageType, InternalImageType>;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer          filter = FilterType::New();
  typename RelabelType::Pointer         relabel = RelabelType::New();

  threshold->SetInput(image);
  threshold->SetInsideValue(itk::NumericTraits<InternalPixelType>::OneValue());
  threshold->SetOutsideValue(itk::NumericTraits<InternalPixelType>::ZeroValue());
  threshold->SetLowerThreshold(0.499);
  threshold->SetUpperThreshold(1.001);
  threshold->Update();

  filter->SetInput(threshold->GetOutput());
  // if (argc > 5)
  {
    int fullyConnected = 1; // std::stoi( argv[5] );
    filter->SetFullyConnected(fullyConnected);
  }
  relabel->SetInput(filter->GetOutput());
  unsigned int minSize = 50;
  std::cout << " min Size " << minSize << std::endl;
  relabel->SetMinimumObjectSize(minSize);
  //    relabel->SetUseHistograms(true);

  try
  {
    relabel->Update();
  }
  catch (const itk::ExceptionObject & excep)
  {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
  }

  typename TImage::Pointer Clusters = AllocImage<TImage>(image, 0);

  using Iterator = itk::ImageRegionIteratorWithIndex<InternalImageType>;
  Iterator vfIter(relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion());

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
      Clusters->SetPixel(vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]);
      //  if ( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      //    maximgval=Clusters->GetPixel( vfIter.GetIndex());
    }
    else
    {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
    }
  }

  float maximgval = 0;
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (Clusters->GetPixel(vfIter.GetIndex()) > maximgval)
    {
      maximgval = Clusters->GetPixel(vfIter.GetIndex());
    }
  }

  std::cout << " max size " << maximgval << std::endl;
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    if (Clusters->GetPixel(vfIter.GetIndex()) >= maximgval)
    {
      Clusters->SetPixel(vfIter.GetIndex(), 1);
    }
    else
    {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
    }
  }

  //  for (int i=0; i<=maximum; i++)
  //  std::cout << " label " << i << " ct is: " << histogram[i] << std::endl;

  return Clusters;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
CheckTopology(std::vector<std::string> args, std::ostream *)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "CheckTopology");
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

  if (argc < 2)
  {
    std::cerr << "Parameter  missing" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << argv[0] << "  image.nii  {g0image.nii}  {threshold}" << std::endl;
    std::cerr << " If you put an arg for g0image then image will be smoothed and thresholded \n until it has genus "
                 "zero or the smoothing kernel gets too large "
              << std::endl;
    return EXIT_FAILURE;
  }

  float thresh = -1; // 0.0001;
  if (argc > 3)
  {
    thresh = atof(argv[3]);
  }
  using PixelType = float;
  constexpr unsigned int ImageDimension = 3; // AvantsImageDimension;
  using ImageType = itk::Image<PixelType, ImageDimension>;

  ImageType::Pointer image = ImageType::New();
  ReadImage<ImageType>(image, argv[1]);
  image = BinaryThreshold<ImageType>(0.5, 1.e9, 1, image);
  float initG = GetImageTopology<ImageType>(image);

  if (initG < 0 && argc > 2)
  {
    std::cout << "smoothing into a Genus Zero image with thresh " << thresh << std::endl;
    float              G = 1;
    float              smooth = 1;
    ImageType::Pointer simage;
    while (!itk::Math::FloatAlmostEqual(G, itk::NumericTraits<float>::ZeroValue()) && smooth < 20)
    {
      simage = SmoothImage<ImageType>(image, smooth);
      NormalizeImage<ImageType>(simage);
      simage = BinaryThreshold<ImageType>(thresh, 1.e9, 1, simage);
      ImageType::Pointer bigimage = GetLargestComponent<ImageType>(simage);
      G = GetImageTopology<ImageType>(bigimage);
      smooth = smooth + 1;
      simage = bigimage;
      std::cout << " G " << G << " at smoothing " << smooth << std::endl;
    }

    std::cout << " Final G " << G << " at smoothing " << smooth << std::endl;

    float        G2 = 0;
    unsigned int mct = 0;
    float        err = 1.e9;
    float        lasterr = 1.e10;
    float        derr = lasterr - err;
    while (itk::Math::FloatAlmostEqual(G2, itk::NumericTraits<float>::ZeroValue()) &&
           derr > itk::NumericTraits<float>::ZeroValue())
    {
      lasterr = err;
      err = 0;
      ImageType::Pointer out = ants::Morphological<ImageType>(simage, 3, 0, 1);
      ImageType::Pointer bigimage = GetLargestComponent<ImageType>(out);
      G2 = GetImageTopology<ImageType>(bigimage);
      using ImageIteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
      ImageIteratorType iter(bigimage, bigimage->GetLargestPossibleRegion());
      iter.GoToBegin();
      while (!iter.IsAtEnd())
      {
        err += fabs(iter.Get() - image->GetPixel(iter.GetIndex()));
        ++iter;
      }

      mct++;
      derr = lasterr - err;
      std::cout << " G2 " << G2 << " at morph " << mct << " err " << err << std::endl;
      if (itk::Math::FloatAlmostEqual(G2, itk::NumericTraits<float>::ZeroValue()) &&
          derr > itk::NumericTraits<float>::ZeroValue())
      {
        simage = GetLargestComponent<ImageType>(out);
      }
    }

    ANTs::WriteImage<ImageType>(simage, argv[2]);
  }
  else if (argc > 2)
  {
    ANTs::WriteImage<ImageType>(image, argv[2]);
  }
  return EXIT_SUCCESS;
}
} // namespace ants
