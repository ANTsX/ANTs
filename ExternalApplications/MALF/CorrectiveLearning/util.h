/*
This file implements some utility functions

Hongzhi Wang 07/19/2010
*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std;

namespace
{
const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType >  WriterType;
typedef float WritePixelType;
typedef itk::Image< WritePixelType, Dimension > WriteImageType;
typedef itk::ImageLinearConstIteratorWithIndex< ImageType >  ConstIteratorType;
typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
typedef itk::ImageRegionIterator< ImageType>        IteratorType;
typedef ImageType::RegionType RegionType;
typedef itk::Image<float, Dimension> PosteriorImage;
typedef PosteriorImage::Pointer PosteriorImagePtr;
typedef std::map<float, PosteriorImagePtr> PosteriorMap;
}

// perform thresholding and dilation to obtain region of interest for some segmentation label
// im: input image
// dmask: the resulting region of interest mask. Value 0 represents background, 1 is the ROI.  
// Tlabel: specifies the label of interst. If Tlabel<0, then all non-background labels are used.
// R: specifies the dilation radius.
void myDilate(ImageType::Pointer im, ImageType::Pointer dmask, int Tlabel, int R){
    int x, y, z, j, k, l;
    ImageType::IndexType idx;
    ImageType::Pointer tmask = ImageType::New();
    tmask->SetRegions(im->GetRequestedRegion());
    tmask->SetSpacing(im->GetSpacing() );
    tmask->SetOrigin(im->GetOrigin() );
    tmask->SetDirection(im->GetDirection());
    tmask->Allocate();

    int Dx=im->GetRequestedRegion().GetSize()[0];
    int Dy=im->GetRequestedRegion().GetSize()[1];
    int Dz=im->GetRequestedRegion().GetSize()[2];

    IndexIteratorType dmaskit(dmask, dmask->GetRequestedRegion());
    IndexIteratorType tmaskit(tmask, tmask->GetRequestedRegion());
    IndexIteratorType imit(im, im->GetRequestedRegion());
    for (dmaskit.GoToBegin(),tmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit,++tmaskit){
      dmaskit.Set(0);    
      tmaskit.Set(0);
    }
    int tc=0;
    for (imit.GoToBegin(); !imit.IsAtEnd(); ++imit){
      idx = imit.GetIndex();
      if ((Tlabel>=0 && imit.Value()==Tlabel) || (Tlabel<0 && imit.Value()>0)){
        tc++;
        x=idx[0];
        for (j=-R;j<R+1;j++){
          idx[0] = x+j;
          if (idx[0]<0 || idx[0]>=Dx)
            continue;

          dmaskit.SetIndex(idx);
          dmaskit.Set(1);
        }
      }
    }
//    cout<<"tc: "<<tc<<endl;
    for (dmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit){
      idx = dmaskit.GetIndex();
      if (dmaskit.Value()==1){
        y=idx[1];
        for (k=-R;k<R+1;k++){
          idx[1] = y+k;
          if (idx[1]<0 || idx[1]>=Dy)
            continue;

          tmaskit.SetIndex(idx);
          tmaskit.Set(1);
        }
      }
    }

    for (tmaskit.GoToBegin(); !tmaskit.IsAtEnd(); ++tmaskit){
      idx = tmaskit.GetIndex();
      if (tmaskit.Value()==1){
        z=idx[2];
        for (l=-R;l<R+1;l++){
          idx[2] = z+l;
          if (idx[2]<0 || idx[2]>=Dz)
            continue;

          dmaskit.SetIndex(idx);
          dmaskit.Set(1);
        }
      }
    }
}

// This function applies MRF optimization.
// ICM (Besag 1986, http://www.stat.duke.edu/~scs/Courses/Stat376/Papers/GibbsFieldEst/BesagDirtyPicsJRSSB1986.pdf)
// Originally implemented by Paul for the c3d software. Adopted by Hongzhi for segadapter.
void MRF_ICM(ImageType::Pointer nseg, PosteriorMap m_PosteriorMap, int ML, vector<int> labelset, double beta, int niter)
{
    // Allocate array for neighborhood histogram
    double *nhist = new double[ML+1];

    // Define an inner region (no boundary voxels)
    RegionType r_inner = nseg->GetBufferedRegion();
    for(size_t d = 0; d < 3; d++)
    {
      r_inner.SetIndex(d, r_inner.GetIndex(d)+1);
      r_inner.SetSize(d, r_inner.GetSize(d)-2);
    }

    // Create neighborhood iterator
    typedef itk::NeighborhoodIterator<ImageType> NeighborIterType;
    NeighborIterType::RadiusType radius;
    radius.Fill(1);
    NeighborIterType nit(radius, nseg, r_inner);
    nit.SetNeedToUseBoundaryCondition(false);

    // Do some iterations
    for(int q = 0; q < niter; q++)
    {
      cout<<"iter: "<<q<<endl;
      // Keep track of number of updates
      size_t n_upd = 0;

      // Iterate over the inner region
      IndexIteratorType rit(nseg, r_inner);

      for(rit.GoToBegin(); !rit.IsAtEnd(); ++rit)
      {
        // Current pixel value
        PixelType x_i = rit.Value();

        // Clear the neighborhood histogram
        for(int j = 0; j < ML; j++)
          nhist[j] = 0.0;

        // Make up for the fact that the current voxel will be counted
        nhist[(int)x_i] = -1;

        // Iterate the neighborhood
        nit.SetLocation(rit.GetIndex()); // Slow , change later
        for(size_t k = 0; k < nit.Size(); k++)
          nhist[(int)nit.GetPixel(k)]++;

        // For each candidate label value, compute conditional posterior
        int j_best = (int)x_i;
        double post_best = 1e100;
        for(size_t j = 0; j < labelset.size(); j++)
        {
          int label=labelset[j];
          double p = m_PosteriorMap[label]->GetPixel(rit.GetIndex());
          if(p > 0)
          {
            double likelihood = -log(p);
            double prior = - beta * nhist[label];
            double post = likelihood + prior;
            if(post_best > post)
            {
              post_best = post;
              j_best=label;
            }
          }
        }
        if(x_i != j_best)
        {
          rit.Set(j_best);
          n_upd++;
        }
      }
    cout<<n_upd<<" voxels modified."<<endl;
    if(n_upd == 0)
      {
      cout << "  Early convergence after " << q << " iterations" << endl;
      break;
      }
    }
    cout<<"MRF optimization is finished."<<endl;
    delete[] nhist;
}

