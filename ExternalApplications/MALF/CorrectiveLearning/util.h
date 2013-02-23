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

const unsigned int Dimension = 3;
typedef float                                             PixelType;
typedef itk::Image<PixelType, Dimension>                  ImageType;
typedef itk::ImageFileReader<ImageType>                   ReaderType;
typedef itk::ImageFileWriter<ImageType>                   WriterType;
typedef float                                             WritePixelType;
typedef itk::Image<WritePixelType, 3>                     WriteImageType;
typedef itk::ImageLinearConstIteratorWithIndex<ImageType> ConstIteratorType;
typedef itk::NeighborhoodIterator<ImageType>              NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex<ImageType>      IndexIteratorType;
typedef itk::ImageRegionIterator<ImageType>               IteratorType;

// perform thresholding and dilation to obtain region of interest for some segmentation label
// im: input image
// dmask: the resulting region of interest mask. Value 0 represents background, 1 is the ROI.
// Tlabel: specifies the label of interst. If Tlabel<0, then all non-background labels are used.
// R: specifies the dilation radius.
void myDilate(ImageType::Pointer im, ImageType::Pointer dmask, int Tlabel, int R)
{
  int                  x, y, z, j, k, l;
  ImageType::IndexType idx;
  ImageType::Pointer   tmask = ImageType::New();

  tmask->SetRegions(im->GetRequestedRegion() );
  tmask->SetSpacing(im->GetSpacing() );
  tmask->SetOrigin(im->GetOrigin() );
  tmask->SetDirection(im->GetDirection() );
  tmask->Allocate();

  int Dx = im->GetRequestedRegion().GetSize()[0];
  int Dy = im->GetRequestedRegion().GetSize()[1];
  int Dz = im->GetRequestedRegion().GetSize()[2];

  IndexIteratorType dmaskit(dmask, dmask->GetRequestedRegion() );
  IndexIteratorType tmaskit(tmask, tmask->GetRequestedRegion() );
  IndexIteratorType imit(im, im->GetRequestedRegion() );
  for( dmaskit.GoToBegin(), tmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit, ++tmaskit )
    {
    dmaskit.Set(0);
    tmaskit.Set(0);
    }
  for( imit.GoToBegin(); !imit.IsAtEnd(); ++imit )
    {
    idx = imit.GetIndex();
    if( (Tlabel >= 0 && imit.Value() == Tlabel) || (Tlabel < 0 && imit.Value() > 0) )
      {
      x = idx[0];
      for( j = -R; j < R + 1; j++ )
        {
        idx[0] = x + j;
        if( idx[0] < 0 || idx[0] >= Dx )
          {
          continue;
          }

        dmaskit.SetIndex(idx);
        dmaskit.Set(1);
        }
      }
    }
  for( dmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit )
    {
    idx = dmaskit.GetIndex();
    if( dmaskit.Value() == 1 )
      {
      y = idx[1];
      for( k = -R; k < R + 1; k++ )
        {
        idx[1] = y + k;
        if( idx[1] < 0 || idx[1] >= Dy )
          {
          continue;
          }

        tmaskit.SetIndex(idx);
        tmaskit.Set(1);
        }
      }
    }
  for( tmaskit.GoToBegin(); !tmaskit.IsAtEnd(); ++tmaskit )
    {
    idx = tmaskit.GetIndex();
    if( tmaskit.Value() == 1 )
      {
      z = idx[2];
      for( l = -R; l < R + 1; l++ )
        {
        idx[2] = z + l;
        if( idx[2] < 0 || idx[2] >= Dz )
          {
          continue;
          }

        dmaskit.SetIndex(idx);
        dmaskit.Set(1);
        }
      }
    }
}
