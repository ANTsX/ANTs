#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"
#include "vnl/algo/vnl_determinant.h"

#include "ReadWriteData.h"

namespace ants
{
template <typename TField, typename TImage>
typename TImage::Pointer
GetVectorComponent(typename TField::Pointer field, unsigned int index)
{
  // Initialize the Moving to the displacement field
  typedef TImage ImageType;

  typename ImageType::Pointer sfield =
    AllocImage<ImageType>(field);

  typedef itk::ImageRegionIteratorWithIndex<TField> Iterator;
  Iterator vfIter( field,  field->GetLargestPossibleRegion() );
  for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    typename TField::PixelType v1 = vfIter.Get();
    sfield->SetPixel(vfIter.GetIndex(), v1[index]);
    }

  return sfield;
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, float sig)
{
// find min value
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator vfIter(image, image->GetLargestPossibleRegion() );
  for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    typename TImage::PixelType v1 = vfIter.Get();
    if( std::isnan(v1) )
      {
      vfIter.Set(0);
      }
    }
  typedef itk::DiscreteGaussianImageFilter<TImage, TImage> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(sig);
  filter->SetUseImageSpacingOn();
  filter->SetMaximumError(.01f);
  filter->SetInput(image);
  filter->Update();
  typename TImage::Pointer out = filter->GetOutput();

  return out;
}

template <typename TImage>
void
SmoothDeformation(typename TImage::Pointer vectorimage, float sig)
{
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<float, 3>  ImageType;
  typename ImageType::Pointer subimgx = GetVectorComponent<TImage, ImageType>(vectorimage, 0);
  subimgx = SmoothImage<ImageType>(subimgx, sig);
  typename ImageType::Pointer subimgy = GetVectorComponent<TImage, ImageType>(vectorimage, 1);
  subimgy = SmoothImage<ImageType>(subimgy, sig);
  typename ImageType::Pointer subimgz = GetVectorComponent<TImage, ImageType>(vectorimage, 2);
  subimgz = SmoothImage<ImageType>(subimgz, sig);

  typedef itk::ImageRegionIteratorWithIndex<TImage> IteratorType;
  IteratorType Iterator( vectorimage, vectorimage->GetLargestPossibleRegion().GetSize() );
  Iterator.GoToBegin();
  while(  !Iterator.IsAtEnd()  )
    {
    VectorType vec;
    vec[0] = subimgx->GetPixel(Iterator.GetIndex() );
    vec[1] = subimgy->GetPixel(Iterator.GetIndex() );
    vec[2] = subimgz->GetPixel(Iterator.GetIndex() );
    Iterator.Set(vec);
    ++Iterator;
    }

  return;
}

template <typename TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground,
             typename TImage::PixelType newval, typename TImage::Pointer input, float distthresh )
{
  std::cout << " Label Surf " << std::endl;
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  // ORIENTATION ALERT -- original code set spacing & origin without
  // also setting orientation.
  typename   ImageType::Pointer     Image =
    AllocImage<ImageType>(input);

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;

  typename iteratorType::RadiusType rad;
  for( int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = static_cast<unsigned int>( distthresh + 0.5f );
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

//  std::cout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    typename TImage::IndexType ind2;
    if( itk::Math::FloatAlmostEqual( p, foreground ) )
      {
      bool atedge = false;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        ind2 = GHood.GetIndex(i);
        float dist = 0.0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          dist += (float)(ind[j] - ind2[j]) * (float)(ind[j] - ind2[j]);
          }
        dist = sqrt(dist);
        if( ! itk::Math::FloatAlmostEqual( GHood.GetPixel(i), foreground ) && dist <  distthresh  )
          {
          atedge = true;
          }
        }
      if( atedge && itk::Math::FloatAlmostEqual( p, foreground ) )
        {
        Image->SetPixel(ind, newval);
        }
      else
        {
        Image->SetPixel(ind, 0);
        }
      }
    ++GHood;
    }

  return Image;
}

template <typename TImage>
typename TImage::Pointer  Morphological( typename TImage::Pointer input, float rad, bool option)
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::PixelType PixelType;

  if( !option )
    {
    std::cout << " eroding the image " << std::endl;
    }
  else
    {
    std::cout << " dilating the image " << std::endl;
    }
  typedef itk::BinaryBallStructuringElement<
      PixelType,
      ImageDimension>             StructuringElementType;

  // Software Guide : BeginCodeSnippet
  typedef itk::BinaryErodeImageFilter<
      TImage,
      TImage,
      StructuringElementType>  ErodeFilterType;

  typedef itk::BinaryDilateImageFilter<
      TImage,
      TImage,
      StructuringElementType>  DilateFilterType;

  typename ErodeFilterType::Pointer  binaryErode  = ErodeFilterType::New();
  typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

  StructuringElementType structuringElement;

  structuringElement.SetRadius( (unsigned long) rad );  // 3x3x3 structuring element

  structuringElement.CreateStructuringElement();

  binaryErode->SetKernel(  structuringElement );
  binaryDilate->SetKernel( structuringElement );

  //  It is necessary to define what could be considered objects on the binary
  //  images. This is specified with the methods \code{SetErodeValue()} and
  //  \code{SetDilateValue()}. The value passed to these methods will be
  //  considered the value over which the dilation and erosion rules will apply
  binaryErode->SetErodeValue( 1 );
  binaryDilate->SetDilateValue( 1 );

  typename TImage::Pointer temp;
  if( option )
    {
    binaryDilate->SetInput( input );
    binaryDilate->Update();
    temp = binaryDilate->GetOutput();
    }
  else
    {
    binaryErode->SetInput( input );  // binaryDilate->GetOutput() );
    binaryErode->Update();
    temp = binaryErode->GetOutput();

    typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;
    ImageIteratorType o_iter( temp, temp->GetLargestPossibleRegion() );
    o_iter.GoToBegin();
    while( !o_iter.IsAtEnd() )
      {
      if( o_iter.Get() > 0.5f && input->GetPixel(o_iter.GetIndex() ) > 0.5f )
        {
        o_iter.Set(1);
        }
      else
        {
        o_iter.Set(0);
        }
      ++o_iter;
      }
    }

  return temp;
}

template <typename TImage, typename TField>
typename TField::Pointer
FMMGrad(typename TImage::Pointer wm, typename TImage::Pointer gm )
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  typename TField::Pointer sfield = AllocImage<TField>(wm);

  typename ImageType::Pointer surf = LabelSurface<ImageType>(1, 1, wm, 1.9 );

  typedef itk::FastMarchingUpwindGradientImageFilter<ImageType, ImageType> FloatFMType;
  typename FloatFMType::Pointer marcher = FloatFMType::New();
  typedef typename FloatFMType::NodeType      NodeType;
  typedef typename FloatFMType::NodeContainer NodeContainer;
  // setup alive points
  typename NodeContainer::Pointer alivePoints = NodeContainer::New();
  typename NodeContainer::Pointer targetPoints = NodeContainer::New();
  typename NodeContainer::Pointer trialPoints = NodeContainer::New();
  typedef itk::ImageRegionIteratorWithIndex<TImage> IteratorType;
  IteratorType thIt( wm, wm->GetLargestPossibleRegion().GetSize() );
  thIt.GoToBegin();
  unsigned long bb = 0, cc = 0, dd = 0;
  while(  !thIt.IsAtEnd()  )
    {
    if( thIt.Get() > 0.1 && surf->GetPixel(thIt.GetIndex() ) == 0 )
      {
      NodeType node;
      node.SetValue( 0 );
      node.SetIndex(thIt.GetIndex() );
      alivePoints->InsertElement(bb, node);
      bb++;
      }
    if( gm->GetPixel(thIt.GetIndex() ) == 0 && wm->GetPixel(thIt.GetIndex() ) == 0  )
      {
      NodeType node;
      node.SetValue( 0 );
      node.SetIndex(thIt.GetIndex() );
      targetPoints->InsertElement(cc, node);
      cc++;
      }
    if( surf->GetPixel(thIt.GetIndex() ) == 1 )
      {
      NodeType node;
      node.SetValue( 0 );
      node.SetIndex(thIt.GetIndex() );
      trialPoints->InsertElement(cc, node);
      dd++;
      }
    ++thIt;
    }

  marcher->SetTargetReachedModeToAllTargets();
  marcher->SetAlivePoints( alivePoints );
  marcher->SetTrialPoints( trialPoints );
  marcher->SetTargetPoints( targetPoints );
  marcher->SetInput( gm );
  double stoppingValue = 1000.0;
  marcher->SetStoppingValue( stoppingValue );
  marcher->GenerateGradientImageOn();
  marcher->Update();
  WriteImage<ImageType>(marcher->GetOutput(), "marcher.nii.gz");

  thIt.GoToBegin();
  while(  !thIt.IsAtEnd()  )
    {
    typename TField::PixelType vec;
    for( dd = 0; dd < ImageDimension; dd++ )
      {
      vec[dd] = marcher->GetGradientImage()->GetPixel(thIt.GetIndex() )[dd];
      }
    ++thIt;
    }

  return sfield;
}

template <typename TImage, typename TField>
typename TField::Pointer
LaplacianGrad(typename TImage::Pointer wm, typename TImage::Pointer gm, float sig, unsigned int numits, float tolerance)
{
  typedef  typename TImage::IndexType IndexType;
  IndexType ind;
  typedef TImage ImageType;
  typedef TField GradientImageType;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
    GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;

  typename TField::Pointer sfield =
    AllocImage<TField>(wm);

  typename TImage::Pointer laplacian = SmoothImage<TImage>(wm, 1);
  laplacian->FillBuffer(0);
  typedef itk::ImageRegionIteratorWithIndex<TImage> IteratorType;
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  Iterator.GoToBegin();

  // initialize L(wm)=1, L(gm)=0.5, else 0
  while(  !Iterator.IsAtEnd()  )
    {
    ind = Iterator.GetIndex();
    if( wm->GetPixel(ind) >= 0.5f )
      {
      laplacian->SetPixel(ind, 1);
      }
    else
      {
      laplacian->SetPixel(ind, 2.);
      }
    ++Iterator;
    }

  // smooth and then reset the values
  float        meanvalue = 0, lastmean = 1;
  unsigned int iterations = 0;
  while( static_cast<float>( std::fabs(meanvalue - lastmean) ) > tolerance  && iterations < numits )
    {
    iterations++;
    std::cout << "  % " << (float) iterations
      / (float)(numits + 1) << " delta-mean " << fabs(meanvalue - lastmean) <<  std::endl;
    laplacian = SmoothImage<TImage>(laplacian, sqrt(sig) );
    Iterator.GoToBegin();
    unsigned int ct = 0;
    lastmean = meanvalue;
    while(  !Iterator.IsAtEnd()  )
      {
      ind = Iterator.GetIndex();
      if( wm->GetPixel(ind) >= 0.5f )
        {
        laplacian->SetPixel(ind, 1);
        }
      else if( gm->GetPixel(ind) < 0.5f  && wm->GetPixel(ind) < 0.5f )
        {
        laplacian->SetPixel(ind, 2.);
        }
      else
        {
        meanvalue += laplacian->GetPixel(ind);  ct++;
        }
      ++Iterator;
      }

    meanvalue /= (float)ct;
    }

  // /  WriteImage<ImageType>(laplacian, "laplacian.hdr");

  GradientImageFilterPointer filter = GradientImageFilterType::New();
  filter->SetInput(  laplacian );
  filter->SetSigma(sig * 0.5f);
  filter->Update();
  return filter->GetOutput();
}

template <typename TImage, typename TField, typename TInterp, typename TInterp2>
float IntegrateLength( typename TImage::Pointer /* gmsurf */,  typename TImage::Pointer /* thickimage */,
                       typename TImage::IndexType velind,  typename TField::Pointer lapgrad,  float itime,
                       float starttime, float /* finishtime */,
                       bool timedone, float deltaTime, typename TInterp::Pointer vinterp,
                       typename TInterp2::Pointer sinterp, unsigned int /* task */,
                       bool /* propagate */, bool domeasure,   unsigned int m_NumberOfTimePoints,
                       typename TImage::SpacingType spacing, float vecsign,
                       float timesign, float gradsign, unsigned int ct, typename TImage::Pointer wm,
                       typename TImage::Pointer gm,
                       float priorthickval,  typename TImage::Pointer smooththick, bool printprobability,
                       typename TImage::Pointer /* sulci */ )
{
  typedef typename TField::PixelType                               VectorType;
  typedef typename TField::PointType                               DPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TField, float> DefaultInterpolatorType;

  VectorType zero;
  zero.Fill(0);
  VectorType disp;
  disp.Fill(0);
  ct = 0;
  DPointType pointIn1;
  DPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::IndexType IndexType;
  for( unsigned int jj = 0; jj < ImageDimension; jj++ )
    {
    IndexType index;
    index[jj] = velind[jj];
    pointIn1[jj] = velind[jj] * lapgrad->GetSpacing()[jj];
    }
  // if( task == 0 )
  //   {
  //   propagate = false;
  //   }
  // else
  //   {
  //   propagate = true;
  //   }
  itime = starttime;
  timedone = false;
  float totalmag = 0;
  if( domeasure )
    {
    while( !timedone )
      {
      float scale = 1; // *m_DT[timeind]/m_DS[timeind];
      //     std::cout << " scale " << scale << std::endl;
      double itimetn1 = static_cast<double>( itime - timesign * deltaTime * scale );
      double itimetn1h = static_cast<double>( itime - timesign * deltaTime * 0.5f * scale );
      if( itimetn1h < 0 )
        {
        itimetn1h = 0;
        }
      if( itimetn1h > m_NumberOfTimePoints - 1 )
        {
        itimetn1h = m_NumberOfTimePoints - 1;
        }
      if( itimetn1 < 0 )
        {
        itimetn1 = 0;
        }
      if( itimetn1 > m_NumberOfTimePoints - 1 )
        {
        itimetn1 = m_NumberOfTimePoints - 1;
        }
      // first get current position of particle
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        IndexType index;
        index[jj] = velind[jj];
        pointIn1[jj] = velind[jj] * lapgrad->GetSpacing()[jj];
        }
      //      std::cout << " ind " << index  << std::endl;
      // now index the time varying field at that position.
      typename DefaultInterpolatorType::OutputType f1;  f1.Fill(0);
      typename DefaultInterpolatorType::OutputType f2;  f2.Fill(0);
      typename DefaultInterpolatorType::OutputType f3;  f3.Fill(0);
      typename DefaultInterpolatorType::OutputType f4;  f4.Fill(0);
      typename DefaultInterpolatorType::ContinuousIndexType  Y1;
      typename DefaultInterpolatorType::ContinuousIndexType  Y2;
      typename DefaultInterpolatorType::ContinuousIndexType  Y3;
      typename DefaultInterpolatorType::ContinuousIndexType  Y4;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        pointIn2[jj] = static_cast<typename DPointType::CoordRepType>( disp[jj] ) + pointIn1[jj];
        vcontind[jj] = pointIn2[jj] / lapgrad->GetSpacing()[jj];
        Y1[jj] = vcontind[jj];
        Y2[jj] = vcontind[jj];
        Y3[jj] = vcontind[jj];
        Y4[jj] = vcontind[jj];
        }
      // Y1[ImageDimension]=itimetn1;
      // Y2[ImageDimension]=itimetn1h;
      // Y3[ImageDimension]=itimetn1h;
      //      Y4[ImageDimension]=itime;

      using ContinuousIndexType = typename DefaultInterpolatorType::ContinuousIndexType;
      using CoordRepType = typename ContinuousIndexType::CoordRepType;

      f1 = vinterp->EvaluateAtContinuousIndex( Y1 );
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        Y2[jj] += static_cast<CoordRepType>( f1[jj] ) * static_cast<CoordRepType>( deltaTime ) * static_cast<CoordRepType>( 0.5 );
        }
      bool isinside = true;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        if( Y2[jj] < 1 || Y2[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2 )
          {
          isinside = false;
          }
        }
      if( isinside )
        {
        f2 = vinterp->EvaluateAtContinuousIndex( Y2 );
        }
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        Y3[jj] += static_cast<CoordRepType>( f2[jj] ) * static_cast<CoordRepType>( deltaTime ) * static_cast<CoordRepType>( 0.5 );
        }
      isinside = true;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        if( Y3[jj] < 1 || Y3[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2 )
          {
          isinside = false;
          }
        }
      if( isinside )
        {
        f3 = vinterp->EvaluateAtContinuousIndex( Y3 );
        }
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        Y4[jj] += static_cast<CoordRepType>( f3[jj] ) * static_cast<CoordRepType>( deltaTime );
        }
      isinside = true;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        if( Y4[jj] < 1 || Y4[jj] > lapgrad->GetLargestPossibleRegion().GetSize()[jj] - 2 )
          {
          isinside = false;
          }
        }
      if( isinside )
        {
        f4 = vinterp->EvaluateAtContinuousIndex( Y4 );
        }
      using DPointCoordRepType = typename DPointType::CoordRepType;
      DPointCoordRepType twoValue = static_cast<DPointCoordRepType>( 2.0 );
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        pointIn3[jj] = pointIn2[jj] + static_cast<DPointCoordRepType>( gradsign * vecsign * deltaTime / 6.0f )
          * ( static_cast<DPointCoordRepType>( f1[jj] ) + twoValue * static_cast<DPointCoordRepType>( f2[jj] )
          + twoValue * static_cast<DPointCoordRepType>( f3[jj] ) + static_cast<DPointCoordRepType>( f4[jj] ) );
        }

      VectorType out;
      float      mag = 0, dmag = 0;
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        out[jj] = pointIn3[jj] - pointIn1[jj];
        mag += static_cast<float>( itk::Math::sqr(pointIn3[jj] - pointIn2[jj]) );
        dmag += static_cast<float>( itk::Math::sqr(pointIn3[jj] - pointIn1[jj]) );
        disp[jj] = out[jj];
        }
      dmag = sqrt(dmag);
      totalmag += static_cast<float>( sqrt(mag) );

      ct++;
      //      if (!propagate) //thislength=dmag;//
//         thislength += totalmag;
      itime = itime + deltaTime * timesign;
      IndexType myind;
      for( unsigned int qq = 0; qq <  ImageDimension; qq++ )
        {
        myind[qq] = (unsigned long)(pointIn3[qq] / spacing[qq] + 0.5);
        }

      if( (gm->GetPixel(myind) < 0.5f && wm->GetPixel(myind) < 0.5f) ||
          (wm->GetPixel(myind) >= 0.5f && gm->GetPixel(myind) < 0.5f) ||
          mag < 1.e-1f * deltaTime )
        {
        timedone = true;
        }
      if( gm->GetPixel(myind) < 0.5f )
        {
        timedone = true;
        }
      if( static_cast<float>( ct ) >  2.0f / deltaTime )
        {
        timedone = true;
        }
      if( totalmag >  priorthickval )
        {
        timedone = true;
        }
      if( smooththick )
        {
        if( (totalmag - smooththick->GetPixel(velind) ) > 1 )
          {
          timedone = true;
          }
        }

      if( printprobability )
        {
        std::cout << " ind " << Y1 << " prob " << sinterp->EvaluateAtContinuousIndex(Y1) << " t " << itime << std::endl;
        }
      }
    }

  return totalmag;
}

template <unsigned int ImageDimension>
int LaplacianThickness(int argc, char *argv[])
{
  float        gradstep = -50.0; // atof(argv[3])*(-1.0);
  unsigned int nsmooth = 2;
  float        smoothparam = 1;
  float        priorthickval = 500;
  double       dT = 0.01;
  std::string  wfn = std::string(argv[1]);
  std::string  gfn = std::string(argv[2]);
  int          argct = 3;
  std::string  outname = std::string(argv[argct]); argct++;

  if( argc > argct )
    {
    smoothparam = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    priorthickval = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    dT = atof(argv[argct]);
    }
  argct++;
  float dosulc = 0;
  if( argc >  argct )
    {
    dosulc = atof(argv[argct]);
    }
  argct++;
  float tolerance = 0.001;
  if( argc >  argct )
    {
    tolerance = atof(argv[argct]);
    }
  argct++;
  std::cout << " using tolerance " << tolerance << std::endl;
  typedef float                                      PixelType;
  typedef itk::Vector<float, ImageDimension>         VectorType;
  typedef itk::Image<VectorType, ImageDimension>     DisplacementFieldType;
  typedef itk::Image<PixelType, ImageDimension>      ImageType;
  typedef typename  ImageType::SpacingType           SpacingType;

  //  typename tvt::Pointer gWarp;
  // ReadImage<tvt>( gWarp, ifn.c_str() );

  typename ImageType::Pointer thickimage;
  ReadImage<ImageType>(thickimage, wfn.c_str() );
  thickimage->FillBuffer(0);
  typename ImageType::Pointer thickimage2;
  ReadImage<ImageType>(thickimage2, wfn.c_str() );
  thickimage2->FillBuffer(0);
  typename ImageType::Pointer wm;
  ReadImage<ImageType>(wm, wfn.c_str() );
  typename ImageType::Pointer gm;
  ReadImage<ImageType>(gm, gfn.c_str() );
  SpacingType spacing = wm->GetSpacing();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  typename ImageType::Pointer wmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, wm);
  typename DisplacementFieldType::Pointer lapgrad = nullptr;
  typename DisplacementFieldType::Pointer lapgrad2 = nullptr;
  typename ImageType::Pointer gmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, gm);

/** get sulcal priors */
  typename ImageType::Pointer sulci = nullptr;
  if( dosulc > 0 )
    {
    std::cout << "  using sulcal prior " << std::endl;
    typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;
    typename  FilterType::Pointer distmap = FilterType::New();
    distmap->InputIsBinaryOn();
    distmap->SetUseImageSpacing(true);
    distmap->SetInput(wmb);
    distmap->Update();
    typename ImageType::Pointer distwm = distmap->GetOutput();

    typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType> dgf;
    typename dgf::Pointer lfilter = dgf::New();
    lfilter->SetSigma(smoothparam);
    lfilter->SetInput(distwm);
    lfilter->Update();
    typename ImageType::Pointer image2 = lfilter->GetOutput();
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(   0 );
    rescaler->SetOutputMaximum( 1 );
    rescaler->SetInput( image2 );
    rescaler->Update();
    sulci =  rescaler->GetOutput();
    WriteImage<ImageType>(sulci, "sulci.nii");

    Iterator.GoToBegin();
    while(  !Iterator.IsAtEnd()  )
      {
//    std::cout << " a good value for use sulcus prior is 0.002  -- in a function :
//  1/(1.+exp(-0.1*(sulcprob-0.275)/use-sulcus-prior)) " << std::endl;
//
      float gmprob = gm->GetPixel(Iterator.GetIndex() );
      if( itk::Math::FloatAlmostEqual( gmprob, 0.0f ) )
        {
        gmprob = 0.05f;
        }
      float sprob = sulci->GetPixel(Iterator.GetIndex() );
      sprob = 1.0f / (1.0f + std::exp(-0.1f * (sprob - 0.5f) / dosulc) );
      sulci->SetPixel(Iterator.GetIndex(), sprob );
//    if (gmprob > 0) std::cout << " gmp " << gmprob << std::endl;
      ++Iterator;
      }

    std::cout << " modified gm prior by sulcus prior " << std::endl;
    WriteImage<ImageType>(sulci, "sulcigm.nii");

    typedef itk::GradientRecursiveGaussianImageFilter<ImageType, DisplacementFieldType>
      GradientImageFilterType;
    typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
    GradientImageFilterPointer filter = GradientImageFilterType::New();
    filter->SetInput(  distwm );
    filter->SetSigma(smoothparam);
    filter->Update();
    lapgrad2 = filter->GetOutput();

//      return 0;
/** sulc priors done */
    }

  lapgrad = LaplacianGrad<ImageType, DisplacementFieldType>(wmb, gmb, smoothparam, 500, tolerance);
  //  lapgrad=FMMGrad<ImageType,DisplacementFieldType>(wmb,gmb);

  //  LabelSurface(typename TImage::PixelType foreground,
  //       typename TImage::PixelType newval, typename TImage::Pointer input, float distthresh )
  float distthresh = 1.9;

  typename ImageType::Pointer wmgrow = Morphological<ImageType>(wmb, 1, true);
  typename ImageType::Pointer surf = LabelSurface<ImageType>(1, 1, wmgrow, distthresh);
  typename ImageType::Pointer gmsurf = LabelSurface<ImageType>(1, 1, gmb, distthresh);
  // now integrate
  //

  double timezero = 0; // 1
  double timeone = 1;  // (s[ImageDimension]-1-timezero);

  //  unsigned int m_NumberOfTimePoints = s[ImageDimension];

  float starttime = timezero; // timezero;
  float finishtime = timeone; // s[ImageDimension]-1;//timeone;
  // std::cout << " MUCKING WITH START FINISH TIME " <<  finishtime <<  std::endl;

  typename DisplacementFieldType::IndexType velind;
  typename ImageType::Pointer smooththick = nullptr;
  float timesign = 1.0;
  if( starttime  >  finishtime )
    {
    timesign = -1.0;
    }
  unsigned int m_NumberOfTimePoints = 2;
  typedef   DisplacementFieldType                                                        TimeVaryingVelocityFieldType;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float> DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer sinterp =  ScalarInterpolatorType::New();
  sinterp->SetInputImage(gm);
  if( sulci )
    {
    sinterp->SetInputImage(sulci);
    }
  VectorType zero;
  zero.Fill(0);

  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> VIteratorType;
  VIteratorType VIterator( lapgrad, lapgrad->GetLargestPossibleRegion().GetSize() );
  VIterator.GoToBegin();
  while(  !VIterator.IsAtEnd()  )
    {
    VectorType vec = VIterator.Get();
    float      mag = 0;
    for( unsigned int qq = 0; qq < ImageDimension; qq++ )
      {
      mag += vec[qq] * vec[qq];
      }
    mag = sqrt(mag);
    if( mag > 0 )
      {
      vec = vec / mag;
      }
    VIterator.Set(vec * gradstep);
    if( lapgrad2 )
      {
      vec = lapgrad2->GetPixel(VIterator.GetIndex() );
      mag = 0;
      for( unsigned int qq = 0; qq < ImageDimension; qq++ )
        {
        mag += vec[qq] * vec[qq];
        }
      mag = sqrt(mag);
      if( mag > 0 )
        {
        vec = vec / mag;
        }
      lapgrad2->SetPixel(VIterator.GetIndex(), vec * gradstep);
      }
    ++VIterator;
    }

  bool propagate = false;
  for( unsigned int smoothit = 0; smoothit < nsmooth; smoothit++ )
    {
    std::cout << " smoothit " << smoothit << std::endl;
    Iterator.GoToBegin();
    unsigned int cter = 0;
    while(  !Iterator.IsAtEnd()  )
      {
      velind = Iterator.GetIndex();
      //      float thislength=0;
      for( unsigned int task = 0; task < 1; task++ )
        {
        float itime = starttime;

        unsigned long ct = 0;
        bool          timedone = false;

        VectorType disp;
        disp.Fill(0.0);
        double deltaTime = dT, vecsign = 1.0;
        bool   domeasure = false;
        float  gradsign = 1.0;
        bool   printprobability = false;
//    std::cout << " wmb " << wmb->GetPixel(velind) << " gm " << gm->GetPixel(velind) << std::endl;
//    if (surf->GetPixel(velind) != 0) printprobability=true;
        if( gm->GetPixel(velind) > 0.25f ) // && wmb->GetPixel(velind) < 1 )
          {
          cter++;
          domeasure = true;
          }
        vinterp->SetInputImage(lapgrad);
        gradsign = -1.0; vecsign = -1.0;
        float len1 = IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
            (gmsurf, thickimage, velind, lapgrad,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
            sinterp, task, propagate, domeasure, m_NumberOfTimePoints, spacing, vecsign, gradsign, timesign, ct, wm, gm,
            priorthickval, smooththick, printprobability,
            sulci );

        gradsign = 1.0;  vecsign = 1;
        float len2 = IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
            (gmsurf, thickimage, velind, lapgrad,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
            sinterp, task, propagate, domeasure, m_NumberOfTimePoints, spacing, vecsign, gradsign, timesign, ct, wm, gm,
            priorthickval - len1, smooththick, printprobability,
            sulci );

        float len3 = 1.e9, len4 = 1.e9;
        if( lapgrad2 )
          {
          vinterp->SetInputImage(lapgrad2);
          gradsign = -1.0; vecsign = -1.0;
          len3 = IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
              (gmsurf, thickimage, velind, lapgrad2,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
              sinterp, task, propagate, domeasure, m_NumberOfTimePoints, spacing, vecsign, gradsign, timesign, ct, wm,
              gm,
              priorthickval, smooththick, printprobability,
              sulci );

          gradsign = 1.0;  vecsign = 1;
          len4 = IntegrateLength<ImageType, DisplacementFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
              (gmsurf, thickimage, velind, lapgrad2,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
              sinterp, task, propagate, domeasure, m_NumberOfTimePoints, spacing, vecsign, gradsign, timesign, ct, wm,
              gm,
              priorthickval - len3, smooththick, printprobability,
              sulci );
          }
        float totalength = len1 + len2;
//    if (totalength > 5 && totalength <  8) std::cout<< " t1 " << len3+len4 << " t2 " << len1+len2 << std::endl;
        if( len3 + len4 < totalength )
          {
          totalength = len3 + len4;
          }

        if( smoothit == 0 )
          {
          if( itk::Math::FloatAlmostEqual( thickimage2->GetPixel(velind), 0.0f )  )
            {
            thickimage2->SetPixel(velind, totalength);
            }
          else if( (totalength) > 0 &&  thickimage2->GetPixel(velind) < (totalength) )
            {
            thickimage2->SetPixel(velind, totalength);
            }
          }
        if( smoothit > 0 && smooththick )
          {
          thickimage2->SetPixel(velind, totalength * 0.5f + smooththick->GetPixel(velind) * 0.5f );
          }

        if( domeasure && (totalength) > 0 && cter % 10000 == 0 )
          {
          std::cout << " len1 " << len1 << " len2 " << len2 << " ind " << velind << std::endl;
          }
        }
      ++Iterator;
      }

    smooththick = SmoothImage<ImageType>(thickimage2, 1.0);

// set non-gm voxels to zero
    IteratorType gIterator( gm, gm->GetLargestPossibleRegion().GetSize() );
    gIterator.GoToBegin();
    while(  !gIterator.IsAtEnd()  )
      {
      if( gm->GetPixel(gIterator.GetIndex() ) < 0.25f )
        {
        thickimage2->SetPixel(gIterator.GetIndex(), 0);
        }
      ++gIterator;
      }

    std::cout << " writing " << outname << std::endl;
    WriteImage<ImageType>(thickimage2, outname.c_str() );
    }
//  WriteImage<ImageType>(thickimage,"turd.hdr");

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int LaplacianThickness( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "LaplacianThickness" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );

  if( argc < 4 )
    {
    std::cout << "Usage:   " << argv[0]
             <<
      " WM.nii GM.nii   Out.nii  {smoothparam=1} {priorthickval=500} {dT=0.01} {sulcus-prior=0} {laplacian-tolerance=0.001}"
             << std::endl;
    std::cout
      <<
      " a good value for sulcus prior (if not 0, which disables its use) is 0.15 -- in a function :  1/(1.+exp(-0.1*(laplacian-img-value-sulcprob)/0.01)) "
      << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  std::string ifn = std::string(argv[1]);
  //  std::cout << " image " << ifn << std::endl;
  // Get the image dimension
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(ifn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(ifn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim =  imageIO->GetNumberOfDimensions();

  //   std::cout << " dim " << dim << std::endl;
  switch( dim )
    {
    case 2:
      {
      return LaplacianThickness<2>(argc, argv);
      }
      break;
    case 3:
      {
      return LaplacianThickness<3>(argc, argv);
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // namespace ants
