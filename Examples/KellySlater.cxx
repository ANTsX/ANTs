// #include "DoSomethingToImage.cxx"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"

#include "itkWarpImageFilter.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkDisplacementFieldFromMultiTransformFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkANTSImageRegistrationOptimizer.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCovariantVector.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "ReadWriteImage.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

#include "itkCentralDifferenceImageFunction.h"
#include "itkSurfaceCurvatureBase.h"
#include "itkSurfaceImageCurvature.h"

template <class TField, class TImage>
typename TImage::Pointer
GetVectorComponent(typename TField::Pointer field, unsigned int index)
{
  // Initialize the Moving to the displacement field
  typedef TField FieldType;
  typedef TImage ImageType;

  typename ImageType::Pointer sfield = ImageType::New();
  sfield->SetSpacing( field->GetSpacing() );
  sfield->SetOrigin( field->GetOrigin() );
  sfield->SetDirection( field->GetDirection() );
  sfield->SetLargestPossibleRegion(field->GetLargestPossibleRegion() );
  sfield->SetRequestedRegion(field->GetRequestedRegion() );
  sfield->SetBufferedRegion( field->GetBufferedRegion() );
  sfield->Allocate();

  typedef itk::ImageRegionIteratorWithIndex<TField> Iterator;
  Iterator vfIter( field,  field->GetLargestPossibleRegion() );
  for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    typename TField::PixelType v1 = vfIter.Get();
    sfield->SetPixel(vfIter.GetIndex(), v1[index]);
    }

  return sfield;
}

template <class TImage>
typename TImage::Pointer BinaryThreshold(
  typename TImage::PixelType low,
  typename TImage::PixelType high,
  typename TImage::PixelType replaceval, typename TImage::Pointer input)
{
  // std::cout << " Binary Thresh " << std::endl;

  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef itk::BinaryThresholdImageFilter<TImage, TImage> InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  inputThresholder->SetInsideValue(  replaceval );
  int outval = 0;
  if( (double) replaceval == (double) -1 )
    {
    outval = 1;
    }
  inputThresholder->SetOutsideValue( outval );

  if( high < low )
    {
    high = 255;
    }
  inputThresholder->SetLowerThreshold( (PixelType) low );
  inputThresholder->SetUpperThreshold( (PixelType) high);
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <class TImage>
typename TImage::Pointer
MaurerDistanceMap(
  typename TImage::PixelType pixlo,
  typename TImage::PixelType pixhi,
  typename TImage::Pointer input)
{
  // std::cout << " DDMap " << std::endl;

  typedef TImage ImageType;

  typedef itk::SignedMaurerDistanceMapImageFilter<
      ImageType, ImageType>  FilterType;

  typename  FilterType::Pointer filter = FilterType::New();
  filter->SetSquaredDistance(false);
  //  filter->InputIsBinaryOn();
  filter->SetUseImageSpacing(true);
  filter->SetBackgroundValue(0);
  filter->SetInput(BinaryThreshold<TImage>(pixlo, pixhi, pixhi, input) );
  filter->Update();

  //  WriteImage<ImageType>(filter->GetOutput(),"temp1.nii");
  return filter->GetOutput();
}

template <class TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer image, double sig)
{
// find min value
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator vfIter(image, image->GetLargestPossibleRegion() );
  for( vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    typename TImage::PixelType v1 = vfIter.Get();
    if( vnl_math_isnan(v1) )
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

template <class TImage>
void
SmoothDeformation(typename TImage::Pointer vectorimage, double sig)
{
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::PixelType            RealType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<RealType, ImageDimension>  ImageType;
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

template <class TImage, class TDisplacementField>
typename TImage::Pointer
CopyImage(TDisplacementField* field )
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  //  unsigned int row=0;
  // unsigned int col=0;
  typedef typename TImage::PixelType            PixelType;
  typedef itk::Image<PixelType, ImageDimension> RealImageType;
  typename RealImageType::RegionType m_JacobianRegion;

  typename RealImageType::Pointer m_RealImage = NULL;
  m_RealImage = RealImageType::New();
  m_RealImage->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  m_RealImage->SetBufferedRegion( field->GetLargestPossibleRegion().GetSize() );
  m_RealImage->SetSpacing(field->GetSpacing() );
  m_RealImage->SetDirection( field->GetDirection() );
  m_RealImage->SetOrigin(field->GetOrigin() );
  m_RealImage->Allocate();
  m_RealImage->FillBuffer(0);

  return m_RealImage;
}

template <class TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground,
             typename TImage::PixelType newval, typename TImage::Pointer input, double distthresh )
{
  std::cout << " Label Surf " << std::endl;

  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typename   ImageType::Pointer     Image = ImageType::New();
  Image->SetLargestPossibleRegion(input->GetLargestPossibleRegion()  );
  Image->SetBufferedRegion(input->GetLargestPossibleRegion() );
  Image->Allocate();
  Image->FillBuffer( 0.0 );
  Image->SetSpacing(input->GetSpacing() );
  Image->SetOrigin(input->GetOrigin() );
  typedef itk::NeighborhoodIterator<ImageType> iteratorType;

  typename iteratorType::RadiusType rad;
  for( int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = (unsigned int)(distthresh + 0.5);
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

//  std::cout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    typename TImage::IndexType ind2;
    if( p == foreground )
      {
      bool atedge = false;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        ind2 = GHood.GetIndex(i);
        double dist = 0.0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          dist += (double)(ind[j] - ind2[j]) * (double)(ind[j] - ind2[j]);
          }
        dist = sqrt(dist);
        if( GHood.GetPixel(i) != foreground && dist <  distthresh  )
          {
          atedge = true;
          }
        }
      if( atedge && p == foreground )
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

template <class TImage>
typename TImage::Pointer  Morphological( typename TImage::Pointer input, double rad, bool option)
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
      if( o_iter.Get() > 0.5 && input->GetPixel(o_iter.GetIndex() ) > 0.5 )
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

template <class TImage, class TField>
typename TField::Pointer
LaplacianGrad(typename TImage::Pointer wm, typename TImage::Pointer gm, double sig)
{
  typedef  typename TImage::IndexType IndexType;
  IndexType ind;
  typedef TImage ImageType;
  typedef TField GradientImageType;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
    GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;

  typename TField::Pointer sfield = TField::New();
  sfield->SetSpacing( wm->GetSpacing() );
  sfield->SetOrigin( wm->GetOrigin() );
  sfield->SetDirection( wm->GetDirection() );
  sfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  sfield->SetRequestedRegion(wm->GetRequestedRegion() );
  sfield->SetBufferedRegion( wm->GetBufferedRegion() );
  sfield->Allocate();

  typename TImage::Pointer laplacian = SmoothImage<TImage>(wm, 3);
  laplacian->FillBuffer(0);
  typedef itk::ImageRegionIteratorWithIndex<TImage> IteratorType;
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  Iterator.GoToBegin();

  // initialize L(wm)=1, L(gm)=0.5, else 0
  while(  !Iterator.IsAtEnd()  )
    {
    ind = Iterator.GetIndex();
    if( wm->GetPixel(ind) )
      {
      laplacian->SetPixel(ind, 1);
      }
    else if( gm->GetPixel(ind) && wm->GetPixel(ind) == 0 )
      {
      laplacian->SetPixel(ind, 0.);
      }
    else
      {
      laplacian->SetPixel(ind, 0.);
      }
    ++Iterator;
    }
  // smooth and then reset the values
  for( unsigned int iterations = 0; iterations < 100; iterations++ )
    {
    laplacian = SmoothImage<TImage>(laplacian, sqrt(sig) );
    Iterator.GoToBegin();
    while(  !Iterator.IsAtEnd()  )
      {
      ind = Iterator.GetIndex();
      if( wm->GetPixel(ind) )
        {
        laplacian->SetPixel(ind, 1);
        }
      else if( gm->GetPixel(ind) == 0 && wm->GetPixel(ind) == 0 )
        {
        laplacian->SetPixel(ind, 0.);
        }
      ++Iterator;
      }
    }

  GradientImageFilterPointer filter = GradientImageFilterType::New();
  filter->SetInput(  laplacian );
  filter->SetSigma(sig);
  filter->Update();
  return filter->GetOutput();
}

template <class TImage, class TField>
typename TField::Pointer
ExpDiffMap(typename TField::Pointer velofield,  typename TImage::Pointer wm,  double sign, unsigned int numtimepoints )
{
  typedef TImage                     ImageType;
  typedef TField                     DisplacementFieldType;
  typedef typename TField::PixelType PixelType;
  typename TField::PixelType zero, disp;
  enum { ImageDimension = TImage::ImageDimension };
  disp.Fill(0);
  zero.Fill(0);

  typename DisplacementFieldType::Pointer incrfield = DisplacementFieldType::New();
  incrfield->SetSpacing( velofield->GetSpacing() );
  incrfield->SetOrigin( velofield->GetOrigin() );
  incrfield->SetDirection( velofield->GetDirection() );
  incrfield->SetLargestPossibleRegion(velofield->GetLargestPossibleRegion() );
  incrfield->SetRequestedRegion(velofield->GetRequestedRegion() );
  incrfield->SetBufferedRegion( velofield->GetBufferedRegion() );
  incrfield->Allocate();
  incrfield->FillBuffer(zero);

  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  Iterator.GoToBegin();
  while(  !Iterator.IsAtEnd()  )
    {
    incrfield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() ) * sign);
    ++Iterator;
    }

  // generate phi
  typedef itk::MatrixOffsetTransformBase<PixelType, ImageDimension, ImageDimension>           AffineTransformType;
  typedef itk::DisplacementFieldFromMultiTransformFilter<TField, TField, AffineTransformType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetOutputSize(velofield->GetLargestPossibleRegion().GetSize() );
  warper->SetOutputSpacing(velofield->GetSpacing() );
  warper->SetOutputOrigin(velofield->GetOrigin() );
  warper->SetOutputDirection(velofield->GetDirection() );
  warper->DetermineFirstDeformNoInterp();

  unsigned int ttiter = 0;
  while( ttiter < numtimepoints )      // 10 time integration points
    {
    ttiter++;
    warper->PushBackDisplacementFieldTransform(incrfield);
    }

  warper->Update();
  return warper->GetOutput();
}

template <class TImage, class TField>
typename TField::Pointer
DiReCTCompose(typename TField::Pointer velofield, typename TField::Pointer diffmap )
{
  typedef TImage                     ImageType;
  typedef TField                     DisplacementFieldType;
  typedef typename TField::PixelType PixelType;
  typename TField::PixelType zero, disp;
  enum { ImageDimension = TImage::ImageDimension };
  disp.Fill(0);
  zero.Fill(0);

  typedef itk::MatrixOffsetTransformBase<PixelType, ImageDimension, ImageDimension>           AffineTransformType;
  typedef itk::DisplacementFieldFromMultiTransformFilter<TField, TField, AffineTransformType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetOutputSize(velofield->GetLargestPossibleRegion().GetSize() );
  warper->SetOutputSpacing(velofield->GetSpacing() );
  warper->SetOutputOrigin(velofield->GetOrigin() );
  warper->SetOutputDirection(velofield->GetDirection() );
  warper->DetermineFirstDeformNoInterp();
  warper->PushBackDisplacementFieldTransform(diffmap);
  warper->PushBackDisplacementFieldTransform(velofield);
  warper->Update();
  return warper->GetOutput();
}

template <class TImage, class TField>
void
InvertField( typename TField::Pointer field,
             typename TField::Pointer inverseFieldIN, double weight = 1.0,
             double toler = 0.1, int maxiter = 20, bool /* print */ = false)
{
  enum { ImageDimension = TImage::ImageDimension };
  typedef TField                     DisplacementFieldType;
  typedef typename TField::Pointer   DisplacementFieldPointer;
  typedef typename TField::PixelType VectorType;
  typedef TImage                     ImageType;
  typedef typename TImage::Pointer   ImagePointer;
  double       mytoler = toler;
  unsigned int mymaxiter = maxiter;

  VectorType zero; zero.Fill(0);
  //  if (this->GetElapsedIterations() < 2 ) maxiter=10;

  ImagePointer realImage = ImageType::New();
  realImage->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  realImage->SetBufferedRegion( field->GetLargestPossibleRegion().GetSize() );
  realImage->SetSpacing(field->GetSpacing() );
  realImage->SetOrigin(field->GetOrigin() );
  realImage->SetDirection(field->GetDirection() );
  realImage->Allocate();

  typedef typename DisplacementFieldType::PixelType                VectorType;
  typedef typename DisplacementFieldType::IndexType                IndexType;
  typedef typename VectorType::ValueType                           ScalarType;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;

  typedef itk::ANTSImageRegistrationOptimizer<ImageDimension, double> ROType;
  typename ROType::Pointer m_MFR = ROType::New();

  DisplacementFieldPointer inverseField = DisplacementFieldType::New();
  inverseField->SetSpacing( field->GetSpacing() );
  inverseField->SetOrigin( field->GetOrigin() );
  inverseField->SetDirection( field->GetDirection() );
  inverseField->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  inverseField->SetRequestedRegion(field->GetRequestedRegion() );
  inverseField->SetBufferedRegion( field->GetLargestPossibleRegion() );
  inverseField->Allocate();
  inverseField->FillBuffer(zero);

  DisplacementFieldPointer lagrangianInitCond = DisplacementFieldType::New();
  lagrangianInitCond->SetSpacing( field->GetSpacing() );
  lagrangianInitCond->SetOrigin( field->GetOrigin() );
  lagrangianInitCond->SetDirection( field->GetDirection() );
  lagrangianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  lagrangianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
  lagrangianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
  lagrangianInitCond->Allocate();
  DisplacementFieldPointer eulerianInitCond = DisplacementFieldType::New();
  eulerianInitCond->SetSpacing( field->GetSpacing() );
  eulerianInitCond->SetOrigin( field->GetOrigin() );
  eulerianInitCond->SetDirection( field->GetDirection() );
  eulerianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  eulerianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
  eulerianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
  eulerianInitCond->Allocate();

  typedef typename DisplacementFieldType::SizeType SizeType;
  SizeType size = field->GetLargestPossibleRegion().GetSize();

  typename ImageType::SpacingType spacing = field->GetSpacing();
  unsigned long npix = 1;
  for( int j = 0; j < ImageDimension; j++ )  // only use in-plane spacing
    {
    npix *= field->GetLargestPossibleRegion().GetSize()[j];
    }

  double   max = 0;
  Iterator iter( field, field->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    IndexType  index = iter.GetIndex();
    VectorType vec1 = iter.Get();
    VectorType newvec = vec1 * weight;
    lagrangianInitCond->SetPixel(index, newvec);
    inverseField->SetPixel(index, inverseFieldIN->GetPixel(index) );

    double mag = 0;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      mag += newvec[jj] * newvec[jj];
      }
    mag = sqrt(mag);
    if( mag > max )
      {
      max = mag;
      }
    }

  eulerianInitCond->FillBuffer(zero);

  double scale = (1.) / max;
  if( scale > 1. )
    {
    scale = 1.0;
    }
//    double initscale=scale;
  Iterator vfIter( inverseField, inverseField->GetLargestPossibleRegion() );

//  int num=10;
//  for (int its=0; its<num; its++)
  double       difmag = 10.0;
  unsigned int ct = 0;
  double       meandif = 1.e8;
//    int badct=0;
//  while (difmag > subpix && meandif > subpix*0.1 && badct < 2 )//&& ct < 20 && denergy > 0)
//    double length=0.0;
  double stepl = 2.;

  double epsilon = (double)size[0] / 256;
  if( epsilon > 1 )
    {
    epsilon = 1;
    }

  while( difmag > mytoler && ct<mymaxiter && meandif> 0.001 )
    {
    meandif = 0.0;

    // this field says what position the eulerian field should contain in the E domain
    m_MFR->ComposeDiffs(inverseField, lagrangianInitCond,    eulerianInitCond, 1);
    difmag = 0.0;
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      IndexType  index = vfIter.GetIndex();
      VectorType update = eulerianInitCond->GetPixel(index);
      double     mag = 0;
      for( int j = 0; j < ImageDimension; j++ )
        {
        update[j] *= (-1.0);
        mag += (update[j] / spacing[j]) * (update[j] / spacing[j]);
        }
      mag = sqrt(mag);
      meandif += mag;
      if( mag > difmag )
        {
        difmag = mag;
        }
      //      if (mag < 1.e-2) update.Fill(0);

      eulerianInitCond->SetPixel(index, update);
      realImage->SetPixel(index, mag);
      }
    meandif /= (double)npix;
    if( ct == 0 )
      {
      epsilon = 0.75;
      }
    else
      {
      epsilon = 0.5;
      }
    stepl = difmag * epsilon;
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      double     val = realImage->GetPixel(vfIter.GetIndex() );
      VectorType update = eulerianInitCond->GetPixel(vfIter.GetIndex() );
      if( val > stepl )
        {
        update = update * (stepl / val);
        }
      VectorType upd = vfIter.Get() + update * (epsilon);
      vfIter.Set(upd);
      }
    ct++;
    }
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    IndexType index = iter.GetIndex();
    inverseFieldIN->SetPixel(index, inverseField->GetPixel(index) );
    }

  //    std::cout <<" difmag " << difmag << ": its " << ct <<  " len " << m_MFR->MeasureDeformation(inverseField ) <<
  // std::endl;

  return;
}

template <unsigned int ImageDimension>
int LaplacianThicknessExpDiff2(int argc, char *argv[])
{
  int         argct = 2;
  std::string segfn = std::string(argv[argct]); argct++;

  std::string  wfn = std::string(argv[argct]); argct++;
  std::string  gfn = std::string(argv[argct]); argct++;
  std::string  outname = std::string(argv[argct]); argct++;
  unsigned int numtimepoints = 10;

  typedef double RealType;
  RealType gradstep = (RealType)(-1.0) * 0.5; // (ImageDimension-1);
  if( argc > argct )
    {
    gradstep = atof(argv[argct]) * (-1.0);
    }
  gradstep *= 1.0 / (RealType)numtimepoints * 10;  argct++;
  unsigned int alltheits = 50;
  if( argc > argct )
    {
    alltheits = atoi(argv[argct]);
    }
  argct++;
  RealType thickprior = 6.0;
  if( argc > argct )
    {
    thickprior = atof(argv[argct]);
    }
  argct++;
  // bool useCurvaturePrior = false;
  // if( argc > argct )
  //   {
  //   useCurvaturePrior = atoi(argv[argct]);
  //   }
  argct++;
  RealType smoothingsigma = 1.5;
  if( argc > argct )
    {
    smoothingsigma = atof(argv[argct]);
    }
  argct++;
  // bool useEuclidean = true;
  // if( argc > argct )
  //   {
  //   useEuclidean = atoi(argv[argct]);
  //   }
  argct++;
  std::cout << " smooth " << smoothingsigma << " thp " << thickprior << " gs " << gradstep << std::endl;
  typedef RealType                                                      PixelType;
  typedef itk::Vector<RealType, ImageDimension>                         VectorType;
  typedef itk::Image<VectorType, ImageDimension>                        DisplacementFieldType;
  typedef itk::Image<PixelType, ImageDimension>                         ImageType;
  typedef itk::ImageFileReader<ImageType>                               readertype;
  typedef itk::ImageFileWriter<ImageType>                               writertype;
  typedef typename  ImageType::IndexType                                IndexType;
  typedef typename  ImageType::SizeType                                 SizeType;
  typedef typename  ImageType::SpacingType                              SpacingType;
  typedef itk::Image<VectorType, ImageDimension + 1>                    tvt;
  typedef itk::ANTSImageRegistrationOptimizer<ImageDimension, RealType> ROType;
  typename ROType::Pointer m_MFR = ROType::New();

  typename ImageType::Pointer segmentationimage;
  ReadImage<ImageType>(segmentationimage, segfn.c_str() );
  typename ImageType::Pointer wm;
  ReadImage<ImageType>(wm, wfn.c_str() );
  typename ImageType::DirectionType omat = wm->GetDirection();
  typename ImageType::DirectionType fmat = wm->GetDirection();
  fmat.SetIdentity();
  std::cout << " Setting Identity Direction  " << fmat << std::endl;
  wm->SetDirection(fmat);
  typename ImageType::Pointer totalimage;
  ReadImage<ImageType>(totalimage, wfn.c_str() );
  totalimage->SetDirection(fmat);
  typename ImageType::Pointer hitimage;
  ReadImage<ImageType>(hitimage, wfn.c_str() );
  hitimage->SetDirection(fmat);
  typename ImageType::Pointer gm;
  ReadImage<ImageType>(gm, gfn.c_str() );
  gm->SetDirection(fmat);
  wm->SetDirection(fmat);
  segmentationimage->SetDirection(fmat);

  typename DisplacementFieldType::Pointer lapgrad;
  typename ImageType::Pointer gmb = BinaryThreshold<ImageType>(2, 2, 1, segmentationimage);  // fixme
  typename ImageType::Pointer wmb = BinaryThreshold<ImageType>(3, 3, 1, segmentationimage);  // fixme
  typename ImageType::Pointer laplacian = SmoothImage<ImageType>(wm, smoothingsigma);
  lapgrad = LaplacianGrad<ImageType, DisplacementFieldType>(wmb, gmb, 1);

  typename DisplacementFieldType::Pointer corrfield = DisplacementFieldType::New();
  corrfield->SetSpacing( wm->GetSpacing() );
  corrfield->SetOrigin( wm->GetOrigin() );
  corrfield->SetDirection( wm->GetDirection() );
  corrfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  corrfield->SetRequestedRegion(wm->GetRequestedRegion() );
  corrfield->SetBufferedRegion( wm->GetBufferedRegion() );
  corrfield->Allocate();
  VectorType zero;
  zero.Fill(0);
  corrfield->FillBuffer(zero);
  typename DisplacementFieldType::Pointer incrfield = DisplacementFieldType::New();
  incrfield->SetSpacing( wm->GetSpacing() );
  incrfield->SetOrigin( wm->GetOrigin() );
  incrfield->SetDirection( wm->GetDirection() );
  incrfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrfield->Allocate();
  incrfield->FillBuffer(zero);

  typename DisplacementFieldType::Pointer invfield = DisplacementFieldType::New();
  invfield->SetSpacing( wm->GetSpacing() );
  invfield->SetOrigin( wm->GetOrigin() );
  invfield->SetDirection( wm->GetDirection() );
  invfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  invfield->SetRequestedRegion(wm->GetRequestedRegion() );
  invfield->SetBufferedRegion( wm->GetBufferedRegion() );
  invfield->Allocate();
  invfield->FillBuffer(zero);

  typename DisplacementFieldType::Pointer incrinvfield = DisplacementFieldType::New();
  incrinvfield->SetSpacing( wm->GetSpacing() );
  incrinvfield->SetOrigin( wm->GetOrigin() );
  incrinvfield->SetDirection( wm->GetDirection() );
  incrinvfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrinvfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrinvfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrinvfield->Allocate();
  incrinvfield->FillBuffer(zero);

  typename DisplacementFieldType::Pointer velofield = DisplacementFieldType::New();
  velofield->SetSpacing( wm->GetSpacing() );
  velofield->SetOrigin( wm->GetOrigin() );
  velofield->SetDirection( wm->GetDirection() );
  velofield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  velofield->SetRequestedRegion(wm->GetRequestedRegion() );
  velofield->SetBufferedRegion( wm->GetBufferedRegion() );
  velofield->Allocate();
  velofield->FillBuffer(zero);

  //  LabelSurface(typename TImage::PixelType foreground,
  //       typename TImage::PixelType newval, typename TImage::Pointer input, RealType distthresh )
  RealType distthresh = 1.5;
  typename ImageType::Pointer wmgrow = Morphological<ImageType>(wmb, 0, true);
  typename ImageType::Pointer bsurf = LabelSurface<ImageType>(1, 1, wmgrow, distthresh); // or wmb ?
  typename ImageType::Pointer speedprior = NULL;
  WriteImage<ImageType>(bsurf, "surf.nii.gz");
  //    typename RealTypeImageType::Pointer distfromboundary =
  //  typename ImageType::Pointer surf=MaurerDistanceMap<ImageType>(0.5,1.e9,bsurf);
  // surf= SmoothImage<ImageType>(surf,3);
  typename ImageType::Pointer finalthickimage = BinaryThreshold<ImageType>(3, 3, 1, segmentationimage); // fixme

  gmb = BinaryThreshold<ImageType>(2, 3, 1, segmentationimage);  // fixme
  typename ImageType::Pointer gmgrow = Morphological<ImageType>(gmb, 1, true);
  typename ImageType::Pointer gmsurf = LabelSurface<ImageType>(1, 1, gmgrow, distthresh); // or wmb ?
  //  WriteImage<ImageType>(gmsurf,"surfdefgm.nii.gz");
  //  WriteImage<ImageType>(bsurf,"surfdefwm.nii.gz");

  typedef   DisplacementFieldType
    TimeVaryingVelocityFieldType;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType>                          FieldIterator;
  typedef typename DisplacementFieldType::IndexType                                         DIndexType;
  typedef typename DisplacementFieldType::PointType                                         DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType                                  VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType                                  VPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, RealType> DefaultInterpolatorType;
  typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, RealType>        DefaultInterpolatorType2;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  vinterp->SetInputImage(lapgrad);
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer ginterp =  ScalarInterpolatorType::New();
  typename ScalarInterpolatorType::Pointer winterp =  ScalarInterpolatorType::New();
  winterp->SetInputImage(wm);
  ginterp->SetInputImage(gm);

  DPointType pointIn1;
  DPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;

  typedef itk::ImageRegionIteratorWithIndex<ImageType>             IteratorType;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> VIteratorType;
  VIteratorType VIterator( lapgrad, lapgrad->GetLargestPossibleRegion().GetSize() );
  VIterator.GoToBegin();
  while(  !VIterator.IsAtEnd()  )
    {
    // the velocity field solution value
    VectorType vec = VIterator.Get();
    RealType   mag = 0;
    for( unsigned dd = 0; dd < ImageDimension; dd++ )
      {
      mag += vec[dd] * vec[dd];
      }
    mag = sqrt(mag);
    //      if (mag > 0) vec=vec/mag;
    VIterator.Set( (vec) * gradstep);
    ++VIterator;
    }

  //  m_MFR->SmoothDisplacementFieldGauss(lapgrad,1.7);
  std::cout << " Scaling done " << std::endl;

  typename ImageType::Pointer thickimage = laplacian;
  VectorType disp;
  VectorType incdisp;
  disp.Fill(0.0);
  incdisp.Fill(0.0);
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  RealType     totalerr = 1.e8, lasterr = 1.e10;
  unsigned     its = 0;
  wmgrow->FillBuffer(0);
  RealType      dmag = 0;
  RealType      thicknesserror = 0;
  unsigned long thickerrct = 0;
  unsigned int  badct = 0;
  RealType      thickoffset = 0;
  bool          checknans = true;

  while( its < alltheits &&  badct < 4 )
    {
    its++;
    if( totalerr > lasterr )
      {
      badct++; std::cout << " badct " << badct << std::endl;
      }
    else
      {
      badct = 0;
      }
    lasterr = totalerr;
    // Optimization Error initialized for this iteration
    totalerr = 0;
    incrfield->FillBuffer(zero);
    incrfield->FillBuffer(zero);
    incrinvfield->FillBuffer(zero);

    // generate phi
    //      corrfield->FillBuffer(zero);
    invfield->FillBuffer(zero);
    unsigned int ttiter = 0;

    thickimage->FillBuffer(0);
    hitimage->FillBuffer(0);
    totalimage->FillBuffer(0);
    thicknesserror = 0;
    thickerrct = 1;
    bool debug = false;
    bool spatprior = false;
    typename ImageType::Pointer priorim = NULL;
    if( speedprior )
      {
      spatprior = true;
      priorim = speedprior;
      }
    typename ImageType::Pointer wpriorim = NULL;
    RealType origthickprior = thickprior;

    while( ttiter < numtimepoints )    // N time integration points
      {
      //      m_MFR->Compose(incrinvfield,invfield,NULL);
      m_MFR->ComposeDiffs(invfield, incrinvfield, invfield, 1);

      if( debug )
        {
        std::cout << " exp " << std::endl;
        }
      // Integrate the negative velocity field to generate diffeomorphism corrfield step 3(a)
      //      corrfield=ExpDiffMap<ImageType,DisplacementFieldType>( velofield,  wm, -1, numtimepoints-ttiter);
      //      std::cout  << " corrf len " << m_MFR->MeasureDeformation( corrfield ) << std::endl;
      if( debug )
        {
        std::cout << " gmdef " << std::endl;
        }
      typename ImageType::Pointer gmdef = gm; // m_MFR->WarpImageBackward(gm,corrfield);
      totalerr = 0;

      typename ImageType::Pointer surfdef = m_MFR->WarpImageBackward(wm, invfield);
      if( debug )
        {
        std::cout << " thkdef " << std::endl;
        }
      typename ImageType::Pointer thkdef = m_MFR->WarpImageBackward(thickimage, invfield);
      if( debug )
        {
        std::cout << " thindef " << std::endl;
        }
      typename ImageType::Pointer thindef = m_MFR->WarpImageBackward(bsurf, invfield);
      if( spatprior )
        {
        wpriorim = m_MFR->WarpImageBackward(priorim, invfield);
        }

      typedef DisplacementFieldType GradientImageType;
      typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
        GradientImageFilterType;
      typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
      GradientImageFilterPointer gfilter = GradientImageFilterType::New();
      gfilter->SetInput(  surfdef );
      gfilter->SetSigma( smoothingsigma );
      gfilter->Update();
      typename DisplacementFieldType::Pointer   lapgrad2 = gfilter->GetOutput();

      // this is the "speed" image
      typename ImageType::Pointer speed_image = CopyImage<ImageType, DisplacementFieldType>(invfield);
      IteratorType xxIterator( speed_image, speed_image->GetLargestPossibleRegion().GetSize() );
      xxIterator.GoToBegin();
      RealType maxlapgrad2mag = 0;
      while(  !xxIterator.IsAtEnd()  )
        {
        typename ImageType::IndexType speedindex = xxIterator.GetIndex();
        if( segmentationimage->GetPixel(speedindex) == 2 )  // fixme
          {
          thickprior = origthickprior;
          VectorType wgradval = lapgrad2->GetPixel(speedindex);
          RealType   wmag = 0;
          for( unsigned kq = 0; kq < ImageDimension; kq++ )
            {
            wmag += wgradval[kq] * wgradval[kq];
            }
          if( fabs(wmag) < 1.e-6 )
            {
            wmag = 0;
            }
          wmag = sqrt(wmag);
          if( checknans )
            {
            if( vnl_math_isnan(wmag) || vnl_math_isinf(wmag) || wmag == 0 )
              {
              wgradval.Fill(0);
              lapgrad2->SetPixel(speedindex, wgradval);
              wmag = 0;
              }
            else
              {
              lapgrad2->SetPixel(speedindex, wgradval / wmag);
              }
            }
          totalerr += fabs(surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex) );
//          RealType thkval=thkdef->GetPixel(speedindex);
//          RealType thkval=finalthickimage->GetPixel(speedindex);
//          RealType fval=1; //(thickprior-thkval);
          //          if ( fval > 0 ) fval=1; else fval=-1;
// speed function here IMPORTANT!!
          RealType dd = (surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex) ) * gradstep;
          dd *= gm->GetPixel(speedindex);
          if( checknans )
            {
            if( vnl_math_isnan(dd) || vnl_math_isinf(dd) )
              {
              dd = 0;
              }
            }
          speed_image->SetPixel(speedindex, dd);
          if( wmag * dd > maxlapgrad2mag )
            {
            maxlapgrad2mag = wmag * dd;
            }
          }
        else
          {
          speed_image->SetPixel(speedindex, 0);
          }
        ++xxIterator;
        }

      if( maxlapgrad2mag < 1.e-4 )
        {
        maxlapgrad2mag = 1.e9;
        }
      if( ttiter == numtimepoints - 1 )
        {
        if( ImageDimension == 2 )
          {
          WriteImage<ImageType>(surfdef, "surfdef.nii.gz");
          }
        if( ImageDimension == 2 )
          {
          WriteImage<ImageType>(thindef, "thindef.nii.gz");
          }
        if( ImageDimension == 2 )
          {
          WriteImage<ImageType>(gmdef, "gmdef.nii.gz");
          }
        if( ImageDimension == 2 )
          {
          WriteImage<ImageType>(thkdef, "thick2.nii.gz");
          }
        }
      /* Now that we have the gradient image, we need to visit each voxel and compute objective function */
      //      std::cout << " maxlapgrad2mag " << maxlapgrad2mag << std::endl;
      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        typename DisplacementFieldType::IndexType velind = Iterator.GetIndex();
        VectorType wgradval = lapgrad2->GetPixel(velind); // *5.0/(maxlapgrad2mag*(RealType)numtimepoints);
        disp = wgradval * speed_image->GetPixel(velind);
        incrfield->SetPixel(velind, incrfield->GetPixel(velind) + disp);

        if( ttiter == 0 )  // make euclidean distance image
          {
          dmag = 0;
          disp = corrfield->GetPixel(velind);
          for( unsigned int jj = 0; jj < ImageDimension; jj++ )
            {
            dmag += disp[jj] * disp[jj];
            }
          RealType bval = bsurf->GetPixel(velind);
          if( checknans )
            {
            if( vnl_math_isnan(dmag) || vnl_math_isinf(dmag) )
              {
              dmag = 0;
              }
            if( vnl_math_isnan(bval) || vnl_math_isinf(bval) )
              {
              bval = 0;
              }
            }
          /** Change 2-26-2010 = incoporate gm prob in length ... */
          dmag = sqrt(dmag) * bval; // *gm->GetPixel(velind);
          thickimage->SetPixel(velind, dmag);
          totalimage->SetPixel(velind, dmag);
          hitimage->SetPixel(velind, bval);
          }
        else if( segmentationimage->GetPixel(velind) == 2 )     // fixme
          {
          RealType thkval = thkdef->GetPixel(velind);
          RealType putval = thindef->GetPixel(velind);
          hitimage->SetPixel(velind, hitimage->GetPixel(velind) + putval);
          totalimage->SetPixel(velind, totalimage->GetPixel(velind) + thkval);
          }

        ++Iterator;
        }

      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        IndexType velind = Iterator.GetIndex();
        bool      shouldbezero = false;
        if( segmentationimage->GetPixel(velind) == 0 )
          {
          shouldbezero = true;
          }
        if( !shouldbezero )
          {
          if( bsurf->GetPixel(velind) == 0 && gmsurf->GetPixel(velind) == 0 && segmentationimage->GetPixel(velind) !=
              2 )
            {
            shouldbezero = true;
            }
          }
        if( shouldbezero   )
          {
          velofield->SetPixel(velind, zero);
          corrfield->SetPixel(velind, zero);
          invfield->SetPixel(velind, zero);
          }
        incrinvfield->SetPixel(velind, velofield->GetPixel(velind) );
        ++Iterator;
        }

      if( ttiter == 0 )
        {
        corrfield->FillBuffer(zero);
        }

      InvertField<ImageType, DisplacementFieldType>( invfield, corrfield, 1.0, 0.1, 20, true);
      InvertField<ImageType, DisplacementFieldType>( corrfield, invfield, 1.0, 0.1, 20, true);
      //      InvertField<ImageType,DisplacementFieldType>( invfield, corrfield, 1.0,0.1,20,true);
      //  InvertField<ImageType,DisplacementFieldType>( corrfield, invfield, 1.0,0.1,20,true);
      ttiter++;
      }

    Iterator.GoToBegin();
    RealType maxth = 0;
    while(  !Iterator.IsAtEnd()  )
      {
      typename DisplacementFieldType::IndexType velind = Iterator.GetIndex();
      // increment velocity field at every voxel v = v + u, step 4
      velofield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() )
                          + incrfield->GetPixel(Iterator.GetIndex() ) );
      RealType hitval = hitimage->GetPixel(velind);
      RealType thkval = 0;
      if( hitval > 0.001 )    /** potential source of problem 2 -- this value could be smaller ... */
        {
        thkval = totalimage->GetPixel(velind) / hitval - thickoffset;
        }
      if( thkval > 10 )
        {
        std::cout << "thkval " << thkval << " hitval " << hitval << " total " << totalimage->GetPixel(velind)
                  << std::endl;
        }
      if( thkval < 0 )
        {
        thkval = 0;
        }
      if( segmentationimage->GetPixel(velind) == 2 )
        {
        finalthickimage->SetPixel(velind, thkval);
        }
      else
        {
        finalthickimage->SetPixel(velind, 0);
        }
      if( thkval > maxth )
        {
        maxth = thkval;
        }
      if( finalthickimage->GetPixel(velind)  > thickprior )
        {
        finalthickimage->SetPixel(velind, thickprior );
        }
      ++Iterator;
      }

    if( debug )
      {
      std::cout << " now smooth " << std::endl;
      }
    m_MFR->SmoothDisplacementFieldGauss(velofield, smoothingsigma);
    WriteImage<DisplacementFieldType>(corrfield, "corrfield.nii.gz");
    WriteImage<DisplacementFieldType>(invfield, "invfield.nii.gz");

    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DisplacementFieldType>(velofield,velofieldname.c_str());
    // std::string incrfieldname = outname + "incrfield";
    // WriteDisplacementField<DisplacementFieldType>(incrfield,incrfieldname.c_str());

    // std::string tname = outname + "dork1.nii.gz";
    // WriteImage<ImageType>(hitimage,tname.c_str());
    // tname = outname + "dork2.nii.gz";
    // WriteImage<ImageType>(totalimage,tname.c_str());
    if( thickerrct == 0 )
      {
      thickerrct = 1;
      }
    std::cout << " error " << totalerr << " at it " << its  << " th-err " << thicknesserror / (RealType)thickerrct
              << " max thick " << maxth << std::endl;
//    std::string sulcthickname =outname + "sulcthick.nii";
    //    if (ImageDimension==2) WriteJpg<ImageType>(finalthickimage,"thick.jpg");
    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DisplacementFieldType>(velofield,velofieldname.c_str());
    if( debug )
      {
      std::cout << "outside it " << its << std::endl;
      }
    // std::cin.get();
    }

  finalthickimage->SetDirection(omat);
  WriteImage<ImageType>(finalthickimage, outname.c_str() );
  finalthickimage->SetDirection(fmat);

  return 0;

  thickimage->FillBuffer(0);
  typename ImageType::Pointer thkdef = m_MFR->WarpImageBackward(finalthickimage, invfield);
  Iterator.GoToBegin();
  while(  !Iterator.IsAtEnd()  )
    {
    RealType tt1 = finalthickimage->GetPixel(Iterator.GetIndex() );
    RealType tt = thkdef->GetPixel(Iterator.GetIndex() );
    if( tt1 > tt )
      {
      tt = tt1;
      }
    thickimage->SetPixel(Iterator.GetIndex(), tt);
    ++Iterator;
    }

  std::string thickname = outname;
  thickimage->SetDirection(omat);
  WriteImage<ImageType>(thickimage, thickname.c_str() );

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 6 )
    {
    std::cout << "Usage:   " << argv[0]
              <<
      " ImageDimension Segmentation.nii.gz WMProb.nii.gz GMProb.nii.gz   Out.nii {GradStep-1-2D,2-3D}   {#Its-~50}  {ThickPriorValue-6} {Bool-use-curvature-prior} {smoothing} {BoolUseEuclidean?}"
              << std::endl;
    std::cout << " this is a kind of binary image registration thing with diffeomorphisms " << std::endl;
    std::cout
      << " Segmentation.nii.gz -- should contain the value 3 where WM exists and the value 2 where GM exists "
      << std::endl;
    return 1;
    }

  unsigned int dim = atoi(argv[1]);
  std::cout << " dim " << dim << std::endl;

  switch( dim )
    {
    case 2:
      {
      LaplacianThicknessExpDiff2<2>(argc, argv);
      }
      break;
    case 3:
      {
      LaplacianThicknessExpDiff2<3>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return EXIT_SUCCESS;

  return 1;
}
