// #include "DoSomethingToImage.cxx"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"

#include "itkWarpImageFilter.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkDeformationFieldFromMultiTransformFilter.h"

#include "itkImageFileWriter.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkImageToHistogramGenerator.h"
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
  if( (float) replaceval == (float) -1 )
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
SmoothImage(typename TImage::Pointer image, float sig)
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
SmoothDeformation(typename TImage::Pointer vectorimage, float sig)
{
  enum { ImageDimension = TImage::ImageDimension };

  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<float, ImageDimension>  ImageType;
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

template <class TImage, class TDeformationField>
typename TImage::Pointer
ComputeJacobian(TDeformationField* field )
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  //  unsigned int row=0;
  // unsigned int col=0;
  typedef itk::Image<float, ImageDimension> FloatImageType;
  typename FloatImageType::RegionType m_JacobianRegion;

  typename TImage::SizeType s = field->GetLargestPossibleRegion().GetSize();
  typename TImage::SpacingType sp = field->GetSpacing();

  typename FloatImageType::Pointer m_FloatImage = NULL;
  m_FloatImage = FloatImageType::New();
  m_FloatImage->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  m_FloatImage->SetBufferedRegion( field->GetLargestPossibleRegion().GetSize() );
  m_FloatImage->SetSpacing(field->GetSpacing() );
  m_FloatImage->SetDirection( field->GetDirection() );
  m_FloatImage->SetOrigin(field->GetOrigin() );
  m_FloatImage->Allocate();
  m_FloatImage->FillBuffer(0);

  typename FloatImageType::SizeType m_FieldSize = field->GetLargestPossibleRegion().GetSize();

  typedef itk::ImageRegionIteratorWithIndex<FloatImageType> Iterator;
  Iterator wimIter( m_FloatImage, m_FloatImage->GetLargestPossibleRegion()  );
  wimIter.GoToBegin();
  for( ; !wimIter.IsAtEnd(); ++wimIter )
    {
    wimIter.Set(1.0);
    }

  typedef  vnl_matrix<double> MatrixType;
  MatrixType jMatrix, idMatrix, avgMatrix;
  jMatrix.set_size(ImageDimension, ImageDimension);
  avgMatrix.set_size(ImageDimension, ImageDimension);
  avgMatrix.fill(0);
  itk::ImageRegionIteratorWithIndex<TDeformationField>
  m_FieldIter( field, field->GetLargestPossibleRegion() );
  typename TImage::IndexType rindex;
  typename TImage::IndexType ddrindex;
  typename TImage::IndexType ddlindex;

  typename TImage::IndexType difIndex[ImageDimension][2];

  double       det = 0.0;
  unsigned int posoff = 1;
  float        difspace = 1.0;
  float        space = 1.0;
  if( posoff == 0 )
    {
    difspace = 1.0;
    }

  typedef itk::Vector<float, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension> FieldType;

  typename FieldType::PixelType dPix;
  typename FieldType::PixelType lpix;
  typename FieldType::PixelType llpix;
  typename FieldType::PixelType rpix;
  typename FieldType::PixelType rrpix;
  typename FieldType::PixelType cpix;

  float volumeelt = 1.0;
  for( int j = 0; j < ImageDimension; j++ )
    {
    volumeelt *= sp[j];
    }
  //  double totaljac=0.0;

  // /the finite difference equations
  float wC, wLL, wL, wR, wRR;
  // 3rd deriv - 4th order
  wC = 0.0;
  wLL = 1.; wL = -2.0; wR =  2.0; wRR = -1.0;
  // 4th deriv - 4th order
  wC = -6.0;
  wLL = 1.; wL = -4.0; wR = -4.0; wRR = 1.0;
  // 2nd deriv - 4th order
  wC = 30.0;
  wLL = -1.0; wL = 16.0; wR = 16.0; wRR = -1.0;
  float total = wC; // wLL + wL + wR + wRR;
  if( total == 0.0 )
    {
    total = 1.0;
    }

  unsigned long ct = 0;
  for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
    {
    rindex = m_FieldIter.GetIndex();
    float mindist = 1.0;
    bool  oktosample = true;
    float dist = 100.0;
    for( unsigned int row = 0; row < ImageDimension; row++ )
      {
      dist = fabs( (float)rindex[row]);
      if( dist < mindist )
        {
        oktosample = false;
        }
      dist = fabs( (float)s[row] - (float)rindex[row]);
      if( dist < mindist )
        {
        oktosample = false;
        }
      }
    if( oktosample )
      {
      ct++;
      typename TImage::IndexType temp = rindex;
      cpix = field->GetPixel(rindex);
      for( unsigned int row = 0; row < ImageDimension; row++ )
        {
        difIndex[row][0] = rindex;
        difIndex[row][1] = rindex;
        ddrindex = rindex;
        ddlindex = rindex;
        if( (int) rindex[row] < (int)(m_FieldSize[row] - 2) )
          {
          difIndex[row][0][row] = rindex[row] + posoff;
          ddrindex[row] = rindex[row] + posoff * 2;
          }
        if( rindex[row] > 1 )
          {
          difIndex[row][1][row] = rindex[row] - 1;
          ddlindex[row] = rindex[row] - 2;
          }

        float h = 0.5;
        space = 1.0; // should use image spacing here?

        rpix = field->GetPixel(difIndex[row][1]);
        rpix = rpix * h + cpix * (1. - h);
        lpix = field->GetPixel(difIndex[row][0]);
        lpix = lpix * h + cpix * (1. - h);
        //    dPix = ( rpix - lpix)*(1.0)/(2.0);

        rrpix = field->GetPixel(ddrindex);
        rrpix = rrpix * h + rpix * (1. - h);
        llpix = field->GetPixel(ddlindex);
        llpix = llpix * h + lpix * (1. - h);
        dPix = ( lpix * (-8.0) + rpix * 8.0 - rrpix + llpix ) * (-1.0) * space / (12.0); // 4th order centered
                                                                                         // difference
        // dPix=( lpix - rpix )*(1.0)*space/(2.0*h); //4th order centered difference
        for( unsigned int col = 0; col < ImageDimension; col++ )
          {
          float val;
          if( row == col )
            {
            val = dPix[col] / sp[col] + 1.0;
            }
          else
            {
            val = dPix[col] / sp[col];
            }
          //        std::cout << " row " << row << " col " << col << " val " << val << std::endl;
          jMatrix.put(col, row, val);
          avgMatrix.put(col, row, avgMatrix.get(col, row) + val);
          }
        }

      // the determinant of the jacobian matrix
      // std::cout << " get det " << std::endl;
      det = vnl_determinant(jMatrix);
      //    float prodval = m_FloatImage->GetPixel(rindex);
      if( det < 0.0 )
        {
        det = 0;
        }

      m_FloatImage->SetPixel(rindex,  det );

      // totaljac+=det;
      } // oktosample if
    }

  return m_FloatImage;
}

template <class TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground,
             typename TImage::PixelType newval, typename TImage::Pointer input, float distthresh )
{
  std::cout << " Label Surf " << std::endl;
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typename   ImageType::Pointer     Image = ImageType::New();
  Image->SetLargestPossibleRegion(input->GetLargestPossibleRegion()  );
  Image->SetBufferedRegion(input->GetLargestPossibleRegion() );
  Image->Allocate();
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
        float dist = 0.0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          dist += (float)(ind[j] - ind2[j]) * (float)(ind[j] - ind2[j]);
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
typename TImage::Pointer
SpeedPrior(typename TImage::Pointer image1, typename TImage::Pointer  wmimage,  typename TImage::Pointer surf )
{
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };

  typedef itk::SurfaceImageCurvature<ImageType> ParamType;
  typename ParamType::Pointer Parameterizer = ParamType::New();

  float sig = 1.5;

  Parameterizer->SetInput(wmimage);
  Parameterizer->SetNeighborhoodRadius( 1. );
  Parameterizer->SetSigma(sig);

  Parameterizer->SetUseLabel(false);
  Parameterizer->SetUseGeodesicNeighborhood(false);
  float sign = 1.0;
  Parameterizer->SetkSign(sign);
  Parameterizer->SetThreshold(0);
  Parameterizer->ComputeFrameOverDomain( 3 );

  typename ImageType::Pointer outimage = Parameterizer->GetFunctionImage();

  float         max = 0;
  float         min = 1.e9;
  float         mean = 0.0;
  unsigned long ct = 0;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator iter( outimage,  outimage->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    float pix = iter.Get();
    mean += pix;
    ct++;
    if( pix > max )
      {
      max = pix;
      }
    if( pix < min )
      {
      min = pix;
      }
    }
  mean /= (float)ct;
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    float pix = iter.Get();
    pix = (pix - min) / (max - min);
    iter.Set(pix);
    }

  return outimage;
}

template <class TImage>
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
LaplacianGrad(typename TImage::Pointer wm, typename TImage::Pointer gm, float sig)
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
ExpDiffMap(typename TField::Pointer velofield,  typename TImage::Pointer wm,  float sign, unsigned int numtimepoints )
{
  typedef TImage ImageType;
  typedef TField DeformationFieldType;
  typename TField::PixelType zero, disp;
  enum { ImageDimension = TImage::ImageDimension };
  disp.Fill(0);
  zero.Fill(0);

  typename DeformationFieldType::Pointer incrfield = DeformationFieldType::New();
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
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>             AffineTransformType;
  typedef itk::DeformationFieldFromMultiTransformFilter<TField, TField, AffineTransformType> WarperType;
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
    warper->PushBackDeformationFieldTransform(incrfield);
    }

  warper->Update();
  return warper->GetOutput();
}

template <class TImage, class TField>
typename TField::Pointer
DiReCTCompose(typename TField::Pointer velofield, typename TField::Pointer diffmap )
{
  typedef TImage ImageType;
  typedef TField DeformationFieldType;
  typename TField::PixelType zero, disp;
  enum { ImageDimension = TImage::ImageDimension };
  disp.Fill(0);
  zero.Fill(0);

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>             AffineTransformType;
  typedef itk::DeformationFieldFromMultiTransformFilter<TField, TField, AffineTransformType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetOutputSize(velofield->GetLargestPossibleRegion().GetSize() );
  warper->SetOutputSpacing(velofield->GetSpacing() );
  warper->SetOutputOrigin(velofield->GetOrigin() );
  warper->SetOutputDirection(velofield->GetDirection() );
  warper->DetermineFirstDeformNoInterp();
  warper->PushBackDeformationFieldTransform(diffmap);
  warper->PushBackDeformationFieldTransform(velofield);
  warper->Update();
  return warper->GetOutput();
}

template <unsigned int ImageDimension>
int LaplacianThicknessExpDiff(int argc, char *argv[])
{
  int         argct = 2;
  std::string wfn = std::string(argv[argct]); argct++;

  std::string  gfn = std::string(argv[argct]); argct++;
  std::string  outname = std::string(argv[argct]); argct++;
  unsigned int numtimepoints = 10;
  float        gradstep = (float)(-1.0) * 0.5; // (ImageDimension-1);
  if( argc > argct )
    {
    gradstep = atof(argv[argct]) * (-1.0) * 1.0 / (float)numtimepoints;
    }
  argct++;
  unsigned int alltheits = 50;
  if( argc > argct )
    {
    alltheits = atoi(argv[argct]);
    }
  argct++;
  float thickprior = 6.0;
  if( argc > argct )
    {
    thickprior = atof(argv[argct]);
    }
  argct++;
  bool useCurvaturePrior = false;
  if( argc > argct )
    {
    useCurvaturePrior = atoi(argv[argct]);
    }
  argct++;
  float smoothingsigma = 1;
  if( argc > argct )
    {
    smoothingsigma = atof(argv[argct]);
    }
  argct++;
  bool useEuclidean = true;
  if( argc > argct )
    {
    useEuclidean = atoi(argv[argct]);
    }
  argct++;
  std::cout << " smooth " << smoothingsigma << " thp " << thickprior << " gs " << gradstep << std::endl;
  typedef float                                                      PixelType;
  typedef itk::Vector<float, ImageDimension>                         VectorType;
  typedef itk::Image<VectorType, ImageDimension>                     DeformationFieldType;
  typedef itk::Image<PixelType, ImageDimension>                      ImageType;
  typedef itk::ImageFileReader<ImageType>                            readertype;
  typedef itk::ImageFileWriter<ImageType>                            writertype;
  typedef typename  ImageType::IndexType                             IndexType;
  typedef typename  ImageType::SizeType                              SizeType;
  typedef typename  ImageType::SpacingType                           SpacingType;
  typedef itk::Image<VectorType, ImageDimension + 1>                 tvt;
  typedef itk::ANTSImageRegistrationOptimizer<ImageDimension, float> ROType;
  typename ROType::Pointer m_MFR = ROType::New();

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
  SpacingType spacing = wm->GetSpacing();
  typename DeformationFieldType::Pointer lapgrad;
  typename ImageType::Pointer gmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, gm);
  typename ImageType::Pointer wmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, wm);
  typename ImageType::Pointer laplacian = SmoothImage<ImageType>(wm, smoothingsigma);
  lapgrad = LaplacianGrad<ImageType, DeformationFieldType>(wmb, gmb, 1);

  typename DeformationFieldType::Pointer corrfield = DeformationFieldType::New();
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
  typename DeformationFieldType::Pointer incrfield = DeformationFieldType::New();
  incrfield->SetSpacing( wm->GetSpacing() );
  incrfield->SetOrigin( wm->GetOrigin() );
  incrfield->SetDirection( wm->GetDirection() );
  incrfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrfield->Allocate();
  incrfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer invfield = DeformationFieldType::New();
  invfield->SetSpacing( wm->GetSpacing() );
  invfield->SetOrigin( wm->GetOrigin() );
  invfield->SetDirection( wm->GetDirection() );
  invfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  invfield->SetRequestedRegion(wm->GetRequestedRegion() );
  invfield->SetBufferedRegion( wm->GetBufferedRegion() );
  invfield->Allocate();
  invfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer incrinvfield = DeformationFieldType::New();
  incrinvfield->SetSpacing( wm->GetSpacing() );
  incrinvfield->SetOrigin( wm->GetOrigin() );
  incrinvfield->SetDirection( wm->GetDirection() );
  incrinvfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrinvfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrinvfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrinvfield->Allocate();
  incrinvfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer velofield = DeformationFieldType::New();
  velofield->SetSpacing( wm->GetSpacing() );
  velofield->SetOrigin( wm->GetOrigin() );
  velofield->SetDirection( wm->GetDirection() );
  velofield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  velofield->SetRequestedRegion(wm->GetRequestedRegion() );
  velofield->SetBufferedRegion( wm->GetBufferedRegion() );
  velofield->Allocate();
  velofield->FillBuffer(zero);

  //  LabelSurface(typename TImage::PixelType foreground,
  //       typename TImage::PixelType newval, typename TImage::Pointer input, float distthresh )
  float distthresh = 1.1;
  typename ImageType::Pointer wmgrow = Morphological<ImageType>(wmb, 1, true);
  typename ImageType::Pointer bsurf = LabelSurface<ImageType>(1, 1, wmgrow, distthresh);
  typename ImageType::Pointer speedprior = NULL;
  if(  useCurvaturePrior )
    {
    speedprior = SpeedPrior<ImageType>(gm, wm, bsurf);
    }
  // WriteImage<ImageType>(bsurf,"surf.hdr");
  //	typename DoubleImageType::Pointer distfromboundary =
  //  typename ImageType::Pointer surf=MaurerDistanceMap<ImageType>(0.5,1.e9,bsurf);
  // surf= SmoothImage<ImageType>(surf,3);
  typename ImageType::Pointer finalthickimage = BinaryThreshold<ImageType>(0.5, 1.e9, 1, wm);

  // typename ImageType::Pointer gmsurf = LabelSurface<ImageType>(1,1,gmb, distthresh);
  //  gmsurf=MaurerDistanceMap<ImageType>(0.5,1.e9,gmsurf);
  //  gmsurf= SmoothImage<ImageType>(gmsurf,3);

  typename ImageType::SizeType s = wm->GetLargestPossibleRegion().GetSize();
  typename DeformationFieldType::IndexType velind;
  typedef   DeformationFieldType                                                         TimeVaryingVelocityFieldType;
  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType>                        FieldIterator;
  typedef typename DeformationFieldType::IndexType                                       DIndexType;
  typedef typename DeformationFieldType::PointType                                       DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType                               VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType                               VPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float> DefaultInterpolatorType;
  typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, float>         DefaultInterpolatorType2;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  vinterp->SetInputImage(lapgrad);
  typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer ginterp =  ScalarInterpolatorType::New();
  typename ScalarInterpolatorType::Pointer winterp =  ScalarInterpolatorType::New();
  winterp->SetInputImage(wm);
  ginterp->SetInputImage(gm);

  DPointType pointIn1;
  DPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;

  typename ImageType::Pointer surfdef;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>            IteratorType;
  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> VIteratorType;
  VIteratorType VIterator( lapgrad, lapgrad->GetLargestPossibleRegion().GetSize() );
  VIterator.GoToBegin();
  while(  !VIterator.IsAtEnd()  )
    {
    VectorType vec = VIterator.Get();
    float      mag = 0;
    for( unsigned dd = 0; dd < ImageDimension; dd++ )
      {
      mag += vec[dd] * vec[dd];
      }
    mag = sqrt(mag);
    if( mag > 0 )
      {
      vec = vec / mag;
      }
    VIterator.Set( (vec) * gradstep);
    ++VIterator;
    }

  //  m_MFR->SmoothDeformationFieldGauss(lapgrad,1.7);
  std::cout << " Scaling done " << std::endl;

  //  float thislength=0;
  unsigned long ct = 1;
  bool          timedone = false;

  typename ImageType::Pointer thickimage = laplacian;
  VectorType disp;
  VectorType incdisp;
  disp.Fill(0.0);
  incdisp.Fill(0.0);
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  timedone = false;
  float    totalerr = 1.e8, lasterr = 1.e10;
  unsigned its = 0;
  wmgrow->FillBuffer(0);
  float         dmag = 0;
  float         thicknesserror = 0;
  unsigned long thickerrct = 0;
  unsigned int  badct = 0;
  float         thickoffset = 0;

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
    ct = 0;
    incrfield->FillBuffer(zero);
    incrfield->FillBuffer(zero);
    incrinvfield->FillBuffer(zero);

    // generate phi
    corrfield->FillBuffer(zero);
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
    float origthickprior = thickprior;

    while( ttiter < numtimepoints )    // N time integration points
      {
      //  void ComposeDiffs(DeformationFieldPointer fieldtowarpby, DeformationFieldPointer field,
      // DeformationFieldPointer fieldout, float sign);
      m_MFR->ComposeDiffs(invfield, incrinvfield, invfield, 1);

      if( debug )
        {
        std::cout << " exp " << std::endl;
        }
      // Integrate the negative velocity field to generate diffeomorphism corrfield step 3(a)
      //	  ExpDiffMap<ImageType,DeformationFieldType>( velofield, corrfield, wm, -1, numtimepoints-ttiter);
      corrfield = ExpDiffMap<ImageType, DeformationFieldType>( velofield,  wm, -1, numtimepoints - ttiter);
      // why integrate velofield, but only compose incrinvfield ???
      // technically, we should warp the gm image by corrfield but this can be avoided
      totalerr = 0;
      //	  typename ImageType::Pointer gmdef =m_MFR->WarpImageBackward(gm,corrfield);
      typename ImageType::Pointer gmdef = m_MFR->WarpMultiTransform(gm, gm, NULL, corrfield, false, NULL );
      //	  typename ImageType::Pointer surfdef=m_MFR->WarpImageBackward(wm,invfield);
      typename ImageType::Pointer surfdef = m_MFR->WarpMultiTransform(wm, wm, NULL, invfield, false, NULL );
      //	  typename ImageType::Pointer thkdef =m_MFR->WarpImageBackward(thickimage,invfield);
      typename ImageType::Pointer thkdef = m_MFR->WarpMultiTransform(thickimage, thickimage, NULL, invfield, false,
                                                                     NULL );
      //	  typename ImageType::Pointer thindef =m_MFR->WarpImageBackward(bsurf,invfield);
      typename ImageType::Pointer thindef = m_MFR->WarpMultiTransform(bsurf, bsurf, NULL, invfield, false, NULL );
      //	  if (spatprior) wpriorim=m_MFR->WarpImageBackward(priorim,invfield);
      if( spatprior )
        {
        wpriorim = m_MFR->WarpMultiTransform(priorim, priorim, NULL, invfield, false, NULL );
        }

      typedef DeformationFieldType GradientImageType;
      typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
      GradientImageFilterType;
      typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
      GradientImageFilterPointer gfilter = GradientImageFilterType::New();
      gfilter->SetInput(  surfdef );
      gfilter->SetSigma( smoothingsigma );
      gfilter->Update();
      typename DeformationFieldType::Pointer   lapgrad2 = gfilter->GetOutput();

      GradientImageFilterPointer gfilter2 = GradientImageFilterType::New();
      gfilter2->SetInput(  gm );
      gfilter2->SetSigma( smoothingsigma );
      //	  gfilter2->Update();
      //	  typename DeformationFieldType::Pointer   lapgrad3=gfilter2->GetOutput();

/** the code below sets up the scalar "speed" function that multiplies
    the gradient driving the registration -- akin to "momentum" */
      typename ImageType::Pointer lapjac = ComputeJacobian<ImageType, DeformationFieldType>(invfield);
      IteratorType xxIterator( lapjac, lapjac->GetLargestPossibleRegion().GetSize() );
      xxIterator.GoToBegin();
      while(  !xxIterator.IsAtEnd()  )
        {
        typename ImageType::IndexType speedindex = xxIterator.GetIndex();
        if( gm->GetPixel(speedindex) >= 0.5 )
          {
//	      float thkval=thkdef->GetPixel(speedindex);
          float thkval = finalthickimage->GetPixel(speedindex);
          float prior = 1;
          if( spatprior )
            {
            float prval = wpriorim->GetPixel(speedindex);
            float partialvol = surfdef->GetPixel(speedindex);
            if( prval > 0.5 && partialvol > 1.e-3 )
              {
              prior = prval / partialvol;                           // 7;//0.5*origthickprior;// prval;
              }
            // if (prior > 1 ) prior=1;
            }
          // else thickprior = origthickprior;
          // } else
          thickprior = origthickprior;

          VectorType wgradval = lapgrad2->GetPixel(speedindex); // velofield->GetPixel(speedindex);//
          //	      VectorType ggradval=lapgrad3->GetPixel(speedindex);
          double dp = 0;
          double gmag = 0, wmag = 0;
          for( unsigned kq = 0; kq < ImageDimension; kq++ )
            {
            //	      gmag+= ggradval[kq]*ggradval[kq];
            wmag += wgradval[kq] * wgradval[kq];
            }
          if( fabs(wmag) < 1.e-9 )
            {
            wmag = 0;
            }
          if( fabs(gmag) < 1.e-9 )
            {
            gmag = 0;
            }
          gmag = sqrt(gmag);
          wmag = sqrt(wmag);
          //	      if (gmag > 0 && wmag > 0)
          //	{
          //	for (unsigned kq=0;kq<ImageDimension; kq++) dp+= ggradval[kq]/gmag*wgradval[kq]/wmag;
          //	}
//	      if (fabs(dp) < 0.6) dp=0;
//	      dp=fabs(dp);
          dp = 1.0; // -dp;
          // tempim->SetPixel(speedindex,dp);

          double fval = (thickprior - thkval);
          double sigmoidf = 1;
//	      double sigmoidf=1.0- thkval/thickprior;
          if( fval >= 0 )
            {
            sigmoidf = 1.0 / (1.0 + exp( (-1.0) * fval * 0.01) );
            }
          if( fval < 0 )
            {
            sigmoidf = -1.0 * (1.0 - thickprior / thkval);
            }
//	      else sigmoidf=1.0- thkval/thickprior;
          float thkscale = thickprior / thkval;
          if( thkscale < 0.99 )
            {
            thkscale = 0.99;               // *+thkscale*0.1;
            }
          if( fval < 0 )
            {
            velofield->SetPixel(speedindex, velofield->GetPixel(speedindex) * thkscale);
            }

          float dd = surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex);
          //	      float gmd=gmdef->GetPixel(speedindex);
          totalerr += fabs(dd);
          if( wm->GetPixel(speedindex) > 0.5 && bsurf->GetPixel(speedindex) < 0.5 )
            {
            dd = 0;
            }
          float stopval = gm->GetPixel(speedindex);
//	      if (wm->GetPixel(speedindex) > stopval) stopval=1;
          float jwt = xxIterator.Get();
          if( jwt < 1 )
            {
            jwt = 1;
            }
          else
            {
            jwt = 1.0 / xxIterator.Get();
            }
//	      dd*=stopval*jwt*thindef->GetPixel(speedindex)*sigmoidf*gradstep*dp*gmd*jwt;
          dd *= stopval * sigmoidf * gradstep * jwt * prior; // speed function here IMPORTANT!!
          lapjac->SetPixel(speedindex, dd);
          }
        else
          {
          lapjac->SetPixel(speedindex, 0);
          }
        ++xxIterator;
        }

/** smooth the momentum image */
      lapjac = SmoothImage<ImageType>(lapjac, 1);
      //	  if (ImageDimension==2) WriteJpg<ImageType>(surfdef,"surfdef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(thindef,"thindef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(gmdef,"gmdef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(lapjac,"diff.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(wpriorim,"prior.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(thkdef,"thick2.jpg");
//	  if (ImageDimension==2) WriteJpg<ImageType>(tempim,"dotp.jpg");
      // exit(0);

      /* Now that we have the gradient image, we need to visit each voxel and compute objective function */
      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        velind = Iterator.GetIndex();
        //	      float currentthickvalue=finalthickimage->GetPixel(velind);
        VectorType wgradval = lapgrad2->GetPixel(velind);

        disp = wgradval * lapjac->GetPixel(velind);

        incrfield->SetPixel(velind, incrfield->GetPixel(velind) + disp);

        if( ttiter == 0 ) // make euclidean distance image
          {
          dmag = 0;
          disp = corrfield->GetPixel(velind);
          for( unsigned int jj = 0; jj < ImageDimension; jj++ )
            {
            dmag += disp[jj] * disp[jj];
            }
          dmag = sqrt(dmag) * bsurf->GetPixel(velind);
          thickimage->SetPixel(velind, dmag);
          totalimage->SetPixel(velind, dmag);
          hitimage->SetPixel(velind, bsurf->GetPixel(velind) );
          }
        else if( gm->GetPixel(velind) >= 0.5 )
          {
          float thkval = thkdef->GetPixel(velind);
          float putval = thindef->GetPixel(velind);
          //		  float getval=hitimage->GetPixel(velind);
          hitimage->SetPixel(velind, hitimage->GetPixel(velind) + putval);
          totalimage->SetPixel(velind, totalimage->GetPixel(velind) + thkval);
          }
        ++Iterator;
        }

      //	  if (ttiter ==0)
      // WriteImage<ImageType>(totalimage,"totalimage.hdr");
      // WriteImage<ImageType>(hitimage,"hitimage.hdr");

      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        incrinvfield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() ) );
        ++Iterator;
        }

      ttiter++;
      }

    Iterator.GoToBegin();
    float maxth = 0;
    while(  !Iterator.IsAtEnd()  )
      {
      velind = Iterator.GetIndex();
      /* increment velocity field at every voxel v = v + u, step 4 */
      velofield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() )
                          + incrfield->GetPixel(Iterator.GetIndex() ) );
      float hitval = hitimage->GetPixel(velind);
      if( hitval == 0 )
        {
        hitval = 1;
        }
      float thkval = totalimage->GetPixel(velind) / hitval - thickoffset;
      if( thkval < 0 )
        {
        thkval = 0;
        }
      finalthickimage->SetPixel(velind, thkval);

      if( thkval > maxth )
        {
        maxth = thkval;
        }
      ++Iterator;
      }

    if( debug )
      {
      std::cout << " now smooth " << std::endl;
      }
    m_MFR->SmoothDeformationFieldGauss(velofield, smoothingsigma);

    if( thickerrct == 0 )
      {
      thickerrct = 1;
      }
    std::cout << " error " << totalerr << " at it " << its  << " th-err " << thicknesserror / (float)thickerrct
              << " max thick " << maxth << std::endl;
    finalthickimage->SetDirection(omat);
    WriteImage<ImageType>(finalthickimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int LaplacianThicknessExpDiff2(int argc, char *argv[])
{
  int         argct = 2;
  std::string wfn = std::string(argv[argct]); argct++;

  std::string  gfn = std::string(argv[argct]); argct++;
  std::string  outname = std::string(argv[argct]); argct++;
  unsigned int numtimepoints = 10;
  float        gradstep = (float)(-1.0) * 0.5; // (ImageDimension-1);
  if( argc > argct )
    {
    gradstep = atof(argv[argct]) * (-1.0);
    }
  gradstep *= 1.0 / (float)numtimepoints * 100;  argct++;
  unsigned int alltheits = 50;
  if( argc > argct )
    {
    alltheits = atoi(argv[argct]);
    }
  argct++;
  float thickprior = 6.0;
  if( argc > argct )
    {
    thickprior = atof(argv[argct]);
    }
  argct++;
  bool useCurvaturePrior = false;
  if( argc > argct )
    {
    useCurvaturePrior = atoi(argv[argct]);
    }
  argct++;
  float smoothingsigma = 1.5;
  if( argc > argct )
    {
    smoothingsigma = atof(argv[argct]);
    }
  argct++;
  bool useEuclidean = true;
  if( argc > argct )
    {
    useEuclidean = atoi(argv[argct]);
    }
  argct++;
  std::cout << " smooth " << smoothingsigma << " thp " << thickprior << " gs " << gradstep << std::endl;
  typedef float                                                      PixelType;
  typedef itk::Vector<float, ImageDimension>                         VectorType;
  typedef itk::Image<VectorType, ImageDimension>                     DeformationFieldType;
  typedef itk::Image<PixelType, ImageDimension>                      ImageType;
  typedef itk::ImageFileReader<ImageType>                            readertype;
  typedef itk::ImageFileWriter<ImageType>                            writertype;
  typedef typename  ImageType::IndexType                             IndexType;
  typedef typename  ImageType::SizeType                              SizeType;
  typedef typename  ImageType::SpacingType                           SpacingType;
  typedef itk::Image<VectorType, ImageDimension + 1>                 tvt;
  typedef itk::ANTSImageRegistrationOptimizer<ImageDimension, float> ROType;
  typename ROType::Pointer m_MFR = ROType::New();

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
  SpacingType spacing = wm->GetSpacing();
  typename DeformationFieldType::Pointer lapgrad;
  typename ImageType::Pointer gmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, gm);
  typename ImageType::Pointer wmb = BinaryThreshold<ImageType>(0.5, 1.e9, 1, wm);
  typename ImageType::Pointer laplacian = SmoothImage<ImageType>(wm, smoothingsigma);
  lapgrad = LaplacianGrad<ImageType, DeformationFieldType>(wmb, gmb, 1);

  typename DeformationFieldType::Pointer corrfield = DeformationFieldType::New();
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
  typename DeformationFieldType::Pointer incrfield = DeformationFieldType::New();
  incrfield->SetSpacing( wm->GetSpacing() );
  incrfield->SetOrigin( wm->GetOrigin() );
  incrfield->SetDirection( wm->GetDirection() );
  incrfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrfield->Allocate();
  incrfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer invfield = DeformationFieldType::New();
  invfield->SetSpacing( wm->GetSpacing() );
  invfield->SetOrigin( wm->GetOrigin() );
  invfield->SetDirection( wm->GetDirection() );
  invfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  invfield->SetRequestedRegion(wm->GetRequestedRegion() );
  invfield->SetBufferedRegion( wm->GetBufferedRegion() );
  invfield->Allocate();
  invfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer incrinvfield = DeformationFieldType::New();
  incrinvfield->SetSpacing( wm->GetSpacing() );
  incrinvfield->SetOrigin( wm->GetOrigin() );
  incrinvfield->SetDirection( wm->GetDirection() );
  incrinvfield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  incrinvfield->SetRequestedRegion(wm->GetRequestedRegion() );
  incrinvfield->SetBufferedRegion( wm->GetBufferedRegion() );
  incrinvfield->Allocate();
  incrinvfield->FillBuffer(zero);

  typename DeformationFieldType::Pointer velofield = DeformationFieldType::New();
  velofield->SetSpacing( wm->GetSpacing() );
  velofield->SetOrigin( wm->GetOrigin() );
  velofield->SetDirection( wm->GetDirection() );
  velofield->SetLargestPossibleRegion(wm->GetLargestPossibleRegion() );
  velofield->SetRequestedRegion(wm->GetRequestedRegion() );
  velofield->SetBufferedRegion( wm->GetBufferedRegion() );
  velofield->Allocate();
  velofield->FillBuffer(zero);

  //  LabelSurface(typename TImage::PixelType foreground,
  //       typename TImage::PixelType newval, typename TImage::Pointer input, float distthresh )
  float distthresh = 1.1;
  typename ImageType::Pointer wmgrow = Morphological<ImageType>(wmb, 1, true);
  typename ImageType::Pointer bsurf = LabelSurface<ImageType>(1, 1, wmgrow, distthresh); // or wmb ?
  typename ImageType::Pointer speedprior = NULL;
  if(  useCurvaturePrior )
    {
    speedprior = SpeedPrior<ImageType>(gm, wm, bsurf);
    }
  // WriteImage<ImageType>(bsurf,"surf.hdr");
  //	typename DoubleImageType::Pointer distfromboundary =
  //  typename ImageType::Pointer surf=MaurerDistanceMap<ImageType>(0.5,1.e9,bsurf);
  // surf= SmoothImage<ImageType>(surf,3);
  typename ImageType::Pointer finalthickimage = BinaryThreshold<ImageType>(0.5, 1.e9, 1, wm);

  // typename ImageType::Pointer gmsurf = LabelSurface<ImageType>(1,1,gmb, distthresh);
  //  gmsurf=MaurerDistanceMap<ImageType>(0.5,1.e9,gmsurf);
  //  gmsurf= SmoothImage<ImageType>(gmsurf,3);

  typename ImageType::SizeType s = wm->GetLargestPossibleRegion().GetSize();
  typename DeformationFieldType::IndexType velind;  velind.Fill(0);
  typedef   DeformationFieldType                                                         TimeVaryingVelocityFieldType;
  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType>                        FieldIterator;
  typedef typename DeformationFieldType::IndexType                                       DIndexType;
  typedef typename DeformationFieldType::PointType                                       DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType                               VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType                               VPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float> DefaultInterpolatorType;
  typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, float>         DefaultInterpolatorType2;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  vinterp->SetInputImage(lapgrad);
  typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer ginterp =  ScalarInterpolatorType::New();
  typename ScalarInterpolatorType::Pointer winterp =  ScalarInterpolatorType::New();
  winterp->SetInputImage(wm);
  ginterp->SetInputImage(gm);

  DPointType pointIn1;
  DPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;

  typename ImageType::Pointer surfdef;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>            IteratorType;
  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> VIteratorType;
  VIteratorType VIterator( lapgrad, lapgrad->GetLargestPossibleRegion().GetSize() );
  VIterator.GoToBegin();
  while(  !VIterator.IsAtEnd()  )
    {
    VectorType vec = VIterator.Get();
    float      mag = 0;
    for( unsigned dd = 0; dd < ImageDimension; dd++ )
      {
      mag += vec[dd] * vec[dd];
      }
    mag = sqrt(mag);
    if( mag > 0 )
      {
      vec = vec / mag;
      }
    VIterator.Set( (vec) * gradstep);
    ++VIterator;
    }

  //  m_MFR->SmoothDeformationFieldGauss(lapgrad,1.7);
  std::cout << " Scaling done " << std::endl;

  //  float thislength=0;
  unsigned long ct = 1;
  bool          timedone = false;

  typename ImageType::Pointer thickimage = laplacian;
  VectorType disp;
  VectorType incdisp;
  disp.Fill(0.0);
  incdisp.Fill(0.0);
  IteratorType Iterator( wm, wm->GetLargestPossibleRegion().GetSize() );
  timedone = false;
  float    totalerr = 1.e8, lasterr = 1.e10;
  unsigned its = 0;
  wmgrow->FillBuffer(0);
  float         dmag = 0;
  float         thicknesserror = 0;
  unsigned long thickerrct = 0;
  unsigned int  badct = 0;
  float         thickoffset = 0;
  bool          checknans = false;

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
    ct = 0;
    incrfield->FillBuffer(zero);
    incrfield->FillBuffer(zero);
    incrinvfield->FillBuffer(zero);

    // generate phi
    corrfield->FillBuffer(zero);
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
    float origthickprior = thickprior;

    while( ttiter < numtimepoints )    // N time integration points
      {
      //	  m_MFR->Compose(incrinvfield,invfield,NULL);
      m_MFR->ComposeDiffs(invfield, incrinvfield, invfield, 1);

      if( debug )
        {
        std::cout << " exp " << std::endl;
        }
      // Integrate the negative velocity field to generate diffeomorphism corrfield step 3(a)
      //	  ExpDiffMap<ImageType,DeformationFieldType>( velofield, corrfield, wm, -1, numtimepoints-ttiter);
      corrfield = ExpDiffMap<ImageType, DeformationFieldType>( velofield,  wm, -1, numtimepoints - ttiter);
      // why integrate velofield, but only compose incrinvfield ???
      // technically, we should warp the gm image by corrfield but this can be avoided
      if( debug )
        {
        std::cout << " gmdef " << std::endl;
        }
      typename ImageType::Pointer gmdef = m_MFR->WarpImageBackward(gm, corrfield);
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

      typedef DeformationFieldType GradientImageType;
      typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
      GradientImageFilterType;
      typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
      GradientImageFilterPointer gfilter = GradientImageFilterType::New();
      gfilter->SetInput(  surfdef );
      gfilter->SetSigma( smoothingsigma );
      gfilter->Update();
      typename DeformationFieldType::Pointer   lapgrad2 = gfilter->GetOutput();

      GradientImageFilterPointer gfilter2 = GradientImageFilterType::New();
      gfilter2->SetInput(  gm );
      gfilter2->SetSigma( smoothingsigma );
      gfilter2->Update();
      typename DeformationFieldType::Pointer   lapgrad3 = gfilter2->GetOutput();

      typename ImageType::Pointer lapjac = ComputeJacobian<ImageType, DeformationFieldType>(invfield);
      IteratorType xxIterator( lapjac, lapjac->GetLargestPossibleRegion().GetSize() );
      xxIterator.GoToBegin();
      float maxlapgrad2mag = 0;
      while(  !xxIterator.IsAtEnd()  )
        {
        typename ImageType::IndexType speedindex = xxIterator.GetIndex();
        if( gm->GetPixel(speedindex) >= 0.5 )
          {
//	      float thkval=thkdef->GetPixel(speedindex);
          float thkval = finalthickimage->GetPixel(speedindex);
          float prior = 1;
          if( spatprior )
            {
            prior = wpriorim->GetPixel(speedindex);
            //		float partialvol=surfdef->GetPixel(speedindex) ;
            // if (prval > 0.5 && partialvol >1.e-3 ) prior = prval/partialvol;//7;//0.5*origthickprior;// prval;
            // if (prior > 100 ) prior=100;  /** Potential cause of problem 1 -- this line added */
            }
          // else thickprior = origthickprior;
          // } else
          thickprior = origthickprior;

          VectorType wgradval = lapgrad2->GetPixel(speedindex); // velofield->GetPixel(speedindex);//
          //	      VectorType ggradval=lapgrad3->GetPixel(speedindex);
          double dp = 0;
          double gmag = 0, wmag = 0;
          for( unsigned kq = 0; kq < ImageDimension; kq++ )
            {
            //	      gmag+= ggradval[kq]*ggradval[kq];
            wmag += wgradval[kq] * wgradval[kq];
            }
          if( fabs(wmag) < 1.e-6 )
            {
            wmag = 0;
            }
          if( fabs(gmag) < 1.e-6 )
            {
            gmag = 0;
            }
          gmag = sqrt(gmag);
          wmag = sqrt(wmag);
          if( checknans )
            {
            if( vnl_math_isnan(wmag) || vnl_math_isinf(wmag) )
              {
              wgradval.Fill(0);
              lapgrad2->SetPixel(speedindex, wgradval);
              wmag = 0;
              }
            }
          //      if (gmag > 0 && wmag > 0)
          //	{
          //		for (unsigned kq=0;kq<ImageDimension; kq++) dp+= ggradval[kq]/gmag*wgradval[kq]/wmag;
          //	}
//	      if (fabs(dp) < 0.6) dp=0;
//	      dp=fabs(dp);
          dp = 1.0; // -dp;
          // tempim->SetPixel(speedindex,dp);

          double fval = (thickprior - thkval);
          double sigmoidf = 1;
//	      double sigmoidf=1.0- thkval/thickprior;
          if( fval >= 0 )
            {
            sigmoidf = 1.0 / (1.0 + exp( (-1.0) * fval * 0.01) );
            }
          if( fval < 0 )
            {
            sigmoidf = -1.0 * (1.0 - thickprior / thkval);
            }
//	      else sigmoidf=1.0- thkval/thickprior;
          float thkscale = thickprior / thkval;
          if( thkscale < 0.99 )
            {
            thkscale = 0.99;               // *+thkscale*0.1;
            }
          if( fval < 0 )
            {
            velofield->SetPixel(speedindex, velofield->GetPixel(speedindex) * thkscale);
            }

          float dd = surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex);
          //	      float gmd=gmdef->GetPixel(speedindex);
          totalerr += fabs(dd);
          if( wm->GetPixel(speedindex) > 0.5 && bsurf->GetPixel(speedindex) < 0.5 )
            {
            dd = 0;
            }
          float stopval = gm->GetPixel(speedindex);
//	      if (wm->GetPixel(speedindex) > stopval) stopval=1;
          float jwt = xxIterator.Get();
          if( jwt < 1 )
            {
            jwt = 1;
            }
          else
            {
            jwt = 1.0 / xxIterator.Get();
            }
//	      dd*=stopval*jwt*thindef->GetPixel(speedindex)*sigmoidf*gradstep*dp*gmd*jwt;
          dd *= stopval * sigmoidf * gradstep * jwt * prior; // speed function here IMPORTANT!!
          if( checknans )
            {
            if( vnl_math_isnan(dd) || vnl_math_isinf(dd) )
              {
              dd = 0;
              }
            }
          if( wmag * dd > 1 )
            {
            dd = stopval * (surfdef->GetPixel(speedindex) - gmdef->GetPixel(speedindex) ) * gradstep;
            }
          lapjac->SetPixel(speedindex, dd);
          //	              std::cout <<" dd " << dd << " prior " << prior << " wmag " << wmag << std::endl;
          if( wmag * dd > maxlapgrad2mag )
            {
            maxlapgrad2mag = wmag * dd;
            }
          }
        else
          {
          lapjac->SetPixel(speedindex, 0);
          }
        ++xxIterator;
        }

      if( maxlapgrad2mag < 1.e-4 )
        {
        maxlapgrad2mag = 1.e9;
        }
      lapjac = SmoothImage<ImageType>(lapjac, 1);
      //	  if (ImageDimension==2) WriteJpg<ImageType>(surfdef,"surfdef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(thindef,"thindef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(gmdef,"gmdef.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(lapjac,"diff.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(wpriorim,"prior.jpg");
      // if (ImageDimension==2) WriteJpg<ImageType>(thkdef,"thick2.jpg");
//	  if (ImageDimension==2) WriteJpg<ImageType>(tempim,"dotp.jpg");
      // exit(0);

      /* Now that we have the gradient image, we need to visit each voxel and compute objective function */
      std::cout << " maxlapgrad2mag " << maxlapgrad2mag << std::endl;
      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        velind = Iterator.GetIndex();
        //	      float currentthickvalue=finalthickimage->GetPixel(velind);
        VectorType wgradval = lapgrad2->GetPixel(velind); // *5.0/(maxlapgrad2mag*(float)numtimepoints);

        disp = wgradval * lapjac->GetPixel(velind);
        incrfield->SetPixel(velind, incrfield->GetPixel(velind) + disp);

        if( ttiter == 0 ) // make euclidean distance image
          {
          dmag = 0;
          disp = corrfield->GetPixel(velind);
          for( unsigned int jj = 0; jj < ImageDimension; jj++ )
            {
            dmag += disp[jj] * disp[jj];
            }
          float bval = bsurf->GetPixel(velind);
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
          dmag = sqrt(dmag) * bval;
          thickimage->SetPixel(velind, dmag);
          totalimage->SetPixel(velind, dmag);
          hitimage->SetPixel(velind, bval);
          }
        else if( gm->GetPixel(velind) >= 0.5 )
          {
          float thkval = thkdef->GetPixel(velind);
          float putval = thindef->GetPixel(velind);
          //	  float getval=hitimage->GetPixel(velind);
          hitimage->SetPixel(velind, hitimage->GetPixel(velind) + putval);
          totalimage->SetPixel(velind, totalimage->GetPixel(velind) + thkval);
          }

        //	      std::cout << "disp " << incrfield->GetPixel(velind) << " hit " << hitimage->GetPixel(velind) << " thk
        // " << totalimage->GetPixel(velind) << std::endl;
        ++Iterator;
        }

      //	  if (ttiter ==0) {
      // WriteImage<ImageType>(totalimage,"Ztotalimage.nii.gz");
      // WriteImage<ImageType>(hitimage,"Zhitimage.nii.gz");
      // WriteImage<ImageType>(lapjac,"Zlapjac.nii.gz"); }

      Iterator.GoToBegin();
      while(  !Iterator.IsAtEnd()  )
        {
        incrinvfield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() ) );
        ++Iterator;
        }

      ttiter++;
      }

    Iterator.GoToBegin();
    float maxth = 0;
    while(  !Iterator.IsAtEnd()  )
      {
      velind = Iterator.GetIndex();
      // increment velocity field at every voxel v = v + u, step 4
      velofield->SetPixel(Iterator.GetIndex(), velofield->GetPixel(Iterator.GetIndex() )
                          + incrfield->GetPixel(Iterator.GetIndex() ) );
      float hitval = hitimage->GetPixel(velind);
      float thkval = 0;
      if( hitval >= 0.001 )  /** potential source of problem 2 -- this value could be smaller ... */
        {
        thkval = totalimage->GetPixel(velind) / hitval - thickoffset;
        }
      if( thkval < 0 )
        {
        thkval = 0;
        }
      finalthickimage->SetPixel(velind, thkval);
      if( thkval > maxth )
        {
        maxth = thkval;
        }
      ++Iterator;
      }

    if( debug )
      {
      std::cout << " now smooth " << std::endl;
      }
    m_MFR->SmoothDeformationFieldGauss(velofield, smoothingsigma);
    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DeformationFieldType>(velofield,velofieldname.c_str());
    // std::string incrfieldname = outname + "incrfield";
    // WriteDisplacementField<DeformationFieldType>(incrfield,incrfieldname.c_str());

    // std::string tname = outname + "dork1.nii.gz";
    // WriteImage<ImageType>(hitimage,tname.c_str());
    // tname = outname + "dork2.nii.gz";
    // WriteImage<ImageType>(totalimage,tname.c_str());
    if( thickerrct == 0 )
      {
      thickerrct = 1;
      }
    std::cout << " error " << totalerr << " at it " << its  << " th-err " << thicknesserror / (float)thickerrct
              << " max thick " << maxth << std::endl;
//    std::string sulcthickname =outname + "sulcthick.nii";
    finalthickimage->SetDirection(omat);
    WriteImage<ImageType>(finalthickimage, outname.c_str() );
    finalthickimage->SetDirection(fmat);
    //    if (ImageDimension==2) WriteJpg<ImageType>(finalthickimage,"thick.jpg");
    //    std::string velofieldname = outname + "velofield";
    // WriteDisplacementField<DeformationFieldType>(velofield,velofieldname.c_str());
    if( debug )
      {
      std::cout << "outside it " << its << std::endl;
      }
    // std::cin.get();
    }

  thickimage->FillBuffer(0);
  typename ImageType::Pointer thkdef = m_MFR->WarpImageBackward(finalthickimage, invfield);
  Iterator.GoToBegin();
  while(  !Iterator.IsAtEnd()  )
    {
    float tt1 = finalthickimage->GetPixel(Iterator.GetIndex() );
    float tt = thkdef->GetPixel(Iterator.GetIndex() );
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
  if( argc < 5 )
    {
    std::cout << "Useage ex:   " << argv[0]
              <<
    " ImageDimension WM.nii GM.nii   Out.nii {GradStep-1-2D,2-3D}   {#Its-~50}  {ThickPriorValue-6} {Bool-use-curvature-prior} {smoothing} {BoolUseEuclidean?}"
              << std::endl;
    std::cout << " this is a kind of binary image registration thing with diffeomorphisms " << std::endl;
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
