



// #include "curvatureapp.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <vcl_fstream.h>
#include "itkSurfaceCurvatureBase.h"
#include "itkSurfaceImageCurvature.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"

#include "itkCenteredAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

namespace itk
{
typedef itk::Image<unsigned char, 3> ImageType;
typedef itk::Image<float, 3>         floatImageType;

// Template function to fill in an image with a circle.
template <class TImage>
typename TImage::Pointer
GenerateEllipse(
  double longradius,
  double shortradius, float rot = 0, int axis = 0)
{
  typedef typename TImage::PixelType PixelType;
  enum { ImageDimension = TImage::ImageDimension };
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef typename ImageType::IndexType         IndexType;
  typedef typename ImageType::SizeType          SizeType;
  typedef typename ImageType::RegionType        RegionType;

  // --------------------------------------------------------
  std::cout << "Generate image ";
  std::cout << std::endl;

  unsigned long sizeArray[ImageDimension];
  double        center[ImageDimension];
  for( int i = 0; i < ImageDimension; i++ )
    {
    float edge = 20.0;
    sizeArray[i] = (unsigned long) (2.0 * longradius + edge * 2.0);
    center[i] = (double)(edge + longradius);
    }
  SizeType size;
  size.SetSize( sizeArray );

  IndexType index;
  index.Fill( 0 );

  RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  typename ImageType::Pointer image = ImageType::New();
  image->SetLargestPossibleRegion( region );
  image->SetBufferedRegion( region );
  image->Allocate();

  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.Begin();

  double radius1;
  double radius2;
  double radius3;

  if( axis == 0 )
    {
    radius2 = shortradius; radius1 = longradius;  radius3 = shortradius;
    }
  if( axis == 1 )
    {
    radius2 = longradius;  radius1 = shortradius; radius3 = shortradius;
    }
  if( axis == 2 )
    {
    radius2 = shortradius; radius1 = shortradius; radius3 = longradius;
    }
  for( ; !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    double distance1 = 0;
    double distance2 = 0;
    double distance3 = 0;
    double distance  = 0;
    for( unsigned int j = 0; j < TImage::ImageDimension; j++ )
      {
      if( j == 0 )
        {
        distance1 += vnl_math_sqr( ( (double) index[j] - center[j]) / radius1);
        }
      if( j == 1 )
        {
        distance2 += vnl_math_sqr( ( (double) index[j] - center[j]) / radius2);
        }
      if( j == ImageDimension - 1 )
        {
        distance3 += vnl_math_sqr( ( (double) index[j] - center[j]) / radius3);
        }
      }
    distance = distance1 + distance2 + distance3;
    if( distance <= 1.0 )
      {
      it.Set( (PixelType) ( (1. - distance + 1.) * 125.0) );
      }
    else
      {
      it.Set( 0 );
      }
    }

  // Create an affine transformation
  typedef itk::CenteredAffineTransform<double, ImageDimension> CenteredAffineTransformType;
  typename CenteredAffineTransformType::Pointer aff = CenteredAffineTransformType::New();
  aff->SetCenter( center );
  aff->Rotate(0, 1, rot, 0);

  // Create a linear interpolation image function

  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetInputImage(image);

  // Create and configure a resampling filter
  typename ResampleImageFilter<ImageType, ImageType>::Pointer resample;
  resample = itk::ResampleImageFilter<ImageType, ImageType>::New();
  resample->SetInput(image);
  resample->SetSize(size);
  resample->SetTransform(aff);
  resample->SetInterpolator(interp);

  index.Fill( 0 );
  resample->SetOutputStartIndex( index );

  typename ImageType::PointType origin;
  origin.Fill( 0.0 );
  resample->SetOutputOrigin( origin );

  typename ImageType::SpacingType spacing;
  spacing.Fill( 1.0 );
  resample->SetOutputSpacing( spacing );

  // Run the resampling filter
  resample->Update();

  // Check if desired results were obtained
  bool passed = true;
  typename ImageType::RegionType region2;
  region2 = resample->GetOutput()->GetRequestedRegion();
  itk::ImageRegionIteratorWithIndex<ImageType>
  iter2(resample->GetOutput(), region2);

  return resample->GetOutput();
//      return image;
}
}

namespace itk
{
template <class TImage>
typename TImage::Pointer BinaryThreshold(
  typename TImage::PixelType low,
  typename TImage::PixelType high,
  typename TImage::PixelType replaceval, typename TImage::Pointer input)
{
  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef BinaryThresholdImageFilter<TImage, TImage>
    InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  inputThresholder->SetInsideValue( 1 );
  inputThresholder->SetOutsideValue( 0 );

  inputThresholder->SetLowerThreshold( (PixelType) low );
  inputThresholder->SetUpperThreshold( (PixelType) high);
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <class TImage>
typename TImage::Pointer  MorphologicalClosing
  (float rad, typename TImage::Pointer input)
{
  typedef TImage ImageType;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::PixelType PixelType;

  std::cout << " Opening the image " << std::endl;
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

  binaryDilate->SetInput( input );
  binaryDilate->Update();
  binaryErode->SetInput( binaryDilate->GetOutput() );
  binaryErode->Update();
  return binaryErode->GetOutput();
}

template <class ImageType>
typename ImageType::Pointer SegmentImage(typename ImageType::Pointer image, float thresh, float ballradius)
{
  std::string fn;

  enum { ImageDimension = ImageType::ImageDimension };

  typename ImageType::Pointer threshimg;
  threshimg =  BinaryThreshold<ImageType>(thresh, 1.e5, 1, image);

  // typename TImage::Pointer  MorphologicalClosing
  // (typename TImage::PixelType foreground,float rad, typename TImage::Pointer input)

    {
    typename itk::ImageFileWriter<ImageType>::Pointer writer;
    writer = itk::ImageFileWriter<ImageType>::New();
    std::string ofn = "thresh1.nii";
    writer->SetFileName(ofn.c_str() );
    writer->SetInput( threshimg );
    writer->Write();
    }

  threshimg = MorphologicalClosing<ImageType>(ballradius, threshimg);
    {
    typename itk::ImageFileWriter<ImageType>::Pointer writer;
    writer = itk::ImageFileWriter<ImageType>::New();
    std::string ofn = "thresh2.nii";
    writer->SetFileName(ofn.c_str() );
    writer->SetInput( threshimg );
    writer->Write();
    }
  threshimg = MorphologicalClosing<ImageType>(ballradius, threshimg);

  typedef itk::ImageRegionIteratorWithIndex<ImageType> iteratorType;
  iteratorType Iter( threshimg, threshimg->GetLargestPossibleRegion() );
  for( Iter.GoToBegin(); !Iter.IsAtEnd(); ++Iter  )
    {
    if( Iter.Get() > 0 )
      {
      Iter.Set(image->GetPixel(Iter.GetIndex() ) );
      }
    else
      {
      Iter.Set(-1.0);
      }
    }

    {
    typename itk::ImageFileWriter<ImageType>::Pointer writer;
    writer = itk::ImageFileWriter<ImageType>::New();
    std::string ofn = "thresh3.nii";
    writer->SetFileName(ofn.c_str() );
    writer->SetInput( threshimg );
    writer->Write();
    }

  return threshimg;
}

void test1()
{
  typedef itk::SurfaceCurvatureBase<ImageType> ParamType;
  ParamType::Pointer Parameterizer = ParamType::New();

//  Parameterizer->TestEstimateTangentPlane(p);
  Parameterizer->FindNeighborhood();
//  Parameterizer->WeightedEstimateTangentPlane(  Parameterizer->GetOrigin() );

  Parameterizer->EstimateTangentPlane(
    Parameterizer->GetAveragePoint() );
  Parameterizer->PrintFrame();
//  Parameterizer->SetOrigin(Parameterizer->GetAveragePoint());
  for( int i = 0; i < 3; i++ )
    {
    Parameterizer->ComputeWeightsAndDirectionalKappaAndAngles
      (Parameterizer->GetOrigin() );
    Parameterizer->ComputeFrame(Parameterizer->GetOrigin() );
    Parameterizer->EstimateCurvature();
    Parameterizer->PrintFrame();
    }

  Parameterizer->ComputeJoshiFrame(Parameterizer->GetOrigin() );  Parameterizer->PrintFrame();
  std::cout << " err 1 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin() )
            << " err 2 " << Parameterizer->ErrorEstimate(Parameterizer->GetOrigin(), -1) << std::endl;
}
}

using namespace itk;
int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << " usage :  SurfaceCurvature FileNameIn FileNameOut sigma option  " << std::endl;
    std::cout << " e.g  :   SurfaceCurvature    BrainIn.nii BrainOut.nii   3  0 " << std::endl;
    std::cout << " option 0 means just compute mean curvature from intensity " << std::endl;
    std::cout << " option 5 means characterize surface from intensity " << std::endl;
    return 0;
    }

//  const unsigned int ImageDimension=ImageType::ImageDimension;
  typedef itk::Image<float, 3> ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef itk::Image<float, 3>                  floatImageType;
  typedef itk::SurfaceImageCurvature<ImageType> ParamType;
  ParamType::Pointer Parameterizer = ParamType::New();

  typedef  itk::ImageFileReader<ImageType> ReaderType;
//  typedef  itk::RawImageIO< unsigned char ,ImageDimension>   RawIOType;
  typedef  ImageType::PixelType PixType;

  int   opt = 0;
  float sig = 1.0;
  if( argc > 3 )
    {
    sig = atof( argv[3]);
    }
  std::cout << " sigma " << sig << std::endl;
  if( argc > 4 )
    {
    opt = (int) atoi(argv[4]);
    }

  if( opt < 0 )
    {
    // if (argc <= 7)  std::cout << " requires addt'l cmd line inputs : longradius  shortradius rot axis " << std::endl;
    // floatImageType::Pointer ell =
    // GenerateEllipse<floatImageType>(atof(argv[4]),atof(argv[5]),atof(argv[6]),atoi(argv[7]));
    // itk::ImageFileWriter<floatImageType>::Pointer writer;
    // writer = itk::ImageFileWriter<floatImageType>::New();
    // std::string ofn="ellipse"+std::string(argv[1]);
    // std::cout << " writing result " << ofn <<  std::endl;
    // writer->SetFileName(ofn.c_str());
    // writer->SetInput( ell );
    // writer->Write();
    return 0;
    }

  ReaderType::Pointer reader  = ReaderType::New();
  std::string         fn = std::string(argv[1]);
  reader->SetFileName(fn.c_str() );

  std::string m_FileExtension = fn;
  m_FileExtension = "." + m_FileExtension.substr(m_FileExtension.find(".") + 1, m_FileExtension.length() );

  reader->Update();
  ImageType::Pointer input = reader->GetOutput();
  std::cout << " done reading ";

  //  float ballradius = 2.0;
  // if (argc >= 6) ballradius = (float) atof(argv[5]);
  // if (ballradius > 0 && thresh > 0) input = SegmentImage<ImageType>(input, thresh, ballradius);

  Parameterizer->SetInput(input);

  //  Parameterizer->ProcessLabelImage();
  Parameterizer->SetNeighborhoodRadius( 1. );
//  std::cout << " set sig " ;  std::cin >> sig;
  if( sig <= 0.5 )
    {
    sig = 1.66;
    }
  std::cout << " sigma " << sig << " option " << opt << std::endl;
  Parameterizer->SetSigma(sig);

  if( opt == 1 )
    {
    Parameterizer->SetUseLabel(true);
    Parameterizer->SetUseGeodesicNeighborhood(true);
    }
  else
    {
    Parameterizer->SetUseLabel(false);
    Parameterizer->SetUseGeodesicNeighborhood(false);
    float sign = 1.0;
    if( opt == 3 )
      {
      sign = -1.0;
      }
    std::cout << " setting outward direction as " << sign;
    Parameterizer->SetkSign(sign);
    Parameterizer->SetThreshold(0);
    }
//  Parameterizer->ComputeSurfaceArea();
//  Parameterizer->IntegrateFunctionOverSurface();
//  Parameterizer->IntegrateFunctionOverSurface(true);

  std::cout << " computing frame " << std::endl;
  if( opt != 5 )
    {
    Parameterizer->ComputeFrameOverDomain( 3 );
    }
  else
    {
    Parameterizer->ComputeFrameOverDomain( opt );
    }

  //   Parameterizer->SetNeighborhoodRadius( 2 );
  //  Parameterizer->LevelSetMeanCurvature();
  //  Parameterizer->SetNeighborhoodRadius( 2.9   );
  //  Parameterizer->IntegrateFunctionOverSurface(false);
  //  Parameterizer->SetNeighborhoodRadius( 1.5  );
  //  Parameterizer->IntegrateFunctionOverSurface(true);
  //   for (int i=0; i<1; i++) Parameterizer->PostProcessGeometry();

  ImageType::Pointer      output = NULL;
  floatImageType::Pointer smooth = NULL;

  //  Parameterizer->GetFunctionImage()->SetSpacing( input->GetSpacing() );
  //  Parameterizer->GetFunctionImage()->SetDirection( input->GetDirection() );
  //  Parameterizer->GetFunctionImage()->SetOrigin( input->GetOrigin() );
  //  smooth->SetSpacing(reader->GetOutput()->GetSpacing());
  // SmoothImage(Parameterizer->GetFunctionImage(),smooth,3);
  // NormalizeImage(smooth,output,mn);
  //  NormalizeImage(Parameterizer->GetFunctionImage(),output,mn);

  typedef unsigned char ptype;
  itk::ImageFileWriter<floatImageType>::Pointer writer;
  writer = itk::ImageFileWriter<floatImageType>::New();
  writer->SetFileName(argv[2]);
  writer->SetInput( Parameterizer->GetFunctionImage() );
  writer->Write();

  std::cout << " done writing ";

  return 1;
}
