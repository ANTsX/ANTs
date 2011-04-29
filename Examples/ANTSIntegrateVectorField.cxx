
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"

#include "itkWarpImageFilter.h"

#include "itkImageFileWriter.h"

// #include "itkScalarImageToHistogramGenerator.h"
// #include "itkImageToHistogramGenerator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "ReadWriteImage.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

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

template <class TImage, class TField, class TInterp, class TInterp2>
float IntegrateLength( typename TImage::Pointer gmsurf,  typename TImage::Pointer thickimage,
                       typename TImage::IndexType velind,  typename TField::Pointer lapgrad,  float itime,
                       float starttime, float finishtime, bool timedone, float deltaTime,
                       typename TInterp::Pointer vinterp, typename TImage::SpacingType spacing, float vecsign,
                       float timesign, float gradsign )
{
  typedef   TField                                                 TimeVaryingVelocityFieldType;
  typedef typename TField::PixelType                               VectorType;
  typedef itk::ImageRegionIteratorWithIndex<TField>                FieldIterator;
  typedef typename TField::IndexType                               DIndexType;
  typedef typename TField::PointType                               DPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TField, float> DefaultInterpolatorType;

  VectorType zero;
  zero.Fill(0);
  VectorType disp;
  disp.Fill(0);
  unsigned int ct = 0;
  DPointType   pointIn1;
  DPointType   pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;
  enum { ImageDimension = TImage::ImageDimension };
  typedef typename TImage::IndexType IndexType;
  unsigned int m_NumberOfTimePoints = 2;
  IndexType    index;
  for( unsigned int jj = 0; jj < ImageDimension; jj++ )
    {
    index[jj] = velind[jj];
    pointIn1[jj] = velind[jj] * lapgrad->GetSpacing()[jj];
    }
  itime = starttime;
  timedone = false;
  float totalmag = 0;
  while( !timedone )
    {
    float scale = 1;  // *m_DT[timeind]/m_DS[timeind];
    //     std::cout << " scale " << scale << std::endl;
    double itimetn1 = itime - timesign * deltaTime * scale;
    double itimetn1h = itime - timesign * deltaTime * 0.5 * scale;
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
    IndexType index;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
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
      pointIn2[jj] = disp[jj] + pointIn1[jj];
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

    f1 = vinterp->EvaluateAtContinuousIndex( Y1 );
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      Y2[jj] += f1[jj] * deltaTime * 0.5;
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
      Y3[jj] += f2[jj] * deltaTime * 0.5;
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
      Y4[jj] += f3[jj] * deltaTime;
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
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      pointIn3[jj] = pointIn2[jj] + gradsign * vecsign * deltaTime / 6.0
        * ( f1[jj] + 2.0 * f2[jj] + 2.0 * f3[jj] + f4[jj] );
      }

    VectorType out;
    float      mag = 0, dmag = 0, voxmag = 0;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      out[jj] = pointIn3[jj] - pointIn1[jj];
      mag += (pointIn3[jj] - pointIn2[jj]) * (pointIn3[jj] - pointIn2[jj]);
      voxmag += (pointIn3[jj] - pointIn2[jj]) / spacing[jj] * (pointIn3[jj] - pointIn2[jj]) / spacing[jj];
      dmag += (pointIn3[jj] - pointIn1[jj]) * (pointIn3[jj] - pointIn1[jj]);
      disp[jj] = out[jj];
      }
    voxmag = sqrt(voxmag);
    dmag = sqrt(dmag);
    totalmag += sqrt(mag);

    ct++;
    //      if (!propagate) //thislength=dmag;//
//         thislength += totalmag;
    itime = itime + deltaTime * timesign;
    IndexType myind;
    for( unsigned int qq = 0; qq <  ImageDimension; qq++ )
      {
      myind[qq] = (unsigned long)(pointIn3[qq] / spacing[qq] + 0.5);
      }

    if( gmsurf->GetPixel(myind) < 1 )
      {
      timedone = true;
      }
    if( ct >  1000 )
      {
      std::cout << " stopping b/c exceed 1000 points " << voxmag <<  std::endl;  timedone = true;
      }
    if( voxmag < 0.1 )
      {
      timedone = true;
      }
    }

  return totalmag;
}

template <unsigned int ImageDimension>
int IntegrateVectorField(int argc, char *argv[])
{
  typedef float                                  PixelType;
  typedef itk::Vector<float, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::Image<PixelType, ImageDimension>  ImageType;
  typedef itk::ImageFileReader<ImageType>        readertype;
  typedef itk::ImageFileWriter<ImageType>        writertype;
  typedef typename  ImageType::IndexType         IndexType;
  typedef typename  ImageType::SizeType          SizeType;
  typedef typename  ImageType::SpacingType       SpacingType;

  double      dT = 0.001;
  float       gradstep = 1. / dT; // atof(argv[3])*(-1.0);
  std::string vectorfn = std::string(argv[1]);
  std::string roifn = std::string(argv[2]);
  int         argct = 3;
  std::string outname = std::string(argv[argct]); argct++;
  std::string lenoutname = std::string("");
  if( argc > argct )
    {
    lenoutname = std::string(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    gradstep *= atof(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer ROIimage;
  ReadImage<ImageType>(ROIimage, roifn.c_str() );
  typename ImageType::Pointer thickimage;
  ReadImage<ImageType>(thickimage, roifn.c_str() );
  thickimage->FillBuffer(0);
  typename DeformationFieldType::Pointer VECimage;
  ReadImage<DeformationFieldType>(VECimage, vectorfn.c_str() );
  SpacingType spacing = ROIimage->GetSpacing();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType Iterator( ROIimage, ROIimage->GetLargestPossibleRegion().GetSize() );

  double timezero = 0; // 1
  typename ImageType::SizeType s = ROIimage->GetLargestPossibleRegion().GetSize();
  double timeone = 1;          // (s[ImageDimension]-1-timezero);
  float  starttime = timezero; // timezero;
  float  finishtime = timeone; // s[ImageDimension]-1;//timeone;

  typename DeformationFieldType::IndexType velind;
  float timesign = 1.0;
  if( starttime  >  finishtime )
    {
    timesign = -1.0;
    }
  typedef   DeformationFieldType                                                         TimeVaryingVelocityFieldType;
  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType>                        FieldIterator;
  typedef typename DeformationFieldType::IndexType                                       DIndexType;
  typedef typename DeformationFieldType::PointType                                       DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType                               VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType                               VPointType;
  typedef itk::VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType, float> DefaultInterpolatorType;
  typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, float>         DefaultInterpolatorType2;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
  VectorType zero;
  zero.Fill(0);

  DPointType pointIn1;
  DPointType pointIn2;
  typename DefaultInterpolatorType::ContinuousIndexType  vcontind;
  DPointType pointIn3;

  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> VIteratorType;
  VIteratorType VIterator( VECimage, VECimage->GetLargestPossibleRegion().GetSize() );
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
    ++VIterator;
    }

  Iterator.GoToBegin();
  while(  !Iterator.IsAtEnd()  )
    {
    velind = Iterator.GetIndex();
    float      itime = starttime;
    bool       timedone = false;
    VectorType disp;
    disp.Fill(0.0);
    double deltaTime = dT, vecsign = 1.0;
    float  gradsign = 1.0;
    if( ROIimage->GetPixel(velind) == 2 )
      {
      vinterp->SetInputImage(VECimage);
      gradsign = -1.0; vecsign = -1.0;
      float len1 = IntegrateLength<ImageType, DeformationFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
          (ROIimage, thickimage, velind, VECimage,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
          spacing, vecsign, gradsign, timesign);

      gradsign = 1.0;  vecsign = 1;
      float len2 = IntegrateLength<ImageType, DeformationFieldType, DefaultInterpolatorType, ScalarInterpolatorType>
          (ROIimage, thickimage, velind, VECimage,  itime, starttime, finishtime,  timedone,  deltaTime,  vinterp,
          spacing, vecsign, gradsign, timesign );

      float totalength = len1 + len2;
      thickimage->SetPixel(velind, totalength);
      if( (totalength) > 0 )
        {
        std::cout << " len1 " << len1 << " len2 " << len2 << " ind " << velind << std::endl;
        }
      }
    ++Iterator;
    }

  WriteImage<ImageType>(thickimage, lenoutname.c_str() );

  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << "Usage:   " << argv[0]
              << "  VecImageIN.nii.gz ROIMaskIN.nii.gz FibersOUT.vtk  LengthImageOUT.nii.gz   " << std::endl;
    std::cout
      <<
    " The vector field should have vectors as voxels , the ROI is an integer image, fibers out will be vtk text files .... "
      << std::endl;
    std::cout << "  ROI-Mask controls where the integration is performed and the start point region ... " << std::endl;
    std::cout << " e.g. the brain will have value 1 , the ROI has value 2 , then all starting seed points "
              << std::endl;
    std::cout
      << " for the integration will start in the region labeled 2 and be constrained to the region labeled 1. "
      << std::endl;
    return 1;
    }

  std::string               ifn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(ifn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(ifn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim =  imageIO->GetNumberOfDimensions();

  switch( dim )
    {
    case 2:
      {
      IntegrateVectorField<2>(argc, argv);
      }
      break;
    case 3:
      {
      IntegrateVectorField<3>(argc, argv);
      }
      break;
    case 4:
      {
      IntegrateVectorField<4>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return EXIT_SUCCESS;

  return 1;
}
