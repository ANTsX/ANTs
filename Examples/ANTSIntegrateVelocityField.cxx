#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkANTSImageRegistrationOptimizer.h"
#include "itkTimeVaryingVelocityFieldIntegrationImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkTimeVaryingVelocityFieldTransform.h"

#include "itkImageFileWriter.h"

// #include "itkScalarImageToHistogramGenerator.h"
// #include "itkImageToHistogramGenerator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vnl/algo/vnl_determinant.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "ReadWriteImage.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

template <unsigned int ImageDimension>
int IntegrateVelocityField(int argc, char *argv[])
{
  int         argct = 1;
  std::string imgfn = std::string(argv[argct]); argct++;

  std::string vectorfn = std::string(argv[argct]); argct++;
  std::string outname = std::string(argv[argct]); argct++;

  typedef float PixelType;
  PixelType timezero = 0;
  PixelType timeone = 1;
  PixelType dT = 0.01;
  if( argc > argct )
    {
    timezero = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    timeone = atof(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    dT = atof(argv[argct]);
    }
  argct++;
  std::cout << " time-0 " << timezero << " dt " << dT << " time-1 " << timeone << std::endl;
  PixelType starttime = timezero;
  PixelType finishtime = timeone;
  typedef float                                       PixelType;
  typedef itk::Vector<PixelType, ImageDimension>      VectorType;
  typedef itk::Image<VectorType, ImageDimension>      DisplacementFieldType;
  typedef itk::Image<VectorType, ImageDimension + 1>  TimeVaryingVelocityFieldType;
  typedef itk::Image<PixelType, ImageDimension>       ImageType;
  typedef typename  ImageType::IndexType              IndexType;
  typedef typename  ImageType::SizeType               SizeType;
  typedef typename  ImageType::SpacingType            SpacingType;
  typedef TimeVaryingVelocityFieldType                tvt;
  typedef itk::ImageFileReader<tvt>                   readertype;
  typedef itk::ImageFileWriter<DisplacementFieldType> writertype;
  typename ImageType::Pointer image;
  ReadImage<ImageType>(image, imgfn.c_str() );
  typename tvt::Pointer timeVaryingVelocity;
  ReadImage<tvt>(timeVaryingVelocity, vectorfn.c_str() );

  typename DisplacementFieldType::Pointer deformation = DisplacementFieldType::New();
  deformation->SetSpacing(    image->GetSpacing() );
  deformation->SetOrigin(     image->GetOrigin() );
  deformation->SetDirection(  image->GetDirection() );
  deformation->SetRegions(    image->GetLargestPossibleRegion() );
  deformation->Allocate();
  VectorType zero;
  zero.Fill(0);
  deformation->FillBuffer(zero);
  if( !timeVaryingVelocity )
    {
    std::cout << " No TV Field " << std::endl;  return EXIT_FAILURE;
    }
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex<tvt>                   TVFieldIterator;
  typedef typename DisplacementFieldType::IndexType                DIndexType;
  typedef typename DisplacementFieldType::PointType                DPointType;
  typedef typename TimeVaryingVelocityFieldType::IndexType         VIndexType;
  typedef typename TimeVaryingVelocityFieldType::PointType         VPointType;

  if( starttime < 0 )
    {
    starttime = 0;
    }
  if( starttime > 1 )
    {
    starttime = 1;
    }
  if( finishtime < 0 )
    {
    finishtime = 0;
    }
  if( finishtime > 1 )
    {
    finishtime = 1;
    }

  typedef itk::TimeVaryingVelocityFieldIntegrationImageFilter
    <TimeVaryingVelocityFieldType, DisplacementFieldType> IntegratorType;
  typename IntegratorType::Pointer integrator = IntegratorType::New();
  integrator->SetInput( timeVaryingVelocity );
  integrator->SetLowerTimeBound( starttime );
  integrator->SetUpperTimeBound( finishtime );
  integrator->SetNumberOfIntegrationSteps( (unsigned int ) 1 / dT );
  integrator->Update();

  WriteImage<DisplacementFieldType>( integrator->GetOutput(), outname.c_str() );
  return 0;
}

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << "Usage:   " << argv[0]
              << " reference_image  VelocityIn.mhd DeformationOut.nii.gz  time0 time1 dT  " << std::endl;
    return 1;
    }
  std::cout << " start " << std::endl;
  std::string               ifn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(ifn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(ifn.c_str() );
  imageIO->ReadImageInformation();
  unsigned int dim =  imageIO->GetNumberOfDimensions();
  std::cout << " dim " << dim << std::endl;

  switch( dim )
    {
    case 2:
      IntegrateVelocityField<2>(argc, argv);
      break;
    case 3:
      IntegrateVelocityField<3>(argc, argv);
      break;
    case 4:
      IntegrateVelocityField<4>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return EXIT_SUCCESS;
}
