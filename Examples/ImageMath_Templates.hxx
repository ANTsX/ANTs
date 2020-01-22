/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>
#include <vnl/vnl_inverse.h>
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "antsAllocImage.h"
#include "antsSCCANObject.h"
#include "itkAlternatingValueDifferenceImageFilter.h"
#include "itkAlternatingValueSimpleSubtractionImageFilter.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkArray.h"
#include "itkAverageOverDimensionImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkBlackTopHatImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCompositeValleyFunction.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConvolutionImageFilter.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkCyclicShiftImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkFastMarchingExtensionImageFilterBase.h"
#include "itkFastMarchingExtensionImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHistogram.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageGaussianModelEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkLabelPerimeterEstimationCalculator.h"
#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkListSample.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkMRIBiasFieldCorrectionFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkMultivariateLegendrePolynomial.h"
#include "itkNeighborhood.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodFirstOrderStatisticsImageFilter.h"
#include "itkNormalVariateGenerator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkSampleToHistogramFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkSize.h"
#include "itkSliceTimingCorrectionImageFilter.h"
#include "itkSphereSpatialFunction.h"
#include "itkSplitAlternatingTimeSeriesImageFilter.h"
#include "itkSTAPLEImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSumProjectionImageFilter.h"
#include "itkTDistribution.h"
#include "itkTileImageFilter.h"
#include "itkTimeProbe.h"
#include "itkTranslationTransform.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkWhiteTopHatImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkTransformFactory.h"
#include "itkSurfaceImageCurvature.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"

#include <fstream>
#include <iostream>
#include <map> // Here I'm using a map but you could choose even other containers
#include <sstream>
#include <string>

#include "ReadWriteData.h"
#include "TensorFunctions.h"
#include "antsMatrixUtilities.h"
#include "antsFastMarchingImageFilter.h"
#include "itkFastMarchingImageFilterBase.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"

namespace ants
{

static inline std::string ANTSGetFilePrefix(const char *str)
{
  const std::string      filename = str;
  const std::string::size_type pos = filename.rfind( "." );
  const std::string            filepre = std::string( filename, 0, pos );

#if 0 // HACK:  This does nothing useful
  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      // extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    }
#endif
  return filepre;
}

static inline std::string ANTSOptionName(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            name = std::string( filename, 0, pos );

  return name;
}

static inline std::string ANTSOptionValue(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            value = std::string( filename, pos + 1, filename.length() );

  return value;
}


// int is the key, string the return value
static inline std::map<unsigned int, std::string> RoiList(std::string file)
{
  unsigned int wordindex = 0;
  std::string  tempstring = "";

  std::map<unsigned int, std::string> RoiList;
  //  RoiList[0]=std::string("Background");
  char         str[2000];
  std::fstream file_op(file.c_str(), std::ios::in);

  while( file_op >> str )
    {
    tempstring = std::string(str);
    RoiList[wordindex] = tempstring;
    wordindex++;
    }

  return RoiList; // returns the maximum index
}

template <typename T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base & (*f)(std::ios_base &) )
{
  std::istringstream iss(s);

  iss >> f >> t;

  // Check to see that there is nothing left over
  if( !iss.eof() )
    {
    return false;
    }

  return true;
}

template <typename T>
std::string ants_to_string(T t)
{
  std::stringstream istream;

  istream << t;
  return istream.str();
}


template <unsigned int ImageDimension>
int FrobeniusNormOfMatrixDifference(int argc, char *argv[])
{
  if( argc < 6 )
    {
    std::cout << " FrobeniusNormOfMatrixDifference: too few options " << std::endl;
    std::cout << " ImageMath 3 out FrobeniusNormOfMatrixDifference aff1.mat aff2.mat " << std::endl;
    return 1;
    }
  int         argct = 4;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension,
                                         ImageDimension> AffineTransformType;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();
  typedef itk::TransformFileReader TransformReaderType;
  typename TransformReaderType::Pointer transformReader1 = TransformReaderType::New();
  transformReader1->SetFileName( fn1.c_str() );
  try
    {
    transformReader1->Update();
    }
  catch( itk::ExceptionObject & /* excp */ )
    {
    std::cout << "no transformation1 that can be read" << fn1 << std::endl;
    return 0;
    }
  typename TransformReaderType::Pointer transformReader2 = TransformReaderType::New();
  transformReader2->SetFileName(  fn2.c_str()  );
  try
    {
    transformReader2->Update();
    }
  catch( itk::ExceptionObject & /* excp */ )
    {
    std::cout << "no transformation2 that can be read" << fn2 << std::endl;
    return 0;
    }
  typename AffineTransformType::Pointer aff1 =
    dynamic_cast<AffineTransformType *>( (transformReader1->GetTransformList() )->front().GetPointer() );
  typename AffineTransformType::Pointer aff2 =
    dynamic_cast<AffineTransformType *>( (transformReader2->GetTransformList() )->front().GetPointer() );

  typename AffineTransformType::MatrixType::InternalMatrixType diffmat = aff2->GetMatrix().GetVnlMatrix()
    - aff1->GetMatrix().GetVnlMatrix();
  std::cout << diffmat.frobenius_norm() << std::endl;
  return 0;
}

template <unsigned int ImageDimension>
void ClosestSimplifiedHeaderMatrix(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << " need more args -- see usage   " << std::endl;
    }
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef double RealType;
  typedef itk::Matrix<RealType, ImageDimension, ImageDimension>                           MatrixType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>( image1, fn1.c_str() );
  MatrixType A_solution = image1->GetDirection();
  //  vnl_svd<RealType>    nearestorth( image1->GetDirection().GetVnlMatrix() );
  //  vnl_matrix<RealType> A_solution = nearestorth.V() * nearestorth.U().transpose();
  MatrixType mydir;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    for ( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      RealType v = A_solution(d,dd);
      long vlong = static_cast<long>( itk::Math::abs ( v ) + 0.5 );
      if ( v < 0 ) vlong = vlong * (-1);
      mydir(d,dd) = static_cast<RealType>( vlong );
      }
  image1->SetDirection( mydir );
  // also do the origin
  typename ImageType::PointType origin = image1->GetOrigin();
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    RealType rndflt = origin[d] * 100.0;
    int rndflti = static_cast<int>( rndflt + 0.5 );
    rndflt = static_cast<RealType>( rndflti ) / 100.0;
    origin[d] = rndflt;
    }
  image1->SetOrigin( origin );
  WriteImage<ImageType>( image1, outname.c_str() );
  return;
}


template <unsigned int ImageDimension>
void ReflectionMatrix(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cout << " need more args -- see usage   " << std::endl;
    }
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int axis = std::stoi(argv[argct]);  argct++;
  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>( image1, fn1.c_str() );
  // compute center of mass
  typedef typename itk::ImageMomentsCalculator<ImageType> ImageCalculatorType;
  typename ImageCalculatorType::Pointer calculator1 = ImageCalculatorType::New();
  calculator1->SetImage(  image1 );
  typename ImageCalculatorType::VectorType fixed_center;
  fixed_center.Fill(0);
  try
    {
    calculator1->Compute();
    fixed_center = calculator1->GetCenterOfGravity();
    }
  catch( ... )
    {
    // std::cout << " zero image2 error ";
    fixed_center.Fill(0);
    }
  // std::cout << " CenterOfMass " << fixed_center << std::endl;
  typename AffineTransformType::Pointer aff = AffineTransformType::New();
  aff->SetIdentity();
  typename AffineTransformType::ParametersType myoff = aff->GetFixedParameters();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    myoff[i] = fixed_center[i];
    }
  typename AffineTransformType::MatrixType mymat = aff->GetMatrix();
  if( axis < ImageDimension )
    {
    mymat[axis][axis] = ( -1.0 );
    }
  aff->SetFixedParameters( myoff );
  aff->SetMatrix( mymat );
  typedef itk::TransformFileWriter TransformWriterType;
  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( aff );
  transformWriter->SetFileName( outname.c_str() );
#if ITK_VERSION_MAJOR >= 5
  transformWriter->SetUseCompression(true);
#endif
  transformWriter->Update();
  return;
}


template <unsigned int ImageDimension>
void MakeAffineTransform(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cout << " need more args -- see usage   " << std::endl;
    }

  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  typename AffineTransformType::Pointer aff = AffineTransformType::New();
  aff->SetIdentity();
  typedef itk::TransformFileWriter TransformWriterType;
  typename TransformWriterType::Pointer transformWriter =
    TransformWriterType::New();
  transformWriter->SetInput( aff );
  transformWriter->SetFileName( outname.c_str() );
#if ITK_VERSION_MAJOR >= 5
  transformWriter->SetUseCompression(true);
#endif
  transformWriter->Update();
  return;
}


template <unsigned int ImageDimension>
int GetLargestComponent(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string   fn1 = std::string(argv[argct]);   argct++;
  unsigned long smallest = 50;
  if( argc > argct )
    {
    smallest = std::stoi(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  // compute the voxel volume
  typename ImageType::SpacingType spacing = image1->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= static_cast<float>( spacing[i] );
    }

  typedef itk::Image<unsigned long, ImageDimension>                          labelimagetype;
  typedef ImageType                                                          InternalImageType;
  typedef itk::BinaryThresholdImageFilter<InternalImageType, labelimagetype> ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype> FilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, ImageType>        RelabelType;

  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  typename FilterType::Pointer filter = FilterType::New();
  typename RelabelType::Pointer relabel = RelabelType::New();

  //  InternalPixelType threshold_low, threshold_hi;

  threshold->SetInput(image1);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(0.25);
  threshold->SetUpperThreshold(1.e9);
  threshold->Update();

  filter->SetInput(threshold->GetOutput() );
  filter->SetFullyConnected( 0 );
  filter->Update();
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( smallest );
  //    relabel->SetUseHistograms(true);

  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & /* excep */ )
    {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
    }

  //  WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
  //  return 0;
  typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(relabel->GetOutput(), 0);
  // typename ImageType::Pointer Clusters=relabel->GetOutput();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  float                     maxtstat = 0;
  std::vector<unsigned int> histogram( (int)maximum + 1);
  std::vector<float>        clustersum( (int)maximum + 1);
  for( int i = 0; i <= maximum; i++ )
    {
    histogram[i] = 0;
    clustersum[i] = 0;
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      float vox = image1->GetPixel(vfIter.GetIndex() );
      histogram[(unsigned int)vfIter.Get()] = histogram[(unsigned int)vfIter.Get()] + 1;
      clustersum[(unsigned int)vfIter.Get()] += vox;
      if( vox > maxtstat )
        {
        maxtstat = vox;
        }
      }
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( vfIter.Get() > 0 )
      {
      Clusters->SetPixel( vfIter.GetIndex(), histogram[(unsigned int)vfIter.Get()]  );
      //  if ( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      //    maximgval=Clusters->GetPixel( vfIter.GetIndex());
      }
    else
      {
      Clusters->SetPixel(vfIter.GetIndex(), 0);
      }
    }

  float maximgval = 0;
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) > maximgval )
      {
      maximgval = Clusters->GetPixel( vfIter.GetIndex() );
      }
    }

  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( Clusters->GetPixel( vfIter.GetIndex() ) >= maximgval )
      {
      image1->SetPixel( vfIter.GetIndex(), 1);
      }
    else
      {
      image1->SetPixel( vfIter.GetIndex(), 0);
      }
    }

  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(image1, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int ClusterThresholdVariate(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  typedef unsigned long                                                    ULPixelType;
  typedef itk::Image<ULPixelType, ImageDimension>                          labelimagetype;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                     fIterator;
  typedef itk::ImageRegionIteratorWithIndex<labelimagetype>                labIterator;
  typedef itk::ConnectedComponentImageFilter<ImageType, labelimagetype>    FilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, labelimagetype> RelabelType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  maskfn = std::string(argv[argct]);   argct++;
  unsigned int minclustersize = 50;
  if( argc > argct )
    {
    minclustersize = std::stoi( argv[argct] );
    }
  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  typename FilterType::Pointer filter = FilterType::New();
  typename RelabelType::Pointer relabel = RelabelType::New();
  filter->SetInput( image );
  filter->SetFullyConnected( 0 );
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( 1 );
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & /* excep */ )
    {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
    }

  labIterator                vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );
  float                      maximum = relabel->GetNumberOfObjects();
  std::vector<unsigned long> histogram( (int)maximum + 1);
  for( int i = 0; i <= maximum; i++ )
    {
    histogram[i] = 0;
    }
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    float vox = vfIter.Get();
    if( vox > 0 )
      {
      if( vox > 0 )
        {
        histogram[(unsigned long)vox] = histogram[(unsigned long)vox] + 1;
        }
      }
    }

  // get the largest component's size
  unsigned long largest_component_size = 0;
  for( int i = 0; i <= maximum; i++ )
    {
    if( largest_component_size < histogram[i] )
      {
      largest_component_size = histogram[i];
      }
    }

  if(  largest_component_size < minclustersize )
    {
    minclustersize = largest_component_size - 1;
    }
//  now create the output vector
// iterate through the image and set the voxels where  countinlabel[(unsigned
// long)(labelimage->GetPixel(vfIter.GetIndex()) - min)]
// is < MinClusterSize
  unsigned long vecind = 0;
  fIterator     mIter( mask,  mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() > 0 )
      {
      float vox = mask->GetPixel(vfIter.GetIndex() );
      if( vox >= 0  )
        {
        const unsigned long clustersize =
          histogram[(unsigned long)(relabel->GetOutput()->GetPixel(mIter.GetIndex() ) )];
        if( clustersize <= minclustersize )
          {
          image->SetPixel( mIter.GetIndex(), 0 );
          }
        vecind++;
        }
      }
    }

  if( outname.length() > 3 )
    {
    WriteImage<ImageType>( image, outname.c_str() );
    }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int ExtractSlice(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>                       OutImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int slice = std::stoi(argv[argct]);   argct++;
  // std::cout << " Extract slice " << slice << " from dimension" << ImageDimension << std::endl;
  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  if( slice >= timedims )
    {
    // std::cout << " max slice number is " << timedims << std::endl;
    return 1;
    }
  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  extractRegion.SetIndex(ImageDimension - 1, slice );

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( image1 );
  //  extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();
  outimage = extractFilter->GetOutput();

  /*  typename ImageType::SpacingType qspc=warpthisimage->GetSpacing();
  typename ImageType::PointType qorg=warpthisimage->GetOrigin();
  typename ImageType::DirectionType qdir=warpthisimage->GetDirection();
  qdir.Fill(0);
  for (unsigned int qq=0; qq<ImageDimension-1; qq++) {
    for (unsigned int pp=0; pp<ImageDimension-1;pp++) {
      qdir[qq][pp]=img_mov->GetDirection()[qq][pp];
    }
      qspc[qq]=img_mov->GetSpacing()[qq];
      qorg[qq]=img_mov->GetOrigin()[qq];
    }
    warpthisimage->SetSpacing(qspc);
    warpthisimage->SetOrigin(qorg);
    warpthisimage->SetDirection(qdir);
  */
  if( outname.length() > 3 )
    {
    WriteImage<OutImageType>(outimage, outname.c_str() );
    }
  else
    {
    return 1;
    }

  return 0;
}

template <unsigned int ImageDimension>
int Finite(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  float replaceValue = 0.0;
  if( argc > 4 )
    {
    replaceValue = atof( argv[4] );
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    float val = vfIter2.Get();
    if( ! std::isfinite( val ) )
      {
      vfIter2.Set( replaceValue );
      }
    }

  WriteImage<ImageType>(image1, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int ThresholdAtMean(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       percentofmean = 1.0;
  if( argc > argct )
    {
    percentofmean = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  float        mean = 0, max = -1.e9, min = 1.e9;
  unsigned long ct = 0;
  Iterator      vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    float val = static_cast<float>( vfIter2.Get() );
    mean += val;
    if( val > max )
      {
      max = val;
      }
    else if( val < min )
      {
      min = val;
      }
    ct++;
    }
  if( ct > 0 )
    {
    mean /= static_cast<float>( ct );
    }

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilterType;
  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();

  threshold->SetInput(image1);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(mean * percentofmean);
  threshold->SetUpperThreshold(max);
  threshold->Update();
  WriteImage<ImageType>(threshold->GetOutput(), outname.c_str() );
  return 0;
}



template <unsigned int ImageDimension>
int SetTimeSpacing(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       timespacing = 1.0;
  if( argc > argct )
    {
    timespacing = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::SpacingType spacing = image1->GetSpacing();
  spacing[ ImageDimension - 1 ] = timespacing;
  image1->SetSpacing( spacing );
  WriteImage<ImageType>( image1 , outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int SetTimeSpacingWarp(int argc, char *argv[])
{
  typedef float                                              RealType;
  typedef itk::Vector<RealType, ImageDimension>              VectorType;
  typedef itk::Image<VectorType, ImageDimension>             ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       timespacing = 1.0;
  if( argc > argct )
    {
    timespacing = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::SpacingType spacing = image1->GetSpacing();
  spacing[ ImageDimension - 1 ] = timespacing;
  image1->SetSpacing( spacing );
  WriteImage<ImageType>( image1 , outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int FlattenImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       percentofmax = 1.0;
  if( argc > argct )
    {
    percentofmax = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  PixelType max = -1.e9f, min = 1.e9f;
  Iterator vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    PixelType val = vfIter2.Get();
    if( val > max )
      {
      max = val;
      }
    else if( val < min )
      {
      min = val;
      }
    }

  typename ImageType::Pointer out = MakeNewImage<ImageType>(image1, 0);
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    PixelType val = vfIter2.Get();
    if( val > max * percentofmax )
      {
      val = (max * percentofmax);
      }
    out->SetPixel(vfIter2.GetIndex(), val);
    }

  // std::cout << " Flattening to :  " << percentofmax << std::endl;
  WriteImage<ImageType>(out, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int TruncateImageIntensity( unsigned int argc, char *argv[] )
{
  typedef int   PixelType;
  typedef float RealType;

  // usage  ImageMath 3 out.nii.gz  TrunateImageIntensity InImage.nii.gz FractionLo(e.g.0.025) FractionHi(e.g.0.975)
  // Bins Mask
  if( argc < 4 )
    {
    std::cout << " need more args -- see usage   " << std::endl
              <<
      " ImageMath 3 outimage.nii.gz  TruncateImageIntensity inputImage  {lowerQuantile=0.025} {upperQuantile=0.975}  {numberOfBins=65}  {binary-maskImage} {copy-image-space-from-input-to-mask}"
              << std::endl;  throw std::exception();
    }

  unsigned int      argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       lo = 0.025;
  if( argc > argct )
    {
    lo = atof(argv[argct]);
    }
  argct++;
  float hi = 0.975;
  if( argc > argct )
    {
    hi = atof(argv[argct]);
    }
  else
    {
    hi = 1.0f - lo;
    } argct++;
  unsigned int numberOfBins = 64;
  if( argc > argct )
    {
    numberOfBins = std::stoi(argv[argct]);
    }
  argct++;

  // std::cout << " bin " << numberOfBins << " lo " << lo << " Hi " << hi << std::endl;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension>  RealImageType;

  typename RealImageType::Pointer image;
  ReadImage<RealImageType>( image, fn1.c_str() );

  typename ImageType::Pointer mask;
  if( argc > argct )
    {
    try
      {
      ReadImage<ImageType>(mask,  argv[argct]  );
      }
    catch( ... )
      {
      // std::cout << " can't read mask " << std::endl;
      mask = nullptr;
      }
    }
  argct++;
  bool copyInputSpaceToMask = false;
  if( argc > argct ) copyInputSpaceToMask = true;


  if( mask.IsNull() )
    {
    mask = AllocImage<ImageType>( image, itk::NumericTraits<PixelType>::OneValue());
    }
  if ( copyInputSpaceToMask )
    {
    mask->CopyInformation( image );
    mask->SetOrigin( image->GetOrigin() );
    mask->SetSpacing( image->GetSpacing() );
    mask->SetDirection( image->GetDirection() );
    }
  // std::cout << " iterate " << std::endl;

  itk::ImageRegionIterator<RealImageType> ItI( image,
                                               image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
                                           mask->GetLargestPossibleRegion() );

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();
  ItM.GoToBegin();
  for( ItI.GoToBegin(); !ItI.IsAtEnd();  ++ItI )
    {
    // std::cout << " ind " << ItI.GetIndex() << std::endl;
    if( ItI.Get() >  0 && ItM.Get() >= 0.5 )
      {
      if( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      ItM.Set( itk::NumericTraits<PixelType>::OneValue() );
      }
    else
      {
      ItM.Set( 0 );
      }
    if( std::isnan( ItI.Get() ) || std::isinf( ItI.Get() ) )
      {
      ItM.Set( 0 );
      }
    ++ItM;
    }
  // std::cout << " label " << std::endl;
  typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( image );
  stats->SetLabelInput( mask );
  stats->SetUseHistograms( true );
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();
  // std::cout << " labeld " << std::endl;
  typedef typename HistogramGeneratorType::HistogramType HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  float lowerQuantile = static_cast<float>( histogram->Quantile( 0, lo ) );
  float upperQuantile = static_cast<float>( histogram->Quantile( 0, hi ) );

  // std::cout << "Lower quantile: " << lowerQuantile << std::endl;
  // std::cout << "Upper quantile: " << upperQuantile << std::endl;
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if( ItI.Get() < lowerQuantile )
      {
      ItI.Set(  lowerQuantile );
      }
    if( ItI.Get() > upperQuantile )
      {
      ItI.Set(  upperQuantile  );
      }
    }

  if( outname.length() > 3 )
    {
    WriteImage<RealImageType>( image, argv[2] );
    }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int TileImages(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  unsigned int nx = std::stoi(argv[argct]);   argct++;

  unsigned int numberofimages = 0;

  typename ImageType::Pointer averageimage = nullptr;
  typename ImageType::Pointer image2 = nullptr;

  typename ImageType::SizeType size;
  size.Fill( 0 );
  typename ImageType::SizeType maxSize;
  maxSize.Fill( 0 );

  double meanval = 1;
  unsigned int bigimage = 0;
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();

    for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
      {
      itk::SizeValueType currentDimensionSize = imageIO->GetDimensions( i );
      size[i] = currentDimensionSize;

      if( currentDimensionSize > maxSize[i] )
        {
        maxSize[i] = currentDimensionSize;
        bigimage = j;
        }
      }
    }
  std::cout << " bigimage " << bigimage << " size " << maxSize << std::endl;

  ReadImage<ImageType>(image2, argv[bigimage]);

  // std::cout << " largest image " << size << std::endl;

/** declare the tiled image */
  unsigned int xsize = maxSize[0];
  unsigned int ysize = maxSize[1];
  typename ImageType::SizeType tilesize;
  unsigned int ny = static_cast<unsigned int>( static_cast<float>( numberofimages ) /
                    static_cast<float>( nx ) + 0.5f);
  if( nx * ny < numberofimages )
    {
    ny++;
    }
  // std::cout << " nx " << nx << " ny " << ny << std::endl;
  tilesize[0] = xsize * nx;
  tilesize[1] = ysize * ny;
  typename ImageType::RegionType region;
  region.SetSize( tilesize );

  bool normalizei = false;
  typename ImageType::Pointer tiledimage =
    AllocImage<ImageType>(region,
                          image2->GetSpacing(),
                          image2->GetOrigin(),
                          image2->GetDirection() );

  unsigned int imagecount = 0, imagexct = 0, imageyct = 0;
  for( unsigned int j = argct; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    ReadImage<ImageType>(image2, fn.c_str() );

    unsigned long ct = 0;
    if( normalizei )
      {
      meanval = 0.0;
      Iterator vfIter2( image2,  image2->GetLargestPossibleRegion() );
      for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
        {
        meanval += static_cast<double>( vfIter2.Get() );
        ct++;
        }
      if( ct > 0 )
        {
        meanval /= static_cast<double>( ct );
        }
      if( meanval <= 0 )
        {
        meanval = 1.0;
        }
      }

    imagexct = imagecount % nx;
    imageyct = imagecount / nx;
    // std::cout << "doing " << fn << "  " << imagecount << " x " << imagexct <<  " y " << imageyct << std::endl;
    imagecount++;
    Iterator vfIter( image2,  image2->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      typename ImageType::IndexType locind = vfIter.GetIndex();
      typename ImageType::IndexType globind;
      globind[0] = size[0] * imagexct + locind[0];
      globind[1] = size[1] * imageyct + locind[1];
      double val =  static_cast<double>( vfIter.Get() ) / meanval;
      tiledimage->SetPixel(globind,   val );
      }
    }

  WriteImage<ImageType>(tiledimage, outname.c_str() );
  return 0;
  typedef itk::Image<unsigned char, ImageDimension>                  ByteImageType;
  typedef itk::RescaleIntensityImageFilter<ImageType, ByteImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( tiledimage );

  // std::cout << " writing output ";
  WriteImage<ByteImageType>( rescaler->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int ConvertLandmarkFile(unsigned int argc, char *argv[])
{
  unsigned int argct = 2;

  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string infn = std::string(argv[argct]); argct++;
  float       pointp = 1;

  typedef itk::PointSet<long, ImageDimension>          PointSetType;
  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( infn.c_str() );
  reader->SetRandomPercentage( 1 );
  if( pointp > 0 && pointp < 1  )
    {
    reader->SetRandomPercentage( pointp );
    }
  reader->Update();

  // std::cout << "Number of labels: " << reader->GetNumberOfLabels() << std::endl;
  // std::cout << "Labels: ";
  for( unsigned int i = 0; i < reader->GetNumberOfLabels(); i++ )
    {
    // std::cout << reader->GetLabelSet()->operator[](i) << " ";

    }
  // std::cout << std::endl;

  typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outname.c_str() );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return 0;
}

template <unsigned int ImageDimension>
int TriPlanarView(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  unsigned int argct = 2;
  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string maskfn = std::string(argv[argct]); argct++;
  // std::cout << " file name " << maskfn << std::endl;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  //  WriteImage<ImageType>(mask,"temp.nii");
  float clamppercent1 = 0.1;
  if( argc > argct )
    {
    clamppercent1 = atof(argv[argct]);
    }
  argct++;
  if( clamppercent1 > 1 )
    {
    clamppercent1 = 1;
    }
  float clamppercent2 = 0.1;
  if( argc > argct )
    {
    clamppercent2 = atof(argv[argct]);
    }
  argct++;
  if( clamppercent2 > 1 )
    {
    clamppercent2 = 1;
    }

  typename ImageType::SizeType size = mask->GetLargestPossibleRegion().GetSize();
  unsigned int xslice = size[0] / 2;
  if( argc > argct )
    {
    xslice = std::stoi(argv[argct]);
    }
  argct++;
  unsigned int yslice = size[1] / 2;
  if( argc > argct )
    {
    yslice = std::stoi(argv[argct]);
    }
  argct++;
  unsigned int zslice = size[2] / 2;
  if( argc > argct )
    {
    zslice = std::stoi(argv[argct]);
    }
  argct++;

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( mask );
  rescaler->Update();
  mask = rescaler->GetOutput();

  //  typedef itk::IntensityWindowingImageFilter<ImageType,ImageType> wFilterType;
  //  typename wFilterType::Pointer wer = wFilterType::New();
  // wer->SetInput(mask);

/** declare the tiled image */
  unsigned long xsize = size[0];
  unsigned long ysize = size[1];
  unsigned long zsize = size[2];
  typename MatrixImageType::SizeType ztilesize;
  ztilesize[0] = xsize;
  ztilesize[1] = ysize;
  typename MatrixImageType::SizeType ytilesize;
  ytilesize[0] = xsize;
  ytilesize[1] = zsize;
  typename MatrixImageType::SizeType xtilesize;
  xtilesize[0] = ysize;
  xtilesize[1] = zsize;
  typename MatrixImageType::SizeType tilesize;
  tilesize[0] = xtilesize[0] + ytilesize[0] + ztilesize[0];
  tilesize[1] = xtilesize[1];
  if( ytilesize[1] > tilesize[1] )
    {
    tilesize[1] = ytilesize[1];
    }
  if( ztilesize[1] > tilesize[1] )
    {
    tilesize[1] = ztilesize[1];
    }
  // std::cout << " allocate matrix " << tilesize << std::endl;
  typename MatrixImageType::RegionType region;
  region.SetSize( tilesize );

  typename MatrixImageType::Pointer matimage =
    AllocImage<MatrixImageType>(region);

  unsigned int lowgetridof = (unsigned int) (clamppercent1 * 256);
  unsigned int higetridof = (unsigned int) (256 - clamppercent2 * 256);
  // std::cout << " get rid of " << getridof << std::endl;
  matimage->FillBuffer(lowgetridof);
  // now loop over each slice and put the pixels in the right place in matimage
  typename MatrixImageType::IndexType index2d;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter2( mask,  mask->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    double val1 = vfIter2.Get();
    if( val1 > higetridof )
      {
      vfIter2.Set(higetridof );
      }
    if( val1 < lowgetridof   )
      {
      vfIter2.Set(lowgetridof );
      }
    // first do z-slice
    if( vfIter2.GetIndex()[2] == (long) zslice )
      {
      double val = vfIter2.Get();
      typename ImageType::IndexType index = vfIter2.GetIndex();
      index2d[0] = index[0] + xtilesize[0] + ytilesize[0];
      index2d[1] = index[1];
      index2d[1] = tilesize[1] - index2d[1] - 1;
      matimage->SetPixel(index2d, val);
      }
    if( vfIter2.GetIndex()[1] == (long)yslice )
      {
      double val = vfIter2.Get();
      typename ImageType::IndexType index = vfIter2.GetIndex();
      index2d[0] = index[0] + xtilesize[0];
      index2d[1] = index[2];
      index2d[1] = tilesize[1] - index2d[1] - 1;
      matimage->SetPixel(index2d, val);
      }
    if( vfIter2.GetIndex()[0] == (long)xslice )
      {
      double val = vfIter2.Get();
      typename ImageType::IndexType index = vfIter2.GetIndex();
      index2d[0] = index[1];
      index2d[1] = index[2];
      index2d[1] = tilesize[1] - index2d[1] - 1;
      matimage->SetPixel(index2d, val);
      }
    }

  typedef itk::Image<unsigned char, 2>                                     ByteImageType;
  typedef itk::RescaleIntensityImageFilter<MatrixImageType, ByteImageType> RescaleFilterType2;
  typename RescaleFilterType2::Pointer rescaler2 = RescaleFilterType2::New();
  rescaler2->SetOutputMinimum(   0 );
  rescaler2->SetOutputMaximum( 255 );
  rescaler2->SetInput( matimage );
  rescaler2->Update();
  // std::cout << " writing output ";
  WriteImage<ByteImageType>( rescaler2->GetOutput(), outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int ConvertVectorToImage(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  typedef itk::ImageRegionIteratorWithIndex<MatrixImageType>              vIterator;

  int argct = 2;
  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string maskfn = std::string(argv[argct]); argct++;
  std::string vecfn = std::string(argv[argct]); argct++;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  typename MatrixImageType::Pointer vecimg = nullptr;
  ReadImage<MatrixImageType>(vecimg, vecfn.c_str() );
  unsigned long voxct = 0, mct = 0;
  Iterator      mIter( mask, mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5f )
      {
      mct++;
      }
    }
  vIterator vIter(vecimg, vecimg->GetLargestPossibleRegion() );
  for(  vIter.GoToBegin(); !vIter.IsAtEnd(); ++vIter )
    {
    voxct++;
    }

  // std::cout << " vct " << voxct << " mct " << mct << std::endl;

  typename ImageType::Pointer outimage = nullptr;
  ReadImage<ImageType>(outimage, maskfn.c_str() );
  outimage->FillBuffer(0);

  vIter.GoToBegin();
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5f )
      {
      outimage->SetPixel(mIter.GetIndex(), vIter.Get() );
      ++vIter;
      }
    }

  WriteImage<ImageType>(outimage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int CorruptImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       noiselevel = 1.0;
  if( argc > argct )
    {
    noiselevel = atof(argv[argct]); argct++;
    }
  float smoothlevel = 1.0;
  if( argc > argct )
    {
    smoothlevel = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    double r = (   (double)rand()  / ( (double)(RAND_MAX) + (double)(1) ) ) - 0.5;
    iter.Set(iter.Get() + static_cast<PixelType>( r ) * static_cast<PixelType>( noiselevel ) );
    }

  if( smoothlevel >  0 )
    {
      {
      typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
      typename dgf::Pointer filter = dgf::New();
      filter->SetVariance(smoothlevel);
      filter->SetUseImageSpacingOff();
      filter->SetMaximumError(.01f);
      filter->SetInput(image1);
      filter->Update();
      image1 = filter->GetOutput();
      }
    }
  WriteImage<ImageType>(image1, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int Where(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 4;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       value = atof(argv[argct]); argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  float tol = 0.0;
  if( argc > argct )
    {
    tol = atof(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }

  unsigned long ct = 0;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    if( !image2 )
      {
      if( static_cast<float>( std::fabs(iter.Get() - value) ) < tol )
        {
        // std::cout << iter.GetIndex() << std::endl;
        ct++;
        }
      }
    else if( image2->GetPixel(iter.GetIndex() ) > 0 && static_cast<float>( fabs(iter.Get() - value) ) < tol )
      {
      // std::cout << iter.GetIndex() << std::endl;
      ct++;
      }
    }
  // std::cout << ct <<  " voxels have the value " << value << std::endl;
  return 0;
}

template <unsigned int ImageDimension>
int SetOrGetPixel(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       value = atof(argv[argct]);
  std::string Get = std::string(argv[argct]); argct++;
  float       indx = atof(argv[argct]); argct++;
  float       indy = 0;
  if( ImageDimension >= 2 )
    {
    indy = atof(argv[argct]); argct++;
    }
  float indz = 0;
  if( ImageDimension >= 3 )
    {
    indz = atof(argv[argct]); argct++;
    }
  bool usephyspace = false;
  if( argc > argct )
    {
    usephyspace = std::stoi(argv[argct]); argct++;
    }
  bool get = false;
  if( strcmp(Get.c_str(), "Get") == 0 )
    {
    get = true;
    }
  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    ReadImage<ImageType>(image2, fn1.c_str() );
    }
  if( !image1 )
    {
    // std::cout << " no image ! " << std::endl; throw std::exception();
    }

  typename ImageType::IndexType index;
  index.Fill(0);
  if( usephyspace == false )
    {
    index[0] = (long int)indx;
    index[1] = (long int)indy;
    if( ImageDimension == 3 )
      {
      index[2] = (long int)indz;
      }
    }
  else
    {
    typename ImageType::PointType porig;
    porig[0] = indx;
    porig[1] = indy;
    if( ImageDimension == 3 )
      {
      porig[2] = indz;
      }
    image1->TransformPhysicalPointToIndex(porig, index);
    }
  // std::cout << " use phy " << usephyspace << " " << indx << " " << indy << " " << indz << std::endl;
  // std::cout << " Ind " << index << std::endl;
  bool isinside = true;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    float shifted = index[i];
    if( shifted < 0 || shifted >  image1->GetLargestPossibleRegion().GetSize()[i] - 1  )
      {
      isinside = false;
      }
    }

  if( isinside == true )
    {
    if( get )
      {
      // std::cout << " GetValue at " << index << " is " << image1->GetPixel(index) << std::endl;
      // std::cout << image1->GetPixel(index) << std::endl;
      }
    else
      {
      // std::cout << " SetValue at " << index << " value " << value << " replaces " <<  image1->GetPixel(index)
      //         << std::endl;
      image2->SetPixel(index, value);
      WriteImage<ImageType>(image2, outname.c_str() );
      }
    }
  else
    {
    // std::cout << "NA" << std::endl;
    }

  return 0;
}

template <unsigned int ImageDimension>
int HistogramMatching(int argc, char * argv[])
{
  typedef float                                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension>                   ImageType;
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> MatchingFilterType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;
  long        bins = 255;
  if( argc > argct )
    {
    bins = std::stoi( argv[argct] );
    }
  argct++;
  long points = 64;
  if( argc > argct )
    {
    points = std::stoi( argv[argct] );
    }
  argct++;
  bool useThresholdAtMeanIntensity = false;
  if( argc > argct )
    {
    useThresholdAtMeanIntensity = static_cast<bool>( std::stoi( argv[argct] ) );
    }
  argct++;

  typename ImageType::Pointer source;
  ReadImage<ImageType>( source, fn1.c_str() );

  typename ImageType::Pointer reference;
  ReadImage<ImageType>( reference, fn2.c_str() );

  typename MatchingFilterType::Pointer match = MatchingFilterType::New();
  match->SetSourceImage( source );
  match->SetReferenceImage( reference );
  match->SetNumberOfHistogramLevels( bins );
  match->SetThresholdAtMeanIntensity( useThresholdAtMeanIntensity );
  match->SetNumberOfMatchPoints( points );
  match->Update();

  WriteImage<ImageType>( match->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int RescaleImage( int argc, char * argv[] )
{
  if( argc < 7 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;
    throw std::exception();
    }

  typedef float                                                  PixelType;
  typedef itk::Image<PixelType, ImageDimension>                  ImageType;

  // Usage:  ImageMath 3 output.nii.gz RescaleImage input.nii.gz min max

  const std::string outputName = std::string( argv[2] );
  const std::string inputName = std::string( argv[4] );
  const PixelType   min = static_cast<PixelType>( atof( argv[5] ) );
  const PixelType   max = static_cast<PixelType>( atof( argv[6] ) );

  typename ImageType::Pointer input;
  ReadImage<ImageType>( input, inputName.c_str() );

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum( min );
  rescaler->SetOutputMaximum( max );
  rescaler->SetInput( input );
  rescaler->Update();

  WriteImage<ImageType>( rescaler->GetOutput(), outputName.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int WindowImage( int argc, char * argv[] )
{
  if( argc < 9 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;
    throw std::exception();
    }

  typedef float                                                  PixelType;
  typedef itk::Image<PixelType, ImageDimension>                  ImageType;

  // Usage:  ImageMath 3 output.nii.gz RescaleImage input.nii.gz min max

  const std::string outputName = std::string( argv[2] );
  const std::string inputName = std::string( argv[4] );
  const PixelType   inputMin = static_cast<PixelType>( atof( argv[5] ) );
  const PixelType   inputMax = static_cast<PixelType>( atof( argv[6] ) );
  const PixelType   outputMin = static_cast<PixelType>( atof( argv[7] ) );
  const PixelType   outputMax = static_cast<PixelType>( atof( argv[8] ) );

  typename ImageType::Pointer input;
  ReadImage<ImageType>( input, inputName.c_str() );

  typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingFilterType;
  typename IntensityWindowingFilterType::Pointer rescaler = IntensityWindowingFilterType::New();
  rescaler->SetWindowMinimum( inputMin );
  rescaler->SetWindowMaximum( inputMax );
  rescaler->SetOutputMinimum( outputMin );
  rescaler->SetOutputMaximum( outputMax );
  rescaler->SetInput( input );
  rescaler->Update();

  WriteImage<ImageType>( rescaler->GetOutput(), outputName.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int NeighborhoodStats( int itkNotUsed( argc ), char * argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::VectorImage<PixelType, ImageDimension>                                              VectorImageType;
  typedef itk::FlatStructuringElement<ImageDimension>                                              KernelType;
  typedef itk::NeighborhoodFirstOrderStatisticsImageFilter<ImageType, VectorImageType, KernelType> TextureFilterType;

  const std::string  outputName = std::string( argv[2] );
  const std::string  inputName = std::string( argv[4] );
  const unsigned int whichStat = static_cast<unsigned int>( std::stoi( argv[5] ) );
  const unsigned int rad = static_cast<unsigned int>( std::stoi( argv[6] ) );

  typename ImageType::Pointer input;
  ReadImage<ImageType>( input, inputName.c_str() );

  typename KernelType::SizeType radius;
  radius.Fill( rad );
  KernelType kernel = KernelType::Box( radius );

  typename TextureFilterType::Pointer filter = TextureFilterType::New();
  filter->SetKernel( kernel );
  filter->SetInput( input );
  filter->Update();

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetInput( filter->GetOutput() );

  switch( whichStat )
    {
    case 0:
      {
      indexSelectionFilter->SetIndex( 0 );
      break;
      }
    case 1:
      {
      indexSelectionFilter->SetIndex( 1 );
      break;
      }
    case 2:
      {
      indexSelectionFilter->SetIndex( 2 );
      break;
      }
    case 3:
      {
      indexSelectionFilter->SetIndex( 3 );
      break;
      }
    case 4:
      {
      indexSelectionFilter->SetIndex( 4 );
      break;
      }
    case 5:
      {
      indexSelectionFilter->SetIndex( 5 );
      break;
      }
    case 6:
      {
      indexSelectionFilter->SetIndex( 6 );
      break;
      }
    case 7:
      {
      indexSelectionFilter->SetIndex( 7 );
      break;
      }
    default:
      {
      std::cerr << "Unrecognized option: " << whichStat << std::endl;
      return EXIT_FAILURE;
      }
    }
  indexSelectionFilter->Update();

  WriteImage<ImageType>( indexSelectionFilter->GetOutput(), outputName.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int PadImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);
  argct++;
  const float padvalue = atof(argv[argct]);
  argct++;

  PixelType padVoxelValue = itk::NumericTraits<PixelType>::Zero;
  if( argc > 6 )
    {
    padVoxelValue = static_cast<PixelType>( atof( argv[argct] ) );
    }

  typename ImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }

  typename ImageType::PointType origin2 = image1->GetOrigin();

  typename ImageType::SizeType size = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::SizeType newsize = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::RegionType newregion;
  // determine new image size
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    float dimsz = (float)size[i];
    newsize[i] = (unsigned int)(dimsz + padvalue * 2);
    }
  newregion.SetSize(newsize);
  newregion.SetIndex(image1->GetLargestPossibleRegion().GetIndex() );

  typename ImageType::Pointer padimage =
    AllocImage<ImageType>(newregion,
                          image1->GetSpacing(),
                          origin2,
                          image1->GetDirection(), padVoxelValue);

  typename ImageType::IndexType index;
  typename ImageType::IndexType index2;
  if( padvalue > 0 )
    {
    index.Fill(0);
    index2.Fill( (unsigned int)fabs(padvalue) );
    }
  else
    {
    index2.Fill(0);
    index.Fill( (unsigned int)fabs(padvalue) );
    }

  typename ImageType::PointType point1, pointpad;
  image1->TransformIndexToPhysicalPoint(index, point1);
  padimage->TransformIndexToPhysicalPoint(index2, pointpad);

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin2[i] += (point1[i] - pointpad[i]);
    }

  padimage->SetOrigin(origin2);

  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    typename ImageType::IndexType oindex = iter.GetIndex();
    typename ImageType::IndexType padindex = iter.GetIndex();

    bool isinside = true;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      float shifted = ( (float)oindex[i] + padvalue);
      if( shifted < 0 || shifted > newsize[i] - 1 )
        {
        isinside = false;
        }
      //      if (shifted < 0) shifted=0;
      // padindex[i]=
      }
    if( isinside )
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        float shifted = ( (float)oindex[i] + padvalue);
        padindex[i] = (unsigned int)shifted;
        }
      padimage->SetPixel(padindex, iter.Get() );
      }
    }

  WriteImage<ImageType>(padimage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int SigmoidImage(int argc, char *argv[])
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  int               argct = 2;
  const std::string outname = std::string( argv[argct++] );
  argct++;
  std::string inputFilename = std::string( argv[argct++] );

  double alpha = 1.0;
  if( argc > argct )
    {
    alpha = atof( argv[argct++] );
    }
  double beta = 0.0;
  if( argc > argct )
    {
    beta = atof( argv[argct++] );
    }

  typename ImageType::Pointer inputImage = nullptr;
  if( inputFilename.length() > 3 )
    {
    ReadImage<ImageType>( inputImage, inputFilename.c_str() );
    }

  typedef itk::SigmoidImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetAlpha( alpha );
  filter->SetBeta( beta );
  filter->SetOutputMinimum( 0.0 );
  filter->SetOutputMaximum( 1.0 );
  filter->Update();

  WriteImage<ImageType>( filter->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int CoordinateComponentImages( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  if( argc < 4 )
    {
    std::cerr << "Error:  Not enough args.  See help.";
    return EXIT_FAILURE;
    }

  int               argct = 2;
  const std::string outPrefix = std::string( argv[argct++] );
  argct++;
  std::string domainImageFile = std::string( argv[argct++] );

  typename ImageType::Pointer domainImage = nullptr;
  if( domainImageFile.length() > 3 )
    {
    ReadImage<ImageType>( domainImage, domainImageFile.c_str() );
    }

  std::vector<typename ImageType::Pointer> coordinateImages( ImageDimension );
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    coordinateImages[d] = ImageType::New();
    coordinateImages[d]->CopyInformation( domainImage );
    coordinateImages[d]->SetRegions( domainImage->GetRequestedRegion() );
    coordinateImages[d]->Allocate();
    }

  itk::ImageRegionIteratorWithIndex<ImageType> It( domainImage,
    domainImage->GetRequestedRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::PointType imagePoint;
    domainImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      coordinateImages[d]->SetPixel( It.GetIndex(), imagePoint[d] );
      }
    }

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    std::string outputFile = outPrefix + std::to_string( d ) + std::string( ".nii.gz" );
    WriteImage<ImageType>( coordinateImages[d], outputFile.c_str() );
    }
  return 0;
}

template <unsigned int ImageDimension>
int SharpenImage(int argc, char *argv[])
{
  if( argc < 5 )
    {
    // std::cout << "Error.  Not enough arguments.  See help menu." << std::endl;
    return EXIT_FAILURE;
    }

  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  const std::string outputFilename = std::string( argv[2] );
  const std::string inputFilename = std::string( argv[4] );

  typename ImageType::Pointer inputImage = nullptr;
  ReadImage<ImageType>( inputImage, inputFilename.c_str() );

  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->Update();

  WriteImage<ImageType>( filter->GetOutput(), outputFilename.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int CenterImage2inImage1(int argc, char *argv[])
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  image1->FillBuffer(0);
  typename ImageType::Pointer image2 = nullptr;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }

  // compute center of mass in image2 in point space
  double iweight = 0;
  typename ImageType::PointType cm_point; cm_point.Fill(0);
  Iterator iter( image2,  image2->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    typename ImageType::PointType point;
    image2->TransformIndexToPhysicalPoint(iter.GetIndex(), point);
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      cm_point[d] += point[d] * static_cast<double>( iter.Get() );
      }
    iweight += static_cast<double>( iter.Get() );
    }

  // center of image1
  typename ImageType::IndexType image_1_center_index;
  image_1_center_index.Fill(0);
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    cm_point[d] /= iweight;
    image_1_center_index[d] = image1->GetLargestPossibleRegion().GetSize()[d] / 2;
    }
  typename ImageType::PointType image_1_center_point;
  image1->TransformIndexToPhysicalPoint(image_1_center_index, image_1_center_point);

  // now we translate the cm_point to the center of image1
  typename ImageType::PointType trans;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    trans[d] = image_1_center_point[d] - cm_point[d];
    }
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    typename ImageType::PointType point;
    image2->TransformIndexToPhysicalPoint(iter.GetIndex(), point);
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = point[d] + trans[d];
      }

    typename ImageType::IndexType newindex;
    newindex.Fill(0);
    bool isinside = image1->TransformPhysicalPointToIndex(point, newindex);
    if( isinside )
      {
      image1->SetPixel(newindex, iter.Get() );
      }
    }

  WriteImage<ImageType>(image1, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesMask( int argc, char *argv[] )
{

  if ( argc <= 5 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                            PixelType;
  typedef itk::Image<PixelType, ImageDimension>            ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>        MaskImageType;
  typedef itk::ImageRegionIteratorWithIndex<MaskImageType> Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer timeseries = ImageType::New();
  ReadImage<ImageType>( timeseries, fn1.c_str() );

  typename MaskImageType::Pointer mask = MaskImageType::New();
  ReadImage<MaskImageType>( mask, fn2.c_str() );

  Iterator it( mask, mask->GetLargestPossibleRegion() );
  while ( ! it.IsAtEnd() )
    {
    if( itk::Math::FloatAlmostEqual( it.Value(), itk::NumericTraits<PixelType>::ZeroValue() ) )// reference
      {
      typename MaskImageType::IndexType maskIdx = it.GetIndex();
      typename ImageType::IndexType timeIdx;

      for (unsigned int i=0; i<(ImageDimension-1); i++)
        {
        timeIdx[i] = maskIdx[i];
        }

      for (unsigned int t=0; t<timeseries->GetLargestPossibleRegion().GetSize()[ImageDimension-1]; t++)
        {
        timeIdx[ImageDimension-1] = t;
        timeseries->SetPixel(timeIdx, 0);
        }

      }

    ++it;
    }

  WriteImage<ImageType>(timeseries, outname.c_str());
  return 0;

}

template <unsigned int ImageDimension>
int TimeSeriesDisassemble(int argc, char *argv[])
{
  if( argc <= 4 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename OutImageType::PointType outOrigin;
  for( unsigned int i = 0; i < (ImageDimension - 1); i++ )
    {
    outOrigin[i] = image1->GetOrigin()[i];
    }

  unsigned int n_sub_vols = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];

  // Extract filename while allowing directory names with '.' in them
  // (cluster temp dirs)
  std::string::size_type idx = outname.find_last_of('/');
  std::string            dirname = outname.substr(0, idx + 1);
  std::string            filename = outname.substr(idx + 1);
  std::string::size_type idx2 = filename.find_first_of('.');
  std::string            tempname = filename.substr(0, idx2);
  std::string            extension = filename.substr(idx2);
  for( unsigned int i = 0; i < n_sub_vols; i++ )
    {
    std::string       s;
    std::stringstream out;
    out << (1000 + i);
    s = out.str();
    std::string kname = dirname + tempname + s + extension;

    typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
    extractRegion.SetSize(ImageDimension - 1, 0);
    extractRegion.SetIndex(ImageDimension - 1, i );

    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput( image1 );
    extractFilter->SetDirectionCollapseToSubmatrix();
    extractFilter->SetExtractionRegion( extractRegion );
    extractFilter->Update();
    outimage = extractFilter->GetOutput();
    outimage->SetOrigin( outOrigin );
    WriteImage<OutImageType>(outimage, kname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesAssemble(int argc, char *argv[])
{
  if( argc <= 6 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension - 1>    ImageType;
  typedef itk::Image<PixelType, ImageDimension>        OutImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  float time = atof( argv[argct] );           argct++;
  float origin = atof( argv[argct] );         argct++;

  typename OutImageType::Pointer outimage = OutImageType::New();

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;

  for( int i = 6; i < argc; i++ )
    {
    typename ImageType::Pointer image1 = nullptr;
    ReadImage<ImageType>(image1, argv[i] );

    if( i == 6 )
      {
      typename OutImageType::SizeType outSize;
      typename OutImageType::SpacingType outSpacing;
      typename OutImageType::PointType outOrigin;
      typename OutImageType::DirectionType outDirection;
      for( unsigned int d = 0; d < (ImageDimension - 1); d++ )
        {
        outSize[d] = image1->GetLargestPossibleRegion().GetSize()[d];
        outSpacing[d] = image1->GetSpacing()[d];
        outOrigin[d] = image1->GetOrigin()[d];
        for( unsigned int e = 0; e < (ImageDimension - 1); e++ )
          {
          outDirection(e, d) = image1->GetDirection() (e, d);
          }
        }
      for( unsigned int d = 0; d < (ImageDimension - 1); d++ )
        {
        outDirection(d, ImageDimension - 1) = 0;
        outDirection(ImageDimension - 1, d) = 0;
        }
      outDirection(ImageDimension - 1, ImageDimension - 1) = 1.0;

      outSize[ImageDimension - 1] = argc - 6;
      outSpacing[ImageDimension - 1] = time;
      outOrigin[ImageDimension - 1] = origin;

      typename OutImageType::RegionType outRegion;
      outRegion.SetSize( outSize );
      outimage->SetRegions( outRegion );
      outimage->SetSpacing( outSpacing );
      outimage->SetOrigin( outOrigin );
      outimage->SetDirection( outDirection );
      outimage->Allocate();
      }

    ImageIt it( image1, image1->GetLargestPossibleRegion() );
    while( !it.IsAtEnd() )
      {
      typename OutImageType::IndexType index;
      for( unsigned int d = 0; d < (ImageDimension - 1); d++ )
        {
        index[d] = it.GetIndex()[d];
        }
      index[ImageDimension - 1] = i - 6;
      outimage->SetPixel(index, it.Value() );
      ++it;
      }
    }

  WriteImage<OutImageType>( outimage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesSubset(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int n_sub_vols = std::stoi(argv[argct]);   argct++;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    // std::cout << "Failed to read input image" << std::endl;
    return 1;
    }

  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  float        step = (float)timedims / (float)n_sub_vols;
  if( n_sub_vols >= timedims )
    {
    n_sub_vols = timedims; step = 1;
    }
  for( unsigned int i = 0; i < n_sub_vols; i++ )
    {
    std::string       s;
    std::stringstream out;
    out << (100 + i);
    s = out.str();
    std::string kname = tempname + s + extension;
    typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
    extractRegion.SetSize(ImageDimension - 1, 0);
    unsigned int sub_vol = (unsigned int)( (float)i * step);
    if( sub_vol >= timedims )
      {
      sub_vol = timedims - 1;
      }
    extractRegion.SetIndex(ImageDimension - 1, sub_vol );

    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput( image1 );
    //    extractFilter->SetDirectionCollapseToIdentity();
    extractFilter->SetDirectionCollapseToSubmatrix();
    extractFilter->SetExtractionRegion( extractRegion );
    extractFilter->Update();
    outimage = extractFilter->GetOutput();
    WriteImage<OutImageType>(outimage, kname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesSimpleSubtraction(int argc, char *argv[])
{
  typedef float                                     PixelType;
  typedef itk::Image<PixelType, ImageDimension>     InputImageType;
  typedef itk::Image<PixelType, ImageDimension - 1> OutputImageType;

  typedef itk::AlternatingValueSimpleSubtractionImageFilter<InputImageType, InputImageType>
    ImageFilterType;
  typedef itk::AverageOverDimensionImageFilter<InputImageType, OutputImageType>
    MeanFilterType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string fn1 = std::string(argv[argct++]);
  bool        mean = false;

  typename ImageFilterType::Pointer filter = ImageFilterType::New();

  if( argc >= 6 )
    {
    if( std::stoi(argv[argct++]) > 0 )
      {
      mean = true;
      }
    }

  typename InputImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<InputImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  filter->SetInput( image1 );
  filter->Update();

  if( mean )
    {
    typename MeanFilterType::Pointer meanFilter = MeanFilterType::New();
    meanFilter->SetInput( filter->GetOutput() );
    meanFilter->SetAveragingDimension( ImageDimension - 1 );
    meanFilter->SetDirectionCollapseToSubmatrix();
    meanFilter->Update();
    WriteImage<OutputImageType>(meanFilter->GetOutput(), outname.c_str() );
    }
  else
    {
    WriteImage<InputImageType>(filter->GetOutput(), outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesInterpolationSubtraction(int argc, char *argv[])
{
  typedef float                                     PixelType;
  typedef itk::Image<PixelType, ImageDimension>     InputImageType;
  typedef itk::Image<PixelType, ImageDimension - 1> OutputImageType;

  typedef itk::AlternatingValueDifferenceImageFilter<InputImageType, InputImageType>
    ImageFilterType;
  typedef itk::AverageOverDimensionImageFilter<InputImageType, OutputImageType>
    MeanFilterType;

  typedef itk::BSplineInterpolateImageFunction<InputImageType, double> BSplineInterpolatorType;
  typedef typename BSplineInterpolatorType::Pointer                    BSplineInterpolatorPointerType;

  typedef itk::WindowedSincInterpolateImageFunction<InputImageType, 4> SincInterpolatorType;
  typedef typename SincInterpolatorType::Pointer                       SincInterpolatorPointerType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string fn1 = std::string(argv[argct++]);

  typename ImageFilterType::Pointer filter = ImageFilterType::New();
  filter->SetSubtractionDimension( ImageDimension - 1 );

  if( argc >= 6 )
    {
    std::string interp = argv[argct++];
    if( strcmp( "sinc", interp.c_str() ) == 0 )
      {
      constexpr unsigned int SincRadius = 4;
      // std::cout << "Using sinc interpolation" << std::endl;
      SincInterpolatorPointerType labelInterp = SincInterpolatorType::New();
      SincInterpolatorPointerType controlInterp = SincInterpolatorType::New();
      filter->SetControlInterpolator( controlInterp );
      filter->SetLabelInterpolator( labelInterp );
      filter->SetIndexPadding( SincRadius );
      }
    else if( strcmp( "bspline", interp.c_str() ) == 0 )
      {
      // std::cout << "Using bspline interpolation of order 3" << std::endl;
      BSplineInterpolatorPointerType labelInterpB = BSplineInterpolatorType::New();
      labelInterpB->SetSplineOrder( 3 );
      BSplineInterpolatorPointerType controlInterpB = BSplineInterpolatorType::New();
      controlInterpB->SetSplineOrder( 3 );
      filter->SetControlInterpolator( controlInterpB );
      filter->SetLabelInterpolator( labelInterpB );
      filter->SetIndexPadding( 1 );
      }
    else
      {
      // std::cout << "Using linear interpolation" << std::endl;
      }

    }

  bool mean = false;
  if( argc >= 7 )
    {
    if( std::stoi(argv[argct++]) > 0 )
      {
      mean = true;
      }
    }

  typename InputImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<InputImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  filter->SetInput( image1 );
  filter->Update();
  if( mean )
    {
    typename MeanFilterType::Pointer meanFilter = MeanFilterType::New();
    meanFilter->SetInput( filter->GetOutput() );
    meanFilter->SetAveragingDimension( ImageDimension - 1 );
    meanFilter->SetDirectionCollapseToSubmatrix();
    meanFilter->Update();
    WriteImage<OutputImageType>(meanFilter->GetOutput(), outname.c_str() );
    }
  else
    {
    WriteImage<InputImageType>(filter->GetOutput(), outname.c_str() );
    }

  if( argc >= 8  )
    {
    std::string control_out = argv[argct++];
    WriteImage<InputImageType>(filter->GetModifiableControlOutputImage(), control_out.c_str() );
    }
  if( argc >= 8  )
    {
    std::string label_out = argv[argct++];
    WriteImage<InputImageType>(filter->GetModifiableLabelOutputImage(), label_out.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int SplitAlternatingTimeSeries(int argc, char *argv[])
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::SplitAlternatingTimeSeriesImageFilter<ImageType, ImageType>
    ImageFilterType;

  if( argc < 5 )
    {
    // std::cout << "Usage: ImageMath 4 split.nii.gz AlternatingTimeSeriesExtraction time.nii.gz" << std::endl;
    return 1;
    }

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  std::string::size_type idx;
  idx = outname.find_first_of('.');

  std::string basename = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  std::string zero( "0" );
  std::string one( "1" );

  std::string outname0 = basename + zero + extension;
  std::string outname1 = basename + one + extension;

  typename ImageFilterType::Pointer filter = ImageFilterType::New();
  typename ImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  filter->SetInput( image1 );
  filter->Update();

  WriteImage<ImageType>(filter->GetOutput(0), outname0.c_str() );
  WriteImage<ImageType>(filter->GetOutput(1), outname1.c_str() );
  return 0;

}

template <unsigned int ImageDimension>
int SliceTimingCorrection(int argc, char *argv[])
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;
  typedef itk::Image<PixelType, ImageDimension> OutputImageType;

  typedef itk::SliceTimingCorrectionImageFilter<InputImageType, OutputImageType> ImageFilterType;

  typedef itk::BSplineInterpolateImageFunction<InputImageType, double> BSplineInterpolatorType;
  typedef typename BSplineInterpolatorType::Pointer                    BSplineInterpolatorPointerType;

  typedef itk::WindowedSincInterpolateImageFunction<InputImageType, ImageDimension - 1> SincInterpolatorType;
  typedef typename SincInterpolatorType::Pointer                                        SincInterpolatorPointerType;

  if( argc < 5 )
    {
    std::cout
      <<
      "Usage: ImageMath 4 out.nii.gz SliceTimingCorrection input.nii.gz [sliceTiming] [sinc/bspline] [sincRadius=4/bsplineOrder=3]"
      << std::endl;
    return 1;
    }

  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string fn1 = std::string(argv[argct++]);

  float sliceTiming = 0;

  if( argc > 5 )
    {
    sliceTiming = atof( argv[argct++] );
    }

  typename ImageFilterType::Pointer filter = ImageFilterType::New();
  filter->SetTimeDimension( ImageDimension - 1 );

  if( argc >= 7 )
    {
    std::string interp = argv[argct++];
    if( strcmp( "sinc", interp.c_str() ) == 0 )
      {
      unsigned int sincRadius = 4;
      if( argc >= 8 )
        {
        sincRadius = std::stoi( argv[argct++] );
        }
      // std::cout << "Using sinc interpolation of radius " << sincRadius << std::endl;

      SincInterpolatorPointerType sInterp = SincInterpolatorType::New();
      filter->SetInterpolator( sInterp );
      filter->SetIndexPadding( sincRadius );
      }
    else if( strcmp( "bspline", interp.c_str() ) == 0 )
      {
      unsigned int order = 3;
      if( argc >= 8 )
        {
        order = std::stoi( argv[argct++] );
        }
      // std::cout << "Using bspline interpolation of order " << order << std::endl;

      BSplineInterpolatorPointerType interpB = BSplineInterpolatorType::New();
      interpB->SetSplineOrder( order );

      filter->SetInterpolator( interpB );
      filter->SetIndexPadding( 1 );
      }
    else
      {
      // std::cout << "Using linear interpolation" << std::endl;
      }
    }

  typename InputImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<InputImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  if( itk::Math::FloatAlmostEqual( sliceTiming, 0.0f ) )
    {
    sliceTiming =
      image1->GetSpacing()[ImageDimension - 1] / image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 2];
    // std::cout << "Using default slice timing = " << sliceTiming << std::endl;
    }

  // FIXME - rounding error hack
  if ( sliceTiming * static_cast<float>( image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 2] ) > static_cast<float>( image1->GetSpacing()[ImageDimension - 1] ) ) {
    sliceTiming = sliceTiming - 1.0e-8f;
    // std::cout << "Corrected timing = " << sliceTiming << std::endl;
  }
  else {
    // std::cout << "Slice timing is valid" << std::endl;
    }

  filter->SetSliceTiming( sliceTiming );
  filter->SetInput( image1 );
  filter->DebugOn();

    try
      {
      filter->Update();
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      return EXIT_FAILURE;
      }
  //filter->Update();

  WriteImage<OutputImageType>(filter->GetOutput(), outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int AverageOverDimension(int argc, char *argv[])
{
  typedef float                                     PixelType;
  typedef itk::Image<PixelType, ImageDimension>     ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1> AverageImageType;

  typedef itk::AverageOverDimensionImageFilter<ImageType, AverageImageType>
    ImageFilterType;

  if( argc < 6 )
    {
    // std::cout << "Usage: ImageMath 4 average.nii.gz AverageOverDimension time.nii.gz dimension" << std::endl;
    return EXIT_FAILURE;
    }
  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string  fn1 = std::string(argv[argct++]);
  unsigned int dim = std::stoi( argv[argct++] );

  typename ImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename ImageFilterType::Pointer filter = ImageFilterType::New();
  filter->SetInput( image1 );
  filter->SetAveragingDimension( dim );
  filter->SetDirectionCollapseToSubmatrix();
  filter->Update();

  WriteImage<AverageImageType>( filter->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesRegionSCCA(int argc, char *argv[])
{
  typedef float                                    PixelType;
  typedef itk::Image<PixelType, ImageDimension>    InputImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::MinimumMaximumImageCalculator<LabelImageType> LabelCalculatorType;
  typedef itk::ants::antsSCCANObject<InputImageType, double> SCCANType;

  typedef typename SCCANType::MatrixType MatrixType;
  typedef typename SCCANType::VectorType VectorType;

  if( argc < 6 )
    {
    return EXIT_FAILURE;
    }
  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string labelName = std::string(argv[argct++]);
  std::string timeName  = std::string(argv[argct++]);

  // FIXME - add option for multi input for combined CCA

  typename LabelImageType::Pointer labels = nullptr;
  ReadImage<LabelImageType>( labels, labelName.c_str() );

  typename InputImageType::Pointer time = nullptr;
  ReadImage<InputImageType>( time, timeName.c_str() );

  typename LabelCalculatorType::Pointer calc = LabelCalculatorType::New();
  calc->SetImage( labels );
  calc->ComputeMaximum();
  unsigned int nLabels = calc->GetMaximum();
  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  // unsigned int labelCounts[nLabels];
  unsigned int *labelCounts = new unsigned int[nLabels];
  for( unsigned int i = 0; i < nLabels; i++ )
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;

    labelCounts[i] = 0;
    for( unsigned int v = 0; v < nVoxels; v++ )
      {
      idx[0] = v;
      if( labels->GetPixel(idx) == (i + 1) )
        {
        ++labelCounts[i];
        }
      }
    }

  typename InputImageType::Pointer connmat = InputImageType::New();
  typename InputImageType::RegionType region;
  region.SetSize(0, nLabels);
  region.SetSize(1, nLabels);
  connmat->SetRegions( region );
  connmat->Allocate();
  connmat->FillBuffer(0);

  // Coorelation parameters
  bool         robust = false;
  unsigned int iterct = 20;
  bool         useL1 = false;
  float        gradstep = itk::Math::abs ( useL1 );
  bool         keepPositive = false;
  float        sparsity = 1.0;
  unsigned int minClusterSize = 1;
  unsigned int minRegionSize = 1;

  // used to rankify matrices if using robust
  typename SCCANType::Pointer cca_rankify = SCCANType::New();
  for( unsigned int i = 0; i < nLabels; i++ )
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;

    MatrixType P(nTimes, labelCounts[i], 0.0);

    unsigned int iCount = 0;
    for( unsigned int v = 0; v < nVoxels; v++ )
      {
      idx[0] = v;
      typename InputImageType::IndexType timeIdx;
      timeIdx[1] = v;

      if( labels->GetPixel(idx) == (i + 1) )
        {
        for( unsigned int t = 0; t < nTimes; t++ )
          {
          timeIdx[0] = t;
          P(t, iCount) = time->GetPixel(timeIdx);
          }
        ++iCount;
        }
      }

    if( robust && ( labelCounts[i] >= minRegionSize)  )
      {
      P = cca_rankify->RankifyMatrixColumns(P);
      }

    if( labelCounts[i] >= minRegionSize )
      {
      for( unsigned int j = i + 1; j < nLabels; j++ )
        {
        MatrixType Q(nTimes, labelCounts[j], 0.0);
        typename LabelImageType::IndexType idx2;
        idx2[1] = 0;

        unsigned int jCount = 0;
        for( unsigned int v2 = 0; v2 < nVoxels; v2++ )
          {
          idx2[0] = v2;
          typename InputImageType::IndexType timeIdx2;
          timeIdx2[1] = v2;

          if( labels->GetPixel(idx2) == (j + 1) )
            {
            for( unsigned int t2 = 0; t2 < nTimes; t2++ )
              {
              timeIdx2[0] = t2;
              Q(t2, jCount) = time->GetPixel(timeIdx2);
              }
            ++jCount;
            }
          }

        if( robust )
          {
          Q = cca_rankify->RankifyMatrixColumns(Q);
          }

        if( labelCounts[j] >= minRegionSize )
          {

          // Correlation magic goes here
          typename SCCANType::Pointer cca = SCCANType::New();
          cca->SetSilent( true );
          cca->SetMaximumNumberOfIterations(iterct);
          cca->SetUseL1( useL1 );
          cca->SetGradStep( gradstep );
          cca->SetKeepPositiveP( keepPositive );
          cca->SetKeepPositiveQ( keepPositive );
          cca->SetFractionNonZeroP( sparsity );
          cca->SetFractionNonZeroQ( sparsity );
          cca->SetMinClusterSizeP( minClusterSize );
          cca->SetMinClusterSizeQ( minClusterSize );
          cca->SetMatrixP( P );
          cca->SetMatrixQ( Q );
          VectorType sccancorrs = cca->GetCanonicalCorrelations();

          typename InputImageType::IndexType connIdx;
          connIdx[0] = i;
          connIdx[1] = j;
          connmat->SetPixel( connIdx, sccancorrs[0] );
          connIdx[0] = j;
          connIdx[1] = i;
          connmat->SetPixel( connIdx, sccancorrs[0] );

          }

        }
      }
    }

  WriteImage<InputImageType>(connmat, outname.c_str() );

  delete []labelCounts;

  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesRegionCorr(int argc, char *argv[])
{
  typedef float                                    PixelType;
  typedef itk::Image<PixelType, ImageDimension>    InputImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::MinimumMaximumImageCalculator<LabelImageType> LabelCalculatorType;
  typedef itk::ants::antsSCCANObject<InputImageType, double> SCCANType;

  typedef typename SCCANType::MatrixType MatrixType;
  typedef typename SCCANType::VectorType VectorType;

  if( argc < 6 )
    {
    return EXIT_FAILURE;
    }

  int         argct = 2;
  std::string outname = std::string(argv[argct++]);
  std::string operation = std::string(argv[argct++]);
  std::string labelName = std::string(argv[argct++]);
  std::string timeName  = std::string(argv[argct++]);

  unsigned int minRegionSize = 3;

  if( argc > 6 )
    {
    minRegionSize = std::stoi( argv[argct++] );
    }

  // FIXME - add option for multi input for combined CCA

  typename LabelImageType::Pointer labels = nullptr;
  ReadImage<LabelImageType>( labels, labelName.c_str() );

  typename InputImageType::Pointer time = nullptr;
  ReadImage<InputImageType>( time, timeName.c_str() );

  typename LabelCalculatorType::Pointer calc = LabelCalculatorType::New();
  calc->SetImage( labels );
  calc->ComputeMaximum();
  unsigned int nLabels = calc->GetMaximum();
  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  VectorType labelCounts( nLabels, 0 );

  typename InputImageType::Pointer connmat = InputImageType::New();
  typename InputImageType::RegionType region;
  region.SetSize(0, nLabels);
  region.SetSize(1, nLabels);
  connmat->SetRegions( region );
  connmat->Allocate();
  connmat->FillBuffer(-1);

  MatrixType timeSig( nLabels, nTimes, 0.0 );
  for( unsigned int i = 0; i < nLabels; i++ )
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;
    for( unsigned int v = 0; v < nVoxels; v++ )
      {
      idx[0] = v;
      if( labels->GetPixel(idx) == (i + 1) )
        {
        labelCounts[i]++;

        typename InputImageType::IndexType timeIdx;
        timeIdx[1] = v;
        for( unsigned int t = 0; t < nTimes; t++ )
          {
          timeIdx[0] = t;
          timeSig(i, t) += static_cast<double>( time->GetPixel(timeIdx) );
          }
        }
      }
    }
  for( unsigned int i = 0; i < nLabels; i++ )
    {
    for( unsigned int j = 0; j < nTimes; j++ )
      {
      timeSig(i, j) /= labelCounts[i];
      }
    }
  for( unsigned int i = 0; i < nLabels; i++ )
    {
    for( unsigned int j = (i + 1); j < nLabels; j++ )
      {

      if( (labelCounts[i] > minRegionSize) && (labelCounts[j] > minRegionSize ) )
        {
        VectorType p = timeSig.get_row(i);
        VectorType q = timeSig.get_row(j);

        double corr = 0.0;
        double xysum = 0;
        for( unsigned int z = 0; z < p.size(); z++ )
          {
          xysum += (p[z] * q[z]);
          }

        double frac = 1.0 / (double)p.size();
        double xsum = p.sum();
        double ysum = q.sum();
        double xsqr = p.squared_magnitude();
        double ysqr = q.squared_magnitude();
        double numer = xysum - frac * xsum * ysum;
        double denom = sqrt( ( xsqr - frac * xsum * xsum) * ( ysqr - frac * ysum * ysum) );
        if( denom > 0 )
          {
          corr = numer / denom;
          }
        if( !std::isfinite( corr ) )
          {
          corr = 0.0;
          }

        typename InputImageType::IndexType connIdx;
        connIdx[0] = i;
        connIdx[1] = j;
        connmat->SetPixel( connIdx, corr );
        connIdx[0] = j;
        connIdx[1] = i;
        connmat->SetPixel( connIdx, corr );

        }

      }
    }

  WriteImage<InputImageType>(connmat, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int PASLQuantifyCBF(int argc, char *argv[])
{
  typedef float                                     PixelType;
  typedef itk::Image<PixelType, ImageDimension>     TimeImageType;
  typedef itk::Image<PixelType, ImageDimension - 1> ImageType;

  typedef itk::PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TimeImageType, ImageType, TimeImageType>
    FilterType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string fn1 = std::string(argv[argct++]);
  std::string m0name = std::string(argv[argct++]);

  typename FilterType::Pointer getCBF = FilterType::New();

  // scan for optional parameters
  while( argct < argc )
    {
    if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "TI1") == 0 )
      {
      getCBF->SetTI1( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "TI2") == 0 )
      {
      getCBF->SetTI2( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "T1blood") == 0 )
      {
      getCBF->SetT1blood( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "Lambda") == 0 )
      {
      getCBF->SetLambda( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "Alpha") == 0 )
      {
      getCBF->SetAlpha( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "SliceDelay") == 0 )
      {
      getCBF->SetSliceDelay( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }

    argct++;
    }

  // read in optional parameters here
  // float m_TI1;
  // float m_TI2;
  // float m_T1blood;
  // float m_lmabda;
  // float m_alpha;
  // float m_sliceDelay;

  typename TimeImageType::Pointer diff = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<TimeImageType>(diff, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename ImageType::Pointer m0 = nullptr;
  if( m0name.length() > 3 )
    {
    ReadImage<ImageType>(m0, m0name.c_str() );
    }
  else
    {
    return 1;
    }

  getCBF->SetDifferenceImage( diff );
  getCBF->SetReferenceImage( m0 );
  getCBF->Update();

  WriteImage<TimeImageType>( getCBF->GetOutput(), outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int PCASLQuantifyCBF(int argc, char *  /*NOT USED argv*/[])
{
  if( argc < 6 )
    {
    return EXIT_FAILURE;
    }
  /*
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> TimeImageType;
  typedef itk::Image<PixelType, ImageDimension-1> ImageType;

  typedef itk::PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TimeImageType, ImageType, TimeImageType> FilterType;
  int         argct = 2;
  std::string outname = std::string(argv[argct++]);
  std::string operation = std::string(argv[argct++]);
  std::string fn1 = std::string(argv[argct++]);
  std::string m0name = std::string(argv[argct++]);

  typename FilterType::Pointer getCBF = FilterType::New();

  // scan for optional parameters
  while( argct < argc )
    {
    if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "TI1") == 0 )
      {
      getCBF->SetTI1( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "TI2") == 0 )
      {
      getCBF->SetTI2( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "T1blood") == 0 )
      {
      getCBF->SetT1blood( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "Lambda") == 0 )
      {
      getCBF->SetLambda( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "Alpha") == 0 )
      {
      getCBF->SetAlpha( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }
    else if( strcmp(ANTSOptionName( argv[argct] ).c_str(), "SliceDelay") == 0 )
      {
      getCBF->SetSliceDelay( atof( ANTSOptionValue( argv[argct] ).c_str() ) );
      }

    argct++;
    }


  typename TimeImageType::Pointer diff = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<TimeImageType>(diff, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename ImageType::Pointer m0 = nullptr;
  if( m0name.length() > 3 )
    {
    ReadImage<ImageType>(m0, m0name.c_str() );
    }
  else
    {
    return 1;
    }

  getCBF->SetDifferenceImage( diff );
  getCBF->SetReferenceImage( m0 );
  getCBF->Update();

  WriteImage<TimeImageType>( getCBF->GetOutput(), outname.c_str() );
  */
  return 0;
}

template <unsigned int ImageDimension>
int ComputeTimeSeriesLeverage(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::IndexType                IndexType;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int k_neighbors = std::stoi(argv[argct]);   argct++;
  typename ImageType::Pointer image1 = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }
  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  unsigned int sub_vol = 0;
  extractRegion.SetIndex(ImageDimension - 1, sub_vol );
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( image1 );
  //    extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();
  typename OutImageType::Pointer outimage = extractFilter->GetOutput();
  outimage->FillBuffer(0);

  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

  // step 1.  compute , for each image in the time series, the effect on the average.
  // step 2.  the effect is defined as the influence of that point on the average or, more simply, the distance of that
  // image from the average ....
  typedef vnl_vector<Scalar> timeVectorType;
  timeVectorType mSample(timedims, 0);
  timeVectorType mLeverage(timedims, 0);
  timeVectorType kDistance(timedims, 0);
  SliceIt        vfIter2( outimage, outimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    IndexType    tind;
    // first collect all samples for that location
    for( unsigned int i = 0; i < ImageDimension - 1; i++ )
      {
      tind[i] = ind[i];
      }
    double all_mean = 0;
    for( unsigned int t = 0; t < timedims; t++ )
      {
      tind[ImageDimension - 1] = t;
      Scalar pix = image1->GetPixel(tind);
      mSample(t) = pix;
      all_mean += pix;
      }
    // compute mean time series value at this voxel
    all_mean /= (double)timedims;
    // second compute the leverage for each time point and add that to the total leverage
    //
    // this is a simple approach --- just the difference from the mean.
    //
    for( unsigned int t = 0; t < timedims; t++ )
      {
      mLeverage(t) += fabs(all_mean - mSample(t) ) / (Scalar)timedims;
      }
    }

  // now use k neighbors to get a distance
  kDistance.fill(0);
  for( unsigned int t = 0; t < timedims; t++ )
    {
    int lo = (int)t - k_neighbors / 2;
    int hi = (int)t + k_neighbors / 2;
    if( lo < (int) 0 )
      {
      lo = 0;
      }
    if( hi > (int)(timedims - 1) )
      {
      hi = timedims - 1;
      }
    unsigned int ct = 0;
    for( int k = lo; k < hi; k++ )
      {
      if( k != (int)t )
        {
        kDistance(t) += fabs(mLeverage(t) - mLeverage(k) );
        ct++;
        }
      }
    kDistance(t) /= (double)ct;
    kDistance(t) = kDistance(t) / mLeverage(t);
    }

  // now write the mLeverage value for each time point ...
  std::ofstream logfile;
  logfile.open(outname.c_str() );
  if( logfile.good() )
    {
    // std::cout << "Raw_Leverage,K_Neighbors_Distance" <<  std::endl;
    logfile << "Raw_Leverage,K_Neighbors_Distance" <<  std::endl;
    for( unsigned int t = 0; t < timedims; t++ )
      {
      // std::cout <<  mLeverage(t) << "," << kDistance(t) << std::endl;
      logfile <<  mLeverage(t) << "," << kDistance(t) << std::endl;
      }
    }
  logfile.close();
  return 0;
}

template <unsigned int ImageDimension>
int TimeSeriesToMatrix(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, 2>                     MatrixImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::IndexType                IndexType;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  bool tomha = true;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       ext = itksys::SystemTools::GetFilenameExtension( outname );
  if( ( strcmp(ext.c_str(), ".csv") != 0 ) && (  strcmp(ext.c_str(), ".mha") != 0 )  )
    {
    // std::cout << " must use .csv or .mha as output file extension " << std::endl;
    return EXIT_FAILURE;
    }
  if( ( strcmp(ext.c_str(), ".csv") == 0 )  ) tomha = false;
  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string maskfn = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer mask = nullptr;
  typename MatrixImageType::Pointer matriximage = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }
  if( maskfn.length() > 3 )
    {
    ReadImage<OutImageType>(mask, maskfn.c_str() );
    }
  else
    {
    return 1;
    }
  unsigned int  timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  unsigned long voxct = 0;
  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;
  SliceIt mIter( mask, mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= static_cast<PixelType>( 0.5 ) )
      {
      voxct++;
      }
    }
  // allocate the matrix image
  typename MatrixImageType::SizeType size;
  size[0] = timedims;
  size[1] = voxct;
  typename MatrixImageType::RegionType newregion;
  newregion.SetSize(size);
  if ( tomha ) matriximage = AllocImage<MatrixImageType>(newregion, 0);

  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  unsigned int sub_vol = 0;
  extractRegion.SetIndex(ImageDimension - 1, sub_vol );
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( image1 );
  //    extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();
  typename OutImageType::Pointer outimage = extractFilter->GetOutput();
  outimage->FillBuffer(0);

  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

  typedef vnl_vector<Scalar> timeVectorType;
  timeVectorType mSample(timedims, 0);
  typedef itk::Array2D<double> MatrixType;
  std::vector<std::string> ColumnHeaders;
  MatrixType               matrix;
  if ( ! tomha )
    {
    matrix.set_size(timedims, voxct);
    matrix.Fill(0);
    }
  SliceIt vfIter2( outimage, outimage->GetLargestPossibleRegion() );
  voxct = 0;
  PixelType meanval = 0;
  unsigned long fullct = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( mask->GetPixel(ind) >= static_cast<PixelType>( 0.5 ) )
      {
      IndexType tind;
      // first collect all samples for that location
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
	typename MatrixImageType::IndexType matind;
	matind.Fill(0);
	matind[1] = t;
	matind[2] = voxct;
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mSample(t) = pix;
        if ( ! tomha ) matrix[t][voxct] = pix;
        if ( tomha )
          {
	  matriximage->SetPixel( matind, pix );
	  }
	meanval += static_cast<PixelType>( pix );
	fullct++;
        }
      std::string colname = std::string("V") + ants_to_string<unsigned int>(voxct);
      ColumnHeaders.push_back( colname );
      voxct++;
      } // check mask
    }
  // std::cout << " Mean " << meanval / fullct << std::endl;
  if ( ! tomha )
    {
    // write out the array2D object
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outname );
    writer->SetInput( &matrix );
    writer->SetColumnHeaders( ColumnHeaders );
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& /* exp */ )
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  if ( tomha )
    {
    WriteImage<MatrixImageType>( matriximage , outname.c_str() );
    }
  return 0;
}

template <unsigned int ImageDimension>
int PASL(int argc, char *argv[])
{
  if( argc <= 3 )
    {
    // std::cout << " too few options " << argv[0] << std::endl;
    // std::cout << argv[0] << " NDImage  Bool_FirstImageIsControl optional-M0mask.nii.gz " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef double                                       RealType;
  typedef vnl_vector<RealType>                         timeVectorType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  bool        firstiscontrol = std::stoi(argv[argct]);   argct++;
  std::string m0fn = "";
  if( argc > argct )
    {
    m0fn = std::string(argv[argct]);
    }
  argct++;
  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;
  typename OutImageType::Pointer M0image = nullptr;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;

  // RealType M0W = 1300; // FIXME
  // RealType TE = 4000;
  // RealType calculatedM0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
  RealType calculatedM0 = 2800; // from "Impact of equilibrium magnetization of blood on ASL quantification" by YChen et
                                // al

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  if( firstiscontrol )
    {
    extractRegion.SetIndex(ImageDimension - 1, 0 );
    }
  else
    {
    extractRegion.SetIndex(ImageDimension - 1, 1 );
    }

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( image1 );
  //  extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();
  M0image = extractFilter->GetOutput();

  outimage = AllocImage<OutImageType>(M0image, 0);

  bool haveM0 = true;
  if( m0fn.length() > 3 )
    {
    ReadImage<OutImageType>( M0image, m0fn.c_str() );
    }
  else
    {
    haveM0 = false;
    std::cout
      <<
      "Warning: using calculatedM0 as reference M0 value --- see see 'Impact of equilibrium magnetization of blood on ASL quantification' "
      << std::endl;
    M0image->FillBuffer( calculatedM0 );
    }

  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator vfIter2( outimage,  outimage->GetLargestPossibleRegion() );

  timeVectorType sample( timedims, 0 );
  timeVectorType cbf( timedims, 0 );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    IndexType    tind;
    for( unsigned int i = 0; i < ImageDimension - 1; i++ )
      {
      tind[i] = ind[i];
      }
    RealType      total = 0;
    unsigned long cbfct = 0;
    float      M_0   = static_cast<float>( M0image->GetPixel( ind ) ); // FIXME can be taken from an input reference image or defined for
                                                    // each tissue
    bool getCBF = true;
    if( haveM0 && itk::Math::FloatAlmostEqual( M_0, 0.0f ) )
      {
      getCBF = false;
      }
    else if( haveM0 )
      {
      M_0 = calculatedM0;
      }
    if( getCBF )
      {
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        RealType pix = image1->GetPixel(tind);
        sample( t ) = pix;
        if( ( t % 2 ) == 1 )
          {                                                /**  the best resource i've found so far equation 1 http://cfn.upenn.edu/perfusion/articles/perfmri_9.pdf
            "Pediatric Perfusion Imaging Using Pulsed Arterial Spin Labeling"
      CBF images calculated as:
          f =  \frac{      \lambda DeltaM        }  {     2 alpha M_0 TI_1 exp( - TI_2 / T_{1a} )  }
                TI_2 = TI_1 + t * slice_delay
             where t is the image index , DeltaM is the difference signal between tag and control acquisitions,
             lambda = 0.9 ml/g is the blood/tissue water partition, T_{1a} = 1200 ms is the longitudinal relaxation time of blood,
             alpha = 0.95 is the inversion (or labeling or tagging?) efciency, TI_1 = 800 millisec is the duration
             between the inversion and saturation pulses, TI_2 = TI_1 + w is the image acquisition time.  M_0 is
             the acquired image.  These parameters were primarily based on experience in healthy adults; potential effects
             of ignoring the difference between adults and children on CBF quantication will be discussed below.
             ...
             also see https://gate.nmr.mgh.harvard.edu/wiki/whynhow/images/e/e2/ASL_whyNhow.pdf
   */
          RealType lambda = 0.9;                           //  grams / mL
          RealType alpha = 0.95;                           // labeling efficiency
          RealType deltaM = sample( t - 1 ) - sample( t ); // if 1st image is control
          if( !firstiscontrol )
            {
            deltaM = sample( t ) - sample( t - 1 );                //  2nd image is control
            }
          bool     is1pt5T = false;
          RealType T_1a = 1650; // 3T
          if( is1pt5T )
            {
            T_1a = 1390;            // 1.5T
            }
          // see "Impact of equilibrium magnetization of blood on ASL quantification"
          RealType TI_1  = 600;        // FIXME milliseconds
          RealType slice_delay = 42.0; // FIXME milliseconds
          // TI2(slice) = TI2 + slice_number * slice_delay (slice_delay = the time taken to acquire each slice)
          RealType TI_2  = TI_1 + t * slice_delay;
          RealType scaling = 2 * alpha * static_cast<RealType>( M_0 ) * TI_1 * static_cast<RealType>( exp( -TI_2 / T_1a ) );
          cbf( t ) = lambda * deltaM / scaling;
          total += cbf( t );
          cbfct++;
          }
        }
      }
    RealType mean = total / (RealType) cbfct;
    vfIter2.Set( mean );
    }

  /** From Quantitative Imaging of Perhsion Using a Single Subtraction (QUIPSS and QUIPSS 11)
  In  a proton density weighted, high-resolution, gradient-echo conventional image
  (TE = 5 ms ,  TR  = 1000 ms ,  a  = l o o ) ,  the measured ratio R  of proton density of
  blood in  the saggital sinus to  that  of white matter was 1.06.
  In a single-shot EPI image (TR = \infty ),  the signal M_{0WM} from  white matter was  measured.
  The fully T_1 relaxed signal from blood was then taken to be

  M_OB  = R M_{0WM} exp [ ( 1 / T_{2WM} - 1 / T_{2B} ) TE ] ,
  where T_{2WM} = 80ms  T_{2B} = 200 ms and TE = 5ms.
  */
  //  RealType M_0B = 1.06 * M0wm * exp( 5 / 80 - 5 / 200 );

  /*
RealType M0W = 1; // FIXME
M_0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
The term, M_{0B} is calculated as follows (Wong1998),
    M0b = A * M_{0WM} * exp(1/T2_{WM}^* - 1/T2_B^*) * TE
    where:
    A is the proton density ratio between blood and white matter (assumed to be 1.06)
    T2^* (GRE echo-planar imaging)
      T2_{WM} is 55 msec  (1.5T), 40 (3.0T), and 30 (4.0T)
      T2_B    is 100 msec (1.5T), 80 (3.0T), and 60 (4.0T)
    M_{0WM} is the mean value in an homogenous white matter region from a image acquired with short TE long TR.
  */

  WriteImage<OutImageType>(outimage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int pCASL(int argc, char *argv[])
{
  if( argc <= 3 )
    {
    // std::cout << " too few options " << argv[0] << std::endl;
    // std::cout << argv[0] << " NDImage  Bool_FirstImageIsControl optional-M0mask.nii.gz " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef double                                       RealType;
  typedef vnl_vector<RealType>                         timeVectorType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  bool        firstiscontrol = std::stoi(argv[argct]);   argct++;
  std::string m0fn = "";
  if( argc > argct )
    {
    m0fn = std::string(argv[argct]);
    }
  argct++;
  typename ImageType::Pointer image1;
  typename OutImageType::Pointer outimage;
  typename OutImageType::Pointer M0image;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  typename ImageType::RegionType extractRegion = image1->GetLargestPossibleRegion();
  extractRegion.SetSize(ImageDimension - 1, 0);
  if( firstiscontrol )
    {
    extractRegion.SetIndex(ImageDimension - 1, 0 );
    }
  else
    {
    extractRegion.SetIndex(ImageDimension - 1, 1 );
    }

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput( image1 );
  //  extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->SetDirectionCollapseToSubmatrix();
  extractFilter->SetExtractionRegion( extractRegion );
  extractFilter->Update();
  M0image = extractFilter->GetOutput();

  outimage = AllocImage<OutImageType>(M0image, 0);

  // RealType M0W = 1300; // FIXME
  // RealType TE = 4000;
  // RealType calculatedM0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
  RealType calculatedM0 = 2800; // from "Impact of equilibrium magnetization of blood on ASL quantification" by YChen et
                                // al

  bool haveM0 = true;
  if( m0fn.length() > 3 )
    {
    ReadImage<OutImageType>( M0image, m0fn.c_str() );
    }
  else
    {
    haveM0 = false;
    std::cout
      <<
      "Warning: using calculated value as reference M0 value --- see see 'Impact of equilibrium magnetization of blood on ASL quantification' "
      << std::endl;
    M0image->FillBuffer( calculatedM0 );
    }

  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator vfIter2( outimage,  outimage->GetLargestPossibleRegion() );

  timeVectorType sample( timedims, 0 );
  timeVectorType cbf( timedims, 0 );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    IndexType    tind;
    for( unsigned int i = 0; i < ImageDimension - 1; i++ )
      {
      tind[i] = ind[i];
      }
    RealType      total = 0;
    unsigned long cbfct = 0;
    float      M_0   = static_cast<float>( M0image->GetPixel( ind ) ); // FIXME can be taken from an input reference image or defined for
                                                    // each tissue
    bool getCBF = true;
    if( haveM0 && itk::Math::FloatAlmostEqual( M_0, 0.0f ) )
      {
      getCBF = false;
      }
    else if( haveM0 )
      {
      M_0 = calculatedM0;
      }
    if( getCBF )
      {
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        RealType pix = image1->GetPixel(tind);
        sample( t ) = pix;
        if( ( t % 2 ) == 1 )
          {
          // see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3049525/?tool=pubmed  Quantitative CBF section
          /**
        Quantitative CBF

        Cerebral blood flow in mL per 100g per minute was calculated pixel-by-pixel using equation (1), where the PLD was taken to be longer than ATT, reducing to (Wang et al, 2002),
        CBF =  \frac{      \lambda DeltaM     }  {     2 alpha M_0 T_1a [ exp( - w / T_1a ) - exp( - ( tau + w ) / T_1a )  ] }
        w of 0.7seconds was determined based on the results of experiment (1) for optimal contrast of the GM. Cerebral blood flow was calculated for a single PLD. Quantitative CBF for the whole brain, GM, and WM were tabulated. The bottom slice was excluded from this analysis because it covered only a small part of the cerebrum.
           */
          // f =  \frac{      \lambda DeltaM     }  {     2 alpha M_0 T_1a [ exp( - w / T_1a ) - exp( - ( tau + w ) /
          // T_1a )  ] }
          RealType lambda = 0.9;                           //  grams / mL
          RealType deltaM = sample( t ) - sample( t - 1 ); //  control - tagged if  control is odd
          RealType alpha = 0.85;                           // labeling efficiency
          // or Tagging efficiency: Magic parameter. Reference values: pCASL 0.85, CASL at 3T 0.68, PASL 0.95.
          bool     is1pt5T = false;
          RealType T_1a = 1650; // 3T  from ASL_whyNhow.pdf
          if( is1pt5T )
            {
            T_1a = 1390;            // 1.5T
            }
          RealType T_1t = 1300; // 3T
          if( is1pt5T )
            {
            T_1t = 900;            // 1.5T
            }
          // from "Impact of equilibrium magnetization of blood on ASL quantification"
          RealType tau  = 2100; // FIXME milliseconds from PMC3049525
          RealType w    = 700;  // FIXME milliseconds from PMC3049525
          // Label width: Not in dicom, but sequence-specific -- magic parameter. Reference values: pCASL 1.5, CASL 1.6,
          // PASL 0.7.
          // from PMC3049525
          const RealType scaling = 4 * alpha * static_cast<RealType>( M_0 ) * T_1t
            * ( exp( -1.0 * ( tau + w ) / T_1a ) - exp( -1.0 * w / T_1t )  );
          cbf( t ) = lambda * deltaM * ( -1.0 )  / scaling;
          total += cbf( t );
          cbfct++;
          }
        }
      }
    RealType mean = total / (RealType) cbfct;
    vfIter2.Set( mean );
    }

  // see "Impact of equilibrium magnetization of blood on ASL quantification"
  /*
RealType M0W = 1300; // FIXME
RealType TE = 4000;
M_0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
The term, M_{0B} is calculated as follows (Wong1998),
    M0b = A * M_{0WM} * exp(1/T2_{WM}^* - 1/T2_B^*) * TE
    where:
    A is the proton density ratio between blood and white matter (assumed to be 1.06)
    T2^* (GRE echo-planar imaging)
      T2_{WM} is 55 msec  (1.5T), 40 (3.0T), and 30 (4.0T)
      T2_B    is 100 msec (1.5T), 80 (3.0T), and 60 (4.0T)
    M_{0WM} is the mean value in an homogenous white matter region from a image acquired with short TE long TR.
  */

  WriteImage<OutImageType>(outimage, outname.c_str() );

  return 0;
}

/* Calculate the Magnetization Transfer Ratio from a reference image
* and an MT image. Truncate values to be in range [0,1]
*/

template <unsigned int ImageDimension>
int MTR(int argc, char *argv[])
{
  if( argc <= 5 )
    {
    // std::cout << " too few options " << argv[0] << std::endl;
    // std::cout << argv[0] << " M0Image.nii.gz M1Image.nii.gz [OptionalMask.nii.gz] " << std::endl;
    return 1;
    }

  typedef itk::Image<float, ImageDimension> ImageType;

  typename ImageType::Pointer M0;
  typename ImageType::Pointer M1;
  ReadImage<ImageType>(M0, argv[4] );
  ReadImage<ImageType>(M1, argv[5] );

  typename ImageType::Pointer MTR =
    AllocImage<ImageType>(M0, 0);

  typename ImageType::Pointer mask;
  if( argc > 6 )
    {
    ReadImage<ImageType>(mask, argv[6]);
    }
  else
    {
    mask = AllocImage<ImageType>(M0, 1);
    }

  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIt;
  ImageIt it( mask, mask->GetLargestPossibleRegion() );

  while( !it.IsAtEnd() )
    {
    if( it.Value() > 0 )
      {
      float m0 = M0->GetPixel( it.GetIndex() );
      float m1 = M1->GetPixel( it.GetIndex() );
      float mtr = ( m0 - m1 ) / m0;

      if( mtr < 0 )
        {
        mtr = 0;
        }
      else if( mtr > 1 )
        {
        mtr = 1;
        }
      MTR->SetPixel( it.GetIndex(), mtr );
      }
    ++it;
    }

  WriteImage<ImageType>( MTR, argv[2] );
  return 0;
}

template <unsigned int ImageDimension>
int CompCorrAuto(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::IndexType                IndexType;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  fn_label = std::string(argv[argct]);   argct++;
  unsigned int n_comp_corr_vecs = 4; // number of eigenvectors to get from high variance voxels
  if( argc > argct )
    {
    n_comp_corr_vecs = std::stoi(argv[argct]);
    }
  argct++;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;
  typename OutImageType::Pointer outimage2 = nullptr;
  typename OutImageType::Pointer label_image = nullptr;
  typename OutImageType::Pointer var_image = nullptr;


  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(label_image, fn_label.c_str() );
    }
  else
    {
    return 1;
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(outimage, fn_label.c_str() );
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(outimage2, fn_label.c_str() );
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(var_image, fn_label.c_str() );
    }
  var_image->FillBuffer(0);
  outimage->FillBuffer(0);
  outimage2->FillBuffer(0);
  // std::cout << " read images " << std::endl;
  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  // std::cout << "timedims " << timedims << " size " << image1->GetLargestPossibleRegion().GetSize() << std::endl;

  // first, count the label numbers
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator   vfIter2( label_image,  label_image->GetLargestPossibleRegion() );
  unsigned long ct_nuis = 0;
  unsigned long ct_vox = 0;
  // std::cout << " verify input " << std::endl;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( itk::Math::FloatAlmostEqual( vfIter2.Get(), itk::NumericTraits<PixelType>::OneValue() ) )      // in brain
      {
      ct_vox++;
      }
    }
  // std::cout << " counted " << ct_vox << " voxels " <<  std::endl;
  if( ct_vox == 0 )
    {
    // std::cout << ct_vox << " not enough voxels labeled as gm (or brain) " << std::endl;
    return 1;
    }
  // step 1.  compute , in label 3 ( the nuisance region ), the representative value of the time series over the region.
  //  at the same time, compute the average value in label 2 ( the reference region ).
  // step 2.  factor out the nuisance region from the activation at each voxel in the ROI (nonzero labels).
  // step 3.  compute the correlation of the reference region with every voxel in the roi.
  typedef vnl_matrix<Scalar> timeMatrixType;
  typedef vnl_vector<Scalar> timeVectorType;
  timeMatrixType mSample(timedims, ct_vox, 0);
  unsigned long  nuis_vox = 0;
  timeVectorType sample(timedims, 0);
  //  FIRST -- get high variance (in time) voxels
  float maxvar = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( vfIter2.Get() > 0 )      // in-brain
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      float total = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        sample(t) = pix;
        total += static_cast<float>( pix );
        }
      float mean = total / (float)timedims;
      float var = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        float samplet = static_cast<float>( sample(t) );
        var += ( samplet - mean) * (samplet - mean);
        }
      var /= (float)(timedims);
      var = sqrt(var);
      if( var > maxvar )
        {
        maxvar = var;
        }
      var_image->SetPixel(ind, var);
      }
    }
  // std::cout << " got var " << std::endl;
  // now build the histogram
  unsigned int   histsize = 50;
  float          binsize = maxvar / histsize;
  timeVectorType varhist(histsize, 0);
  float          varhistsum = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    float        var = var_image->GetPixel(ind);
    if( vfIter2.Get() > 0 && var > 0 )      // in-brain
      {
      int bin = (int)( var / binsize ) - 1;
      if( bin < 0 )
        {
        bin = 0;
        }
      varhist[bin] += 1;
      varhistsum += 1;
      }
    }
  varhist = varhist / varhistsum;
  // std::cout << " got var hist " << std::endl;
  float temp = 0;
  float varval_csf = 0;
  for( unsigned int j = 0; j < histsize; j++ )
    {
    temp += static_cast<float>( varhist(j) );
    if( temp >= 0.95f && varval_csf <=  0 )
      {
      varval_csf = (float)j * binsize;
      }
    }

  // std::cout << " maxvar " << maxvar << " varval_csf " << varval_csf << std::endl;
  ct_nuis = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( var_image->GetPixel(ind) > varval_csf  )      // nuisance
      {
      ct_nuis++;
      }
    }
  timeMatrixType mNuisance(timedims, ct_nuis, 0);
  nuis_vox = 0;
  unsigned long brain_vox = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( var_image->GetPixel(ind) > varval_csf  )      // nuisance
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mNuisance(t, nuis_vox) = pix;
        }
      nuis_vox++;
      }
    if( vfIter2.Get() > 0  )      // in brain
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mSample(t, brain_vox) = pix;
        }
      brain_vox++;
      }
    }
  // factor out the nuisance variables by OLS
  timeVectorType vGlobal = matrixOps->AverageColumns(mSample);
  //  typedef itk::Array2D<double> csvMatrixType;
  timeMatrixType           reducedNuisance(timedims, n_comp_corr_vecs + 1);
  std::vector<std::string> ColumnHeaders;
  std::string              colname = std::string("GlobalSignal");
  ColumnHeaders.push_back( colname );
  reducedNuisance.set_column(0, vGlobal);
  if( ct_nuis == 0 )
    {
    n_comp_corr_vecs = 1;
    }
  for( unsigned int i = 0; i < n_comp_corr_vecs; i++ )
    {
    timeVectorType nuisi = matrixOps->GetCovMatEigenvector(mNuisance, i);
    reducedNuisance.set_column(i + 1, nuisi);
    colname = std::string("CompCorrVec") + ants_to_string<unsigned int>(i + 1);
    ColumnHeaders.push_back( colname );
    }

  // write out these nuisance variables
  // write out the array2D object
  typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
  WriterType::Pointer writer = WriterType::New();
  std::string         kname = tempname + std::string("_compcorr.csv");
  writer->SetFileName( kname );
  writer->SetInput( &reducedNuisance );
  writer->SetColumnHeaders( ColumnHeaders );
  try
    {
    writer->Write();
    }
  catch( itk::ExceptionObject& /* exp */ )
    {
    // std::cout << "Exception caught!" << std::endl;
    // std::cout << exp << std::endl;
    return EXIT_FAILURE;
    }

  timeMatrixType RRt = matrixOps->ProjectionMatrix(reducedNuisance);
  mSample = matrixOps->NormalizeMatrix(mSample);
  mSample = mSample - RRt * mSample;
  brain_vox = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( vfIter2.Get() > 0 )
      {
      timeVectorType samp = mSample.get_column(brain_vox);
// correct the original image
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        image1->SetPixel(tind, samp[t]);
        }
      brain_vox++;
      }
    }
  kname = tempname + std::string("_corrected") + extension;
  WriteImage<ImageType>(image1, kname.c_str() );
  kname = tempname + std::string("_variance") + extension;
  WriteImage<OutImageType>(var_image, kname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int ThreeTissueConfounds(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::IndexType                IndexType;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  fn_label = std::string(argv[argct]);   argct++;
  unsigned int wmlabel = 1;
  unsigned int csflabel = 3;
  if( argc > argct )
    {
    csflabel = std::stoi(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    wmlabel = std::stoi(argv[argct]);
    }
  argct++;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = nullptr;
  typename OutImageType::Pointer outimage = nullptr;
  typename OutImageType::Pointer outimage2 = nullptr;
  typename OutImageType::Pointer label_image = nullptr;
  typename OutImageType::Pointer var_image = nullptr;


  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    return 1;
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(label_image, fn_label.c_str() );
    }
  else
    {
    return 1;
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(outimage, fn_label.c_str() );
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(outimage2, fn_label.c_str() );
    }
  if( fn_label.length() > 3 )
    {
    ReadImage<OutImageType>(var_image, fn_label.c_str() );
    }
  var_image->FillBuffer(0);
  outimage->FillBuffer(0);
  outimage2->FillBuffer(0);
  // std::cout << " read images " << std::endl;
  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  // std::cout << "timedims " << timedims << " size " << image1->GetLargestPossibleRegion().GetSize() << std::endl;

  // first, count the label numbers
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator   vfIter2( label_image,  label_image->GetLargestPossibleRegion() );
  unsigned long ct_nuis = 0;
  unsigned long ct_ref = 0;
  unsigned long ct_gm = 0;
  // std::cout << " verify input " << std::endl;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( itk::Math::FloatAlmostEqual( vfIter2.Get(), static_cast<PixelType>( csflabel ) ) )      // nuisance
      {
      ct_nuis++;
      }
    if( itk::Math::FloatAlmostEqual( vfIter2.Get(), static_cast<PixelType>( wmlabel ) ) )      // reference
      {
      ct_ref++;
      }
    if( vfIter2.Get() > 0 )      // gm roi
      {
      ct_gm++;
      }
    }
  // std::cout << " counted " << ct_gm << " gm voxels " << ct_ref << " reference region voxels " << std::endl;
  if( ct_gm == 0 )
    {
    // std::cout << ct_gm << " not enough voxels labeled as gm (or brain) " << ct_gm << std::endl;
    return 1;
    }
  if( ct_ref == 0 )
    {
    // std::cout << ct_ref << " not enough voxels labeled as reference region " << std::endl;
    return 1;
    }
  // step 1.  compute , in label 3 ( the nuisance region ), the representative value of the time series over the region.
  //  at the same time, compute the average value in label 2 ( the reference region ).
  // step 2.  factor out the nuisance region from the activation at each voxel in the ROI (nonzero labels).
  // step 3.  compute the correlation of the reference region with every voxel in the roi.
  typedef vnl_matrix<Scalar> timeMatrixType;
  typedef vnl_vector<Scalar> timeVectorType;
  timeMatrixType mReference(timedims, ct_ref, 0);
  timeMatrixType mSample(timedims, ct_gm, 0);
  unsigned long  nuis_vox = 0;
  unsigned long  ref_vox = 0;
  unsigned long  gm_vox = 0;
  timeVectorType smoother(timedims, 0);
  timeVectorType smoother_out(timedims, 0);
  //  FIRST -- get high variance (in time) voxels
  float maxvar = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( vfIter2.Get() > 0 )      // in-brain
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      float total = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        float pix = static_cast<float>( image1->GetPixel(tind) );
        smoother(t) = pix;
        total += pix;
        }
      float mean = total / (float)timedims;
      float var = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        float smoothert = static_cast<float>( smoother(t) );
        var += ( (smoothert - mean) * (smoothert - mean) );
        }
      var /= (float)(timedims);
      var = sqrt(var);
      if( var > maxvar )
        {
        maxvar = var;
        }
      var_image->SetPixel(ind, var);
      }
    }
  // std::cout << " got var " << std::endl;
  // now build the histogram
  unsigned int   histsize = 50;
  float          binsize = maxvar / histsize;
  timeVectorType varhist(histsize, 0);
  float          varhistsum = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    float        var = var_image->GetPixel(ind);
    if( vfIter2.Get() > 0 && var > 0 )      // in-brain
      {
      int bin = (int)( var / binsize ) - 1;
      if( bin < 0 )
        {
        bin = 0;
        }
      varhist[bin] += 1;
      varhistsum += 1;
      }
    }
  varhist = varhist / varhistsum;
  // std::cout << " got var hist " << std::endl;
  float temp = 0;
  float varval_csf = 0;
  for( unsigned int j = 0; j < histsize; j++ )
    {
    temp += static_cast<float>( varhist(j) );
    if( temp >= 0.95f && varval_csf <=  0 )
      {
      varval_csf = (float)j * binsize;
      }
    }

  // std::cout << " maxvar " << maxvar << " varval_csf " << varval_csf << std::endl;
  //  WriteImage<OutImageType>(var_image,"varimage.nii.gz");
  //
  ct_nuis = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    // OutIndexType ind = vfIter2.GetIndex();
    if ( itk::Math::FloatAlmostEqual( vfIter2.Get(), static_cast<PixelType>( csflabel ) ) )// reference
    //    if( var_image->GetPixel(ind) > varval_csf  )      // nuisance
      {
      ct_nuis++;
      }
    }
  timeMatrixType mNuisance(timedims, ct_nuis, 0);
  ref_vox = 0; nuis_vox = 0; gm_vox = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    //      if ( vfIter2.Get() == 3 ) { // nuisance
    //    if( var_image->GetPixel(ind) > varval_csf  )      // nuisance
    if( itk::Math::FloatAlmostEqual( vfIter2.Get(), static_cast<PixelType>( csflabel ) ) )// reference
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mNuisance(t, nuis_vox) = pix;
        }
      nuis_vox++;
      }
    if( itk::Math::FloatAlmostEqual( vfIter2.Get(), static_cast<PixelType>( wmlabel ) ) )// reference
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mReference(t, ref_vox) = pix;
        }
      ref_vox++;
      }
    if( vfIter2.Get() > 0  )      // in brain
      {
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mSample(t, gm_vox) = pix;
        }
      gm_vox++;
      }
    }

  // factor out the nuisance variables by OLS
  unsigned int nnuis = 3; // global , csf , wm
  if( ct_nuis == 0 )
    {
    nnuis = 1;
    }
  timeMatrixType reducedNuisance(timedims, nnuis);
  timeVectorType vGlobal = matrixOps->AverageColumns( mSample );
  reducedNuisance.set_column( 0, vGlobal);
  vGlobal = matrixOps->AverageColumns( mNuisance ); // csf
  reducedNuisance.set_column( 1, vGlobal);
  vGlobal = matrixOps->AverageColumns( mReference ); // wm
  reducedNuisance.set_column( 2, vGlobal);

  std::vector<std::string> ColumnHeaders;
  std::string              colname = std::string("GlobalSignal");
  ColumnHeaders.push_back( colname );
  colname = std::string("CSF");
  ColumnHeaders.push_back( colname );
  colname = std::string("WM");
  ColumnHeaders.push_back( colname );

  // write out these nuisance variables
  // write out the array2D object
  typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
  WriterType::Pointer writer = WriterType::New();
  std::string         kname = tempname + std::string("_compcorr.csv");
  writer->SetFileName( kname );
  writer->SetInput( &reducedNuisance );
  writer->SetColumnHeaders( ColumnHeaders );
  try
    {
    writer->Write();
    }
  catch( itk::ExceptionObject& /* exp */ )
    {
    // std::cout << "Exception caught!" << std::endl;
    // std::cout << exp << std::endl;
    return EXIT_FAILURE;
    }

  return 0;

  timeMatrixType RRt = matrixOps->ProjectionMatrix(reducedNuisance);
  mReference = matrixOps->NormalizeMatrix(mReference);
  mReference = mReference - RRt * mReference;
  mSample = matrixOps->NormalizeMatrix(mSample);
  mSample = mSample - RRt * mSample;
  // reduce your reference region to the first & second eigenvector
  timeVectorType vReference = matrixOps->GetCovMatEigenvector(mReference, 0);
  timeVectorType vReference2 = matrixOps->AverageColumns(mReference);
  Scalar         testcorr = matrixOps->PearsonCorr(vReference, vReference2);
  if( testcorr < 0 )
    {
    vReference = vReference * (-1);
    }
  if( vReference.size() != timedims )
    {
    // std::cout << " CompCorr Error exiting " << std::endl; throw std::exception();
    }
  gm_vox = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( vfIter2.Get() > 0 )
      {
      timeVectorType samp = mSample.get_column(gm_vox);
// correct the original image
      IndexType tind;
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        image1->SetPixel(tind, samp[t]);
        }
// compute the gm-reference correlation
      Scalar corr = matrixOps->PearsonCorr(samp, vReference);
      //    Scalar corr2=matrixOps->PearsonCorr(samp,vReference2);
      outimage->SetPixel(ind, corr);
      //    outimage2->SetPixel(ind,corr2);
      outimage2->SetPixel(ind, samp.two_norm() ); // the power of the time series
      gm_vox++;
      }
    }
  // std::cout << "write results" << std::endl;
  kname = tempname + std::string("first_evec") + extension;
  WriteImage<OutImageType>(outimage, kname.c_str() );
  //  kname=tempname+std::string("second_evec")+extension;
  kname = tempname + std::string("power") + extension;
  WriteImage<OutImageType>(outimage2, kname.c_str() );
  kname = tempname + std::string("_corrected") + extension;
  WriteImage<ImageType>(image1, kname.c_str() );
  kname = tempname + std::string("_variance") + extension;
  WriteImage<OutImageType>(var_image, kname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int StackImage(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                          PixelType;
  typedef itk::Image<PixelType, ImageDimension>          ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>   Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;

  std::string fn1 = std::string(argv[argct]);

  //unsigned int nImages = argc - argct;
  //// std::cout << "Stacking " << nImages << " images" << std::endl;

  // Reference image is the first image passed
  // All other input images must match in all dimensions except the stack dimension
  typename ImageType::Pointer refImage = nullptr;

  // Re-used image pointer for all images to be stacked
  typename ImageType::Pointer image1 = nullptr;

  ReadImage<ImageType>(refImage, fn1.c_str());

  typename ImageType::SpacingType refSpacing = refImage->GetSpacing();
  typename ImageType::DirectionType refDirection = refImage->GetDirection();
  typename ImageType::PointType refOrigin = refImage->GetOrigin();
  typename ImageType::SizeType refSize = refImage->GetLargestPossibleRegion().GetSize();

  unsigned int nDims = refImage->GetImageDimension();


  if ( nDims != ImageDimension )
  {
    // std::cout << "Image dimensions not consistent with passed parameters" << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int stackLength = refSize[nDims-1];

  // Check headers align to within this tolerance. Set eps relatively large, aims to catch gross errors
  // without getting hung up on floating point precision issues, which can be sizeable when dealing with
  // nifti I/O
  float eps = 1E-4;

  for ( int i=(argct+1); i<argc; i++)
  {

    std::string fn2 = std::string(argv[i]);

    ReadImage<ImageType>(image1, fn2.c_str());

    typename ImageType::SpacingType im1Spacing = image1->GetSpacing();
    typename ImageType::DirectionType im1Direction = image1->GetDirection();
    typename ImageType::PointType im1Origin = image1->GetOrigin();
    typename ImageType::SizeType im1Size = image1->GetLargestPossibleRegion().GetSize();


    if ( image1->GetImageDimension() != nDims )
    {
      // std::cout << "Inconsistent image dimension in " << argv[i] << std::endl;
      return EXIT_FAILURE;
    }

    // Check input images
    for ( unsigned int d=0; d<nDims; d++)
    {
      if ( static_cast<float>( std::fabs( im1Spacing[d] - refSpacing[d] ) ) > eps )
      {
        // std::cout << "Inconsistent image spacing not allowed" << std::endl;
        return EXIT_FAILURE;
      }

      for ( unsigned int d2=0; d2<nDims; d2++)
      {
        if ( static_cast<float>( std::fabs(im1Direction(d,d2) - refDirection(d,d2)) ) > eps )
        {
          // std::cout << "Inconsistent image direction not allowed" << std::endl;
          return EXIT_FAILURE;
        }
      }

      // Only check origin and size up to nDims - 1
      // Thus we allow stacked images to have different origins and dimension along
      // the stack dimension. For example, allow stacking of two time series with
      // 100 volumes and 50 volumes respectively
      if ( d < (nDims-1) )
      {
        if ( static_cast<float>( std::fabs(im1Origin[d] - refOrigin[d] ) ) > eps )
        {
          // std::cout << "Inconsistent image origins not allowed" << std::endl;
          return EXIT_FAILURE;
        }

        if ( im1Size[d] != refSize[d] )
        {
          // std::cout << "Size variation in stacking dimension only" << std::endl;
          return EXIT_FAILURE;
        }
      }


    }

    stackLength += im1Size[nDims-1];
  }

  typename ImageType::SizeType stackSize  = refSize;
  stackSize[nDims-1] = stackLength;

  typename ImageType::RegionType region = refImage->GetLargestPossibleRegion();
  region.SetSize( stackSize );
  typename ImageType::Pointer stackImage = ImageType::New();
  stackImage->SetRegions( region );
  stackImage->SetDirection( refDirection );
  stackImage->SetOrigin( refOrigin );
  stackImage->SetSpacing( refSpacing );
  stackImage->Allocate();

  unsigned int offset = 0;
  while( argc > argct )
    {

    std::string fn2 = std::string(argv[argct++]);
    ReadImage<ImageType>(image1, fn2.c_str());

    Iterator iter( image1,  image1->GetLargestPossibleRegion() );
    for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      typename ImageType::IndexType oIndex = iter.GetIndex();
      oIndex[nDims-1] = oIndex[nDims-1] + offset;
      stackImage->SetPixel(oIndex, iter.Value() );
      }

    offset += image1->GetLargestPossibleRegion().GetSize()[nDims-1];
    }

  WriteImage<ImageType>(stackImage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int Stack2Images(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);  argct++;
  std::string fn2 = std::string(argv[argct]);  argct++;
  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer tiledimage = nullptr;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }

  itk::FixedArray< unsigned int, ImageDimension > layout;
  for ( unsigned int i = 0; i < (ImageDimension-1); i++ )
    {
    layout[i]=1;
    }
  layout[ ImageDimension - 1 ] = 0;
  typedef itk::TileImageFilter <ImageType, ImageType >
    TileImageFilterType;
  typename TileImageFilterType::Pointer tileFilter
    = TileImageFilterType::New ();
  tileFilter->SetLayout( layout );
  unsigned int inputImageNumber = 0;
  tileFilter->SetInput( inputImageNumber++, image1 );
  tileFilter->SetInput( inputImageNumber++, image2 );
  tileFilter->SetDefaultPixelValue( 0 );
  tiledimage = tileFilter->GetOutput();
  WriteImage<ImageType>(tiledimage, outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int MakeImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  unsigned int sizevalx = std::stoi(argv[argct]);   argct++;
  unsigned int sizevaly = 0;
  if( argc > argct )
    {
    sizevaly = std::stoi(argv[argct]); argct++;
    }
  unsigned int sizevalz = 0;
  if( argc > argct && ImageDimension > 2 )
    {
    sizevalz = std::stoi(argv[argct]); argct++;
    }

  unsigned int sizevalt = 0;
  if( argc > argct && ImageDimension > 3 )
    {
    sizevalt = std::stoi(argv[argct]); argct++;
    }

  typename ImageType::SizeType size;
  size[0] = sizevalx;
  size[1] = sizevaly;
  if( ImageDimension > 2 )
    {
    size[2] = sizevalz;
    }
  if( ImageDimension > 3 )
    {
    size[3] = sizevalt;
    }
  typename ImageType::RegionType newregion;
  newregion.SetSize(size);

  typename ImageType::Pointer padimage = AllocImage<ImageType>(newregion, 0);
  WriteImage<ImageType>(padimage, outname.c_str() );

  return 0;
}

template <typename TImage>
typename TImage::Pointer
LabelSurface(typename TImage::Pointer input, typename TImage::Pointer input2  )
{
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typename   ImageType::Pointer Image =
    AllocImage<ImageType>(input);

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );
  GHood.GoToBegin();
  while( !GHood.IsAtEnd() )
    {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    if( p >= 0.5f )
      {
      bool atedge = false;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        const typename TImage::IndexType & ind2 = GHood.GetIndex(i);
        float dist = 0.0;
        for( int j = 0; j < ImageDimension; j++ )
          {
          dist += (float)(ind[j] - ind2[j]) * (float)(ind[j] - ind2[j]);
          }
        dist = sqrt(dist);
        bool secondval = true;
        if( input2 )
          {
          if( input2->GetPixel(ind2) >= 0.5f )
            {
            secondval = true;
            }
          }
        if( GHood.GetPixel(i) < 0.5f && dist <  2.0f && secondval  )
          {
          atedge = true;
          }
        }
      if( atedge && p >=  0.5f )
        {
        Image->SetPixel(ind, 1);
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


template <unsigned int ImageDimension>
int LabelSurfaceArea(int argc, char *argv[])
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;

  typedef double                                            Scalar;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer input = nullptr;
  ReadImage<ImageType>( input, fn1.c_str() );
  typename   ImageType::Pointer areaImage = AllocImage<ImageType>( input );
  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  rad.Fill( 1 );
  if ( argc > argct )  rad.Fill( std::stoi( argv[argct] ) );
  Scalar voxspc = 0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    voxspc += input->GetSpacing()[i];
  voxspc = voxspc / static_cast<Scalar>( ImageDimension );
  Scalar voxspc2 = voxspc * voxspc;
  Scalar dm1 = static_cast<Scalar>( ImageDimension - 1 );
  Scalar refarea = std::pow( static_cast<Scalar>( rad[1] ) , dm1 );
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );
  GHood.GoToBegin();
  while( !GHood.IsAtEnd() )
    {
    typename ImageType::IndexType ind = GHood.GetIndex();
    Scalar area = 0.0;
    Scalar locrefarea = refarea;
    if (  GHood.GetCenterPixel() > 0.1f )
    {
    for( unsigned int i = 0; i < GHood.Size(); i++ )
      {
      const typename ImageType::IndexType & ind2 = GHood.GetIndex(i);
      Scalar dist = 0.0;
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        dist += static_cast<Scalar>( ( ind[j] - ind2[j] ) * (ind[j] - ind2[j]) );
        }
      dist = sqrt( dist );
      if( dist <  (2.0*voxspc) )
        {
        area += ( ( static_cast<Scalar>( GHood.GetPixel( i ) ) * voxspc2 ) / locrefarea );
        }
      }
    }
    areaImage->SetPixel( ind, area );
    ++GHood;
    }
  WriteImage<ImageType>( areaImage, outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int FitSphere(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    std::cout << " too few options " << std::string(argv[1]) << std::endl;
    return 1;
    }
  /*
  typedef float  PixelType;
  typedef itk::Vector<float,ImageDimension>         VectorType;
  typedef itk::Image<VectorType,ImageDimension>     FieldType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> readertype;
  typedef itk::ImageFileWriter<ImageType> writertype;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int argct=2;
  std::string outname=std::string(argv[argct]); argct++;
  std::string operation = std::string(argv[argct]);  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if (argc > argct) fn2=std::string(argv[argct]);   argct++;
  float MaxRad=5;
  if (argc > argct) MaxRad = atof(argv[argct]);   argct++;

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer radimage = nullptr;
  typename ImageType::Pointer radimage2 = nullptr;
  typename ImageType::Pointer priorimage = nullptr;
  typename ImageType::Pointer wmimage = nullptr;
  if (fn2.length() > 3)   ReadImage<ImageType>(wmimage, fn2.c_str());
  // std::cout <<"  read " << fn1 << " MXR " << MaxRad << std::endl;
  ReadImage<ImageType>(image1, fn1.c_str());
  ReadImage<ImageType>(radimage, fn1.c_str());
  ReadImage<ImageType>(radimage2, fn1.c_str());
  ReadImage<ImageType>(priorimage, fn1.c_str());
  radimage->FillBuffer(0);
  radimage2->FillBuffer(0);
  priorimage->FillBuffer(0);
  typename ImageType::SpacingType spacing=image1->GetSpacing();

  typename ImageType::Pointer surf = LabelSurface<ImageType>(image1,wmimage);

  typedef itk::LinearInterpolateImageFunction<ImageType,float> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer ginterp =  ScalarInterpolatorType::New();
  ginterp->SetInputImage(image1);
  typename ScalarInterpolatorType::ContinuousIndexType  Y1;
  typename ScalarInterpolatorType::ContinuousIndexType  Y2;
  typename ScalarInterpolatorType::ContinuousIndexType  GMx;
  typename ScalarInterpolatorType::ContinuousIndexType  WMx;
  typename ScalarInterpolatorType::Pointer winterp= nullptr;
  if (wmimage)
    {
      winterp=ScalarInterpolatorType::New();
      winterp->SetInputImage(wmimage);
    }
  //  float x=0,y=0,z=0;
  //  float xc=0,yc=0,zc=0;
  float globalbestrad=0;
  typename ImageType::IndexType bestind;
  bestind.Fill(0);
  typename ImageType::IndexType ind2;
  ind2.Fill(0);
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  // std::cout <<"  Begin " << std::endl;
  unsigned long npx=0;
//  float pi=3.141;
 unsigned long numpx=image1->GetBufferedRegion().GetNumberOfPixels();
  //unsigned int prog=0;

//  float gmtotal,wmtotal,gvol,wvol,warea,garea;
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
      npx++;
      typename ImageType::IndexType ind=iter.GetIndex();
//      float val=surf->GetPixel(ind);//iter.Get();
//      float minrad=1.e9;
*/
  /*
  if (val > 0.5 && fabs((float)ind[2]-80.) < 6 )
{
  //float  bestrad=0;
  //parameterize sphere at this index
  float epspi=pi*0.1;
  float epsrad=0.2;
  gvol=0;wvol=0;warea=0;garea=0;
  float svol=0;//4./3.*pi*MaxRad*MaxRad*MaxRad;
  //      float sarea=0;//4.*pi*MaxRad*MaxRad;
  //float bestvol=0;
  //      float wdiff=0;
  //      for (float theta=0; theta<=pi; theta+=epspi)
  float theta=0;
  gmtotal=0;
  wmtotal=0;
  float glength=0,wlength=0;
  while ( theta < pi )
    {
  garea=0;warea=0;
  float psi=0;
  while ( psi < pi )
    {
      glength=0;wlength=0;
      float rr=0;
      bool raddone=false;
      GMx.Fill(0);  WMx.Fill(0);
      while ( rr <= MaxRad && !raddone)
    {
      // integrate the wm/gm probability along this radius, at these angles
      Y1[0]=(float)ind[0]/spacing[0]+rr*cos(psi)*sin(theta);
      Y1[1]=(float)ind[1]/spacing[1]+rr*sin(psi)*sin(theta);
      Y1[2]=(float)ind[2]/spacing[2]+rr*cos(theta);
      Y2[0]=(float)ind[0]/spacing[0]+rr*cos(psi+pi)*sin(theta);
      Y2[1]=(float)ind[1]/spacing[1]+rr*sin(psi+pi)*sin(theta);
      Y2[2]=(float)ind[2]/spacing[2]+rr*cos(theta);
      float gval1 = ginterp->EvaluateAtContinuousIndex( Y1 );
      float gval2 = ginterp->EvaluateAtContinuousIndex( Y2 );
      float wval1=0,wval2=0;
      if (wmimage) wval1 = winterp->EvaluateAtContinuousIndex( Y1 );
      if (wmimage) wval2 = winterp->EvaluateAtContinuousIndex( Y2 );
      glength+=(gval1+gval2)*epsrad*2;
      wlength+=(wval2+wval2)*epsrad*2;
      gmtotal+=(gval1+gval2);
      wmtotal+=(wval1+wval2);
      for (unsigned int dd=0; dd<ImageDimension; dd++)
        {
          GMx[dd]+=(Y1[dd]*gval1+Y2[dd]*gval2);
          WMx[dd]+=(Y1[dd]*wval1+Y2[dd]*wval2);
        }
      rr+=epsrad;
    }//update rr

      garea+=glength*epspi/pi;
      warea+=wlength*epspi/pi;
      psi+=epspi;
    }//update psi
  gvol+=garea;//epspi/pi;
  wvol+=warea;//epspi/pi;
  svol+=(gvol+wvol);
  theta+=epspi;
    }//update theta

      for (unsigned int dd=0; dd<ImageDimension; dd++)
        {
          GMx[dd]/=gmtotal;
          WMx[dd]/=wmtotal;
        }
      float cmdist=0;
      for (unsigned int dd=0; dd<ImageDimension; dd++)
        {
          cmdist+=(GMx[dd]-WMx[dd])*(GMx[dd]-WMx[dd]);
        }
      cmdist=sqrt(cmdist);
      //          // std::cout << " GMT " << gmtotal << " WMT " << wmtotal << " dist " << cmdist << std::endl;
  float gmrad=std::pow( 3.*gvol/(4.*pi) , 1./3.);
  float gwrat=0,gvrat=0;
  if (warea > 0) gwrat=garea/warea;
  if (wvol > 0) gvrat=gvol/wvol;
  priorimage->SetPixel(ind,gmtotal/(wmtotal+gmtotal));
  radimage->SetPixel(ind,cmdist);

}
*/
/*
      //      radimage2->SetPixel(ind,gvrat);
      if (image1->GetPixel(ind) >= 0.5)
    {
      bool okfit=true;
      float dorad=1;
      float bestrad=1;
      while (okfit && dorad <= MaxRad )
        {
          typedef itk::NeighborhoodIterator<ImageType>  iteratorType;
          typename iteratorType::RadiusType rad,rad2;
          rad2.Fill(0);
          for (unsigned int j=0; j<ImageDimension; j++) rad[j]=(long unsigned int) dorad;
          float tardist=(dorad*sqrt((double)2));
          iteratorType GHood(rad, image1,image1->GetLargestPossibleRegion());
          GHood.SetLocation(ind);
          typename ImageType::PixelType p = GHood.GetCenterPixel();
          unsigned int goodct=0;
          unsigned int possct=0;
          if ( p > 0 && radimage->GetPixel(ind) <= 0  )
        {
          for (unsigned int i = 0; i < GHood.Size(); i++)
            {
              ind2=GHood.GetIndex(i);
              float dist=sqrt(((float)ind2[0]-(float)ind[0])*((float)ind2[0]-(float)ind[0])+
                      ((float)ind2[1]-(float)ind[1])*((float)ind2[1]-(float)ind[1])+
                      ((float)ind2[2]-(float)ind[2])*((float)ind2[2]-(float)ind[2]));
              if ( GHood.GetPixel(i) == p && dist <= tardist)
            {
              goodct++;
              possct++;
            }
              else if ( dist <= tardist ) possct++;
              //          // std::cout << " Ind " <<  ind << " : " <<  bestrad << " tardist " << tardist << " gct " << goodct <<" pos " << possct << " dist " << dist << " ind2 " << ind2 << std::endl;
            }
          if (goodct==possct)
            {
              bestrad=dorad; radimage->SetPixel(ind,bestrad*(-1.0));
            }
          else { okfit=false; radimage->SetPixel(ind,radimage->GetPixel(ind)*(-1.0)); }
        }
          dorad=dorad+1;
        }
      if (bestrad >= globalbestrad)
        {
          globalbestrad=bestrad;
          bestind=ind;
        }

      if (npx % 10000 == 0)
        {
          // std::cout <<" prog " << (float)npx/(float)numpx << std::endl;
          //          WriteImage<ImageType>(radimage,outname.c_str());
          //          WriteImage<ImageType>(radimage2,(std::string("Sphere")+outname).c_str());
          //WriteImage<ImageType>(priorimage,(std::string("Prior")+outname).c_str());
        }
    }

    }


  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
      float val=iter.Get();
      typename ImageType::IndexType ind=iter.GetIndex();
      if (val > 0)
    {
      unsigned int dorad=(unsigned int)fabs(radimage->GetPixel(ind));
      typedef itk::NeighborhoodIterator<ImageType>  iteratorType;
      typename iteratorType::RadiusType rad;
      for (unsigned int j=0; j<ImageDimension; j++) rad[j]= dorad;
      float tardist=(dorad*sqrt((double)2));
      float diameter=2.0*(float)dorad*sqrt((double)2);
      iteratorType GHood(rad, image1,image1->GetLargestPossibleRegion());
      GHood.SetLocation(ind);
      typename ImageType::PixelType p = GHood.GetCenterPixel();
      for (unsigned int i = 0; i < GHood.Size(); i++)
        {
          ind2=GHood.GetIndex(i);
          float dist=sqrt(((float)ind2[0]-(float)ind[0])*((float)ind2[0]-(float)ind[0])+
                  ((float)ind2[1]-(float)ind[1])*((float)ind2[1]-(float)ind[1])+
                  ((float)ind2[2]-(float)ind[2])*((float)ind2[2]-(float)ind[2]));
          if ( GHood.GetPixel(i) == p && dist <= tardist && diameter > priorimage->GetPixel(ind2))
        {
          priorimage->SetPixel(ind2,diameter);
        }
        }
    }
    }




  // now, make rad image
  // std::cout << " Best " << bestind << " gbr " << globalbestrad << std::endl;
  typedef itk::NeighborhoodIterator<ImageType>  iteratorType;
  typename iteratorType::RadiusType rad;
  for (unsigned int j=0; j<ImageDimension; j++) rad[j]= (long unsigned int)globalbestrad;
  iteratorType GHood(rad, image1,image1->GetLargestPossibleRegion());
  GHood.SetLocation(bestind);
  typename ImageType::PixelType p = GHood.GetCenterPixel();
  for (unsigned int i = 0; i < GHood.Size(); i++)
    {
      ind2=GHood.GetIndex(i);
      float dist=0;
      dist+=sqrt((float)(ind2[0]-(float)bestind[0])*(float)(ind2[0]-(float)bestind[0])+
         (float)(ind2[1]-(float)bestind[1])*(float)(ind2[1]-(float)bestind[1])+
         (float)(ind2[2]-(float)bestind[2])*(float)(ind2[2]-(float)bestind[2]));
      if ( dist <= (globalbestrad*sqrt((double)2)))
    {
      radimage2->SetPixel(ind2,p);
    }
    }


  //  WriteImage<ImageType>(radimage,outname.c_str());
  WriteImage<ImageType>(radimage2,outname.c_str());
  //WriteImage<ImageType>(priorimage,(std::string("Prior")+outname).c_str());
 */
  return 0;
}

template <unsigned int ImageDimension>
int ImageMath(int argc, char *argv[])
{
  typedef float PixelType;
  //  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       operation = std::string(argv[argct]);  argct++;
  std::string       fn1 = std::string(argv[argct]);   argct++;
  std::string       fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer varimage = nullptr;

  bool  isfloat = false;
  float floatval = 1.0;
  if( from_string<float>(floatval, fn2, std::dec) )
    {
    isfloat = true;
    }
  else
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }
  ReadImage<ImageType>(image1, fn1.c_str() );
  varimage = AllocImage<ImageType>(image1);

  if( strcmp(operation.c_str(), "mresample") == 0 && !isfloat )
    {
    typename ImageType::SpacingType spc = image2->GetSpacing();
    typedef itk::TranslationTransform<double, ImageDimension> TransformType0;
    typename TransformType0::Pointer             m_Transform0 = TransformType0::New();
    typename TransformType0::ParametersType trans = m_Transform0->GetParameters();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      trans[i] = 0;
      spc[i] = image1->GetSpacing()[i] * image1->GetLargestPossibleRegion().GetSize()[i]
        / image2->GetLargestPossibleRegion().GetSize()[i];
      }
    image2->SetSpacing(spc);
    image2->SetOrigin(image1->GetOrigin() );
    image2->SetDirection(image1->GetDirection() );
    m_Transform0->SetParameters(trans);
    // std::cout << " trans " << m_Transform0->GetParameters() << " Nspc " << image2->GetSpacing() << std::endl;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform( m_Transform0 );
    resample->SetInput( image2 );
    resample->SetOutputParametersFromImage(  image1 );
    typename ImageType::IndexType zeroind;  zeroind.Fill(0);
    resample->SetDefaultPixelValue( image1->GetPixel(zeroind) );
    resample->UpdateLargestPossibleRegion();
    image2 = resample->GetOutput();

    WriteImage<ImageType>(image2, outname.c_str() );
    return 0;
    }

  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < ImageDimension; i++ )
    {
    volumeelement *= static_cast<float>( varimage->GetSpacing()[i] );
    }

  float         result = 0;
  unsigned long ct = 0;
  Iterator      vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    IndexType ind = vfIter2.GetIndex();
    float     pix2;
    if( isfloat )
      {
      pix2 = floatval;
      }
    else
      {
      pix2 = image2->GetPixel(ind);
      }
    float pix1 = image1->GetPixel(ind);

    if( strcmp(operation.c_str(), "m") == 0 )
      {
      result = pix1 * pix2;
      }
    else if( strcmp(operation.c_str(), "+") == 0 )
      {
      result = pix1 + pix2;
      }
    else if( strcmp(operation.c_str(), "-") == 0 )
      {
      result = pix1 - pix2;
      }
    else if( strcmp(operation.c_str(), "/") == 0 )
      {
      if( pix2 > 0 )
        {
        result = pix1 / pix2;
        }
      }
    else if( strcmp(operation.c_str(), "^") == 0 )
      {
      result = std::pow(pix1, pix2);
      }
    else if( strcmp(operation.c_str(), "exp") == 0 )
      {
      result = exp(pix1 * pix2);
      }
    else if( strcmp(operation.c_str(), "max") == 0 )
      {
      result = std::max( pix1, pix2 );
      }
    else if( strcmp(operation.c_str(), "abs") == 0 )
      {
      result = static_cast<float>( fabs(pix1) );
      }
    else if( strcmp(operation.c_str(), "addtozero") == 0 && itk::Math::FloatAlmostEqual( pix1, 0.0f ) )
      {
      result = pix1 + pix2;
      }
    else if( strcmp(operation.c_str(), "addtozero") == 0 && ! itk::Math::FloatAlmostEqual( pix1, 0.0f ) )
      {
      result = pix1;
      }
    else if( strcmp(operation.c_str(), "overadd") == 0 && ! itk::Math::FloatAlmostEqual( pix2, 0.0f ) )
      {
      result = pix2;
      }
    else if( strcmp(operation.c_str(), "overadd") == 0 )
      {
      result = pix1;
      }
    else if( strcmp(operation.c_str(), "Decision") == 0 )
      {
      result = 1.0f / (1.0f + static_cast<float>( exp(-1.0f * ( pix1 - 0.25f) / pix2) ) );
      }
    else if( strcmp(operation.c_str(), "total") == 0 )
      {
      result += pix1 * pix2;
      }
    else if( strcmp(operation.c_str(), "mean") == 0 )
      {
      result += pix1 * pix2;
      ct++;
      }

    vfIter2.Set(result);
    }
  if( strcmp(operation.c_str(), "total") == 0 )
    {
    std::cout << "total: " << result << " total-volume: " << result * volumeelement << std::endl;
    }
  else if( strcmp(operation.c_str(), "mean") == 0 )
    {
    std::cout << result / ct << std::endl;
    }
  else
    {
    // std::cout << "operation " << operation << std::endl;
    }
  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(varimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int FuseNImagesIntoNDVectorField(int argc, char *argv[])
{
  typedef itk::Vector<float, ImageDimension>       VectorType;
  typedef itk::Image<VectorType, ImageDimension>   FieldType;
  typedef itk::Image<float, ImageDimension>        ImageType;
  typedef  typename ImageType::IndexType           IndexType;
  typedef itk::ImageRegionIteratorWithIndex<FieldType>   Iterator;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       operation = std::string(argv[argct]);  argct++;
  std::string       fn1 = std::string(argv[argct]);   argct++;
  std::string       fn2 = "";
  std::string       fn3 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct++]);
    }
  if( argc > argct )
    {
    fn3 = std::string(argv[argct++]);
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer image3 = nullptr;
  if ( fn1.length() > 3 ) ReadImage<ImageType>(image1, fn1.c_str() );
  if ( fn2.length() > 3 ) ReadImage<ImageType>(image2, fn2.c_str() );
  if ( fn3.length() > 3 ) ReadImage<ImageType>(image3, fn3.c_str() );
  typename FieldType::Pointer vimage =
    AllocImage<FieldType>(image1->GetLargestPossibleRegion(),
      image1->GetSpacing(),
      image1->GetOrigin(),
      image1->GetDirection() );

  VectorType vec;
  vec.Fill(0);
  Iterator vfIter2( vimage,  vimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    IndexType ind = vfIter2.GetIndex();
    vec[0] = image1->GetPixel(ind);
    vec[1] = image2->GetPixel(ind);
    if ( ImageDimension >  2 ) vec[2] = image3->GetPixel(ind);
    vfIter2.Set(vec);
    }
  WriteImage<FieldType>(vimage, outname.c_str() );
  return EXIT_SUCCESS;
}



template <unsigned int ImageDimension>
int VImageMath(int argc, char *argv[])
{
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef FieldType                                                       ImageType;
  typedef typename ImageType::PixelType                                   PixelType;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       operation = std::string(argv[argct]);  argct++;
  std::string       fn1 = std::string(argv[argct]);   argct++;
  std::string       fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct++]);
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer varimage = nullptr;

  bool  isfloat = false;
  float floatval = 1.0;
  if( from_string<float>(floatval, fn2, std::dec) )
    {
    isfloat = true;
    }
  else
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }
  ReadImage<ImageType>(image1, fn1.c_str() );
  varimage = AllocImage<ImageType>(image1);

  if( strcmp(operation.c_str(), "mresample") == 0 && !isfloat )
    {
    typename ImageType::SpacingType spc = image2->GetSpacing();
    typedef itk::TranslationTransform<double, ImageDimension> TransformType0;
    typename TransformType0::Pointer             m_Transform0 = TransformType0::New();
    typename TransformType0::ParametersType trans = m_Transform0->GetParameters();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      trans[i] = 0;
      spc[i] = image1->GetSpacing()[i] * image1->GetLargestPossibleRegion().GetSize()[i]
        / image2->GetLargestPossibleRegion().GetSize()[i];
      }
    image2->SetSpacing(spc);
    image2->SetOrigin(image1->GetOrigin() );
    image2->SetDirection(image1->GetDirection() );
    m_Transform0->SetParameters(trans);
    // std::cout << " trans " << m_Transform0->GetParameters() << " Nspc " << image2->GetSpacing() << std::endl;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform( m_Transform0 );
    resample->SetInput( image2 );
    resample->SetOutputParametersFromImage(  image1 );
    typename ImageType::IndexType zeroind;  zeroind.Fill(0);
    resample->SetDefaultPixelValue( image1->GetPixel(zeroind) );
    resample->UpdateLargestPossibleRegion();
    image2 = resample->GetOutput();

    WriteImage<ImageType>(image2, outname.c_str() );
    return 0;
    }

  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < ImageDimension; i++ )
    {
    volumeelement *= static_cast<float>( varimage->GetSpacing()[i] );
    }

  PixelType result;
  result.Fill( 0 );
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    IndexType ind = vfIter2.GetIndex();
    PixelType pix2;
    if( isfloat )
      {
      pix2 = floatval;
      }
    else
      {
      pix2 = image2->GetPixel(ind);
      }
    PixelType pix1 = image1->GetPixel(ind);

    if( strcmp(operation.c_str(), "vm") == 0 )
      {
      result = pix1 * pix2;
      }
    else if( strcmp(operation.c_str(), "v+") == 0 )
      {
      result = pix1 + pix2;
      }
    else if( strcmp(operation.c_str(), "v-") == 0 )
      {
      result = pix1 - pix2;
      }
    else if( strcmp(operation.c_str(), "v/") == 0 )
      {
      // std::cout << "Operation v/ not implemented" << std::endl;
      //	result = // pix1 / pix2;
      }
    else if( strcmp(operation.c_str(), "vtotal") == 0 )
      {
      result += pix1 * pix2;
      }
    vfIter2.Set(result);
    }
  if( strcmp(operation.c_str(), "vtotal") == 0 )
    {
    std::cout << "total: " << result << " total-volume: " << result * volumeelement << std::endl;
    }
  else
    {
    // std::cout << "operation " << operation << std::endl;
    }
  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(varimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int SmoothTensorImage(int argc, char *argv[])
{
  typedef float                                  ValueType;
  typedef itk::DiffusionTensor3D<ValueType>      TensorType;
  typedef itk::Image<TensorType, ImageDimension> TensorImageType;

  typedef itk::RecursiveGaussianImageFilter<TensorImageType, TensorImageType>
    GaussianFilterType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct++]);
  std::string       operation = std::string(argv[argct++]);
  std::string       fn1 = std::string(argv[argct++]);
  float             sigma = atof(argv[argct++]);
  std::string       fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct++]);
    }

  typename TensorImageType::Pointer inDT = nullptr;
  ReadTensorImage<TensorImageType>(inDT, fn1.c_str() );

  typename GaussianFilterType::Pointer gFilter = GaussianFilterType::New();
  gFilter->SetInput( inDT );
  gFilter->SetSigma( sigma );
  gFilter->Update();

  WriteTensorImage<TensorImageType>(gFilter->GetOutput(), outname.c_str() );
  return 0;

}

template <unsigned int ImageDimension>
int TensorFunctions(int argc, char *argv[])
{
  typedef float PixelType;
  // Tensor4DType may contain b0
  typedef float                                              D4TensorType;
  typedef itk::SymmetricSecondRankTensor<float, 3>           TensorType;
  typedef typename itk::RGBPixel<unsigned char>              RGBType;
  typedef itk::Image<TensorType, ImageDimension>             TensorImageType;
  typedef itk::Image<D4TensorType, 4>                        D4TensorImageType;
  typedef typename TensorImageType::IndexType                IndexType;
  typedef itk::Image<PixelType, ImageDimension>              ImageType;
  typedef itk::Image<RGBType, ImageDimension>                ColorImageType;
  typedef itk::ImageFileWriter<ColorImageType>               ColorWriterType;
  typedef itk::ImageRegionIteratorWithIndex<TensorImageType> Iterator;

  typedef itk::Vector<float, ImageDimension>                 VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       operation = std::string(argv[argct]);  argct++;
  std::string       fn1 = std::string(argv[argct]);   argct++;
  std::string       fn2 = ""; // used for whichvec and mask file name below
  float             backgroundMD = 0.0; // mean diffusivity of isotropic background voxels

  typename TensorImageType::Pointer timage = nullptr;    // input tensor image
  typename ImageType::Pointer       vimage = nullptr;    // output scalar image
  typename ColorImageType::Pointer  cimage = nullptr;    // output color image
  typename VectorImageType::Pointer vecimage = nullptr;  // output vector image
  typename TensorImageType::Pointer toimage = nullptr;   // output tensor image
  typename ImageType::Pointer       mimage = nullptr;    // mask image

  if( strcmp(operation.c_str(), "4DTensorTo3DTensor") == 0 )
    {
    std::cout
      <<
      " Convert a 4D tensor to a 3D tensor --- if there are 7 components to the tensor, we throw away the first component b/c its probably b0 "
      << std::endl;
    itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn1.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn1.c_str() );
    imageIO->ReadImageInformation();
    unsigned int dim = imageIO->GetNumberOfDimensions();
    if( dim == 4 )
      {
      typename D4TensorImageType::Pointer d4img = nullptr;
      ReadImage<D4TensorImageType>(d4img, fn1.c_str() );
      unsigned int d4size = d4img->GetLargestPossibleRegion().GetSize()[3];
      if( d4size != 6 && d4size != 7 )
        {
        std::cout << " you should not be using this function if the input data is not a tensor. " << std::endl;
        std::cout
          <<
          " there is no way for us to really check if your use of this function is correct right now except checking the size of the 4th dimension which should be 6 or 7 (the latter if you store b0 in the first component) --- you should really store tensors not as 4D images but as 3D images with tensor voxel types. "
          << std::endl;
        throw std::exception();
        }
      typename TensorImageType::SizeType size;
      typename TensorImageType::RegionType tensorregion;
      typename TensorImageType::SpacingType spacing;
      typename TensorImageType::PointType origin;
      typename TensorImageType::DirectionType direction;
      for( unsigned int dd = 0; dd < ImageDimension; dd++ )
        {
        size[dd] = d4img->GetLargestPossibleRegion().GetSize()[dd];
        origin[dd] = d4img->GetOrigin()[dd];
        spacing[dd] = d4img->GetSpacing()[dd];
        for( unsigned int ee = 0; ee < ImageDimension; ee++ )
          {
          direction[dd][ee] = d4img->GetDirection()[dd][ee];
          }
        }
      tensorregion.SetSize(size);
      timage =
        AllocImage<TensorImageType>(tensorregion,
                                    spacing,
                                    origin,
                                    direction);

      // now iterate through & set the values of the tensors.
      Iterator tIter(timage, timage->GetLargestPossibleRegion() );
      for(  tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter )
        {
        typename TensorImageType::IndexType ind = tIter.GetIndex();
        typename D4TensorImageType::IndexType ind2;
        for( unsigned int dd = 0; dd < ImageDimension; dd++ )
          {
          ind2[dd] = ind[dd];
          }
        TensorType pix6 = tIter.Get();
        if( d4size == 6 )
          {
          for( unsigned int ee = 0; ee < d4size; ee++ )
            {
            ind2[3] = ee;
            if( ee ==  2 )
              {
              pix6[2] = d4img->GetPixel(ind2);
              }
            else if( ee ==  3 )
              {
              pix6[3] = d4img->GetPixel(ind2);
              }
            else
              {
              pix6[ee] = d4img->GetPixel(ind2);
              }
            }
          }
        // ITK-way
        //  xx, xy, yy, xz , yz, zz
        // VTK/other-way
        //  xx, xy, xz, yy , yz, zz
        else if( d4size == 7 )
          {
          for( unsigned int ee = 1; ee < d4size; ee++ )
            {
            ind2[3] = ee;
            pix6[ee - 1] = d4img->GetPixel(ind2);
            }
          }
        timage->SetPixel(ind, pix6);
        }
      WriteTensorImage<TensorImageType>(timage, outname.c_str(), false);
      return 0;
      }
    // std::cout << " cannot convert --- input image not 4D --- " << fn1 << std::endl;
    return 0;
    }

  if( strcmp(operation.c_str(), "ComponentTo3DTensor") == 0 )
    {
    std::string extension = std::string( ".nii.gz" );
    if( argc > argct )
      {
      extension = std::string( argv[argct++] );
      }

    typename ImageType::Pointer xx, xy, xz, yy, yz, zz;

    std::string fn1xx = fn1 + std::string( "xx" ) + extension;
    std::string fn1xy = fn1 + std::string( "xy" ) + extension;
    std::string fn1xz = fn1 + std::string( "xz" ) + extension;
    std::string fn1yy = fn1 + std::string( "yy" ) + extension;
    std::string fn1yz = fn1 + std::string( "yz" ) + extension;
    std::string fn1zz = fn1 + std::string( "zz" ) + extension;

    ReadImage<ImageType>( xx, fn1xx.c_str() );
    ReadImage<ImageType>( xy, fn1xy.c_str() );
    ReadImage<ImageType>( xz, fn1xz.c_str() );
    ReadImage<ImageType>( yy, fn1yy.c_str() );
    ReadImage<ImageType>( yz, fn1yz.c_str() );
    ReadImage<ImageType>( zz, fn1zz.c_str() );

    timage = AllocImage<TensorImageType>(xx);

    Iterator tIter( timage, timage->GetLargestPossibleRegion() );
    for( tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter )
      {
      typename TensorImageType::IndexType ind = tIter.GetIndex();
      TensorType pix6 = tIter.Get();

      pix6[0] = xx->GetPixel( ind );
      pix6[1] = xy->GetPixel( ind );
      pix6[2] = xz->GetPixel( ind );
      pix6[3] = yy->GetPixel( ind );
      pix6[4] = yz->GetPixel( ind );
      pix6[5] = zz->GetPixel( ind );

      tIter.Set( pix6 );
      }
    WriteTensorImage<TensorImageType>( timage, outname.c_str(), false );
    return 0;
    }

  if( strcmp(operation.c_str(), "ExtractComponentFrom3DTensor") == 0 )
    {
    ReadImage<TensorImageType>( timage, fn1.c_str() );

    unsigned int which = 0;
    if( argc > argct )
      {
      std::string component = std::string( argv[argct++] );
      if( component.find( "xx" ) != std::string::npos )
        {
        which = 0;
        }
      else if( component.find( "xy" ) != std::string::npos )
        {
        which = 1;
        }
      else if( component.find( "xz" ) != std::string::npos )
        {
        which = 2;
        }
      else if( component.find( "yy" ) != std::string::npos )
        {
        which = 3;
        }
      else if( component.find( "yz" ) != std::string::npos )
        {
        which = 4;
        }
      else if( component.find( "zz" ) != std::string::npos )
        {
        which = 5;
        }
      else
        {
        std::cout << "Unrecognized component.  Need to specify "
                  << "xx, xy, xz, yy, yz, or zz";
        return EXIT_FAILURE;
        }
      }
    else
      {
      // std::cout << "Error:  need to specify component (xx, xy, xz, yy, yz, zz)";
      return EXIT_FAILURE;
      }

    typename ImageType::Pointer componentImage =
      AllocImage<ImageType>(timage, 0.0);

    Iterator tIter( timage, timage->GetLargestPossibleRegion() );
    for( tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter )
      {
      typename TensorImageType::IndexType ind = tIter.GetIndex();
      TensorType pix6 = tIter.Get();

      componentImage->SetPixel( ind, pix6[which] );
      }
    WriteTensorImage<TensorImageType>( timage, outname.c_str(), false );
    return 0;
    }

  unsigned int whichvec = ImageDimension - 1;
  if( argc > argct )
    {
    // Might be a file or a number for whichvec
    fn2 = std::string(argv[argct]);
    try
      {
      whichvec = std::stoi(fn2.c_str() );
      }
    catch(std::invalid_argument& itkNotUsed(e))
      {
      // arg is not whichvec
      }
    argct++;
    }

  if ( argc > argct )
    {
    // Throws exception if arg is not a valid float
    backgroundMD = std::stof( std::string(argv[argct]).c_str() );
    argct++;
    }

  ReadTensorImage<TensorImageType>(timage, fn1.c_str(), false);

  if( strcmp(operation.c_str(), "TensorIOTest") == 0 )
    {
    // std::cout << " test function for tensor I/O " << std::endl;
    WriteTensorImage<TensorImageType>(timage, outname.c_str(), false);
    return 0;
    }
  // std::cout << " imagedir " << timage->GetDirection() << std::endl;

  if( strcmp(operation.c_str(), "TensorColor") == 0 )
    {
    cimage =
      AllocImage<ColorImageType>(timage);

    if( argc > 5 )
      {
      // std::cout << "Using mask image: " << fn2 << std::endl;
      ReadImage<ImageType>(mimage, fn2.c_str() );
      }

    }
  else if( strcmp(operation.c_str(), "TensorMask") == 0 )
    {
    // std::cout << "Using mask image: " << fn2 << std::endl;
    ReadImage<ImageType>(mimage, fn2.c_str() );

    typename TensorImageType::PixelType zero;

    zero.Fill(0);

    toimage = AllocImage<TensorImageType>(timage, zero);
    }
  else if( strcmp(operation.c_str(), "TensorToVector") == 0 )
    {
    VectorType zero;  zero.Fill(0);
    vecimage = AllocImage<VectorImageType>(timage, zero);
    }
  else if( (strcmp(operation.c_str(), "TensorToPhysicalSpace") == 0) ||
           (strcmp(operation.c_str(), "TensorToLocalSpace") == 0)  ||
           (strcmp(operation.c_str(), "ValidTensor") == 0)  )
    {
    typename TensorImageType::PixelType zero;
    zero.Fill(0);
    toimage = AllocImage<TensorImageType>(timage, zero);
    }
  else
    {
    vimage = AllocImage<ImageType>(timage);
    }

  TensorType backgroundTensor;  // for masking background tensors
  for( unsigned int i = 0; i < 6; i++ )
    {
    backgroundTensor[i] = 0.0;
    }

  if ( backgroundMD > 0.0f )
    {
    // std::cout << "Setting background voxels to isotropic tensors with mean diffusivity " << backgroundMD << std::endl;
    backgroundTensor[0] = backgroundMD;
    backgroundTensor[3] = backgroundMD;
    backgroundTensor[5] = backgroundMD;
    }


  RGBType rgbZero;

  rgbZero[0] = 0;
  rgbZero[1] = 0;
  rgbZero[2] = 0;

  Iterator tIter(timage, timage->GetLargestPossibleRegion() );
  for(  tIter.GoToBegin(); !tIter.IsAtEnd(); ++tIter )
    {
    IndexType ind = tIter.GetIndex();
    float     result = 0;

    if( strcmp(operation.c_str(), "TensorFA") == 0 )
      {
      result = 0.0;
      if( IsRealTensor<TensorType>( tIter.Value() ) )
        {
        result = GetTensorFA<TensorType>(tIter.Value() );
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorMeanDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 0);
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorRadialDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 2);
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorEigenvalue") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 3 + whichvec);
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorAxialDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 5);
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorFANumerator") == 0 )
      {
      result = GetTensorFANumerator<TensorType>(tIter.Value() );
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorFADenominator") == 0 )
      {
      result = GetTensorFADenominator<TensorType>(tIter.Value() );
      if( std::isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorColor") == 0 )
      {
      if( argc > 5 )
        {
        if( mimage->GetPixel( tIter.GetIndex() ) > 0 )
          {
          cimage->SetPixel(ind, GetTensorRGB<TensorType>(tIter.Value() ) );
          }
        else
          {
          cimage->SetPixel(ind, rgbZero);
          }
        }
      else
        {
        RGBType rgb = GetTensorRGB<TensorType>(tIter.Value() );
        cimage->SetPixel(ind, rgb);
        }
      }
    else if( strcmp(operation.c_str(), "TensorMask") == 0 )
      {
      float maskVal = mimage->GetPixel(ind);

      if( maskVal > 0.0f )
        {
        toimage->SetPixel( ind, tIter.Value() );
        }
      else
        {
        toimage->SetPixel( ind, backgroundTensor );
        }

      }
    else if( strcmp(operation.c_str(), "TensorToVector") == 0 )
      {
      VectorType vv = GetTensorPrincipalEigenvector<TensorType>(tIter.Value(), whichvec);
      vecimage->SetPixel(ind, vv);
      }

    else if( strcmp(operation.c_str(), "TensorToVectorComponent") == 0 )
      {
      if( whichvec <= 2 )
        {
        VectorType vv = GetTensorPrincipalEigenvector<TensorType>(tIter.Value(), 2);
        vimage->SetPixel(ind, vv[whichvec]);
        }
      else if( whichvec > 2 && whichvec < 9 )
        {
        vimage->SetPixel(ind, tIter.Value()[whichvec - 3]);
        }
      }
    else if( strcmp(operation.c_str(), "TensorToPhysicalSpace") == 0 )
      {
      typename TensorType::EigenValuesArrayType eigenValues;
      typename TensorType::EigenVectorsMatrixType eigenVectors;
      typename TensorType::EigenVectorsMatrixType eigenVectorsPhysical;
      typename TensorType::EigenVectorsMatrixType eigenValuesMatrix;
      eigenValuesMatrix.Fill( 0.0 );
      tIter.Value().ComputeEigenAnalysis( eigenValues, eigenVectors );
      for( unsigned int i = 0; i < 3; i++ )
        {
        eigenValuesMatrix(i, i) = fabs( eigenValues[i] );

        itk::Vector<float, 3> ev;
        for( unsigned int j = 0; j < 3; j++ )
          {
          ev[j] = eigenVectors(j, i);
          }

        itk::Vector<float, 3> evp;
        timage->TransformLocalVectorToPhysicalVector( ev, evp );

        itk::Vector<float, 3> evl;
        timage->TransformPhysicalVectorToLocalVector( evp, evl );
        for( unsigned int j = 0; j < 3; j++ )
          {
          eigenVectorsPhysical(j, i) = evp[j];
          }
        }

      typename TensorType::MatrixType::InternalMatrixType phyTensor
        = eigenVectorsPhysical.GetTranspose() * eigenValuesMatrix.GetVnlMatrix() * eigenVectorsPhysical.GetVnlMatrix();

      TensorType oTensor = Matrix2Vector<TensorType, TensorType::MatrixType::InternalMatrixType>( phyTensor );
      toimage->SetPixel( tIter.GetIndex(), oTensor );
      }
    else if( strcmp(operation.c_str(), "TensorToLocalSpace") == 0 )
      {
      typename TensorType::EigenValuesArrayType eigenValues;
      typename TensorType::EigenVectorsMatrixType eigenVectors;
      typename TensorType::EigenVectorsMatrixType eigenValuesMatrix;
      eigenValuesMatrix.Fill( 0.0 );
      tIter.Value().ComputeEigenAnalysis( eigenValues, eigenVectors );
      for( unsigned int i = 0; i < 3; i++ )
        {
        eigenValuesMatrix(i, i) = fabs( eigenValues[i] );

        itk::Vector<float, 3> ev;
        for( unsigned int j = 0; j < 3; j++ )
          {
          ev[j] = eigenVectors(j, i);
          }

        itk::Vector<float, 3> evp;
        timage->TransformPhysicalVectorToLocalVector( ev, evp );
        for( unsigned int j = 0; j < 3; j++ )
          {
          eigenVectors(j, i) = evp[j];
          }
        }

      typename TensorType::MatrixType::InternalMatrixType lclTensor
        = eigenVectors.GetTranspose() * eigenValuesMatrix.GetVnlMatrix() * eigenVectors.GetVnlMatrix();

      TensorType oTensor = Matrix2Vector<TensorType, TensorType::MatrixType::InternalMatrixType>( lclTensor );
      toimage->SetPixel( tIter.GetIndex(), oTensor );
      }
    else if( strcmp(operation.c_str(), "ValidTensor") == 0 )
      {
      typename TensorType::EigenValuesArrayType eigenValues;
      typename TensorType::EigenVectorsMatrixType eigenVectors;
      typename TensorType::EigenVectorsMatrixType eigenValuesMatrix;
      eigenValuesMatrix.Fill( 0.0 );
      tIter.Value().ComputeEigenAnalysis( eigenValues, eigenVectors );
      // NOT USED bool hasNeg = false;

      typename TensorType::MatrixType::InternalMatrixType lclTensor
        = eigenVectors.GetTranspose() * eigenValuesMatrix.GetVnlMatrix() * eigenVectors.GetVnlMatrix();

      TensorType oTensor = Matrix2Vector<TensorType, TensorType::MatrixType::InternalMatrixType>( lclTensor );
      toimage->SetPixel( tIter.GetIndex(), oTensor );
      }
    }
  if( strcmp(operation.c_str(), "TensorColor") == 0 )
    {
    typename ColorWriterType::Pointer cwrite = ColorWriterType::New();
    cwrite->SetInput(cimage);
    cwrite->SetFileName(outname.c_str() );
    cwrite->Update();
    }
  else if( strcmp(operation.c_str(), "TensorToVector") == 0 )
    {
    WriteImage<VectorImageType>(vecimage, outname.c_str() );
    }
  else if( (strcmp(operation.c_str(), "TensorToPhysicalSpace") == 0) ||
           (strcmp(operation.c_str(), "TensorToLocalSpace") == 0 ) ||
           (strcmp(operation.c_str(), "TensorMask") == 0 ) ||
           (strcmp(operation.c_str(), "ValidTensor") == 0 ) )
    {
    WriteTensorImage<TensorImageType>(toimage, outname.c_str(), false );
    }
  else
    {
    // std::cout << "Writing scalar image" << std::endl;
    WriteImage<ImageType>(vimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int CompareHeadersAndImages(int argc, char *argv[])
{
  typedef float PixelType;
  //  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  bool isfloat = false;
  try
    {
    ReadImage<ImageType>( image2, fn2.c_str() );
    }
  catch( ... )
    {
    // std::cout << " Error reading " << fn2 << std::endl;
    isfloat = true;
    }

  float floatval = 1.0;
  if( isfloat )
    {
    floatval = atof(argv[argct]);
    }
  else
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }

  try
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  catch( ... )
    {
    // std::cout << " read 1 error ";
    }

  // compute error in spacing, in orientation and in offset
  unsigned int failure = 0;
  float        sperr = 0, merr = 0, operr = 0, orsignerr = 0;

  typename ImageType::SpacingType sp1, sp2;
  sp1 = image1->GetSpacing();
  sp2 = image2->GetSpacing();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    float temp = sp1[i] - sp2[i];
    sperr += temp * temp;
    }
  // std::cout << " SpacingError: " << sqrt(sperr) << std::endl;

  typename ImageType::PointType op1, op2;
  op1 = image1->GetOrigin();
  op2 = image2->GetOrigin();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    float temp = op1[i] - op2[i];
    operr += temp * temp;
    if( op1[i] > 0 && op2[i] <= 0 )
      {
      orsignerr += 1;
      }
    else if( op1[i] < 0 && op2[i] >= 0 )
      {
      orsignerr += 1;
      }
    }
  // std::cout << " OriginError: " << sqrt(operr) << std::endl;
  // std::cout << " OriginSignError: " << orsignerr << std::endl;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      float temp = image1->GetDirection()[i][j] - image2->GetDirection()[i][j];
      merr += temp * temp;
      }
    }
  // std::cout << " OrientError: " << sqrt(merr) << std::endl;

  bool samesize = true;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( image1->GetLargestPossibleRegion().GetSize()[i] != image2->GetLargestPossibleRegion().GetSize()[i] )
      {
      samesize = false;
      }
    }

  if( samesize )
    {
    // run a quick registration
    if( false )
      {
      typedef typename itk::ImageMomentsCalculator<ImageType> ImageCalculatorType;
      typename ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();
      calculator->SetImage(  image1 );

      typename ImageCalculatorType::VectorType fixed_center;
      fixed_center.Fill(0);
      typename ImageCalculatorType::VectorType moving_center;
      moving_center.Fill(0);
      try
        {
        calculator->Compute();
        fixed_center = calculator->GetCenterOfGravity();
        calculator->SetImage(  image2 );
        try
          {
          calculator->Compute();
          moving_center = calculator->GetCenterOfGravity();
          }
        catch( ... )
          {
          // std::cout << " zero image2 error ";
          fixed_center.Fill(0);
          }
        }
      catch( ... )
        {
        // std::cout << " zero image1 error ";
        }

      typedef itk::TranslationTransform<double, ImageDimension> TransformType0;
      typename TransformType0::Pointer             m_Transform0 = TransformType0::New();
      typename TransformType0::ParametersType trans = m_Transform0->GetParameters();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        trans[i] = moving_center[i] - fixed_center[i];
        }
      m_Transform0->SetParameters(trans);
      // std::cout << " trans " << m_Transform0->GetParameters() << std::endl;
      typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
      typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
      resample->SetTransform( m_Transform0 );
      resample->SetInput( image2 );
      resample->SetOutputParametersFromImage(  image1 );
      typename ImageType::IndexType zeroind;  zeroind.Fill(0);
      resample->SetDefaultPixelValue( image1->GetPixel(zeroind) );
      resample->UpdateLargestPossibleRegion();
      typename ImageType::Pointer varimage = resample->GetOutput();
      }
    float         i1norm = 0, i2norm = 0, i1i2norm = 0; // i3norm=0,i1i3norm=0;
    unsigned long ct1 = 1, ct2 = 1, ct12 = 1;           // ct3=1,ct13=1;
    Iterator      vfIter2( image1,  image1->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      IndexType ind = vfIter2.GetIndex();
      float     pix2;
      if( isfloat )
        {
        pix2 = floatval;
        }
      else
        {
        pix2 = image2->GetPixel(ind);
        }
      if( std::isnan(pix2) || std::isinf(pix2) )
        {
        pix2 = 0; image2->SetPixel(ind, 0);
        }
      //      float pix3=varimage->GetPixel(ind);
      float pix1 = image1->GetPixel(ind);
      if( pix1 > 0  )
        {
        i1norm += static_cast<float>( fabs(pix1) ); ct1++;
        }
      if( pix2 > 0  )
        {
        i2norm += static_cast<float>( fabs(pix2) ); ct2++;
        }
      // if (pix3 > 0  )  {  i3norm+=fabs(pix3); ct3++; }
      }
    float mean1 = i1norm / ct1;
    if( itk::Math::FloatAlmostEqual( mean1, 0.0f ) )
      {
      mean1 = 1;
      }
    float mean2 = i2norm / ct2;
    if( itk::Math::FloatAlmostEqual( mean2, 0.0f ) )
      {
      mean2 = 1;
      }
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      IndexType ind = vfIter2.GetIndex();
      float     pix2;
      if( isfloat )
        {
        pix2 = floatval;
        }
      else
        {
        pix2 = image2->GetPixel(ind) / mean2;
        }
      //      float pix3=vfIter2.Get()/mean2;
      float pix1 = image1->GetPixel(ind) / mean1;

      if( pix1 > 0 || pix2 > 0 )
        {
        i1i2norm += static_cast<float>( std::fabs(pix1 - pix2) );
        ct12++;
        }
      //      if ( pix1 > 0 || pix3 > 0)
      // {
      //  i1i3norm+=fabs(pix1-pix3);
      //  ct13++;
      // }
      }
    //  float idice1 = 1.0 - 2.0*i1i3norm/ct13 /  (  i1norm/ct1 + i3norm / ct3 );
    // std::cout << " DiceImageDifference: " << idice0 << " IntensityDifference: " << i1i2norm << std::endl;
    // std::cout << " CenterOfMassTransImageDifference: " << idice1 << " and " << i1i3norm << std::endl;
    }

  typename ImageType::PointType fixedorig = image2->GetOrigin();
  if( orsignerr > 0 )
    {
    failure = 1;
    // now fix the error
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( (op1[i] > 0 && op2[i] < 0 ) || (op1[i] < 0 && op2[i] > 0) )
        {
        fixedorig[i] *= (-1.0);
        }
      else if( (op1[i] > 0 && itk::Math::FloatAlmostEqual( op2[i], 0.0 ) ) || ( op1[i] < 0 && itk::Math::FloatAlmostEqual( op2[i], 0.0 ) ) )
        {
        fixedorig[i] = image1->GetOrigin()[i];
        }
      }
    }
  image2->SetOrigin(fixedorig);
  //  if ( sqrt(operr) > 100) failure=
  if( sqrt(merr) >= 0.7 && failure == 0 )
    {
    failure = 2;
    }
  if( sqrt(merr) >= 0.7 && failure != 0 )
    {
    failure = 3;
    }
  if( failure > 1 )
    {
    image2->SetDirection(image1->GetDirection() );
    }
  // write repaired images
  WriteImage<ImageType>( image2, outname.c_str() );
  // std::cout << "  FailureState: " << failure << " for " << fn2  << std::endl;
  return failure;
}

// template<typename TImage>
// typename TImage::Pointer
// SegmentKMeans(typename TImage::Pointer image , unsigned int nclasses)
// {
//
//   typedef TImage ImageType;
//   typedef typename TImage::PixelType PixelType;
//   enum { ImageDimension = ImageType::ImageDimension };
//   typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
//
//   typedef itk::Statistics::ScalarImageToListAdaptor< ImageType >   AdaptorType;
//
//   typename AdaptorType::Pointer adaptor = AdaptorType::New();
//
//   adaptor->SetImage( image );
//
//   // Define the Measurement vector type from the AdaptorType
//   typedef typename AdaptorType::MeasurementVectorType  MeasurementVectorType;
//
//   // Create the K-d tree structure
//   typedef itk::Statistics::WeightedCentroidKdTreeGenerator<
//                                                       AdaptorType >
//                                                               TreeGeneratorType;
//
//   typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
//
//   treeGenerator->SetSample( adaptor );
//   treeGenerator->SetBucketSize( 16 );
//   treeGenerator->Update();
//
//   typedef typename TreeGeneratorType::KdTreeType TreeType;
//   typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
//
//   typename EstimatorType::Pointer estimator = EstimatorType::New();
//
//   typename EstimatorType::ParametersType initialMeans( nclasses );
//
//   Iterator vfIter2( image,  image->GetLargestPossibleRegion() );
//   double mx =-1.e12, mn=1.e12;
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//       double px = vfIter2.Get();
//       if (px > mx) mx=px;
//       else if (px < mn) mn=px;
//     }
//   float range=(mx-mn);
//
//   float bins=1.0/((float)nclasses+1.);
//   for (unsigned int i=0;  i<nclasses; i++)
//     {
//     initialMeans[i]=mn+(0+i*bins)*range;
//     // std::cout << " Initial Means " << initialMeans[i] << " ";
//     }
//   // std::cout << std::endl;
//   estimator->SetParameters( initialMeans );
//
//   estimator->SetKdTree( treeGenerator->GetOutput() );
//   estimator->SetMaximumIteration( 200 );
//   estimator->SetCentroidPositionChangesThreshold(0.0);
//   estimator->StartOptimization();
//
//   typename EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
//
//
//   typename ImageType::Pointer varimage=ImageType::New();
//   varimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   varimage->SetBufferedRegion( image->GetLargestPossibleRegion() );
//   varimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   varimage->Allocate();
//   varimage->SetSpacing(image->GetSpacing());
//   varimage->SetOrigin(image->GetOrigin());
//   varimage->SetDirection(image->GetDirection());
//
// //  float var=sqrt(range);
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//     double px = vfIter2.Get();
//       unsigned int best=0;
//       float mindist=1.e9;
//       for ( unsigned int i = 0 ; i < nclasses ; ++i )
//     {
//     float dist=fabs(px-estimatedMeans[i]);
//     if (dist < mindist) { mindist=dist;  best=i; }
//     //vec[i]=exp(-1.0*dist*dist/var);
//     }
//       varimage->SetPixel(vfIter2.GetIndex(),best+1);
// //      vecimage->SetPixel(vfIter2.GetIndex(),vec);
//     }
//
//   return varimage;
//
//
// }

// template<typename TImage>
// typename TImage::Pointer
// BayesianSegmentation(typename TImage::Pointer image , unsigned int nclasses,  std::string priorfn ,  unsigned int
// nsmooth = 2 )
// {
//
//   typedef TImage ImageType;
//   typedef typename TImage::PixelType PixelType;
//   enum { ImageDimension = ImageType::ImageDimension };
//   typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
//
//   //  const unsigned int ImageDimension = AvantsImageDimension;
//   typedef itk::Vector<float,ImageDimension>         VectorType;
//   typedef itk::Image<VectorType,ImageDimension>     FieldType;
//   typedef itk::ImageFileReader<ImageType> readertype;
//   typedef itk::ImageFileWriter<ImageType> writertype;
//   typedef  typename ImageType::IndexType IndexType;
//   typedef  typename ImageType::SizeType SizeType;
//   typedef  typename ImageType::SpacingType SpacingType;
//   typedef itk::AffineTransform<double,ImageDimension>   AffineTransformType;
//   typedef itk::LinearInterpolateImageFunction<ImageType,double>  InterpolatorType1;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>  InterpolatorType2;
//   typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
//
//
//
//   typedef itk::Statistics::ScalarImageToListAdaptor< ImageType >   AdaptorType;
//   typename AdaptorType::Pointer adaptor = AdaptorType::New();
//   adaptor->SetImage( image );
//   // Define the Measurement vector type from the AdaptorType
//   typedef typename AdaptorType::MeasurementVectorType  MeasurementVectorType;
//   // Create the K-d tree structure
//   typedef itk::Statistics::WeightedCentroidKdTreeGenerator<
//                                                       AdaptorType >
//                                                               TreeGeneratorType;
//   typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
//   treeGenerator->SetSample( adaptor );
//   treeGenerator->SetBucketSize( 16 );
//   treeGenerator->Update();
//   typedef typename TreeGeneratorType::KdTreeType TreeType;
//   typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
//   typename EstimatorType::Pointer estimator = EstimatorType::New();
//   typename EstimatorType::ParametersType initialMeans( nclasses );
//   Iterator vfIter2( image,  image->GetLargestPossibleRegion() );
//   double mx =-1.e12, mn=1.e12;
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//       double px = vfIter2.Get();
//       if (px > mx) mx=px;
//       else if (px < mn) mn=px;
//     }
//   float range=(mx-mn);
//   float bins=1.0/((float)nclasses+1.);
//   for (unsigned int i=0;  i<nclasses; i++)
//     {
//     initialMeans[i]=mn+(0+i*bins)*range;
//     }
//   estimator->SetParameters( initialMeans );
//   estimator->SetKdTree( treeGenerator->GetOutput() );
//   estimator->SetMaximumIteration( 200 );
//   estimator->SetCentroidPositionChangesThreshold(0.0);
//   estimator->StartOptimization();
//   typename EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
//
//   typename ImageType::Pointer varimage=ImageType::New();
//   varimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   varimage->SetBufferedRegion( image->GetLargestPossibleRegion() );
//   varimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   varimage->SetSpacing(image->GetSpacing());
//   varimage->SetOrigin(image->GetOrigin());
//   varimage->SetDirection(image->GetDirection());
//   varimage->Allocate();
//
//   std::vector<double>  estimatedVar(nclasses,0);
//   std::vector<double>  estimatedCounts(nclasses,0);
//
// //  float var=sqrt(range);
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//     double px = vfIter2.Get();
//       unsigned int best=0;
//       float mindist=1.e9;
//       for ( unsigned int i = 0 ; i < nclasses ; ++i )
//     {
//     float dist=fabs(px-estimatedMeans[i]);
//     if (dist < mindist) { mindist=dist;  best=i; }
//     //vec[i]=exp(-1.0*dist*dist/var);
//     }
//       varimage->SetPixel(vfIter2.GetIndex(),best);
//     }
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//     double px = vfIter2.Get();
//     unsigned int i = (unsigned int) varimage->GetPixel(vfIter2.GetIndex());
//     estimatedCounts[i]=estimatedCounts[i]+1;
//     float dist=(px -estimatedMeans[i]);
//     estimatedVar[i]+=dist*dist;
//     }
//
//   for (unsigned int i=0; i<nclasses; i++)
//     {
//     float   ct= estimatedCounts[i];
//     if (ct > 0) estimatedVar[i]=(estimatedVar[i])/ct;
//     else estimatedVar[i]=0;
//     // std::cout << " Sample SD Ests " << sqrt(estimatedVar[i]) << " Mean " << estimatedMeans[i] <<  std::endl;
//     }
//
//   typedef float InputPixelType;
//   typedef itk::VectorImage< InputPixelType, ImageDimension > InputImageType;
//   typedef  float  LabelType;
//   typedef float          PriorType;
//   typedef float          PosteriorType;
//   typedef itk::BayesianClassifierImageFilter<
//                               InputImageType,LabelType,
//                               PosteriorType,PriorType >   ClassifierFilterType;
//
//   typename InputImageType::Pointer vecImage = InputImageType::New();
//   typedef typename InputImageType::PixelType VecPixelType;
//   vecImage->SetSpacing(image->GetSpacing());
//   vecImage->SetOrigin(image->GetOrigin());
//   vecImage->SetRegions( image->GetLargestPossibleRegion() );
//   vecImage->SetVectorLength(nclasses);
//   vecImage->Allocate();
//   VecPixelType vvv(nclasses);
//   vvv.Fill(0);
//   vecImage->FillBuffer(vvv);
//
//   for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )  {
//      double px = vfIter2.Get();
//      VecPixelType  probs(nclasses);
//      float total=0;
//      for (unsigned int i=0; i<nclasses; i++)
//        {
//        probs[i]=exp(-1.0*(estimatedMeans[i]-px)*(estimatedMeans[i]-px)/(2.0*estimatedVar[i]));
//        total+=probs[i];
//        }
//      for (unsigned int i=0; i<nclasses; i++)
//        {
//        if (total>0) probs[i]/=total;
//        }
//      vecImage->SetPixel( vfIter2.GetIndex(), probs);
//   }
//
//   typename ClassifierFilterType::Pointer filter = ClassifierFilterType::New();
//
//
//   typedef itk::ImageFileReader< InputImageType >  ReaderType;
//   typename ReaderType::Pointer reader =  ReaderType::New();
//   typename InputImageType::Pointer priors= nullptr;
//
//   if (priorfn.length()  > 3 )
//     {
//     // std::cout << " Setting Priors " << priorfn << std::endl;
//     bool geometric=false;
//     if ( strcmp(priorfn.c_str(),"Geometric") == 0) geometric=true;
//     if (geometric)
//       {
//       // std::cout <<" Using a geometric thickness prior to aid cortical segmentation " << std::endl;
//       typename ImageType::Pointer outbrainmask = BinaryThreshold<TImage>(0,nclasses-3,1,varimage);
//       typename ImageType::Pointer inwmask = BinaryThreshold<TImage>(nclasses-1,nclasses,1,varimage);
//       typename ImageType::Pointer outwmask = BinaryThreshold<TImage>(0,nclasses-2,1,varimage);
//       typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType >  FilterType;
//       typename  FilterType::Pointer distmap = FilterType::New();
//       distmap->InputIsBinaryOn();
//       distmap->SetUseImageSpacing(true);
//       distmap->SetInput(outbrainmask);
//       distmap->Update();
//       typename ImageType::Pointer distcortex=distmap->GetOutput();
//
//       typename  FilterType::Pointer distmap2 = FilterType::New();
//       distmap2->InputIsBinaryOn();
//       distmap2->SetUseImageSpacing(true);
//       distmap2->SetInput(inwmask);
//       distmap2->Update();
//       typename ImageType::Pointer distwm=distmap2->GetOutput();
//
//       typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType,ImageType >  dgf;
//       typename dgf::Pointer lfilter = dgf::New();
//       lfilter->SetSigma(1.3);
//       lfilter->SetInput(distwm);
//       lfilter->Update();
//       typename ImageType::Pointer image2=lfilter->GetOutput();
//       typedef itk::RescaleIntensityImageFilter<ImageType,ImageType > RescaleFilterType;
//       typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
//       rescaler->SetOutputMinimum(   0 );
//       rescaler->SetOutputMaximum( 1 );
//       rescaler->SetInput( image2 );
//       rescaler->Update();
//       typename ImageType::Pointer sulci=  rescaler->GetOutput();
//
//       priors= InputImageType::New();
//       typedef typename InputImageType::PixelType VecPixelType;
//       priors->SetSpacing(image->GetSpacing());
//       priors->SetOrigin(image->GetOrigin());
//       priors->SetDirection(image->GetDirection());
//       priors->SetRegions( image->GetLargestPossibleRegion() );
//       priors->SetVectorLength(nclasses);
//       priors->Allocate();
//       // std::cout <<" Allocated " << std::endl;
//
//       for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
// //    // std::cout <<" ind " <<vfIter2.GetIndex() << std::endl;
// //    float outbrain = outbrainmask->GetPixel( vfIter2.GetIndex());
// //    float inw = inwmask->GetPixel( vfIter2.GetIndex());
//     float distance = distcortex->GetPixel( vfIter2.GetIndex());
//     float wdistance = distwm->GetPixel( vfIter2.GetIndex());
//     VecPixelType  probs(nclasses);
//     probs.Fill(1.0/(float)nclasses);
//     VecPixelType  posteriors=vecImage->GetPixel(vfIter2.GetIndex());
//     if (nclasses > 2)
//     {
//     float dvar=2;
//     float basedist=4;
//     float distmag=basedist-distance;
// //    if (distmag > 0) distmag=0;
//     distmag*=distmag;
//     float gdistprob=1.0-exp(-1.0*distmag/dvar);
//
//     float wdistprob=1.0/(1.0+exp(-1.0*distance/5));
//
//     float wdistmag=basedist-wdistance;
//     if (wdistmag > 0) wdistmag=0;
//     wdistmag*=wdistmag;
//     float gdistprob2=exp(-1.0*wdistmag/dvar);
//
//     float sulcprob=sulci->GetPixel( vfIter2.GetIndex());
//     float sdiff=(0.25-sulcprob);
//     if (sdiff > 0) sdiff=0;
//     sdiff*=sdiff;
//     sulcprob=exp(-1.0*sdiff/0.5);
// //    // std::cout << " Sulc " << sulcprob << std::endl;
// //    bool test = (outbrain < 1 && inw > 1);
//     if (  true  )
//       {
//       for (unsigned int i=0; i<nclasses; i++)
//         {
//         if ( i == (nclasses-2)) probs[i]=gdistprob*(gdistprob2);
//         else if ( i == (nclasses-1)) probs[i]=wdistprob*wdistprob;
//         //else
//           if ( i == (nclasses-3)) probs[i]=sulcprob;
// //        else if (i > 0)  probs[i]=sulcprob;
//         }
//       }
//     else
//       {
//       for (unsigned int i=0; i<nclasses; i++) probs[i]=posteriors[i];//1.0/(float)nclasses;
//       }
//     for (unsigned int i=0; i<nclasses; i++) posteriors[i]*=probs[i];
//     }
//       float prtotal=0;
//       float pototal=0;
//       for (unsigned int i=0; i<nclasses; i++) {  prtotal+=probs[i];  pototal+=posteriors[i]; }
//
//       for (unsigned int i=0; i<nclasses; i++)
//     {
//     if (prtotal>0) probs[i]/=prtotal; else probs[i]=1.0/(float)nclasses;
//     if (pototal>0) posteriors[i]/=pototal; else posteriors[i]=1.0/(float)nclasses;
//     }
//       priors->SetPixel( vfIter2.GetIndex(), probs);
//       vecImage->SetPixel( vfIter2.GetIndex(), posteriors);
//       }
//       // std::cout << " ok " << std::endl;
//
//
//       }
//     else
//       {
//       reader->SetFileName( priorfn.c_str() );
//       reader->Update();
//       priors=reader->GetOutput();
//
//       for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
//     VecPixelType  posteriors=vecImage->GetPixel(vfIter2.GetIndex());
//     VecPixelType  probs=priors->GetPixel(vfIter2.GetIndex());
//     for (unsigned int i=0; i<nclasses; i++) posteriors[i]*=probs[i];
//
//     float prtotal=0;
//     float pototal=0;
//     for (unsigned int i=0; i<nclasses; i++) {  prtotal+=probs[i];  pototal+=posteriors[i]; }
//
//     for (unsigned int i=0; i<nclasses; i++)
//       {
//       if (prtotal>0) probs[i]/=prtotal; else probs[i]=1.0/(float)nclasses;
//       if (pototal>0) posteriors[i]/=pototal; else posteriors[i]=1.0/(float)nclasses;
//       }
//     vecImage->SetPixel( vfIter2.GetIndex(), posteriors);
//     }
//       }
// //    if (priors) filter->SetInput( 1,  priors ); // Bug --
// //    classification filter does not actually use priors
//     } else // std::cout << " No Priors " << std::endl;
//
//   filter->SetInput(  vecImage );
//
//
//   if( nsmooth >= 1  )
//     {
//     // std::cout << " Smoothing Iterations:  " << nsmooth << std::endl;
//     filter->SetNumberOfSmoothingIterations( nsmooth );
//     typedef typename ClassifierFilterType::ExtractedComponentImageType ExtractedComponentImageType;
//     typedef itk::DiscreteGaussianImageFilter<
//       ExtractedComponentImageType, ExtractedComponentImageType >  SmoothingFilterType;
//     typedef itk::BilateralImageFilter<
//       ExtractedComponentImageType, ExtractedComponentImageType >  SmoothingFilterType2;
//     typename SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
//     smoother->SetVariance(1.0);
//     smoother->SetUseImageSpacingOff();
//     //smoother->SetDomainSigma(1);
//     //smoother->SetRangeSigma(1);
//     filter->SetSmoothingFilter( smoother );
//     }
//
//   // SET FILTER'S PRIOR PARAMETERS
//   // do nothing here to default to uniform priors
//   // otherwise set the priors to some user provided values
//
//   //
//   // Setup writer.. Rescale the label map to the dynamic range of the
//   // datatype and write it
//   //
//   filter->Update();
//
// return  filter->GetOutput();
//
// }

template <unsigned int ImageDimension>
int NegativeImage(int /*argc */, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>( image1, fn1.c_str() );
  Iterator vfIter2( image1,  image1->GetLargestPossibleRegion() );
  double   mx = -1.e12, mn = 1.e12;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    double px = vfIter2.Get();
    if( px > mx )
      {
      mx = px;
      }
    else if( px < mn )
      {
      mn = px;
      }
    }
  if( itk::Math::FloatAlmostEqual( mx, mn ) )
    {
    mx = 1; mn = 0;
    }
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    IndexType ind = vfIter2.GetIndex();
    double    pix = image1->GetPixel(ind);
    pix = (pix - mn) / (mx - mn);
    pix = (1.0 - pix) * (mx - mn);

    vfIter2.Set(pix);
    }

  WriteImage<ImageType>( image1, outname.c_str() );

  return 0;
}

// template<typename TImage>
// typename TImage::Pointer
// //void
// SegmentMRF(typename TImage::Pointer image ,
//        typename TImage::Pointer labelimage, unsigned int nclasses, float smf, unsigned int maxits)
// {
//
//   typedef TImage ImageType;
//   typedef typename TImage::PixelType PixelType;
//   enum { ImageDimension = ImageType::ImageDimension };
//   enum { Dimension = ImageType::ImageDimension };
//
//
//   //  const unsigned int NUMBANDS=1;
//   typedef itk::Image<itk::Vector<double,1>,ImageDimension> VecImageType;
//   typedef typename VecImageType::PixelType  VecPixelType;
//
//   /** copy the input image into this vector image.  stupid.  */
//
//   typename VecImageType::Pointer vecImage = VecImageType::New();
//   vecImage->SetSpacing(image->GetSpacing());
//   vecImage->SetOrigin(image->GetOrigin());
//   vecImage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   vecImage->SetBufferedRegion(  image->GetLargestPossibleRegion() );
//   vecImage->Allocate();
//   VecPixelType vvv;
//   vvv.Fill(0);
//   vecImage->FillBuffer(vvv);
//
//   // setup the iterators
//   //  typedef VecImageType::PixelType::VectorType VecPixelType;
//
//   enum { VecImageDimension = VecImageType::ImageDimension };
//   typedef itk::ImageRegionIterator< VecImageType > VecIterator;
//
//   VecIterator outIt( vecImage, vecImage->GetBufferedRegion() );
//   for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
//     {
//       vvv[0]=image->GetPixel(outIt.GetIndex());
//       outIt.Set(vvv);
//     }
//
//
//   namespace stat = itk::Statistics;
//   typedef itk::Image<PixelType,ImageDimension> ClassImageType;
//   typedef stat::MahalanobisDistanceMembershipFunction< VecPixelType >
//     MembershipFunctionType ;
//   typedef typename MembershipFunctionType::Pointer MembershipFunctionPointer ;
//   typedef std::vector< MembershipFunctionPointer >   MembershipFunctionPointerVector;
//
//   typedef itk::ImageGaussianModelEstimator<VecImageType,
//     MembershipFunctionType, ClassImageType>
//     ImageGaussianModelEstimatorType;
//
//   typename ImageGaussianModelEstimatorType::Pointer  applyEstimateModel =
//     ImageGaussianModelEstimatorType::New();
//   /*
//   typedef itk::Statistics::WeightedCentroidKdTreeGenerator<
//                                                       AdaptorType >
//                                                               TreeGeneratorType;
//   typedef typename TreeGeneratorType::KdTreeType TreeType;
//   typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
//
//   */
//   applyEstimateModel->SetNumberOfModels(nclasses);
//   applyEstimateModel->SetInputImage(vecImage);
//   applyEstimateModel->SetTrainingImage(labelimage);
//
//   //Run the gaussian classifier algorithm
//   applyEstimateModel->Update();
//   applyEstimateModel->Print(// std::cout);
//
//   MembershipFunctionPointerVector membershipFunctions =
//     applyEstimateModel->GetMembershipFunctions();
//
//   //----------------------------------------------------------------------
//   //Set the decision rule
//   //----------------------------------------------------------------------
//   typedef typename itk::DecisionRuleBase::Pointer DecisionRuleBasePointer;
//
//   typedef itk::MinimumDecisionRule DecisionRuleType;
//   typename DecisionRuleType::Pointer
//     myDecisionRule = DecisionRuleType::New();
//
//   //----------------------------------------------------------------------
//   // Set the classifier to be used and assigne the parameters for the
//   // supervised classifier algorithm except the input image which is
//   // grabbed from the MRF application pipeline.
//   //----------------------------------------------------------------------
//   //---------------------------------------------------------------------
//   typedef PixelType MeasurementVectorType;
//
//   typedef itk::ImageClassifierBase< VecImageType,
//     ClassImageType > ClassifierType;
//
//   typedef typename itk::ClassifierBase<VecImageType>::Pointer
//     ClassifierBasePointer;
//
//   typedef typename ClassifierType::Pointer ClassifierPointer;
//   ClassifierPointer myClassifier = ClassifierType::New();
//   // Set the Classifier parameters
//   myClassifier->SetNumberOfClasses(nclasses);
//
//   // Set the decison rule
//   myClassifier->
//     SetDecisionRule((DecisionRuleBasePointer) myDecisionRule );
//
//   //Add the membership functions
//   double meanDistance=0;
//   for( unsigned int i=0; i<nclasses; i++ )
//     {
//     myClassifier->AddMembershipFunction( membershipFunctions[i] );
//     meanDistance+=membershipFunctions[i]->GetMean()[0];
//     }
//   meanDistance/=(float)nclasses;
//   // std::cout << " mean dist " << meanDistance << std::endl;
//
//
//   //----------------------------------------------------------------------
//   // Set the MRF labeller and populate the parameters
//   //----------------------------------------------------------------------
//
//   //Set the MRF labeller
//   typedef itk::MRFImageFilter<VecImageType,ClassImageType> MRFImageFilterType;
//   typename MRFImageFilterType::Pointer applyMRFImageFilter = MRFImageFilterType::New();
//
//   // Set the MRF labeller parameters
//   applyMRFImageFilter->SetNumberOfClasses( nclasses );
//   applyMRFImageFilter->SetMaximumNumberOfIterations( maxits );
//   applyMRFImageFilter->SetErrorTolerance( 1.e-4 );
//   applyMRFImageFilter->SetSmoothingFactor( smf );
//   applyMRFImageFilter->SetInput(vecImage);
//   applyMRFImageFilter->SetClassifier( myClassifier );
//
//   //For setting up a square/cubic or hypercubic neighborhood
//   applyMRFImageFilter->SetNeighborhoodRadius( 1 );
//   std::vector<double> weights =
//     applyMRFImageFilter->GetMRFNeighborhoodWeight();
//   std::vector<double> testNewNeighborhoodWeight( weights.size(), 1);
//   double totalWeight = 0;
//   for(std::vector< double >::const_iterator wcIt = weights.begin();
//       wcIt != weights.end(); ++wcIt )
//     {
//     totalWeight += *wcIt;
//     }
//   unsigned int jj = 0;
//   for(std::vector< double >::iterator wIt = weights.begin();
//       wIt != weights.end(); wIt++ )
//     {
//      testNewNeighborhoodWeight[jj] = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
//      //   // std::cout << " ow " << weights[jj] << " nw " <<  testNewNeighborhoodWeight[jj] << std::endl;
//     jj++;
//     }
//
// //  applyMRFImageFilter->SetMRFNeighborhoodWeight( testNewNeighborhoodWeight );
//   //  applyMRFImageFilter->SetMRFNeighborhoodWeight( weights );
//
//   //Kick off the MRF labeller function
//   applyMRFImageFilter->Update();
//
//   applyMRFImageFilter->Print(// std::cout);
//   // std::cout << "Number of Iterations : " << applyMRFImageFilter->GetNumberOfIterations()
//     << std::endl;
//   // std::cout << "Stop condition: (1) Maximum number of iterations (2) Error tolerance:  "
//     << applyMRFImageFilter->GetStopCondition() << std::endl;
//
//   typename ClassImageType::Pointer  outClassImage = applyMRFImageFilter->GetOutput();
//
//   //Testing of different parameter access functions in the filter
//   // std::cout << "The number of classes labelled was: " <<
//     applyMRFImageFilter->GetNumberOfClasses() << std::endl;
//   // std::cout << "The maximum number of iterations were: " <<
//     applyMRFImageFilter->GetMaximumNumberOfIterations() << std::endl;
//   // std::cout << "The error tolerace threshold was: " <<
//     applyMRFImageFilter->GetErrorTolerance() << std::endl;
//   // std::cout << "The smoothing MRF parameter used was: " <<
//     applyMRFImageFilter->GetSmoothingFactor() << std::endl;
//   // std::cout << "The MRF neighborhood weights are: " << std::endl;
//
//
//   return  outClassImage;
//
// }

template <typename TImage>
typename TImage::Pointer
// void
itkMRIBiasFieldCorrectionFilter(typename TImage::Pointer image,
                                typename TImage::Pointer labelimage, unsigned int sd = 2)
{
  // std::cout << "doing Bias corr " << std::endl;

  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;

//    bool SaveImages = false;

  //  WriteImage<ImageType>(image,"temp1.nii");
  //  WriteImage<ImageType>(labelimage,"temp2.nii");

  // class statistics for two classes: a bright sphere and background
  // 3 or 1 classes ? SD
  // unsigned int numclasses=1;
  unsigned int              numclasses = 4;
  itk::Array<double>        classMeans(numclasses);
  itk::Array<double>        classSigmas(numclasses);
  itk::Array<unsigned long> classCounts(numclasses);
  classMeans.Fill(0.0);
  classSigmas.Fill(0.0);
  classCounts.Fill(0);

  // ******** get output mask from image SD
  typename ImageType::Pointer outputmask = BinaryThreshold<TImage>(0.1, 1.e9, 1, image);
  // ******** get output mask from image SD

  //  WriteImage<TImage>(mask,"outMK1.nii");
  // SD play with mask to see if we can move the artifact line
  //  mask=MorphologicalErosion<TImage>(1.5, mask);
  // mask=MorphologicalErosion<TImage>(3.5, mask);
  // WriteImage<TImage>(mask,"outMK2.nii");

  ImageIteratorType o_iter( labelimage, labelimage->GetLargestPossibleRegion() );
  o_iter.GoToBegin();
  while( !o_iter.IsAtEnd() )
    {
    unsigned int label = (unsigned int) o_iter.Get();
    float        pix = image->GetPixel(o_iter.GetIndex() );
    // uncomment for using one tissue class only for input mask SD
    /*
        if ( mask->GetPixel(o_iter.GetIndex()) == 1 )
     {
            label=label-1;
            unsigned int ind=0;
            classCounts[ind] +=1;
            float n = classCounts[ind];
            classMeans[ind] = (n-1.0)/( n )*  classMeans[ind] + 1.0 / n * pix;
            double sqrdf =  ( pix -  classMeans[ind])*( pix -  classMeans[ind]);
            if (n > 1)classSigmas[ind]=  (n-1)/n * classSigmas[ind] +  1.0 / (n-1) * sqrdf;
     }
     */

    // uncomment for using 3 tissue classes SD
    if( outputmask->GetPixel(o_iter.GetIndex() ) == 1 )
      {
      // label=label-1;
      classCounts[label] += 1;
      double n = static_cast<double>( classCounts[label] );
      classMeans[label] = (n - 1.0) / ( n ) *  classMeans[label] + 1.0 / n * static_cast<double>( pix );
      double sqrdf =  ( static_cast<double>( pix ) -  classMeans[label]) * ( static_cast<double>( pix ) -  classMeans[label]);
      if( n > 1 )
        {
        classSigmas[label] =  (n - 1.0) / n * classSigmas[label] +  1.0 / (n - 1.0) * sqrdf;
        }
      }

    ++o_iter;
    }
  for( unsigned int k = 0; k < numclasses; k++ )
    {
    classSigmas[k] = sqrt(classSigmas[k]);
    // std::cout << " Initial Means pre-bias " << classMeans[k] << " sig " <<  classSigmas[k] << std::endl;
    }

  // creats a normal random variate generator
  // itk::Statistics::NormalVariateGenerator::Pointer randomGenerator =
  //  itk::Statistics::NormalVariateGenerator::New() ;

  // creates a bias correction filter and run it.
  typedef itk::MRIBiasFieldCorrectionFilter<ImageType, ImageType, ImageType> FilterType;

  // std::cout << "before new filter" << std::endl;
  typename FilterType::Pointer filter = FilterType::New();
  // std::cout << "after new filter" << std::endl;

  //  typename FilterType::BiasFieldType::CoefficientArrayType

  filter->SetInput( image.GetPointer() );
  // filter->SetInput( image ) ;
  filter->IsBiasFieldMultiplicative( true );    // correct with multiplicative bias
  unsigned int biasDegree = 4;
  filter->SetBiasFieldDegree( biasDegree );    // default value = 3
  filter->SetTissueClassStatistics( classMeans, classSigmas );
  // filter->SetOptimizerGrowthFactor( 1.01 ) ; // default value
  // filter->SetOptimizerInitialRadius( 0.02 ) ; // default value

  // SD debug don't do interslice correction
  // filter->SetUsingInterSliceIntensityCorrection( true ) ; // default value
  filter->SetUsingInterSliceIntensityCorrection( false );    // default value

  filter->SetVolumeCorrectionMaximumIteration( 200);      // default value = 100
  filter->SetInterSliceCorrectionMaximumIteration( 100 ); // default value = 100
  filter->SetUsingSlabIdentification( false );            // default value = false
                                                          // filter->SetSlabBackgroundMinimumThreshold( 0 ) ; // default
                                                          // value
                                                          // filter->SetSlabNumberOfSamples( 10 ) ; // default value
                                                          // filter->SetSlabTolerance(0.0) ; // default value
  filter->SetSlicingDirection(sd);                        // default value
  filter->SetUsingBiasFieldCorrection( true );            // default value
  filter->SetGeneratingOutput( true );                    // default value

  filter->SetInputMask( outputmask );

  // ******** Try different output mask SD
  //  filter->SetOutputMask( labelimage ) ;
  filter->SetOutputMask( outputmask );
  //  filter->SetDebug( true );
  //    WriteImage<ImageType>(mask,"mask.nii");
  //    WriteImage<ImageType>(outputmask,"outputmask.nii");
  //    WriteImage<ImageType>(image,"image.nii");
  //    WriteImage<ImageType>(labelimage,"labelimage.nii");
  // ******** Try different output mask SD

  // default schedule is 2 2 2 - 1 1 1, let's change this
  bool setnewsched = false;
  if( setnewsched )
    {
    unsigned int                      nlev = 1;
    typename FilterType::ScheduleType schedule( nlev, ImageDimension);
    schedule.Fill( 4 );
    //      for (unsigned int jj=0; jj<ImageDimension; jj++) schedule[0][jj]=4;
    // for (unsigned int jj=0; jj<ImageDimension; jj++) schedule[1][jj]=2;
    // for (unsigned int jj=0; jj<ImageDimension; jj++) schedule[2][jj]=1;
    filter->SetNumberOfLevels( nlev );      // Important to set this first, otherwise the filter rejects the new
                                            // schedule
    filter->SetSchedule( schedule );
    }
  //  filter->SetInitialBiasFieldCoefficients(initCoefficients);
  filter->SetVolumeCorrectionMaximumIteration( 200 );     // default value = 100
  filter->SetInterSliceCorrectionMaximumIteration( 100 ); // default value = 100
                                                          // filter->SetOptimizerInitialRadius( 0.02 ) ; // default
                                                          // value
                                                          // timing
  // long int t1 = time(nullptr);
  filter->Update();
  // long int t2 = time(nullptr);
  // std::cout << "Run time (in s)" << t2 - t1  << std::endl;

  return filter->GetOutput();
}

// template<typename TImage>
// typename TImage::Pointer
// //void
// SegmentMRFKM(typename TImage::Pointer image ,
//          typename TImage::Pointer labelimage, unsigned int nclasses, float smf, unsigned int maxit)
// {
//
//   typedef TImage ImageType;
//   typedef typename TImage::PixelType PixelType;
//   enum { ImageDimension = ImageType::ImageDimension };
//   enum { Dimension = ImageType::ImageDimension };
//
//
//   //  const unsigned int NUMBANDS=1;
//   typedef itk::Image<itk::Vector<double,1>,ImageDimension> VecImageType;
//   typedef typename VecImageType::PixelType  VecPixelType;
//
//   /** copy the input image into this vector image.  stupid.  */
//
//   typename VecImageType::Pointer vecImage = VecImageType::New();
//   vecImage->SetSpacing(image->GetSpacing());
//   vecImage->SetOrigin(image->GetOrigin());
//   vecImage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//   vecImage->SetBufferedRegion(  image->GetLargestPossibleRegion() );
//   vecImage->Allocate();
//   VecPixelType vvv;
//   vvv.Fill(0);
//   vecImage->FillBuffer(vvv);
//
//   // setup the iterators
//   //  typedef VecImageType::PixelType::VectorType VecPixelType;
//
//   enum { VecImageDimension = VecImageType::ImageDimension };
//   typedef itk::ImageRegionIterator< VecImageType > VecIterator;
//
//   VecIterator outIt( vecImage, vecImage->GetBufferedRegion() );
//   for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
//     {
//       vvv[0]=image->GetPixel(outIt.GetIndex());
//       outIt.Set(vvv);
//     }
//
//
//   namespace stat = itk::Statistics;
//
//   typedef itk::Image<PixelType,ImageDimension> ClassImageType;
//
//
//   //----------------------------------------------------------------------
//   //Set membership function (Using the statistics objects)
//   //----------------------------------------------------------------------
//
//   typedef itk::Statistics::DistanceToCentroidMembershipFunction< VecPixelType >
//     MembershipFunctionType ;
//
//   typedef typename MembershipFunctionType::Pointer MembershipFunctionPointer ;
//
//   typedef std::vector< MembershipFunctionPointer >
//     MembershipFunctionPointerVector;
//
//   //----------------------------------------------------------------------
//   //Set the image model estimator
//   //----------------------------------------------------------------------
//   typedef itk::ImageKmeansModelEstimator< VecImageType,
//     MembershipFunctionType> ImageKmeansModelEstimatorType;
//
//   typename ImageKmeansModelEstimatorType::Pointer
//     applyKmeansModelEstimator = ImageKmeansModelEstimatorType::New();
//
//   //----------------------------------------------------------------------
//   //Set the parameters of the clusterer
//   //----------------------------------------------------------------------
//
//   // std::cout << "Starting to build the K-means model ....." << std::endl;
//
//   applyKmeansModelEstimator->SetInputImage( vecImage );
//   applyKmeansModelEstimator->SetNumberOfModels(nclasses);
//   applyKmeansModelEstimator->SetThreshold(0.0001);
//   applyKmeansModelEstimator->Update();
//
//   MembershipFunctionPointerVector membershipFunctions =
//     applyKmeansModelEstimator->GetMembershipFunctions();
//
//   typedef std::vector<double> TempVectorType;
//   typedef TempVectorType::iterator TempVectorIterator;
//   TempVectorIterator  start, end;
//
//   std::vector<double> kmeansResultForClass(membershipFunctions.size());
//
//
//   // std::cout << "Result of K-Means clustering" << std::endl;
//
//   double meanDistance=0;
//   for(unsigned int classIndex=0; classIndex < membershipFunctions.size();
//     classIndex++ )
//     {
//     kmeansResultForClass[classIndex] =
//       (double) (membershipFunctions[classIndex]->GetCentroid())[0];
//     meanDistance+=kmeansResultForClass[classIndex];//membershipFunctions[i]->GetMean()[0];
//     }
//   meanDistance/=(float)nclasses;
//   // std::cout << " mean dist " << meanDistance << std::endl;
//
//
//   start = kmeansResultForClass.begin();
//   end   = kmeansResultForClass.end();
//
//   std::sort( start, end );
//
//   vnl_vector<double> temp =  membershipFunctions[0]->GetCentroid();
//   for(unsigned int classIndex=0; classIndex < membershipFunctions.size();
//     classIndex++ )
//     {
//     temp[0] = (double) kmeansResultForClass[classIndex];
//     membershipFunctions[classIndex]->SetCentroid(temp);
//     }
//
//   for(unsigned int classIndex=0; classIndex < membershipFunctions.size();
//     classIndex++ )
//     {
//     // std::cout <<  (membershipFunctions[classIndex]->GetCentroid())[0] << std::endl;
//     }
//
//   //----------------------------------------------------------------------
//   //Set the decision rule
//   //----------------------------------------------------------------------
//   typedef itk::DecisionRuleBase::Pointer DecisionRuleBasePointer;
//
//   typedef itk::MinimumDecisionRule DecisionRuleType;
//   DecisionRuleType::Pointer
//     classifierDecisionRule = DecisionRuleType::New();
//
//   //------------------------------------------------------
//   //Instantiate the classifier model (as the input image is in right format)
//   //------------------------------------------------------
//
//   //Assign a class label image type
//
//   typedef itk::Image<PixelType,ImageDimension> ClassImageType;
//
//   typedef itk::ImageClassifierBase< VecImageType,ClassImageType >
//     SupervisedClassifierType;
//
//   typename SupervisedClassifierType::Pointer
//     classifierPointer = SupervisedClassifierType::New();
//
//
//   //------------------------------------------------------
//   // Set the Classifier parameters
//   //------------------------------------------------------
//   classifierPointer->SetNumberOfClasses( nclasses );
//   classifierPointer->SetInputImage( vecImage );
//
//   // Set the decison rule
//   classifierPointer->
//     SetDecisionRule( (DecisionRuleBasePointer) classifierDecisionRule );
//
//   MembershipFunctionPointer membershipFunction;
//   //------------------------------------------------------
//   //Set the classifier membership functions
//   //------------------------------------------------------
//   for( unsigned int i=0; i<nclasses; i++ )
//     {
//     classifierPointer->AddMembershipFunction( membershipFunctions[i] );
//     }
//
//   //Do the classification
//   //Run the kmeans classifier algorithm
//   classifierPointer->Update();
//
//   //Get the classified image
//   typedef typename ClassImageType::Pointer ClassifiedImagePointer;
//   ClassifiedImagePointer outClassImage =
//     classifierPointer->GetClassifiedImage();
//
//   //------------------------------------------------------
//   //Mask the output of the classifier
//   //------------------------------------------------------
//
//   // Declare the type for the MaskInput filter
//
//   typedef itk::MaskImageFilter< ClassImageType,
//                            ClassImageType,
//                            ClassImageType  >   MaskFilterType;
//
//   typedef typename ClassImageType::Pointer   MaskedOutputImagePointer;
//   typedef typename MaskFilterType::Pointer        MaskFilterTypePointer;
//
//   // Create an ADD Filter
//   MaskFilterTypePointer maskfilter = MaskFilterType::New();
//
//   // Connect the input images
//   maskfilter->SetInput1( outClassImage );
//   maskfilter->SetInput2( labelimage );
//
//   // Execute the filter
//   maskfilter->Update();
//
//   // Get the Smart Pointer to the Filter Output
//   MaskedOutputImagePointer maskedOutputImage = maskfilter->GetOutput();
//
//   //  this->SetClassifiedImage( maskedOutputImage );
//
//   //------------------------------------------------------
//   //Set the MRF labeller and populate the parameters
//   //------------------------------------------------------
//   //Set the MRF labeller
//   typedef itk::MRFImageFilter<VecImageType,ClassImageType>
//     MRFFilterType;
//
//   typename MRFFilterType::Pointer applyMRFFilter = MRFFilterType::New();
//
//   // Set the MRF labeller parameters
//   applyMRFFilter->SetNumberOfClasses(nclasses);
//   unsigned int m_MaximumNumberOfIterations=maxit;
//   applyMRFFilter->SetMaximumNumberOfIterations(m_MaximumNumberOfIterations);
//   float m_ErrorTolerance=1.e-5;
//   applyMRFFilter->SetErrorTolerance(m_ErrorTolerance);
//   float m_SmoothingFactor=smf;
//   applyMRFFilter->SetSmoothingFactor( m_SmoothingFactor );
//
//   //For setting up a square/cubic or hypercubic neighborhood
//   applyMRFFilter->SetNeighborhoodRadius( 1 );
//   std::vector<double> weights =
//     applyMRFFilter->GetMRFNeighborhoodWeight();
//   std::vector<double> testNewNeighborhoodWeight( weights.size(), 1);
//   double totalWeight = 0;
//   for(std::vector< double >::const_iterator wcIt = weights.begin();
//       wcIt != weights.end(); ++wcIt )
//     {
//     totalWeight += *wcIt;
//     }
//   unsigned int jj = 0;
//   for(std::vector< double >::iterator wIt = weights.begin();
//       wIt != weights.end(); wIt++ )
//     {
//      testNewNeighborhoodWeight[jj] = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
//      //// std::cout << " ow " << weights[jj] << " nw " <<  testNewNeighborhoodWeight[jj] << std::endl;
//     jj++;
//     }
//
//   applyMRFFilter->SetMRFNeighborhoodWeight( testNewNeighborhoodWeight );
//
//   applyMRFFilter->SetInput(vecImage);
//   applyMRFFilter->SetClassifier( classifierPointer );
//
//   //Kick off the MRF labeller function
//   applyMRFFilter->Update();
//
//   applyMRFFilter->Print(// std::cout);
//   outClassImage = applyMRFFilter->GetOutput();
//
//   //------------------------------------------------------
//   //Mask the output of the classifier
//   //------------------------------------------------------
//
//   // Declare the type for the MaskInput filter
//
//   // Create an ADD Filter
//   MaskFilterTypePointer maskfilter2 = MaskFilterType::New();
//
//   // Connect the input images
//   maskfilter2->SetInput1( outClassImage );
//   maskfilter2->SetInput2( labelimage );
//
//   // Execute the filter
//   maskfilter2->Update();
//
//   // Get the Smart Pointer to the Filter Output
//   maskedOutputImage = maskfilter2->GetOutput();
//
//   //  this->SetClassifiedImage( maskedOutputImage );
//
//   return outClassImage;//maskedOutputImage;
//
// }

template <unsigned int ImageDimension>
int SmoothImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  std::vector<float> sigmaVector;
  if( argc > argct )
    {
    sigmaVector = ConvertVector<float>( argv[argct] );
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer varimage = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();

  if( sigmaVector.size() == 1 )
    {
    filter->SetVariance( itk::Math::sqr ( sigmaVector[0] ) );
    }
  else if( sigmaVector.size() == ImageDimension )
    {
    typename dgf::ArrayType varianceArray;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      varianceArray[d] = itk::Math::sqr ( sigmaVector[d] );
      }
    filter->SetVariance( varianceArray );
    }
  else
    {
    std::cerr << "Incorrect sigma vector size.  Must either be of size 1 or ImageDimension." << std::endl;
    }
  bool usespacing = true;
  if( !usespacing )
    {
    filter->SetUseImageSpacingOff();
    }
  else
    {
    filter->SetUseImageSpacingOn();
    }
  filter->SetMaximumError(.01f);
  filter->SetInput(image1);
  filter->Update();
  varimage = filter->GetOutput();
  WriteImage<ImageType>( varimage, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int MorphImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string       operation = std::string(argv[argct]);  argct++;
  std::string       fn1 = std::string(argv[argct]);   argct++;
  float             sigma = 1.0;
  if( argc > argct )
    {
    sigma = atof(argv[argct]);
    }

  unsigned int morphopt = 1;
  if( strcmp(operation.c_str(), "ME") == 0 )
    {
    morphopt = 0;
    }
  else if( strcmp(operation.c_str(), "MO") == 0 )
    {
    morphopt = 2;
    }
  else if( strcmp(operation.c_str(), "MC") == 0 )
    {
    morphopt = 3;
    }
  else if( strcmp(operation.c_str(), "GE") == 0 )
    {
    morphopt = 4;
    }
  else if( strcmp(operation.c_str(), "GD") == 0 )
    {
    morphopt = 5;
    }
  else if( strcmp(operation.c_str(), "GO") == 0 )
    {
    morphopt = 6;
    }
  else if( strcmp(operation.c_str(), "GC") == 0 )
    {
    morphopt = 7;
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );

  // SRD set dilatevalue
  float dilateval = 1;
  if( argc > argct + 1 )
    {
    dilateval = atof(argv[argct + 1]);
    }

  image1 = ants::Morphological<ImageType>(image1, sigma, morphopt, dilateval);

  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(image1, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int FastMarchingSegmentation( unsigned int argc, char *argv[] )
{
  unsigned int      argct = 2;
  const std::string outname = std::string(argv[argct]);
  std::cout << outname << argc << " This function has been disabled --- see PropagateLabelsThroughMask " << std::endl;
  return 0;
}


template <unsigned int ImageDimension>
int PropagateLabelsThroughMask(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  float             thresh = 0.5;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  float       stopval = 100.0;
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  else
    {
    // std::cout << " not enough parameters -- need label image " << std::endl;  return 0;
    }
  if(  argc > argct )
    {
    stopval = atof(argv[argct]);   argct++;
    }
  unsigned int topocheck = 0;
  if(  argc > argct )
    {
    topocheck = std::stoi(argv[argct]);   argct++;
    }

  typename ImageType::Pointer speedimage = nullptr;
  ReadImage<ImageType>(speedimage, fn1.c_str() );
  typename ImageType::Pointer labimage = nullptr;
  ReadImage<ImageType>(labimage, fn2.c_str() );
  typename ImageType::Pointer fastimage = nullptr;
  ReadImage<ImageType>(fastimage, fn1.c_str() );
  typename ImageType::Pointer outlabimage = nullptr;
  ReadImage<ImageType>(outlabimage, fn2.c_str() );
  fastimage->FillBuffer(1.e9);
  // compute max label
  double   maxlabel = 0;
  Iterator vfIter2( labimage,  labimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool   isinside = true;
    double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
    double labval = labimage->GetPixel(vfIter2.GetIndex() );
    if( speedval < static_cast<double>( thresh ) )
      {
      isinside = false;
      }
    if( isinside )
      {
      if( labval > maxlabel )
        {
        maxlabel = labval;
        }
      }
    }
  typedef  itk::FMarchingImageFilter<ImageType, ImageType> FastMarchingFilterType;
  typedef  typename FastMarchingFilterType::LabelImageType LabelImageType;
  typename FastMarchingFilterType::Pointer  fastMarching;
  for( unsigned int lab = 1; lab <= (unsigned int)maxlabel; lab++ )
    {

    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( labimage );
    thresholder->SetLowerThreshold( lab );
    thresholder->SetUpperThreshold( lab );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::LabelContourImageFilter<ImageType, ImageType>
      ContourFilterType;
    typename ContourFilterType::Pointer contour = ContourFilterType::New();
    contour->SetInput( thresholder->GetOutput() );
    contour->FullyConnectedOff();
    contour->SetBackgroundValue( itk::NumericTraits<typename LabelImageType::PixelType>::ZeroValue() );
    contour->Update();
    typename ImageType::Pointer contourimage = contour->GetOutput();

    fastMarching = FastMarchingFilterType::New();
    fastMarching->SetInput( speedimage );
    fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
    if( topocheck == 1 )  // Strict
      {
      // std::cout << " strict " << std::endl;
      fastMarching->SetTopologyCheck( FastMarchingFilterType::Strict );
      }
    if( topocheck == 2 )  // No handles
      {
      // std::cout << " no handles " << std::endl;
      fastMarching->SetTopologyCheck( FastMarchingFilterType::NoHandles );
      }
    typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
    typedef typename FastMarchingFilterType::NodeType      NodeType;
    typename NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();
    typename NodeContainer::Pointer alivePoints = NodeContainer::New();
    alivePoints->Initialize();
    unsigned long aliveCount = 0;
    unsigned long ct = 0;
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      double labval = labimage->GetPixel( vfIter2.GetIndex() );
      double contourval = contourimage->GetPixel( vfIter2.GetIndex() );
      if( ( (unsigned int) contourval == 1 ) && ( (unsigned int) labval == lab ) )
        {
	NodeType     node;
	const double seedValue = 0.0;
	node.SetValue( seedValue );
	node.SetIndex( vfIter2.GetIndex() );
	seeds->InsertElement( ct, node );
	ct++;
        }
      if( ( (unsigned int) contourval == 0 ) && ( (unsigned int) labval == lab ) )
        {
	NodeType     node;
	const double seedValue = 0.0;
	node.SetValue( seedValue );
	node.SetIndex( vfIter2.GetIndex() );
	alivePoints->InsertElement( aliveCount, node );
	aliveCount++;
        }
      }
    fastMarching->SetTrialPoints(  seeds  );
    fastMarching->SetAlivePoints( alivePoints );
    fastMarching->SetStoppingValue(  stopval );
    fastMarching->Update();
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      bool   isinside = true;
      double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
      double labval = labimage->GetPixel(vfIter2.GetIndex() );
      if( speedval < static_cast<double>( thresh ) )
        {
        isinside = false;
        }
      if( isinside && itk::Math::FloatAlmostEqual( labval, 0.0 ) )
        {
        double fmarrivaltime = fastMarching->GetOutput()->GetPixel( vfIter2.GetIndex() );
        double mmm = fastimage->GetPixel(vfIter2.GetIndex() );
        if( fmarrivaltime < mmm )
          {
          fastimage->SetPixel(vfIter2.GetIndex(),  fmarrivaltime );
          outlabimage->SetPixel(vfIter2.GetIndex(), lab );
          }
        }
      else if( !isinside )
        {
        outlabimage->SetPixel(vfIter2.GetIndex(), 0 );
        }
      }
    }
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );
  std::string kname = tempname + std::string("_speed") + extension;
  std::string lname = tempname + std::string("_label") + extension;
  WriteImage<ImageType>(fastimage, kname.c_str() );
  WriteImage<ImageType>(outlabimage, outname.c_str() );

// this nonsense fixes a type error
  typedef itk::CastImageFilter<LabelImageType, LabelImageType>                 CastFilterType;
  typename CastFilterType::Pointer castRegions = CastFilterType::New();
  castRegions->SetInput( fastMarching->GetLabelImage() );
  WriteImage<LabelImageType>( castRegions->GetOutput(), lname.c_str() );
  return 0;
}



template <unsigned int ImageDimension>
int FastMarchingExtension(int argc, char *argv[])
{
  typedef float                                     PixelType;
  typedef itk::Image<PixelType, ImageDimension>     ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++; // speed image
  std::string fn2 = ""; // label image
  std::string fn3 = ""; // value image
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  else
    {
    // std::cout << " not enough parameters -- need label image " << std::endl;
    return 0;
    }
  if(  argc > argct )
    {
    fn3 = std::string(argv[argct]);   argct++;
    }
  else
    {
    // std::cout << " not enough parameters -- need value image " << std::endl;
    return 0;
    }
  typename ImageType::Pointer speedimage = nullptr;
  ReadImage<ImageType>(speedimage, fn1.c_str() );
  typename ImageType::Pointer labimage = nullptr;
  ReadImage<ImageType>(labimage, fn2.c_str() );
  typename ImageType::Pointer valimage = nullptr;
  ReadImage<ImageType>(valimage, fn3.c_str() );
  typedef itk::FastMarchingThresholdStoppingCriterion< ImageType, ImageType >
    CriterionType;
  typedef typename CriterionType::Pointer CriterionPointer;
  CriterionPointer criterion = CriterionType::New();
  criterion->SetThreshold( 1.e9 ); // something large
  typedef  itk::FastMarchingExtensionImageFilterBase<ImageType, ImageType,
    PixelType,1>  MarcherBaseType;
  typedef  typename MarcherBaseType::LabelImageType LabelImageType;
  typename MarcherBaseType::Pointer  fastMarching;
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( labimage );
  thresholder->SetLowerThreshold( 0.5 );
  thresholder->SetUpperThreshold( 1.001 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::LabelContourImageFilter<ImageType, ImageType> ContourFilterType;
  typename ContourFilterType::Pointer contour = ContourFilterType::New();
  contour->SetInput( thresholder->GetOutput() );
  contour->FullyConnectedOff();
  contour->SetBackgroundValue( itk::NumericTraits<typename LabelImageType::PixelType>::ZeroValue() );
  contour->Update();
  typename ImageType::Pointer contourimage = contour->GetOutput();
  // contour defines starting points

  fastMarching = MarcherBaseType::New();
  fastMarching->SetInput( speedimage );
  typedef typename MarcherBaseType::NodePairType           NodePairType;
  typedef typename MarcherBaseType::NodePairContainerType  NodePairContainerType;
  typedef typename MarcherBaseType::AuxValueVectorType     AuxValueVectorType;
  typedef typename MarcherBaseType::AuxValueContainerType  AuxValueContainerType;
  typename AuxValueContainerType::Pointer auxAliveValues = AuxValueContainerType::New();
  typename AuxValueContainerType::Pointer auxTrialValues = AuxValueContainerType::New();
  typename NodePairContainerType::Pointer seeds = NodePairContainerType::New();
  seeds->Initialize();
  typename NodePairContainerType::Pointer alivePoints = NodePairContainerType::New();
  alivePoints->Initialize();
  unsigned int seedct = 0, alivect = 0;
  Iterator vfIter2( labimage,  labimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    typename ImageType::IndexType ind = vfIter2.GetIndex();
    double labval = labimage->GetPixel( ind );
    double contourval = contourimage->GetPixel( ind );
    if ( ( (unsigned int) contourval == 1 )  )
      {
      seeds->push_back( NodePairType(  ind, 0. ) );
      AuxValueVectorType vector;
      vector[0] = valimage->GetPixel( ind  );
      auxTrialValues->push_back( vector );
      seedct++;
      }
    if (  ( labval > 0  ) && ((unsigned int) contourval != 1 ) )
      {
      alivePoints->push_back( NodePairType(  ind, 0. ) );
      AuxValueVectorType vector;
      vector[0] = valimage->GetPixel( ind  );
      auxAliveValues->push_back( vector );
      alivect++;
      }
  }
  fastMarching->SetTrialPoints(  seeds  );
  fastMarching->SetAuxiliaryTrialValues( auxTrialValues );
  fastMarching->SetAlivePoints( alivePoints );
  fastMarching->SetAuxiliaryAliveValues( auxAliveValues );
  fastMarching->SetStoppingCriterion( criterion );
  fastMarching->Update();
  WriteImage<ImageType>( fastMarching->GetAuxiliaryImage(0), outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int itkPropagateLabelsThroughMask(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  float             thresh = 1.e-9;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  float       stopval = 100.0;
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  else
    {
    // std::cout << " not enough parameters -- need label image " << std::endl;  return 0;
    }
  if(  argc > argct )
    {
    stopval = atof(argv[argct]);   argct++;
    }
  unsigned int topocheck = 0;
  if(  argc > argct )
    {
    topocheck = std::stoi(argv[argct]);   argct++;
    }

  typename ImageType::Pointer speedimage = nullptr;
  ReadImage<ImageType>(speedimage, fn1.c_str() );
  typename ImageType::Pointer labimage = nullptr;
  ReadImage<ImageType>(labimage, fn2.c_str() );
  typename ImageType::Pointer fastimage = nullptr;
  ReadImage<ImageType>(fastimage, fn1.c_str() );
  typename ImageType::Pointer outlabimage = nullptr;
  ReadImage<ImageType>(outlabimage, fn2.c_str() );
  fastimage->FillBuffer(1.e9);
  // compute max label
  double   maxlabel = 0;
  Iterator vfIter2( labimage,  labimage->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    bool   isinside = true;
    double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
    double labval = labimage->GetPixel(vfIter2.GetIndex() );
    if( speedval < static_cast<double>( thresh ) )
      {
      isinside = false;
      }
    if( isinside )
      {
      if( labval > maxlabel )
        {
        maxlabel = labval;
        }
      }
    }
  typedef itk::FastMarchingThresholdStoppingCriterion< ImageType, ImageType >
      CriterionType;
  typedef typename CriterionType::Pointer CriterionPointer;
  CriterionPointer criterion = CriterionType::New();
  criterion->SetThreshold( stopval );
  typedef  itk::FastMarchingImageFilterBase<ImageType, ImageType> FastMarchingFilterType;
  typedef  typename FastMarchingFilterType::LabelImageType LabelImageType;
  typename FastMarchingFilterType::Pointer  fastMarching;
  for( unsigned int lab = 1; lab <= (unsigned int)maxlabel; lab++ )
    {

    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( labimage );
    thresholder->SetLowerThreshold( lab );
    thresholder->SetUpperThreshold( lab );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::LabelContourImageFilter<ImageType, ImageType>
      ContourFilterType;
    typename ContourFilterType::Pointer contour = ContourFilterType::New();
    contour->SetInput( thresholder->GetOutput() );
    contour->FullyConnectedOff();
    contour->SetBackgroundValue( itk::NumericTraits<typename LabelImageType::PixelType>::ZeroValue() );
    contour->Update();
    typename ImageType::Pointer contourimage = contour->GetOutput();

    fastMarching = FastMarchingFilterType::New();
    fastMarching->SetInput( speedimage );
    fastMarching->SetStoppingCriterion( criterion );

    if( topocheck == 1 )  // Strict
      {
      // std::cout << " strict " << std::endl;
      fastMarching->SetTopologyCheck( FastMarchingFilterType::Strict );
      }
    if( topocheck == 2 )  // No handles
      {
      // std::cout << " no handles " << std::endl;
      fastMarching->SetTopologyCheck( FastMarchingFilterType::NoHandles );
      }
    typedef typename FastMarchingFilterType::NodePairContainerType NodeContainer;
    typedef typename FastMarchingFilterType::NodePairType      NodePairType;
    typename NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();
    typename NodeContainer::Pointer alivePoints = NodeContainer::New();
    alivePoints->Initialize();
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      double labval = labimage->GetPixel( vfIter2.GetIndex() );
      double contourval = contourimage->GetPixel( vfIter2.GetIndex() );
      if( ( (unsigned int) contourval == 1 ) && ( (unsigned int) labval == lab ) )
        {
	seeds->push_back( NodePairType(  vfIter2.GetIndex(), 0. ) );
        }
      if( ( (unsigned int) contourval == 0 ) && ( (unsigned int) labval == lab ) )
        {
	alivePoints->push_back( NodePairType(  vfIter2.GetIndex(), 0. ) );
        }
      }
    fastMarching->SetTrialPoints(  seeds  );
    fastMarching->SetAlivePoints( alivePoints );
    fastMarching->Update();
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      bool   isinside = true;
      double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
      double labval = labimage->GetPixel(vfIter2.GetIndex() );
      if( speedval < static_cast<double>( thresh ) )
        {
        isinside = false;
        }
      if( isinside && itk::Math::FloatAlmostEqual( labval, 0.0 ) )
        {
        double fmarrivaltime = fastMarching->GetOutput()->GetPixel( vfIter2.GetIndex() );
        double mmm = fastimage->GetPixel(vfIter2.GetIndex() );
        if( fmarrivaltime < mmm )
          {
          fastimage->SetPixel(vfIter2.GetIndex(),  fmarrivaltime );
          outlabimage->SetPixel(vfIter2.GetIndex(), lab );
          }
        }
      else if( !isinside )
        {
        outlabimage->SetPixel(vfIter2.GetIndex(), 0 );
        }
      }
    }

  WriteImage<ImageType>(outlabimage, outname.c_str() );
/*
  // Write debug images, but string manipulation breaks on any path containing periods

  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );
  std::string kname = tempname + std::string("_speed") + extension;
  std::string lname = tempname + std::string("_label") + extension;
  WriteImage<ImageType>(fastimage, kname.c_str() );

  // this nonsense fixes a type error
  typedef itk::CastImageFilter<LabelImageType, LabelImageType>                 CastFilterType;
  typename CastFilterType::Pointer castRegions = CastFilterType::New();
  castRegions->SetInput( fastMarching->GetLabelImage() );
  WriteImage<LabelImageType>( castRegions->GetOutput(), lname.c_str() );
*/

  return 0;
}





template <unsigned int ImageDimension>
int DistanceMap(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int argct = 2;
  if( argc < 5 )
    {
    // std::cout << "Missing required arguments ( output name, operation & fn1)" << std::endl;
    throw;
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );

  typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;

  typename  FilterType::Pointer filter = FilterType::New();
  filter->InputIsBinaryOff();
  filter->SetUseImageSpacing(true);
  filter->SetInput(image1);
  filter->Update();

  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(filter->GetOutput(), outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int GenerateMaurerDistanceImage( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  // Usage:  ImageMath 3 output.nii.gz MaurerDistance input.nii.gz {foreground=1}

  const std::string outputName = std::string( argv[2] );
  const std::string inputName = std::string( argv[4] );

  typename ImageType::Pointer input = nullptr;
  ReadImage<ImageType>( input, inputName.c_str() );

  PixelType foreground = 1;
  if( argc > 5 )
    {
    foreground = static_cast<PixelType>( atof( argv[5] ) );
    }

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( input );
  thresholder->SetLowerThreshold( foreground );
  thresholder->SetUpperThreshold( foreground );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( thresholder->GetOutput() );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false );
  filter->Update();

  WriteImage<ImageType>( filter->GetOutput(), outputName.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int FillHoles(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<int, ImageDimension>                                 LabelImageType;
  typedef itk::CastImageFilter<ImageType, LabelImageType>                 CastFilterType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       holeparam = 2.0;
  if( argc > argct )
    {
    holeparam = atof(argv[argct]);
    }

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::Pointer imageout = nullptr;
  ReadImage<ImageType>(imageout, fn1.c_str() );
  image1 = BinaryThreshold<ImageType>(0.5, 1.e9, 1, image1);

  typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;
  typename  FilterType::Pointer filter = FilterType::New();
  filter->InputIsBinaryOff();
  filter->SetUseImageSpacing(false);
  filter->SetInput(image1);
  filter->Update();
  // algorithm :
  // 1. get distance map of object
  // 2. threshold
  // 3. label connected components
  // 4. label surface
  // 5. if everywhere on surface is next to object then it's a hole
  // 6. make sure it's not the background
  typedef itk::Image<int, ImageDimension> labelimagetype;
  typename ImageType::Pointer dist = filter->GetOutput();
  typename ImageType::Pointer regions = BinaryThreshold<ImageType>(0.001, 1.e9, 1, dist);

  typedef itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype> ccFilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, ImageType>        RelabelType;
  typename RelabelType::Pointer relabel = RelabelType::New();
  typename ccFilterType::Pointer ccfilter = ccFilterType::New();

  typename CastFilterType::Pointer castRegions = CastFilterType::New();
  castRegions->SetInput( regions );

  ccfilter->SetInput(castRegions->GetOutput() );
  ccfilter->SetFullyConnected( 0 );
  ccfilter->Update();
  relabel->SetInput( ccfilter->GetOutput() );

  relabel->SetMinimumObjectSize( 0 );
  //    relabel->SetUseHistograms(true);
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & /* excep */ )
    {
    // std::cout << "Relabel: exception caught !" << std::endl;
    // std::cout << excep << std::endl;
    }

  // WriteImage<ImageType>(relabel->GetOutput(),"test.nii");

  if( itk::Math::FloatAlmostEqual( holeparam, 2.0f ) )
    {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> RelabelIterator;
    RelabelIterator vfIter( relabel->GetOutput(),
                            relabel->GetOutput()->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      if( vfIter.Get() > 1 )
        {
        imageout->SetPixel(vfIter.GetIndex(), 1);
        }
      }

    WriteImage<ImageType>(imageout, outname.c_str() );

    return 0;
    }

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  // now we have the exact number of objects labeled independently
  for( unsigned int lab = 2; lab <= maximum; lab++ )
    {
    float erat = 2;
    if( holeparam  <= 1 )
      {
      GHood.GoToBegin();

      unsigned long objectedge = 0;
      unsigned long backgroundedge = 0;
      unsigned long totaledge = 0;
      unsigned long volume = 0;

      while( !GHood.IsAtEnd() )
        {
        typename ImageType::PixelType p = GHood.GetCenterPixel();
        typename ImageType::IndexType ind2;
        if( itk::Math::FloatAlmostEqual( p, static_cast<float>( lab ) ) )
          {
          volume++;
          for( unsigned int i = 0; i < GHood.Size(); i++ )
            {
            ind2 = GHood.GetIndex(i);
            float val2 = image1->GetPixel(ind2);
            if( val2 >= 0.5f && ! itk::Math::FloatAlmostEqual( GHood.GetPixel(i), static_cast<float>( lab ) ) )
              {
              objectedge++;
              totaledge++;
              }
            else if( val2 < 0.5f && ! itk::Math::FloatAlmostEqual( GHood.GetPixel(i), static_cast<float>( lab ) ) )
              {
              backgroundedge++;
              totaledge++;
              }
            }
          }
        ++GHood;
        }

      erat = (float)objectedge / (float)totaledge;
      // std::cout << " Lab " << lab << " volume " << volume << " v-rat " << vrat << " edge " << erat << std::endl;
      }

    if( erat > holeparam ) // fill the hole
      {
      // std::cout << " Filling " << lab << " of " << maximum <<  std::endl;
      typedef itk::ImageRegionIteratorWithIndex<ImageType> RelabelIterator;
      RelabelIterator vfIter( relabel->GetOutput(),
                              relabel->GetOutput()->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        if( itk::Math::FloatAlmostEqual( vfIter.Get(), static_cast<float>( lab ) ) )
          {
          imageout->SetPixel(vfIter.GetIndex(), 1);
          }
        }
      }
    }

  WriteImage<ImageType>(imageout, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int NormalizeImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       option = 0;

  bool useMaskImage = false;
  typename ImageType::Pointer maskImage = nullptr;
  if( argc > argct )
    {

    std::string maskFileName = std::string( argv[argct] );
    ReadImage<ImageType>( maskImage, maskFileName.c_str() );
    if( maskImage.IsNotNull() )
      {
      useMaskImage = true;
      }
    else
      {
      option = atof( argv[argct] );
      }
    }
  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );

  if( useMaskImage )
    {
    itk::ImageRegionIterator<ImageType> It( maskImage,
      maskImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> It2( image,
      image->GetLargestPossibleRegion() );

    float roiMean = 0.0;
    float count = 0.0;
    for( It.GoToBegin(), It2.GoToBegin(); !It.IsAtEnd(); ++It, ++It2 )
      {
      if( ! itk::Math::FloatAlmostEqual( It.Get(), 0.0f ) )
        {
        roiMean += It2.Get();
        count += 1.0f;
        }
      }
    roiMean /= count;

    for( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
      {
      It2.Set( It2.Get() / roiMean );
      }

    if( outname.length() > 3 )
      {
      WriteImage<ImageType>( image, outname.c_str() );
      }
    }
  else
    {
    float         max = 0;
    float         min = 1.e9;
    float         mean = 0.0;
    unsigned long ct = 0;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    Iterator iter( image,  image->GetLargestPossibleRegion() );
    for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      float pix = iter.Get();
      //      if (option == 0) if (pix < 0) pix=0;
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
      if( itk::Math::FloatAlmostEqual( option, 0.0f ) )
        {
        pix = (pix - min) / (max - min);
        iter.Set(pix);
        }
      else
        {
        iter.Set(pix / mean);
        }
      }

    if( outname.length() > 3 )
      {
      WriteImage<ImageType>( image, outname.c_str() );
      }
    }

  return 0;
}

template <unsigned int ImageDimension>
int PrintHeader(int argc, char *argv[])
{
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;

  int         argct = 4;
  std::string fn1 = std::string(argv[argct]);
  if( argc > 20 )
    {
    // std::cout << " k " << std::endl;
    }

  typename ImageType::Pointer image;
  ReadImage<ImageType>(image, fn1.c_str() );
  // std::cout << " Spacing " << image->GetSpacing() << std::endl;
  // std::cout << " Origin " << image->GetOrigin() << std::endl;
  // std::cout << " Direction " << std::endl << image->GetDirection() << std::endl;
  // std::cout << " Size " << std::endl << image->GetLargestPossibleRegion().GetSize() << std::endl;

  //  if (strcmp(operation.c_str(),"n_last_dim") == 0){
  // unsigned int lastdim=image->GetLargestPossibleRegion().GetSize()[ImageDimension-1];
  //   std::ofstream logfile;
  // logfile.open(outname.c_str() );
  // if (logfile.good()  )
  // {
  //  logfile << lastdim << std::endl;
  // }
  // cd // std::cout << lastdim << std::endl;
  // }
  return 1;
}

template <unsigned int ImageDimension>
int GradientImage(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sig = 1;
  if( argc > argct )
    {
    sig = atof(argv[argct]);
    }
  argct++;
  if( sig <= 0 )
    {
    sig = 0.5;
    }
  bool normalize = false;
  if( argc > argct )
    {
    normalize = std::stoi(argv[argct]);
    }

  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetSigma(sig);
  filter->SetInput(image);
  filter->Update();
  image2 = filter->GetOutput();

  if( !normalize )
    {
    WriteImage<ImageType>( image2, outname.c_str() );
    return 0;
    }

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 1 );
  rescaler->SetInput( image2 );
  rescaler->Update();
  WriteImage<ImageType>( rescaler->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int LaplacianImage(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sig = 1;
  if( argc > argct )
    {
    sig = atof(argv[argct]);
    }
  if( sig <= 0 )
    {
    sig = 0.5;
    }
  bool normalize = false;
  if( argc > argct )
    {
    normalize = std::stoi(argv[argct]);
    }

  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );

  typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetSigma(sig);
  filter->SetInput(image);
  filter->Update();
  image2 = filter->GetOutput();
  if( !normalize )
    {
    WriteImage<ImageType>( image2, outname.c_str() );
    return 0;
    }

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 1 );
  rescaler->SetInput( image2 );
  rescaler->Update();
  WriteImage<ImageType>( rescaler->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int CannyImage(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sig = 1;
  if( argc > argct )
    {
    sig = atof(argv[argct]); argct++;
    }
  if( sig <= 0 )
    {
    sig = 0.5;
    }
  PixelType lowerThreshold = -1.0;
  if( argc > argct )
    {
    lowerThreshold = atof(argv[argct]); argct++;
    }
  PixelType upperThreshold = -1.0;
  if( argc > argct )
    {
    upperThreshold = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  typedef itk::CannyEdgeDetectionImageFilter< ImageType, ImageType >
  FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetVariance( sig );
  if ( ! itk::Math::FloatAlmostEqual( upperThreshold, -1.0f ) ) {
  filter->SetUpperThreshold( (PixelType) (upperThreshold) );
  }
  if ( ! itk::Math::FloatAlmostEqual( lowerThreshold, -1.0f ) ) {
  filter->SetLowerThreshold( (PixelType) (lowerThreshold) );
  }
  filter->Update();
  WriteImage<ImageType>( filter->GetOutput(), outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int PoissonDiffusion( int argc, char *argv[])
{
  if( argc < 6 )
    {
    std::cout << "Usage error---not enough arguments.   See help menu."
              << std::endl;
    throw std::exception();
    }

  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension>       LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[4] );
  reader->Update();

  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( reader->GetOutput() );
  duplicator->Update();

  typename ImageType::Pointer output = duplicator->GetOutput();
  output->DisconnectPipeline();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[5] );
  labelReader->Update();

  float label = 1.0;
  if( argc > 7 )
    {
    label = atof( argv[7] );
    }
  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( labelReader->GetOutput() );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetLowerThreshold( label );
  thresholder->SetUpperThreshold( label );
  thresholder->Update();

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType2;
  typename ThresholderType2::Pointer thresholder2 = ThresholderType2::New();
  thresholder2->SetInput( reader->GetOutput() );
  thresholder2->SetOutsideValue( 0 );
  thresholder2->SetInsideValue( 1 );
  thresholder2->SetLowerThreshold( 0.2 );
  thresholder2->SetUpperThreshold( 1.e9 );
  thresholder2->Update();

  float sigma = 1.0;
  if( argc > 6 )
    {
    sigma = atof( argv[6] );
    }

  float convergence = itk::NumericTraits<float>::max();
  float convergenceThreshold = 1e-10;
  if( argc > 9 )
    {
    convergenceThreshold = atof( argv[9] );
    }
  unsigned int maximumNumberOfIterations = 500;
  if( argc > 8 )
    {
    maximumNumberOfIterations = std::stoi( argv[8] );
    }
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  float        lastmean = 0;
  unsigned int iterations = 0;
  while( iterations++ < maximumNumberOfIterations && convergence >= convergenceThreshold )
    {
    // std::cout << "  Iteration " << iterations << ": " << convergence << std::endl;
    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( itk::Math::sqr ( sigma ) );
    smoother->SetMaximumError( 0.01f );
    smoother->SetInput( output );

    /*
    typedef itk::MaximumImageFilter<ImageType, ImageType, ImageType>
      MaximumFilterType;
    typename MaximumFilterType::Pointer maximumFilter = MaximumFilterType::New();
    maximumFilter->SetInput1( smoother->GetOutput() );
    maximumFilter->SetInput2( reader->GetOutput() );

    typedef itk::MultiplyImageFilter<ImageType, LabelImageType, ImageType>
      MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput1( maximumFilter->GetOutput() );
    multiplier->SetInput2( thresholder->GetOutput() );

    typedef itk::AddImageFilter<ImageType, ImageType, ImageType>
     AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( reader->GetOutput() );
    adder->SetInput2( multiplier->GetOutput() );

    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType>
      SubtracterType;
    typename SubtracterType::Pointer subtracter = SubtracterType::New();
    subtracter->SetInput1( adder->GetOutput() );
    subtracter->SetInput2( output );
    */
    typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType>
      StatsFilterType;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetInput( smoother->GetOutput() );
    stats->SetLabelInput( thresholder->GetOutput() );
    stats->Update();

    convergence = static_cast<float>( stats->GetMean( 1 ) ) - lastmean;
    lastmean = static_cast<float>( stats->GetMean( 1 ) );
    output =  smoother->GetOutput();
    output->DisconnectPipeline();

    Iterator vfIter( output,  output->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      if( itk::Math::FloatAlmostEqual( thresholder2->GetOutput()->GetPixel(vfIter.GetIndex() ), 1.0f ) )
        {
        vfIter.Set(reader->GetOutput()->GetPixel(vfIter.GetIndex() ) );
        }
      }
    }

  Iterator vfIter( output,  output->GetLargestPossibleRegion() );
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if( itk::Math::FloatAlmostEqual( thresholder2->GetOutput()->GetPixel(vfIter.GetIndex() ), 1.0f ) )
      {
      vfIter.Set(reader->GetOutput()->GetPixel(vfIter.GetIndex() ) );
      }
    if(  thresholder->GetOutput()->GetPixel(vfIter.GetIndex() ) == 0 )
      {
      vfIter.Set(0);
      }
    }

  WriteImage<ImageType>( output, argv[2] );

  return 0;
}

template <unsigned int ImageDimension>
void
RemoveLabelInterfaces(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl; return;
    }
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer input = nullptr;
  ReadImage<ImageType>(input, fn1.c_str() );

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

// std::cout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename ImageType::PixelType p = GHood.GetCenterPixel();
    typename ImageType::IndexType ind = GHood.GetIndex();
    if( p > 0  )
      {
      bool atedge = false;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        if( GHood.GetPixel(i) > 0 && ! itk::Math::FloatAlmostEqual( GHood.GetPixel(i), p ) )
          {
          atedge = true;
          }
        }
      if( atedge )
        {
        input->SetPixel(ind, 0);
        }
      }
    ++GHood;
    }

  WriteImage<ImageType>(input, outname.c_str() );

  return;
}

template <unsigned int ImageDimension>
void
ReplaceVoxelValue( int argc, char *argv[] )
{
  // Usage:  ImageMath 3 output.nii.gz ReplaceVoxelValue input.nii.gz a b c

  if( argc < 8 )
    {
    // std::cout << " too few options " << std::endl; return;
    }

  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  const std::string outputFile = std::string( argv[2] );
  const std::string inputFile = std::string( argv[4] );
  const PixelType   lowThreshold = atof( argv[5] );
  const PixelType   highThreshold = atof( argv[6] );
  const PixelType   replacementValue = atof( argv[7] );

  typename ImageType::Pointer inputImage = nullptr;
  ReadImage<ImageType>( inputImage, inputFile.c_str() );

  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  typename itk::ImageRegionIterator<ImageType> ItI( inputImage, inputImage->GetRequestedRegion() );
  typename itk::ImageRegionIterator<ImageType> ItO( outputImage, outputImage->GetRequestedRegion() );
  for( ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO )
    {
    PixelType inputVoxel = ItI.Get();
    if( inputVoxel >= lowThreshold && inputVoxel <= highThreshold )
      {
      ItO.Set( replacementValue );
      }
    else
      {
      ItO.Set( inputVoxel );
      }
    }

  WriteImage<ImageType>( outputImage, outputFile.c_str() );
}

template <unsigned int ImageDimension>
void
EnumerateLabelInterfaces(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string outname2 = "";
  if( argc  > argct )
    {
    outname2 = std::string(argv[argct]);
    }
  argct++;
  float nfrac = 0.1;
  if( argc  > argct )
    {
    nfrac = atof(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer input = nullptr;
  ReadImage<ImageType>(input, fn1.c_str() );
  typename ImageType::Pointer output = nullptr;
  ReadImage<ImageType>(output, fn1.c_str() );
  output->FillBuffer(0);
  typename ImageType::Pointer colored = nullptr;
  ReadImage<ImageType>(colored, fn1.c_str() );
  colored->FillBuffer(0);

  unsigned int max = 0;
  Iterator     o_iter( input, input->GetLargestPossibleRegion() );
  o_iter.GoToBegin();
  while( !o_iter.IsAtEnd() )
    {
    unsigned int label = (unsigned int) o_iter.Get();
    if( label > max )
      {
      max = label;
      }
    ++o_iter;
    }

  // std::cout << " Max Label " << max << std::endl;
  typedef itk::Image<float, 2> myInterfaceImageType;
  typename myInterfaceImageType::SizeType size;
  size[0] = max + 1;
  size[1] = max + 1;
  typename myInterfaceImageType::RegionType region;
  typename myInterfaceImageType::SpacingType spacing;
  spacing.Fill(1);
  typename myInterfaceImageType::PointType origin;
  origin.Fill(max / 2);
  region.SetSize(size);
  typename myInterfaceImageType::DirectionType direction;
  direction.SetIdentity();
  typename myInterfaceImageType::Pointer faceimage =
    AllocImage<myInterfaceImageType>(region, spacing, origin, direction, 0);

  typename myInterfaceImageType::Pointer colorimage =
    AllocImage<myInterfaceImageType>(region, spacing, origin, direction, 0);


// we can use this to compute a 4-coloring of the brain

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

// std::cout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    const typename ImageType::PixelType & p = GHood.GetCenterPixel();
    if( p > 0  )
      {
      bool atedge = false;

      unsigned long linearinda = 0, linearindb = 0, linearind = 0;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        if( GHood.GetPixel(i) > 0 && ! itk::Math::FloatAlmostEqual( GHood.GetPixel(i), p ) )
          {
          atedge = true;
          typename myInterfaceImageType::IndexType inda, indb;
          inda[0] = (long int)p;  inda[1] = (long int)GHood.GetPixel(i);
          indb[1] = (long int)p;  indb[0] = (long int)GHood.GetPixel(i);
          faceimage->SetPixel(inda, faceimage->GetPixel(inda) + 1);
          faceimage->SetPixel(indb, faceimage->GetPixel(indb) + 1);
          linearinda = inda[1] * size[0] + inda[0];
          linearindb = indb[1] * size[0] + indb[0];
          if( linearinda < linearindb )
            {
            linearind = linearinda;
            }
          else
            {
            linearind = linearindb;
            }
          }
        }
      if( atedge )
        {
        const typename ImageType::IndexType & ind = GHood.GetIndex();
        output->SetPixel(ind, linearind);
        }
      }
    ++GHood;
    }

// first normalize the interfaces
  typename myInterfaceImageType::IndexType find;
  typename myInterfaceImageType::IndexType find2;
  find2.Fill(0);
  for( unsigned int j = 0; j <= max; j++ )
    {
    find[0] = j;
    float total = 0;
    for( unsigned int i = 0; i <= max; i++ )
      {
      find[1] = i;
      total += faceimage->GetPixel(find);
//      if ( faceimage->GetPixel(find) > 50 )
      // std::cout << i <<"  :: " << faceimage->GetPixel(find)  << std::endl;
      }
    // std::cout << " total interfaces for label :  " << j << " are " << total << std::endl;
    for( unsigned int i = 0; i <= max; i++ )
      {
      find[1] = i;
      if( total > 0 )
        {
        faceimage->SetPixel(find, faceimage->GetPixel(find) / total);
        }
      if( faceimage->GetPixel(find) >  0.01f )
        {
        // std::cout << i << "  :: " << faceimage->GetPixel(find)  << std::endl;
        }
      }
    }

// colors 1 , 2 , 3 , 4
// set all to color 1.
  typedef std::vector<unsigned int> ColorSetType;
  colorimage->FillBuffer(0);
//  for (unsigned int j=1; j<=max; j++)
  for( unsigned int j = max; j >= 1; j-- )
    {
    ColorSetType myColorSet1;
    ColorSetType myColorSet2;
    ColorSetType myColorSet3;
    ColorSetType myColorSet4;
    ColorSetType myColorSet5;
    find[0] = j;
// list all indices that have an interface with j
//    for (unsigned int i=1; i<=max; i++)
    for( unsigned int i = max; i >= 1; i-- )
      {
      if( i !=  j )
        {
        find[1] = i;
        find2[0] = i;
        unsigned int color = (unsigned int) colorimage->GetPixel(find2);
        if( faceimage->GetPixel(find) >  nfrac ) // then color
          {
          if( color == 1 )
            {
            myColorSet1.push_back( (unsigned int) i );
            }
          if( color == 2 )
            {
            myColorSet2.push_back( (unsigned int) i );
            }
          if( color == 3 )
            {
            myColorSet3.push_back( (unsigned int) i );
            }
          if( color == 4 )
            {
            myColorSet4.push_back( (unsigned int) i );
            }
          if( color == 5 )
            {
            myColorSet5.push_back( (unsigned int) i );
            }
// now you know all of j's neighbors and their colors
          }
        }
      }

// we have to recolor j so that his color is different
// than any in the label set
    unsigned int okcolor = 6;
    find[0] = j;
    find[1] = 0;
    if( myColorSet1.empty()  )
      {
      okcolor = 1;
      }
    else if( myColorSet2.empty()  )
      {
      okcolor = 2;
      }
    else if( myColorSet3.empty()  )
      {
      okcolor = 3;
      }
    else if( myColorSet4.empty()  )
      {
      okcolor = 4;
      }
    else if( myColorSet5.empty()  )
      {
      okcolor = 5;
      }

    colorimage->SetPixel(find, okcolor);
    // std::cout << " Label " << j << " color " << okcolor << std::endl;
    }

  o_iter.GoToBegin();
  colored->FillBuffer(0);
  while( !o_iter.IsAtEnd() )
    {
    unsigned int label = (unsigned int) o_iter.Get();
    if( label > 0 )
      {
      find[0] = label;
      find[1] = 0;
      unsigned int color = (unsigned int)colorimage->GetPixel(find);
      colored->SetPixel(o_iter.GetIndex(), color);
      }
    ++o_iter;
    }

  WriteImage<ImageType>(output, outname.c_str() );
  WriteImage<myInterfaceImageType>(faceimage, (std::string("face") + outname).c_str() );
  if( outname2.length() > 3 )
    {
    WriteImage<ImageType>(colored, outname2.c_str() );
    }

  return;
}

template <unsigned int ImageDimension>
int CountVoxelDifference(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  std::string maskfn = "";
  if(  argc > argct )
    {
    maskfn = std::string(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  ReadImage<ImageType>(image2, fn2.c_str() );
  ReadImage<ImageType>(mask, maskfn.c_str() );

  unsigned long maskct = 0;
  float         err = 0, negerr = 0;

  Iterator It( mask, mask->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() > 0  )
      {
      float p1 = image1->GetPixel(It.GetIndex() );
      float p2 = image2->GetPixel(It.GetIndex() );
      float locerr = p1 - p2;
      err += static_cast<float>( fabs(locerr) );
      if( locerr < 0 )
        {
        negerr += locerr;
        }
      maskct++;
      }
    }

  if( maskct == 0 )
    {
    maskct = 1;
    }
  std::ofstream logfile;
  logfile.open(outname.c_str() );
  if( logfile.good()  )
    {
    logfile << " Err " <<  " : " << err << " %ER " <<  " : " << err / (maskct) * 100.0f << " NER " <<  " : " << negerr
      / maskct * 100.0f << std::endl;
    }
  return 0;
}

template <unsigned int ImageDimension>
int DiceAndMinDistSum(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  std::string outdistfn = "";
  if(  argc > argct )
    {
    outdistfn = std::string(argv[argct]);   argct++;
    }
  std::string diceimagename = ANTSGetFilePrefix(outname.c_str() ) + std::string("dice.nii.gz");
  std::string mdsimagename = ANTSGetFilePrefix(outname.c_str() ) + std::string("mds.nii.gz");

  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::Pointer outdist = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }
  if( outdistfn.length() > 3 )
    {
    ReadImage<ImageType>(outdist, fn1.c_str() );
    outdist->FillBuffer(0);
    }
  // 1: for each label, and each image, compute the distance map
  // 2: sum min-dist(i) over every point in each image
  // 3: take average over all points
  typedef float                  PixelType;
  typedef std::vector<PixelType> LabelSetType;
  LabelSetType myLabelSet1;
    {
/** count the labels in the image */
    Iterator It( image1, image1->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( fabs(label) > 0 )
        {
        if( find( myLabelSet1.begin(), myLabelSet1.end(), label )
            == myLabelSet1.end() )
          {
          myLabelSet1.push_back( label );
          }
        }
      }
    }

  LabelSetType myLabelSet2;
  unsigned int labct = 0;
    {
    Iterator It( image2, image2->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( fabs(label) > 0 )
        {
        if( find( myLabelSet2.begin(), myLabelSet2.end(), label )
            == myLabelSet2.end()   &&
            find( myLabelSet1.begin(), myLabelSet1.end(), label )
            != myLabelSet1.end() )
          {
          myLabelSet2.push_back( label );
          labct++;
          }
        }
      }
    }

  vnl_vector<double> distances(labct, 0.0);
  vnl_vector<double> dicevals(labct, 0.0);
  vnl_vector<double> rovals(labct, 0.0);
  vnl_vector<double> tpvals(labct, 0.0);
  vnl_vector<double> tpvals2(labct, 0.0);

  /** now we have the common labels */
  std::ofstream logfile;
  logfile.open(outname.c_str() );
  labct = 0;
  typename LabelSetType::const_iterator it;
  for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )
    {
    typename ImageType::Pointer mask1 = BinaryThreshold<ImageType>(*it, *it, 1, image1);
    typename ImageType::Pointer surf = nullptr;
    typename ImageType::Pointer d1 = nullptr;
    typename ImageType::Pointer d2 = nullptr;

    float count1 = 0; // count vox in label 1
    float count2 = 0; // count vox in label 2
    float counti = 0; // count vox intersection in label 1 and label 2
    float countu = 0; // count vox union label 1 and label 2
    //    float surfdist = 0, surfct = 0;
    //    float truepos=0;
    if( outdist )
      {
      surf = LabelSurface<ImageType>(mask1, mask1);
      //    WriteImage<ImageType>(surf,outdistfn.c_str());
      // throw std::exception();
      typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;
      typename  FilterType::Pointer dfilter1 = FilterType::New();
      dfilter1->InputIsBinaryOn();
      dfilter1->SetUseImageSpacing(true);
      dfilter1->SetInput(mask1);

      typename ImageType::Pointer mask2 = BinaryThreshold<ImageType>(*it, *it, 1, image2);
      typename  FilterType::Pointer dfilter2 = FilterType::New();
      dfilter2->InputIsBinaryOn();
      dfilter2->SetUseImageSpacing(true);
      dfilter2->SetInput(mask2);

      dfilter1->Update();
      d1 = dfilter1->GetOutput();

      dfilter2->Update();
      d2 = dfilter2->GetOutput();
      }

    float    dist1 = 0;
    float    dist2 = 0;
    Iterator It2( image2, image2->GetLargestPossibleRegion() );
    for( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
      {
      if( itk::Math::FloatAlmostEqual( It2.Get(), *it ) || itk::Math::FloatAlmostEqual( image1->GetPixel( It2.GetIndex() ), *it ) )
        {
        countu++;
        }
      if( itk::Math::FloatAlmostEqual( It2.Get(), *it ) && itk::Math::FloatAlmostEqual( image1->GetPixel( It2.GetIndex() ), *it ) )
        {
        counti++;
        }

      if( itk::Math::FloatAlmostEqual( It2.Get(), *it ) )
        {
        count2++;
        if( d1 )
          {
          dist2 += d1->GetPixel(It2.GetIndex() );
          }
        }
      if( itk::Math::FloatAlmostEqual( image1->GetPixel(It2.GetIndex() ), *it ) )
        {
        count1++;
        if( d2 )
          {
          dist1 += d2->GetPixel(It2.GetIndex() );
          }
        }
      if( outdist )
        {
        if( surf->GetPixel(It2.GetIndex() ) > 0 )
          {
          float sdist = d2->GetPixel(It2.GetIndex() );
          outdist->SetPixel(It2.GetIndex(), sdist);
          // surfdist += sdist;
          // surfct += 1;
          }
        }
      }
    //    // std::cout << " sdist " << surfdist << " sct " << surfct << std::endl

    if( outdist )
      {
      WriteImage<ImageType>(outdist, outdistfn.c_str() );
      }

    if( count2 + count1 > 0 )
      {
      distances[labct] += static_cast<double>( (dist2 + dist1) / (count2 + count1) );
      dicevals[labct] = 2.0f * counti / (count2 + count1);
      rovals[labct] = counti / (countu);
      tpvals[labct] = counti / count1;
      tpvals2[labct] = counti / count2;
      }
    labct++;
    }

  labct = 0;
  float sum = 0, sumdice = 0, sumro = 0;
  //  std::sort(myLabelSet2.begin(),myLabelSet2.end());
  unsigned long labelcount = 0;
  for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )
    {
    labelcount++;
    }
  float        temp = sqrt( (float)labelcount);
  unsigned int squareimagesize = (unsigned int)(temp + 1);
  typedef itk::Image<float, 2> TwoDImageType;
  typename TwoDImageType::RegionType newregion;
  typename TwoDImageType::SizeType size;
  size[0] = size[1] = squareimagesize;
  newregion.SetSize(size);
  typename TwoDImageType::Pointer squareimage =
    AllocImage<TwoDImageType>(newregion, 0.0);

  typename TwoDImageType::Pointer squareimage2 =
    AllocImage<TwoDImageType>(newregion, 0.0);

  std::vector<std::string> ColumnHeaders;
  std::vector<std::string> RowHeaders;
  unsigned int             NumberOfLabels = labelcount;

  if( outdist )
    {
    vnl_matrix<double> OutputValues(NumberOfLabels, 2);

    ColumnHeaders.emplace_back("Label Name");
    ColumnHeaders.emplace_back("Min_Distance");
    ColumnHeaders.emplace_back("Dice");
    for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )
      {
      OutputValues(labct, 0) = distances[labct];
      OutputValues(labct, 1) = dicevals[labct];
      std::string LabelName = "Label_";
      int         LabelNumber;
      LabelNumber = *it;
      char LabelNumberAsString[50];
      sprintf(LabelNumberAsString, "%.2d", LabelNumber);
      LabelName = LabelName + LabelNumberAsString;
      RowHeaders.push_back(LabelName);
      labct++;
      }
    typedef  itk::CSVNumericObjectFileWriter<double, 1, 1> CSVType;
    CSVType::Pointer OutputCSV = CSVType::New();
    OutputCSV->SetInput(&OutputValues);
    OutputCSV->SetFileName(std::string(outname.c_str() ) + ".csv");
    OutputCSV->SetColumnHeaders(ColumnHeaders);
    OutputCSV->SetRowHeaders(RowHeaders);
    try
      {
      OutputCSV->Write();
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      }
    }
  else
    {
    vnl_matrix<double> OutputValues(NumberOfLabels, 4);
    ColumnHeaders.emplace_back("Label Name");
    ColumnHeaders.emplace_back("Dice");
    ColumnHeaders.emplace_back("RO");
    ColumnHeaders.emplace_back("Percent_of_Region_1_In_Overlap");
    ColumnHeaders.emplace_back("Percent_of_Region_2_In_Overlap");
    for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )
      {
      OutputValues(labct, 0) = dicevals[labct];
      OutputValues(labct, 1) = rovals[labct];
      OutputValues(labct, 2) = tpvals[labct];
      OutputValues(labct, 3) = tpvals2[labct];

      std::string LabelName = "Label_";
      int         LabelNumber;
      LabelNumber = *it;
      char LabelNumberAsString[50];
      sprintf(LabelNumberAsString, "%.2d", LabelNumber);
      LabelName = LabelName + LabelNumberAsString;
      RowHeaders.push_back(LabelName);
      labct++;
      }

    typedef  itk::CSVNumericObjectFileWriter<double, 1, 1> CSVType;
    CSVType::Pointer OutputCSV = CSVType::New();
    OutputCSV->SetInput(&OutputValues);
    OutputCSV->SetFileName(std::string(outname.c_str() ) + ".csv");
    OutputCSV->SetColumnHeaders(ColumnHeaders);
    OutputCSV->SetRowHeaders(RowHeaders);
    try
      {
      OutputCSV->Write();
      // std::cout << "Output written to " << outname.c_str() << ".csv." << std::endl;
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      }
    }

  labelcount = 0;
  labct = 0;
  for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )

    {
    sum += static_cast<float>( distances[labct] );
    sumdice += static_cast<float>( dicevals[labct] );
    sumro += static_cast<float>( rovals[labct] );
// square image
    squareimage->GetBufferPointer()[labct] = distances[labct];
    squareimage2->GetBufferPointer()[labct] = dicevals[labct];
    labct++;
    }

  if( labct == 0 )
    {
    labct = 1;
    }
  if( logfile.good()  )
    {
    if( outdist )
      {
      logfile << " AvgMinDist " <<  " : " << sum / (float)labct << " AvgDice " <<  " : " << sumdice / (float)labct
              << " AvgRO " << " : " << sumro / (float)labct << std::endl;
      }
    else
      {
      logfile << " AvgDice " <<  " : " << sumdice / (float)labct << std::endl;
      }
    }

  WriteImage<TwoDImageType>(squareimage, mdsimagename.c_str() );
  WriteImage<TwoDImageType>(squareimage2, diceimagename.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int Lipschitz( int argc, char *argv[] )
{
  if( argc > 20 )
    {
    // std::cout << " k " << std::endl;
    }
  // std::cout << " Compute Lipschitz continuity of the mapping " << std::endl;

  typedef float                                              RealType;
  typedef itk::Image<RealType, ImageDimension>               RealImageType;
  typedef itk::Vector<RealType, ImageDimension>              VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string vecname = std::string(argv[argct]);   argct++;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( vecname.c_str() );
  // reader->SetUseAvantsNamingConvention( true );
  reader->Update();
  typename VectorImageType::Pointer vecimage = reader->GetOutput();

  typename RealImageType::Pointer lipcon =
    AllocImage<RealImageType>(vecimage, 0);

  itk::TimeProbe timer;
  timer.Start();
/** global Lipschitz points */
  typename VectorImageType::PointType gx; gx.Fill(0);
  typename VectorImageType::PointType gxt; gxt.Fill(0);
  typename VectorImageType::PointType gy; gy.Fill(0);
  typename VectorImageType::PointType gyt; gyt.Fill(0);
  typename VectorImageType::PointType lx;  lx.Fill(0);
  typename VectorImageType::PointType lxt;  lxt.Fill(0);
  typename VectorImageType::PointType ly;  ly.Fill(0);
  typename VectorImageType::PointType lyt;  lyt.Fill(0);
  float         globalmaxval = 0;
  unsigned long ct1 = 0;
//  unsigned long numpx=vecimage->GetBufferedRegion().GetNumberOfPixels();
  Iterator It1( vecimage, vecimage->GetLargestPossibleRegion() );
  for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
    {
    ct1++;
    typename VectorImageType::PointType x;
    typename VectorImageType::PointType xt;
    vecimage->TransformIndexToPhysicalPoint( It1.GetIndex(), x);
    typename VectorImageType::PixelType vecx = It1.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      xt[i] = x[i] + static_cast<double>( vecx[i] );
      }

    float         localmaxval = 0;
    unsigned long ct2 = 0;
    Iterator      It2( vecimage, vecimage->GetLargestPossibleRegion() );
    for( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
      {
      ct2++;
      if( ct2 != ct1 )
        {
        typename VectorImageType::PointType y;
        typename VectorImageType::PointType yt;
        vecimage->TransformIndexToPhysicalPoint( It2.GetIndex(), y);
        typename VectorImageType::PixelType vecy = It2.Get();

        double numer = 0.0;
        double denom = 0.0;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          yt[i] = y[i] + static_cast<double>( vecy[i] );
          numer += (yt[i] - xt[i]) * (yt[i] - xt[i]);
          denom += (y[i] - x[i]) * (y[i] - x[i]);
          }
        numer = sqrt(numer);
        denom = sqrt(denom);
        float localval = numer / denom;
        if( localval > localmaxval )
          {
          lx = x;
          ly = y;
          lxt = xt;
          lyt = yt;
          localmaxval = localval;
          }
        }
      }
    if( localmaxval > globalmaxval )
      {
      gx = lx;
      gy = ly;
      gxt = lxt;
      gyt = lyt;
      globalmaxval = localmaxval;
      }
    lipcon->SetPixel(It1.GetIndex(), localmaxval);
//      if (ct1 % 1000 == 0) // std::cout << " Progress : " << (float ) ct1 / (float) numpx *100.0 << " val " <<
// localmaxval << std::endl;
    }

  // std::cout << " Lipschitz continuity related to: " << globalmaxval << std::endl;
  // std::cout << " Tx :  " << gxt << "  Ty: " << gyt << std::endl;
  // std::cout << " x :  " << gx << "  y: " << gy << std::endl;
  timer.Stop();
//    // std::cout << "Elapsed time: " << timer.GetMeanTime()  << std::endl;

  if( outname.length() > 3 )
    {
    WriteImage<RealImageType>( lipcon, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int ExtractVectorComponent( int argc, char *argv[] )
{
  if( argc <= 2 )
    {
    // std::cout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                       PixelType;
  typedef itk::VectorImage<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PixelType, ImageDimension>       RealImageType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  inname = std::string(argv[argct]);   argct++;
  unsigned int whichvec = std::stoi(argv[argct]);   argct++;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( inname.c_str() );
  reader1->Update();
  typename ImageType::Pointer vecimage = reader1->GetOutput();
  if( whichvec >= vecimage->GetVectorLength() )
    {
    std::cout << " input image " << inname << " only has " << vecimage->GetVectorLength() << " components "
              << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    typename RealImageType::Pointer component
      = AllocImage<RealImageType>(vecimage, 0);
    typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
    Iterator It1( vecimage, vecimage->GetLargestPossibleRegion() );
    for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
      {
      component->SetPixel(It1.GetIndex(), It1.Get()[whichvec]);
      }
    WriteImage<RealImageType>( component, outname.c_str() );
    }
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int InvId( int argc, char *argv[] )
{
  if( argc > 2 )
    {
    // std::cout << " Compute  phi(  phi^{-1}(x)) " << std::endl;
    }
  else
    {
    return 1;
    }
  typedef float                                              RealType;
  typedef itk::Image<RealType, ImageDimension>               RealImageType;
  typedef itk::Vector<RealType, ImageDimension>              VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string vecname1 = std::string(argv[argct]);   argct++;
  std::string vecname2 = std::string(argv[argct]);   argct++;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( vecname1.c_str() );
  //  reader1->SetUseAvantsNamingConvention( true );
  reader1->Update();
  typename VectorImageType::Pointer vecimage1 = reader1->GetOutput();
  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( vecname2.c_str() );
  //  reader2->SetUseAvantsNamingConvention( true );
  reader2->Update();
  typename VectorImageType::Pointer vecimage2 = reader2->GetOutput();

  typename RealImageType::Pointer invid =
    AllocImage<RealImageType>(vecimage1, 0);

  itk::TimeProbe timer;
  timer.Start();

  typedef itk::VectorLinearInterpolateImageFunction<VectorImageType, float> DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
  vinterp->SetInputImage(vecimage2);

/** global Lipschitz points */
  typename VectorImageType::PointType gx;  gx.Fill(0);
  double    globalmaxval = 0;
  Iterator It1( vecimage1, vecimage1->GetLargestPossibleRegion() );
  for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
    {
    typename VectorImageType::PointType x;
    typename DefaultInterpolatorType::PointType xt, yt;
    vecimage1->TransformIndexToPhysicalPoint( It1.GetIndex(), x);
    typename VectorImageType::PixelType vecx = It1.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      xt[i] = x[i] + static_cast<double>( vecx[i] );
      }
// above, we warped the point -- now we index the other field with this
// point

    typename VectorImageType::PixelType disp;
    if( vinterp->IsInsideBuffer(xt) )
      {
      disp = vinterp->Evaluate(xt);
      }
    else
      {
      disp.Fill(0);
      }
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      yt[jj] = xt[jj] + disp[jj];
      }

    double error = 0.0;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      error += static_cast<double>( itk::Math::sqr( static_cast<double>( yt[jj] ) - static_cast<double>( x[jj] ) ) );
      }
    error = sqrt(error);
    if( error >  globalmaxval )
      {
      globalmaxval = error; gx = x;
      }
    invid->SetPixel(It1.GetIndex(), static_cast<float>( error ) );
    }
  // std::cout << " Max error " << globalmaxval << " at " << gx << std::endl;
  timer.Stop();
//    // std::cout << "Elapsed time: " << timer.GetMeanTime()  << std::endl;

  WriteImage<RealImageType>( invid, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int ReplicateDisplacement( int argc, char *argv[] )
{
  if( argc > 6 )
    {
    // std::cout << " Replicate a ND displacement to ND+1 dimensions " << std::endl;
    }
  else
    {
    // std::cout << "ImageMath 3 out4DWarp.nii.gz ReplicateDisplacement in3DWarp.nii.gz nreplications time-spacing time-origin" << std::endl;
    // std::cout << "ImageMath 3 out4DWarp.nii.gz ReplicateDisplacement in3DWarp.nii.gz 10 2.5 0.0" << std::endl;
    return 1;
    }
  typedef float                                              RealType;
  typedef itk::Vector<RealType, ImageDimension>              VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> Iterator;

  typedef itk::Vector<RealType, ImageDimension+1>            VectorRType;
  typedef itk::Image<VectorRType, ImageDimension+1>          VectorRImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string vecname1 = std::string(argv[argct]);   argct++;
  unsigned int timedims = std::stoi(argv[argct]);  argct++;
  RealType tr = atof(argv[argct]);  argct++;
  RealType torigin = atof(argv[argct]);  argct++;
  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( vecname1.c_str() );
  reader1->Update();
  typename VectorImageType::Pointer vecimage1 = reader1->GetOutput();
  typename VectorRImageType::Pointer outputImage = VectorRImageType::New();
  typename VectorRImageType::RegionType outRegion;
  typename VectorRImageType::SizeType outSize;
  typename VectorRImageType::SpacingType outSpacing;
  typename VectorRImageType::PointType outOrigin;
  typename VectorRImageType::DirectionType outDirection;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    outSize[d] = vecimage1->GetLargestPossibleRegion().GetSize()[d];
    outSpacing[d] = vecimage1->GetSpacing()[d];
    outOrigin[d] = vecimage1->GetOrigin()[d];
    for( unsigned int e = 0; e < ImageDimension; e++ )
      {
      outDirection(e, d) = vecimage1->GetDirection() (e, d);
      }
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    outDirection(d, ImageDimension) = 0;
    outDirection(ImageDimension, d) = 0;
    }
  outDirection(ImageDimension, ImageDimension) = 1.0;
  outSize[ImageDimension] = timedims;
  outSpacing[ImageDimension] = tr;
  outOrigin[ImageDimension] = torigin;
  outRegion.SetSize( outSize );
  outputImage->SetRegions( outRegion );
  outputImage->SetSpacing( outSpacing );
  outputImage->SetOrigin( outOrigin );
  outputImage->SetDirection( outDirection );
  outputImage->Allocate();
  VectorRType vec;
  vec.Fill( 0 );
  outputImage->FillBuffer( vec );
  // perform the replication
  typename VectorImageType::IndexType ind;
  typename VectorRImageType::IndexType indp1;
  Iterator It1( vecimage1, vecimage1->GetLargestPossibleRegion() );
  for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
    {
    ind = It1.GetIndex();
    typename VectorImageType::PixelType vecx = It1.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      vec[i] = vecx[i];
      indp1[i] = ind[i];
      }
    for( unsigned int i = 0; i < timedims; i++ )
      {
      indp1[ImageDimension] = i;
      outputImage->SetPixel( indp1, vec );
      }
    }
  WriteImage<VectorRImageType>( outputImage, outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int ShiftImageSlicesInTime( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    return 1;
    }
  typedef float                                              RealType;
  typedef itk::Image<RealType, ImageDimension>               ImageType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  unsigned int shiftamount = 1;
  if(  argc > argct )
    {
    shiftamount = std::stoi( argv[argct] );  argct++;
    }
  unsigned int shiftdim = ImageDimension - 1;
  if(  argc > argct )
    {
    shiftdim = std::stoi( argv[argct] );  argct++;
    }
  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  typedef itk::CyclicShiftImageFilter<ImageType, ImageType> FilterType;
  typename  FilterType::Pointer cyfilter = FilterType::New();
  cyfilter->SetInput(image);
  typename FilterType::OffsetType myshift;
  myshift.Fill( 0 );
  myshift[shiftdim]=shiftamount;
  cyfilter->SetShift( myshift );
  WriteImage<ImageType>(cyfilter->GetOutput(), outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int ReplicateImage( int argc, char *argv[] )
{
  if( argc > 6 )
    {
    // std::cout << " Replicate a ND image to ND+1 dimensions " << std::endl;
    }
  else
    {
    // std::cout << "ImageMath 3 out4D.nii.gz ReplicateImage in3D.nii.gz nreplications time-spacing time-origin" << std::endl;
    // std::cout << "ImageMath 3 out4D.nii.gz ReplicateImage in3D.nii.gz 10 2.5 0.0" << std::endl;
    return 1;
    }
  typedef float                                              RealType;
  typedef itk::Image<RealType, ImageDimension>               ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>       Iterator;
  typedef itk::Image<RealType, ImageDimension+1>             RImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string vecname1 = std::string(argv[argct]);   argct++;
  unsigned int timedims = std::stoi(argv[argct]);  argct++;
  float tr = atof(argv[argct]);  argct++;
  float torigin = atof(argv[argct]);  argct++;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( vecname1.c_str() );
  reader1->Update();
  typename ImageType::Pointer vecimage1 = reader1->GetOutput();
  typename RImageType::Pointer outputImage = RImageType::New();
  typename RImageType::RegionType outRegion;
  typename RImageType::SizeType outSize;
  typename RImageType::SpacingType outSpacing;
  typename RImageType::PointType outOrigin;
  typename RImageType::DirectionType outDirection;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    outSize[d] = vecimage1->GetLargestPossibleRegion().GetSize()[d];
    outSpacing[d] = vecimage1->GetSpacing()[d];
    outOrigin[d] = vecimage1->GetOrigin()[d];
    for( unsigned int e = 0; e < ImageDimension; e++ )
      {
      outDirection(e, d) = vecimage1->GetDirection() (e, d);
      }
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    outDirection(d, ImageDimension) = 0;
    outDirection(ImageDimension, d) = 0;
    }
  outDirection(ImageDimension, ImageDimension) = 1.0;
  outSize[ImageDimension] = timedims;
  outSpacing[ImageDimension] = tr;
  outOrigin[ImageDimension] = torigin;
  outRegion.SetSize( outSize );
  outputImage->SetRegions( outRegion );
  outputImage->SetSpacing( outSpacing );
  outputImage->SetOrigin( outOrigin );
  outputImage->SetDirection( outDirection );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );
  // perform the replication
  typename ImageType::IndexType ind;
  typename RImageType::IndexType indp1;
  ind.Fill(0);
  indp1.Fill( 0 );

  Iterator It1( vecimage1, vecimage1->GetLargestPossibleRegion() );
  for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
    {
    ind = It1.GetIndex();
    typename ImageType::PixelType vecx = It1.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      indp1[i] = ind[i];
      }
    for( unsigned int i = 0; i < timedims; i++ )
      {
      indp1[ImageDimension] = i;
      outputImage->SetPixel( indp1, vecx );
      }
    }
  WriteImage<RImageType>( outputImage, outname.c_str() );
  return 0;
}


template <unsigned int ImageDimension>
int LabelStats(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  std::string       imagename = ANTSGetFilePrefix(outname.c_str() ) + std::string("_square.nii.gz");
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer valimage = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(valimage, fn2.c_str() );
    }

  typedef float                  PixelType;
  typedef std::vector<PixelType> LabelSetType;
  LabelSetType myLabelSet;
/** count the labels in the image */
  Iterator It( image, image->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PixelType label = It.Get();
    if( fabs(label) > 0 )
      {
      if( find( myLabelSet.begin(), myLabelSet.end(), label )
          == myLabelSet.end() )
        {
        myLabelSet.push_back( label );
        }
      }
    }

  // compute the voxel volume
  typename ImageType::SpacingType spacing = image->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= static_cast<float>( spacing[i] );
    }

  std::ofstream logfile;
  logfile.open(outname.c_str() );
  logfile << "x,y,z,t,label,mass,volume,count" << std::endl;

  std::sort(myLabelSet.begin(), myLabelSet.end() );
  typename LabelSetType::const_iterator it;
  unsigned long labelcount = 0;
  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
    {
    labelcount++;
    }

  float        temp = sqrt( (float)labelcount);
  unsigned int squareimagesize = (unsigned int)(temp + 1);

  typedef itk::Image<float, 2> TwoDImageType;
  typename TwoDImageType::RegionType newregion;
  typename TwoDImageType::SizeType size;
  size[0] = size[1] = squareimagesize;
  newregion.SetSize(size);
  typename TwoDImageType::Pointer squareimage =
    AllocImage<TwoDImageType>(newregion, 0);

  labelcount = 0;
  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
    {
    float currentlabel = *it;
    float totalvolume = 0;
    float totalmass = 0;
    float totalct = 0;
    typename ImageType::PointType myCenterOfMass;
    myCenterOfMass.Fill(0);
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( itk::Math::FloatAlmostEqual( label, currentlabel ) )
        {
        totalct += 1;
        if( valimage )
          {
          totalmass += valimage->GetPixel( It.GetIndex() );
          }

        // compute center of mass
        typename ImageType::PointType point;
        image->TransformIndexToPhysicalPoint( It.GetIndex(), point );
        for( unsigned int i = 0; i < spacing.Size(); i++ )
          {
          myCenterOfMass[i] += point[i];
          }
        }
      }

    totalvolume = volumeelement * totalct;

    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= static_cast<double>( totalct );
      }

    if( !valimage )
      {
//      // std::cout << " Volume Of Label " << *it << " is " << totalvolume <<   "  Avg-Location " << myCenterOfMass
          //      << std::endl;
      }
    else // if ( totalvolume > 500 &&  totalmass/totalct > 1/500 )  {
      {
//      // std::cout << " Volume Of Label " << *it << " is " << totalvolume <<   "  Avg-Location " << myCenterOfMass
//                << " mass is " << totalmass << " average-val is " << totalmass / totalct << std::endl;
//      //      // std::cout << *it << "  " <<  totalvolume <<  " & " <<  totalmass/totalct   << " \ " << std::endl;
      }

// square image
    squareimage->GetBufferPointer()[labelcount] = totalmass / totalct;
    if( ImageDimension == 2 )
      {
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << ",0,0," << currentlabel << "," << totalmass << "," << totalvolume << "," << totalct << std::endl;
      }
    if( ImageDimension == 3 )
      {
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << "," << myCenterOfMass[2] << ",0," << currentlabel << "," << totalmass << "," << totalvolume << "," << totalct << std::endl;
      }
    if( ImageDimension == 4 )
      {
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << "," << myCenterOfMass[2] << ","
              << myCenterOfMass[3] << "," << currentlabel << "," << totalmass << "," << totalvolume << "," << totalct << std::endl;
      }
    labelcount++;
    }

  logfile.close();

  //WriteImage<TwoDImageType>(squareimage, imagename.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int LabelThickness(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return EXIT_FAILURE;
    }
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer eimage = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  ReadImage<ImageType>(eimage, fn1.c_str() );
  // typename TImage::Pointer  Morphological( typename TImage::Pointer input, float rad, unsigned int option,
  //                                       float dilateval)
  eimage = ants::Morphological<ImageType>(image, 1, 4, 1);
  typedef float                  PixelType;
  typedef std::vector<PixelType> LabelSetType;
  LabelSetType myLabelSet;
  /** count the labels in the image */
  unsigned long maxlab = 0;
  Iterator      iIt( image, image->GetLargestPossibleRegion() );
  for( iIt.GoToBegin(); !iIt.IsAtEnd(); ++iIt )
    {
    PixelType label = iIt.Get();
    if( fabs(label) > 0 )
      {
      if( find( myLabelSet.begin(), myLabelSet.end(), label )
          == myLabelSet.end() )
        {
        myLabelSet.push_back( label );
        }
      if( label > maxlab )
        {
        maxlab = (unsigned long)label;
        }
      }
    }

  std::sort(myLabelSet.begin(), myLabelSet.end() );
  typename LabelSetType::const_iterator it;
  unsigned long labelcount = 0;
  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
    {
    labelcount++;
    }

  typename ImageType::SpacingType spacing = image->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= static_cast<float>( spacing[i] );
    }
  volumeelement = std::pow( static_cast<double>( volumeelement ), static_cast<double>( 0.3333 ) );

  vnl_vector<double> surface(maxlab + 1, 0);
  vnl_vector<double> volume(maxlab + 1, 0);

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  rad.Fill( 1 );
  iteratorType GHood(rad, image, image->GetLargestPossibleRegion() );
  float        Gsz = (float)GHood.Size();
  // iterate over the label image and index into the stat image
  for( iIt.GoToBegin(); !iIt.IsAtEnd(); ++iIt )
    {
    PixelType label = iIt.Get();
    if(  label > 0 )
      {
      volume[(unsigned long) label] = volume[(unsigned long) label] + 1;
      }
    label = label - eimage->GetPixel( iIt.GetIndex() );
    if(  label > 0 )
      {
      GHood.SetLocation( iIt.GetIndex() );
      bool isedge = false;
      for( unsigned int ii = 0; ii < Gsz; ii++ )
        {
        if( itk::Math::FloatAlmostEqual( GHood.GetPixel( ii ), 0.0f ) )
          {
          isedge = true;
          }
        }
      if( isedge )
        {
        surface[(unsigned long) label] = surface[(unsigned long) label] + 1;
        }
      }
    }
  for(double i : surface)
    {
    if( i > 0 )
      {
      // std::cout << " S " << surface[i] << " V " << volume[i] << " T " << volume[i] / surface[i] * 2.0 <<  std::endl;
      }
    }
  for( iIt.GoToBegin(); !iIt.IsAtEnd(); ++iIt )
    {
    PixelType label = iIt.Get();
    if( label > 0 )
      {
      PixelType thicknessprior = volumeelement;
      if( surface[label] > 0 )
        {
        thicknessprior = static_cast<float>( volume[label] / surface[label] ) * 2.0f * volumeelement;
        }
      eimage->SetPixel( iIt.GetIndex(), thicknessprior );
      }
    }
  WriteImage<ImageType>( eimage, outname.c_str() );
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int LabelThickness2( int argc, char *argv[] )
{
  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef float                                 RealType;
  typedef itk::Image<RealType, ImageDimension>  RealImageType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return EXIT_FAILURE;
    }
  int               argct = 2;
  const std::string outname = std::string( argv[argct] );
  argct += 2;

  std::string fn = std::string( argv[argct++] );

  typename LabelImageType::Pointer labelImage = nullptr;
  ReadImage<LabelImageType>( labelImage, fn.c_str() );

  typename LabelImageType::SpacingType spacing = labelImage->GetSpacing();
  float volumeElement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeElement *= static_cast<float>( spacing[i] );
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, RealImageType> GeometryFilterType;
  typename GeometryFilterType::Pointer geometryFilter = GeometryFilterType::New();
  geometryFilter->SetInput( labelImage );
  geometryFilter->CalculatePixelIndicesOff();
  geometryFilter->CalculateOrientedBoundingBoxOff();
  geometryFilter->CalculateOrientedLabelRegionsOff();
  geometryFilter->Update();

  typedef itk::LabelPerimeterEstimationCalculator<LabelImageType> AreaFilterType;
  typename AreaFilterType::Pointer areaFilter = AreaFilterType::New();
  areaFilter->SetImage( labelImage );
  areaFilter->SetFullyConnected( true );
  areaFilter->Compute();

  typename GeometryFilterType::LabelsType allLabels = geometryFilter->GetLabels();
  typename GeometryFilterType::LabelsType::iterator allLabelsIt;

  std::sort( allLabels.begin(), allLabels.end() );
  for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
    {
    if( *allLabelsIt == 0 )
      {
      continue;
      }
    //    RealType volume = geometryFilter->GetVolume( *allLabelsIt ) * volumeElement;
    // RealType perimeter = areaFilter->GetPerimeter( *allLabelsIt );
    // RealType thicknessPrior = volume / perimeter;
    // std::cout << "Label "  << *allLabelsIt << ": ";
    // std::cout << "volume = " << volume << ", ";
    // std::cout << "area = " << perimeter << ", ";
    // std::cout << "thickness = " << thicknessPrior << std::endl;
    }

  typename RealImageType::Pointer thicknessPriorImage = RealImageType::New();
  thicknessPriorImage->CopyInformation( labelImage );
  thicknessPriorImage->SetRegions( labelImage->GetLargestPossibleRegion() );
  thicknessPriorImage->Allocate();
  thicknessPriorImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<LabelImageType> It( labelImage, labelImage->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    LabelType label = It.Get();
    if( label == 0 )
      {
      continue;
      }
    RealType volume = geometryFilter->GetVolume( label ) * volumeElement;
    RealType perimeter = areaFilter->GetPerimeter( label );

    RealType thicknessPrior = volume / perimeter;
    thicknessPriorImage->SetPixel( It.GetIndex(), thicknessPrior );
    }

  WriteImage<RealImageType>( thicknessPriorImage, outname.c_str() );
  return EXIT_SUCCESS;
}

// now words can be accessed like this WordList[n]; where 'n' is the index
template <unsigned int ImageDimension>
int ROIStatistics(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  //  if(grade_list.find("Tim") == grade_list.end()) {  // std::cout<<"Tim is not in the map!"<<endl; }
  // mymap.find('a')->second
  int argct = 2;
  if( argc < 6 )
    {
    std::cout << " not enough parameters --- usage example 1 :" << "" << std::endl;
    std::cout << argv[0]
              << " ImageMath  3 output.csv ROIStatistics roinames.txt LabelImage.nii.gz ValueImage.nii.gz  "
              << std::endl;
    throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  typedef vnl_matrix<double> MatrixType;
  argct += 2;
  std::string fn0 = std::string(argv[argct]);   argct++;
  // std::cout << "  fn0 " << fn0 << std::endl;
  std::map<unsigned int, std::string> roimap = RoiList(fn0);
  std::string                         fn1 = std::string(argv[argct]);   argct++;
  std::string                         fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  std::string fn3 = "";
  if(  argc > argct )
    {
    fn3 = std::string(argv[argct]);   argct++;
    }
  std::string fn4 = "";
  if(  argc > argct )
    {
    fn4 = std::string(argv[argct]);   argct++;
    }
  std::string fn5 = "";
  if(  argc > argct )
    {
    fn5 = std::string(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image = nullptr;
  typename ImageType::Pointer valimage = nullptr;
  typename ImageType::Pointer valimage3 = nullptr;
  typename ImageType::Pointer valimage4 = nullptr;
  typename ImageType::Pointer valimage5 = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(valimage, fn2.c_str() );
    }
  if( fn3.length() > 3 )
    {
    ReadImage<ImageType>(valimage3, fn3.c_str() );
    }
  if( fn4.length() > 3 )
    {
    ReadImage<ImageType>(valimage4, fn4.c_str() );
    }
  if( fn5.length() > 3 )
    {
    ReadImage<ImageType>(valimage5, fn5.c_str() );
    }

  typedef float                  PixelType;
  typedef std::vector<PixelType> LabelSetType;
  LabelSetType myLabelSet;
/** count the labels in the image */
  unsigned long maxlab = 0;
  Iterator      It( image, image->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PixelType label = It.Get();
    if( fabs(label) > 0 )
      {
      if( find( myLabelSet.begin(), myLabelSet.end(), label )
          == myLabelSet.end() )
        {
        myLabelSet.push_back( label );
        }
      if( label > maxlab )
        {
        maxlab = (unsigned long)label;
        }
      }
    }
  // define the output matrix
  // for a csv file of  n_regions for rows
  // with columns of   Volume,CenterOfMassX,CenterOfMassY,CenterOfMassZ,CenterOfMassTime,
  MatrixType mCSVMatrix(maxlab, ImageDimension + 1, 0);

  //  maxlab=32; // for cortical analysis
  // compute the voxel volume
  typename ImageType::SpacingType spacing = image->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= static_cast<float>( spacing[i] );
    }

  vnl_vector<double> pvals(maxlab + 1, 1.0);
  vnl_vector<double> pvals3(maxlab + 1, 1.0);
  vnl_vector<double> pvals4(maxlab + 1, 1.0);
  vnl_vector<double> pvals5(maxlab + 1, 1.0);
  vnl_vector<double> clusters(maxlab + 1, 0.0);
  vnl_vector<double> masses(maxlab + 1, 0.0);
  //  typename ImageType::PointType mycomlist[33];
  std::ofstream logfile;
  logfile.open(outname.c_str() );

  std::sort(myLabelSet.begin(), myLabelSet.end() );
  typename LabelSetType::const_iterator it;
  unsigned long labelcount = 0;
  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
    {
    labelcount++;
    }

  typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer interp =  ScalarInterpolatorType::New();
  if( valimage )
    {
    interp->SetInputImage(valimage);
    }
  else
    {
    interp->SetInputImage(image);
    }
  // iterate over the label image and index into the stat image
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PixelType label = It.Get();
    typename ImageType::PointType point;
    image->TransformIndexToPhysicalPoint(It.GetIndex(), point);
    float vv = interp->Evaluate( point );
    //    float vv = valimage->GetPixel(It.GetIndex());
    if(  label > 0  && fabs(vv) > 1.e-9 )
      {
      clusters[(unsigned long) label] += 1;
      masses[(unsigned long) label] += static_cast<double>( vv );
      }
    }
  logfile << "ROIName,ROINumber,ClusterSize,Mass,Mean,comX,comY,comZ,comT" << std::endl;
  //  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
  for( unsigned int mylabel = 1; mylabel < maxlab + 1; mylabel++ )
    {
    unsigned int roi = mylabel - 1;
    /* first count which roi it is by iterating through the label sets
     it = myLabelSet.begin();
    while (  (*it) != mylabel && it != myLabelSet.end() )
  {
    // std::cout << " it " << *it << " roi " << roi << " mylabel " << mylabel << std::endl;
    roi++;
    ++it;
    }*/

    typename ImageType::PointType myCenterOfMass;
    myCenterOfMass.Fill(0);
    unsigned long totalct = 0;
    if( clusters[mylabel] > 0 )
      {
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        PixelType label = It.Get();
        if( itk::Math::FloatAlmostEqual( label, static_cast<PixelType>( mylabel ) ) )
          {
          totalct += 1;
          // compute center of mass
          typename ImageType::PointType point;
          image->TransformIndexToPhysicalPoint(It.GetIndex(), point);
          for( unsigned int i = 0; i < spacing.Size(); i++ )
            {
            myCenterOfMass[i] += point[i];
            }
          }
        }
      } // if clusters big enough
    if( totalct == 0 )
      {
      totalct = 1;
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= static_cast<double>( totalct );
      }
    double comy = 0, comz = 0, comt = 0;
    double comx = myCenterOfMass[0];
    if( ImageDimension >= 2 )
      {
      comy = myCenterOfMass[1];
      }
    if( ImageDimension >= 3 )
      {
      comz = myCenterOfMass[2];
      }
    if( ImageDimension == 4 )
      {
      comt = myCenterOfMass[3];
      }
    if( roi < roimap.size() )
      {
      logfile << roimap[roi] << "," << roi + 1 << ","
              << clusters[roi
                  + 1] << ","
              << masses[roi
                + 1] << ","
              << masses[roi
                + 1] / clusters[roi + 1] << "," << comx << "," << comy << "," << comz << "," << comt << std::endl;
      }
    }
  logfile.close();

  return 0;

  float        temp = sqrt( (float)labelcount);
  unsigned int squareimagesize = (unsigned int)(temp + 1);

  typedef itk::Image<float, 2> TwoDImageType;
  typename TwoDImageType::RegionType newregion;
  typename TwoDImageType::SizeType size;
  size[0] = size[1] = squareimagesize;
  newregion.SetSize(size);
  typename TwoDImageType::Pointer squareimage =
    AllocImage<TwoDImageType>(newregion, 0);

  labelcount = 0;
  typename ImageType::PointType myCenterOfMass;
  myCenterOfMass.Fill(0);
  //  for (unsigned int i=0; i<=maxlab; i++) mycomlist[i]=myCenterOfMass;
  for( it = myLabelSet.begin(); it != myLabelSet.end(); ++it )
    {
    float currentlabel = *it;
    float totalvolume = 0;
    float totalmass = 0;
    float totalct = 0;
    float maxoneminuspval = 0.0;
    float maxoneminuspval3 = 0.0;
    float maxoneminuspval4 = 0.0;
    float maxoneminuspval5 = 0.0;
    myCenterOfMass.Fill(0);
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( itk::Math::FloatAlmostEqual( label, static_cast<PixelType>( currentlabel ) ) )
        {
        totalvolume += volumeelement;
        totalct += 1;
        float vv = valimage->GetPixel(It.GetIndex() );
        if( valimage )
          {
          totalmass += vv;
          }
        if(  vv > maxoneminuspval )
          {
          maxoneminuspval = vv;
          }
        if( valimage3 )
          {
          vv = valimage3->GetPixel(It.GetIndex() );
          if(  fabs(vv) > fabs(maxoneminuspval3) )
            {
            maxoneminuspval3 = vv;
            }
          }
        if( valimage4 )
          {
          vv = valimage4->GetPixel(It.GetIndex() );
          if(  vv > maxoneminuspval4 )
            {
            maxoneminuspval4 = vv;
            }
          }
        if( valimage5 )
          {
          vv = valimage5->GetPixel(It.GetIndex() );
          if(  vv > maxoneminuspval5 )
            {
            maxoneminuspval5 = vv;
            }
          }

        // compute center of mass
        typename ImageType::PointType point;
        image->TransformIndexToPhysicalPoint(It.GetIndex(), point);
        for( unsigned int i = 0; i < spacing.Size(); i++ )
          {
          myCenterOfMass[i] += point[i];
          }
        }
      }
    if( itk::Math::FloatAlmostEqual( totalct, 0.0f ) )
      {
      totalct = 1.0f;
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= static_cast<double>( totalct );
      }

    clusters[(unsigned long)*it] = totalvolume;
    masses[(unsigned long)*it] = totalmass / totalct;
    pvals[(unsigned long)*it] = 1.0f - maxoneminuspval;
    pvals3[(unsigned long)*it] = maxoneminuspval3;
    pvals4[(unsigned long)*it] = 1.0f - maxoneminuspval4;
    pvals5[(unsigned long)*it] = 1.0f - maxoneminuspval5;
    // mycomlist[(unsigned long)*it]=myCenterOfMass;
// square image
    squareimage->GetBufferPointer()[labelcount] = totalmass / totalct;
    labelcount++;
    }
  for( unsigned int roi = 1; roi <= maxlab; roi++ )
    {
    logfile << roimap.find(roi)->second << ",";
    }
  for( unsigned int roi = maxlab + 1; roi <= maxlab * 2; roi++ )
    {
    if( roi < maxlab * 2 )
      {
      logfile << roimap.find(roi - maxlab)->second + std::string("mass") << ",";
      }
    else
      {
      logfile << roimap.find(roi - maxlab)->second + std::string("mass") << std::endl;
      }
    }
  for( unsigned int roi = 1; roi <= maxlab; roi++ )
    {
    logfile << clusters[roi] << ",";
    }
  for( unsigned int roi = maxlab + 1; roi <= maxlab * 2; roi++ )
    {
    if( roi < maxlab * 2 )
      {
      logfile << masses[roi - maxlab] << ",";
      }
    else
      {
      logfile << masses[roi - maxlab] << std::endl;
      }
    }
  logfile.close();
  return 0;
}

template <unsigned int ImageDimension>
int ByteImage(      int /*argc */, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);

  typename ImageType::Pointer image = nullptr;
  typename ByteImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );

  typedef itk::RescaleIntensityImageFilter<ImageType, ByteImageType> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( image );

  WriteImage<ByteImageType>( rescaler->GetOutput(), outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int PValueImage(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int dof = 1;
  if( argc > argct )
    {
    dof = std::stoi(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );

  // std::cout << " read Image" << fn1 << " dof " << dof << std::endl;
  typedef itk::Statistics::TDistribution DistributionType;
  typename DistributionType::Pointer distributionFunction = DistributionType::New();
  distributionFunction->SetDegreesOfFreedom(  dof );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( image,  image->GetLargestPossibleRegion() );
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    float val = image->GetPixel(vfIter.GetIndex() );
    if( !std::isnan(val) && !std::isinf(val) )
      {
      if( fabs( val ) > 0 )
        {
        float pval = 0;
        if( val >= 1 )
          {
          pval = 1.0 - distributionFunction->EvaluateCDF( val, dof );
          }
        else if( val < 1 )
          {
          pval = distributionFunction->EvaluateCDF( val, dof );
          }
        image->SetPixel(vfIter.GetIndex(), pval );
        }
      }
    else
      {
      image->SetPixel(vfIter.GetIndex(), 0);
      }
    }

  WriteImage<ImageType>(image, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int ConvertImageSetToMatrix(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;
    throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  std::string       ext = itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int rowcoloption = std::stoi(argv[argct]);   argct++;
  std::string  maskfn = std::string(argv[argct]); argct++;
  unsigned int numberofimages = 0;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  unsigned long voxct = 0;
  Iterator      mIter( mask, mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5f )
      {
      voxct++;
      }
    }

  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::SizeType size;
  size.Fill(0);
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        // std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  // std::cout << " largest image " << size << " num images " << numberofimages << " voxct " << voxct << std::endl;
  unsigned long xx1 = 0, yy1 = 0;
  if( rowcoloption == 0 )
    {
    // std::cout << " row option " << std::endl;
    xx1 = voxct;
    yy1 = numberofimages;
    }
  if( rowcoloption == 1 )
    {
    // std::cout << " col option " << std::endl;
    yy1 = voxct;
    xx1 = numberofimages;
    }
  unsigned long xsize = xx1;
  unsigned long ysize = yy1;

  if( strcmp(ext.c_str(), ".csv") == 0 )
    {
    typedef itk::Array2D<double> MatrixType;
    std::vector<std::string> ColumnHeaders;
    MatrixType               matrix(xsize, ysize);
    matrix.Fill(0);
    unsigned int imagecount = 0;
    for( unsigned int j = argct; j < argc; j++ )
      {
      std::string fn = std::string(argv[j]);
      ReadImage<ImageType>(image2, fn.c_str() );
      // std::cout << " image " << j << " is "  << fn << std::endl;
      unsigned long xx = 0, yy = 0, tvoxct = 0;
      if( rowcoloption == 0 )
        {
        yy = imagecount;
        }
      if( rowcoloption == 1 )
        {
        xx = imagecount;
        }
      for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
        {
        if( mIter.Get() >= 0.5f )
          {
          if( rowcoloption == 0 )
            {
            xx = tvoxct;
            }
          if( rowcoloption == 1 )
            {
            yy = tvoxct;
            }
          matrix[xx][yy] = image2->GetPixel(mIter.GetIndex() );
          std::string colname = std::string("V") + ants_to_string<unsigned int>(tvoxct);
          if( imagecount == 0 )
            {
            ColumnHeaders.push_back( colname );
            }
          tvoxct++;
          }
        }
      imagecount++;
      }

    // write out the array2D object
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outname );
    writer->SetInput( &matrix );
    if( rowcoloption == 0 )
      {
      writer->SetRowHeaders( ColumnHeaders );
      }
    if( rowcoloption == 1 )
      {
      writer->SetColumnHeaders( ColumnHeaders );
      }
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
/** declare the tiled image */
    typename MatrixImageType::SizeType tilesize;
    tilesize[0] = xsize;
    tilesize[1] = ysize;
    // std::cout << " allocate matrix " << tilesize << std::endl;
    typename MatrixImageType::RegionType region;
    region.SetSize( tilesize );

    typename MatrixImageType::Pointer matimage =
      AllocImage<MatrixImageType>(region);

    unsigned int imagecount = 0;
    for( unsigned int j = argct; j < argc; j++ )
      {
      std::string fn = std::string(argv[j]);
      ReadImage<ImageType>(image2, fn.c_str() );
      // std::cout << " image " << j << " is "  << fn << std::endl;
      unsigned long xx = 0, yy = 0, tvoxct = 0;
      if( rowcoloption == 0 )
        {
        yy = imagecount;
        }
      if( rowcoloption == 1 )
        {
        xx = imagecount;
        }
      typename MatrixImageType::IndexType mind;
      for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
        {
        if( mIter.Get() >= 0.5f )
          {
          if( rowcoloption == 0 )
            {
            xx = tvoxct;
            }
          if( rowcoloption == 1 )
            {
            yy = tvoxct;
            }
          mind[0] = xx;
          mind[1] = yy;
          //          // std::cout << " Mind " << mind << std::endl;
          matimage->SetPixel(mind, image2->GetPixel(mIter.GetIndex() ) );
          tvoxct++;
          }
        }
      imagecount++;
      }

    // std::cout << " mat size " << matimage->GetLargestPossibleRegion().GetSize() << std::endl;
    WriteImage<MatrixImageType>(matimage, outname.c_str() );
    }
  return 0;
}

template <unsigned int ImageDimension>
int RandomlySampleImageSetToCSV(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRandomConstIteratorWithIndex<ImageType>               Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  std::string       ext = itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int n_samples = std::stoi(argv[argct]);   argct++;
  /* std::string maskfn=std::string(argv[argct]); argct++;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask,maskfn.c_str());
  Iterator mIter( mask,mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    if (mIter.Get() >= 0.5) voxct++;
  */
  unsigned int numberofimages = 0;

  typename ImageType::Pointer image2 = nullptr;
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    }
  unsigned long xsize = numberofimages;
  unsigned long ysize = n_samples;

  if( strcmp(ext.c_str(), ".csv") == 0 )
    {
    typedef itk::Array2D<double> MatrixType;
    std::vector<std::string> ColumnHeaders;
    MatrixType               matrix(xsize, ysize);
    matrix.Fill(0);
    unsigned int imagecount = 0;
    for( unsigned int j = argct; j < argc; j++ )
      {
      std::string fn = std::string(argv[j]);
      ReadImage<ImageType>(image2, fn.c_str() );
      Iterator mIter( image2, image2->GetLargestPossibleRegion() );
      mIter.SetNumberOfSamples(n_samples);
      // std::cout << " image " << j << " is "  << fn << std::endl;
      unsigned long voxct = 0;
      for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
        {
        matrix[imagecount][voxct] = image2->GetPixel(mIter.GetIndex() );
        if( imagecount == 0 )
          {
          std::string colname = std::string("V") + ants_to_string<unsigned int>(voxct);
          ColumnHeaders.push_back( colname );
          }
        voxct++;
        }
      imagecount++;
      }

    // write out the array2D object
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outname );
    writer->SetInput( &matrix );
    writer->SetColumnHeaders( ColumnHeaders );
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    // std::cout << " need a csv file as output type , you tried " << outname << std::endl;
    }
  return 0;
}

template <typename T>
inline std::string to_string(const T& t)
{
  std::stringstream ss;

  ss << t;
  return ss.str();
}

template <unsigned int ImageDimension>
int ConvertImageSetToEigenvectors(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    // std::cout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  std::string       ext = std::string(".csv"); // itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int n_evecs = std::stoi(argv[argct]);   argct++;
  unsigned int rowcoloption = 1;
  std::string  maskfn = std::string(argv[argct]); argct++;
  unsigned int numberofimages = 0;
  typename ImageType::Pointer mask = nullptr;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  /** 1. compute max value in mask */
  unsigned long maxval = 0;
  for( Iterator mIter( mask, mask->GetLargestPossibleRegion() );
       !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() > maxval )
      {
      maxval = (unsigned long) mIter.Get();
      }
    }
  if( maxval == 0 )
    {
    // std::cout << " Max value in mask is <= 0, aborting. " << maxval << std::endl;
    throw std::exception();
    }

  typedef itk::Array2D<double> MatrixType;
  typename ImageType::Pointer image2 = nullptr;
  typename ImageType::SizeType size;
  size.Fill(0);
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);

        // std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  // std::cout << " largest image " << size << " num images " << numberofimages << std::endl;
  MatrixType avg_matrix(numberofimages, maxval);
  avg_matrix.Fill(0);
  for( unsigned long mv = 1; mv <= maxval; mv++ )
    {
    /** 2. count the voxels in this label */
    unsigned long voxct = 0;
    for( Iterator mIter( mask, mask->GetLargestPossibleRegion() );
         !mIter.IsAtEnd(); ++mIter )
      {
      if( itk::Math::FloatAlmostEqual( mIter.Get(), static_cast<float>( mv ) ) )
        {
        voxct++;
        }
      }

    unsigned long xx1 = 0, yy1 = 0;
    if( rowcoloption == 0 )
      {
      // std::cout << " row option " << std::endl;  xx1 = voxct;  yy1 = numberofimages;
      }
    if( rowcoloption == 1 )
      {
      // std::cout << " col option " << std::endl;  yy1 = voxct;  xx1 = numberofimages;
      }
    unsigned long xsize = xx1;
    unsigned long ysize = yy1;

    if( strcmp(ext.c_str(), ".csv") == 0 )
      {
      MatrixType matrix(xsize, ysize);
      matrix.Fill(0);
      unsigned int imagecount = 0;
      for( unsigned int j = argct; j < argc; j++ )
        {
        std::string fn = std::string(argv[j]);
        ReadImage<ImageType>(image2, fn.c_str() );
        // std::cout << " image " << j << " is "  << fn << std::endl;
        unsigned long xx = 0, yy = 0, tvoxct = 0;
        if( rowcoloption == 0 )
          {
          yy = imagecount;
          }
        if( rowcoloption == 1 )
          {
          xx = imagecount;
          }
        for( Iterator mIter( mask, mask->GetLargestPossibleRegion() );
             !mIter.IsAtEnd(); ++mIter )
          {
          if( itk::Math::FloatAlmostEqual( mIter.Get(), static_cast<PixelType>( mv ) ) )
            {
            if( rowcoloption == 0 )
              {
              xx = tvoxct;
              }
            if( rowcoloption == 1 )
              {
              yy = tvoxct;
              }
            matrix[xx][yy] = image2->GetPixel(mIter.GetIndex() );
            tvoxct++;
            }
          }
        imagecount++;
        }

      typedef double                                            Scalar;
      typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
      typename matrixOpType::Pointer matrixOps = matrixOpType::New();
      typedef vnl_vector<Scalar> ScalarVectorType;
      MatrixType evec_matrix(xsize, n_evecs + 1);
      evec_matrix.Fill(0);
      ScalarVectorType avg = matrixOps->AverageColumns(matrix);
      avg_matrix.set_column(mv, avg);
      evec_matrix.set_column(0, avg);
      MatrixType tempm = matrixOps->GetCovMatEigenvectors(matrix);
      for( unsigned int i = 0; i < n_evecs; i++ )
        {
        /** the GetCovMatEigenvector function normalizes the matrix & computes the covariance matrix internally */
        //      VectorType evec=matrixOps->GetCovMatEigenvector(matrix,i);
        evec_matrix.set_column(i + 1, tempm.get_column(i) );
        }

      // write out the array2D object
      typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      std::string         num = std::string("evecs_for_label") + to_string<unsigned long>(mv);
      std::string         eoutname = outname + num + ext;
      writer->SetFileName( eoutname );
      writer->SetInput( &evec_matrix );
      try
        {
        writer->Write();
        }
      catch( itk::ExceptionObject& /* exp */)
        {
        // std::cout << "Exception caught!" << std::endl;
        // std::cout << exp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      // std::cout << " can only write out csv files " << std::endl;
      }
    } // end loop over mv variable

    {
    // write out the array2D object
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    std::string         num = std::string("avg_for_labels");
    std::string         eoutname = outname + num + ext;
    writer->SetFileName( eoutname );
    writer->SetInput( &avg_matrix );
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& /* exp */)
      {
      // std::cout << "Exception caught!" << std::endl;
      // std::cout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }

  return 0;
}

template <unsigned int ImageDimension>
int ConvertImageToFile(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, fn1.c_str() );
  typename ImageType::Pointer mask = nullptr;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(mask, fn2.c_str() );
    }

  // std::cout << " read Image" << fn1 << " mask? " << fn2 << std::endl;
  std::ofstream logfile;
  logfile.open(outname.c_str() );
  if( logfile.good() )
    {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
    ImageIterator vfIter( image,  image->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      bool getval = true;
      if( mask->GetPixel(vfIter.GetIndex() ) < 0.5f )
        {
        getval = false;
        }
      if( getval )
        {
        float val = image->GetPixel(vfIter.GetIndex() );
        logfile << " " << val;
        }
      }
    logfile << std::endl;
    }
  logfile.close();

  return 0;
}

template <unsigned int ImageDimension>
int CorrelationUpdate(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]); argct++;
    }
  else
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }
  unsigned int radius = 2;
  if( argc > argct )
    {
    radius = std::stoi(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::Pointer imageout = nullptr;
  ReadImage<ImageType>(imageout, fn1.c_str() );
  imageout->FillBuffer(0);
  typename ImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image2, fn2.c_str() );

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = radius;
    }
  iteratorType GHood(rad, image1, image1->GetLargestPossibleRegion() );
  float        Gsz = (float)GHood.Size();
  GHood.GoToBegin();
  while( !GHood.IsAtEnd() )
    {
    typename ImageType::IndexType ind = GHood.GetIndex();
    bool isinside = true;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      float shifted = ind[j];
      if( shifted < (radius + 1) || shifted >  image1->GetLargestPossibleRegion().GetSize()[j] - radius - 1  )
        {
        isinside = false;
        }
      }
    if( isinside )
      {
      if( image1->GetPixel(ind) > 0 || image2->GetPixel(ind) > 0 )
        {
        // compute mean difference
        float diff = 0.0;
        for( unsigned int i = 0; i < GHood.Size(); i++ )
          {
          const typename ImageType::IndexType & ind2 = GHood.GetIndex(i);
          diff += (image1->GetPixel(ind2) - image2->GetPixel(ind2) );
          }
        diff /= Gsz;
        float upd = (image1->GetPixel(ind) - image2->GetPixel(ind) ) - diff;
        imageout->SetPixel(ind, upd);
        }
      }
    ++GHood;
    }

  WriteImage<ImageType>(imageout, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int MajorityVoting( int argc, char *argv[] )
{
  typedef int                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>         ImageType;
  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  IteratorType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );

  // Check for usage 1, e.g.,
  //     ImageMath 3 outputLabels.nii.gz MajorityVoting labelImage1.nii.gz ... labelImageN.nii.gz
  //  vs.
  // usage 2 (hidden from casual user), e.g.,
  //     ImageMath 3 outputLabels.nii.gz MajorityVoting 0.9 labelImage1.nii.gz ... labelImageN.nii.gz

  typename ImageType::Pointer image4 = nullptr;
  bool doUsage1 = ReadImage<ImageType>( image4, argv[4] );

  unsigned long numberOfImages = 0;
  unsigned int startIndex = 0;
  float percentage = 0.0;

  if( doUsage1 )
    {
    numberOfImages = argc - 4;
    startIndex = 4;
    }
  else
    {
    numberOfImages = argc - 5;
    startIndex = 5;
    percentage = atof( argv[4] );
    }

  if( numberOfImages <= 1 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  typename std::vector<typename ImageType::Pointer> images;
  images.resize( numberOfImages );

  for( int i = startIndex; i < argc; i++ )
    {
    if( i == 4 )
      {
      images[0] = image4;
      }
    else
      {
      ReadImage<ImageType>( images[i - startIndex], argv[i] );
      }
    }

  typename ImageType::Pointer output = nullptr;

  // Find maximum label
  int maxLabel = 0;
  for( unsigned int i = 0; i < numberOfImages; i++ )
    {
    typename CalculatorType::Pointer calc = CalculatorType::New();
    calc->SetImage( images[i] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }
  unsigned long nLabels = maxLabel + 1; // account for label=0

  output = AllocImage<ImageType>( images[0], 0 );

  typename ImageType::Pointer outputMask = nullptr;
  outputMask = AllocImage<ImageType>( images[0], 0 );

  IteratorType it( output, output->GetLargestPossibleRegion() );
  IteratorType itM( outputMask, outputMask->GetLargestPossibleRegion() );

  itk::Array<unsigned long> votes;
  votes.SetSize( nLabels );

  while( !it.IsAtEnd() )
    {
    votes.Fill( 0 );
    unsigned long maxVotes = 0;
    unsigned long votedLabel = 0;
    for( unsigned long i = 0; i < numberOfImages; i++ )
      {
      unsigned long label = images[i]->GetPixel( it.GetIndex() );
      votes.SetElement( label, votes.GetElement( label ) + 1 );

      if( votes.GetElement( label ) > maxVotes )
        {
        maxVotes = votes.GetElement( label );
        votedLabel = label;
        }
      }

    if( doUsage1 )
      {
      it.Set( votedLabel );
      }
    else
      {
      if( static_cast<float>( maxVotes ) / static_cast<float>( numberOfImages ) >= percentage )
        {
        it.Set( votedLabel );
        }
      else
        {
        itM.Set( 1 );
        }
      }

    ++it;
    ++itM;
    }

  WriteImage<ImageType>( output, outputName.c_str() );

  if( ! doUsage1 )
    {
    std::vector<std::string> fileExtensions;
    fileExtensions.emplace_back(".nii" );
    fileExtensions.emplace_back(".mha" );
    fileExtensions.emplace_back(".nrrd" );

    std::string outputMaskName;
    for(auto & fileExtension : fileExtensions)
      {
      if( outputName.find_last_of( fileExtension ) != std::string::npos )
        {
        outputMaskName = outputName.insert(
          outputName.find_last_of( fileExtension ) - fileExtension.length(), "_Mask" );
        break;
        }
      }
    if( outputMaskName.empty() )
      {
      // std::cout << " Unrecognized output file extension. " << std::endl;
      return 1;
      }
    WriteImage<ImageType>( outputMask, outputMaskName.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int MostLikely( int argc, char *argv[] )
{
  typedef float                                               PixelType;
  typedef itk::Image<PixelType, ImageDimension>               ImageType;
  typedef itk::Image<int, ImageDimension>                     LabeledImageType;
  typedef itk::ImageRegionIteratorWithIndex<LabeledImageType> IteratorType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );
  float       sigma = atof( argv[4] );

  typename LabeledImageType::Pointer output = LabeledImageType::New();
  typename ImageType::Pointer prob = ImageType::New();
  // Read input segmentations
  for( int i = 5; i < argc; i++ )
    {
    typename ImageType::Pointer iLabel = ImageType::New();
    ReadImage<ImageType>( iLabel, argv[i] );

    if( i == 5 )
      {
      output->SetRegions( iLabel->GetLargestPossibleRegion() );
      output->SetSpacing( iLabel->GetSpacing() );
      output->SetOrigin( iLabel->GetOrigin() );
      output->SetDirection( iLabel->GetDirection() );
      output->Allocate();
      output->FillBuffer( 0 );

      prob->SetRegions( iLabel->GetLargestPossibleRegion() );
      prob->SetSpacing( iLabel->GetSpacing() );
      prob->SetOrigin( iLabel->GetOrigin() );
      prob->SetDirection( iLabel->GetDirection() );
      prob->Allocate();
      prob->FillBuffer( 0 );
      }

    IteratorType it( output, output->GetLargestPossibleRegion() );

    while( !it.IsAtEnd() )
      {
      if( ( iLabel->GetPixel( it.GetIndex() ) > sigma ) &&
          ( iLabel->GetPixel( it.GetIndex() ) > prob->GetPixel( it.GetIndex() ) ) )
        {
        prob->SetPixel( it.GetIndex(), iLabel->GetPixel( it.GetIndex() ) );
        output->SetPixel( it.GetIndex(), i - 4 );
        }

      ++it;
      }
    }

  WriteImage<LabeledImageType>( output, outputName.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int STAPLE( int argc, char *argv[] )
{
  typedef int                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension>     OutputImageType;

  typedef itk::MinimumMaximumImageCalculator<ImageType>      CalculatorType;
  typedef itk::STAPLEImageFilter<ImageType, OutputImageType> StapleFilterType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );
  typename StapleFilterType::Pointer stapler = StapleFilterType::New();
  std::string::size_type idx;
  idx = outputName.find_first_of('.');
  std::string tempname = outputName.substr(0, idx);
  std::string extension = outputName.substr(idx, outputName.length() );
  float       confidence = atof( argv[4] ); // = 0.5

  // stapler->SetForegroundValue( foreground );
  stapler->SetConfidenceWeight( confidence );

  // Read input segmentations
  typename ImageType::Pointer              nullimage( nullptr );
  std::vector<typename ImageType::Pointer> images(argc - 5, nullimage );
  typename CalculatorType::Pointer calc = CalculatorType::New();
  int maxLabel = 0;
  for( int i = 5; i < argc; i++ )
    {
    images[i - 5] = ImageType::New();
    ReadImage<ImageType>( images[i - 5], argv[i] );
    stapler->SetInput( i - 5, images[i - 5] );
    // std::cout << "Input image " << i - 5 << " " << argv[i] << std::endl;

    calc->SetImage( images[i - 5] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }

  // std::cout << "Examining " << maxLabel << " labels" << std::endl;
  for( int label = 1; label <= maxLabel; label++ )
    {
    char num[16];
    sprintf( num, "%04d", label ); //NOTE: %04d is between 4 and 10 bytes

    std::string oname = tempname + num + extension;
    stapler->SetForegroundValue( label );
    stapler->Update();
    WriteImage<OutputImageType>( stapler->GetOutput(), oname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int AverageLabels( int argc, char *argv[] )
{
  typedef int                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension>     OutputImageType;

  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  IteratorType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string            outputName = std::string( argv[2] );
  std::string::size_type idx;
  idx = outputName.find_first_of('.');
  std::string tempname = outputName.substr(0, idx);
  std::string extension = outputName.substr(idx, outputName.length() );

  // Read input segmentations
  typename ImageType::Pointer              nullimage( nullptr );
  std::vector<typename ImageType::Pointer> images(argc - 4, nullimage );

  typename CalculatorType::Pointer calc = CalculatorType::New();
  int maxLabel = 0;
  for( int i = 4; i < argc; i++ )
    {
    images[i - 4] = ImageType::New();
    ReadImage<ImageType>( images[i - 4], argv[i] );
    // std::cout << "Input image " << i - 4 << " " << argv[i] << std::endl;

    calc->SetImage( images[i - 4] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }

  // std::cout << "Examining " << maxLabel << " labels" << std::endl;
  typename OutputImageType::Pointer              nullout( nullptr );
  std::vector<typename OutputImageType::Pointer> outimages;
  for( int label = 0; label < maxLabel; label++ )
    {
    typename OutputImageType::Pointer img = OutputImageType::New();
    img->SetRegions( images[0]->GetLargestPossibleRegion() );
    img->SetDirection( images[0]->GetDirection() );
    img->SetOrigin( images[0]->GetOrigin() );
    img->SetSpacing( images[0]->GetSpacing() );
    img->Allocate();
    outimages.push_back( img );
    }
  for( unsigned int i = 0; i < images.size(); i++ )
    {
    IteratorType it( images[i], images[i]->GetLargestPossibleRegion() );
    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      if( it.Value() > 0 )
        {
        outimages[it.Value() - 1]->SetPixel( it.GetIndex(),
                                             outimages[it.Value() - 1]->GetPixel( it.GetIndex() )
                                             + 1.0f / static_cast<float>( images.size() ) );
        }
      }
    }
  for( int label = 1; label <= maxLabel; label++ )
    {
    char num[16];
    sprintf( num, "%04d", label ); //NOTE: %04d is between 4 and 10 bytes

    std::string oname = tempname + num + extension;
    WriteImage<OutputImageType>( outimages[label - 1], oname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int CorrelationVoting( int argc, char *argv[] )
{
  typedef float                                              PixelType;
  typedef int                                                LabelType;
  typedef itk::Image<PixelType, ImageDimension>              ImageType;
  typedef itk::Image<LabelType, ImageDimension>              LabelImageType;
  typedef itk::MinimumMaximumImageCalculator<LabelImageType> CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<LabelImageType>  IteratorType;
  typedef itk::NeighborhoodIterator<ImageType>               NeighborhoodIteratorType;

  if( argc < 6 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );

  // Read input images and segmentations
  const int nImages = (argc - 5) / 2;

  int radius = 5;
  if( argc > ( 5 + 2 * nImages ) )
    {
    radius = std::stoi( argv[5 + 2 * nImages] );
    }

  typename ImageType::Pointer target;
  ReadImage<ImageType>( target, argv[4] );

  typename std::vector<typename ImageType::Pointer>      images(nImages);
  typename std::vector<typename LabelImageType::Pointer> labels(nImages);
  for( int i = 5; i < (5 + nImages); i++ )
    {
    ReadImage<ImageType>( images[i - 5], argv[i] );
    }
  for( int i = 5 + nImages; i < (5 + 2 * nImages); i++ )
    {
    ReadImage<LabelImageType>( labels[i - 5 - nImages], argv[i] );
    }

  // Find maximum label
  int maxLabel = 0;
  for( int i = 0; i < nImages; i++ )
    {
    typename CalculatorType::Pointer calc = CalculatorType::New();
    calc->SetImage( labels[i] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }
  unsigned long nLabels = maxLabel + 1; // account for label=0

  typename LabelImageType::Pointer output =
    AllocImage<LabelImageType>(images[0], 0);

  IteratorType      it( output, output->GetLargestPossibleRegion() );
  itk::Array<float> votes;
  votes.SetSize( nLabels );

  typename NeighborhoodIteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = radius;
    }
  NeighborhoodIteratorType metricIt(rad, target, target->GetLargestPossibleRegion() );

  while( !it.IsAtEnd() )
    {
    votes.Fill(0);
    float maxVotes = 0;
    int   votedLabel = 0;
    for( int i = 0; i < nImages; i++ )
      {
      int label = labels[i]->GetPixel( it.GetIndex() );
      votes.SetElement(label, votes.GetElement(label) + 1 );

      if( votes.GetElement(label) > maxVotes )
        {
        maxVotes = votes.GetElement(label);
        votedLabel = label;
        }
      }

    // If all agree, assign label immediately
    if( static_cast<int>( maxVotes ) == nImages )
      {
      it.Set( votedLabel );
      }
    else
      {
      votes.Fill( 0.0 );
      maxVotes = 0.0;
      votedLabel = 0;

      itk::VariableLengthVector<float> weights( nImages );
      for( int i = 0; i < nImages; i++ )
        {
        metricIt.SetLocation( it.GetIndex() );
        float targetMean = 0.0;
        float targetVar = 0.0;
        float imageMean = 0.0;
        float imageVar = 0.0;
        float k = 0;
        float product = 0.0;
        for( unsigned int j = 0; j < metricIt.Size(); j++ )
          {
          typename NeighborhoodIteratorType::OffsetType internal;
          typename NeighborhoodIteratorType::OffsetType offset;
          if( metricIt.IndexInBounds( j, internal, offset ) )
            {
            k++;
            typename ImageType::IndexType idx = metricIt.GetIndex( j );
            if( itk::Math::FloatAlmostEqual( k, 1.0f ) )
              {
              targetMean = target->GetPixel( idx );
              imageMean = images[i]->GetPixel( idx );
              targetVar = 0.0;
              imageVar = 0.0;
              }
            else
              {
              float oldMean = targetMean;
              float value = target->GetPixel( idx );
              targetMean = targetMean + (value - targetMean) / k;
              targetVar = targetVar + (value - oldMean) * ( value - targetMean );

              oldMean = imageMean;
              float iValue = images[i]->GetPixel( idx );
              imageMean = imageMean + ( iValue - imageMean ) / k;
              imageVar = imageVar + ( iValue - oldMean) * ( iValue - imageMean );

              product += value * iValue;
              }
            } // metricIt.IndexInBounds()
          }   // j >= metricIt.Size()

        targetVar /= (k - 1);
        imageVar /= (k - 1);
        float pearson =
          ( product - k * targetMean * imageMean ) / ( (k - 1) * std::sqrt(targetVar) * std::sqrt(imageVar) );
        weights.SetElement( i, std::fabs(pearson) );
        } // i >= nImages
      for( int i = 0; i < nImages; i++ )
        {
        int label = labels[i]->GetPixel( it.GetIndex() );
        votes.SetElement(label, votes.GetElement(label) + weights.GetElement(i) );

        if( votes.GetElement(label) > maxVotes )
          {
          maxVotes = votes.GetElement(label);
          votedLabel = label;
          }
        }

      it.Set(votedLabel);
      } // else find weighted vote winner

    ++it;
    } // iterate over output image

  WriteImage<LabelImageType>( output, outputName.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int ImageMetrics( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  // typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4
  //  <ImageType, ImageType, ImageType>                         MetricType;

  if( argc < 5 )
    {
    // std::cout << "ERROR: Not enough inputs " << std::endl;
    return 1;
    }

  typename ImageType::Pointer img1;
  ReadImage<ImageType>( img1, argv[4] );

  typename ImageType::Pointer img2;
  ReadImage<ImageType>( img2, argv[5] );

  float value = 0.0;

  if( strcmp(argv[3], "NeighborhoodCorrelation") == 0 )
    {
    typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4
      <ImageType, ImageType, ImageType>                          MetricType;

    typedef typename itk::ImageMaskSpatialObject<ImageDimension> MaskType;
    typedef typename MaskType::ImageType                         MaskImageType;

    int r = 5;
    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );

    if( argc > 6 )
      {
      r = std::stoi( argv[6] );
      }
    if( argc > 7 )
      {
      typename MaskType::Pointer mask = MaskType::New();
      typename MaskImageType::Pointer maskImage = MaskImageType::New();
      ReadImage<MaskImageType>( maskImage, argv[7] );
      mask->SetImage( maskImage );
      metric->SetFixedImageMask( mask );
      metric->SetMovingImageMask( mask );
      }

    typename MetricType::RadiusType radius;
    radius.Fill( r );
    metric->SetRadius( radius );

    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "NormalizedCorrelation") == 0 )
    {
    typedef itk::CorrelationImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    typedef typename itk::ImageMaskSpatialObject<ImageDimension> MaskType;
    typedef typename MaskType::ImageType                         MaskImageType;

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );

    if( argc > 6 )
      {
      typename MaskType::Pointer mask = MaskType::New();
      typename MaskImageType::Pointer maskImage = MaskImageType::New();
      ReadImage<MaskImageType>( maskImage, argv[6] );
      mask->SetImage( maskImage );
      metric->SetFixedImageMask( mask );
      metric->SetMovingImageMask( mask );
      }

    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "Demons") == 0 )
    {
    typedef itk::DemonsImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    typedef typename itk::ImageMaskSpatialObject<ImageDimension> MaskType;
    typedef typename MaskType::ImageType                         MaskImageType;

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );

    if( argc > 6 )
      {
      typename MaskType::Pointer mask = MaskType::New();
      typename MaskImageType::Pointer maskImage = MaskImageType::New();
      ReadImage<MaskImageType>( maskImage, argv[6] );
      mask->SetImage( maskImage );
      metric->SetFixedImageMask( mask );
      metric->SetMovingImageMask( mask );
      }

    // FIXME - Calling initialize on demons causes seg fault
    // std::cout << "Demons is currently broken" << std::endl;
    return 1;

    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "Mattes") == 0 )
    {
    typedef itk::MattesMutualInformationImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    typedef typename itk::ImageMaskSpatialObject<ImageDimension> MaskType;
    typedef typename MaskType::ImageType                         MaskImageType;

    int bins = 32;
    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );

    if( argc > 6 )
      {
      bins = std::stoi( argv[6] );
      }

    if( argc > 7 )
      {
      typename MaskType::Pointer mask = MaskType::New();
      typename MaskImageType::Pointer maskImage = MaskImageType::New();
      ReadImage<MaskImageType>( maskImage, argv[7] );
      mask->SetImage( maskImage );
      metric->SetFixedImageMask( mask );
      metric->SetMovingImageMask( mask );
      }

    metric->SetNumberOfHistogramBins( bins );
    metric->Initialize();
    value = metric->GetValue();
    }

  std::cout << value << std::endl;

  return 0;
}

template <unsigned int ImageDimension>
int PearsonCorrelation( int argc, char *argv[] )
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

  if( argc < 5 )
    {
    // std::cout << "ERROR: Not enough inputs " << std::endl;
    return 1;
    }

  typename ImageType::Pointer img1;
  ReadImage<ImageType>( img1, argv[4] );

  typename ImageType::Pointer img2;
  ReadImage<ImageType>( img2, argv[5] );

  typename ImageType::Pointer mask;

  if( argc > 6 )
    {
    ReadImage<ImageType>( mask, argv[6] );
    }
  else
    {
    mask = AllocImage<ImageType>(img1, 1);
    }

  float xmean = 0.0;
  float ymean = 0.0;
  float xsquared = 0.0;
  float ysquared = 0.0;
  float product = 0.0;
  float k = 0;

  IteratorType it( mask, mask->GetLargestPossibleRegion() );
  while( !it.IsAtEnd() )
    {
    if( it.Value() > 0 )
      {
      k++;
      typename ImageType::IndexType idx = it.GetIndex();

      float value = img1->GetPixel( idx );
      xmean += value;
      xsquared += (value * value);

      value = img2->GetPixel( idx );
      ymean += value;
      ysquared += (value * value);

      product += img1->GetPixel( idx ) *  img2->GetPixel( idx );
      }
    ++it;
    }

  xmean /= k;
  ymean /= k;

  float pearson = ( product - k * xmean * ymean ) /
                  ( std::sqrt(xsquared - (k * xmean * xmean)) * std::sqrt(ysquared - (k * ymean * ymean)) );
  std::cout << pearson << std::endl;

  return 0;
}

template <unsigned int ImageDimension>
int Project( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  if( argc < 5 )
    {
    // std::cout << "Not enough input parameters" << std::endl;
    return EXIT_FAILURE;
    }
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int axis = std::stoi(argv[argct]);   argct++;
  unsigned int which = std::stoi(argv[argct]);   argct++;
  typename ImageType::Pointer imagein;
  typename ImageType::Pointer imageout;
  ReadImage<ImageType>( imagein, fn1.c_str() );
  if( which == 0 )
    {
    typedef itk::SumProjectionImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetProjectionDimension( axis );
    filter->SetInput( imagein );
    filter->Update();
    imageout = filter->GetOutput();
    WriteImage<ImageType>( imageout, outname.c_str() );
    }
  if( which == 1 )
    {
    typedef itk::MaximumProjectionImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetProjectionDimension( axis );
    filter->SetInput( imagein );
    filter->Update();
    imageout = filter->GetOutput();
    WriteImage<ImageType>( imageout, outname.c_str() );
    }
  if( which == 2 )
    {
    typedef itk::MinimumProjectionImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetProjectionDimension( axis );
    filter->SetInput( imagein );
    filter->Update();
    imageout = filter->GetOutput();
    WriteImage<ImageType>( imageout, outname.c_str() );
    }
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int Translate( int argc, char *argv[] )
{
  typedef float                                             PixelType;
  typedef itk::Image<PixelType, ImageDimension>             ImageType;
  typedef itk::ResampleImageFilter<ImageType, ImageType>    FilterType;
  typedef itk::TranslationTransform<double, ImageDimension> TransformType;
  typedef typename TransformType::OutputVectorType          OffsetType;

  if( argc < (int)(ImageDimension + 4 ) )
    {
    // std::cout << "Not enough input parameters" << std::endl;
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer input = ImageType::New();
  ReadImage<ImageType>( input, argv[4] );
  // std::cout << "Input image: " << argv[4] << std::endl;

  typename TransformType::Pointer translate = TransformType::New();
  OffsetType offset;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    offset[i] = atof( argv[i + 5] );
    }
  translate->SetOffset( offset );
  // std::cout << "Translation: " << offset << std::endl;

  typename FilterType::Pointer resample = FilterType::New();
  resample->SetInput( input );
  resample->SetTransform( translate );
  resample->SetOutputParametersFromImage( input );
  resample->SetDefaultPixelValue( 0 );
  resample->Update();

  WriteImage<ImageType>( resample->GetOutput(), argv[2] );

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int MinMaxMean( int argc, char *argv[] )
{
  typedef float                                         PixelType;
  typedef itk::Image<PixelType, ImageDimension>         ImageType;
  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::ImageMomentsCalculator<ImageType>        MomentsCalculatorType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }

  typename ImageType::Pointer image = ImageType::New();
  ReadImage<ImageType>( image, argv[4] );

  typename CalculatorType::Pointer calc = CalculatorType::New();
  calc->SetImage( image );
  calc->ComputeMaximum();
  calc->ComputeMinimum();

  typename MomentsCalculatorType::Pointer calc2 = MomentsCalculatorType::New();
  calc2->SetImage( image );
  calc2->Compute();

  float mean = calc2->GetTotalMass();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    mean /= image->GetLargestPossibleRegion().GetSize()[i];
    }

  // std::cout << calc->GetMinimum() << " " << calc->GetMaximum() << " " << mean << std::endl;

  return 0;
}

template <unsigned int ImageDimension, typename TRealType, typename TImageType, typename TGImageType, typename TInterp>
TRealType PatchCorrelation(  itk::NeighborhoodIterator<TImageType> GHood,  itk::NeighborhoodIterator<TImageType> GHood2,
                             std::vector<unsigned int> activeindex, std::vector<TRealType> weight,
                             typename TGImageType::Pointer gimage,
                             typename TGImageType::Pointer gimage2,
                             TInterp interp2 )
{
  typedef TRealType                                      RealType;
  typedef typename TImageType::PointType                 PointType;
  typedef itk::CovariantVector<RealType, ImageDimension> GradientPixelType;
  typedef vnl_vector<RealType>                           VectorType;
  typedef typename TImageType::IndexType                 IndexType;
  unsigned int           Gsz = activeindex.size();
  std::vector<PointType> imagepatch1;
  std::vector<PointType> imagepatch2;
  VectorType             sample1( Gsz, 0);
  VectorType             sample2( Gsz, 0);
  vnl_matrix<RealType>   grad1mat( Gsz, ImageDimension); grad1mat.fill( 0 );
  vnl_matrix<RealType>   grad2mat( Gsz, ImageDimension); grad2mat.fill( 0 );
  // compute mean difference
  PointType avgpoint1;
  PointType avgpoint2;
  avgpoint1.Fill( 0 );
  avgpoint2.Fill( 0 );
  RealType wt = itk::NumericTraits<RealType>::OneValue() / static_cast<RealType>( Gsz );
  for( unsigned int ii = 0; ii < Gsz; ii++ )
    {
    sample1[ii] = GHood.GetPixel( activeindex[ii] );
    sample2[ii] = GHood2.GetPixel(  activeindex[ii] );
    IndexType gind = GHood.GetIndex(  activeindex[ii] );
    IndexType gind2 = GHood2.GetIndex(  activeindex[ii] );
    if( ( IsInside<TGImageType>( gimage, gind ) ) && (  IsInside<TGImageType>( gimage2, gind2 )  )  )
      {
      GradientPixelType grad1 = gimage->GetPixel( gind ) * weight[ii];
      GradientPixelType grad2 = gimage2->GetPixel( gind2 ) * weight[ii];
      for( unsigned int jj = 0; jj < ImageDimension; jj++ )
        {
        grad1mat( ii, jj ) = grad1[jj]; grad2mat( ii, jj ) = grad2[jj];
        }
      PointType point1;
      PointType point2;
      gimage->TransformIndexToPhysicalPoint(   gind, point1 );
      gimage2->TransformIndexToPhysicalPoint( gind2, point2 );
      for( unsigned int dd = 0; dd < ImageDimension; dd++ )
        {
        avgpoint1[dd] = avgpoint1[dd] + point1[dd] * static_cast<double>( wt );
        avgpoint2[dd] = avgpoint2[dd] + point2[dd] * static_cast<double>( wt );
        }
      imagepatch1.push_back( point1 );
      imagepatch2.push_back( point2 );
      }
    else
      {
      return 0;
      }
    }
  RealType mean1 = sample1.mean();
  RealType mean2 = sample2.mean();
  sample1 = ( sample1 - mean1 );
  sample2 = ( sample2 - mean2 );
  RealType sd1 = sqrt( sample1.squared_magnitude() );
  RealType sd2 = sqrt( sample2.squared_magnitude() );
  RealType correlation = inner_product( sample1, sample2 ) / ( sd1 * sd2 );

  bool ok = true;
  /** compute patch orientation */
  vnl_matrix<RealType>                cov1 = grad1mat.transpose() * grad1mat;
  vnl_matrix<RealType>                cov2 = grad2mat.transpose() * grad2mat;
  vnl_symmetric_eigensystem<RealType> eig1( cov1 );
  vnl_symmetric_eigensystem<RealType> eig2( cov2 );
  unsigned int                        eigind0 = 1; // biggest eigenvalue
  unsigned int                        eigind1 = 0;
  if( ImageDimension == 3 )
    {
    eigind0 = 2; // biggest eigenvalue
    eigind1 = 1;
    }
  vnl_vector<RealType> evec1_2ndary = eig1.get_eigenvector( eigind1 );
  vnl_vector<RealType> evec1_primary = eig1.get_eigenvector( eigind0 );
  vnl_vector<RealType> evec2_2ndary = eig2.get_eigenvector( eigind1 );
  vnl_vector<RealType> evec2_primary = eig2.get_eigenvector( eigind0 );
  /** Solve Wahba's problem --- http://en.wikipedia.org/wiki/Wahba%27s_problem */
  // NOT USED RealType wt0 = fabs( eig1.get_eigenvalue( eigind0 ) ) + fabs( eig2.get_eigenvalue( eigind0 ) );
  // NOT USED RealType wt1 = fabs( eig1.get_eigenvalue( eigind1 ) ) + fabs( eig2.get_eigenvalue( eigind1 ) );
  // NOT USED RealType evsum = wt0 + wt1;
  vnl_matrix<RealType> B = outer_product( evec2_primary, evec1_primary );
  if( ImageDimension == 3 )
    {
    //    B = outer_product( evec2_2ndary , evec1_2ndary ) * wt1 / evsum +
    //        outer_product( evec2_primary , evec1_primary ) * wt0 / evsum;
    B = outer_product( evec2_2ndary, evec1_2ndary )
      + outer_product( evec2_primary, evec1_primary );
    }
  vnl_svd<RealType>    wahba( B );
  vnl_matrix<RealType> A_solution = wahba.V() * wahba.U().transpose();
  //   " << vnl_determinant<RealType>(  wahba.U()) << std::endl;
  // now rotate the points to the same frame and sample the neighborhoods again
  for( unsigned int ii = 0; ii < Gsz; ii++ )
    {
    PointType            ptran = imagepatch2[ii];
    vnl_vector<RealType> ptran2( ptran.Size(), 0 );
    for( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      ptran[dd] -= avgpoint2[dd]; ptran2[dd] = ptran[dd];
      }
    // rotate ptran
    ptran2 = ( A_solution ) * ptran2;
    for( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      ptran[dd] = ptran2[dd] + static_cast<RealType>( avgpoint2[dd] );
      }
    if( interp2->IsInsideBuffer( ptran ) )
      {
      sample2[ii] = interp2->Evaluate( ptran );
      }
    else
      {
      ok = false;
      }
    }
  if( ok )
    {
    mean2 = sample2.mean();
    sample2 = ( sample2 - mean2 );
    sd2 = sqrt( sample2.squared_magnitude() );
    correlation = inner_product( sample1, sample2 ) / ( sd1 * sd2 );
    } // done applying wahba solution
  else
    {
    correlation = 0;
    }

  if( std::isnan( correlation ) || std::isinf( correlation )  )
    {
    return 0;
    }
  else
    {
    return correlation;
    }
}

template <typename TRealType>
void Sinkhorn( vnl_matrix<TRealType>&  correspondencematrix  )
{
  // std::cout << " SH begin " << std::endl;

  typedef TRealType RealType;
  // now that we have the correspondence matrix we convert it to a "probability" matrix
  for( unsigned int loop = 0; loop < 2; loop++ )
    {
    for( unsigned int ii = 0; ii < correspondencematrix.cols(); ii++ )
      {
      vnl_vector<RealType> correspondencematrixcol = correspondencematrix.get_column( ii );
      RealType             csum = correspondencematrixcol.sum();
      if( csum > 0 )
        {
        correspondencematrix.set_column( ii, correspondencematrixcol / csum );
        }
      }
    for( unsigned int ii = 0; ii < correspondencematrix.rows(); ii++ )
      {
      vnl_vector<RealType> correspondencematrixrow = correspondencematrix.get_row( ii );
      RealType             rsum = correspondencematrixrow.sum();
      if( rsum > 0 )
        {
        correspondencematrix.set_row( ii, correspondencematrixrow / rsum );
        }
      }
    }
  // std::cout << " SH done " << std::endl;
}

template <typename TImage>
typename TImage::Pointer ComputeLaplacianImage( typename TImage::Pointer input )
{
  typedef itk::LaplacianRecursiveGaussianImageFilter<TImage, TImage> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetSigma( 1.0 );
  filter->SetInput( input );
  filter->Update();
  return filter->GetOutput();
}

template <unsigned int ImageDimension>
int PureTissueN4WeightMask( int argc, char *argv[] )
{
  // This function is used to create the weight mask for N4.
  // The basic idea is that we want a weight mask which emphasizes
  // those voxels which are comprised of "pure" tissue types, e.g.
  // 100% prob. that it is gray matter.  In words, suppose I want a
  // pure tissue weight mask comprised of gray matter and white matter.
  // To do this, I calculate the probability that something is gray
  // matter and not white matter or that something is white matter and
  // not gray matter.  Generalized, this formula is
  //
  //   weightMask = \sum_{i=1}^N P_i ( \Prod_{j=1,j!=i}^N (1 - P_j ) )
  //
  // where P_i is the i^th probability

  if( argc < 5 )
    {
    // std::cout << "Usage: " << argv[0] << " ImageDimension";
    // std::cout << " outputWeightImage PureTissueN4WeightMask";
    // std::cout << " probImage1 probImage2 ... probImageN" << std::endl;
    return EXIT_FAILURE;
    }

  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  std::vector<typename ImageType::Pointer> images;
  for( int n = 4; n < argc; n++ )
    {
    typename ImageType::Pointer image;
    ReadImage<ImageType>( image, argv[n] );
    images.push_back( image );
    }

  typename ImageType::Pointer output = ImageType::New();
  output->CopyInformation( images[0] );
  output->SetRegions( images[0]->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<ImageType> ItO( output, output->GetLargestPossibleRegion() );
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    PixelType probability = 0.0f;
    for( unsigned int i = 0; i < images.size(); i++ )
      {
      PixelType negation = 1.0f;
      for( unsigned int j = 0; j < images.size(); j++ )
        {
        if( i == j )
          {
          continue;
          }
        negation *= ( 1.0f - images[j]->GetPixel( ItO.GetIndex() ) );
        }
      probability += negation * images[i]->GetPixel( ItO.GetIndex() );
      }
    ItO.Set( probability );
    }

  WriteImage<ImageType>( output, argv[2] );

  return EXIT_SUCCESS;
}


template <unsigned int ImageDimension>
int PMSmoothImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;

  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sigma = 1.0;
  if( argc > argct )
    {
    sigma = atof(argv[argct]); argct++;
    }
  float       conductance = 0.25;
  if( argc > argct )
    {
    conductance = atof(argv[argct]); argct++;
    }
  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer varimage = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  double     spacingsize = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    double sp = static_cast<double>( image1->GetSpacing()[d] );
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );
  typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType,
							ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image1 );
  filter->SetNumberOfIterations( sigma );
  double mytimestep = spacingsize / std::pow( 2.0, static_cast<double>(ImageDimension+1) );
  double reftimestep = 0.4 / std::pow( 2.0, static_cast<double>(ImageDimension+1) );
  if ( mytimestep > reftimestep ) mytimestep = reftimestep;
  filter->SetTimeStep( mytimestep );
  filter->SetConductanceParameter( conductance ); // might need to change this
  filter->Update();
  varimage = filter->GetOutput();
  WriteImage<ImageType>( varimage, outname.c_str() );
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int ConvolveImage( int argc, char *argv[] )
{
  typedef float                                  PixelType;
  typedef itk::Image<PixelType, ImageDimension>  ImageType;

  int argct = 2;
  const std::string outputFileName = std::string( argv[argct] );
  argct += 2;

  std::string inputFileName = std::string( argv[argct] );   argct++;
  std::string kernelFileName = std::string( argv[argct] );  argct++;

  bool normalize = true;
  if( argc > argct )
    {
    normalize = static_cast<bool>( std::stoi( argv[argct] ) );
    }

  typename ImageType::Pointer inputImage = nullptr;
  typename ImageType::Pointer kernelImage = nullptr;

  ReadImage<ImageType>( inputImage, inputFileName.c_str() );
  ReadImage<ImageType>( kernelImage, kernelFileName.c_str() );

  typedef itk::ConvolutionImageFilter<ImageType> FilterType;
  typename FilterType::Pointer convoluter = FilterType::New();
  convoluter->SetInput( inputImage );
  convoluter->SetKernelImage( kernelImage );
  convoluter->SetNormalize( normalize );
  convoluter->Update();

  typename ImageType::Pointer outputImage = convoluter->GetOutput();
  outputImage->Update();
  outputImage->DisconnectPipeline();

  WriteImage<ImageType>( outputImage, outputFileName.c_str() );
  return EXIT_SUCCESS;
}


template <unsigned int ImageDimension>
int InPaint(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sigma = 1.0;
  if( argc > argct )
    {
    sigma = atof(argv[argct]); argct++;
    }
  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  PixelType     spacingsize = 0;
  PixelType     minsp = image1->GetSpacing()[0];
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    PixelType sp = image1->GetSpacing()[d];
    if ( sp < minsp ) minsp = sp;
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );
  float       inpow = spacingsize * 2.0f;
  if( argc > argct )
    {
    inpow = atof(argv[argct]); argct++;
    }
  // # 1 job - create a kernel
  typename ImageType::Pointer kernel = ImageType::New();
  typename ImageType::IndexType start;
  typename ImageType::SizeType size;
  typename ImageType::RegionType region;
  start.Fill(0);
  size.Fill(3);
  region.SetSize(size);
  region.SetIndex(start);
  kernel->SetRegions(region);
  kernel->Allocate();
  kernel->SetSpacing( image1->GetSpacing() );
  unsigned long kernelsize = region.GetNumberOfPixels();
  itk::ImageRegionIterator<ImageType> imageIterator(kernel, region);
  unsigned int ct = 0;
  typename ImageType::PointType centerPoint;
  centerPoint.Fill( 0 );
  typename ImageType::PointType locPoint;
  locPoint.Fill( 0 );
  while(!imageIterator.IsAtEnd())
    {
    if ( ct == static_cast<unsigned int>( std::floor(static_cast<float>( kernelsize ) / 2.0f ) ) )
      {
      kernel->TransformIndexToPhysicalPoint(  imageIterator.GetIndex(), centerPoint );
      }
    ++ct;
    ++imageIterator;
    }
  unsigned int ct2 = 0;
  PixelType totalval = 0;
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    kernel->TransformIndexToPhysicalPoint(  imageIterator.GetIndex(), locPoint );
    PixelType val = 0;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      PixelType delt = ( locPoint[d] - centerPoint[d] );
      val += delt * delt;
      }
    PixelType distval = sqrt( val );
    if ( distval > inpow ) val = 1.e8;
    if ( val > 0 )
      {
      imageIterator.Set( 1.0f / val );
      totalval += imageIterator.Get( );
      }

    if ( ct2 == static_cast<unsigned int>( std::floor( static_cast<float>( ct ) / 2.0f ) ) ) imageIterator.Set( 1.e-8 );
    ++ct2;
    ++imageIterator;
    }
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set( imageIterator.Get() / totalval );
    ++imageIterator;
    }
  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( image1 );
  duplicator->Update();
  typename ImageType::Pointer varimage = duplicator->GetOutput();
  for ( unsigned int i = 0; i < sigma; i++ )
    {
    typedef itk::ConvolutionImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( varimage );
    filter->SetKernelImage(kernel);
    filter->Update();
    varimage = filter->GetOutput();
    Iterator vfIter( varimage,  varimage->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      PixelType pixval = image1->GetPixel( vfIter.GetIndex() );
      if ( pixval > 0 ) vfIter.Set( pixval );
      }
    }
  WriteImage<ImageType>( varimage, outname.c_str() );
  return EXIT_SUCCESS;
}



template <unsigned int ImageDimension>
int InPaint2(int argc, char *argv[])
{
  //  I_t = \nabla ( \Delta I ) \cdot N   where N is \perpto the normal direction
  typedef float                                                           PixelType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sigma = 1.0;
  if( argc > argct )
    {
    sigma = atof(argv[argct]); argct++;
    }
  typename ImageType::Pointer image1 = nullptr;
  typename ImageType::Pointer varimage = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str() );
  PixelType     spacingsize = 0;
  PixelType     minsp = image1->GetSpacing()[0];
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    PixelType sp = image1->GetSpacing()[d];
    if ( sp < minsp ) minsp = sp;
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );
  float       inpow = spacingsize * 2.0f;
  if( argc > argct )
    {
    inpow = atof(argv[argct]); argct++;
    }
  // # 1 job - create a kernel
  typename ImageType::Pointer kernel = ImageType::New();
  typename ImageType::IndexType start;
  typename ImageType::SizeType size;
  typename ImageType::RegionType region;
  start.Fill(0);
  size.Fill(3);
  region.SetSize(size);
  region.SetIndex(start);
  kernel->SetRegions(region);
  kernel->Allocate();
  kernel->SetSpacing( image1->GetSpacing() );
  unsigned long kernelsize = region.GetNumberOfPixels();
  itk::ImageRegionIterator<ImageType> imageIterator(kernel, region);
  unsigned int ct = 0;
  typename ImageType::PointType centerPoint;
  typename ImageType::PointType locPoint;
  while(!imageIterator.IsAtEnd())
    {
    if ( ct == static_cast<unsigned int>( std::floor( static_cast<float>( kernelsize ) / 2.0f ) ) )
      {
      kernel->TransformIndexToPhysicalPoint(  imageIterator.GetIndex(), centerPoint );
      }
    ++ct;
    ++imageIterator;
    }
  unsigned int ct2 = 0;
  PixelType totalval = 0;
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    kernel->TransformIndexToPhysicalPoint(  imageIterator.GetIndex(), locPoint );
    PixelType val = 0;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      PixelType delt = ( locPoint[d] - centerPoint[d] );
      val += delt * delt;
      }
    PixelType distval = sqrt( val );
    if ( distval > inpow ) val = 1.e8;
    if ( val > 0 )
      {
      imageIterator.Set( 1.0f / val );
      totalval += imageIterator.Get( );
      }
    if ( ct2 == static_cast<unsigned int>( std::floor( static_cast<float>( ct ) / 2.0f ) ) ) imageIterator.Set( 1.e-8 );
    ++ct2;
    ++imageIterator;
    }
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set( imageIterator.Get() / totalval );
    ++imageIterator;
    }
  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( image1 );
  duplicator->Update();
  varimage =  duplicator->GetOutput();
  for ( unsigned int i = 0; i < sigma; i++ )
    {
    typedef itk::ConvolutionImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( varimage );
    filter->SetKernelImage(kernel);
    filter->Update();
    varimage = filter->GetOutput();
    Iterator vfIter( varimage,  varimage->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      PixelType pixval = image1->GetPixel( vfIter.GetIndex() );
      if ( pixval > 0 ) vfIter.Set( pixval );
      }
    }
  WriteImage<ImageType>( varimage, outname.c_str() );
  return EXIT_SUCCESS;
}



template <unsigned int ImageDimension>
int Check3TissueLabeling( int argc, char *argv[] )
{
  // This function is used for quality control in the abp.sh pipeline.
  // On the UVa cluster, the labels after the segmentation step on
  // random data sets would be of some odd permutation.  Under expected
  // conditions the labels should be  CSF -> 1, GM -> 2, WM -> 3 but for
  // some reason which we can't reproduce, they'd be some other ordering.
  // The warped priors are correct so we use those images and the segmentation
  // to reorder the labels and move the posteriors where appropriate.

  if( argc < 8 )
    {
    // std::cout << "Usage: " << argv[0] << " ImageDimension";
    // std::cout << " segmentationImage Check3TissueLabeling";
    // std::cout << " priorWarped1 priorWarped2 priorWarped3";
    // std::cout << " posteriorWarped1 posteriorWarped2 posteriorWarped3" << std::endl;
    return EXIT_FAILURE;
    }

  typedef double       PixelType;
  typedef unsigned int LabelType;

  constexpr unsigned int NumberOfLabels = 3;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typename LabelImageType::Pointer labelImage;
  ReadImage<LabelImageType>( labelImage, argv[2] );

  typename ImageType::Pointer priors[NumberOfLabels];
  typename ImageType::Pointer posteriors[NumberOfLabels];
  for( unsigned int d = 0; d < NumberOfLabels; d++ )
    {
    ReadImage<ImageType>( priors[d], argv[d + 4] );
    ReadImage<ImageType>( posteriors[d], argv[d + 7] );
    }

  typename LabelImageType::Pointer maxPriorLabelImage = LabelImageType::New();
  maxPriorLabelImage->CopyInformation( labelImage );
  maxPriorLabelImage->SetRegions( labelImage->GetRequestedRegion() );
  maxPriorLabelImage->Allocate();
  maxPriorLabelImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( labelImage, labelImage->GetRequestedRegion() );
  itk::ImageRegionIterator<LabelImageType>          ItM( maxPriorLabelImage, maxPriorLabelImage->GetRequestedRegion() );
  for( ItL.GoToBegin(), ItM.GoToBegin(); !ItL.IsAtEnd(); ++ItM, ++ItL )
    {
    if( ItL.Get() == 0 )
      {
      continue;
      }

    LabelType maxLabel = 1;
    PixelType maxPrior = priors[0]->GetPixel( ItL.GetIndex() );
    for( LabelType d = 2; d <= 3; d++ )
      {
      PixelType prior = priors[d - 1]->GetPixel( ItL.GetIndex() );
      if( prior > maxPrior )
        {
        maxPrior = prior;
        maxLabel = d;
        }
      }
    ItM.Set( maxLabel );
    }

  itk::Matrix<LabelType, 6, 3> permutations;

  unsigned int which = 0;

  permutations(which, 0) = 1;
  permutations(which, 1) = 2;
  permutations(which, 2) = 3;

  permutations(++which, 0) = 1;
  permutations(which, 1) = 3;
  permutations(which, 2) = 2;

  permutations(++which, 0) = 2;
  permutations(which, 1) = 1;
  permutations(which, 2) = 3;

  permutations(++which, 0) = 2;
  permutations(which, 1) = 3;
  permutations(which, 2) = 1;

  permutations(++which, 0) = 3;
  permutations(which, 1) = 1;
  permutations(which, 2) = 2;

  permutations(++which, 0) = 3;
  permutations(which, 1) = 2;
  permutations(which, 2) = 1;

  PixelType maxDice = 0.0;
  int       maxPermutationRow = -1;
  for( unsigned r = 0; r < 6; r++ )
    {
    typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( labelImage );
    duplicator->Update();

    typename LabelImageType::Pointer permutedLabelImage = duplicator->GetOutput();

    itk::ImageRegionIterator<LabelImageType> ItP( permutedLabelImage, permutedLabelImage->GetRequestedRegion() );
    for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
      {
      LabelType permutedLabel = ItP.Get();
      if( permutedLabel != 0 )
        {
        unsigned int whichColumn = permutedLabel - 1;
        ItP.Set( permutations( r, whichColumn ) );
        }
      }

    typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetSourceImage( permutedLabelImage );
    filter->SetTargetImage( maxPriorLabelImage );
    filter->Update();

    PixelType dice = filter->GetMeanOverlap();
    // std::cout << r << ": " << dice << std::endl;
    if( dice > maxDice )
      {
      maxPermutationRow = r;
      maxDice = dice;
      }
    }

  if( maxPermutationRow == -1 )
    {
    std::cerr << "Unexpected problem." << std::endl;
    return EXIT_FAILURE;
    }

  LabelType movingLabels[3];
  LabelType fixedLabels[3];
  for( unsigned int d = 0; d < 3; d++ )
    {
    // std::cout << d + 1 << " -> " << permutations( maxPermutationRow, d ) << std::endl;
    movingLabels[d] = permutations( maxPermutationRow, d );
    fixedLabels[d] = d + 1;
    }

  if( maxPermutationRow == 0 )
    {
    // std::cout << "No need to change labels/posteriors." << std::endl;
    }
  else
    {
    for( unsigned int d = 0; d < 3; d++ )
      {
      LabelType movingLabel = movingLabels[d];
      LabelType fixedLabel = fixedLabels[d];

      if( movingLabel == fixedLabel )
        {
        continue;
        }
      else
        {
        WriteImage<ImageType>( posteriors[movingLabels[d] - 1], argv[7 + d] );
        WriteImage<ImageType>( posteriors[movingLabels[movingLabels[d] - 1] - 1], argv[7 + movingLabels[d] - 1] );

        LabelType tmp = movingLabels[d];
        movingLabels[d] = movingLabels[tmp - 1];
        movingLabels[tmp - 1] = tmp;
        }
      if( d == 0 )
        {
        typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
        typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage( labelImage );
        duplicator->Update();

        typename LabelImageType::Pointer permutedLabelImage = duplicator->GetOutput();

        itk::ImageRegionIteratorWithIndex<LabelImageType> ItP( permutedLabelImage,
                                                               permutedLabelImage->GetRequestedRegion() );
        for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
          {
          LabelType permutedLabel = ItP.Get();
          if( permutedLabel != 0 )
            {
            unsigned int whichColumn = permutedLabel - 1;
            ItP.Set( permutations( maxPermutationRow, whichColumn ) );
            }
          }

        // std::cout << "Relabeling segmentation image." << std::endl;
        WriteImage<LabelImageType>( permutedLabelImage, argv[2] );
        }
      }
    }
  return EXIT_SUCCESS;
}

template <typename T>
struct blob_index_cmp
  {
  blob_index_cmp(const T arr) : barr(arr)
  {
  }

  bool operator()(const size_t a, const size_t b) const
  {
    return barr[a] < barr[b];
  }

  const T barr;
  };

template <unsigned int ImageDimension, typename TImage, typename BlobsListType>
void getBlobCorrespondenceMatrix( unsigned int radval, typename TImage::Pointer image,  typename TImage::Pointer image2,
                                  vnl_matrix<float>& correspondencematrix, BlobsListType blobs1,  BlobsListType blobs2,
                                  float gradsig, bool dosinkhorn )
{
  typedef float                                                                   PixelType;
  typedef float                                                                   RealType;
  typedef itk::Image<PixelType, ImageDimension>                                   ImageType;
  typedef itk::CovariantVector<RealType, ImageDimension>                          GradientPixelType;
  typedef itk::Image<GradientPixelType, ImageDimension>                           GradientImageType;
  typedef itk::SmartPointer<GradientImageType>                                    GradientImagePointer;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer                               GradientImageFilterPointer;
  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType>              BlobFilterType;
  typedef typename BlobFilterType::BlobPointer                                    BlobPointer;

  GradientImageFilterPointer gfilter = GradientImageFilterType::New();
  gfilter->SetInput( image );
  gfilter->SetSigma( gradsig );
  gfilter->Update();
  GradientImagePointer gimage = gfilter->GetOutput();

  GradientImageFilterPointer gfilter2 = GradientImageFilterType::New();
  gfilter2->SetInput( image2 );
  gfilter2->SetSigma( gradsig );
  gfilter2->Update();
  GradientImagePointer gimage2 = gfilter2->GetOutput();

  // now compute some feature characteristics in each blob
  correspondencematrix.set_size( blobs1.size(), blobs2.size() );
  correspondencematrix.fill( 0 );
  typedef typename ImageType::IndexType        IndexType;
  typedef itk::NeighborhoodIterator<ImageType> niteratorType;
  typename niteratorType::RadiusType rad;
  rad.Fill( radval );
  niteratorType GHood(  rad, image, image->GetLargestPossibleRegion() );
  niteratorType GHood2( rad, image2, image2->GetLargestPossibleRegion() );
  IndexType     zeroind;  zeroind.Fill( radval );
  GHood.SetLocation( zeroind );
  // get indices within a ND-sphere
  std::vector<unsigned int> activeindex;
  std::vector<RealType>     weights;
  RealType                  weightsum = 0;
  for( unsigned int ii = 0; ii < GHood.Size(); ii++ )
    {
    IndexType ind = GHood.GetIndex( ii );
    RealType  dist = 0;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      dist += ( ind[jj] - zeroind[jj] ) * ( ind[jj] - zeroind[jj] );
      }
    dist = sqrt( dist );
    if( dist <= radval )
      {
      activeindex.push_back( ii );
      RealType wt =  exp( -1.0f * dist / ( radval * radval ) );
      weights.push_back( wt );
      weightsum += ( wt );
      }
    }
  for(float & weight : weights)
    {
    weight = weight / weightsum;
    }
  BlobPointer bestblob = nullptr;
  if( ( !blobs2.empty() ) && ( !blobs1.empty() ) )
    {
    unsigned int count2;
    RealType     smallval = 1.e-4;
    typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
    typedef typename ScalarInterpolatorType::Pointer              InterpPointer;
    InterpPointer interp1 =  ScalarInterpolatorType::New();
    interp1->SetInputImage(image);
    InterpPointer interp2 =  ScalarInterpolatorType::New();
    interp2->SetInputImage(image2);
    unsigned int count1 = 0;
    for( unsigned int i = 0; i < blobs1.size(); i++ )
      {
      BlobPointer blob1 = blobs1[i];
      IndexType   indexi = blob1->GetCenter();
      if( image->GetPixel( indexi ) > smallval )
        {
        GHood.SetLocation( indexi );
        bestblob = nullptr;
        count2 = 0;
        for( unsigned int j = 0; j < blobs2.size(); j++ )
          {
          const BlobPointer & blob2 = blobs2[j];
          const IndexType &   indexj = blob2->GetCenter();
          if( image2->GetPixel( indexj ) > smallval )
            {
            GHood2.SetLocation( indexj );
            RealType correlation =
              PatchCorrelation<ImageDimension, RealType, ImageType, GradientImageType, InterpPointer>( GHood, GHood2,
                                                                                                       activeindex,
                                                                                                       weights, gimage,
                                                                                                       gimage2,
                                                                                                       interp2 );
            if( correlation < 0 )
              {
              correlation = 0;
              }
            correspondencematrix( i, j ) = correlation;
            count2++;
            }
          }
        count1++;
        } // imagei GetPixel gt 0
      if( i % 100 == 0 )
        {
        // std::cout << " Progress : " << (float ) i / (float) blobs1.size() * 100.0 << std::endl;
        }
      }
    if( dosinkhorn )
      {
      Sinkhorn<RealType>( correspondencematrix );
      }
    }
  return;
}

template <unsigned int ImageDimension>
int BlobDetector( int argc, char *argv[] )
{
  typedef float                                                                   PixelType;
  typedef float                                                                   RealType;
  typedef itk::Image<PixelType, ImageDimension>                                   ImageType;

  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }
  // sensitive parameters are set here - begin
  RealType     gradsig = 1.0;      // sigma for gradient filter
  unsigned int stepsperoctave = 10; // number of steps between doubling of scale
  RealType     minscale = std::pow( 1.0, 1.0 );
  RealType     maxscale = std::pow( 2.0, 10.0 );
  RealType     uniqfeat_thresh = 0.01;
  RealType     smallval = 1.e-2; // assumes images are normalizes in [ 0, 1 ]
  bool         dosinkhorn = false;
  RealType     maxradiusdiffallowed = 0.25; // IMPORTANT feature size difference
  RealType     kneighborhoodval = 3;        // IMPORTANT - defines how many nhood nodes to use in k-hood definition
  unsigned int radval = 20;                 // IMPORTANT radius for correlation
  RealType     dthresh = 0.02;              // IMPORTANT distance preservation threshold
  // sensitive parameters are set here - end
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  std::string       outname2 = std::string("temp.nii.gz");
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int nblobs = std::stoi( argv[argct] );   argct++;
  std::string  fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]); argct++;
    }
  if( argc > argct )
    {
    outname2 = std::string(argv[argct]); argct++;
    }
  RealType corrthresh = 0;
  if( argc > argct )
    {
    corrthresh = atof(argv[argct]); argct++;
    }
  if( argc > argct )
    {
    radval = std::stoi(argv[argct]); argct++;
    }
  if( argc > argct )
    {
    dthresh = atof(argv[argct]); argct++;
    }
  typename ImageType::Pointer image;
  typename ImageType::Pointer image2;
  ReadImage<ImageType>( image, fn1.c_str() );
  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typename BlobFilterType::Pointer blobFilter = BlobFilterType::New();
  typedef typename BlobFilterType::BlobPointer BlobPointer;
  blobFilter->SetStartT( minscale );
  blobFilter->SetEndT( maxscale );
  blobFilter->SetStepsPerOctave( stepsperoctave );
  blobFilter->SetNumberOfBlobs( nblobs );
  blobFilter->SetInput( image ); /*ComputeLaplacianImage<ImageType>( image ) ); */
  blobFilter->Update();
  typedef typename BlobFilterType::BlobRadiusImageType BlobRadiusImageType;
  typename BlobRadiusImageType::Pointer labimg = blobFilter->GetBlobRadiusImage();
  typename BlobRadiusImageType::Pointer labimg2;
  WriteImage<BlobRadiusImageType>( labimg, outname.c_str() );
  typedef typename BlobFilterType::BlobsListType BlobsListType;
  BlobsListType blobs1 =  blobFilter->GetBlobs();

  vnl_matrix<RealType> correspondencematrix1;
  correspondencematrix1.set_size( blobs1.size(), blobs1.size() );
  correspondencematrix1.fill( 1 );
  vnl_matrix<RealType> correspondencematrix2;
  vnl_matrix<RealType> correspondencematrix;
  // getBlobCorrespondenceMatrix<ImageDimension,ImageType,BlobsListType >( radval, image, image, correspondencematrix1,
  // blobs1, blobs1, gradsig , dosinkhorn );

  BlobsListType blobs2;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>( image2, fn2.c_str() );
    typename BlobFilterType::Pointer blobFilter2 = BlobFilterType::New();
    blobFilter2->SetStartT( minscale );
    blobFilter2->SetEndT( maxscale );
    blobFilter2->SetStepsPerOctave( stepsperoctave );
    blobFilter2->SetNumberOfBlobs( nblobs );
    blobFilter2->SetInput( image2 );
    blobFilter2->Update();
    labimg2 = blobFilter2->GetBlobRadiusImage();
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    labimg->FillBuffer( 0 );
    labimg2->FillBuffer( 0 );
    blobs2 =  blobFilter2->GetBlobs();
    correspondencematrix2.set_size( blobs2.size(), blobs2.size() );
    correspondencematrix2.fill( 1 );
    //    getBlobCorrespondenceMatrix<ImageDimension,ImageType,BlobsListType > ( radval, image2, image2,
    // correspondencematrix2, blobs2, blobs2, gradsig, dosinkhorn );
    }
  else
    {
    return EXIT_SUCCESS;
    }
  // std::cout << " Blob1Length " << blobs1.size() << " Blob2Length " << blobs2.size() << std::endl;
  // now compute some feature characteristics in each blob
  typedef typename ImageType::IndexType IndexType;
  IndexType   zeroind;  zeroind.Fill( radval );
  BlobPointer bestblob = nullptr;
  if( ( !blobs2.empty() ) && ( !blobs1.empty() ) )
    {
    getBlobCorrespondenceMatrix<ImageDimension, ImageType, BlobsListType>
      ( radval, image, image2, correspondencematrix, blobs1, blobs2, gradsig, dosinkhorn );
    // vnl_matrix<RealType> diagcorr = outer_product( correspondencematrix1.get_diagonal(),
    //  correspondencematrix2.get_diagonal() );
    // Sinkhorn<RealType>( diagcorr );
    // for ( unsigned int row = 0; row < correspondencematrix.rows(); row++ )
    //  for ( unsigned int col = 0; col < correspondencematrix.cols(); col++ )
    //	correspondencematrix( row, col ) *= diagcorr( row, col );
    //    Sinkhorn<RealType>( correspondencematrix );
    unsigned int matchpt = 1;
    std::cout << " now compute pairwise matching " << correspondencematrix.max_value() << " reducing to "
              << corrthresh << std::endl;
    unsigned int count1 = 0;
    typedef std::pair<BlobPointer, BlobPointer> BlobPairType;
    std::vector<BlobPairType> blobpairs;
    while( ( matchpt < ( corrthresh + 1 ) ) && ( count1 <  blobs1.size() ) )
      {
      unsigned int maxpair = correspondencematrix.arg_max();
      unsigned int maxrow = ( unsigned int )  maxpair / correspondencematrix.cols();
      unsigned int maxcol = maxpair - maxrow * correspondencematrix.cols();
      BlobPointer  blob1 = blobs1[maxrow];
      bestblob = blobs2[maxcol];
      if( bestblob &&  bestblob->GetObjectRadius() > 1 )
        {
        if( static_cast<RealType>( std::fabs( bestblob->GetObjectRadius() - blob1->GetObjectRadius() ) ) < maxradiusdiffallowed )
          {
          if( bestblob && ( image->GetPixel( blob1->GetCenter() ) > smallval )  &&
              ( image2->GetPixel( bestblob->GetCenter() )  > smallval ) &&
              ( correspondencematrix1(maxrow, maxrow) > uniqfeat_thresh ) &&
              ( correspondencematrix2(maxcol, maxcol) > uniqfeat_thresh  )
              )
            {
            BlobPairType blobpairing = std::make_pair( blob1, bestblob );
            blobpairs.push_back( blobpairing );
            std::cout << " best correlation " << correspondencematrix.absolute_value_max() << " rad1 "
                      << blob1->GetObjectRadius() << " rad2 " << bestblob->GetObjectRadius() << " : " << matchpt
                      << std::endl;
            labimg->SetPixel(     blob1->GetCenter(), matchpt ); // ( int ) ( 0.5 +   ( *i )->GetObjectRadius() ) );
            labimg2->SetPixel( bestblob->GetCenter(), matchpt ); // ( int ) ( 0.5 + bestblob->GetObjectRadius() ) );
            matchpt++;
            }
          }
        }
      correspondencematrix.set_row( maxrow, correspondencematrix.get_row( 0 ).fill( 0 ) );
      correspondencematrix.set_column( maxcol, correspondencematrix.get_column( 0 ).fill( 0 ) );
      count1++;
      }

    /** For every blob, compute the distance to its neighbors before and after matching */
    vnl_matrix<RealType> distmatpre( blobpairs.size(), blobpairs.size() );
    distmatpre.fill( 0 );
    vnl_matrix<RealType> distmatpost( blobpairs.size(), blobpairs.size() );
    distmatpost.fill( 0 );
    vnl_matrix<RealType> distratiomat( blobpairs.size(), blobpairs.size() );
    distratiomat.fill( 0 );
    if( true )
      {
      for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
        {
        IndexType             blobind = blobpairs[bp].first->GetCenter();
        IndexType             blobpairind  = blobpairs[bp].second->GetCenter();
        std::vector<RealType> distspre;
        std::vector<RealType> distspost;
        std::vector<size_t>   distspreind;
        std::vector<size_t>   distspostind;
        for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
          {
          IndexType blobneighborind = blobpairs[bp2].first->GetCenter();
          IndexType blobpairneighborind  = blobpairs[bp2].second->GetCenter();
          RealType  dist1 = 0;
          RealType  dist2 = 0;
          for( unsigned int dim = 0; dim < ImageDimension; dim++ )
            {
            RealType delta1 = blobind[dim] - blobneighborind[dim];
            RealType delta2 = blobpairind[dim]  - blobpairneighborind[dim];
            dist1 += delta1 * delta1;
            dist2 += delta2 * delta2;
            }
          RealType drat = 0;
          if( dist1 > 0 )
            {
            drat = dist2 / dist1;
            }
          distspre.push_back( dist1 );
          distspost.push_back( dist2 );
          distspreind.push_back(  bp2 );
          distspostind.push_back( bp2 );
          //        // std::cout << " blob " << bp << " vs " << bp2 << sqrt( dist1 ) << " v " << sqrt( dist2 ) <<
          // std::endl;
          distmatpre(   bp, bp2 ) = distmatpre(   bp2, bp ) = dist1;
          distmatpost(  bp, bp2 ) = distmatpost(  bp2, bp ) = dist2;
          distratiomat( bp, bp2 ) = distratiomat( bp2, bp ) = drat;
          }
        }
      // now we have the distance ratio matrix - let's find a cluster of nodes with values near 1
      // count the k neighborhood for each blobpair possibility
      for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
        {
        IndexType    blobind = blobpairs[bp].first->GetCenter();
        IndexType    blobpairind  = blobpairs[bp].second->GetCenter();
        unsigned int kct = 0;
        for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
          {
          if( ( bp2 != bp ) && ( itk::Math::abs ( distratiomat( bp2, bp ) - 1 ) <  dthresh ) )
            {
            kct++;
            }
          }
        if( kct < kneighborhoodval )
          {
          labimg->SetPixel(  blobind, 0 );     // ( int ) ( 0.5 +   ( *i )->GetObjectRadius() ) );
          labimg2->SetPixel( blobpairind, 0 ); // ( int ) ( 0.5 + bestblob->GetObjectRadius() ) );
          }
        else
          {
          // std::cout << " blob " << bp << " keep " <<  distratiomat.get_row( bp ) << std::endl;
          }
        }
      } // if false
    if(  ( correspondencematrix1.rows() > 0 ) && ( false ) )
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "corrmat1.csv" );
      writer->SetInput( &correspondencematrix1 );
      writer->SetInput( &distmatpre );
      writer->Write();
      }
    if(  ( correspondencematrix2.rows() > 0 ) && ( false ) )
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "corrmat2.csv" );
      writer->SetInput( &correspondencematrix2 );
      writer->SetInput( &distratiomat );
      writer->Write();
      }
    WriteImage<BlobRadiusImageType>( labimg, outname.c_str() );
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    // std::cout << " Matched " << matchpt << " blobs " << std::endl;
    }
  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int MatchBlobs( int argc, char *argv[] )
{
  typedef float                                                      PixelType;
  typedef float                                                      RealType;
  typedef itk::Image<PixelType, ImageDimension>                      ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>               IteratorType;
  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typedef typename BlobFilterType::BlobPointer                       BlobPointer;
  typedef typename BlobFilterType::BlobRadiusImageType               BlobRadiusImageType;
  typedef typename BlobFilterType::BlobsListType                     BlobsListType;
  typedef typename BlobFilterType::BlobType                          BlobType;
  typedef typename ImageType::IndexType                              IndexType;
  if( argc < 5 )
    {
    // std::cout << " Not enough inputs " << std::endl;
    return 1;
    }
  // sensitive parameters are set here - begin
  RealType     gradsig = 1.0;    // sigma for gradient filter
  RealType     smallval = 1.e-2; // assumes images are normalizes in [ 0, 1 ]
  bool         dosinkhorn = false;
  RealType     maxradiusdiffallowed = 0.25; // IMPORTANT feature size difference
  unsigned int radval = 5;                  // IMPORTANT radius for correlation
  RealType     dthresh = 0.02;              // IMPORTANT distance preservation threshold
  // sensitive parameters are set here - end
  int               argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string landmarks1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer image;
  typename ImageType::Pointer imageLM;
  typename ImageType::Pointer image2;
  ReadImage<ImageType>( image, fn1.c_str() );
  ReadImage<ImageType>( imageLM, landmarks1.c_str() );
  ReadImage<ImageType>( image2, fn2.c_str() );
  vnl_matrix<RealType> correspondencematrix;
  typedef itk::MinimumMaximumImageCalculator<ImageType> LabelCalculatorType;
  typename LabelCalculatorType::Pointer calc = LabelCalculatorType::New();
  calc->SetImage( imageLM );
  calc->ComputeMaximum();
  typename ImageType::Pointer labimg  = MakeNewImage<ImageType>( image, 0 );
  typename ImageType::Pointer labimg2 = MakeNewImage<ImageType>( image2, 0 );
  typename ImageType::Pointer confimg2 = MakeNewImage<ImageType>( image2, 0 );
  float        maximum  = calc->GetMaximum();
  unsigned int corrthresh =  ( (unsigned int) maximum ) * 100;
  if( argc > argct )
    {
    corrthresh = atof(argv[argct]); argct++;
    }
  if( argc > argct )
    {
    radval = atof(argv[argct]); argct++;
    }
  if( argc > argct )
    {
    dthresh = atof(argv[argct]); argct++;
    }
  RealType kneighborhoodval = maximum;                                      // IMPORTANT - defines how many nhood nodes
                                                                            // to use in k-hood definition
  IteratorType                                  cIter( imageLM, imageLM->GetLargestPossibleRegion() );
  BlobsListType                                 blobs1;
  itk::Point<double, ImageType::ImageDimension> zeroPoint;
  zeroPoint.Fill(0);
  // std::cout << " N-Landmarks " << maximum << std::endl;
  for(  cIter.GoToBegin(); !cIter.IsAtEnd(); ++cIter )
    {
    RealType val = imageLM->GetPixel( cIter.GetIndex() );
    if( val > 0 )
      {
      typename BlobType::PointType centerPoint;
      image->TransformIndexToPhysicalPoint(  cIter.GetIndex(), centerPoint );
      BlobPointer blob = BlobType::New();
      blob->SetSigmaInObjectSpace( 1 );
      blob->SetScaleSpaceValue( val );
      blob->SetCenter(  cIter.GetIndex() );
      const typename BlobType::VectorType centerVector = centerPoint - zeroPoint;
      blob->GetModifiableObjectToParentTransform()->SetOffset(centerVector);
      blob->Update();
      blobs1.push_back( blob );
      }
    }

  typedef itk::ImageRandomConstIteratorWithIndex<ImageType> randIterator;
  randIterator  mIter( image2, image2->GetLargestPossibleRegion() );
  unsigned long numpx = image2->GetBufferedRegion().GetNumberOfPixels();
  unsigned int  n_samples = ( unsigned int ) ( ( float ) numpx ) * 0.1;
  mIter.SetNumberOfSamples( n_samples );
  BlobsListType blobs2;
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    // transform the center index into offset vector
    typename BlobType::PointType centerPoint;
    image2->TransformIndexToPhysicalPoint(  mIter.GetIndex(), centerPoint );
    BlobPointer blob = BlobType::New();
    blob->SetSigmaInObjectSpace( 1 );
    blob->SetScaleSpaceValue( image2->GetPixel( mIter.GetIndex() ) );
    blob->SetCenter(  mIter.GetIndex() );
    const typename BlobType::VectorType centerVector = centerPoint - zeroPoint;
    blob->GetModifiableObjectToParentTransform()->SetOffset(centerVector);
    blob->Update();
    blobs2.push_back( blob );
    }

  getBlobCorrespondenceMatrix<ImageDimension, ImageType, BlobsListType>
    ( radval, image, image2, correspondencematrix, blobs1, blobs2, gradsig, dosinkhorn );

  unsigned int matchpt = 1;
  std::cout << " now compute pairwise matching " << correspondencematrix.max_value() << " reducing to " << corrthresh
            << std::endl;
  unsigned int count1 = 0;
  typedef std::pair<BlobPointer, BlobPointer> BlobPairType;
  std::vector<BlobPairType> blobpairs;
  vnl_matrix<RealType>      correspondencematrix_hard( correspondencematrix );
  vnl_matrix<RealType>      correspondencematrix_soft( correspondencematrix );
  while( ( matchpt < ( corrthresh + 1 ) )  )
    {
    unsigned int maxpair = correspondencematrix_hard.arg_max();
    if( maxpair < 1.e-9 )
      {
      correspondencematrix_hard.update(  correspondencematrix_soft );
      maxpair = correspondencematrix_hard.arg_max();
      }
    unsigned int maxrow = ( unsigned int )  maxpair / correspondencematrix.cols();
    unsigned int maxcol = maxpair - maxrow * correspondencematrix.cols();
    BlobPointer  blob1 = blobs1[maxrow];
    BlobPointer  bestblob( blobs2[maxcol] );
    if( bestblob &&  bestblob->GetObjectRadius() > 1 )
      {
      if( static_cast<RealType>( std::fabs( bestblob->GetObjectRadius() - blob1->GetObjectRadius() ) ) < maxradiusdiffallowed )
        {
        if( bestblob && ( image->GetPixel( blob1->GetCenter() ) > smallval )  &&
            ( image2->GetPixel( bestblob->GetCenter() )  > smallval ) )
          {
          bestblob->SetScaleSpaceValue( correspondencematrix( maxrow, maxcol ) );
          BlobPairType blobpairing = std::make_pair( blob1, bestblob );
          blobpairs.push_back( blobpairing );
          matchpt++;
          }
        }
      }
    correspondencematrix_hard.set_row( maxrow, correspondencematrix_hard.get_row( maxrow ).fill( 0 ) );
    correspondencematrix_hard.set_column( maxcol, correspondencematrix_hard.get_column( maxcol ).fill( 0 ) );
    correspondencematrix_soft( maxrow, maxcol ) = 0;
    count1++;
    }

  /** For every blob, compute the distance to its neighbors before and after matching */
  vnl_matrix<RealType> distmatpre( blobpairs.size(), blobpairs.size() );
  distmatpre.fill( 0 );
  vnl_matrix<RealType> distmatpost( blobpairs.size(), blobpairs.size() );
  distmatpost.fill( 0 );
  vnl_matrix<RealType> distratiomat( blobpairs.size(), blobpairs.size() );
  distratiomat.fill( 0 );
  if( true )
    {
    for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
      {
      IndexType             blobind = blobpairs[bp].first->GetCenter();
      IndexType             blobpairind  = blobpairs[bp].second->GetCenter();
      std::vector<RealType> distspre;
      std::vector<RealType> distspost;
      std::vector<size_t>   distspreind;
      std::vector<size_t>   distspostind;
      for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
        {
        IndexType blobneighborind = blobpairs[bp2].first->GetCenter();
        IndexType blobpairneighborind  = blobpairs[bp2].second->GetCenter();
        RealType  dist1 = 0;
        RealType  dist2 = 0;
        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
          {
          RealType delta1 = blobind[dim] - blobneighborind[dim];
          RealType delta2 = blobpairind[dim]  - blobpairneighborind[dim];
          dist1 += delta1 * delta1;
          dist2 += delta2 * delta2;
          }
        RealType drat = 0;
        if( dist1 > 0 )
          {
          drat = dist2 / dist1;
          }
        distspre.push_back( dist1 );
        distspost.push_back( dist2 );
        distspreind.push_back(  bp2 );
        distspostind.push_back( bp2 );
        distmatpre(   bp, bp2 ) = distmatpre(   bp2, bp ) = dist1;
        distmatpost(  bp, bp2 ) = distmatpost(  bp2, bp ) = dist2;
        distratiomat( bp, bp2 ) = distratiomat( bp2, bp ) = drat;
        }
      }
    // now we have the distance ratio matrix - let's find a cluster of nodes with values near 1
    // count the k neighborhood for each blobpair possibility
    for( unsigned int bp = 0; bp < blobpairs.size(); bp++ )
      {
      IndexType    blobind = blobpairs[bp].first->GetCenter();
      IndexType    blobpairind  = blobpairs[bp].second->GetCenter();
      unsigned int kct = 0;
      typedef vnl_vector<RealType> kVectorType;
      kVectorType kLog1( kneighborhoodval, 0 );
      for( unsigned int bp2 = 0; bp2 < blobpairs.size(); bp2++ )
        {
        if( ( bp2 != bp ) && ( itk::Math::abs ( distratiomat( bp2, bp ) - 1 ) <  dthresh )
            //	  &&   ( blobpairs[bp2].first->GetScaleSpaceValue() != blobpairs[bp].first->GetScaleSpaceValue() )
            //   &&   ( blobpairs[bp2].second->GetScaleSpaceValue() != blobpairs[bp].second->GetScaleSpaceValue() )
            )
          {
          kct++;
          kLog1(  blobpairs[bp2].first->GetScaleSpaceValue() - 1 ) = 1;
          }
        }
      // if ( ( kLog1.sum() >= ( kneighborhoodval / 2 ) ) && ( kct >= kneighborhoodval ) )
      if( ( kct >= kneighborhoodval ) )
        {
        labimg->SetPixel(  blobind, blobpairs[bp].first->GetScaleSpaceValue() ); // ( int ) ( 0.5 +   ( *i
                                                                                 // )->GetObjectRadius() ) );
        labimg2->SetPixel( blobpairind, blobpairs[bp].first->GetScaleSpaceValue() );
        confimg2->SetPixel( blobpairind, blobpairs[bp].second->GetScaleSpaceValue() );
        //	// std::cout << " blob " << bp << " keep " <<  distratiomat.get_row( bp ) << " LM " <<
        // blobpairs[bp].first->GetScaleSpaceValue() <<  std::endl;
        }
      }
    } // if false
  if( true )
    {
    typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( "temp_corrmat.csv" );
    writer->SetInput( &correspondencematrix );
    writer->Write();
    }
  if( true )
    {
    typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( "temp_distmat2.csv" );
    writer->SetInput( &distratiomat );
    writer->Write();
    }
  std::string outname1 = outname + std::string("lm1.nii.gz");
  WriteImage<BlobRadiusImageType>( labimg, outname1.c_str() );
  std::string outname2 = outname + std::string("lm2.nii.gz");
  WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
  std::string outname3 = outname + std::string("conf2.nii.gz");
  WriteImage<BlobRadiusImageType>( confimg2, outname3.c_str() );
  return EXIT_SUCCESS;
}

//
// ImageMath was a gigantic switch statement that had 3 duplicated
// lists of 'if (operation == <op>)' clauses for 2d, 3d, and 4d. I
// figured out which functions were 2D only, 3D only and 4D Only,
// which were valid for all dimensions,  which were 2d and 3d, and
// which were 3d and 4d.
// So there's a template method for each case, and they're assembled
// for each dimension in an Explicit Template Function below.
template<unsigned DIM>
int
ImageMathHelper2DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "TileImages" )
    {
    TileImages<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesRegionCorr" )
    {
    TimeSeriesRegionCorr<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesRegionSCCA" )
    {
    TimeSeriesRegionSCCA<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

template <unsigned DIM>
int
ImageMathHelper2DOr3D(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "AverageLabels" )
    {
    AverageLabels<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Check3TissueLabeling" )
    {
    Check3TissueLabeling<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "oldPropagateLabelsThroughMask" )
    {
    PropagateLabelsThroughMask<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SetOrGetPixel" )
    {
    SetOrGetPixel<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "STAPLE" )
    {
    STAPLE<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

template <unsigned DIM>
int
ImageMathHelper3DOr4D(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "ConvertLandmarkFile" )
    {
    ConvertLandmarkFile<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PASLQuantifyCBF" )
    {
    PASLQuantifyCBF<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PASL" )
    {
    PASL<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "pCASL" )
    {
    pCASL<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TriPlanarView" )
    {
    TriPlanarView<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

template <unsigned DIM>
int
ImageMathHelper3DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "4DTensorTo3DTensor" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ComponentTo3DTensor" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "FuseNImagesIntoNDVectorField" )
    {
    FuseNImagesIntoNDVectorField<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ExtractComponentFrom3DTensor" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MTR" )
    {
    MTR<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SmoothTensorImage" )
    {
    SmoothTensorImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorAxialDiffusion" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorColor" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorEigenvalue" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorFA" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorFANumerator" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorFADenominator" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorIOTest" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorMask" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorMeanDiffusion" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorRadialDiffusion" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorToLocalSpace" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorToPhysicalSpace" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorToVectorComponent" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TensorToVector" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ValidTensor" )
    {
    TensorFunctions<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

template <unsigned DIM>
int
ImageMathHelper4DOnly(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "AverageOverDimension" )
    {
    AverageOverDimension<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CompCorrAuto" )
    {
    CompCorrAuto<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ComputeTimeSeriesLeverage" )
    {
    ComputeTimeSeriesLeverage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "nvols" )
    {
    PrintHeader<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PCASLQuantifyCBF" )
    {
    // PCASLQuantifyCBF<4>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SliceTimingCorrection" )
    {
    SliceTimingCorrection<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SplitAlternatingTimeSeries" )
    {
    SplitAlternatingTimeSeries<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ThreeTissueConfounds" )
    {
    ThreeTissueConfounds<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesMask" )
    {
    TimeSeriesMask<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesAssemble" )
    {
    TimeSeriesAssemble<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesDisassemble" )
    {
    TimeSeriesDisassemble<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesInterpolationSubtraction" )
    {
    TimeSeriesInterpolationSubtraction<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesSimpleSubtraction" )
    {
    TimeSeriesSimpleSubtraction<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesSubset" )
    {
    TimeSeriesSubset<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TimeSeriesToMatrix" )
    {
    TimeSeriesToMatrix<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

template <unsigned DIM>
int
ImageMathHelperAll(int argc, char **argv)
{
  std::string operation = std::string(argv[3]);
  if( operation == "m")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "mresample")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "+")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "-")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "vm")
    {
    VImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "vmresample")
    {
    VImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "v+")
    {
    VImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "v-")
    {
    VImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "/")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "^")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "exp")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "max")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "abs")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "addtozero")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "overadd")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "total")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "vtotal")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "mean")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Decision")
    {
    ImageMath<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Neg")
    {
    NegativeImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "G")
    {
    SmoothImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Convolve")
    {
    ConvolveImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PeronaMalik")
    {
    PMSmoothImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "InPaint")
    {
    InPaint<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MD" || operation == "ME" )
    {
    MorphImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MO" || operation == "MC" )
    {
    MorphImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "GD" || operation == "GE" )
    {
    MorphImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "GO" || operation == "GC" )
    {
    MorphImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "D" )
    {
    DistanceMap<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MaurerDistance"  )
    {
    GenerateMaurerDistanceImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Normalize" )
    {
    NormalizeImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Grad" )
    {
    GradientImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Laplacian" )
    {
    LaplacianImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Canny" )
    {
      CannyImage<DIM>(argc, argv);
      return EXIT_SUCCESS;
    }
  if( operation == "LabelSurfaceArea" )
    {
    LabelSurfaceArea<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PH" )
    {
    PrintHeader<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CenterImage2inImage1" )
    {
    CenterImage2inImage1<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Byte" )
    {
    ByteImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ReflectionMatrix" )
    {
    ReflectionMatrix<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MakeAffineTransform" )
    {
    MakeAffineTransform<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ClosestSimplifiedHeaderMatrix" )
    {
    ClosestSimplifiedHeaderMatrix<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "LabelStats" )
    {
    LabelStats<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ROIStatistics" )
    {
    ROIStatistics<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "LabelThickness" )
    {
    LabelThickness<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "LabelThickness2" )
    {
    LabelThickness2<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "DiceAndMinDistSum" )
    {
    DiceAndMinDistSum<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Lipschitz" )
    {
    Lipschitz<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "InvId" )
    {
    InvId<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ShiftImageSlicesInTime" )
    {
    ShiftImageSlicesInTime<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ReplicateImage" )
    {
    ReplicateImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ReplicateDisplacement" )
    {
    ReplicateDisplacement<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "GetLargestComponent" )
    {
    GetLargestComponent<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ExtractVectorComponent" )
    {
    ExtractVectorComponent<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ThresholdAtMean" )
    {
    ThresholdAtMean<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SetTimeSpacing" )
    {
    SetTimeSpacing<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SetTimeSpacingWarp" )
    {
    SetTimeSpacingWarp<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "FlattenImage" )
    {
    FlattenImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CorruptImage" )
    {
    CorruptImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Where" )
    {
    Where<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if(  operation == "Finite" )
    {
    Finite<DIM>( argc, argv );
    return EXIT_SUCCESS;
    }
  if( operation == "FillHoles" )
    {
    FillHoles<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "HistogramMatch" )
    {
    HistogramMatching<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "RescaleImage" )
    {
    RescaleImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "WindowImage" )
    {
    WindowImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "WindowImage" )
    {
    WindowImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "NeighborhoodStats" )
    {
    NeighborhoodStats<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PadImage" )
    {
    PadImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "SigmoidImage" )
    {
    SigmoidImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CoordinateComponentImages" )
    {
    CoordinateComponentImages<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Sharpen" )
    {
    SharpenImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MakeImage" )
    {
    MakeImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "stack" )
    {
    StackImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "stack2" )
    {
    Stack2Images<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CompareHeadersAndImages" )
    {
    CompareHeadersAndImages<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CountVoxelDifference" )
    {
    CountVoxelDifference<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "RemoveLabelInterfaces" )
    {
    RemoveLabelInterfaces<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ReplaceVoxelValue" )
    {
    ReplaceVoxelValue<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PoissonDiffusion" )
    {
    PoissonDiffusion<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "EnumerateLabelInterfaces" )
    {
    EnumerateLabelInterfaces<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ConvertImageToFile" )
    {
    ConvertImageToFile<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PValueImage" )
    {
    PValueImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CorrelationUpdate" )
    {
    CorrelationUpdate<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ConvertImageSetToMatrix" )
    {
    ConvertImageSetToMatrix<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "RandomlySampleImageSetToCSV" )
    {
    RandomlySampleImageSetToCSV<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ConvertImageSetToEigenvectors" )
    {
    ConvertImageSetToEigenvectors<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ConvertVectorToImage" )
    {
    ConvertVectorToImage<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PropagateLabelsThroughMask" )
    {
    itkPropagateLabelsThroughMask<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "FastMarchingExtension" )
    {
    FastMarchingExtension<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "FastMarchingSegmentation" )
    {
    FastMarchingSegmentation<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "TruncateImageIntensity" )
    {
    TruncateImageIntensity<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ExtractSlice" )
    {
    ExtractSlice<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "ClusterThresholdVariate" )
    {
    ClusterThresholdVariate<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MajorityVoting" )
    {
    MajorityVoting<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MostLikely" )
    {
    MostLikely<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "CorrelationVoting" )
    {
    CorrelationVoting<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PearsonCorrelation" )
    {
    PearsonCorrelation<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Translate" )
    {
    Translate<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "NeighborhoodCorrelation" )
    {
    ImageMetrics<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "NormalizedCorrelation" )
    {
    ImageMetrics<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Demons" )
    {
    ImageMetrics<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Mattes" )
    {
    return ImageMetrics<DIM>(argc, argv);
    }
  if( operation == "MinMaxMean" )
    {
    MinMaxMean<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "PureTissueN4WeightMask" )
    {
    PureTissueN4WeightMask<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "BlobDetector" )
    {
    BlobDetector<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "MatchBlobs" )
    {
    MatchBlobs<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  if( operation == "Project" )
    {
    Project<DIM>(argc, argv);
    return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}


} // namespace ants
