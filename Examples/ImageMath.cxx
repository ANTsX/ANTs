/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RSfile: ImageMath.cxx,v $
  Language:  C++
  Date:      $Date: 2009/06/02 21:51:08 $
  Version:   $Revision: 1.103 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>
#include <vnl/vnl_inverse.h>
#include "antsAllocImage.h"
#include "antsSCCANObject.h"
#include "itkAlternatingValueDifferenceImageFilter.h"
#include "itkAlternatingValueSimpleSubtractionImageFilter.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkArray.h"
#include "itkAverageOverDimensionImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkCompositeValleyFunction.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkDiffusionTensor3D.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHistogram.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImage.h"
#include "itkImageClassifierBase.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageGaussianModelEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkListSample.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
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
#include "itkNormalVariateGenerator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSampleToHistogramFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkSize.h"
#include "itkSphereSpatialFunction.h"
#include "itkSplitAlternatingTimeSeriesImageFilter.h"
#include "itkSTAPLEImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTDistribution.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTranslationTransform.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
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

#include "ReadWriteImage.h"
#include "TensorFunctions.h"
#include "antsMatrixUtilities.h"
#include "../Temporary/itkFastMarchingImageFilter.h"

namespace ants
{
template <class T>
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

template <class T>
std::string ants_to_string(T t)
{
  std::stringstream istream;

  istream << t;
  return istream.str();
}

std::string ANTSOptionName(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            name = std::string( filename, 0, pos );

  return name;
}

std::string ANTSOptionValue(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            value = std::string( filename, pos + 1, filename.length() );

  return value;
}

std::string ANTSGetFilePrefix(const char *str)
{
  const std::string            filename = str;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );
#if 0 //HACK:  This does nothing useful
  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      //extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    }
#endif
  return filepre;
}

template <unsigned int ImageDimension>
int FrobeniusNormOfMatrixDifference(int argc, char *argv[])
{
  if( argc < 6 )
    {
    antscout << " FrobeniusNormOfMatrixDifference: too few options " << std::endl;
    antscout << " ImageMath 3 out FrobeniusNormOfMatrixDifference aff1.mat aff2.mat " << std::endl;
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
  catch( itk::ExceptionObject & excp )
    {
    ::ants::antscout << "no transformation1 that can be read" << fn1 << std::endl;
    return 0;
    }
  typename TransformReaderType::Pointer transformReader2 = TransformReaderType::New();
  transformReader2->SetFileName(  fn2.c_str()  );
  try
    {
    transformReader2->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    ::ants::antscout << "no transformation2 that can be read" << fn2 << std::endl;
    return 0;
    }
  typename AffineTransformType::Pointer aff1 =
    dynamic_cast<AffineTransformType *>( (transformReader1->GetTransformList() )->front().GetPointer() );
  typename AffineTransformType::Pointer aff2 =
    dynamic_cast<AffineTransformType *>( (transformReader2->GetTransformList() )->front().GetPointer() );

  typename AffineTransformType::MatrixType::InternalMatrixType diffmat = aff2->GetMatrix().GetVnlMatrix()
    - aff1->GetMatrix().GetVnlMatrix();
  ::ants::antscout << diffmat.frobenius_norm() << std::endl;
  return 0;
}

template <unsigned int ImageDimension>
int GetLargestComponent(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int           argct = 2;
  const std::string   outname = std::string(argv[argct]);
  argct += 2;
  std::string   fn1 = std::string(argv[argct]);   argct++;
  unsigned long smallest = 50;
  if( argc > argct )
    {
    smallest = atoi(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  // compute the voxel volume
  typename ImageType::SpacingType spacing = image1->GetSpacing();
  float volumeelement = 1.0;
  for( unsigned int i = 0;  i < spacing.Size(); i++ )
    {
    volumeelement *= spacing[i];
    }

  typedef float InternalPixelType;
  //  typedef unsigned long PixelType;
  //  typedef Image<PixelType,ImageDimension>  labelimagetype;
  typedef itk::Image<unsigned long, ImageDimension>                          labelimagetype;
  typedef ImageType                                                          InternalImageType;
  typedef ImageType                                                          OutputImageType;
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
  catch( itk::ExceptionObject & excep )
    {
    antscout << "Relabel: exception caught !" << std::endl;
    antscout << excep << std::endl;
    }

  //  WriteImage<ImageType>(relabel->GetOutput(),outname.c_str());
  //  return 0;
  typename ImageType::Pointer Clusters = MakeNewImage<ImageType>(relabel->GetOutput(), 0);
  // typename ImageType::Pointer Clusters=relabel->GetOutput();
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( relabel->GetOutput(),  relabel->GetOutput()->GetLargestPossibleRegion() );

  float maximum = relabel->GetNumberOfObjects();
  antscout << " #ob " << maximum << std::endl;
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

  antscout << " max float size "
           <<  (maximgval
       * volumeelement) << " long-size: " << (unsigned long) (maximgval * volumeelement)  << std::endl;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  typedef unsigned long                                                    ULPixelType;
  typedef itk::Image<ULPixelType, ImageDimension>                          labelimagetype;
  typedef ImageType                                                        InternalImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                     fIterator;
  typedef itk::ImageRegionIteratorWithIndex<labelimagetype>                labIterator;
  typedef itk::ConnectedComponentImageFilter<ImageType, labelimagetype>    FilterType;
  typedef itk::RelabelComponentImageFilter<labelimagetype, labelimagetype> RelabelType;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  maskfn = std::string(argv[argct]);   argct++;
  unsigned int minclustersize = 50;
  if( argc > argct )
    {
    minclustersize = atoi( argv[argct] );
    }
  typename ImageType::Pointer image = NULL;
  ReadImage<ImageType>(image, fn1.c_str() );
  typename ImageType::Pointer mask = NULL;
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
  catch( itk::ExceptionObject & excep )
    {
    ::ants::antscout << "Relabel: exception caught !" << std::endl;
    ::ants::antscout << excep << std::endl;
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
      float         vox = mask->GetPixel(vfIter.GetIndex() );
      if( vox >= 0  )
        {
        const unsigned long clustersize = histogram[(unsigned long)(relabel->GetOutput()->GetPixel(mIter.GetIndex() ) )];
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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>                       OutImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int slice = atoi(argv[argct]);   argct++;
  antscout << " Extract slice " << slice << " from dimension" << ImageDimension << std::endl;
  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>     ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

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
    antscout << " max slice number is " << timedims << std::endl;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  float replaceValue = 0.0;
  if( argc > 4 )
    {
    replaceValue = atof( argv[4] );
    }

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    float val = vfIter2.Get();
    if( val != val )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       percentofmean = 1.0;
  if( argc > argct )
    {
    percentofmean = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  double        mean = 0, max = -1.e9, min = 1.e9;
  unsigned long ct = 0;
  Iterator      vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    double val = vfIter2.Get();
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
    mean /= (float)ct;
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
int FlattenImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       percentofmax = 1.0;
  if( argc > argct )
    {
    percentofmax = atof(argv[argct]); argct++;
    }

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  double        max = -1.e9, min = 1.e9;
  Iterator      vfIter2( image1,  image1->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    double val = vfIter2.Get();
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
    double val = vfIter2.Get();
    if( val > max * percentofmax )
      {
      val = (max * percentofmax);
      }
    out->SetPixel(vfIter2.GetIndex(), val);
    }

  antscout << " Flattening to :  " << percentofmax << std::endl;
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
    antscout << " need more args -- see usage   " << std::endl
             <<
      " ImageMath 3 outimage.nii.gz  TruncateImageIntensity inputImage  {lowerQuantile=0.025} {upperQuantile=0.975}  {numberOfBins=65}  {binary-maskImage} "
             << std::endl;  throw std::exception();
    }

  unsigned int argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  float        lo = 0.025;
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
    hi = 1.0 - lo;
    } argct++;
  unsigned int numberOfBins = 64;
  if( argc > argct )
    {
    numberOfBins = atoi(argv[argct]);
    }
  argct++;

  //  antscout << " bin " << numberOfBins << " lo " << lo << " Hi " << hi << std::endl;

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
      antscout << " can't read mask " << std::endl;
      mask = NULL;
      }
    ;
    }
  //  antscout << " Mask " << std::endl;
  if( mask.IsNull() )
    {
    mask = AllocImage<ImageType>( image,
                                  itk::NumericTraits<PixelType>::One);
    }

  //  antscout << " iterate " << std::endl;

  itk::ImageRegionIterator<RealImageType> ItI( image,
                                               image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
                                           mask->GetLargestPossibleRegion() );

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();
  ItM.GoToBegin();
  for( ItI.GoToBegin(); !ItI.IsAtEnd();  ++ItI )
    {
    //  antscout << " ind " << ItI.GetIndex() << std::endl;
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
      ItM.Set( itk::NumericTraits<PixelType>::One );
      }
    else
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    if( vnl_math_isnan( ItI.Get() ) || vnl_math_isinf( ItI.Get() ) )
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    ++ItM;
    }
  //  antscout << " label " << std::endl;
  typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( image );
  stats->SetLabelInput( mask );
  stats->SetUseHistograms( true );
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();
  //  antscout << " labeld " << std::endl;
  typedef typename HistogramGeneratorType::HistogramType HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  double lowerQuantile = histogram->Quantile( 0, lo );
  double upperQuantile = histogram->Quantile( 0, hi );

  antscout << "Lower quantile: " << lowerQuantile << std::endl;
  antscout << "Upper quantile: " << upperQuantile << std::endl;
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if( ItI.Get() <  lowerQuantile )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  unsigned int nx = atoi(argv[argct]);   argct++;

  unsigned int numberofimages = 0;

  typename ImageType::Pointer averageimage = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::SizeType size;
  double meanval = 1;
  size.Fill(0);
  unsigned int bigimage = 0;
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        bigimage = j;
        antscout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  ReadImage<ImageType>(image2, argv[bigimage]);

  antscout << " largest image " << size << std::endl;

/** declare the tiled image */
  unsigned int xsize = size[0];
  unsigned int ysize = size[1];
  typename ImageType::SizeType tilesize;
  unsigned int ny = (unsigned int)( (float)numberofimages / (float)nx + 0.5);
  if( nx * ny < numberofimages )
    {
    ny++;
    }
  antscout << " nx " << nx << " ny " << ny << std::endl;
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
        meanval += vfIter2.Get();
        ct++;
        }
      if( ct > 0 )
        {
        meanval /= (float)ct;
        }
      if( meanval <= 0 )
        {
        meanval = 1.0;
        }
      }

    imagexct = imagecount % nx;
    imageyct = imagecount / nx;
    antscout << "doing " << fn << "  " << imagecount << " x " << imagexct <<  " y " << imageyct << std::endl;
    imagecount++;
    Iterator vfIter( image2,  image2->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      typename ImageType::IndexType locind = vfIter.GetIndex();
      typename ImageType::IndexType globind;
      globind[0] = size[0] * imagexct + locind[0];
      globind[1] = size[1] * imageyct + locind[1];
      double val =  vfIter.Get() / meanval;
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

  antscout << " writing output ";
  WriteImage<ByteImageType>( rescaler->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int ConvertLandmarkFile(unsigned int argc, char *argv[])
{
  unsigned int argct = 2;

  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
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

  antscout << "Number of labels: " << reader->GetNumberOfLabels() << std::endl;
  antscout << "Labels: ";
  for( unsigned int i = 0; i < reader->GetNumberOfLabels(); i++ )
    {
    antscout << reader->GetLabelSet()->operator[](i) << " ";
    }
  antscout << std::endl;

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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  unsigned int argct = 2;
  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string maskfn = std::string(argv[argct]); argct++;
  antscout << " file name " << maskfn << std::endl;
  typename ImageType::Pointer mask = NULL;
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
    xslice = atoi(argv[argct]);
    }
  argct++;
  unsigned int yslice = size[1] / 2;
  if( argc > argct )
    {
    yslice = atoi(argv[argct]);
    }
  argct++;
  unsigned int zslice = size[2] / 2;
  if( argc > argct )
    {
    zslice = atoi(argv[argct]);
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
  antscout << " allocate matrix " << tilesize << std::endl;
  typename MatrixImageType::RegionType region;
  region.SetSize( tilesize );

  typename MatrixImageType::Pointer matimage =
    AllocImage<MatrixImageType>(region);

  unsigned int lowgetridof = (unsigned int) (clamppercent1 * 256);
  unsigned int higetridof = (unsigned int) (256 - clamppercent2 * 256);
  //  antscout << " get rid of " << getridof << std::endl;
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
  antscout << " writing output ";
  WriteImage<ByteImageType>( rescaler2->GetOutput(), outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int ConvertVectorToImage(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  typedef itk::ImageRegionIteratorWithIndex<MatrixImageType>              vIterator;

  int argct = 2;
  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string maskfn = std::string(argv[argct]); argct++;
  std::string vecfn = std::string(argv[argct]); argct++;
  typename ImageType::Pointer mask = NULL;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  typename MatrixImageType::Pointer vecimg = NULL;
  ReadImage<MatrixImageType>(vecimg, vecfn.c_str() );
  unsigned long voxct = 0, mct = 0;
  Iterator      mIter( mask, mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5 )
      {
      mct++;
      }
    }
  vIterator vIter(vecimg, vecimg->GetLargestPossibleRegion() );
  for(  vIter.GoToBegin(); !vIter.IsAtEnd(); ++vIter )
    {
    voxct++;
    }

  antscout << " vct " << voxct << " mct " << mct << std::endl;

  typename ImageType::Pointer outimage = NULL;
  ReadImage<ImageType>(outimage, maskfn.c_str() );
  outimage->FillBuffer(0);

  vIter.GoToBegin();
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5 )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator iter( image1,  image1->GetLargestPossibleRegion() );
  for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    double r = (   (double)rand()  / ( (double)(RAND_MAX) + (double)(1) ) ) - 0.5;
    iter.Set(iter.Get() + r * noiselevel);
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
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

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
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
      if( fabs(iter.Get() - value) < tol )
        {
        antscout << iter.GetIndex() << std::endl;
        ct++;
        }
      }
    else if( image2->GetPixel(iter.GetIndex() ) > 0 &&  fabs(iter.Get() - value) < tol )
      {
      antscout << iter.GetIndex() << std::endl;
      ct++;
      }
    }
  antscout << ct <<  " voxels have the value " << value << std::endl;
  return 0;
}

template <unsigned int ImageDimension>
int SetOrGetPixel(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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
    usephyspace = atoi(argv[argct]); argct++;
    }
  bool get = false;
  if( strcmp(Get.c_str(), "Get") == 0 )
    {
    get = true;
    }
  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    ReadImage<ImageType>(image2, fn1.c_str() );
    }
  if( !image1 )
    {
    antscout << " no image ! " << std::endl; throw std::exception();
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
  //  antscout << " use phy " << usephyspace << " " << indx << " " << indy << " " << indz << std::endl;
  //  antscout << " Ind " << index << std::endl;
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
      //antscout << " GetValue at " << index << " is " << image1->GetPixel(index) << std::endl;
      antscout << image1->GetPixel(index) << std::endl;
      }
    else
      {
      //  antscout << " SetValue at " << index << " value " << value << " replaces " <<  image1->GetPixel(index)
      //         << std::endl;
      image2->SetPixel(index, value);
      WriteImage<ImageType>(image2, outname.c_str() );
      }
    }
  else
    {
    antscout << "NA" << std::endl;
    }

  return 0;
}

template <unsigned int ImageDimension>
int HistogramMatching(int argc, char * argv[])
{
  typedef float                                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension>                   ImageType;
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> MatchingFilterType;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = std::string(argv[argct]);   argct++;
  long        bins = 255;
  if( argc > argct )
    {
    bins = atoi(argv[argct]);
    }
  argct++;
  long points = 64;
  if( argc > argct )
    {
    points = atoi(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer source;
  ReadImage<ImageType>(source, fn1.c_str() );

  typename ImageType::Pointer reference;
  ReadImage<ImageType>(reference, fn2.c_str() );

  typename MatchingFilterType::Pointer match = MatchingFilterType::New();
  match->SetSourceImage(source);
  match->SetReferenceImage(reference);
  match->SetNumberOfHistogramLevels(bins);
  match->SetNumberOfMatchPoints(points);
  match->Update();

  WriteImage<ImageType>(match->GetOutput(), outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int PadImage(int /*argc */, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  const float padvalue = atof(argv[argct]);
  argct +=2;

  typename ImageType::Pointer image1 = NULL;
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
  antscout << " oldsize " << size <<  " newsize " << newsize << std::endl;
  newregion.SetSize(newsize);
  newregion.SetIndex(image1->GetLargestPossibleRegion().GetIndex() );

  typename ImageType::Pointer padimage =
    AllocImage<ImageType>(newregion,
                          image1->GetSpacing(),
                          origin2,
                          image1->GetDirection(), 0);

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

  antscout << " pre " << point1 << " pad " << pointpad << std::endl;
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
  typedef itk::ImageFileReader<ImageType>       ReaderType;
  typedef itk::ImageFileWriter<ImageType>       WriterType;

  int         argct = 2;
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

  typename ImageType::Pointer inputImage = NULL;
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
int CenterImage2inImage1(int argc, char *argv[])
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = NULL;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  image1->FillBuffer(0);
  typename ImageType::Pointer image2 = NULL;
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
      cm_point[d] += point[d] * iter.Get();
      }
    iweight += iter.Get();
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
int TimeSeriesDisassemble(int argc, char *argv[])
{
  if( argc <= 4 )
    {
    antscout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>     ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

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
  antscout << " Extract " << n_sub_vols << " subvolumes " << std::endl;

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
    out << (100 + i);
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
    antscout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension - 1>    ImageType;
  typedef itk::Image<PixelType, ImageDimension>        OutImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  float       time = atof( argv[argct] );           argct++;
  float       origin = atof( argv[argct] );         argct++;

  typename OutImageType::Pointer outimage = OutImageType::New();

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

  antscout << " Merging " << argc - 6 << " subvolumes " << std::endl;
  antscout << " time spacing: " << time << std::endl;
  antscout << " time origin: " << origin << std::endl;
  for( int i = 6; i < argc; i++ )
    {
    typename ImageType::Pointer image1 = NULL;
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
    antscout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int n_sub_vols = atoi(argv[argct]);   argct++;
  antscout << " Extract " << n_sub_vols << " subvolumes " << std::endl;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>     ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  else
    {
    antscout << "Failed to read input image" << std::endl;
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

  int         argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string fn1 = std::string(argv[argct++]);
  bool        mean = false;

  typename ImageFilterType::Pointer filter = ImageFilterType::New();

  if ( argc >= 6 )
    {
    if ( atoi(argv[argct++]) > 0 )
      {
      mean = true;
      }
    }

  typename InputImageType::Pointer image1 = NULL;
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

  if ( mean )
    {
    typename MeanFilterType::Pointer meanFilter = MeanFilterType::New();
    meanFilter->SetInput( filter->GetOutput() );
    meanFilter->SetAveragingDimension( ImageDimension-1 );
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

  typedef itk::BSplineInterpolateImageFunction<InputImageType, double > BSplineInterpolatorType;
  typedef typename BSplineInterpolatorType::Pointer                     BSplineInterpolatorPointerType;

  typedef itk::WindowedSincInterpolateImageFunction<InputImageType, 4>   SincInterpolatorType;
  typedef typename SincInterpolatorType::Pointer                         SincInterpolatorPointerType;

  int         argct = 2;
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
      const unsigned int SincRadius = 4;
      antscout << "Using sinc interpolation" << std::endl;
      SincInterpolatorPointerType labelInterp = SincInterpolatorType::New();
      SincInterpolatorPointerType controlInterp = SincInterpolatorType::New();
      filter->SetControlInterpolator( controlInterp );
      filter->SetLabelInterpolator( labelInterp );
      filter->SetIndexPadding( SincRadius );
      }
    else if ( strcmp( "bspline", interp.c_str() ) == 0 )
      {
      antscout << "Using bspline interpolation of order 3" << std::endl;
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
      antscout << "Using linear interpolation" << std::endl;
      }

    }


  bool mean = false;
  if( argc >= 7 )
    {
    if( atoi(argv[argct++]) > 0 )
      {
      mean = true;
      }
    }

  typename InputImageType::Pointer image1 = NULL;
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
  if ( mean )
    {
    typename MeanFilterType::Pointer meanFilter = MeanFilterType::New();
    meanFilter->SetInput( filter->GetOutput() );
    meanFilter->SetAveragingDimension( ImageDimension-1 );
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
  typedef float                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension >  ImageType;

  typedef itk::SplitAlternatingTimeSeriesImageFilter<ImageType, ImageType>
    ImageFilterType;

  if ( argc < 5 )
    {
    antscout << "Usage: ImageMath 4 split.nii.gz AlternatingTimeSeriesExtraction time.nii.gz" << std::endl;
    return 1;
    }

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;

  std::string::size_type idx;
  idx = outname.find_first_of('.');

  std::string basename = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  std::string zero( "0" );
  std::string one( "1" );

  std::string outname0 = basename + zero + extension;
  std::string outname1 = basename + one + extension;

  typename ImageFilterType::Pointer filter = ImageFilterType::New();
  typename ImageType::Pointer image1 = NULL;
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
int AverageOverDimension(int argc, char *argv[])
{
  typedef float                                   PixelType;
  typedef itk::Image<PixelType, ImageDimension >  ImageType;
  typedef itk::Image<PixelType, ImageDimension-1> AverageImageType;

  typedef itk::AverageOverDimensionImageFilter<ImageType, AverageImageType>
    ImageFilterType;

  if ( argc < 6 )
    {
    antscout << "Usage: ImageMath 4 average.nii.gz AverageOverDimension time.nii.gz dimension" << std::endl;
    return EXIT_FAILURE;
    }
  int          argct = 2;
  const std::string  outname = std::string(argv[argct++]);
  argct++;
  std::string  fn1 = std::string(argv[argct++]);
  unsigned int dim = atoi( argv[argct++] );

  typename ImageType::Pointer image1 = NULL;
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
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        InputImageType;
  typedef itk::Image<unsigned int, ImageDimension>     LabelImageType;

  typedef itk::MinimumMaximumImageCalculator<LabelImageType>  LabelCalculatorType;
  typedef itk::ants::antsSCCANObject<InputImageType, double>  SCCANType;

  typedef typename SCCANType::MatrixType         MatrixType;
  typedef typename SCCANType::VectorType         VectorType;

  if(argc < 6)
    {
    return EXIT_FAILURE;
    }
  int         argct = 2;
  const std::string outname = std::string(argv[argct++]);
  argct++;
  std::string labelName = std::string(argv[argct++]);
  std::string timeName  = std::string(argv[argct++]);

  // FIXME - add option for multi input for combined CCA

  typename LabelImageType::Pointer labels = NULL;
  ReadImage<LabelImageType>( labels, labelName.c_str() );

  typename InputImageType::Pointer time = NULL;
  ReadImage<InputImageType>( time, timeName.c_str() );

  typename LabelCalculatorType::Pointer calc = LabelCalculatorType::New();
  calc->SetImage( labels );
  calc->ComputeMaximum();
  unsigned int nLabels = calc->GetMaximum();
  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  std::cout << "Examining " << nLabels << " regions, covering "
            << nVoxels << " voxels " << std::endl;

  //unsigned int labelCounts[nLabels];
  unsigned int *labelCounts = new unsigned int [nLabels] ;

  for (unsigned int i=0; i<nLabels; i++)
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;

    labelCounts[i] = 0;
    for ( unsigned int v=0; v<nVoxels; v++)
      {
      idx[0] = v;
      if ( labels->GetPixel(idx) == (i+1) )
        {
        ++labelCounts[i];
        }
      }
    }

  typename InputImageType::Pointer connmat = InputImageType::New();
  typename InputImageType::RegionType region;
  region.SetSize(0,nLabels);
  region.SetSize(1,nLabels);
  connmat->SetRegions( region );
  connmat->Allocate();
  connmat->FillBuffer(0);

  // Coorelation parameters
  bool robust = false;
  unsigned int iterct = 20;
  bool useL1 = false;
  float gradstep = vnl_math_abs( useL1 );
  bool keepPositive = false;
  float sparsity = 1.0;
  unsigned int minClusterSize = 1;
  unsigned int minRegionSize = 1;

  // used to rankify matrices if using robust
  typename SCCANType::Pointer cca_rankify = SCCANType::New();

  for (unsigned int i=0; i<nLabels; i++)
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;

    MatrixType P(nTimes, labelCounts[i], 0.0);

    unsigned int iCount = 0;
    for ( unsigned int v=0; v<nVoxels; v++)
      {
      idx[0] = v;
      typename InputImageType::IndexType timeIdx;
      timeIdx[1] = v;

      if ( labels->GetPixel(idx) == (i+1) )
        {
        for ( unsigned int t=0; t<nTimes; t++)
          {
          timeIdx[0] = t;
          P(t,iCount) = time->GetPixel(timeIdx);
          }
        ++iCount;
        }
      }

    if ( robust && ( labelCounts[i] >= minRegionSize)  )
      {
      P = cca_rankify->RankifyMatrixColumns(P);
      }


    if ( labelCounts[i] >= minRegionSize )
      {
      for ( unsigned int j=i+1; j<nLabels; j++)
        {
        MatrixType Q(nTimes, labelCounts[j], 0.0);
        typename LabelImageType::IndexType idx2;
        idx2[1] = 0;

        unsigned int jCount = 0;
        for (unsigned int v2=0; v2<nVoxels; v2++)
          {
          idx2[0] = v2;
          typename InputImageType::IndexType timeIdx2;
          timeIdx2[1] = v2;

          if ( labels->GetPixel(idx2) == (j+1) )
            {
            for ( unsigned int t2=0; t2<nTimes; t2++)
              {
              timeIdx2[0] = t2;
              Q(t2,jCount) = time->GetPixel(timeIdx2);
              }
            ++jCount;
            }
          }

        if ( robust )
          {
          Q = cca_rankify->RankifyMatrixColumns(Q);
          }

        if ( labelCounts[j] >= minRegionSize)
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
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        InputImageType;
  typedef itk::Image<unsigned int, ImageDimension>     LabelImageType;

  typedef itk::MinimumMaximumImageCalculator<LabelImageType>  LabelCalculatorType;
  typedef itk::ants::antsSCCANObject<InputImageType, double>  SCCANType;

  typedef typename SCCANType::MatrixType         MatrixType;
  typedef typename SCCANType::VectorType         VectorType;

  if(argc < 6)
    {
    return EXIT_FAILURE;
    }

  int         argct = 2;
  std::string outname = std::string(argv[argct++]);
  std::string operation = std::string(argv[argct++]);
  std::string labelName = std::string(argv[argct++]);
  std::string timeName  = std::string(argv[argct++]);

  unsigned int minRegionSize = 3;
  
  if ( argc > 6 )
    minRegionSize = atoi( argv[argct++] );

  std::cout << "min region = " << minRegionSize << std::endl;

  // FIXME - add option for multi input for combined CCA
  
  typename LabelImageType::Pointer labels = NULL;
  ReadImage<LabelImageType>( labels, labelName.c_str() );

  typename InputImageType::Pointer time = NULL;
  ReadImage<InputImageType>( time, timeName.c_str() );
  
  typename LabelCalculatorType::Pointer calc = LabelCalculatorType::New();
  calc->SetImage( labels );
  calc->ComputeMaximum();
  unsigned int nLabels = calc->GetMaximum();
  unsigned int nVoxels = labels->GetLargestPossibleRegion().GetSize()[0];
  unsigned int nTimes = time->GetLargestPossibleRegion().GetSize()[0];

  std::cout << "Examining " << nLabels << " regions, covering " 
            << nVoxels << " voxels " << std::endl;

  VectorType labelCounts( nLabels, 0 );
  

  typename InputImageType::Pointer connmat = InputImageType::New();
  typename InputImageType::RegionType region;
  region.SetSize(0,nLabels);
  region.SetSize(1,nLabels);
  connmat->SetRegions( region );
  connmat->Allocate();
  connmat->FillBuffer(-1);

  MatrixType timeSig( nLabels, nTimes, 0.0 );
  for ( unsigned int i=0; i<nLabels; i++)
    {
    typename LabelImageType::IndexType idx;
    idx[1] = 0;
    
    for ( unsigned int v=0; v<nVoxels; v++)
      {      
      idx[0] = v;
      if ( labels->GetPixel(idx) == (i+1) )
        {
        labelCounts[i]++;

        typename InputImageType::IndexType timeIdx;
        timeIdx[1] = v;
        for ( unsigned int t=0; t<nTimes; t++)
          {
          timeIdx[0] = t;
          timeSig(i,t) += time->GetPixel(timeIdx);
          }
        }
      }
    }

  for ( unsigned int i=0; i<nLabels; i++ )
    {
    for ( unsigned int j=0; j<nTimes; j++ )
      {
      timeSig(i,j) /= labelCounts[i];
      }
    }

  std::cout << "Building matrices..." << std::endl;

  for (unsigned int i=0; i<nLabels; i++)
    {
    for ( unsigned int j=(i+1); j<nLabels; j++ )      {
      
      if ( (labelCounts[i] > minRegionSize) && (labelCounts[j] > minRegionSize ) )
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
        if ( ! vnl_math_isfinite( corr ) )
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

std::cout << "Write output " << outname << std::endl;
WriteImage<InputImageType>(connmat, outname.c_str() );
  
  return 0;
}


template <unsigned int ImageDimension>
int PASLQuantifyCBF(int argc, char *argv[])
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> TimeImageType;
  typedef itk::Image<PixelType, ImageDimension-1> ImageType;

  typedef itk::PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TimeImageType, ImageType, TimeImageType>
    FilterType;

  int         argct = 2;
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

  typename TimeImageType::Pointer diff = NULL;
  if( fn1.length() > 3 )
    {
    ReadImage<TimeImageType>(diff, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename ImageType::Pointer m0 = NULL;
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
  if ( argc < 6 )
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


  typename TimeImageType::Pointer diff = NULL;
  if( fn1.length() > 3 )
    {
    ReadImage<TimeImageType>(diff, fn1.c_str() );
    }
  else
    {
    return 1;
    }

  typename ImageType::Pointer m0 = NULL;
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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int k_neighbors = atoi(argv[argct]);   argct++;
  typename ImageType::Pointer image1 = NULL;
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

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;
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
    antscout << "Raw_Leverage,K_Neighbors_Distance" <<  std::endl;
    logfile << "Raw_Leverage,K_Neighbors_Distance" <<  std::endl;
    for( unsigned int t = 0; t < timedims; t++ )
      {
      antscout <<  mLeverage(t) << "," << kDistance(t) << std::endl;
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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int         argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string ext = itksys::SystemTools::GetFilenameExtension( outname );
  if( strcmp(ext.c_str(), ".csv") != 0 )
    {
    antscout << " must use .csv as output file extension " << std::endl;
    return EXIT_FAILURE;
    }
  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string maskfn = std::string(argv[argct]);   argct++;
  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer mask = NULL;
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
    if( mIter.Get() >= 0.5 )
      {
      voxct++;
      }
    }

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

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

  typedef vnl_vector<Scalar> timeVectorType;
  timeVectorType mSample(timedims, 0);
  typedef itk::Array2D<double> MatrixType;
  std::vector<std::string> ColumnHeaders;
  MatrixType               matrix(timedims, voxct);
  matrix.Fill(0);
  SliceIt vfIter2( outimage, outimage->GetLargestPossibleRegion() );
  voxct = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    OutIndexType ind = vfIter2.GetIndex();
    if( mask->GetPixel(ind) >= 0.5 )
      {
      IndexType tind;
      // first collect all samples for that location
      for( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        tind[i] = ind[i];
        }
      for( unsigned int t = 0; t < timedims; t++ )
        {
        tind[ImageDimension - 1] = t;
        Scalar pix = image1->GetPixel(tind);
        mSample(t) = pix;
        matrix[t][voxct] = pix;
        }
      std::string colname = std::string("V") + ants_to_string<unsigned int>(voxct);
      ColumnHeaders.push_back( colname );
      voxct++;
      } // check mask
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
  catch( itk::ExceptionObject& exp )
    {
    antscout << "Exception caught!" << std::endl;
    antscout << exp << std::endl;
    return EXIT_FAILURE;
    }

  return 0;
}

template <unsigned int ImageDimension>
int PASL(int argc, char *argv[])
{
  if( argc <= 3 )
    {
    antscout << " too few options " << argv[0] << std::endl;
    antscout << argv[0] << " NDImage  Bool_FirstImageIsControl optional-M0mask.nii.gz " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef double                                       RealType;
  typedef vnl_vector<RealType>                         timeVectorType;
  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  bool        firstiscontrol = atoi(argv[argct]);   argct++;
  std::string m0fn = "";
  if( argc > argct )
    {
    m0fn = std::string(argv[argct]);
    }
  argct++;
  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;
  typename OutImageType::Pointer M0image = NULL;

  typedef itk::ExtractImageFilter<ImageType, OutImageType> ExtractFilterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>     ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

  //RealType M0W = 1300; // FIXME
  //RealType TE = 4000;
  //RealType calculatedM0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
  RealType calculatedM0 = 2800; // from "Impact of equilibrium magnetization of blood on ASL quantification" by YChen et al

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
    antscout
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
    RealType      M_0   = M0image->GetPixel( ind ); // FIXME can be taken from an input reference image or defined for
                                                    // each tissue
    bool getCBF = true;
    if( haveM0 && M_0 == 0 )
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
          f =  \frac{      \lambda DeltaM        }  {     2 \alpha M_0 TI_1 exp( - TI_2 / T_{1a} )  }
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
          RealType scaling = 2 * alpha * M_0 * TI_1 * exp( -TI_2 / T_1a );
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
    antscout << " too few options " << argv[0] << std::endl;
    antscout << argv[0] << " NDImage  Bool_FirstImageIsControl optional-M0mask.nii.gz " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef double                                       RealType;
  typedef vnl_vector<RealType>                         timeVectorType;
  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  bool        firstiscontrol = atoi(argv[argct]);   argct++;
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
  typedef itk::ImageRegionIteratorWithIndex<ImageType>     ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType>  SliceIt;

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

  //RealType M0W = 1300; // FIXME
  //RealType TE = 4000;
  //RealType calculatedM0 = 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;
  RealType calculatedM0 = 2800; // from "Impact of equilibrium magnetization of blood on ASL quantification" by YChen et al

  bool haveM0 = true;
  if( m0fn.length() > 3 )
    {
    ReadImage<OutImageType>( M0image, m0fn.c_str() );
    }
  else
    {
    haveM0 = false;
    antscout
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
    RealType      M_0   = M0image->GetPixel( ind ); // FIXME can be taken from an input reference image or defined for
                                                    // each tissue
    bool getCBF = true;
    if( haveM0 && M_0 == 0 )
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
        CBF =  \frac{      \lambda DeltaM     }  {     2 \alpha M_0 T_1a [ exp( - w / T_1a ) - exp( - ( tau + w ) / T_1a )  ] }
        w of 0.7seconds was determined based on the results of experiment (1) for optimal contrast of the GM. Cerebral blood flow was calculated for a single PLD. Quantitative CBF for the whole brain, GM, and WM were tabulated. The bottom slice was excluded from this analysis because it covered only a small part of the cerebrum.
           */
          // f =  \frac{      \lambda DeltaM     }  {     2 \alpha M_0 T_1a [ exp( - w / T_1a ) - exp( - ( tau + w ) /
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
          const RealType scaling = 4 * alpha * M_0 * T_1t
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
    antscout << " too few options " << argv[0] << std::endl;
    antscout << argv[0] << " M0Image.nii.gz M1Image.nii.gz [OptionalMask.nii.gz] " << std::endl;
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
    antscout << " too few options " << std::endl;
    return 1;
    }

  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  fn_label = std::string(argv[argct]);   argct++;
  unsigned int n_comp_corr_vecs = 4; // number of eigenvectors to get from high variance voxels
  if( argc > argct )
    {
    n_comp_corr_vecs = atoi(argv[argct]);
    }
  argct++;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;
  typename OutImageType::Pointer outimage2 = NULL;
  typename OutImageType::Pointer label_image = NULL;
  typename OutImageType::Pointer var_image = NULL;

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

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
  antscout << " read images " << std::endl;
  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  antscout << "timedims " << timedims << " size " << image1->GetLargestPossibleRegion().GetSize() << std::endl;

  // first, count the label numbers
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator   vfIter2( label_image,  label_image->GetLargestPossibleRegion() );
  unsigned long ct_nuis = 0;
  unsigned long ct_vox = 0;
  antscout << " verify input " << std::endl;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( vfIter2.Get() == 1 )      // in brain
      {
      ct_vox++;
      }
    }
  antscout << " counted " << ct_vox << " voxels " <<  std::endl;
  if( ct_vox == 0 )
    {
    antscout << ct_vox << " not enough voxels labeled as gm (or brain) " << std::endl;
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
        total += pix;
        }
      float mean = total / (float)timedims;
      float var = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        var += ( (sample(t) - mean) * (sample(t) - mean) );
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
  antscout << " got var " << std::endl;
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
  antscout << " got var hist " << std::endl;
  float temp = 0;
  float varval_csf = 0;
  for( unsigned int j = 0; j < histsize; j++ )
    {
    temp += varhist(j);
    float varth = (float)j / (float)histsize * maxvar;
    antscout << " j " << j << " temp " << temp << " varth " << varth << std::endl;
    if( temp >= 0.95 && varval_csf <=  0 )
      {
      varval_csf = (float)j * binsize;
      }
    }

  antscout << " maxvar " << maxvar << " varval_csf " << varval_csf << std::endl;
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
  catch( itk::ExceptionObject& exp )
    {
    antscout << "Exception caught!" << std::endl;
    antscout << exp << std::endl;
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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                        PixelType;
  typedef itk::Vector<float, ImageDimension>           VectorType;
  typedef itk::Image<VectorType, ImageDimension>       FieldType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1>    OutImageType;
  typedef typename OutImageType::IndexType             OutIndexType;
  typedef itk::ImageFileReader<ImageType>              readertype;
  typedef itk::ImageFileWriter<ImageType>              writertype;
  typedef typename ImageType::IndexType                IndexType;
  typedef typename ImageType::SizeType                 SizeType;
  typedef typename ImageType::SpacingType              SpacingType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

  typedef double                                            Scalar;
  typedef itk::ants::antsMatrixUtilities<ImageType, Scalar> matrixOpType;
  typename matrixOpType::Pointer matrixOps = matrixOpType::New();

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  fn_label = std::string(argv[argct]);   argct++;
  unsigned int wmlabel = 1;
  unsigned int csflabel = 3;
  if( argc > argct )
    {
    csflabel = atoi(argv[argct]);
    }
  argct++;
  if( argc > argct )
    {
    wmlabel = atoi(argv[argct]);
    }
  argct++;
  std::string::size_type idx;
  idx = outname.find_first_of('.');
  std::string tempname = outname.substr(0, idx);
  std::string extension = outname.substr(idx, outname.length() );

  typename ImageType::Pointer image1 = NULL;
  typename OutImageType::Pointer outimage = NULL;
  typename OutImageType::Pointer outimage2 = NULL;
  typename OutImageType::Pointer label_image = NULL;
  typename OutImageType::Pointer var_image = NULL;

  typedef itk::ImageRegionIteratorWithIndex<ImageType>    ImageIt;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> SliceIt;

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
  antscout << " read images " << std::endl;
  unsigned int timedims = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  antscout << "timedims " << timedims << " size " << image1->GetLargestPossibleRegion().GetSize() << std::endl;

  // first, count the label numbers
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> labIterator;
  labIterator   vfIter2( label_image,  label_image->GetLargestPossibleRegion() );
  unsigned long ct_nuis = 0;
  unsigned long ct_ref = 0;
  unsigned long ct_gm = 0;
  antscout << " verify input " << std::endl;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    if( vfIter2.Get() == csflabel )      // nuisance
      {
      ct_nuis++;
      }
    if( vfIter2.Get() == wmlabel )      // reference
      {
      ct_ref++;
      }
    if( vfIter2.Get() > 0 )      // gm roi
      {
      ct_gm++;
      }
    }
  antscout << " counted " << ct_gm << " gm voxels " << ct_ref << " reference region voxels " << std::endl;
  if( ct_gm == 0 )
    {
    antscout << ct_gm << " not enough voxels labeled as gm (or brain) " << ct_gm << std::endl;
    return 1;
    }
  if( ct_ref == 0 )
    {
    antscout << ct_ref << " not enough voxels labeled as reference region " << std::endl;
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
        Scalar pix = image1->GetPixel(tind);
        smoother(t) = pix;
        total += pix;
        }
      float mean = total / (float)timedims;
      float var = 0;
      for( unsigned int t = 0; t < timedims; t++ )
        {
        var += ( (smoother(t) - mean) * (smoother(t) - mean) );
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
  antscout << " got var " << std::endl;
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
  antscout << " got var hist " << std::endl;
  float temp = 0;
  float varval_csf = 0;
  for( unsigned int j = 0; j < histsize; j++ )
    {
    temp += varhist(j);
    float varth = (float)j / (float)histsize * maxvar;
    antscout << " j " << j << " temp " << temp << " varth " << varth << std::endl;
    if( temp >= 0.95 && varval_csf <=  0 )
      {
      varval_csf = (float)j * binsize;
      }
    }

  antscout << " maxvar " << maxvar << " varval_csf " << varval_csf << std::endl;
  //  WriteImage<OutImageType>(var_image,"varimage.nii.gz");
  //
  ct_nuis = 0;
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    // OutIndexType ind = vfIter2.GetIndex();
    if( vfIter2.Get() == csflabel )      // reference
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
    if( vfIter2.Get() == csflabel )      // reference
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
    if( vfIter2.Get() == wmlabel )      // reference
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
  catch( itk::ExceptionObject& exp )
    {
    antscout << "Exception caught!" << std::endl;
    antscout << exp << std::endl;
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
    antscout << " CompCorr Error exiting " << std::endl; throw std::exception();
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
  antscout << "write results" << std::endl;
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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = NULL;
  if( fn1.length() > 3 )
    {
    ReadImage<ImageType>(image1, fn1.c_str() );
    }
  typename ImageType::Pointer image2 = NULL;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(image2, fn2.c_str() );
    }

  typename ImageType::PointType origin2 = image1->GetOrigin();

  typename ImageType::SizeType size = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::SizeType newsize = image1->GetLargestPossibleRegion().GetSize();
  typename ImageType::RegionType newregion;
  // determine new image size

  newsize[ImageDimension
          - 1] =
    (unsigned int)newsize[ImageDimension - 1] + image2->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  antscout << " oldsize " << size <<  " newsize " << newsize << std::endl;
  newregion.SetSize(newsize);
  newregion.SetIndex(image1->GetLargestPossibleRegion().GetIndex() );

  typename ImageType::Pointer padimage =
    AllocImage<ImageType>(newregion,
                          image1->GetSpacing(),
                          origin2,
                          image1->GetDirection(),
                          0);

  typename ImageType::IndexType index; index.Fill(0);
  typename ImageType::IndexType index2; index2.Fill(0);
  typename ImageType::PointType point1, pointpad;
  image1->TransformIndexToPhysicalPoint(index, point1);
  padimage->TransformIndexToPhysicalPoint(index2, pointpad);
  //  antscout << " pre " << point1 << " pad " << pointpad << std::endl;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin2[i] += (point1[i] - pointpad[i]);
    }
  padimage->SetOrigin(origin2);
  const float padvalue = image1->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  Iterator iter( padimage,  padimage->GetLargestPossibleRegion() );
  for( iter.GoToBegin(); !iter.IsAtEnd();
       ++iter )
    {
    const typename ImageType::IndexType & oindex = iter.GetIndex();
    typename ImageType::IndexType padindex = iter.GetIndex();
    bool isinside = false;
    if( oindex[ImageDimension - 1]  >= padvalue )
      {
      isinside = true;
      }
    if( isinside )
      {
      float shifted = ( (float)oindex[ImageDimension - 1] - padvalue);
      padindex[ImageDimension - 1] = (unsigned int)shifted;
      padimage->SetPixel(oindex, image2->GetPixel(padindex) );
      }
    else
      {
      padimage->SetPixel(oindex, image1->GetPixel(oindex) );
      }
    }

  WriteImage<ImageType>(padimage, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int MakeImage(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  unsigned int sizevalx = atoi(argv[argct]);   argct++;
  unsigned int sizevaly = 0;
  if( argc > argct )
    {
    sizevaly = atoi(argv[argct]); argct++;
    }
  unsigned int sizevalz = 0;
  if( argc > argct && ImageDimension == 3 )
    {
    sizevalz = atoi(argv[argct]); argct++;
    }

  typename ImageType::SizeType size;
  size[0] = sizevalx;
  size[1] = sizevaly;
  if( ImageDimension == 3 )
    {
    size[2] = sizevalz;
    }
  typename ImageType::RegionType newregion;
  antscout << " size " << size << std::endl;
  newregion.SetSize(size);

  typename ImageType::Pointer padimage = AllocImage<ImageType>(newregion, 0);
  WriteImage<ImageType>(padimage, outname.c_str() );

  return 0;
}

template <class TImage>
typename TImage::Pointer
LabelSurface(typename TImage::Pointer input, typename TImage::Pointer input2  )
{
  antscout << " Label Surf " << std::endl;
  typedef TImage ImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  // ORIENTATION ALERT -- original code set size,spacing,and origin
  // without setting orientation
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
    typename TImage::IndexType ind2;
    if( p >= 0.5 )
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
        bool secondval = true;
        if( input2 )
          {
          if( input2->GetPixel(ind2) >= 0.5 )
            {
            secondval = true;
            }
          }
        if( GHood.GetPixel(i) < 0.5 && dist <  2. && secondval  )
          {
          atedge = true;
          }
        }
      if( atedge && p >=  0.5 )
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
int FitSphere(int argc, char *argv[])
{
  if( argc <= 2 )
    {
    antscout << " too few options " << std::string(argv[1]) << std::endl;
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

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer radimage = NULL;
  typename ImageType::Pointer radimage2 = NULL;
  typename ImageType::Pointer priorimage = NULL;
  typename ImageType::Pointer wmimage = NULL;
  if (fn2.length() > 3)   ReadImage<ImageType>(wmimage, fn2.c_str());
  antscout <<"  read " << fn1 << " MXR " << MaxRad << std::endl;
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
  typename ScalarInterpolatorType::Pointer winterp=NULL;
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
  antscout <<"  Begin " << std::endl;
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
      //          antscout << " GMT " << gmtotal << " WMT " << wmtotal << " dist " << cmdist << std::endl;
  float gmrad=pow( 3.*gvol/(4.*pi) , 1./3.);
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
              //          antscout << " Ind " <<  ind << " : " <<  bestrad << " tardist " << tardist << " gct " << goodct <<" pos " << possct << " dist " << dist << " ind2 " << ind2 << std::endl;
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
          antscout <<" prog " << (float)npx/(float)numpx << std::endl;
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
  antscout << " Best " << bestind << " gbr " << globalbestrad << std::endl;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string operation = std::string(argv[argct]);  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::Pointer varimage = NULL;

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
    antscout << " trans " << m_Transform0->GetParameters() << " Nspc " << image2->GetSpacing() << std::endl;
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
    volumeelement *= varimage->GetSpacing()[i];
    }

  float    result = 0;
  unsigned long ct = 0;
  Iterator vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
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
      result = pow(pix1, pix2);
      }
    else if( strcmp(operation.c_str(), "exp") == 0 )
      {
      result = exp(pix1 * pix2);
      }
    else if( strcmp(operation.c_str(), "abs") == 0 )
      {
      result = fabs(pix1);
      }
    else if( strcmp(operation.c_str(), "addtozero") == 0 && pix1 == 0 )
      {
      result = pix1 + pix2;
      }
    else if( strcmp(operation.c_str(), "addtozero") == 0 && pix1 != 0 )
      {
      result = pix1;
      }
    else if( strcmp(operation.c_str(), "overadd") == 0 && pix2 != 0 )
      {
      result = pix2;
      }
    else if( strcmp(operation.c_str(), "overadd") == 0 )
      {
      result = pix1;
      }
    else if( strcmp(operation.c_str(), "Decision") == 0 )
      {
      result = 1. / (1. + exp(-1.0 * ( pix1 - 0.25) / pix2) );
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
    antscout << "total: " << result << " total-volume: " << result * volumeelement << std::endl;
    }
  else if( strcmp(operation.c_str(), "mean") == 0 )
    {
    antscout << result/ct << std::endl;
    }
  else
    {
    antscout << "operation " << operation << std::endl;
    }
  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(varimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int VImageMath(int argc, char *argv[])
{
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef FieldType                                                       ImageType;
  typedef typename ImageType::PixelType                                   PixelType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string operation = std::string(argv[argct]);  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct++]);
    }

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::Pointer varimage = NULL;

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
    antscout << " trans " << m_Transform0->GetParameters() << " Nspc " << image2->GetSpacing() << std::endl;
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
    volumeelement *= varimage->GetSpacing()[i];
    }

  PixelType result = itk::NumericTraits<PixelType>::Zero;
  Iterator  vfIter2( varimage,  varimage->GetLargestPossibleRegion() );
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
      ::ants::antscout << "Operation v/ not implemented" << std::endl;
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
    antscout << "total: " << result << " total-volume: " << result * volumeelement << std::endl;
    }
  else
    {
    antscout << "operation " << operation << std::endl;
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
  typedef float                                              ValueType;
  typedef itk::DiffusionTensor3D<ValueType>                  TensorType;
  typedef itk::Image<TensorType,ImageDimension>              TensorImageType;

  typedef itk::RecursiveGaussianImageFilter<TensorImageType, TensorImageType>
    GaussianFilterType;

  int         argct = 2;
  const std::string outname = std::string(argv[argct++]);
  std::string operation = std::string(argv[argct++]);  
  std::string fn1 = std::string(argv[argct++]);
  float sigma = atof(argv[argct++]);
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct++]);
    }

  typename TensorImageType::Pointer inDT = NULL;
  ReadTensorImage<TensorImageType>(inDT, fn1.c_str());

  typename GaussianFilterType::Pointer gFilter = GaussianFilterType::New();
  gFilter->SetInput( inDT );
  gFilter->SetSigma( sigma );
  gFilter->Update();

  WriteTensorImage<TensorImageType>(gFilter->GetOutput(), outname.c_str());
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
  typedef itk::ImageFileReader<ImageType>                    readertype;
  typedef itk::ImageFileWriter<ImageType>                    writertype;
  typedef itk::ImageFileWriter<ColorImageType>               ColorWriterType;
  typedef itk::ImageRegionIteratorWithIndex<TensorImageType> Iterator;

  typedef itk::Vector<float, ImageDimension>                 VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> vecIterator;
  typedef itk::Matrix<float, 3, 3>                           MatrixType;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string operation = std::string(argv[argct]);  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = ""; // used for whichvec and mask file name below
  
  typename TensorImageType::Pointer timage = NULL;    // input tensor image
  typename ImageType::Pointer       vimage = NULL;    // output scalar image
  typename ColorImageType::Pointer  cimage = NULL;    // output color image
  typename VectorImageType::Pointer vecimage = NULL;  // output vector image
  typename TensorImageType::Pointer toimage = NULL;   // output tensor image
  typename ImageType::Pointer       mimage = NULL;    // mask image

  if( strcmp(operation.c_str(), "4DTensorTo3DTensor") == 0 )
    {
    antscout
      <<
      " Convert a 4D tensor to a 3D tensor --- if there are 7 components to the tensor, we throw away the first component b/c its probably b0 "
      << std::endl;
    itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn1.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn1.c_str() );
    imageIO->ReadImageInformation();
    unsigned int dim = imageIO->GetNumberOfDimensions();
    if( dim == 4 )
      {
      typename D4TensorImageType::Pointer d4img = NULL;
      ReadImage<D4TensorImageType>(d4img, fn1.c_str() );
      unsigned int d4size = d4img->GetLargestPossibleRegion().GetSize()[3];
      if( d4size != 6 && d4size != 7 )
        {
        antscout << " you should not be using this function if the input data is not a tensor. " << std::endl;
        antscout
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
    antscout << " cannot convert --- input image not 4D --- " << fn1 << std::endl;
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
        antscout << "Unrecognized component.  Need to specify "
                 << "xx, xy, xz, yy, yz, or zz";
        return EXIT_FAILURE;
        }
      }
    else
      {
      antscout << "Error:  need to specify component (xx, xy, xz, yy, yz, zz)";
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
    fn2 = std::string(argv[argct]);
    whichvec = atoi(fn2.c_str());
    argct++;
    }

  ReadTensorImage<TensorImageType>(timage, fn1.c_str(), false);

  if( strcmp(operation.c_str(), "TensorIOTest") == 0 )
    {
    antscout << " test function for tensor I/O " << std::endl;
    WriteTensorImage<TensorImageType>(timage, outname.c_str(), false);
    return 0;
    }
  antscout << " imagedir " << timage->GetDirection() << std::endl;

  if( strcmp(operation.c_str(), "TensorColor") == 0 )
    {
    cimage =
      AllocImage<ColorImageType>(timage);
    

    if( argc > 5 )
      {
      antscout << "Using mask image: " << fn2 << std::endl;
      ReadImage<ImageType>(mimage, fn2.c_str() );
      }

    }
  else if( strcmp(operation.c_str(), "TensorMask") == 0 )
    {
    antscout << "Using mask image: " << fn2 << std::endl;

    ReadImage<ImageType>(mimage, fn2.c_str());

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

  
  TensorType zeroTensor;  // for masking background tensors 

  for ( unsigned int i = 0; i < 6; i++ ) 
    {
    zeroTensor[i] = 0.0;
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
      result = GetTensorFA<TensorType>(tIter.Value() );
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorMeanDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 0);
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorRadialDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 2);
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorEigenvalue") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 3 + whichvec);
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorAxialDiffusion") == 0 )
      {
      result = GetTensorADC<TensorType>(tIter.Value(), 5);
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorFANumerator") == 0 )
      {
      result = GetTensorFANumerator<TensorType>(tIter.Value() );
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorFADenominator") == 0 )
      {
      result = GetTensorFADenominator<TensorType>(tIter.Value() );
      if( vnl_math_isnan(result) )
        {
        result = 0;
        }
      vimage->SetPixel(ind, result);
      }
    else if( strcmp(operation.c_str(), "TensorColor") == 0 )
      {
      if ( argc > 5 )
        {
        if ( mimage->GetPixel( tIter.GetIndex() ) > 0 )
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
    
      if (maskVal > 0.0) 
        {
          toimage->SetPixel( ind, tIter.Value() );
        }
      else 
        {
          toimage->SetPixel( ind, zeroTensor );
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
        vimage->SetPixel(ind, tIter.Value()[whichvec]);
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
      bool hasNeg = false;
      for( unsigned int i = 0; i < 3; i++ )
        {
        if ( eigenValues[i] < 0 )
          {
          hasNeg = true;
          }
        eigenValuesMatrix(i, i) = fabs( eigenValues[i] );
        }
      if ( hasNeg )
        {
        std::cout << tIter.Value() << std::endl;
        std::cout << eigenValues << std::endl;
        }

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
    antscout << "Writing scalar image" << std::endl;
    WriteImage<ImageType>(vimage, outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int CompareHeadersAndImages(int argc, char *argv[])
{
  typedef float PixelType;
  //  const unsigned int ImageDimension = AvantsImageDimension;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  bool isfloat = false;
  try
    {
    ReadImage<ImageType>( image2, fn2.c_str() );
    }
  catch( ... )
    {
    antscout << " Error reading " << fn2 << std::endl;
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
    antscout << " read 1 error ";
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
  antscout << " SpacingError: " << sqrt(sperr) << std::endl;

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
  antscout << " OriginError: " << sqrt(operr) << std::endl;
  antscout << " OriginSignError: " << orsignerr << std::endl;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      float temp = image1->GetDirection()[i][j] - image2->GetDirection()[i][j];
      merr += temp * temp;
      }
    }
  antscout << " OrientError: " << sqrt(merr) << std::endl;

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
          antscout << " zero image2 error ";
          fixed_center.Fill(0);
          }
        }
      catch( ... )
        {
        antscout << " zero image1 error ";
        }

      typedef itk::TranslationTransform<double, ImageDimension> TransformType0;
      typename TransformType0::Pointer             m_Transform0 = TransformType0::New();
      typename TransformType0::ParametersType trans = m_Transform0->GetParameters();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        trans[i] = moving_center[i] - fixed_center[i];
        }
      m_Transform0->SetParameters(trans);
      antscout << " trans " << m_Transform0->GetParameters() << std::endl;
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
      if( vnl_math_isnan(pix2) || vnl_math_isinf(pix2) )
        {
        pix2 = 0; image2->SetPixel(ind, 0);
        }
      //      float pix3=varimage->GetPixel(ind);
      float pix1 = image1->GetPixel(ind);
      if( pix1 > 0  )
        {
        i1norm += fabs(pix1); ct1++;
        }
      if( pix2 > 0  )
        {
        i2norm += fabs(pix2); ct2++;
        }
      // if (pix3 > 0  )  {  i3norm+=fabs(pix3); ct3++; }
      }
    float mean1 = i1norm / ct1;
    if( mean1 == 0 )
      {
      mean1 = 1;
      }
    float mean2 = i2norm / ct2;
    if( mean2 == 0 )
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
        i1i2norm += fabs(pix1 - pix2);
        ct12++;
        }
      //      if ( pix1 > 0 || pix3 > 0)
      // {
      //  i1i3norm+=fabs(pix1-pix3);
      //  ct13++;
      // }
      }
    float idice0 = 1.0 - 2.0 * i1i2norm / ct12 /  (  i1norm / ct1 + i2norm / ct2 );
    //  float idice1 = 1.0 - 2.0*i1i3norm/ct13 /  (  i1norm/ct1 + i3norm / ct3 );
    antscout << " DiceImageDifference: " << idice0 << " IntensityDifference: " << i1i2norm << std::endl;
    // antscout << " CenterOfMassTransImageDifference: " << idice1 << " and " << i1i3norm << std::endl;
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
      else if( (op1[i] > 0 && op2[i] == 0 ) || (op1[i] < 0 && op2[i] == 0) )
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
  antscout << "  FailureState: " << failure << " for " << fn2  << std::endl;
  return failure;
}

// template<class TImage>
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
//     antscout << " Initial Means " << initialMeans[i] << " ";
//     }
//   antscout << std::endl;
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

// template<class TImage>
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
//     antscout << " Sample SD Ests " << sqrt(estimatedVar[i]) << " Mean " << estimatedMeans[i] <<  std::endl;
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
//   typename InputImageType::Pointer priors=NULL;
//
//   if (priorfn.length()  > 3 )
//     {
//     antscout << " Setting Priors " << priorfn << std::endl;
//     bool geometric=false;
//     if ( strcmp(priorfn.c_str(),"Geometric") == 0) geometric=true;
//     if (geometric)
//       {
//       antscout <<" Using a geometric thickness prior to aid cortical segmentation " << std::endl;
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
//       antscout <<" Allocated " << std::endl;
//
//       for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
//     {
// //    antscout <<" ind " <<vfIter2.GetIndex() << std::endl;
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
// //    antscout << " Sulc " << sulcprob << std::endl;
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
//       antscout << " ok " << std::endl;
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
//     } else antscout << " No Priors " << std::endl;
//
//   filter->SetInput(  vecImage );
//
//
//   if( nsmooth >= 1  )
//     {
//     antscout << " Smoothing Iterations:  " << nsmooth << std::endl;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);

  typename ImageType::Pointer image1 = NULL;
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
  if( mx == mn )
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

// template<class TImage>
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
//   applyEstimateModel->Print(antscout);
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
//   antscout << " mean dist " << meanDistance << std::endl;
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
//      //   antscout << " ow " << weights[jj] << " nw " <<  testNewNeighborhoodWeight[jj] << std::endl;
//     jj++;
//     }
//
// //  applyMRFImageFilter->SetMRFNeighborhoodWeight( testNewNeighborhoodWeight );
//   //  applyMRFImageFilter->SetMRFNeighborhoodWeight( weights );
//
//   //Kick off the MRF labeller function
//   applyMRFImageFilter->Update();
//
//   applyMRFImageFilter->Print(antscout);
//   antscout << "Number of Iterations : " << applyMRFImageFilter->GetNumberOfIterations()
//     << std::endl;
//   antscout << "Stop condition: (1) Maximum number of iterations (2) Error tolerance:  "
//     << applyMRFImageFilter->GetStopCondition() << std::endl;
//
//   typename ClassImageType::Pointer  outClassImage = applyMRFImageFilter->GetOutput();
//
//   //Testing of different parameter access functions in the filter
//   antscout << "The number of classes labelled was: " <<
//     applyMRFImageFilter->GetNumberOfClasses() << std::endl;
//   antscout << "The maximum number of iterations were: " <<
//     applyMRFImageFilter->GetMaximumNumberOfIterations() << std::endl;
//   antscout << "The error tolerace threshold was: " <<
//     applyMRFImageFilter->GetErrorTolerance() << std::endl;
//   antscout << "The smoothing MRF parameter used was: " <<
//     applyMRFImageFilter->GetSmoothingFactor() << std::endl;
//   antscout << "The MRF neighborhood weights are: " << std::endl;
//
//
//   return  outClassImage;
//
// }

template <class TImage>
typename TImage::Pointer
// void
itkMRIBiasFieldCorrectionFilter(typename TImage::Pointer image,
                                typename TImage::Pointer labelimage, unsigned int sd = 2)
{
  antscout << "doing Bias corr " << std::endl;
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
      float n = classCounts[label];
      classMeans[label] = (n - 1.0) / ( n ) *  classMeans[label] + 1.0 / n * pix;
      double sqrdf =  ( pix -  classMeans[label]) * ( pix -  classMeans[label]);
      if( n > 1 )
        {
        classSigmas[label] =  (n - 1) / n * classSigmas[label] +  1.0 / (n - 1) * sqrdf;
        }
      }

    ++o_iter;
    }
  for( unsigned int k = 0; k < numclasses; k++ )
    {
    classSigmas[k] = sqrt(classSigmas[k]);
    antscout << " Initial Means pre-bias " << classMeans[k] << " sig " <<  classSigmas[k] << std::endl;
    }

  // creats a normal random variate generator
  // itk::Statistics::NormalVariateGenerator::Pointer randomGenerator =
  //  itk::Statistics::NormalVariateGenerator::New() ;

  // creates a bias correction filter and run it.
  typedef itk::MRIBiasFieldCorrectionFilter<ImageType, ImageType, ImageType> FilterType;

  antscout << "before new filter" << std::endl;
  typename FilterType::Pointer filter = FilterType::New();
  antscout << "after new filter" << std::endl;

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
  long int t1 = time(NULL);
  filter->Update();
  long int t2 = time(NULL);
  antscout << "Run time (in s)" << t2 - t1  << std::endl;

  return filter->GetOutput();
}

// template<class TImage>
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
//   antscout << "Starting to build the K-means model ....." << std::endl;
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
//   antscout << "Result of K-Means clustering" << std::endl;
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
//   antscout << " mean dist " << meanDistance << std::endl;
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
//     antscout <<  (membershipFunctions[classIndex]->GetCentroid())[0] << std::endl;
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
//      //antscout << " ow " << weights[jj] << " nw " <<  testNewNeighborhoodWeight[jj] << std::endl;
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
//   applyMRFFilter->Print(antscout);
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sigma = 1.0;
  if( argc > argct )
    {
    sigma = atof(argv[argct]);
    }

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer varimage = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(sigma * sigma);
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]); argct++;
  std::string operation = std::string(argv[argct]);  argct++;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       sigma = 1.0;
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

  typename ImageType::Pointer image1 = NULL;
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
  typedef float                                         InternalPixelType;
  typedef itk::Image<InternalPixelType, ImageDimension> InternalImageType;
  typedef itk::Image<InternalPixelType, ImageDimension> ImageType;

  typedef unsigned char                               OutputPixelType;
  typedef itk::Image<OutputPixelType, ImageDimension> OutputImageType;

  //  option->SetUsageOption( 0, "[speedImage,seedImage,<stoppingValue=max>,<topologyCheck=0>]" );
  unsigned int argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  std::string  fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }
  else
    {
    antscout << " not enough parameters -- need label image " << std::endl;  return 0;
    }
  float stoppingValue = 100.0;
  if(  argc > argct )
    {
    stoppingValue = atof(argv[argct]);   argct++;
    }
  int topocheck = 0;
  if(  argc > argct )
    {
    topocheck = atoi(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image1;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typedef itk::FastMarchingImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image1 );

  typedef typename FilterType::NodeContainer  NodeContainer;
  typedef typename FilterType::NodeType       NodeType;
  typedef typename FilterType::LabelImageType LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
  typename LabelImageReaderType::Pointer labelImageReader =
    LabelImageReaderType::New();
  labelImageReader->SetFileName( fn2.c_str() );
  labelImageReader->Update();

  typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType>
    ContourFilterType;
  typename ContourFilterType::Pointer contour = ContourFilterType::New();
  contour->SetInput( labelImageReader->GetOutput() );
  contour->FullyConnectedOff();
  contour->SetBackgroundValue(
    itk::NumericTraits<typename LabelImageType::PixelType>::Zero );
  contour->Update();

  typename NodeContainer::Pointer alivePoints = NodeContainer::New();
  alivePoints->Initialize();
  unsigned long aliveCount = 0;
  typename NodeContainer::Pointer trialPoints = NodeContainer::New();
  trialPoints->Initialize();
  unsigned long trialCount = 0;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL(
    labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItC( contour->GetOutput(),
                                                         contour->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(), ItC.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItC )
    {
    if( ItC.Get() !=
        itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItC.GetIndex();

      NodeType     node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      trialPoints->InsertElement( trialCount++, node );
      }
    else if( ItL.Get() !=
             itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItL.GetIndex();

      NodeType     node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      alivePoints->InsertElement( aliveCount++, node );
      }
    }
  filter->SetTrialPoints( trialPoints );
  filter->SetAlivePoints( alivePoints );

  filter->SetStoppingValue( stoppingValue );
  filter->SetTopologyCheck( FilterType::None );
  if( topocheck == 1 )  // Strict
    {
    antscout << " strict " << std::endl;
    filter->SetTopologyCheck( FilterType::Strict );
    }
  if( topocheck == 2 )  // No handles
    {
    antscout << " no handles " << std::endl;
    filter->SetTopologyCheck( FilterType::NoHandles );
    }

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    antscout << "Exception caught !" << std::endl;
    antscout << excep << std::endl;
    }

  itk::ImageRegionIteratorWithIndex<ImageType> ItF(
    filter->GetOutput(),
    filter->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(), ItF.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItF )
    {
    if( ItL.Get() !=
        itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      ItF.Set( -ItF.Get() );
      }
    }

  if( outname.length() > 3 )
    {
    WriteImage<ImageType>(filter->GetOutput(), outname.c_str() );
    }

  return 0;
}

template <unsigned int ImageDimension>
int PropagateLabelsThroughMask(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  float       thresh = 0.5;
  int         argct = 2;
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
    antscout << " not enough parameters -- need label image " << std::endl;  return 0;
    }
  if(  argc > argct )
    {
    stopval = atof(argv[argct]);   argct++;
    }

  typename ImageType::Pointer speedimage = NULL;
  ReadImage<ImageType>(speedimage, fn1.c_str() );
  typename ImageType::Pointer labimage = NULL;
  ReadImage<ImageType>(labimage, fn2.c_str() );
  typename ImageType::Pointer fastimage = NULL;
  ReadImage<ImageType>(fastimage, fn1.c_str() );
  typename ImageType::Pointer outlabimage = NULL;
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
    if( speedval < thresh )
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
  for( unsigned int lab = 1; lab <= (unsigned int)maxlabel; lab++ )
    {
    typedef  itk::FastMarchingImageFilter<ImageType, ImageType> FastMarchingFilterType;
    typename FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
    fastMarching->SetInput( speedimage );
    fastMarching->SetTopologyCheck( FastMarchingFilterType::None );

    typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
    typedef typename FastMarchingFilterType::NodeType      NodeType;
    typename NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();
    unsigned long ct = 0;
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      bool   isinside = true;
      double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
      double labval = labimage->GetPixel(vfIter2.GetIndex() );
      if( speedval < thresh )
        {
        isinside = false;
        }
      if( isinside && (unsigned int) labval == lab )
        {
        NodeType     node;
        const double seedValue = 0.0;
        node.SetValue( seedValue );
        node.SetIndex( vfIter2.GetIndex() );
        seeds->InsertElement( ct, node );
        ct++;
        }
      }
    fastMarching->SetTrialPoints(  seeds  );
    fastMarching->SetStoppingValue(  stopval );
    fastMarching->Update();
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      bool   isinside = true;
      double speedval = speedimage->GetPixel(vfIter2.GetIndex() );
      double labval = labimage->GetPixel(vfIter2.GetIndex() );
      if( speedval < thresh )
        {
        isinside = false;
        }
      if( isinside && labval == 0 )
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
  WriteImage<ImageType>(fastimage, kname.c_str() );
  WriteImage<ImageType>(outlabimage, outname.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int DistanceMap(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    antscout << "Missing required arguments ( output name, operation & fn1)" << std::endl;
    throw;
    }
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer image1 = NULL;
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
int FillHoles(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<int, ImageDimension>                                 LabelImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  typedef itk::CastImageFilter<ImageType, LabelImageType>                 CastFilterType;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       holeparam = 2.0;
  if( argc > argct )
    {
    holeparam = atof(argv[argct]);
    }

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::Pointer imageout = NULL;
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
  catch( itk::ExceptionObject & excep )
    {
    antscout << "Relabel: exception caught !" << std::endl;
    antscout << excep << std::endl;
    }

  // WriteImage<ImageType>(relabel->GetOutput(),"test.nii");

  if( holeparam == 2 )
    {
    antscout << " Filling all holes " <<  std::endl;
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
        if( p == lab )
          {
          volume++;
          for( unsigned int i = 0; i < GHood.Size(); i++ )
            {
            ind2 = GHood.GetIndex(i);
            float val2 = image1->GetPixel(ind2);
            if( val2 >= 0.5 && GHood.GetPixel(i) != lab )
              {
              objectedge++;
              totaledge++;
              }
            else if( val2 < 0.5 && GHood.GetPixel(i) != lab )
              {
              backgroundedge++;
              totaledge++;
              }
            }
          }
        ++GHood;
        }

      float vrat = (float)totaledge / (float)volume;
      erat = (float)objectedge / (float)totaledge;
      antscout << " Lab " << lab << " volume " << volume << " v-rat " << vrat << " edge " << erat << std::endl;
      }

    if( erat > holeparam ) // fill the hole
      {
      antscout << " Filling " << lab << " of " << maximum <<  std::endl;
      typedef itk::ImageRegionIteratorWithIndex<ImageType> RelabelIterator;
      RelabelIterator vfIter( relabel->GetOutput(),
                              relabel->GetOutput()->GetLargestPossibleRegion() );
      for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
        {
        if( vfIter.Get() == lab )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  float       option = 0;
  if( argc > argct )
    {
    option = atof(argv[argct]);
    }

  typename ImageType::Pointer image = NULL;
  ReadImage<ImageType>(image, fn1.c_str() );

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
    if( option == 0 )
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

  return 0;
}

template <unsigned int ImageDimension>
int PrintHeader(int argc, char *argv[])
{
  typedef  float                                     outPixelType;
  typedef  float                                     floatPixelType;
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::Image<floatPixelType, ImageDimension> IntermediateType;
  typedef itk::Image<outPixelType, ImageDimension>   OutImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;
  typedef itk::ImageFileWriter<OutImageType>         writertype;

  int         argct = 4;
  std::string fn1 = std::string(argv[argct]);
  if( argc > 20 )
    {
    antscout << " k " << std::endl;
    }

  typename ImageType::Pointer image;
  ReadImage<ImageType>(image, fn1.c_str() );
  antscout << " Spacing " << image->GetSpacing() << std::endl;
  antscout << " Origin " << image->GetOrigin() << std::endl;
  antscout << " Direction " << std::endl << image->GetDirection() << std::endl;
  antscout << " Size " << std::endl << image->GetLargestPossibleRegion().GetSize() << std::endl;

  //  if (strcmp(operation.c_str(),"n_last_dim") == 0){
  // unsigned int lastdim=image->GetLargestPossibleRegion().GetSize()[ImageDimension-1];
  //   std::ofstream logfile;
  // logfile.open(outname.c_str() );
  // if (logfile.good()  )
  // {
  //  logfile << lastdim << std::endl;
  // }
  // cd antscout << lastdim << std::endl;
  // }
  return 1;
}

template <unsigned int ImageDimension>
int GradientImage(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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
    normalize = atoi(argv[argct]);
    }

  typename ImageType::Pointer image = NULL;
  typename ImageType::Pointer image2 = NULL;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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
    normalize = atoi(argv[argct]);
    }

  typename ImageType::Pointer image = NULL;
  typename ImageType::Pointer image2 = NULL;
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
int PoissonDiffusion( int argc, char *argv[])
{
  if( argc < 6 )
    {
    antscout << "Usage error---not enough arguments.   See help menu."
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

  typename ImageType::Pointer output = duplicator->GetModifiableOutput();
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
    maximumNumberOfIterations = atoi( argv[8] );
    }
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  float        lastmean = 0;
  unsigned int iterations = 0;
  while( iterations++ < maximumNumberOfIterations && convergence >= convergenceThreshold )
    {
    antscout << "  Iteration " << iterations << ": " << convergence << std::endl;
    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( vnl_math_sqr( sigma ) );
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

    convergence = stats->GetMean( 1 ) - lastmean;
    lastmean = stats->GetMean( 1 );
    output =  smoother->GetOutput();
    output->DisconnectPipeline();

    Iterator vfIter( output,  output->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      if(  thresholder2->GetOutput()->GetPixel(vfIter.GetIndex() ) == 1 )
        {
        vfIter.Set(reader->GetOutput()->GetPixel(vfIter.GetIndex() ) );
        }
      }
    }

  Iterator vfIter( output,  output->GetLargestPossibleRegion() );
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    if(  thresholder2->GetOutput()->GetPixel(vfIter.GetIndex() ) == 1 )
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
    antscout << " too few options " << std::endl; return;
    }
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;

  typename ImageType::Pointer input = NULL;
  ReadImage<ImageType>(input, fn1.c_str() );

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

//  antscout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename ImageType::PixelType p = GHood.GetCenterPixel();
    typename ImageType::IndexType ind = GHood.GetIndex();
    typename ImageType::IndexType ind2;
    if( p > 0  )
      {
      bool atedge = false;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        ind2 = GHood.GetIndex(i);
        if( GHood.GetPixel(i) > 0 && GHood.GetPixel(i)  != p )
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
EnumerateLabelInterfaces(int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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

  typename ImageType::Pointer input = NULL;
  ReadImage<ImageType>(input, fn1.c_str() );
  typename ImageType::Pointer output = NULL;
  ReadImage<ImageType>(output, fn1.c_str() );
  output->FillBuffer(0);
  typename ImageType::Pointer colored = NULL;
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

  antscout << " Max Label " << max << std::endl;
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

  typedef itk::ImageRegionIteratorWithIndex<myInterfaceImageType> FIterator;

// we can use this to compute a 4-coloring of the brain

  typedef itk::NeighborhoodIterator<ImageType> iteratorType;
  typename iteratorType::RadiusType rad;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    rad[j] = 1;
    }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion() );

  GHood.GoToBegin();

//  antscout << " foreg " << (int) foreground;
  while( !GHood.IsAtEnd() )
    {
    typename ImageType::PixelType p = GHood.GetCenterPixel();
    typename ImageType::IndexType ind = GHood.GetIndex();
    typename ImageType::IndexType ind2;
    if( p > 0  )
      {
      bool atedge = false;

      unsigned long linearinda = 0, linearindb = 0, linearind = 0;
      for( unsigned int i = 0; i < GHood.Size(); i++ )
        {
        ind2 = GHood.GetIndex(i);
        if( GHood.GetPixel(i) > 0 && GHood.GetPixel(i)  != p )
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
      // antscout << i <<"  :: " << faceimage->GetPixel(find)  << std::endl;
      }
    antscout << " total interfaces for label :  " << j << " are " << total << std::endl;
    for( unsigned int i = 0; i <= max; i++ )
      {
      find[1] = i;
      if( total > 0 )
        {
        faceimage->SetPixel(find, faceimage->GetPixel(find) / total);
        }
      if( faceimage->GetPixel(find) >  0.01 )
        {
        antscout << i << "  :: " << faceimage->GetPixel(find)  << std::endl;
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
    antscout << " Label " << j << " color " << okcolor << std::endl;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ByteImageType>                             writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::Pointer mask = NULL;
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
      err += fabs(locerr);
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
    logfile << " Err " <<  " : " << err << " %ER " <<  " : " << err / (maskct) * 100. << " NER " <<  " : " << negerr
      / maskct * 100.0 << std::endl;
    }
  return 0;
}

template <unsigned int ImageDimension>
int DiceAndMinDistSum(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ByteImageType>                             writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
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

  typename ImageType::Pointer image1 = NULL;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::Pointer outdist = NULL;
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
    typename ImageType::Pointer surf = NULL;
    typename ImageType::Pointer d1 = NULL;
    typename ImageType::Pointer d2 = NULL;

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
      if( It2.Get() == *it || image1->GetPixel(It2.GetIndex() ) == *it )
        {
        countu++;
        }
      if( It2.Get() == *it && image1->GetPixel(It2.GetIndex() ) == *it )
        {
        counti++;
        }

      if( It2.Get() == *it )
        {
        count2++;
        if( d1 )
          {
          dist2 += d1->GetPixel(It2.GetIndex() );
          }
        }
      if( image1->GetPixel(It2.GetIndex() ) == *it )
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
    //    antscout << " sdist " << surfdist << " sct " << surfct << std::endl

    if( outdist )
      {
      WriteImage<ImageType>(outdist, outdistfn.c_str() );
      }

    if( count2 + count1 > 0 )
      {
      distances[labct] += (dist2 + dist1) / (count2 + count1);
      dicevals[labct] = 2.0 * counti / (count2 + count1);
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

    ColumnHeaders.push_back("Label Name");
    ColumnHeaders.push_back("Min_Distance");
    ColumnHeaders.push_back("Dice");
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
    catch( itk::ExceptionObject& exp )
      {
      antscout << "Exception caught!" << std::endl;
      antscout << exp << std::endl;
      }
    }
  else
    {
    vnl_matrix<double> OutputValues(NumberOfLabels, 4);
    ColumnHeaders.push_back("Label Name");
    ColumnHeaders.push_back("Dice");
    ColumnHeaders.push_back("RO");
    ColumnHeaders.push_back("Percent_of_Region_1_In_Overlap");
    ColumnHeaders.push_back("Percent_of_Region_2_In_Overlap");
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
      antscout << "Output written to " << outname.c_str() << ".csv." << std::endl;
      }
    catch( itk::ExceptionObject& exp )
      {
      antscout << "Exception caught!" << std::endl;
      antscout << exp << std::endl;
      }
    }

  labelcount = 0;
  labct = 0;
  for( it = myLabelSet2.begin(); it != myLabelSet2.end(); ++it )

    {
    sum += distances[labct];
    sumdice += dicevals[labct];
    sumro += rovals[labct];
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
    antscout << " k " << std::endl;
    }
  antscout << " Compute Lipschitz continuity of the mapping " << std::endl;

  typedef float                                              RealType;
  typedef itk::Image<RealType, ImageDimension>               RealImageType;
  typedef itk::Vector<RealType, ImageDimension>              VectorType;
  typedef itk::Image<VectorType, ImageDimension>             VectorImageType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> Iterator;

  int         argct = 2;
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
      xt[i] = x[i] + vecx[i];
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

        float numer = 0;
        float denom = 0;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          yt[i] = y[i] + vecy[i];
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
//      if (ct1 % 1000 == 0) antscout << " Progress : " << (float ) ct1 / (float) numpx *100.0 << " val " <<
// localmaxval << std::endl;
    }

  antscout << " Lipschitz continuity related to: " << globalmaxval << std::endl;
  antscout << " Tx :  " << gxt << "  Ty: " << gyt << std::endl;
  antscout << " x :  " << gx << "  y: " << gy << std::endl;
  timer.Stop();
//    antscout << "Elapsed time: " << timer.GetMeanTime()  << std::endl;

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
    antscout << " too few options " << std::endl;
    return 1;
    }
  typedef float                                       PixelType;
  typedef itk::VectorImage<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PixelType, ImageDimension>       RealImageType;
  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  inname = std::string(argv[argct]);   argct++;
  unsigned int whichvec = atoi(argv[argct]);   argct++;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( inname.c_str() );
  reader1->Update();
  typename ImageType::Pointer vecimage = reader1->GetOutput();
  if( whichvec >= vecimage->GetVectorLength() )
    {
    antscout << " input image " << inname << " only has " << vecimage->GetVectorLength() << " components "
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
    antscout << " Compute  phi(  phi^{-1}(x)) " << std::endl;
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

  int         argct = 2;
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
  float    globalmaxval = 0;
  Iterator It1( vecimage1, vecimage1->GetLargestPossibleRegion() );
  for( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
    {
    typename VectorImageType::PointType x;
    typename DefaultInterpolatorType::PointType xt, yt;
    vecimage1->TransformIndexToPhysicalPoint( It1.GetIndex(), x);
    typename VectorImageType::PixelType vecx = It1.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      xt[i] = x[i] + vecx[i];
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

    float error = 0;
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
      {
      error += (yt[jj] - x[jj]) * (yt[jj] - x[jj]);
      }
    error = sqrt(error);
    if( error >  globalmaxval )
      {
      globalmaxval = error; gx = x;
      }
    invid->SetPixel(It1.GetIndex(), error);
    }
  antscout << " Max error " << globalmaxval << " at " << gx << std::endl;
  timer.Stop();
//    antscout << "Elapsed time: " << timer.GetMeanTime()  << std::endl;

  WriteImage<RealImageType>( invid, outname.c_str() );

  return 0;
}

template <unsigned int ImageDimension>
int LabelStats(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ByteImageType>                             writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  std::string imagename = ANTSGetFilePrefix(outname.c_str() ) + std::string("_square.nii.gz");
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if(  argc > argct )
    {
    fn2 = std::string(argv[argct]);   argct++;
    }

  typename ImageType::Pointer image = NULL;
  typename ImageType::Pointer valimage = NULL;
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
    volumeelement *= spacing[i];
    }

  std::ofstream logfile;
  logfile.open(outname.c_str() );
  logfile << "x,y,z,t,label" << std::endl;

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
      if(  label == currentlabel  )
        {
        totalvolume += volumeelement;
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
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= (float)totalct;
      }

    if( !valimage )
      {
      antscout << " Volume Of Label " << *it << " is " << totalvolume <<   "  Avg-Location " << myCenterOfMass
               << std::endl;
      }
    else // if ( totalvolume > 500 &&  totalmass/totalct > 1/500 )  {
      {
      antscout << " Volume Of Label " << *it << " is " << totalvolume <<   "  Avg-Location " << myCenterOfMass
               << " mass is " << totalmass << " average-val is " << totalmass / totalct << std::endl;
      //      antscout << *it << "  " <<  totalvolume <<  " & " <<  totalmass/totalct   << " \ " << std::endl;
      }

// square image
    squareimage->GetBufferPointer()[labelcount] = totalmass / totalct;
    if ( ImageDimension == 2 ) 
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << ",0,0," << currentlabel << std::endl;
    if ( ImageDimension == 3 ) 
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << "," << myCenterOfMass[2] << ",0," << currentlabel << std::endl;
    if ( ImageDimension == 4 ) 
      logfile << myCenterOfMass[0] << "," << myCenterOfMass[1] << "," << myCenterOfMass[2] << "," << myCenterOfMass[3] << "," << currentlabel << std::endl;
    labelcount++;
    }

  logfile.close();

  WriteImage<TwoDImageType>(squareimage, imagename.c_str() );

  return 0;
}

// int is the key, string the return value
std::map<unsigned int, std::string> RoiList(std::string file)
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

// now words can be accessed like this WordList[n]; where 'n' is the index

template <unsigned int ImageDimension>
int ROIStatistics(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ByteImageType>                             writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  //  if(grade_list.find("Tim") == grade_list.end()) {  antscout<<"Tim is not in the map!"<<endl; }
  // mymap.find('a')->second
  int argct = 2;
  if( argc < 6 )
    {
    antscout << " not enough parameters --- usage example 1 :" << "" << std::endl;
    antscout << argv[0]
             << " ImageMath  3 output.csv ROIStatistics roinames.txt LabelImage.nii.gz ValueImage.nii.gz  "
             << std::endl;
    throw std::exception();
    }
  const std::string outname = std::string(argv[argct]);
  typedef vnl_matrix<double> MatrixType;
  argct += 2;
  std::string fn0 = std::string(argv[argct]);   argct++;
  antscout << "  fn0 " << fn0 << std::endl;
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

  typename ImageType::Pointer image = NULL;
  typename ImageType::Pointer valimage = NULL;
  typename ImageType::Pointer valimage3 = NULL;
  typename ImageType::Pointer valimage4 = NULL;
  typename ImageType::Pointer valimage5 = NULL;
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
    volumeelement *= spacing[i];
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
      clusters[(unsigned long) label] = clusters[(unsigned long) label] + 1;
      masses[(unsigned long) label] = masses[(unsigned long) label] + vv;
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
    antscout << " it " << *it << " roi " << roi << " mylabel " << mylabel << std::endl;
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
        if(  label == mylabel )
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
      myCenterOfMass[i] /= (float)totalct;
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
              << masses[roi + 1] << "," << masses[roi + 1]/clusters[roi+1] << "," << comx << "," << comy << "," << comz << "," << comt << std::endl;
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
      if(  label == currentlabel  )
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
    if( totalct == 0 )
      {
      totalct = 1;
      }
    for( unsigned int i = 0; i < spacing.Size(); i++ )
      {
      myCenterOfMass[i] /= (float)totalct;
      }

    clusters[(unsigned long)*it] = totalvolume;
    masses[(unsigned long)*it] = totalmass / totalct;
    pvals[(unsigned long)*it] = 1.0 - maxoneminuspval;
    pvals3[(unsigned long)*it] = maxoneminuspval3;
    pvals4[(unsigned long)*it] = 1.0 - maxoneminuspval4;
    pvals5[(unsigned long)*it] = 1.0 - maxoneminuspval5;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ByteImageType>                             writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);

  typename ImageType::Pointer image = NULL;
  typename ByteImageType::Pointer image2 = NULL;
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int dof = 1;
  if( argc > argct )
    {
    dof = atoi(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image = NULL;
  ReadImage<ImageType>(image, fn1.c_str() );

  antscout << " read Image" << fn1 << " dof " << dof << std::endl;
  typedef itk::Statistics::TDistribution DistributionType;
  typename DistributionType::Pointer distributionFunction = DistributionType::New();
  distributionFunction->SetDegreesOfFreedom(  dof );
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator vfIter( image,  image->GetLargestPossibleRegion() );
  for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
    float val = image->GetPixel(vfIter.GetIndex() );
    if( !vnl_math_isnan(val) && !vnl_math_isinf(val) )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string  outname = std::string(argv[argct]);
  std::string  ext = itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int rowcoloption = atoi(argv[argct]);   argct++;
  std::string  maskfn = std::string(argv[argct]); argct++;
  unsigned int numberofimages = 0;
  typename ImageType::Pointer mask = NULL;
  ReadImage<ImageType>(mask, maskfn.c_str() );
  unsigned long voxct = 0;
  Iterator      mIter( mask, mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    {
    if( mIter.Get() >= 0.5 )
      {
      voxct++;
      }
    }

  typename ImageType::Pointer image2 = NULL;
  typename ImageType::SizeType size;
  size.Fill(0);
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);
        antscout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  antscout << " largest image " << size << " num images " << numberofimages << " voxct " << voxct << std::endl;
  unsigned long xx1 = 0, yy1 = 0;
  if( rowcoloption == 0 )
    {
    antscout << " row option " << std::endl;  xx1 = voxct;  yy1 = numberofimages;
    }
  if( rowcoloption == 1 )
    {
    antscout << " col option " << std::endl;  yy1 = voxct;  xx1 = numberofimages;
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
      antscout << " image " << j << " is "  << fn << std::endl;
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
        if( mIter.Get() >= 0.5 )
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
    writer->SetColumnHeaders( ColumnHeaders );
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& exp )
      {
      antscout << "Exception caught!" << std::endl;
      antscout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
/** declare the tiled image */
    typename MatrixImageType::SizeType tilesize;
    tilesize[0] = xsize;
    tilesize[1] = ysize;
    antscout << " allocate matrix " << tilesize << std::endl;
    typename MatrixImageType::RegionType region;
    region.SetSize( tilesize );

    typename MatrixImageType::Pointer matimage =
      AllocImage<MatrixImageType>(region);

    unsigned int imagecount = 0;
    for( unsigned int j = argct; j < argc; j++ )
      {
      std::string fn = std::string(argv[j]);
      ReadImage<ImageType>(image2, fn.c_str() );
      antscout << " image " << j << " is "  << fn << std::endl;
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
        if( mIter.Get() >= 0.5 )
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
          //          antscout << " Mind " << mind << std::endl;
          matimage->SetPixel(mind, image2->GetPixel(mIter.GetIndex() ) );
          tvoxct++;
          }
        }
      imagecount++;
      }

    antscout << " mat size " << matimage->GetLargestPossibleRegion().GetSize() << std::endl;
    WriteImage<MatrixImageType>(matimage, outname.c_str() );
    }
  return 0;
}

template <unsigned int ImageDimension>
int RandomlySampleImageSetToCSV(unsigned int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRandomConstIteratorWithIndex<ImageType>               Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string  outname = std::string(argv[argct]);
  std::string  ext = itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int n_samples = atoi(argv[argct]);   argct++;
  /* std::string maskfn=std::string(argv[argct]); argct++;
  typename ImageType::Pointer mask = NULL;
  ReadImage<ImageType>(mask,maskfn.c_str());
  Iterator mIter( mask,mask->GetLargestPossibleRegion() );
  for(  mIter.GoToBegin(); !mIter.IsAtEnd(); ++mIter )
    if (mIter.Get() >= 0.5) voxct++;
  */
  unsigned int  numberofimages = 0;

  typename ImageType::Pointer image2 = NULL;
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
      antscout << " image " << j << " is "  << fn << std::endl;
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
    catch( itk::ExceptionObject& exp )
      {
      antscout << "Exception caught!" << std::endl;
      antscout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    antscout << " need a csv file as output type , you tried " << outname << std::endl;
    }
  return 0;
}

template <class T>
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<PixelType, 2>                                        MatrixImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef typename ImageType::IndexType                                   IndexType;
  typedef typename ImageType::SizeType                                    SizeType;
  typedef typename ImageType::SpacingType                                 SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;

  int argct = 2;
  if( argc < 5 )
    {
    antscout << " need more args -- see usage   " << std::endl;  throw std::exception();
    }
  const std::string  outname = std::string(argv[argct]);
  std::string  ext = std::string(".csv"); // itksys::SystemTools::GetFilenameExtension( outname );
  argct += 2;
  unsigned int n_evecs = atoi(argv[argct]);   argct++;
  unsigned int rowcoloption = 1;
  std::string  maskfn = std::string(argv[argct]); argct++;
  unsigned int numberofimages = 0;
  typename ImageType::Pointer mask = NULL;
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
    antscout << " Max value in mask is <= 0, aborting. " << maxval << std::endl;
    throw std::exception();
    }

  typedef itk::Array2D<double> MatrixType;
  typename ImageType::Pointer image2 = NULL;
  typename ImageType::SizeType size;
  size.Fill(0);
  for( unsigned int j = argct; j < argc; j++ )
    {
    numberofimages++;
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();
    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( imageIO->GetDimensions(i) > size[i] )
        {
        size[i] = imageIO->GetDimensions(i);

        antscout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  antscout << " largest image " << size << " num images " << numberofimages << std::endl;
  MatrixType avg_matrix(numberofimages, maxval);
  avg_matrix.Fill(0);
  for( unsigned long mv = 1; mv <= maxval; mv++ )
    {
    /** 2. count the voxels in this label */
    unsigned long voxct = 0;
    for( Iterator mIter( mask, mask->GetLargestPossibleRegion() );
         !mIter.IsAtEnd(); ++mIter )
      {
      if( mIter.Get() == mv )
        {
        voxct++;
        }
      }

    unsigned long xx1 = 0, yy1 = 0;
    if( rowcoloption == 0 )
      {
      antscout << " row option " << std::endl;  xx1 = voxct;  yy1 = numberofimages;
      }
    if( rowcoloption == 1 )
      {
      antscout << " col option " << std::endl;  yy1 = voxct;  xx1 = numberofimages;
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
        antscout << " image " << j << " is "  << fn << std::endl;
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
          if( mIter.Get() == mv )
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
      catch( itk::ExceptionObject& exp )
        {
        antscout << "Exception caught!" << std::endl;
        antscout << exp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      antscout << " can only write out csv files " << std::endl;
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
    catch( itk::ExceptionObject& exp )
      {
      antscout << "Exception caught!" << std::endl;
      antscout << exp << std::endl;
      return EXIT_FAILURE;
      }
    }

  return 0;
}

template <unsigned int ImageDimension>
int ConvertImageToFile(      int argc, char *argv[])
{
  typedef float                                                           PixelType;
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  int         argct = 2;
  const std::string outname = std::string(argv[argct]);
  argct += 2;
  std::string fn1 = std::string(argv[argct]);   argct++;
  std::string fn2 = "";
  if( argc > argct )
    {
    fn2 = std::string(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image = NULL;
  ReadImage<ImageType>(image, fn1.c_str() );
  typename ImageType::Pointer mask = NULL;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>(mask, fn2.c_str() );
    }

  antscout << " read Image" << fn1 << " mask? " << fn2 << std::endl;
  std::ofstream logfile;
  logfile.open(outname.c_str() );
  if( logfile.good() )
    {
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
    ImageIterator vfIter( image,  image->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      bool getval = true;
      if( mask->GetPixel(vfIter.GetIndex() ) < 0.5 )
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
  typedef itk::Vector<float, ImageDimension>                              VectorType;
  typedef itk::Image<VectorType, ImageDimension>                          FieldType;
  typedef itk::Image<PixelType, ImageDimension>                           ImageType;
  typedef itk::Image<unsigned char, ImageDimension>                       ByteImageType;
  typedef itk::ImageFileReader<ImageType>                                 readertype;
  typedef itk::ImageFileWriter<ImageType>                                 writertype;
  typedef  typename ImageType::IndexType                                  IndexType;
  typedef  typename ImageType::SizeType                                   SizeType;
  typedef  typename ImageType::SpacingType                                SpacingType;
  typedef itk::AffineTransform<double, ImageDimension>                    AffineTransformType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double>          InterpolatorType1;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType2;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                    Iterator;
  int         argct = 2;
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
    antscout << " Not enough inputs " << std::endl;
    return 1;
    }
  unsigned int radius = 2;
  if( argc > argct )
    {
    radius = atoi(argv[argct]);
    }
  argct++;

  typename ImageType::Pointer image1 = NULL;
  ReadImage<ImageType>(image1, fn1.c_str() );
  typename ImageType::Pointer imageout = NULL;
  ReadImage<ImageType>(imageout, fn1.c_str() );
  imageout->FillBuffer(0);
  typename ImageType::Pointer image2 = NULL;
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
        typename ImageType::IndexType ind2;
        // compute mean difference
        float diff = 0.0;
        for( unsigned int i = 0; i < GHood.Size(); i++ )
          {
          ind2 = GHood.GetIndex(i);
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
  typedef itk::ImageFileReader<ImageType>               ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>               ImageWriterType;
  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  IteratorType;

  if( argc < 5 )
    {
    antscout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );

  // Read input segmentations
  const unsigned long                               nImages = argc - 4;
  typename std::vector<typename ImageType::Pointer> images(argc - 4);
  for( int i = 4; i < argc; i++ )
    {
    ReadImage<ImageType>( images[i - 4], argv[i] );
    }

  // Find maximum label
  int maxLabel = 0;
  for( unsigned int i = 0; i < nImages; i++ )
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

  typename ImageType::Pointer output =
    AllocImage<ImageType>(images[0], 0);

  IteratorType              it( output, output->GetLargestPossibleRegion() );
  itk::Array<unsigned long> votes;
  votes.SetSize( nLabels );

  while( !it.IsAtEnd() )
    {
    votes.Fill(0);
    unsigned long maxVotes = 0;
    unsigned long votedLabel = 0;
    for( unsigned long i = 0; i < nImages; i++ )
      {
      unsigned long label = images[i]->GetPixel( it.GetIndex() );
      votes.SetElement(label, votes.GetElement(label) + 1 );

      if( votes.GetElement(label) > maxVotes )
        {
        maxVotes = votes.GetElement(label);
        votedLabel = label;
        }
      }

    it.Set( votedLabel );
    ++it;
    }

  WriteImage<ImageType>( output, outputName.c_str() );
  return 0;
}

template <unsigned int ImageDimension>
int MostLikely( int argc, char *argv[] )
{
  typedef float                                               PixelType;
  typedef itk::Image<PixelType, ImageDimension>               ImageType;
  typedef itk::Image<int, ImageDimension>                     LabeledImageType;
  typedef itk::MinimumMaximumImageCalculator<ImageType>       CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<LabeledImageType> IteratorType;

  if( argc < 5 )
    {
    antscout << " Not enough inputs " << std::endl;
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
  typedef itk::ImageRegionIteratorWithIndex<ImageType>       IteratorType;
  typedef itk::STAPLEImageFilter<ImageType, OutputImageType> StapleFilterType;

  if( argc < 5 )
    {
    antscout << " Not enough inputs " << std::endl;
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
  typename ImageType::Pointer              nullimage( NULL );
  std::vector<typename ImageType::Pointer> images(argc - 5, nullimage );
  typename CalculatorType::Pointer calc = CalculatorType::New();
  int maxLabel = 0;
  for( int i = 5; i < argc; i++ )
    {
    images[i - 5] = ImageType::New();
    ReadImage<ImageType>( images[i - 5], argv[i] );
    stapler->SetInput( i - 5, images[i - 5] );
    antscout << "Input image " << i - 5 << " " << argv[i] << std::endl;

    calc->SetImage( images[i - 5] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }

  antscout << "Examining " << maxLabel << " labels" << std::endl;
  for( int label = 1; label <= maxLabel; label++ )
    {
    char              num[5];
    sprintf( num, "%04d", label );

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
    antscout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string            outputName = std::string( argv[2] );
  std::string::size_type idx;
  idx = outputName.find_first_of('.');
  std::string tempname = outputName.substr(0, idx);
  std::string extension = outputName.substr(idx, outputName.length() );

  // Read input segmentations
  typename ImageType::Pointer              nullimage( NULL );
  std::vector<typename ImageType::Pointer> images(argc - 4, nullimage );

  typename CalculatorType::Pointer calc = CalculatorType::New();
  int maxLabel = 0;
  for( int i = 4; i < argc; i++ )
    {
    images[i - 4] = ImageType::New();
    ReadImage<ImageType>( images[i - 4], argv[i] );
    antscout << "Input image " << i - 4 << " " << argv[i] << std::endl;

    calc->SetImage( images[i - 4] );
    calc->ComputeMaximum();
    if( calc->GetMaximum() > maxLabel )
      {
      maxLabel = calc->GetMaximum();
      }
    }

  antscout << "Examining " << maxLabel << " labels" << std::endl;
  typename OutputImageType::Pointer              nullout( NULL );
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
                                             + 1.0 / images.size() );
        }
      }
    }
  for( int label = 1; label <= maxLabel; label++ )
    {
    char              num[5];
    sprintf( num, "%04d", label );

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
  typedef itk::ImageFileReader<ImageType>                    ImageReaderType;
  typedef itk::ImageFileReader<LabelImageType>               LabelImageReaderType;
  typedef itk::ImageFileWriter<LabelImageType>               LabelImageWriterType;
  typedef itk::MinimumMaximumImageCalculator<LabelImageType> CalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<LabelImageType>  IteratorType;
  typedef itk::NeighborhoodIterator<ImageType>               NeighborhoodIteratorType;

  if( argc < 6 )
    {
    antscout << " Not enough inputs " << std::endl;
    return 1;
    }

  std::string outputName = std::string( argv[2] );

  // Read input images and segmentations
  const int nImages = (argc - 5) / 2;

  int radius = 5;
  if( argc > ( 5 + 2 * nImages ) )
    {
    radius = atoi( argv[5 + 2 * nImages] );
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
    if( maxVotes == nImages )
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
            if( k == 1 )
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
          ( product - k * targetMean * imageMean ) / ( (k - 1) * vcl_sqrt(targetVar) * vcl_sqrt(imageVar) );
        weights.SetElement( i, vcl_fabs(pearson) );
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
  typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4
    <ImageType, ImageType, ImageType>                                          MetricType;

  if( argc < 5 )
    {
    antscout << "ERROR: Not enough inputs " << std::endl;
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
      <ImageType, ImageType, ImageType>                                          MetricType;

    int r = 5;
    if( argc > 6 )
      {
      r = atoi( argv[6] );
      }

    typename MetricType::RadiusType radius;
    radius.Fill( r );

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetRadius( radius );
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );
    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "NormalizedCorrelation") == 0 )
    {
    typedef itk::CorrelationImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );
    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "Demons") == 0 )
    {
    typedef itk::DemonsImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );

    // FIXME - Calling initialize on demons causes seg fault
    antscout << "Demons is currently broken" << std::endl;
    return 1;

    metric->Initialize();
    value = metric->GetValue();
    }
  else if( strcmp(argv[3], "Mattes") == 0 )
    {
    typedef itk::MattesMutualInformationImageToImageMetricv4
      <ImageType, ImageType, ImageType> MetricType;

    int bins = 32;
    if( argc > 6 )
      {
      bins = atoi( argv[6] );
      }

    typename MetricType::Pointer metric = MetricType::New();
    metric->SetFixedImage( img1 );
    metric->SetMovingImage( img2 );
    metric->SetNumberOfHistogramBins( bins );
    metric->Initialize();
    value = metric->GetValue();
    }

  antscout << value << std::endl;

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
    antscout << "ERROR: Not enough inputs " << std::endl;
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

  float mean1 = 0.0;
  float var1 = 0.0;
  float mean2 = 0.0;
  float var2 = 0.0;
  float product = 0.0;
  float k = 0;

  IteratorType it( mask, mask->GetLargestPossibleRegion() );
  while( !it.IsAtEnd() )
    {
    if( it.Value() > 0 )
      {
      k++;
      typename ImageType::IndexType idx = it.GetIndex();
      if( k == 1 )
        {
        mean1 = img1->GetPixel( idx );
        mean2 = img2->GetPixel( idx );
        var1 = 0.0;
        var2 = 0.0;
        }
      else
        {
        float oldMean = mean1;
        float value = img1->GetPixel( idx );
        mean1 = mean1 + (value - mean1) / k;
        var1 = var1 + (value - oldMean) * ( value - mean1 );

        oldMean = mean2;
        value = img2->GetPixel( idx );
        mean2 = mean2 + ( value - mean2 ) / k;
        var2 = var2 + ( value - oldMean) * ( value - mean2 );

        product += img1->GetPixel( idx ) *  img2->GetPixel( idx );
        }
      }
    ++it;
    }

  var1 /= (k - 1);
  var2 /= (k - 1);

  float pearson = ( product - k * mean1 * mean2 ) / ( (k - 1) * vcl_sqrt(var1) * vcl_sqrt(var2) );
  antscout << pearson << std::endl;

  return 0;
}

template <unsigned int ImageDimension>
int MinMaxMean( int argc, char *argv[] )
{
  typedef float                                         PixelType;
  typedef itk::Image<PixelType, ImageDimension>         ImageType;
  typedef itk::ImageFileReader<ImageType>               ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>               ImageWriterType;
  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::ImageMomentsCalculator<ImageType>        MomentsCalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>  IteratorType;

  if( argc < 5 )
    {
    antscout << " Not enough inputs " << std::endl;
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

  antscout << calc->GetMinimum() << " " << calc->GetMaximum() << " " << mean << std::endl;

  return 0;
}

template <unsigned int ImageDimension, class TRealType, class TImageType, class TGImageType, class TInterp>
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
  RealType wt = 1.0 / ( RealType ) Gsz;
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
        avgpoint1[dd] = avgpoint1[dd] + point1[dd] * wt;
        avgpoint2[dd] = avgpoint2[dd] + point2[dd] * wt;
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
  RealType                            evsum = 0, wt0 = 0, wt1 = 0;
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
  wt0 = fabs( eig1.get_eigenvalue( eigind0 ) ) + fabs( eig2.get_eigenvalue( eigind0 ) );
  wt1 = fabs( eig1.get_eigenvalue( eigind1 ) ) + fabs( eig2.get_eigenvalue( eigind1 ) );
  evsum = wt0 + wt1;
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
  //  std::cout <<" V1 " << evec1_primary << std::endl;
  //  std::cout <<" V2 " << evec2_primary << std::endl;
  //  std::cout <<" R(V2) " << A_solution * evec2_primary << " dets " << vnl_determinant<RealType>(  wahba.V()) <<  "
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
      ptran[dd] = ptran2[dd] + avgpoint2[dd];
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

  if( vnl_math_isnan( correlation ) || vnl_math_isinf( correlation )  )
    {
    return 0;
    }
  else
    {
    return correlation;
    }
}

template <class TRealType>
void Sinkhorn( vnl_matrix<TRealType>&  correspondencematrix  )
{
  antscout << " SH begin " << std::endl;
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
  antscout << " SH done " << std::endl;
}

template <class TImage>
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
    antscout << "Usage: " << argv[0] << " ImageDimension";
    antscout << " outputWeightImage PureTissueN4WeightMask";
    antscout << " probImage1 probImage2 ... probImageN" << std::endl;
    return EXIT_FAILURE;
    }

  typedef float       PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>       ReaderType;

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
    PixelType probability = 0.0;
    for( unsigned int i = 0; i < images.size(); i++ )
      {
      PixelType negation = 1.0;
      for( unsigned int j = 0; j < images.size(); j++ )
        {
        if( i == j )
          {
          continue;
          }
        negation *= ( 1.0 - images[j]->GetPixel( ItO.GetIndex() ) );
        }
      probability += negation * images[i]->GetPixel( ItO.GetIndex() );
      }
    ItO.Set( probability );
    }

  WriteImage<ImageType>( output, argv[2] );

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
    antscout << "Usage: " << argv[0] << " ImageDimension";
    antscout << " segmentationImage Check3TissueLabeling";
    antscout << " priorWarped1 priorWarped2 priorWarped3";
    antscout << " posteriorWarped1 posteriorWarped2 posteriorWarped3" << std::endl;
    return EXIT_FAILURE;
    }

  typedef double       PixelType;
  typedef unsigned int LabelType;

  const unsigned int NumberOfLabels = 3;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>       ReaderType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType>  LabelReaderType;

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
  itk::ImageRegionIterator<LabelImageType> ItM( maxPriorLabelImage, maxPriorLabelImage->GetRequestedRegion() );
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
      PixelType prior = priors[d-1]->GetPixel( ItL.GetIndex() );
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
  permutations(which  , 1) = 3;
  permutations(which  , 2) = 2;

  permutations(++which, 0) = 2;
  permutations(which  , 1) = 1;
  permutations(which  , 2) = 3;

  permutations(++which, 0) = 2;
  permutations(which  , 1) = 3;
  permutations(which  , 2) = 1;

  permutations(++which, 0) = 3;
  permutations(which  , 1) = 1;
  permutations(which  , 2) = 2;

  permutations(++which, 0) = 3;
  permutations(which  , 1) = 2;
  permutations(which  , 2) = 1;

  PixelType maxDice = 0.0;
  int maxPermutationRow = -1;
  for( unsigned r = 0; r < 6; r++ )
    {
    typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( labelImage );
    duplicator->Update();

    typename LabelImageType::Pointer permutedLabelImage = duplicator->GetModifiableOutput();

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
    antscout << r << ": " << dice << std::endl;
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
    antscout << d+1 << " -> " << permutations( maxPermutationRow, d ) << std::endl;
    movingLabels[d] = permutations( maxPermutationRow, d );
    fixedLabels[d] = d + 1;
    }

  if( maxPermutationRow == 0 )
    {
    antscout << "No need to change labels/posteriors." << std::endl;
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
        antscout << "Writing posteriors " << movingLabels[d] << " " << argv[7 + d] << std::endl;
        antscout << "Writing posteriors "
                 << movingLabels[movingLabels[d] - 1] << " " << argv[7 + movingLabels[d] - 1] << std::endl;

        WriteImage<ImageType>( posteriors[movingLabels[d] - 1], argv[7 + d] );
        WriteImage<ImageType>( posteriors[movingLabels[movingLabels[d] - 1] - 1], argv[7 + movingLabels[d] - 1] );

        LabelType tmp = movingLabels[d];
        movingLabels[d] = movingLabels[tmp-1];
        movingLabels[tmp-1] = tmp;
        }
      if( d == 0 )
        {
        typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
        typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage( labelImage );
        duplicator->Update();

        typename LabelImageType::Pointer permutedLabelImage = duplicator->GetModifiableOutput();

        itk::ImageRegionIterator<LabelImageType> ItP( permutedLabelImage, permutedLabelImage->GetRequestedRegion() );
        for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
          {
          LabelType permutedLabel = ItP.Get();
          if( permutedLabel != 0 )
            {
            unsigned int whichColumn = permutedLabel - 1;
            ItP.Set( permutations( maxPermutationRow, whichColumn ) );
            }
          }

        antscout << "Relabeling segmentation image." << std::endl;
        WriteImage<LabelImageType>( labelImage, argv[2] );
        }
      }
    }
  return EXIT_SUCCESS;
}

template <class T>
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

template <unsigned int ImageDimension>
int BlobDetector( int argc, char *argv[] )
{
  typedef float                                                                   PixelType;
  typedef float                                                                   RealType;
  typedef itk::Image<PixelType, ImageDimension>                                   ImageType;
  typedef itk::ImageFileReader<ImageType>                                         ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>                                         ImageWriterType;
  typedef itk::MinimumMaximumImageCalculator<ImageType>                           CalculatorType;
  typedef itk::ImageMomentsCalculator<ImageType>                                  MomentsCalculatorType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>                            IteratorType;
  typedef itk::SurfaceImageCurvature<ImageType>                                   ParamType;
  typedef itk::CovariantVector<RealType, ImageDimension>                          GradientPixelType;
  typedef itk::Image<GradientPixelType, ImageDimension>                           GradientImageType;
  typedef itk::SmartPointer<GradientImageType>                                    GradientImagePointer;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer                               GradientImageFilterPointer;
  antscout << " Nearest neighbors should remain nearest neighbors under transformation! - not done " << std::endl;
  if( argc < 5 )
    {
    antscout << " Not enough inputs " << std::endl;
    return 1;
    }
  RealType     gradsig = 1.0;       // sigma for gradient filter
  unsigned int radval = 10;         // radius for correlation
  unsigned int stepsperoctave = 16; // number of steps between doubling of scale
  RealType     minscale = vcl_pow( 1.0, 1 );
  RealType     maxscale = vcl_pow( 2.0, 6 );
  int          argct = 2;
  const std::string  outname = std::string(argv[argct]);
  std::string  outname2 = std::string("temp.nii.gz");
  argct += 2;
  std::string  fn1 = std::string(argv[argct]);   argct++;
  unsigned int nblobs = atoi( argv[argct] );   argct++;
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
  typename ImageType::Pointer image;
  typename ImageType::Pointer image2;
  ReadImage<ImageType>( image, fn1.c_str() );
  GradientImageFilterPointer gfilter = GradientImageFilterType::New();
  gfilter->SetInput( image );
  gfilter->SetSigma( gradsig );
  gfilter->Update();
  GradientImagePointer gimage = gfilter->GetOutput();
  GradientImagePointer gimage2;
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
  BlobsListType blobs2;
  if( fn2.length() > 3 )
    {
    ReadImage<ImageType>( image2, fn2.c_str() );
    GradientImageFilterPointer gfilter2 = GradientImageFilterType::New();
    gfilter2->SetInput( image2 );
    gfilter2->SetSigma( gradsig );
    gfilter2->Update();
    gimage2 = gfilter2->GetOutput();
    typename BlobFilterType::Pointer blobFilter2 = BlobFilterType::New();
    blobFilter2->SetStartT( minscale );
    blobFilter2->SetEndT( maxscale );
    blobFilter2->SetStepsPerOctave( stepsperoctave );
    blobFilter2->SetNumberOfBlobs( nblobs );
    blobFilter2->SetInput( image2 );
    //  blobFilter2->SetInput( ComputeLaplacianImage<ImageType>( image2 ) );
    blobFilter2->Update();
    labimg2 = blobFilter2->GetBlobRadiusImage();
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    labimg->FillBuffer( 0 );
    labimg2->FillBuffer( 0 );
    blobs2 =  blobFilter2->GetBlobs();
    }
  else
    {
    return EXIT_SUCCESS;
    }
  antscout << " Blob1Length " << blobs1.size() << " Blob2Length " << blobs2.size() << std::endl;
  // now compute some feature characteristics in each blob
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
      RealType wt =  exp( -1.0 * dist / radval );
      weights.push_back( wt );
      weightsum += ( wt );
      }
    }
  for( unsigned int ii = 0; ii < weights.size(); ii++ )
    {
    weights[ii] = weights[ii] / weightsum;
    }
  BlobPointer  bestblob = NULL;
  if( ( !blobs2.empty() ) && ( !blobs1.empty() ) )
    {
    unsigned int matchpt = 1;
    unsigned int count2;
    RealType smallval = 1.e-4;
    typedef itk::LinearInterpolateImageFunction<ImageType, float> ScalarInterpolatorType;
    typedef typename ScalarInterpolatorType::Pointer              InterpPointer;
    InterpPointer interp1 =  ScalarInterpolatorType::New();
    interp1->SetInputImage(image);
    InterpPointer interp2 =  ScalarInterpolatorType::New();
    interp2->SetInputImage(image2);
    vnl_matrix<RealType> correspondencematrix( blobs1.size(), blobs2.size() );
    correspondencematrix.fill( 0 );
    unsigned int count1 = 0;
    for( unsigned int i = 0; i < blobs1.size(); i++ )
      {
      BlobPointer blob1 = blobs1[i];
      IndexType   indexi = blob1->GetCenter();
      if( image->GetPixel( indexi ) > smallval )
        {
        GHood.SetLocation( indexi );
        bestblob = NULL;
        count2 = 0;
        for( unsigned int j = 0; j < blobs2.size(); j++ )
          {
          const BlobPointer & blob2 = blobs2[j];
          const IndexType   & indexj = blob2->GetCenter();
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
        antscout << " Progress : " << (float ) i / (float) blobs1.size() * 100.0 << std::endl;
        }
      }
    antscout << " now compute pairwise matching " << correspondencematrix.max_value() << " reducing to " << corrthresh
             << " count1 " << count1 << " count2 " << count2 << std::endl;
    count1 = 0;
    typedef std::pair<BlobPointer, BlobPointer> BlobPairType;
    std::vector<BlobPairType> blobpairs;
    while( ( matchpt < ( corrthresh + 1 ) ) && ( count1 <  blobs1.size() ) )
      {
      unsigned int maxpair = correspondencematrix.arg_max();
      unsigned int maxrow = ( unsigned int )  maxpair / correspondencematrix.cols();
      unsigned int maxcol = maxpair - maxrow * correspondencematrix.cols();
      BlobPointer blob1 = blobs1[maxrow];
      bestblob = blobs2[maxcol];
      if( bestblob )
        {
        if( fabs( bestblob->GetObjectRadius() - blob1->GetObjectRadius() ) < 0.5 )
          {
          if( bestblob && ( image->GetPixel( blob1->GetCenter() ) > smallval )  &&
              ( image2->GetPixel( bestblob->GetCenter() )  > smallval ) )
            {
            BlobPairType blobpairing = std::make_pair( blob1, bestblob );
            blobpairs.push_back( blobpairing );
            antscout << " best correlation " << correspondencematrix.absolute_value_max() << " rad1 "
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
        distspre.push_back( dist1 );
        distspost.push_back( dist2 );
        distspreind.push_back(  bp2 );
        distspostind.push_back( bp2 );
        //        antscout << " blob " << bp << " vs " << bp2 << sqrt( dist1 ) << " v " << sqrt( dist2 ) << std::endl;
        distmatpre( bp, bp2 )  = distmatpre(  bp2, bp ) = dist1;
        distmatpost( bp, bp2 ) = distmatpost( bp2, bp ) = dist2;
        }
      bool kneighborhoodequal = false;
      std::sort( distspreind.begin(), distspreind.end(), blob_index_cmp<std::vector<RealType> &>( distspre  ) );
      std::sort( distspostind.begin(), distspostind.end(), blob_index_cmp<std::vector<RealType> &>( distspost ) );
      //      std::cout << distspreind[0] << "  " << distspreind[1] << " dist0 " <<   distspreind[  distspreind[0] ]  <<
      // std::endl;
      //      std::cout << distspostind[0] << "  " << distspostind[1] << " dist0 " << distspostind[  distspostind[0] ]
      //  <<  std::endl;
      if( ( distspostind[1] == distspreind[1] )   ||
          ( distspostind[1] == distspreind[2] )   ||
          ( distspostind[1] == distspreind[3] )  )
        {
        kneighborhoodequal = true;
        }
      if( !kneighborhoodequal )
        {
        labimg->SetPixel(  blobind, 0 );     // ( int ) ( 0.5 +   ( *i )->GetObjectRadius() ) );
        labimg2->SetPixel( blobpairind, 0 ); // ( int ) ( 0.5 + bestblob->GetObjectRadius() ) );
        }
      antscout << " blob " << bp << " keep " << kneighborhoodequal << std::endl;
      }
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "distpre.csv" );
      writer->SetInput( &distmatpre );
      writer->Write();
      }
      {
      typedef itk::CSVNumericObjectFileWriter<RealType, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( "distpost.csv" );
      writer->SetInput( &distmatpost );
      writer->Write();
      }
    WriteImage<BlobRadiusImageType>( labimg, outname.c_str() );
    WriteImage<BlobRadiusImageType>( labimg2, outname2.c_str() );
    antscout << " Matched " << matchpt << " blobs " << std::endl;
    }
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ImageMath( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ImageMath" );

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
  argv[argc] = 0;
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
  antscout->set_stream( out_stream );

  if( argc < 5 )
    {
    antscout << "\nUsage: " << argv[0]
             << " ImageDimension <OutputImage.ext> [operations and inputs] <Image1.ext> <Image2.ext>" << std::endl;

    antscout << "\nUsage Information " << std::endl;
    antscout << " ImageDimension: 2 or 3 (for 2 or 3 dimensional operations)." << std::endl;
    antscout << " ImageDimension: 4 (for operations on 4D file, e.g. time-series data)." << std::endl;
    antscout << " Operator: See list of valid operators below." << std::endl;
    antscout << " The last two arguments can be an image or float value " << std::endl;
    antscout << " NB: Some options output text files" << std::endl;

    antscout << "\nMathematical Operations:" << std::endl;
    antscout << "  m            : Multiply ---  use vm for vector multiply " << std::endl;
    antscout << "  +             : Add ---  use v+ for vector add " << std::endl;
    antscout << "  -             : Subtract ---  use v- for vector subtract " << std::endl;
    antscout << "  /             : Divide" << std::endl;
    antscout << "  ^            : Power" << std::endl;
    antscout << "  exp            : Take exponent exp(imagevalue*value)" << std::endl;
    antscout << "  addtozero        : add image-b to image-a only over points where image-a has zero values"
             << std::endl;
    antscout << "  overadd        : replace image-a pixel with image-b pixel if image-b pixel is non-zero"
             << std::endl;
    antscout << "  abs            : absolute value " << std::endl;
    antscout
      << "  total            : Sums up values in an image or in image1*image2 (img2 is the probability mask)"
      << std::endl;
    antscout
      << "  mean            :  Average of values in an image or in image1*image2 (img2 is the probability mask)"
      << std::endl;
    antscout
      << "  vtotal            : Sums up volumetrically weighted values in an image or in image1*image2 (img2 is the probability mask)"
      << std::endl;
    antscout << "  Decision        : Computes result=1./(1.+exp(-1.0*( pix1-0.25)/pix2))" << std::endl;
    antscout << "  Neg            : Produce image negative" << std::endl;

    antscout << "\nSpatial Filtering:" <<  std::endl;
    antscout << "  G Image1.ext s    : Smooth with Gaussian of sigma = s" << std::endl;
    antscout << "  MD Image1.ext s    : Morphological Dilation with radius s" << std::endl;
    antscout << "  ME Image1.ext s    : Morphological Erosion with radius s" << std::endl;
    antscout << "  MO Image1.ext s    : Morphological Opening with radius s" << std::endl;
    antscout << "  MC Image1.ext s    : Morphological Closing with radius s" << std::endl;
    antscout << "  GD Image1.ext s    : Grayscale Dilation with radius s" << std::endl;
    antscout << "  GE Image1.ext s    : Grayscale Erosion with radius s" << std::endl;
    antscout << "  GO Image1.ext s    : Grayscale Opening with radius s" << std::endl;
    antscout << "  GC Image1.ext s    : Grayscale Closing with radius s" << std::endl;
    antscout
      <<
      "  BlobDetector Image1.ext NumberOfBlobs  Optional-Input-Image2 Blob-2-out.nii.gz N-Blobs-To-Match  :  blob detection by searching for local extrema of the Laplacian of the Gassian (LoG) "
      << std::endl;

    antscout << "\nTime Series Operations:" << std::endl;
    antscout
      <<
      " CompCorrAuto : Outputs a csv file containing global signal vector and N comp-corr eigenvectors determined from PCA of the high-variance voxels.  Also outputs a comp-corr + global signal corrected 4D image as well as a 3D image measuring the time series variance.  Requires a label image with label 1 identifying voxels in the brain."
      << std::endl;
    antscout
      <<
      "   ImageMath 4 ${out}compcorr.nii.gz ThreeTissueConfounds ${out}.nii.gz  ${out}seg.nii.gz 1 3  "
      << " : Outputs average global, CSF and WM signals.  Requires a label image with 3 labels , csf, gm , wm ."
      << std::endl;
    antscout << "    Usage        : ThreeTissueConfounds 4D_TimeSeries.nii.gz LabeLimage.nii.gz  csf-label wm-label "
             << std::endl;
    antscout
      << " TimeSeriesSubset : Outputs n 3D image sub-volumes extracted uniformly from the input time-series 4D image."
      << std::endl;
    antscout << "    Usage        : TimeSeriesSubset 4D_TimeSeries.nii.gz n " << std::endl;
    antscout
      << " TimeSeriesDisassemble : Outputs n 3D image volumes for each time-point in time-series 4D image."
      << std::endl;
    antscout << "    Usage        : TimeSeriesDisassemble 4D_TimeSeries.nii.gz " << std::endl
             << std::endl;
    antscout
      << " TimeSeriesAssemble : Outputs a 4D time-series image from a list of 3D volumes."
      << std::endl;
    antscout << "    Usage        : TimeSeriesAssemble time_spacing time_origin *images.nii.gz " << std::endl;
    antscout
      <<
      " TimeSeriesToMatrix : Converts a 4D image + mask to matrix (stored as csv file) where rows are time and columns are space ."
      << std::endl;
    antscout << "    Usage        : TimeSeriesToMatrix 4D_TimeSeries.nii.gz mask " << std::endl;
    antscout
      << " TimeSeriesSimpleSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes."
      << std::endl;
    antscout << "    Usage        : TimeSeriesSimpleSubtraction image.nii.gz " << std::endl;
    antscout
      << " TimeSeriesSurroundSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes."
      << std::endl;
    antscout << "    Usage        : TimeSeriesSurroundSubtraction image.nii.gz " << std::endl;
    antscout
      << " TimeSeriesSincSubtraction : Outputs a 3D mean pair-wise difference list of 3D volumes."
      << std::endl;
    antscout << "    Usage        : TimeSeriesSincSubtraction image.nii.gz " << std::endl;
    antscout
      << " SplitAlternatingTimeSeries : Outputs 2 3D time series"
      << std::endl;
    antscout << "    Usage        : SplitAlternatingTimeSeries image.nii.gz " << std::endl;


    antscout
      <<
      " ComputeTimeSeriesLeverage : Outputs a csv file that identifies the raw leverage and normalized leverage for each time point in the 4D image.  leverage, here, is the difference of the time-point image from the average of the n images.  the normalized leverage is =  average( sum_k abs(Leverage(t)-Leverage(k)) )/Leverage(t). "
      << std::endl;
    antscout << "    Usage        : ComputeTimeSeriesLeverage 4D_TimeSeries.nii.gz k_neighbors " << std::endl;

    antscout
      <<
      " PASL : computes the PASL model of CBF  "     << std::endl <<  "f =  \frac{      lambda DeltaM        } "
      << std::endl
      << " {     2 \alpha M_0 TI_1 exp( - TI_2 / T_{1a} )  } " << std::endl;
    antscout
      << "    Usage        : PASL 3D/4D_TimeSeries.nii.gz BoolFirstImageIsControl M0Image parameter_list.txt "
      << std::endl;

    antscout
      <<
      " pCASL : computes the pCASL model of CBF  "     << std::endl
      << " f =  \frac{      lambda DeltaM R_{1a}        }  " << std::endl
      << "  {     2 \alpha M_0 [ exp( - w R_{1a} ) - exp( -w ( \tau + w ) R_{1a}) ]     } " << std::endl;
    antscout << "    Usage        : pCASL 3D/4D_TimeSeries.nii.gz parameter_list.txt " << std::endl;
    antscout
      << " PASLQuantifyCBF : Outputs a 3D CBF image in ml/100g/min from a magnetization ratio image"
      << std::endl;
    antscout
      <<
      "    Usage        : PASLQuantifyCBF mag_ratio.nii.gz [TI1=700] [TI2=1900] [T1blood=1664] [Lambda=0.9] [Alpha=0.95] [SliceDelay-45] "
      << std::endl;

    antscout << "\nTensor Operations:" << std::endl;
    antscout << "  4DTensorTo3DTensor    : Outputs a 3D_DT_Image with the same information. " << std::endl;
    antscout << "    Usage        : 4DTensorTo3DTensor 4D_DTImage.ext" << std::endl;
    antscout << "  ComponentTo3DTensor    : Outputs a 3D_DT_Image with the same information as component images. "
             << std::endl;
    antscout << "    Usage        : ComponentTo3DTensor component_image_prefix[xx,xy,xz,yy,yz,zz] extension"
             << std::endl;
    antscout << "  ExtractComponentFrom3DTensor    : Outputs a component images. " << std::endl;
    antscout << "    Usage        : ExtractComponentFrom3DTensor dtImage.ext which={xx,xy,xz,yy,yz,zz}" << std::endl;
    antscout << "  ExtractVectorComponent: Produces the WhichVec component of the vector " << std::endl;
    antscout << "    Usage        : ExtractVectorComponent VecImage WhichVec" << std::endl;
    antscout << "  TensorColor        : Produces RGB values identifying principal directions " << std::endl;
    antscout << "    Usage        : TensorColor DTImage.ext" << std::endl;
    antscout << "  TensorFA        : " << std::endl;
    antscout << "    Usage        : TensorFA DTImage.ext" << std::endl;
    antscout << "  TensorFADenominator    : " << std::endl;
    antscout << "    Usage        : TensorFADenominator DTImage.ext" << std::endl;
    antscout << "  TensorFANumerator    : " << std::endl;
    antscout << "    Usage        : TensorFANumerator DTImage.ext" << std::endl;
    antscout << "  TensorIOTest    : Will write the DT image back out ... tests I/O processes for consistency. "
             << std::endl;
    antscout << "    Usage        : TensorIOTest DTImage.ext" << std::endl;
    antscout << "  TensorMeanDiffusion      : Mean of the eigenvalues" << std::endl;
    antscout << "    Usage        : TensorMeanDiffusion DTImage.ext" << std::endl;
    antscout << "  TensorRadialDiffusion    : Mean of the two smallest eigenvalues" << std::endl;
    antscout << "    Usage        : TensorRadialDiffusion DTImage.ext" << std::endl;
    antscout << "  TensorAxialDiffusion     : Largest eigenvalue, equivalent to TensorEigenvalue DTImage.ext 2"
             << std::endl;
    antscout << "    Usage        : TensorAxialDiffusion DTImage.ext" << std::endl;
    antscout << "  TensorEigenvalue         : Gets a single eigenvalue 0-2, where 0 = smallest, 2 = largest"
             << std::endl;
    antscout << "    Usage        : TensorEigenvalue DTImage.ext WhichInd" << std::endl;
    antscout
      <<
      "  TensorToVector    : Produces vector field identifying one of the principal directions, 2 = largest eigenvalue"
      << std::endl;
    antscout << "    Usage        : TensorToVector DTImage.ext WhichVec" << std::endl;
    antscout
      <<
      "  TensorToVectorComponent: 0 => 2 produces component of the principal vector field (largest eigenvalue). 3 = 8 => gets values from the tensor "
      << std::endl;
    antscout << "    Usage        : TensorToVectorComponent DTImage.ext WhichVec" << std::endl;
    antscout << "  TensorMask     : Mask a tensor image, sets background tensors to zero " << std::endl;
    antscout << "    Usage        : TensorMask DTImage.ext mask.ext" << std::endl;

 
    antscout << "\nLabel Fusion:" << std::endl;
    antscout << "  MajorityVoting : Select label with most votes from candidates" << std::endl;
    antscout << "    Usage: MajorityVoting LabelImage1.nii.gz .. LabelImageN.nii.gz" << std::endl;
    antscout << "  CorrelationVoting : Select label with local correlation weights" << std::endl;
    antscout << "    Usage: CorrelationVoting Template.ext IntenistyImages* LabelImages* {Optional-Radius=5}"
             << std::endl;
    antscout << "  STAPLE : Select label using STAPLE method" << std::endl;
    antscout << "    Usage: STAPLE confidence-weighting LabelImages*" << std::endl;
    antscout << "    Note:  Gives probabilistic output (float)" << std::endl;
    antscout << "  MostLikely : Select label from from maximum probabilistic segmentations" << std::endl;
    antscout << "    Usage: MostLikely ProbabilityImages*" << std::endl;
    antscout << "  AverageLabels : Select label using STAPLE method" << std::endl;
    antscout << "    Usage: STAPLE LabelImages*" << std::endl;
    antscout << "    Note:  Gives probabilistic output (float)" << std::endl;

    antscout << "\nImage Metrics & Info:" <<  std::endl;
    antscout << "  PearsonCorrelation: r-value from intesities of two images" << std::endl;
    antscout << "    Usage: PearsonCorrelation image1.ext image2.ext {Optional-mask.ext}" << std::endl;
    antscout << "  NeighborhoodCorrelation: local correlations" << std::endl;
    antscout << "    Usage: NeighborhoodCorrelation image1.ext image2.ext {Optional-radius=5}" << std::endl;
    antscout << "  NormalizedCorrelation: r-value from intesities of two images" << std::endl;
    antscout << "    Usage: NormalizedCorrelation image1.ext image2.ext" << std::endl;
    antscout << "  Demons: " << std::endl;
    antscout << "    Usage: Demons image1.ext image2.ext" << std::endl;
    antscout << "  Mattes: mutual information" << std::endl;
    antscout << "    Usage: Mattes image1.ext image2.ext {Optional-number-bins=32}" << std::endl;

    antscout << "\nUnclassified Operators:" << std::endl;

    antscout << "  Byte            : Convert to Byte image in [0,255]" << std::endl;

    antscout
      << "\n  CompareHeadersAndImages: Tries to find and fix header errors. Outputs a repaired image with new header. "
      << std::endl;
    antscout << "                Never use this if you trust your header information. " << std::endl;
    antscout << "      Usage        : CompareHeadersAndImages Image1 Image2" << std::endl;

    antscout
      <<
      "\n  ConvertImageSetToMatrix: Each row/column contains image content extracted from mask applied to images in *img.nii "
      << std::endl;
    antscout << "      Usage        : ConvertImageSetToMatrix rowcoloption Mask.nii *images.nii" << std::endl;
    antscout << " ConvertImageSetToMatrix output can be an image type or csv file type." << std::endl;

    antscout << "\n  RandomlySampleImageSetToCSV: N random samples are selected from each image in a list "
             << std::endl;
    antscout << "      Usage        : RandomlySampleImageSetToCSV N_samples *images.nii" << std::endl;
    antscout << " RandomlySampleImageSetToCSV outputs a csv file type." << std::endl;

    antscout
      <<
      "\n  FrobeniusNormOfMatrixDifference: take the difference between two itk-transform matrices and then compute the frobenius norm"
      << std::endl;
    antscout << "      Usage        : FrobeniusNormOfMatrixDifference mat1 mat2 " << std::endl;
    antscout
      <<
      "\n  ConvertImageSetToEigenvectors: Each row/column contains image content extracted from mask applied to images in *img.nii "
      << std::endl;
    antscout << "      Usage        : ConvertImageSetToEigenvectors N_Evecs Mask.nii *images.nii" << std::endl;
    antscout << " ConvertImageSetToEigenvectors output will be a csv file for each label value > 0 in the mask."
             << std::endl;

    antscout << "\n  ConvertImageToFile    : Writes voxel values to a file  " << std::endl;
    antscout << "      Usage        : ConvertImageToFile imagevalues.nii {Optional-ImageMask.nii}" << std::endl;

    antscout
      << "\n  ConvertLandmarkFile    : Converts landmark file between formats. See ANTS.pdf for description of formats."
      << std::endl;
    antscout << "      Usage        : ConvertLandmarkFile InFile.txt" << std::endl;
    antscout << "      Example 1        : ImageMath 3  outfile.vtk  ConvertLandmarkFile  infile.txt" << std::endl;

    antscout << "\n  ConvertToGaussian    : " << std::endl;
    antscout << "      Usage        : ConvertToGaussian  TValueImage  sigma-float" << std::endl;

    antscout
      <<
      "\n  ConvertVectorToImage    : The vector contains image content extracted from a mask. Here the vector is returned to its spatial origins as image content "
      << std::endl;
    antscout << "      Usage        : ConvertVectorToImage Mask.nii vector.nii" << std::endl;

    antscout << "\n  CorrelationUpdate    : In voxels, compute update that makes Image2 more like Image1."
             << std::endl;
    antscout << "      Usage        : CorrelationUpdate Image1.ext Image2.ext RegionRadius" << std::endl;

    antscout << "\n  CountVoxelDifference    : The where function from IDL " << std::endl;
    antscout << "      Usage        : CountVoxelDifference Image1 Image2 Mask" << std::endl;

    antscout << "\n  CorruptImage        : " << std::endl;
    antscout << "      Usage        : CorruptImage Image NoiseLevel Smoothing" << std::endl;

    antscout << "\n  D            : DistanceTransform" << std::endl;

    antscout
      <<
      "\n  DiceAndMinDistSum    : Outputs DiceAndMinDistSum and Dice Overlap to text log file + optional distance image"
      << std::endl;
    antscout << "      Usage        : DiceAndMinDistSum LabelImage1.ext LabelImage2.ext OptionalDistImage"
             << std::endl;

    antscout << "\n  EnumerateLabelInterfaces: " << std::endl;
    antscout
      << "      Usage        : EnumerateLabelInterfaces ImageIn ColoredImageOutname NeighborFractionToIgnore"
      << std::endl;

    antscout
      << "\n  ClusterThresholdVariate        :  for sparse estimation "
      << std::endl;
    antscout << "      Usage        : ClusterThresholdVariate image mask  MinClusterSize" << std::endl;

    antscout
      << "\n  ExtractSlice        : Extracts slice number from last dimension of volume (2,3,4) dimensions "
      << std::endl;
    antscout << "      Usage        : ExtractSlice volume.nii.gz slicetoextract" << std::endl;

    antscout
      <<
      "\n  FastMarchingSegmentation: final output is the propagated label image. Optional stopping value: higher values allow more distant propagation "
      << std::endl;
    antscout
      <<
      "      Usage        : FastMarchingSegmentation speed/binaryimagemask.ext initiallabelimage.ext Optional-Stopping-Value"
      << std::endl;

    antscout << "\n  FillHoles        : Parameter = ratio of edge at object to edge at background;  --  " << std::endl;
    antscout << "                Parameter = 1 is a definite hole bounded by object only, 0.99 is close" << std::endl;
    antscout << "                Default of parameter > 1 will fill all holes" << std::endl;
    antscout << "      Usage        : FillHoles Image.ext parameter" << std::endl;

    antscout << "  Finite            : replace non-finite values with finite-value (default = 0)" << std::endl;
    antscout << "      Usage        : Finite Image.exdt {replace-value=0}" << std::endl;

    antscout << "\n  FitSphere        : " << std::endl;
    antscout << "      Usage        : FitSphere GM-ImageIn {WM-Image} {MaxRad-Default=5}" << std::endl;

    antscout << "\n  FlattenImage        : Replaces values greater than %ofMax*Max to the value %ofMax*Max "
             << std::endl;
    antscout << "      Usage        : FlattenImage Image %ofMax" << std::endl;

    antscout << "\n  GetLargestComponent    : Get the largest object in an image" << std::endl;
    antscout << "      Usage        : GetLargestComponent InputImage {MinObjectSize}" << std::endl;

    antscout << "\n  Grad            : Gradient magnitude with sigma s (if normalize, then output in range [0, 1])"
             << std::endl;
    antscout << "      Usage        : Grad Image.ext s normalize?" << std::endl;

    antscout << "\n  HistogramMatch    : " << std::endl;
    antscout
      <<
      "      Usage        : HistogramMatch SourceImage ReferenceImage {NumberBins-Default=255} {NumberPoints-Default=64}"
      << std::endl;

    antscout
      <<
      "\n  InvId            : computes the inverse-consistency of two deformations and write the inverse consistency error image "
      << std::endl;
    antscout << "      Usage        : InvId VectorFieldName VectorFieldName" << std::endl;

    antscout << "\n  LabelStats        : Compute volumes / masses of objects in a label image. Writes to text file"
             << std::endl;
    antscout << "      Usage        : LabelStats labelimage.ext valueimage.nii" << std::endl;

    antscout
      << "\n  Laplacian        : Laplacian computed with sigma s (if normalize, then output in range [0, 1])"
      << std::endl;
    antscout << "      Usage        : Laplacian Image.ext s normalize?" << std::endl;

    antscout << "\n  Lipschitz        : Computes the Lipschitz norm of a vector field " << std::endl;
    antscout << "      Usage        : Lipschitz VectorFieldName" << std::endl;

    antscout << "\n  MakeImage        : " << std::endl;
    antscout << "      Usage        : MakeImage SizeX  SizeY {SizeZ};" << std::endl;

    antscout
      << "\n  MTR        : Computes the magnetization transfer ratio ( (M0-M1)/M0 ) and truncates values to [0,1]"
      << std::endl;
    antscout << "      Usage        : MTR M0Image M1Image [MaskImage];" << std::endl;

    antscout << "\n  Normalize        : Normalize to [0,1]. Option instead divides by average value" << std::endl;
    antscout << "      Usage        : Normalize Image.ext opt" << std::endl;

    antscout << "\n  PadImage       : If Pad-Number is negative, de-Padding occurs" << std::endl;
    antscout << "      Usage        : PadImage ImageIn Pad-Number" << std::endl;

    antscout << "\n  SigmoidImage   : " << std::endl;
    antscout << "      Usage        : SigmoidImage ImageIn [alpha=1.0] [beta=0.0]" << std::endl;

    antscout << "\n  CenterImage2inImage1        : " << std::endl;
    antscout << "      Usage       : ReferenceImageSpace ImageToCenter " << std::endl;

    antscout << "\n  PH            : Print Header" << std::endl;

    antscout
      << "\n  PoissonDiffusion        : Solves Poisson's equation in a designated region using non-zero sources"
      << std::endl;
    antscout
      <<
      "      Usage        : PoissonDiffusion inputImage labelImage [sigma=1.0] [regionLabel=1] [numberOfIterations=500] [convergenceThreshold=1e-10]"
      << std::endl;

    antscout
      <<
      "\n  PropagateLabelsThroughMask: Final output is the propagated label image. Optional stopping value: higher values allow more distant propagation"
      << std::endl;
    antscout
      <<
      "      Usage        : PropagateLabelsThroughMask speed/binaryimagemask.nii.gz initiallabelimage.nii.gz Optional-Stopping-Value"
      << std::endl;

    antscout << "\n  PValueImage        : " << std::endl;
    antscout << "      Usage        : PValueImage TValueImage dof" << std::endl;

    antscout << "\n  RemoveLabelInterfaces: " << std::endl;
    antscout << "      Usage        : RemoveLabelInterfaces ImageIn" << std::endl;

    antscout
      <<
      "\n  ROIStatistics        : computes anatomical locations, cluster size and mass of a stat image which should be in the same physical space (but not nec same resolution) as the label image."
      << std::endl;
    antscout << "      Usage        : ROIStatistics LabelNames.txt labelimage.ext valueimage.nii" << std::endl;

    antscout << "\n  SetOrGetPixel    : "  << std::endl;
    antscout << "      Usage        : SetOrGetPixel ImageIn Get/Set-Value IndexX IndexY {IndexZ}" << std::endl;
    antscout
      << "      Example 1        : ImageMath 2 outimage.nii SetOrGetPixel Image Get 24 34; Gets the value at 24, 34"
      << std::endl;
    antscout
      <<
      "      Example 2        : ImageMath 2 outimage.nii SetOrGetPixel Image 1.e9 24 34; This sets 1.e9 as the value at 23 34"
      << std::endl;
    antscout << "                You can also pass a boolean at the end to force the physical space to be used"
             << std::endl;

    antscout
      << "\n  Segment        : Segment an Image  with option of Priors, weight 1 => maximally local/prior-based )"
      << std::endl;
    antscout
      <<
      "      Usage        : Segment Image1.ext N-Classes LocalityVsGlobalityWeight-In-ZeroToOneRange OptionalPriorImages"
      << std::endl;

    antscout << "\n  stack            : Will put 2 images in the same volume" << std::endl;
    antscout << "      Usage        : Stack Image1.ext Image2.ext" << std::endl;

    antscout << "\n  ThresholdAtMean    : See the code" << std::endl;
    antscout << "      Usage        : ThresholdAtMean Image %ofMean" << std::endl;

    antscout << "\n  TileImages    : " << std::endl;
    antscout << "      Usage        : TileImages NumColumns ImageList*" << std::endl;

    antscout << "\n  TriPlanarView    : " << std::endl;
    antscout
      <<
      "      Usage        : TriPlanarView  ImageIn.nii.gz PercentageToClampLowIntensity PercentageToClampHiIntensity x-slice y-slice z-slice"
      << std::endl;

    antscout << "\n  TruncateImageIntensity: " << std::endl;
    antscout
      <<
      "      Usage        : TruncateImageIntensity InputImage.ext {lowerQuantile=0.05} {upperQuantile=0.95} {numberOfBins=65} {binary-maskImage}"
      << std::endl;

    antscout << "\n  Where            : The where function from IDL" << std::endl;
    antscout << "      Usage        : Where Image ValueToLookFor maskImage-option tolerance" << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  std::string operation = std::string(argv[3]);

  switch( atoi(argv[1]) )
    {
    case 2:
      {
      if( strcmp(operation.c_str(), "m") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mresample") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "+") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "-") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vm") == 0 )
        {
        VImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vmresample") == 0 )
        {
        VImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v+") == 0 )
        {
        VImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v-") == 0 )
        {
        VImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "/") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "^") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "exp") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "abs") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "addtozero") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "overadd") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "total") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vtotal") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mean") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Decision") == 0 )
        {
        ImageMath<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Neg") == 0 )
        {
        NegativeImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "G") == 0 )
        {
        SmoothImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MD") == 0 || strcmp(operation.c_str(), "ME") == 0 )
        {
        MorphImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MO") == 0 || strcmp(operation.c_str(), "MC") == 0 )
        {
        MorphImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GD") == 0 || strcmp(operation.c_str(), "GE") == 0 )
        {
        MorphImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GO") == 0 || strcmp(operation.c_str(), "GC") == 0 )
        {
        MorphImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "D") == 0 )
        {
        DistanceMap<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Normalize") == 0 )
        {
        NormalizeImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Grad") == 0 )
        {
        GradientImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Laplacian") == 0 )
        {
        LaplacianImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PH") == 0 )
        {
        PrintHeader<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CenterImage2inImage1") == 0 )
        {
        CenterImage2inImage1<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Byte") == 0 )
        {
        ByteImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "LabelStats") == 0 )
        {
        LabelStats<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ROIStatistics") == 0 )
        {
        ROIStatistics<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "DiceAndMinDistSum") == 0 )
        {
        DiceAndMinDistSum<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Lipschitz") == 0 )
        {
        Lipschitz<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "InvId") == 0 )
        {
        InvId<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GetLargestComponent") == 0 )
        {
        GetLargestComponent<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractVectorComponent") == 0 )
        {
        ExtractVectorComponent<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ThresholdAtMean") == 0 )
        {
        ThresholdAtMean<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FlattenImage") == 0 )
        {
        FlattenImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorruptImage") == 0 )
        {
        CorruptImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TileImages") == 0 )
        {
        TileImages<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Where") == 0 )
        {
        Where<2>(argc, argv);
        }
      else if( strcmp( operation.c_str(), "Finite") == 0 )
        {
        Finite<2>( argc, argv );
        }
      else if( strcmp(operation.c_str(), "FillHoles") == 0 )
        {
        FillHoles<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "HistogramMatch") == 0 )
        {
        HistogramMatching<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PadImage") == 0 )
        {
        PadImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SigmoidImage") == 0 )
        {
        SigmoidImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SetOrGetPixel") == 0 )
        {
        SetOrGetPixel<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MakeImage") == 0 )
        {
        MakeImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "stack") == 0 )
        {
        StackImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CompareHeadersAndImages") == 0 )
        {
        CompareHeadersAndImages<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CountVoxelDifference") == 0 )
        {
        CountVoxelDifference<2>(argc, argv);
        }
      //     else if (strcmp(operation.c_str(),"AddToZero") == 0 )  AddToZero<2>(argc,argv);
      else if( strcmp(operation.c_str(), "RemoveLabelInterfaces") == 0 )
        {
        RemoveLabelInterfaces<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PoissonDiffusion") == 0 )
        {
        PoissonDiffusion<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "EnumerateLabelInterfaces") == 0 )
        {
        EnumerateLabelInterfaces<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageToFile") == 0 )
        {
        ConvertImageToFile<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PValueImage") == 0 )
        {
        PValueImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationUpdate") == 0 )
        {
        CorrelationUpdate<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToMatrix") == 0 )
        {
        ConvertImageSetToMatrix<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "RandomlySampleImageSetToCSV") == 0 )
        {
        RandomlySampleImageSetToCSV<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToEigenvectors") == 0 )
        {
        ConvertImageSetToEigenvectors<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertVectorToImage") == 0 )
        {
        ConvertVectorToImage<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PropagateLabelsThroughMask") == 0 )
        {
        PropagateLabelsThroughMask<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FastMarchingSegmentation") == 0 )
        {
        FastMarchingSegmentation<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TruncateImageIntensity") == 0 )
        {
        TruncateImageIntensity<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractSlice") == 0 )
        {
        ExtractSlice<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ClusterThresholdVariate") == 0 )
        {
        ClusterThresholdVariate<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "STAPLE") == 0 )
        {
        STAPLE<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "AverageLabels") == 0 )
        {
        AverageLabels<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MajorityVoting") == 0 )
        {
        MajorityVoting<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MostLikely") == 0 )
        {
        MostLikely<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationVoting") == 0 )
        {
        CorrelationVoting<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PearsonCorrelation") == 0 )
        {
        PearsonCorrelation<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NeighborhoodCorrelation") == 0 )
        {
        ImageMetrics<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NormalizedCorrelation") == 0 )
        {
        ImageMetrics<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Demons") == 0 )
        {
        ImageMetrics<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Mattes") == 0 )
        {
        ImageMetrics<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MinMaxMean") == 0 )
        {
        MinMaxMean<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Check3TissueLabeling") == 0 )
        {
        Check3TissueLabeling<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PureTissueN4WeightMask") == 0 )
        {
        PureTissueN4WeightMask<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "BlobDetector") == 0 )
        {
        BlobDetector<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesRegionSCCA") == 0 )
        {
        TimeSeriesRegionSCCA<2>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesRegionCorr") == 0 )
        {
        TimeSeriesRegionCorr<2>(argc, argv);
        }
      //     else if (strcmp(operation.c_str(),"ConvertLandmarkFile") == 0)  ConvertLandmarkFile<2>(argc,argv);
      else
        {
        antscout << " cannot find operation : " << operation << std::endl;
        }
      }
      break;
    case 3:
      {
      if( strcmp(operation.c_str(), "m") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mresample") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "+") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "-") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vm") == 0 )
        {
        VImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vmresample") == 0 )
        {
        VImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v+") == 0 )
        {
        VImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v-") == 0 )
        {
        VImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "/") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "^") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "exp") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "abs") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "addtozero") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "overadd") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "total") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vtotal") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mean") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Decision") == 0 )
        {
        ImageMath<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Neg") == 0 )
        {
        NegativeImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "G") == 0 )
        {
        SmoothImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MD") == 0 || strcmp(operation.c_str(), "ME") == 0 )
        {
        MorphImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MO") == 0 || strcmp(operation.c_str(), "MC") == 0 )
        {
        MorphImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GD") == 0 || strcmp(operation.c_str(), "GE") == 0 )
        {
        MorphImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GO") == 0 || strcmp(operation.c_str(), "GC") == 0 )
        {
        MorphImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "D") == 0 )
        {
        DistanceMap<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Normalize") == 0 )
        {
        NormalizeImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Grad") == 0 )
        {
        GradientImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Laplacian") == 0 )
        {
        LaplacianImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PH") == 0 )
        {
        PrintHeader<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CenterImage2inImage1") == 0 )
        {
        CenterImage2inImage1<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Byte") == 0 )
        {
        ByteImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "LabelStats") == 0 )
        {
        LabelStats<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ROIStatistics") == 0 )
        {
        ROIStatistics<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "DiceAndMinDistSum") == 0 )
        {
        DiceAndMinDistSum<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Lipschitz") == 0 )
        {
        Lipschitz<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "InvId") == 0 )
        {
        InvId<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GetLargestComponent") == 0 )
        {
        GetLargestComponent<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractVectorComponent") == 0 )
        {
        ExtractVectorComponent<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ThresholdAtMean") == 0 )
        {
        ThresholdAtMean<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FlattenImage") == 0 )
        {
        FlattenImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorruptImage") == 0 )
        {
        CorruptImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "RemoveLabelInterfaces") == 0 )
        {
        RemoveLabelInterfaces<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "EnumerateLabelInterfaces") == 0 )
        {
        EnumerateLabelInterfaces<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FitSphere") == 0 )
        {
        FitSphere<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Where") == 0 )
        {
        Where<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SmoothTensorImage") == 0 )
        {
        SmoothTensorImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorFA") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorFANumerator") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorFADenominator") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "4DTensorTo3DTensor") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ComponentTo3DTensor") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractComponentFrom3DTensor") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorIOTest") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorMeanDiffusion") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorRadialDiffusion") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorAxialDiffusion") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorEigenvalue") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorColor") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorMask") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorToVector") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorToPhysicalSpace") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorToLocalSpace") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ValidTensor") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TensorToVectorComponent") == 0 )
        {
        TensorFunctions<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PASLQuantifyCBF") == 0 )
        {
        PASLQuantifyCBF<3>(argc, argv);
        }
      else if( strcmp( operation.c_str(), "Finite") == 0 )
        {
        Finite<3>( argc, argv );
        }
      else if( strcmp(operation.c_str(), "FillHoles") == 0 )
        {
        FillHoles<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "HistogramMatch") == 0 )
        {
        HistogramMatching<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PadImage") == 0 )
        {
        PadImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SigmoidImage") == 0 )
        {
        SigmoidImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SetOrGetPixel") == 0 )
        {
        SetOrGetPixel<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MakeImage") == 0 )
        {
        MakeImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "stack") == 0 )
        {
        StackImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CompareHeadersAndImages") == 0 )
        {
        CompareHeadersAndImages<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CountVoxelDifference") == 0 )
        {
        CountVoxelDifference<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageToFile") == 0 )
        {
        ConvertImageToFile<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PValueImage") == 0 )
        {
        PValueImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationUpdate") == 0 )
        {
        CorrelationUpdate<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PoissonDiffusion") == 0 )
        {
        PoissonDiffusion<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToMatrix") == 0 )
        {
        ConvertImageSetToMatrix<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "RandomlySampleImageSetToCSV") == 0 )
        {
        RandomlySampleImageSetToCSV<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToEigenvectors") == 0 )
        {
        ConvertImageSetToEigenvectors<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertVectorToImage") == 0 )
        {
        ConvertVectorToImage<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PropagateLabelsThroughMask") == 0 )
        {
        PropagateLabelsThroughMask<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FastMarchingSegmentation") == 0 )
        {
        FastMarchingSegmentation<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TriPlanarView") == 0 )
        {
        TriPlanarView<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TruncateImageIntensity") == 0 )
        {
        TruncateImageIntensity<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractSlice") == 0 )
        {
        ExtractSlice<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ClusterThresholdVariate") == 0 )
        {
        ClusterThresholdVariate<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertLandmarkFile") == 0 )
        {
        ConvertLandmarkFile<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PASL") == 0 )
        {
        PASL<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "pCASL") == 0 )
        {
        pCASL<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MTR") == 0 )
        {
        MTR<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MajorityVoting") == 0 )
        {
        MajorityVoting<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MostLikely") == 0 )
        {
        MostLikely<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationVoting") == 0 )
        {
        CorrelationVoting<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "STAPLE") == 0 )
        {
        STAPLE<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "AverageLabels") == 0 )
        {
        AverageLabels<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PearsonCorrelation") == 0 )
        {
        PearsonCorrelation<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NeighborhoodCorrelation") == 0 )
        {
        ImageMetrics<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NormalizedCorrelation") == 0 )
        {
        ImageMetrics<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Demons") == 0 )
        {
        ImageMetrics<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Mattes") == 0 )
        {
        ImageMetrics<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MinMaxMean") == 0 )
        {
        MinMaxMean<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "BlobDetector") == 0 )
        {
        BlobDetector<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Check3TissueLabeling") == 0 )
        {
        Check3TissueLabeling<3>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PureTissueN4WeightMask") == 0 )
        {
        PureTissueN4WeightMask<3>(argc, argv);
        }
      else
        {
        antscout << " cannot find operation : " << operation << std::endl;
        }
      }
      break;
    case 4:
      {
      if( strcmp(operation.c_str(), "m") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mresample") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "+") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "-") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vm") == 0 )
        {
        VImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vmresample") == 0 )
        {
        VImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v+") == 0 )
        {
        VImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "v-") == 0 )
        {
        VImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "/") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "^") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "exp") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "abs") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "addtozero") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "overadd") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "total") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "vtotal") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "mean") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Decision") == 0 )
        {
        ImageMath<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Neg") == 0 )
        {
        NegativeImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "G") == 0 )
        {
        SmoothImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MD") == 0 || strcmp(operation.c_str(), "ME") == 0 )
        {
        MorphImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MO") == 0 || strcmp(operation.c_str(), "MC") == 0 )
        {
        MorphImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GD") == 0 || strcmp(operation.c_str(), "GE") == 0 )
        {
        MorphImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GO") == 0 || strcmp(operation.c_str(), "GC") == 0 )
        {
        MorphImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "D") == 0 )
        {
        DistanceMap<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Normalize") == 0 )
        {
        NormalizeImage<4>(argc, argv);
        }

      else if( strcmp(operation.c_str(), "Grad") == 0 )
        {
        GradientImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Laplacian") == 0 )
        {
        LaplacianImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PH") == 0 )
        {
        PrintHeader<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "nvols") == 0 )
        {
        PrintHeader<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CenterImage2inImage1") == 0 )
        {
        CenterImage2inImage1<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Byte") == 0 )
        {
        ByteImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "LabelStats") == 0 )
        {
        LabelStats<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ROIStatistics") == 0 )
        {
        ROIStatistics<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "DiceAndMinDistSum") == 0 )
        {
        DiceAndMinDistSum<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Lipschitz") == 0 )
        {
        Lipschitz<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "InvId") == 0 )
        {
        InvId<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "GetLargestComponent") == 0 )
        {
        GetLargestComponent<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractVectorComponent") == 0 )
        {
        ExtractVectorComponent<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ThresholdAtMean") == 0 )
        {
        ThresholdAtMean<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FlattenImage") == 0 )
        {
        FlattenImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorruptImage") == 0 )
        {
        CorruptImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "RemoveLabelInterfaces") == 0 )
        {
        RemoveLabelInterfaces<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "EnumerateLabelInterfaces") == 0 )
        {
        EnumerateLabelInterfaces<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FitSphere") == 0 )
        {
        FitSphere<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Where") == 0 )
        {
        Where<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PureTissueN4WeightMask") == 0 )
        {
        PureTissueN4WeightMask<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "BlobDetector") == 0 )
        {
        BlobDetector<4>(argc, argv);
        }
      //    else if (strcmp(operation.c_str(),"TensorFA") == 0 )  TensorFunctions<4>(argc,argv);
      // else if (strcmp(operation.c_str(),"TensorIOTest") == 0 )  TensorFunctions<4>(argc,argv);
      // else if (strcmp(operation.c_str(),"TensorMeanDiffusion") == 0 )  TensorFunctions<4>(argc,argv);
      // else if (strcmp(operation.c_str(),"TensorColor") == 0) TensorFunctions<4>(argc,argv);
      // else if (strcmp(operation.c_str(),"TensorToVector") == 0) TensorFunctions<4>(argc,argv);
      // else if (strcmp(operation.c_str(),"TensorToVectorComponent")
      // == 0) TensorFunctions<4>(argc,argv);
      else if( strcmp( operation.c_str(), "Finite") == 0 )
        {
        Finite<4>( argc, argv );
        }
      else if( strcmp(operation.c_str(), "FillHoles") == 0 )
        {
        FillHoles<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "HistogramMatch") == 0 )
        {
        HistogramMatching<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PadImage") == 0 )
        {
        PadImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "SigmoidImage") == 0 )
        {
        SigmoidImage<4>(argc, argv);
        }
      //  else if (strcmp(operation.c_str(),"SetOrGetPixel") == 0 )  SetOrGetPixel<4>(argc,argv);
      else if( strcmp(operation.c_str(), "MakeImage") == 0 )
        {
        MakeImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "stack") == 0 )
        {
        StackImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CompareHeadersAndImages") == 0 )
        {
        CompareHeadersAndImages<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CountVoxelDifference") == 0 )
        {
        CountVoxelDifference<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageToFile") == 0 )
        {
        ConvertImageToFile<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PValueImage") == 0 )
        {
        PValueImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PoissonDiffusion") == 0 )
        {
        PoissonDiffusion<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationUpdate") == 0 )
        {
        CorrelationUpdate<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToMatrix") == 0 )
        {
        ConvertImageSetToMatrix<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "RandomlySampleImageSetToCSV") == 0 )
        {
        RandomlySampleImageSetToCSV<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertImageSetToEigenvectors") == 0 )
        {
        ConvertImageSetToEigenvectors<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertVectorToImage") == 0 )
        {
        ConvertVectorToImage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PropagateLabelsThroughMask") == 0 )
        {
        PropagateLabelsThroughMask<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "FastMarchingSegmentation") == 0 )
        {
        FastMarchingSegmentation<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TriPlanarView") == 0 )
        {
        TriPlanarView<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TruncateImageIntensity") == 0 )
        {
        TruncateImageIntensity<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ExtractSlice") == 0 )
        {
        ExtractSlice<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ClusterThresholdVariate") == 0 )
        {
        ClusterThresholdVariate<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ConvertLandmarkFile") == 0 )
        {
        ConvertLandmarkFile<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesDisassemble") == 0 )
        {
        TimeSeriesDisassemble<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesAssemble") == 0 )
        {
        TimeSeriesAssemble<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesSubset") == 0 )
        {
        TimeSeriesSubset<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesToMatrix") == 0 )
        {
        TimeSeriesToMatrix<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesSimpleSubtraction") == 0 )
        {
        TimeSeriesSimpleSubtraction<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "TimeSeriesInterpolationSubtraction") == 0 )
        {
        TimeSeriesInterpolationSubtraction<4>(argc, argv);
        }
      else if ( strcmp(operation.c_str(), "SplitAlternatingTimeSeries") == 0 )
        {
        SplitAlternatingTimeSeries<4>(argc, argv);
        }
      else if ( strcmp(operation.c_str(), "AverageOverDimension") == 0 )
        {
        AverageOverDimension<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ThreeTissueConfounds") == 0 )
        {
        ThreeTissueConfounds<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CompCorrAuto") == 0 )
        {
        CompCorrAuto<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "ComputeTimeSeriesLeverage") == 0 )
        {
        ComputeTimeSeriesLeverage<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PASLQuantifyCBF") == 0 )
        {
        PASLQuantifyCBF<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PCASLQuantifyCBF") == 0 )
        {
        //PCASLQuantifyCBF<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PASL") == 0 )
        {
        PASL<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "pCASL") == 0 )
        {
        pCASL<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MajorityVoting") == 0 )
        {
        MajorityVoting<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MostLikely") == 0 )
        {
        MostLikely<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "CorrelationVoting") == 0 )
        {
        CorrelationVoting<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "PearsonCorrelation") == 0 )
        {
        PearsonCorrelation<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NeighborhoodCorrelation") == 0 )
        {
        ImageMetrics<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "NormalizedCorrelation") == 0 )
        {
        ImageMetrics<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Demons") == 0 )
        {
        ImageMetrics<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "Mattes") == 0 )
        {
        ImageMetrics<4>(argc, argv);
        }
      else if( strcmp(operation.c_str(), "MinMaxMean") == 0 )
        {
        MinMaxMean<4>(argc, argv);
        }
      else
        {
        antscout << " cannot find 4D operation : " << operation << std::endl;
        }
      }
      break;

    default:
      antscout << " Dimension Not supported " << atoi(argv[1]) << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
