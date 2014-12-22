#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkNeighborhoodIterator.h>


#include "antsUtilities.h"
#include <algorithm>



#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

namespace ants
{
template <unsigned int ImageDimension>
int LesionFilling( int argc, char * argv[] )
{
  typedef int                                   LesionType;
  typedef itk::Image<LesionType, ImageDimension> LesionImageType;
  typedef itk::ImageFileReader<LesionImageType>  LesionReaderType;

  typedef int                                   T1Type;
  typedef itk::Image<LesionType, ImageDimension> T1ImageType;
  typedef itk::ImageFileReader<LesionImageType>  T1ReaderType;
  
  typedef double           RealPixelType;  //  Operations

  const int * ImageDimension = argv[1];
  const char * method = argv[2];
  const char * LesionMapFileName = argv[3];
  if (method == "-neighbourVoxels")
  {
         const char * T1FileName = argv[4]
  }
  else
  {
         const char * WMFileName = argv[4]
         const char * T1FileName = argv[5]
  }


  typename LesionReaderType::Pointer LesionReader = LesionReaderType::New();
  LesionReader->SetFileName( LesionMapFileName );
  LesionReader->Update();
  
  if (method == "-neighbourVoxels")
  {

         //Neighbouring voxel
         //filling lesions with the voxels surrounding them
         //first finding the edges of lesions
         
         typedef itk::BinaryDilateImageFilter<
                                   InputImageType,
                                   OutputImageType,
                                   StructuringElementType >  DilateFilterType;
         
         DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
         
         binaryDilate->SetKernel( structuringElement );

         float variance = 1.0;
         float upperThreshold = 0.0;
         float lowerThreshold = 0.0;
         typedef itk::CastImageFilter< LesionImageType, RealImageType>
                                                                CastToRealFilterType;
         typedef itk::CannyEdgeDetectionImageFilter<RealImageType, RealImageType> CannyFilter;

         CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
         toReal->SetInput( LesionReader->GetOutput() );
         CannyFilter::Pointer cannyFilter = CannyFilter::New();
       
         cannyFilter->SetInput( toReal->GetOutput() );
         
         cannyFilter->SetVariance( variance );
         cannyFilter->SetUpperThreshold( upperThreshold );
         cannyFilter->SetLowerThreshold( lowerThreshold );
  
  /* LesionFilling d -neighbourVoxels binaryLeisonMap T1
   * LesionFilling d -mean binaryLesionMap AtroposWM T1
   * LesionFilling d -median binaryLesionmap AtroposWM T1
   * LesionFilling d -randomSample binaryLesionMap AtroposWM T1
   */
