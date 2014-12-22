#include <algorithm>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "antsUtilities.h"
#include "itkMaskImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

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
         typedef itk::BinaryBallStructuringElement<
                             InputPixelType,
                             Dimension  >             StructuringElementType;
         //Neighbouring voxel
         //filling lesions with the voxels surrounding them
         //first finding the edges of lesions
         //by subtracting dilated lesion map from lesion map itself
         typedef itk::BinaryDilateImageFilter<
                                   InputImageType,
                                   OutputImageType,
                                   StructuringElementType >  DilateFilterType;
         DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
         structuringElement.SetRadius( 1 );  // 3x3 structuring element
         structuringElement.CreateStructuringElement();
         binaryDilate->SetKernel( structuringElement );
         typedef itk::CastImageFilter< LesionImageType, RealImageType>
                                                                CastToRealFilterType;
         CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
         toReal->SetInput( LesionReader->GetOutput() );
         binaryDilate->SetInput( toReal->GetOutput() );
         binaryDilate->SetDilateValue( 1 );
         // subtract dilated image form non-dilated one
         typedef itk::SubtractImageFilter <ImageType, ImageType >
                     SubtractImageFilterType;
         SubtractImageFilterType::Pointer subtractFilter
                       = SubtractImageFilterType::New ();
         //output = image1 - image2
         subtractFilter->SetInput1(toReal->GetOutput() );
         subtractFilter->SetInput2(binaryDilate->GetOutput());
         subtractFilter->Update();
         //multiply the outer lesion mask with T1 to get only the neighbouring voxels
         typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
         MaskFilterType::Pointer maskFilter = MaskFilterType::New();
         maskFilter->SetInput( toReal->GetOutput() );
         maskFilter->SetMaskImage( subtractFilter->GetOutput() );
         //replacing lesions with extracted 
         
  }

  /* LesionFilling d -neighbourVoxels binaryLeisonMap T1
   * LesionFilling d -mean binaryLesionMap AtroposWM T1
   * LesionFilling d -median binaryLesionmap AtroposWM T1
   * LesionFilling d -randomSample binaryLesionMap AtroposWM T1
   */
