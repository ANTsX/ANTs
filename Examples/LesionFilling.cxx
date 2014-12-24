#include <iostream>
#include <ostream>
#include <vector>

#include "antsUtilities.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMaskImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
//LesionFilling dimension t1.nii.gz lesionmask output.nii.gz
namespace ants
{
template <unsigned int Dimension>
int LesionFilling( int argc, char * argv[] )
{
  typedef unsigned char PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ImageReaderType;
  typedef itk::ImageRegionIterator< LesionImageType> IteratorType;
  const int * Dimension = argv[1];
  const char * T1FileName = argv[2]
  const char * LesionMapFileName = argv[3];
  const char * OutputFileName = argv[4]

  typename ImageReaderType::Pointer LesionReader = ImageReaderType::New();
  LesionReader->SetFileName( LesionMapFileName );
  LesionReader->Update();

  typename ImageReaderType::Pointer T1Reader = ImageReaderType::New();
  T1Reader->SetFileName( T1FileName); 
  T1Reader->Update();
  
  typedef double           RealPixelType;  //  Operations
  typedef itk::Image<RealPixelType, Dimension> RealImageType;
  typedef itk::CastImageFilter< LesionImageType, RealImageType>
                                                         CastToRealFilterType;
  CastToRealFilterType::Pointer LesiontoReal = CastToRealFilterType::New();
  LesiontoReal->SetInput( LesionReader->GetOutput() );
  CastToRealFilterType::Pointer T1toReal = CastToRealFilterType::New();
  T1toReal->SetInput( T1Reader->GetOutput() );

  typedef itk::BinaryThresholdImageFilter <LesionImageType, LesionImageType>
                     BinaryThresholdImageFilterType;
  typedef itk::BinaryBallStructuringElement<
                               LesionImageType,
                               Dimension  >             StructuringElementType;
  typedef itk::BinaryDilateImageFilter<
                               LesionImageType,
                               LesionImageType,
                               StructuringElementType >  DilateFilterType;
  typedef itk::SubtractImageFilter <LesionImageType, LesionImageType>
                               SubtractImageFilterType;
  //finding connected components, we assume each component is one lesion
  typedef itk::ConnectedComponentImageFilter <LesionImageType, LesionImageType>
              ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connected =
              ConnectedComponentImageFilterType::New ();
  connected->SetInput( LesiontoReal->GetOutput() ) ;
  connected->Update();
  const int LesionNumber = connected->GetObjectCount() ;
  std::cout << "Number of lesions: " << LesionNumber << std::endl;
  for ( int i = 1; i < LesionNumber; i++)
  {
     std::vector<double> outervoxels;
     BinaryThresholdImageFilterType::Pointer thresholdFilter
                = BinaryThresholdImageFilterType::New();
     thresholdFilter->SetInput(connected->GetOutput());
     thresholdFilter->SetLowerThreshold( (double) i - 0.1);
     thresholdFilter->SetUpperThreshold( (double) i + 0.1);
     thresholdFilter->SetInsideValue  ( 1 );
     thresholdFilter->SetOutsideValue ( 0 );
     //Neighbouring voxel
     //filling lesions with the voxels surrounding them
     //first finding the edges of lesions
     //by subtracting dilated lesion map from lesion map itself
     DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
     structuringElement.SetRadius( 1 );  // 3x3 structuring element
     structuringElement.CreateStructuringElement();
     binaryDilate->SetKernel( structuringElement );
     binaryDilate->SetInput( thresholdFilter->GetOutput() );
     binaryDilate->SetDilateValue( 1 );
     // subtract dilated image form non-dilated one
     SubtractImageFilterType::Pointer subtractFilter
                   = SubtractImageFilterType::New ();
     //output = image1 - image2
     subtractFilter->SetInput1( binaryDilate->GetOutput() );
     subtractFilter->SetInput2( thresholdFilter->GetOutput() );
     subtractFilter->Update();
     //multiply the outer lesion mask with T1 to get only the neighbouring voxels
     typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
     MaskFilterType::Pointer maskFilter = MaskFilterType::New();
     maskFilter->SetInput(T1toReal->GetOutput() );
     maskFilter->SetMaskImage( subtractFilter->GetOutput() );
     //collecting non-zero voxels
     IteratorType it( maskFilter->GetOutput(),
                       maskFilter->GetOutput()->GetLargestPossibleRegion() );
     it.GoToBegin();
     /** Walk over the image. */
     while ( !it.IsAtEnd() )
       {
         if( it.Value() )
         {
           outervoxels.push_back ( it.Get() ); 
         }
         ++it;
       } // end while
     //calculating mean lesion intesity
     //Note: lesions should not be filled with values
     //less than their originial values, this is a
     //trick to exclude any CSF voxels in the outer mask (if any)
     MaskFilterType::Pointer maskFilterLesion = MaskFilterType::New();
     maskFilterLesion->SetInput( T1toReal->GetOutput() );
     maskFilterLesion->SetMaskImage( LesiontoReal->GetOutput() );
     IteratorType it( maskFilterLesion->GetOutput(),
                       maskFilterLesion->GetOutput()->GetLargestPossibleRegion() );
     it.GoToBegin();
     /** Walk over the image. */
     int counter  = 0;
     double meanInsideLesion = 0;
     while ( !it.IsAtEnd() )
       {
         if( it.Value() )
         {
           //coutning number of voxels inside lesion
           counter++;
           meanInsideLesion += it.Get();
         }
         ++it;
       }
     meanInsideLesion /= counter;
     //check that all outer voxels are more than the mean 
     //intensity of the lesion, i.e. not including CSF voxels
     
     IteratorType it( maskFilter->GetOutput(),
                       maskFilter->GetOutput()->GetLargestPossibleRegion() );
     it.GoToBegin();
     std::vector<double> outerWMVoxels;
     while ( !it.IsAtEnd() )
     {
       if ( it.Get() > meanInsideLesion )
       {
         outerWMVoxels.push_back( it.Get() );
       }//end if
       ++it;
    }//end while
    //walk through original T1
    //and chagne inside the lesion with a random pick from
    //collected normal appearing WM voxels (outerWMVoxels)
    IteratorType it( T1Reader->GetOutput(),
                     T1Reader->GetOutput()->GetLargestPossibleRegion() );
    IteratorType itL(thresholdFilter->GetOutput(),
                     thresholdFilter->GetOutput()->GetLargestPossibleRegion() );
    int max = outerWMVoxels.size();
    int min = 0;
    int index = min + (rand() % (int)(max - min + 1)) ;
    it.GoToBegin();
    itL.GoToBegin();
    while ( !it.IsAtEnd() )
    {
      if (itL.Get() == 1)
      {
        it.Set(outerWMVoxels[ index ];
      }
    }
  }//loop for lesions
}//main int
}//namespace ants
