#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>

#include "antsUtilities.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
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
  template <unsigned int ImageDimension>
  int LesionFilling( int argc, char * argv[] )
  {
    const char * T1FileName = argv[2];
    const char * LesionMapFileName = argv[3];
    const char * OutputFileName = argv[4];
    if( argc < 3 )
      {
      std::cout << "Missing arguments, see usage." << std::endl;
      throw;
      }
    if( argc < 4 )
      {
      std::cout << "2 more arguments are necessary, see usage." << std::endl;
      throw;
      }
    if( argc < 5 )
      {
      std::cout << "Missing output filename." << std::endl;
      throw;
      }
    typedef itk::Image< double, ImageDimension> T1ImageType;
    typedef itk::Image< unsigned char, ImageDimension> LesionImageType;
    typedef itk::ImageFileReader<T1ImageType>  T1ImageReaderType;
    typedef itk::ImageFileReader<LesionImageType>  LesionImageReaderType;
    typename LesionImageReaderType::Pointer LesionReader = LesionImageReaderType::New();
    LesionReader->SetFileName( LesionMapFileName );
    try
      {
      LesionReader->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cout << "no lesion mask that can be read" << std::endl;
      return 0;
      }
    typename T1ImageReaderType::Pointer T1Reader = T1ImageReaderType::New();
    T1Reader->SetFileName( T1FileName); 
    try
      {
      T1Reader->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cout << "no T1 image that can be read" << std::endl;
      return 0;
      }
    typename T1ImageType::Pointer outImage = NULL ;
    outImage = T1Reader->GetOutput() ;
    typedef itk::ImageRegionIterator< T1ImageType> IteratorType;
    typedef itk::BinaryThresholdImageFilter <T1ImageType, T1ImageType>
                                 BinaryThresholdImageFilterType;
    typedef itk::BinaryBallStructuringElement<
                                 double,
                                 ImageDimension> StructuringElementType;
    typedef itk::BinaryDilateImageFilter<
                                 T1ImageType,
                                 T1ImageType,
                                 StructuringElementType >  DilateFilterType;
    typedef itk::ConnectedComponentImageFilter <LesionImageType, LesionImageType>
                ConnectedComponentFilterType;
    typename ConnectedComponentFilterType::Pointer connected =
                ConnectedComponentFilterType::New();
    connected->SetInput( LesionReader->GetOutput() );
    connected->Update();
    const int LesionNumber = connected->GetObjectCount() ;
    std::cout << "Number of lesions: " << LesionNumber << std::endl;
    for ( int i = 1;  i < LesionNumber + 1 ; i++)
    {
       typedef itk::CastImageFilter< LesionImageType, T1ImageType> FilterType;
       typename FilterType::Pointer filter = FilterType::New();
       filter->SetInput( connected->GetOutput() );
       filter->Update() ;
       typename BinaryThresholdImageFilterType::Pointer thresholdFilter
                  = BinaryThresholdImageFilterType::New();
       thresholdFilter->SetInput(filter->GetOutput());
       thresholdFilter->SetLowerThreshold( (double) i - 0.1);
       thresholdFilter->SetUpperThreshold( (double) i + 0.1);
       thresholdFilter->SetInsideValue  ( 1 );
       thresholdFilter->SetOutsideValue ( 0 );
       thresholdFilter->Update() ;
       //Neighbouring voxel
       //filling lesions with the voxels surrounding them
       //first finding the edges of lesions
       //by subtracting dilated lesion map from lesion map itself
       typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
  
       StructuringElementType structuringElement;
       structuringElement.SetRadius( 1 );  // 3x3 structuring element
       structuringElement.CreateStructuringElement();
       binaryDilate->SetKernel( structuringElement );
       binaryDilate->SetInput( thresholdFilter->GetOutput() );
       binaryDilate->SetDilateValue( 1 );
       binaryDilate->Update() ;
       // subtract dilated image form non-dilated one
       typedef itk::SubtractImageFilter <T1ImageType, T1ImageType>
                                    SubtractImageFilterType;
       typename SubtractImageFilterType::Pointer subtractFilter
                     = SubtractImageFilterType::New ();
       //output = image1 - image2
       subtractFilter->SetInput1( binaryDilate->GetOutput() );
       subtractFilter->SetInput2( thresholdFilter->GetOutput() );
       subtractFilter->Update();
       //multiply the outer lesion mask with T1 to get only the neighbouring voxels
       typedef itk::MaskImageFilter< T1ImageType, T1ImageType > MaskFilterType;
       typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
       maskFilter->SetInput( outImage ) ;
       maskFilter->SetMaskImage( subtractFilter->GetOutput() );
       maskFilter->Update() ;
       typename T1ImageType::Pointer LesionEdge= NULL ;
       LesionEdge = maskFilter->GetOutput() ;
       
       //calculating mean lesion intesity
       //Note: lesions should not be filled with values
       //less than their originial values, this is a
       //trick to exclude any CSF voxels in surronding voxels (if any)
       typename MaskFilterType::Pointer maskFilterLesion = MaskFilterType::New();
       maskFilterLesion->SetInput( T1Reader->GetOutput() );
       //do we need to pass a double lesion
       maskFilterLesion->SetMaskImage(thresholdFilter->GetOutput() );
       maskFilterLesion->Update() ;
       IteratorType it( maskFilterLesion->GetOutput(),
                         maskFilterLesion->GetOutput()->GetLargestPossibleRegion() );
       it.GoToBegin();
       /** Walk over the image. */
       int counter  = 0;
       double meanInsideLesion = 0;
       while ( !it.IsAtEnd() )
         {
           if( it.Value() != 0)
           {
             //coutning number of voxels inside lesion
             counter++;
             meanInsideLesion += it.Get();
           }
           ++it;
         }
       meanInsideLesion /= (double) counter;
       
       //check that all outer voxels are more than the mean 
       //intensity of the lesion, i.e. not including CSF voxels
       IteratorType itNoCSF( maskFilter->GetOutput(),
                         maskFilter->GetOutput()->GetLargestPossibleRegion() );
       itNoCSF.GoToBegin();
       std::vector<double> outerWMVoxels;
       while ( !itNoCSF.IsAtEnd() )
       {
         if ( itNoCSF.Get() >= meanInsideLesion )
         {
           outerWMVoxels.push_back( itNoCSF.Get() );
         }//end if
         ++itNoCSF;
      }//end while
      //walk through original T1
      //and change inside the lesion with a random pick from
      //collected normal appearing WM voxels (outerWMVoxels)
      IteratorType it4( outImage, 
                        outImage->GetLargestPossibleRegion() );
      IteratorType itL( thresholdFilter->GetOutput(),
                        thresholdFilter->GetOutput()->GetLargestPossibleRegion() );
      int max = outerWMVoxels.size();
      int min = 0;
      it4.GoToBegin();
      itL.GoToBegin();
      while ( !it4.IsAtEnd() )
      {
        if (itL.Get() == 1)
        {
          int index = min  + (rand() % (int)(max - min )) ;
          it4.Set( outerWMVoxels[ index ] );
        }//end if
        ++it4;
        ++itL;
      }//end while
    }//end of loop for lesions
  typedef itk::ImageFileWriter< T1ImageType>  WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( outImage);
  writer->SetFileName( OutputFileName );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    } 
  return EXIT_SUCCESS;
  }//main int
  
  int LesionFilling( std::vector<std::string> args, std::ostream* itkNotUsed( out_stream ) )
  {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert( args.begin(), "LesionFilling" );
  
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
  
    // antscout->set_stream( out_stream );
  
    //LesionFilling dimension t1.nii.gz lesionmask output.nii.gz
    if( argc < 3 )
      {
      std::cout << "Example usage: " << argv[0] << " imageDimension T1_image.nii.gz lesion_mask.nii.gz output_lesion_filled.nii.gz" << std::endl;
  
      if( argc >= 2 &&
          ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
        return EXIT_SUCCESS;
        }
      return EXIT_FAILURE;
      }
  
      switch( atoi( argv[1] ) )
        {
        case 2:
          {
          return LesionFilling<2>( argc, argv );
          }
          break;
        case 3:
          {
          return LesionFilling<3>( argc, argv );
          }
          break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
        }
     return EXIT_SUCCESS;    
  }//int LesionFilling std::vector
}//namespace ants
