#include "stdio.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "ReadWriteImage.h"
#include "TensorFunctions.h"

template <unsigned int ImageDimension>
int AverageTensorImages(unsigned int argc, char *argv[])
{
  typedef itk::Vector<float, 6>                        TensorType;
  typedef itk::Image<TensorType, ImageDimension>       ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

  char * outputName = argv[2];
  int    mathtype = atoi(argv[3]);
  float  numberofimages = (float)argc - 4.0;

  std::cout << "Averaging " << numberofimages << " images " << std::endl;

  typename ImageType::Pointer averageimage = NULL;
  typename ImageType::Pointer image2 = NULL;

  typename ImageType::SizeType size;
  size.Fill(0);
  unsigned int bigimage = 0;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    std::cout << " fn " << fn << std::endl;
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
        std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  std::cout << " largest image " << size << std::endl;

  bool logeuc = true;
  if( mathtype == 1 )
    {
    logeuc = false;
    }

  TensorType nullTensor;
  nullTensor[0] = nullTensor[1] = nullTensor[2] = nullTensor[3]
          = nullTensor[4] = nullTensor[5] = 0;

  ReadTensorImage<ImageType>(averageimage, argv[bigimage], logeuc);
  averageimage->FillBuffer(nullTensor);
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::string fn = std::string(argv[j]);
    ReadTensorImage<ImageType>(image2, fn.c_str(), logeuc);

    IteratorType vfIter( image2,  image2->GetLargestPossibleRegion() );
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
      TensorType val =  vfIter.Get() / numberofimages;
      averageimage->SetPixel(vfIter.GetIndex(), val + averageimage->GetPixel(vfIter.GetIndex() ) );
      }
    }

  WriteTensorImage<ImageType>(averageimage, outputName, logeuc);

  return EXIT_SUCCESS;
}

// Main Program
int main( int argc, char *argv[] )
{
  try
    {
    int dim = atoi(argv[1]);
    // char * outputName = argv[2];
    // int mathtype = atoi(argv[3]);
    int numberofimages = argc - 4;

    if( numberofimages < 1 )
      {
      std::cout << "Basic useage ex: " << std::endl;
      std::cout << argv[0] << " ImageDimension  average.nii mathtype list-of-files-via-wildcard " << std::endl;
      std::cout << " e.g. \n   AverageTensorImages 3  average.nii  1  *registered.nii " << std::endl;
      std::cout << " mathtype=[0=log-euclidean, 1=euclidean] " << std::endl;
      exit( EXIT_FAILURE );
      }

    // Get the image dimension
    switch( dim )
      {
      case 2:
        AverageTensorImages<2>(argc, argv);
        break;
      case 3:
        AverageTensorImages<3>(argc, argv);
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
      }

    return 0;;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
}
