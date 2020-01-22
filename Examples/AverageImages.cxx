/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

// Note: could easily add variance computation
// http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Variance.htm

#include "antsUtilities.h"

#include "itkArray.h"
#include "itkVariableLengthVector.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkOptimalSharpeningImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkResampleImageFilter.h"
#include <algorithm>

namespace ants
{
template <unsigned int ImageDimension, unsigned int NVectorComponents>
int AverageImages1(unsigned int argc, char *argv[])
{
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageFileReader<ImageType>              ImageFileReader;
  typedef itk::ImageFileWriter<ImageType>              writertype;

    {
    const std::string temp(argv[1]);
    if( !( ( temp == "2" ) || ( temp == "3" ) || ( temp == "4" ) ) )
      {
      std::cerr << "ERROR:  Dimension option must be 2 or 3 or 4, " << temp << "given" << std::endl;
      return EXIT_FAILURE;
      }
    }
    {
    const std::string temp(argv[3]);
    if( !( ( temp == "0" ) || ( temp == "1" )  ) )
      {
      std::cerr << "ERROR:  Normalize option must be 0 or 1, " << temp << "given" << std::endl;
      return EXIT_FAILURE;
      }
    }

  const bool  normalizei = std::stoi(argv[3]);
  const float numberofimages = static_cast<float>( argc ) - 4.0f;

  typename ImageType::SizeType maxSize;
  maxSize.Fill( 0 );
  unsigned int bigimage = 0;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    const std::string fn = std::string(argv[j]);
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName(fn.c_str() );
    imageIO->ReadImageInformation();

    for( unsigned int i = 0; i < ImageType::ImageDimension; i++ )
      {
      itk::SizeValueType currentDimensionSize = imageIO->GetDimensions( i );

      if( currentDimensionSize > maxSize[i] )
        {
        maxSize[i] = currentDimensionSize;
        bigimage = j;
        }
      }
    }
  std::cout << " bigimage " << bigimage << " maxSize " << maxSize << std::endl;

  typename ImageFileReader::Pointer reader = ImageFileReader::New();
  reader->SetFileName(argv[bigimage]);
  reader->Update();
  typename ImageType::Pointer averageimage = reader->GetOutput();
  std::cout << " Setting physcal space of output average image based on largest image " << std::endl;
  unsigned int vectorlength = reader->GetImageIO()->GetNumberOfComponents();
  std::cout << " Averaging " << numberofimages << " images with dim = " << ImageDimension << " vector components "
           << vectorlength << std::endl;
  PixelType meanval = 0;
  averageimage->FillBuffer(meanval);  // Reset all images to a mean of zero on the accumulator buffer.
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::cout << " reading " << std::string(argv[j]) << std::endl;
    typename ImageFileReader::Pointer rdr = ImageFileReader::New();
    rdr->SetFileName(argv[j]);
    rdr->Update();
    typedef itk::ResampleImageFilter<ImageType, ImageType, float> ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    // default to identity resampler->SetTransform( transform );
    // default to linearinterp resampler->SetInterpolator( interpolator );
    resampler->SetInput( rdr->GetOutput() );
    resampler->SetOutputParametersFromImage( averageimage );
    resampler->Update();

    typename ImageType::Pointer image2 = resampler->GetOutput();
    Iterator      vfIter2( image2,  image2->GetLargestPossibleRegion() );
    unsigned long ct = 0;
    if( normalizei )
      {
      meanval = 0;
      for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
        {
        const PixelType & localp = image2->GetPixel( vfIter2.GetIndex() );
        meanval = meanval + localp;
        ct++;
        }
      if( ct > 0 )
        {
        meanval = meanval / (float)ct;
        }
      if( meanval <= 0 )
        {
        meanval = (1);
        }
      }
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      PixelType val = vfIter2.Get();
      if( normalizei )
        {
        val /= meanval;
        }
      val = val / (float)numberofimages;
      const PixelType & oldval = averageimage->GetPixel(vfIter2.GetIndex() );
      averageimage->SetPixel(vfIter2.GetIndex(), val + oldval );
      }
    }

  //  typedef itk::OptimalSharpeningImageFilter<ImageType,ImageType > sharpeningFilter;
  typedef itk::LaplacianSharpeningImageFilter<ImageType, ImageType> sharpeningFilter;
  typename sharpeningFilter::Pointer shFilter = sharpeningFilter::New();
  if( normalizei && argc > 3 && vectorlength == 1 )
    {
    shFilter->SetInput( averageimage );
    //    shFilter->SetSValue(0.5);
    averageimage =  shFilter->GetOutput();
    }

  std::cout << " writing output ";
    {
    typename writertype::Pointer writer = writertype::New();
    writer->SetFileName(argv[2]);
    writer->SetInput( averageimage );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension, unsigned int NVectorComponents>
int AverageImages(unsigned int argc, char *argv[])
{
  typedef itk::Vector<float, NVectorComponents>        PixelType;
  typedef itk::Image<PixelType, ImageDimension>        ImageType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageFileReader<ImageType>              ImageFileReader;
  typedef itk::ImageFileWriter<ImageType>              writertype;

  //  bool  normalizei = std::stoi(argv[3]);
  float numberofimages = static_cast<float>( argc ) - 4.0f;
  typename ImageType::Pointer averageimage = nullptr;
  typename ImageType::Pointer image2 = nullptr;

  typename ImageType::SizeType size;
  size.Fill( 0 );
  typename ImageType::SizeType maxSize;
  maxSize.Fill( 0 );

  unsigned int bigimage = 4;
  for( unsigned int j = 4; j < argc; j++ )
    {
    // Get the image dimension
    std::string fn = std::string(argv[j]);
    std::cout << " fn " << fn << " " << ImageDimension << " " << NVectorComponents << std::endl;
    typename itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
    imageIO->SetFileName( fn.c_str() );
    imageIO->ReadImageInformation();

    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      size[i] = imageIO->GetDimensions( i );
      }

    for( unsigned int i = 0; i < imageIO->GetNumberOfDimensions(); i++ )
      {
      if( size[i] > maxSize[i] )
        {
        maxSize[i] = size[i];
        bigimage = j;
        std::cout << " bigimage " << j << " size " << size << std::endl;
        }
      }
    }

  std::cout << " largest image " << size << std::endl;
  typename ImageFileReader::Pointer reader = ImageFileReader::New();
  reader->SetFileName(argv[bigimage]);
  reader->Update();
  averageimage = reader->GetOutput();
  unsigned int vectorlength = reader->GetImageIO()->GetNumberOfComponents();
  std::cout << " Averaging " << numberofimages << " images with dim = " << ImageDimension << " vector components "
           << vectorlength << std::endl;
  typename ImageType::IndexType zindex; zindex.Fill(0);
  PixelType meanval = reader->GetOutput()->GetPixel(zindex);
  meanval.Fill(0);
  averageimage->FillBuffer(meanval);
  for( unsigned int j = 4; j < argc; j++ )
    {
    std::cout << " reading " << std::string(argv[j]) << " for average " << std::endl;
    typename ImageFileReader::Pointer rdr = ImageFileReader::New();
    rdr->SetFileName(argv[j]);
    rdr->Update();
    image2 = rdr->GetOutput();
    Iterator vfIter2( image2,  image2->GetLargestPossibleRegion() );
    for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
      {
      PixelType val = vfIter2.Get();
      double    valnorm = val.GetNorm();
      if( !std::isnan( valnorm  ) &&  !std::isinf( valnorm  )   )
        {
        val = val / (float)numberofimages;
        PixelType oldval = averageimage->GetPixel( vfIter2.GetIndex() );
        averageimage->SetPixel(vfIter2.GetIndex(), val + oldval );
        }
      }
    }

    {
    typename writertype::Pointer writer = writertype::New();
    writer->SetFileName(argv[2]);
    writer->SetInput( averageimage );
    writer->Update();
    }
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int AverageImages( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "AverageImages" );
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
  argv[argc] = nullptr;
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

  if( argc < 5 )
    {
    std::cout << "\n" << std::endl;
    std::cout << "Usage: \n" << std::endl;
    std::cout << argv[0] << " ImageDimension Outputfname.nii.gz Normalize <images> \n" << std::endl;
    std::cout << " Compulsory arguments: \n" << std::endl;
    std::cout << " ImageDimension: 2 or 3 (for 2 or 3 dimensional input).\n " << std::endl;
    std::cout << " Outputfname.nii.gz: the name of the resulting image.\n" << std::endl;
    std::cout
      <<
      " Normalize: 0 (false) or 1 (true); if true, the 2nd image is divided by its mean. This will select the largest image to average into.\n"
      << std::endl;
    std::cout << " Example Usage:\n" << std::endl;
    std::cout << argv[0] << " 3 average.nii.gz  1  *.nii.gz \n" << std::endl;
    std::cout << " \n" << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  const int                 dim = std::stoi( argv[1] );
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(argv[4], itk::ImageIOFactory::FileModeEnum::ReadMode);
  imageIO->SetFileName(argv[4]);
  imageIO->ReadImageInformation();
  unsigned int ncomponents = imageIO->GetNumberOfComponents();

  // Get the image dimension
  switch( dim )
    {
    case 2:
      {
      switch( ncomponents )
        {
        case 2:
          {
          return AverageImages<2, 2>(argc, argv);
          }
          break;
        default:
          {
          return AverageImages1<2, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 3:
      {
      switch( ncomponents )
        {
        case 7:
          {
          return AverageImages<3, 7>(argc, argv);
          }
          break;
        case 6:
          {
          return AverageImages<3, 6>(argc, argv);
          }
          break;
        case 3:
          {
          return AverageImages<3, 3>(argc, argv);
          }
          break;
        case 2:
          {
          return AverageImages<3, 2>(argc, argv);
          }
          break;
        default:
          {
          return AverageImages1<3, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    case 4:
      {
      switch( ncomponents )
        {
        case 7:
          {
          return AverageImages<4, 7>(argc, argv);
          }
          break;
        case 6:
          {
          return AverageImages<4, 6>(argc, argv);
          }
          break;
        case 4:
          {
          return AverageImages<4, 4>(argc, argv);
          }
          break;
        case 3:
          {
          return AverageImages<4, 3>(argc, argv);
          }
          break;
        default:
          {
          return AverageImages1<4, 1>(argc, argv);
          }
          break;
        }
      }
      break;
    default:
      std::cout << " You passed ImageDimension: " << dim << " . Please use only image domains of 2, 3 or 4  "
               << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
