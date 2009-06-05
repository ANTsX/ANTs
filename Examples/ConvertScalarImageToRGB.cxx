#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"

#include "itkRedColormapFunctor.h"
#include "itkGreenColormapFunctor.h"
#include "itkBlueColormapFunctor.h"
#include "itkGreyColormapFunctor.h"
#include "itkHotColormapFunctor.h"
#include "itkCoolColormapFunctor.h"
#include "itkSpringColormapFunctor.h"
#include "itkSummerColormapFunctor.h"
#include "itkAutumnColormapFunctor.h"
#include "itkWinterColormapFunctor.h"
#include "itkCopperColormapFunctor.h"
#include "itkHSVColormapFunctor.h"
#include "itkJetColormapFunctor.h"
#include "itkCustomColormapFunctor.h"
#include "itkOverUnderColormapFunctor.h"

#include "itkScalarToRGBColormapImageFilter.h"

#include <fstream.h>
#include <sstream>
#include <string>

template <unsigned int ImageDimension>
int ConvertScalarImageToRGB( int argc, char *argv[] )
{
  typedef unsigned int                 PixelType;
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
//  typedef itk::RGBAPixel<unsigned char> RGBPixelType;

  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension>    ImageType;
  typedef itk::Image<float, ImageDimension>        RealImageType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  std::string colormapString( argv[5] );

  typedef itk::ScalarToRGBColormapImageFilter<RealImageType,
                                              RGBImageType> RGBFilterType;
  typename RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput( reader->GetOutput() );

  if( colormapString == "red" )
    {
    rgbfilter->SetColormap( RGBFilterType::Red );
    }
  else if( colormapString == "green"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Green );
    }
  else if( colormapString == "blue"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Blue );
    }
  else if( colormapString == "grey"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Grey );
    }
  else if( colormapString == "cool"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Cool );
    }
  else if( colormapString == "hot"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Hot );
    }
  else if( colormapString == "spring"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Spring );
    }
  else if( colormapString == "autumn"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Autumn );
    }
  else if( colormapString == "winter"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Winter );
    }
  else if( colormapString == "copper"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Copper );
    }
  else if( colormapString == "summer"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Summer );
    }
  else if( colormapString == "jet"  )
    {
//    rgbfilter->SetColormap( RGBFilterType::Jet );
    typedef itk::Functor::JetColormapFunctor<typename RealImageType::PixelType,
                                             typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if( colormapString == "hsv"  )
    {
//    rgbfilter->SetColormap( RGBFilterType::HSV );
    typedef itk::Functor::HSVColormapFunctor<typename RealImageType::PixelType,
                                             typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if( colormapString == "overunder"  )
    {
    rgbfilter->SetColormap( RGBFilterType::OverUnder );
    }
  else if( colormapString == "custom"  )
    {
    typedef itk::Functor::CustomColormapFunctor<typename RealImageType::PixelType,
                                                typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();

    ifstream    str( argv[6] );
    std::string line;

    // Get red values
      {
      std::getline( str, line );
      std::istringstream iss( line );
      float              value;
      typename ColormapType::ChannelType channel;
      while( iss >> value )
        {
        channel.push_back( value );
        }

      colormap->SetRedChannel( channel );
      }

    // Get green values
      {
      std::getline( str, line );
      std::istringstream iss( line );
      float              value;
      typename ColormapType::ChannelType channel;
      while( iss >> value )
        {
        channel.push_back( value );
        }

      colormap->SetGreenChannel( channel );
      }
    // Get blue values
      {
      std::getline( str, line );
      std::istringstream iss( line );
      float              value;
      typename ColormapType::ChannelType channel;
      while( iss >> value )
        {
        channel.push_back( value );
        }

      colormap->SetBlueChannel( channel );
      }
    rgbfilter->SetColormap( colormap );
    }

  rgbfilter->GetColormap()->SetMinimumRGBComponentValue( 0 );
  rgbfilter->GetColormap()->SetMaximumRGBComponentValue( 255 );

  try
    {
    rgbfilter->Update();
    }
  catch( ... )
    {
    return EXIT_FAILURE;
    }

  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  typedef itk::ImageFileReader<MaskImageType>       MaskReaderType;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[4] );
  try
    {
    maskreader->Update();

    itk::ImageRegionIterator<MaskImageType> ItM( maskreader->GetOutput(),
                                                 maskreader->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RGBImageType> ItC( rgbfilter->GetOutput(),
                                                rgbfilter->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RealImageType> ItS( reader->GetOutput(),
                                                 reader->GetOutput()->GetLargestPossibleRegion() );

    ItM.GoToBegin();
    ItC.GoToBegin();
    ItS.GoToBegin();

    while( !ItM.IsAtEnd() )
      {
      if( ItM.Get() == 0 )
        {
        RGBPixelType rgbpixel;

        RealType minimumValue = rgbfilter->GetColormap()->GetMinimumInputValue();
        RealType maximumValue = rgbfilter->GetColormap()->GetMaximumInputValue();

        RealType minimumRGBValue
          = rgbfilter->GetColormap()->GetMinimumRGBComponentValue();
        RealType maximumRGBValue
          = rgbfilter->GetColormap()->GetMaximumRGBComponentValue();

        RealType ratio = ( ItS.Get() - minimumValue ) / ( maximumValue - minimumValue );

        rgbpixel.Fill( ratio * ( maximumRGBValue - minimumRGBValue )
                       + minimumRGBValue );

        ItC.Set( rgbpixel );
        }
      ++ItM;
      ++ItC;
      ++ItS;
      }
    }
  catch( ... )
    {
    }

  typedef itk::ImageFileWriter<RGBImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rgbfilter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage mask colormap [customColormapFile]"
              << std::endl;
    std::cout << "  Possible colormaps: grey, red, green, blue, copper, jet, hsv, ";
    std::cout << "spring, summer, autumn, winter, hot, cool, overunder, custom" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      ConvertScalarImageToRGB<2>( argc, argv );
      break;
    case 3:
      ConvertScalarImageToRGB<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
