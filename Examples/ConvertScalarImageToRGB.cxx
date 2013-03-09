
#include "antsUtilities.h"
#include <algorithm>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"

#include "itkRedColormapFunction.h"
#include "itkGreenColormapFunction.h"
#include "itkBlueColormapFunction.h"
#include "itkGreyColormapFunction.h"
#include "itkHotColormapFunction.h"
#include "itkCoolColormapFunction.h"
#include "itkSpringColormapFunction.h"
#include "itkSummerColormapFunction.h"
#include "itkAutumnColormapFunction.h"
#include "itkWinterColormapFunction.h"
#include "itkCopperColormapFunction.h"
#include "itkHSVColormapFunction.h"
#include "itkJetColormapFunction.h"
#include "itkCustomColormapFunction.h"
#include "itkOverUnderColormapFunction.h"

#include "itkScalarToRGBColormapImageFilter.h"

namespace ants
{
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

  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[4] );
  try
    {
    maskreader->Update();
    maskImage = maskreader->GetOutput();
    }
  catch( ... )
    {
    maskImage = NULL;
    }
  ;

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
    rgbfilter->SetColormap( RGBFilterType::Jet );
//    typedef itk::Function::JetColormapFunction<typename RealImageType::PixelType,
//      typename RGBImageType::PixelType> ColormapType;
//    typename ColormapType::Pointer colormap = ColormapType::New();
//    rgbfilter->SetColormap( colormap );
    }
  else if( colormapString == "hsv"  )
    {
    rgbfilter->SetColormap( RGBFilterType::HSV );
//    typedef itk::Function::HSVColormapFunction<typename RealImageType::PixelType,
//      typename RGBImageType::PixelType> ColormapType;
//    typename ColormapType::Pointer colormap = ColormapType::New();
//    rgbfilter->SetColormap( colormap );
    }
  else if( colormapString == "overunder"  )
    {
    rgbfilter->SetColormap( RGBFilterType::OverUnder );
    }
  else if( colormapString == "custom"  )
    {
    typedef itk::Function::CustomColormapFunction<typename RealImageType::PixelType,
                                                  typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();

    std::ifstream str( argv[6] );
    std::string   line;

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
//    rgbfilter->SetColormap( colormap );
    }

  if( maskImage )
    {
    RealType maskMinimumValue = itk::NumericTraits<RealType>::max();
    RealType maskMaximumValue = itk::NumericTraits<RealType>::NonpositiveMin();

    itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
                                                 maskImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RealImageType> ItS( reader->GetOutput(),
                                                 reader->GetOutput()->GetLargestPossibleRegion() );
    for( ItM.GoToBegin(), ItS.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItS )
      {
      if( ItM.Get() != 0 )
        {
        if( maskMinimumValue >= ItS.Get() )
          {
          maskMinimumValue = ItS.Get();
          }
        if( maskMaximumValue <= ItS.Get() )
          {
          maskMaximumValue = ItS.Get();
          }
        }
      }

    rgbfilter->SetUseInputImageExtremaForScaling( false );
    rgbfilter->GetModifiableColormap()->SetMinimumInputValue( maskMinimumValue );
    rgbfilter->GetModifiableColormap()->SetMaximumInputValue( maskMaximumValue );
    }

  rgbfilter->GetModifiableColormap()->SetMinimumRGBComponentValue(
    ( argc > 9 ) ? static_cast<
      typename RGBPixelType::ComponentType>( atof( argv[9] ) ) : 0 );
  rgbfilter->GetModifiableColormap()->SetMaximumRGBComponentValue(
    ( argc > 10 ) ? static_cast<
      typename RGBPixelType::ComponentType>( atof( argv[10] ) ) : 255 );

  if( argc > 8 )
    {
    rgbfilter->SetUseInputImageExtremaForScaling( false );
    rgbfilter->GetModifiableColormap()->SetMinimumInputValue(
      static_cast<RealType>( atof( argv[7] ) ) );
    rgbfilter->GetModifiableColormap()->SetMaximumInputValue(
      static_cast<RealType>( atof( argv[8] ) ) );
    }

  try
    {
    rgbfilter->Update();
    }
  catch( ... )
    {
    return EXIT_FAILURE;
    }

  if( maskImage )
    {
    itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
                                                 maskImage->GetLargestPossibleRegion() );
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

//        RealType minimumValue = rgbfilter->GetModifiableColormap()->GetMinimumInputValue();
//        RealType maximumValue = rgbfilter->GetModifiableColormap()->GetMaximumInputValue();
//
//        RealType minimumRGBValue
//          = rgbfilter->GetModifiableColormap()->GetMinimumRGBComponentValue();
//        RealType maximumRGBValue
//          = rgbfilter->GetModifiableColormap()->GetMaximumRGBComponentValue();
//
//        RealType ratio = ( ItS.Get() - minimumValue ) / ( maximumValue - minimumValue );
//
//        rgbpixel.Fill( ratio * ( maximumRGBValue - minimumRGBValue )
//          + minimumRGBValue );
        rgbpixel.Fill( itk::NumericTraits<typename RGBPixelType::ComponentType>::Zero );

        ItC.Set( rgbpixel );
        }
      ++ItM;
      ++ItC;
      ++ItS;
      }
    }

  typedef itk::ImageFileWriter<RGBImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rgbfilter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ConvertScalarImageToRGB( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ConvertScalarImageToRGB" );
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

  if( argc < 6 )
    {
    antscout << "Usage: " << argv[0] << " imageDimension inputImage outputImage "
             << "mask colormap [customColormapFile] [minimumInput] [maximumInput] "
             << "[minimumRGBOutput] [maximumRGBOutput]" << std::endl;
    antscout << "  Possible colormaps: grey, red, green, blue, copper, jet, hsv, ";
    antscout << "spring, summer, autumn, winter, hot, cool, overunder, custom" << std::endl;
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
      ConvertScalarImageToRGB<2>( argc, argv );
      }
      break;
    case 3:
      {
      ConvertScalarImageToRGB<3>( argc, argv );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
