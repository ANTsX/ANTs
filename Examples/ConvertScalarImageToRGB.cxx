
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
int
ConvertScalarImageToRGB(int argc, char * argv[])
{
  using RGBPixelType = itk::RGBPixel<unsigned char>;
  //  typedef itk::RGBAPixel<unsigned char> RGBPixelType;

  using RealType = float;

  using RealImageType = itk::Image<float, ImageDimension>;
  using RGBImageType = itk::Image<RGBPixelType, ImageDimension>;

  using ReaderType = itk::ImageFileReader<RealImageType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  reader->Update();

  using MaskImageType = itk::Image<unsigned char, ImageDimension>;
  typename MaskImageType::Pointer maskImage = nullptr;
  using MaskReaderType = itk::ImageFileReader<MaskImageType>;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName(argv[4]);
  try
  {
    maskreader->Update();
    maskImage = maskreader->GetOutput();
  }
  catch (...)
  {
    maskImage = nullptr;
  };

  std::string colormapString(argv[5]);

  using RGBFilterType = itk::ScalarToRGBColormapImageFilter<RealImageType, RGBImageType>;
  typename RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput(reader->GetOutput());

  if (colormapString == "red")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Red);
  }
  else if (colormapString == "green")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Green);
  }
  else if (colormapString == "blue")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Blue);
  }
  else if (colormapString == "grey")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Grey);
  }
  else if (colormapString == "cool")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Cool);
  }
  else if (colormapString == "hot")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Hot);
  }
  else if (colormapString == "spring")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Spring);
  }
  else if (colormapString == "autumn")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Autumn);
  }
  else if (colormapString == "winter")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Winter);
  }
  else if (colormapString == "copper")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Copper);
  }
  else if (colormapString == "summer")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Summer);
  }
  else if (colormapString == "jet")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::Jet);
    //    typedef itk::Function::JetColormapFunction<typename RealImageType::PixelType,
    //      typename RGBImageType::PixelType> ColormapType;
    //    typename ColormapType::Pointer colormap = ColormapType::New();
    //    rgbfilter->SetColormap( colormap );
  }
  else if (colormapString == "hsv")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::HSV);
    //    typedef itk::Function::HSVColormapFunction<typename RealImageType::PixelType,
    //      typename RGBImageType::PixelType> ColormapType;
    //    typename ColormapType::Pointer colormap = ColormapType::New();
    //    rgbfilter->SetColormap( colormap );
  }
  else if (colormapString == "overunder")
  {
    rgbfilter->SetColormap(RGBFilterType::RGBColormapFilterEnum::OverUnder);
  }
  else if (colormapString == "custom")
  {
    using ColormapType =
      itk::Function::CustomColormapFunction<typename RealImageType::PixelType, typename RGBImageType::PixelType>;
    typename ColormapType::Pointer colormap = ColormapType::New();

    std::ifstream str(argv[6]);
    std::string   line;

    // Get red values
    {
      std::getline(str, line);
      std::istringstream                 iss(line);
      float                              value;
      typename ColormapType::ChannelType channel;
      while (iss >> value)
      {
        channel.push_back(value);
      }

      colormap->SetRedChannel(channel);
    }

    // Get green values
    {
      std::getline(str, line);
      std::istringstream                 iss(line);
      float                              value;
      typename ColormapType::ChannelType channel;
      while (iss >> value)
      {
        channel.push_back(value);
      }

      colormap->SetGreenChannel(channel);
    }
    // Get blue values
    {
      std::getline(str, line);
      std::istringstream                 iss(line);
      float                              value;
      typename ColormapType::ChannelType channel;
      while (iss >> value)
      {
        channel.push_back(value);
      }

      colormap->SetBlueChannel(channel);
    }
    rgbfilter->SetColormap(colormap);
  }

  RealType minimumValue = itk::NumericTraits<RealType>::max();
  RealType maximumValue = itk::NumericTraits<RealType>::NonpositiveMin();

  itk::ImageRegionIteratorWithIndex<RealImageType> ItS(reader->GetOutput(),
                                                       reader->GetOutput()->GetLargestPossibleRegion());
  for (ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItS)
  {
    if (!maskImage || maskImage->GetPixel(ItS.GetIndex()) != 0)
    {
      if (minimumValue >= ItS.Get())
      {
        minimumValue = ItS.Get();
      }
      if (maximumValue <= ItS.Get())
      {
        maximumValue = ItS.Get();
      }
    }
  }

  rgbfilter->SetUseInputImageExtremaForScaling(false);
  rgbfilter->GetModifiableColormap()->SetMinimumInputValue(minimumValue);
  rgbfilter->GetModifiableColormap()->SetMaximumInputValue(maximumValue);

  rgbfilter->GetModifiableColormap()->SetMinimumRGBComponentValue(
    (argc > 9) ? static_cast<typename RGBPixelType::ComponentType>(atof(argv[9])) : 0);
  rgbfilter->GetModifiableColormap()->SetMaximumRGBComponentValue(
    (argc > 10) ? static_cast<typename RGBPixelType::ComponentType>(atof(argv[10])) : 255);

  if (argc > 7)
  {
    std::string argvString(argv[7]);
    if (argvString != "min" && argvString != "minimum")
    {
      rgbfilter->GetModifiableColormap()->SetMinimumInputValue(static_cast<RealType>(atof(argv[7])));
    }
  }

  if (argc > 8)
  {
    std::string argvString(argv[8]);
    if (argvString != "max" && argvString != "maximum")
    {
      rgbfilter->GetModifiableColormap()->SetMaximumInputValue(static_cast<RealType>(atof(argv[8])));
    }
  }

  try
  {
    rgbfilter->Update();
  }
  catch (...)
  {
    return EXIT_FAILURE;
  }

  if (maskImage)
  {
    itk::ImageRegionIterator<MaskImageType> ItM(maskImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<RGBImageType>  ItC(rgbfilter->GetOutput(),
                                               rgbfilter->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionIterator<RealImageType> ItS2(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

    ItM.GoToBegin();
    ItC.GoToBegin();
    ItS2.GoToBegin();

    while (!ItM.IsAtEnd())
    {
      if (ItM.Get() == 0)
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
        //        RealType ratio = ( ItS2.Get() - minimumValue ) / ( maximumValue - minimumValue );
        //
        //        rgbpixel.Fill( ratio * ( maximumRGBValue - minimumRGBValue )
        //          + minimumRGBValue );
        rgbpixel.Fill(itk::NumericTraits<typename RGBPixelType::ComponentType>::ZeroValue());

        ItC.Set(rgbpixel);
      }
      ++ItM;
      ++ItC;
      ++ItS2;
    }
  }

  using WriterType = itk::ImageFileWriter<RGBImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(rgbfilter->GetOutput());
  writer->SetFileName(argv[3]);
  writer->Update();


  if (argc > 11)
  {
    // Let's arbitrarily choose 256 colors to populate the look up table.

    std::ofstream str(argv[11]);

    RealType minimumValue2 = rgbfilter->GetModifiableColormap()->GetMinimumInputValue();
    RealType maximumValue2 = rgbfilter->GetModifiableColormap()->GetMaximumInputValue();

    RealType deltaValue = (maximumValue2 - minimumValue2) / 255.0f;

    for (unsigned int d = 0; d < 256; d++)
    {
      RealType value = minimumValue2 + d * deltaValue;

      RGBPixelType rgbPixel = rgbfilter->GetModifiableColormap()->operator()(value);

      str << value << "," << static_cast<int>(rgbPixel[0]) << "," << static_cast<int>(rgbPixel[1]) << ","
          << static_cast<int>(rgbPixel[2]) << ",1" << std::endl;
    }
    str.close();
  }


  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ConvertScalarImageToRGB(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ConvertScalarImageToRGB");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  // antscout->set_stream( out_stream );

  std::string usage =
    "  Usage: ConvertScalarImageToRGB imageDimension inputImage outputImage mask colormap [customColormapFile] [minimumInput=min] [maximumInput=max] [minimumRGBOutput=0] [maximumRGBOutput=255] [vtkLookupTable]\n"
    "\n"
    "  Arguments:\n"
    "    imageDimension\n"
    "        The dimension of the input image. Allowed values:\n"
    "           2   — for a 2D image\n"
    "           3   — for a 3D image\n"
    "\n"
    "    inputImage\n"
    "        Path to the input scalar image.\n"
    "\n"
    "    outputImage\n"
    "        Path to the output RGB image that will be generated.\n"
    "\n"
    "    mask\n"
    "        Path to a mask image. This sets pixels outside the mask to zero. Use 'none' for no mask.\n"
    "\n"
    "    colormap\n"
    "        A string specifying which colormap to use. Allowed values include:\n"
    "           grey, red, green, blue, copper, jet, hsv, spring, summer, autumn, winter, hot, cool, overunder, custom\n"
    "        If \"custom\" is chosen, the next argument (customColormapFile) must be provided.\n"
    "\n"
    "    customColormapFile\n"
    "        (Required if colormap is 'custom'.)\n"
    "        A text file defining the custom colormap. The colormap is a 3xN matrix, columns separated by spaces, where each\n"
    "        column defines an (R,G,B) color triplet in the range (0,1). The limits of the colormap are defined by the\n"
    "        minimumInput and maximumInput arguments. Intermediate values are interpolated linearly.\n"
    "\n"
    "    minimumInput (default = min)\n"
    "        Specifies the minimum scalar value for colormap scaling.\n"
    "        If this argument is not provided or is set to 'min' or 'minimum', the minimum from the input image is used\n"
    "        (considering only pixels within the mask if one is provided).\n"
    "\n"
    "    maximumInput (default = max)\n"
    "        Specifies the minimum scalar value for colormap scaling.\n"
    "        If this argument is not provided or is set to 'min' or 'minimum', the minimum from the input image is used\n"
    "        (considering only pixels within the mask if one is provided).\n"
    "\n"
    "    minimumRGBOutput (default = 0.0)\n"
    "        The minimum value for the RGB output components.\n"
    "\n"
    "    maximumRGBOutput (default = 255)\n"
    "        The maximum value for the RGB output components.\n"
    "\n"
    "    vtkLookupTable (Optional)\n"
    "        If provided, the program exports the color mapping to a VTK lookup table with 256 entries.\n"
    "        Each entry is formatted as:\n"
    "           scalarValue, red, green, blue, alpha\n"
    "        This maps scalar values between minimumInput and maximumInput to colors, as specified by the color map. This\n"
    "        is useful to export the mapping for use with programs that support VTK lookup tables.\n";

  if (argc < 6)
  {
    std::cout << usage << std::endl;

    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      ConvertScalarImageToRGB<2>(argc, argv);
    }
    break;
    case 3:
    {
      ConvertScalarImageToRGB<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
