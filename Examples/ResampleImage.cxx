
#include "antsUtilities.h"
#include "antsCommandLineParser.h"
#include <algorithm>

#include <cstdio>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"

#include "itkConstantBoundaryCondition.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "ReadWriteData.h"
#include <string>
#include <vector>

namespace ants
{

namespace
{
enum class SincWindow
{
  Cosine,
  Welch,
  Blackman,
  Lanczos,
  Hamming
};

template <typename ImageType>
typename itk::InterpolateImageFunction<ImageType, double>::Pointer
MakeWindowedSincInterpolator(SincWindow window)
{
  constexpr unsigned int Radius = 3;
  using BaseInterpolatorType = itk::InterpolateImageFunction<ImageType, double>;
  typename BaseInterpolatorType::Pointer interpolator;
  switch (window)
  {
    case SincWindow::Cosine: {
      using InterpolatorType =
        itk::WindowedSincInterpolateImageFunction<ImageType, Radius, itk::Function::CosineWindowFunction<Radius>>;
      interpolator = InterpolatorType::New();
      break;
    }
    case SincWindow::Welch: {
      using InterpolatorType =
        itk::WindowedSincInterpolateImageFunction<ImageType, Radius, itk::Function::WelchWindowFunction<Radius>>;
      interpolator = InterpolatorType::New();
      break;
    }
    case SincWindow::Blackman: {
      using InterpolatorType =
        itk::WindowedSincInterpolateImageFunction<ImageType, Radius, itk::Function::BlackmanWindowFunction<Radius>>;
      interpolator = InterpolatorType::New();
      break;
    }
    case SincWindow::Lanczos: {
      using InterpolatorType =
        itk::WindowedSincInterpolateImageFunction<ImageType, Radius, itk::Function::LanczosWindowFunction<Radius>>;
      interpolator = InterpolatorType::New();
      break;
    }
    case SincWindow::Hamming:
    default: {
      using InterpolatorType = itk::WindowedSincInterpolateImageFunction<ImageType, Radius>;
      interpolator = InterpolatorType::New();
      break;
    }
  }
  return interpolator;
}
} // namespace

template <unsigned int ImageDimension, typename PixelType>
int
ResampleImage(int argc, char * argv[])
{
  using RealType = double;
  using ImageType = itk::Image<PixelType, ImageDimension>;

  typename ImageType::Pointer image = nullptr;
  ReadImage<ImageType>(image, argv[2]);

  using TransformType = itk::IdentityTransform<RealType, ImageDimension>;
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  using ResamplerType = itk::ResampleImageFilter<ImageType, ImageType, RealType>;
  typename ResamplerType::Pointer     resampler = ResamplerType::New();
  typename ResamplerType::SpacingType spacing;
  typename ResamplerType::SizeType    size;
  typename ImageType::IndexType       oldStartIndex = image->GetLargestPossibleRegion().GetIndex();
  typename ImageType::IndexType       newStartIndex;
  newStartIndex.Fill(0); // should be "same" as original start index but in new physical space

  std::vector<RealType> sp = ConvertVector<RealType>(std::string(argv[4]));

  if (argc <= 5 || std::stoi(argv[5]) == 0)
  {
    if (sp.size() == 1)
    {
      spacing.Fill(sp[0]);
    }
    else if (sp.size() == ImageDimension)
    {
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        spacing[d] = sp[d];
      }
    }
    else
    {
      std::cout << "Invalid spacing." << std::endl;
    }
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      RealType spacing_old = image->GetSpacing()[i];
      RealType size_old = image->GetLargestPossibleRegion().GetSize()[i];
      size[i] = static_cast<int>((spacing_old * size_old) / spacing[i] + 0.5);
      RealType oldstart = static_cast<float>(oldStartIndex[i]);
      newStartIndex[i] = static_cast<int>((spacing_old * oldstart) / spacing[i] + 0.5);
    }
  }
  else
  {
    if (sp.size() == 1)
    {
      size.Fill(static_cast<unsigned int>(sp[0]));
    }
    else if (sp.size() == ImageDimension)
    {
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        size[d] = static_cast<unsigned int>(sp[d]);
      }
    }
    else
    {
      std::cout << "Invalid size." << std::endl;
    }
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      auto     spacing_old = static_cast<RealType>(image->GetSpacing()[i]);
      auto     size_old = static_cast<RealType>(image->GetLargestPossibleRegion().GetSize()[i]);
      RealType ratio = (static_cast<RealType>(size_old) - itk::NumericTraits<RealType>::OneValue()) /
                       (static_cast<RealType>(size[i]) - itk::NumericTraits<RealType>::OneValue());
      spacing[i] = spacing_old * ratio;
      RealType oldstart = static_cast<float>(oldStartIndex[i]);
      newStartIndex[i] = static_cast<int>(oldstart * ratio + static_cast<RealType>(0.5));
    }
  }

  const std::string interpolationArgument = argc > 6 ? argv[6] : "Linear";
  using ParserType = itk::ants::CommandLineParser;
  auto parser = ParserType::New();
  auto interpolationOption = ParserType::OptionType::New();
  interpolationOption->AddFunction(interpolationArgument);
  const auto  interpolationFunction = interpolationOption->GetFunction();
  std::string interpolationName = interpolationFunction->GetName();
  ConvertToLowerCase(interpolationName);
  const auto  parameters = interpolationFunction->GetParameters();
  std::string interpolationError;

  // Numeric interpolation choices retain their default parameters. The next
  // positional argument remains the pixel type, including for Gaussian and BSpline.
  const std::vector<std::string> numericNames{ "linear", "nearestneighbor", "gaussian", "windowedsinc", "bspline" };
  if (interpolationName.size() == 1 && interpolationName[0] >= '0' && interpolationName[0] <= '4')
  {
    if (!parameters.empty())
    {
      std::cerr << "Deprecated numeric interpolation options do not accept parameters" << std::endl;
      return EXIT_FAILURE;
    }
    interpolationName = numericNames[interpolationName[0] - '0'];
  }

  using BaseInterpolatorType = itk::InterpolateImageFunction<ImageType, RealType>;
  typename BaseInterpolatorType::Pointer selectedInterpolator;

  try
  {
    if (interpolationName == "linear")
    {
      if (!parameters.empty())
      {
        interpolationError = "Linear does not accept parameters";
      }
      else
      {
        using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, RealType>;
        selectedInterpolator = InterpolatorType::New();
      }
    }
    else if (interpolationName == "nearestneighbor")
    {
      if (!parameters.empty())
      {
        interpolationError = "NearestNeighbor does not accept parameters";
      }
      else
      {
        using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, RealType>;
        selectedInterpolator = InterpolatorType::New();
      }
    }
    else if (interpolationName == "gaussian")
    {
      using InterpolatorType = itk::GaussianInterpolateImageFunction<ImageType, RealType>;

      double sigma[ImageDimension];
      for (unsigned int d = 0; d < ImageDimension; ++d)
      {
        sigma[d] = image->GetSpacing()[d];
      }
      double alpha = 1.0;

      if (parameters.size() > 2)
      {
        interpolationError = "Gaussian accepts at most two parameters: sigma and alpha";
      }
      else if (!parameters.empty())
      {
        std::string sigmaParameter = parameters[0];
        ConvertToLowerCase(sigmaParameter);
        if (sigmaParameter != "spacing")
        {
          const auto sigmaValues = parser->ConvertVector<RealType>(sigmaParameter);
          if (sigmaValues.size() == 1)
          {
            for (unsigned int d = 0; d < ImageDimension; ++d)
            {
              sigma[d] = sigmaValues[0];
            }
          }
          else if (sigmaValues.size() == ImageDimension)
          {
            for (unsigned int d = 0; d < ImageDimension; ++d)
            {
              sigma[d] = sigmaValues[d];
            }
          }
          else
          {
            interpolationError = "Gaussian sigma must be one value or one value per image dimension";
          }
        }
      }
      if (interpolationError.empty() && parameters.size() > 1)
      {
        alpha = parser->Convert<RealType>(parameters[1]);
      }

      if (interpolationError.empty())
      {
        typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetParameters(sigma, alpha);
        selectedInterpolator = interpolator;
      }
    }
    else if (interpolationName == "windowedsinc")
    {
      SincWindow window = SincWindow::Hamming;

      if (parameters.size() > 1)
      {
        interpolationError = "WindowedSinc accepts at most one parameter: type";
      }
      else if (!parameters.empty())
      {
        std::string windowName = parameters[0];
        ConvertToLowerCase(windowName);
        if (windowName == "cosine" || windowName == "c")
        {
          window = SincWindow::Cosine;
        }
        else if (windowName == "welch" || windowName == "w")
        {
          window = SincWindow::Welch;
        }
        else if (windowName == "blackman" || windowName == "b")
        {
          window = SincWindow::Blackman;
        }
        else if (windowName == "lanczos" || windowName == "l")
        {
          window = SincWindow::Lanczos;
        }
        else if (windowName == "hamming" || windowName == "h")
        {
          window = SincWindow::Hamming;
        }
        else
        {
          interpolationError = "unknown WindowedSinc type '" + parameters[0] + "'";
        }
      }

      if (interpolationError.empty())
      {
        selectedInterpolator = MakeWindowedSincInterpolator<ImageType>(window);
      }
    }
    else if (interpolationName == "bspline")
    {
      unsigned int order = 3;
      if (parameters.size() > 1)
      {
        interpolationError = "BSpline accepts at most one parameter: order";
      }
      else if (!parameters.empty())
      {
        order = parser->Convert<unsigned int>(parameters[0]);
        if (order > 5)
        {
          interpolationError = "BSpline order must be an integer from 0 through 5";
        }
      }

      if (interpolationError.empty())
      {
        using InterpolatorType = itk::BSplineInterpolateImageFunction<ImageType, RealType>;
        typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetSplineOrder(order);
        selectedInterpolator = interpolator;
      }
    }
    else
    {
      interpolationError = "unknown interpolator '" + interpolationName + "'";
    }
  }
  catch (const itk::ExceptionObject & error)
  {
    interpolationError = error.GetDescription();
  }

  if (!interpolationError.empty() || selectedInterpolator.IsNull())
  {
    std::cerr << "Invalid interpolation specification '" << interpolationArgument << "': " << interpolationError
              << std::endl;
    return EXIT_FAILURE;
  }

  selectedInterpolator->SetInputImage(image);
  resampler->SetTransform(transform);
  resampler->SetInterpolator(selectedInterpolator);
  resampler->SetInput(image);
  resampler->SetSize(size);
  resampler->SetOutputOrigin(image->GetOrigin());
  resampler->SetOutputDirection(image->GetDirection());
  resampler->SetOutputSpacing(spacing);
  //  resampler->SetOutputStartIndex( newStartIndex );
  resampler->SetDefaultPixelValue(0);
  resampler->Update();
  typename ImageType::Pointer outimage = resampler->GetOutput();
  //  typename ImageType::RegionType region = outimage->GetLargestPossibleRegion();
  //  region.SetIndex( newStartIndex );
  //  outimage->SetLargestPossibleRegion( region );
  ANTs::WriteImage<ImageType>(outimage, argv[3]);
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ResampleImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ResampleImage");

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

  if (argc < 5)
  {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
              << "outputImage MxNxO [size=1,spacing=0] [interpolation=Linear] [pixeltype=6]" << std::endl;
    std::cout << "  imageDimension: 2, 3, or 4" << std::endl;
    std::cout << "  inputImage: input image file name" << std::endl;
    std::cout << "  outputImage: output image file name" << std::endl;
    std::cout << "  MxNxO: either the new size (number of pixels) in each dimension, or the new spacing (physical size of a pixel) in each dimension" << std::endl;
    std::cout << "  size=1,spacing=0: Either 1 or 0, if the MxNxO argument is interpreted as size (1), or spacing (0)" << std::endl;
    std::cout << "  interpolation (parameters are positional; parameter names below document their defaults):"
              << std::endl;
    std::cout << "    Linear (default)" << std::endl;
    std::cout << "    NearestNeighbor" << std::endl;
    std::cout << "    Gaussian[<sigma=spacing>,<alpha=1>]" << std::endl;
    std::cout << "      sigma may be one value or one value per dimension separated by 'x'." << std::endl;
    std::cout << "    WindowedSinc[<type=hamming>] (fixed radius 3)" << std::endl;
    std::cout << "      type: cosine, welch, blackman, lanczos, or hamming." << std::endl;
    std::cout << "    BSpline[<order=3>]" << std::endl;
    std::cout << "    Deprecated numeric interpolation options (default parameters only):" << std::endl;
    std::cout << "      0: Linear,1: NearestNeighbor, 2: Gaussian, 3: WindowedSinc, 4: BSpline" << std::endl;
    std::cout << "  pixeltype: TYPE" << std::endl;
    std::cout << "  0  :  char   " << std::endl;
    std::cout << "  1  :  unsigned char   " << std::endl;
    std::cout << "  2  :  short   " << std::endl;
    std::cout << "  3  :  unsigned short   " << std::endl;
    std::cout << "  4  :  int   " << std::endl;
    std::cout << "  5  :  unsigned int   " << std::endl;
    std::cout << "  6  :  float (default)  " << std::endl;
    std::cout << "  7  :  double  " << std::endl;
    std::cout << "  Note: both input and output images will be processed as pixeltype, casting may cause " << std::endl;
    std::cout << "  rounding before resampling" << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  unsigned int typeoption = 6;
  if (argc > 7)
  {
    typeoption = std::stoi(argv[7]);
  }

  switch (typeoption)
  {
    case 0:
      switch (std::stoi(argv[1]))
      {
        case 2: {
          return ResampleImage<2, char>(argc, argv);
        }
        break;
        case 3: {
          return ResampleImage<3, char>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, char>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
      break;
    case 1:
      switch (std::stoi(argv[1]))
      {
        case 2: {
          return ResampleImage<2, unsigned char>(argc, argv);
        }
        break;
        case 3: {
          return ResampleImage<3, unsigned char>(argc, argv);
        }
        break;
        case 4: {
          return ResampleImage<4, unsigned char>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
      break;
    case 2:
      switch (std::stoi(argv[1]))
      {
        case 2: {
          return ResampleImage<2, short>(argc, argv);
        }
        break;
        case 3: {
          return ResampleImage<3, short>(argc, argv);
        }
        break;
        case 4: {
          return ResampleImage<4, short>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    case 3:
      switch (std::stoi(argv[1]))
      {
        case 2: {
          return ResampleImage<2, unsigned short>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, unsigned short>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, unsigned short>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    case 4:
      switch (std::stoi(argv[1]))
      {
        case 2:
        {
          return ResampleImage<2, int>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, int>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, int>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    case 5:
      switch (std::stoi(argv[1]))
      {
        case 2:
        {
          return ResampleImage<2, unsigned int>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, unsigned int>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, unsigned int>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    case 6:
      switch (std::stoi(argv[1]))
      {
        case 2:
        {
          return ResampleImage<2, float>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, float>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, float>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    case 7:
      switch (std::stoi(argv[1]))
      {
        case 2:
        {
          return ResampleImage<2, double>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, double>(argc, argv);
        }
        break;
        case 4:
        {
          return ResampleImage<4, double>(argc, argv);
        }
        break;
        default:
          std::cout << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    default:
      std::cout << "Unsupported pixel type" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // end of function
} // namespace ants
