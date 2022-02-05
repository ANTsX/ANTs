
#include "antsUtilities.h"
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

  using LinearInterpolatorType = itk::LinearInterpolateImageFunction<ImageType, RealType>;
  typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
  interpolator->SetInputImage(image);

  using NearestNeighborInterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, RealType>;
  typename NearestNeighborInterpolatorType::Pointer nn_interpolator = NearestNeighborInterpolatorType::New();
  nn_interpolator->SetInputImage(image);

  using BSplineInterpolatorType = itk::BSplineInterpolateImageFunction<ImageType, RealType>;
  typename BSplineInterpolatorType::Pointer bs_interpolator = BSplineInterpolatorType::New();
  bs_interpolator->SetInputImage(image);

  using GaussianInterpolatorType = itk::GaussianInterpolateImageFunction<ImageType, RealType>;
  typename GaussianInterpolatorType::Pointer g_interpolator = GaussianInterpolatorType::New();
  g_interpolator->SetInputImage(image);

  using HammingInterpolatorType = itk::WindowedSincInterpolateImageFunction<ImageType, 3>;
  typename HammingInterpolatorType::Pointer sh_interpolator = HammingInterpolatorType::New();
  sh_interpolator->SetInputImage(image);

  using Sinc1InterpolatorType =
    itk::WindowedSincInterpolateImageFunction<ImageType, 3, itk::Function::CosineWindowFunction<3>>;
  typename Sinc1InterpolatorType::Pointer sc_interpolator = Sinc1InterpolatorType::New();
  sc_interpolator->SetInputImage(image);

  using Sinc2InterpolatorType =
    itk::WindowedSincInterpolateImageFunction<ImageType, 3, itk::Function::WelchWindowFunction<3>>;
  typename Sinc2InterpolatorType::Pointer sw_interpolator = Sinc2InterpolatorType::New();
  sw_interpolator->SetInputImage(image);

  using Sinc3InterpolatorType =
    itk::WindowedSincInterpolateImageFunction<ImageType, 3, itk::Function::LanczosWindowFunction<3>>;
  typename Sinc3InterpolatorType::Pointer sl_interpolator = Sinc3InterpolatorType::New();
  sl_interpolator->SetInputImage(image);

  typename Sinc3InterpolatorType::Pointer sb_interpolator = Sinc3InterpolatorType::New();
  sb_interpolator->SetInputImage(image);

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

  char arg7 = '\0';
  if (argc > 7)
  {
    arg7 = *argv[7];
  }

  resampler->SetTransform(transform);
  resampler->SetInterpolator(interpolator);
  if (argc > 6 && std::stoi(argv[6]))
  {
    switch (std::stoi(argv[6]))
    {
      case 0:
      default:
      {
        resampler->SetInterpolator(interpolator);
      }
      break;
      case 1:
      {
        resampler->SetInterpolator(nn_interpolator);
      }
      break;
      case 2:
      {
        double sigma[ImageDimension];
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          sigma[d] = image->GetSpacing()[d];
        }
        double alpha = 1.0;

        if (argc > 7)
        {
          std::vector<RealType> sg = ConvertVector<RealType>(std::string(argv[7]));
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            sigma[d] = sg[d];
          }
        }
        if (argc > 8)
        {
          alpha = static_cast<double>(atof(argv[8]));
        }
        g_interpolator->SetParameters(sigma, alpha);

        resampler->SetInterpolator(g_interpolator);
      }
      break;
      case 3:
      {
        switch (arg7)
        {
          case 'h':
          default:
          {
            resampler->SetInterpolator(sh_interpolator);
          }
          break;
          case 'c':
          {
            resampler->SetInterpolator(sc_interpolator);
          }
          break;
          case 'l':
          {
            resampler->SetInterpolator(sl_interpolator);
          }
          break;
          case 'w':
          {
            resampler->SetInterpolator(sw_interpolator);
          }
          break;
          case 'b':
          {
            resampler->SetInterpolator(sb_interpolator);
          }
          break;
        }
        break;
      }
      case 4:
      {
        if (argc > 7 && std::stoi(argv[7]) >= 0 && std::stoi(argv[7]) <= 5)
        {
          bs_interpolator->SetSplineOrder(std::stoi(argv[7]));
        }
        else
        {
          bs_interpolator->SetSplineOrder(3);
        }
        resampler->SetInterpolator(bs_interpolator);
        break;
      }
    }
  }
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
              << "outputImage MxNxO [size=1,spacing=0] [interpolate type] [pixeltype]" << std::endl;
    std::cout << "  Interpolation type: " << std::endl;
    std::cout << "    0. linear (default)" << std::endl;
    std::cout << "    1. nn " << std::endl;
    std::cout << "    2. gaussian [sigma=imageSpacing] [alpha=1.0]" << std::endl;
    std::cout << "    3. windowedSinc [type = 'c'osine, 'w'elch, 'b'lackman, 'l'anczos, 'h'amming]" << std::endl;
    std::cout << "    4. B-Spline [order=3]" << std::endl;
    std::cout << " pixeltype  :  TYPE " << std::endl;
    std::cout << "  0  :  char   " << std::endl;
    std::cout << "  1  :  unsigned char   " << std::endl;
    std::cout << "  2  :  short   " << std::endl;
    std::cout << "  3  :  unsigned short   " << std::endl;
    std::cout << "  4  :  int   " << std::endl;
    std::cout << "  5  :  unsigned int   " << std::endl;
    std::cout << "  6  :  float (default)  " << std::endl;
    std::cout << "  7  :  double  " << std::endl;
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
        case 2:
        {
          return ResampleImage<2, char>(argc, argv);
        }
        break;
        case 3:
        {
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
        case 2:
        {
          return ResampleImage<2, unsigned char>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, unsigned char>(argc, argv);
        }
        break;
        case 4:
        {
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
        case 2:
        {
          return ResampleImage<2, short>(argc, argv);
        }
        break;
        case 3:
        {
          return ResampleImage<3, short>(argc, argv);
        }
        break;
        case 4:
        {
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
        case 2:
        {
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
