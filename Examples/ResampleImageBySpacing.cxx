/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "antsUtilities.h"
#include <algorithm>
#include "itkImage.h"
#include "ReadWriteData.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ResampleImageBySpacing(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  args.insert(args.begin(), "ResampleImageBySpacing");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
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
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0]
              << "  ImageDimension inputImageFile  outputImageFile outxspc outyspc {outzspacing}  {dosmooth?}  "
                 "{addvox} {nn-interp?}"
              << std::endl;
    std::cout << " addvox pads each dimension by addvox " << std::endl;
    std::cout << "  " << std::endl;
    //    std::cout << " interp 0 = linear, 1 = nn " << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  unsigned int Dimension = std::stoi(argv[1]);

  if (Dimension == 2)
  {
    using InputPixelType = float;
    using InternalPixelType = float;
    using OutputPixelType = float;

    using InputImageType = itk::Image<InputPixelType, 2>;
    using InternalImageType = itk::Image<InternalPixelType, 2>;
    using OutputImageType = itk::Image<OutputPixelType, 2>;

    InputImageType::Pointer inputImage;
    ReadImage<InputImageType>(inputImage, argv[2]);

    const InputImageType::SpacingType & inputSpacing = inputImage->GetSpacing();

    OutputImageType::SpacingType spacing;
    for (int i = 0; i < 2; i++)
    {
      spacing[i] = inputSpacing[i];
    }

    std::cout << " spacing " << spacing << " dim " << 2 << std::endl;

    bool dosmooth = true;
    if (argc > 4)
    {
      spacing[0] = atof(argv[4]);
    }
    if (argc > 5)
    {
      spacing[1] = atof(argv[5]);
    }
    if (argc > 6)
    {
      dosmooth = std::stoi(argv[6]);
    }
    int addvox = 0;
    if (argc > 7)
    {
      addvox = std::stoi(argv[7]);
    }
    bool nn = false;
    if (argc > 8)
    {
      nn = std::stoi(argv[7]);
    }

    std::cout << " spacing2 " << spacing << std::endl;

    InternalImageType::Pointer smoothedImage = inputImage;
    if (dosmooth)
    {
      for (int sm = 0; sm < 2; sm++)
      {
        using GaussianFilterType = itk::RecursiveGaussianImageFilter<OutputImageType, OutputImageType>;

        GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
        smootherX->SetInput(smoothedImage);
        const float sig = atof(argv[4 + sm]) / inputSpacing[sm] - 1.0;
        std::cout << " smoothing by : " << sig << " dir " << sm << std::endl;
        smootherX->SetSigma(sig);
        smootherX->SetDirection(sm);
        smootherX->SetNormalizeAcrossScale(false);
        if (sig > 0 && dosmooth)
        {
          try
          {
            smootherX->Update();
          }
          catch (const itk::ExceptionObject & excep)
          {
            std::cout << "Exception catched !" << std::endl;
            std::cout << excep << std::endl;
          }
          smoothedImage = smootherX->GetOutput();
        }
      }
    }

    // InternalImageType::ConstPointer smoothedImage = smootherY->GetOutput();

    // InternalImageType::ConstPointer smoothedImage = reader->GetOutput();
    // smoothedImage =SmoothImage<ImageType>(reader->GetOutput() , );

    using ResampleFilterType = itk::ResampleImageFilter<InternalImageType, OutputImageType>;

    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    using TransformType = itk::IdentityTransform<double, 2>;

    using InterpolatorType = itk::LinearInterpolateImageFunction<InternalImageType, double>;
    using InterpolatorType2 = itk::NearestNeighborInterpolateImageFunction<InternalImageType, double>;

    InterpolatorType::Pointer  interpolator = InterpolatorType::New();
    InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();

    resampler->SetInterpolator(interpolator);
    if (nn == 1)
    {
      resampler->SetInterpolator(interpolator2);
    }

    InternalImageType::IndexType ind;
    ind.Fill(1);
    resampler->SetDefaultPixelValue(inputImage->GetPixel(ind)); // zero regions without source

    // Use the inputImage as initial template
    resampler->SetOutputParametersFromImage(inputImage);
    // Reset spacing by explicit specification
    std::cout << " out space " << spacing << std::endl;
    resampler->SetOutputSpacing(spacing);

    InputImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
    using SizeValueType = InputImageType::SizeType::SizeValueType;
    InputImageType::SizeType size;
    for (int i = 0; i < 2; i++)
    {
      size[i] = static_cast<SizeValueType>(inputSize[i] * inputSpacing[i] / spacing[i] + addvox);
    }

    std::cout << " output size " << size << " spc " << spacing << std::endl;
    resampler->SetSize(size);
    resampler->SetInput(smoothedImage);
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampler->SetTransform(transform);
    resampler->Update();
    ANTs::WriteImage<OutputImageType>(resampler->GetOutput(), argv[3]);
  }

  if (Dimension == 3)
  {
    using InputPixelType = float;
    using InternalPixelType = float;
    using OutputPixelType = float;

    using InputImageType = itk::Image<InputPixelType, 3>;
    using InternalImageType = itk::Image<InternalPixelType, 3>;
    using OutputImageType = itk::Image<OutputPixelType, 3>;
    InputImageType::Pointer inputImage;
    ReadImage<InputImageType>(inputImage, argv[2]);

    const InputImageType::SpacingType & inputSpacing = inputImage->GetSpacing();

    OutputImageType::SpacingType spacing;
    for (int i = 0; i < 3; i++)
    {
      spacing[i] = inputSpacing[i];
    }

    std::cout << " spacing " << spacing << " dim " << 3 << std::endl;

    bool dosmooth = true;
    if (argc > 4)
    {
      spacing[0] = atof(argv[4]);
    }
    if (argc > 5)
    {
      spacing[1] = atof(argv[5]);
    }
    if (argc > 6)
    {
      spacing[2] = atof(argv[6]);
    }
    if (argc > 7)
    {
      dosmooth = std::stoi(argv[7]);
    }
    int addvox = 0;
    if (argc > 8)
    {
      addvox = std::stoi(argv[8]);
    }
    bool nn = false;
    if (argc > 9)
    {
      nn = std::stoi(argv[9]);
    }

    std::cout << " spacing2 " << spacing << std::endl;

    InternalImageType::Pointer smoothedImage = inputImage;
    if (dosmooth)
    {
      for (int sm = 0; sm < 3; sm++)
      {
        using GaussianFilterType = itk::RecursiveGaussianImageFilter<OutputImageType, OutputImageType>;

        GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
        smootherX->SetInput(smoothedImage);
        const float sig = atof(argv[4 + sm]) / inputSpacing[sm] - 1.0;
        std::cout << " smoothing by : " << sig << " dir " << sm << std::endl;
        smootherX->SetSigma(sig);
        smootherX->SetDirection(sm);
        smootherX->SetNormalizeAcrossScale(false);
        if (sig > 0 && dosmooth)
        {
          try
          {
            smootherX->Update();
          }
          catch (const itk::ExceptionObject & excep)
          {
            std::cout << "Exception catched !" << std::endl;
            std::cout << excep << std::endl;
          }
          smoothedImage = smootherX->GetOutput();
        }
      }
    }

    // InternalImageType::ConstPointer smoothedImage = smootherY->GetOutput();

    // InternalImageType::ConstPointer smoothedImage = reader->GetOutput();
    // smoothedImage =SmoothImage<ImageType>(reader->GetOutput() , );

    using ResampleFilterType = itk::ResampleImageFilter<InternalImageType, OutputImageType>;

    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    using TransformType = itk::IdentityTransform<double, 3>;

    using InterpolatorType = itk::LinearInterpolateImageFunction<InternalImageType, double>;
    using InterpolatorType2 = itk::NearestNeighborInterpolateImageFunction<InternalImageType, double>;

    InterpolatorType::Pointer  interpolator = InterpolatorType::New();
    InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();

    resampler->SetInterpolator(interpolator);
    if (nn == 1)
    {
      resampler->SetInterpolator(interpolator2);
    }

    InternalImageType::IndexType ind;
    ind.Fill(1);
    resampler->SetDefaultPixelValue(inputImage->GetPixel(ind)); // zero regions without source

    // Use the inputImage as initial template
    resampler->SetOutputParametersFromImage(inputImage);
    // Reset spacing by explicit specification
    std::cout << " out space " << spacing << std::endl;
    resampler->SetOutputSpacing(spacing);

    InputImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
    using SizeValueType = InputImageType::SizeType::SizeValueType;
    InputImageType::SizeType size;
    for (int i = 0; i < 3; i++)
    {
      size[i] = static_cast<SizeValueType>(inputSize[i] * inputSpacing[i] / spacing[i] + addvox);
    }

    std::cout << " output size " << size << " spc " << spacing << std::endl;
    resampler->SetSize(size);
    resampler->SetInput(smoothedImage);
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampler->SetTransform(transform);
    resampler->Update();
    ANTs::WriteImage<OutputImageType>(resampler->GetOutput(), argv[3]);
  }
  /*ADDING 4-dimensional images */
  if (Dimension == 4)
  {
    using InputPixelType = float;
    using InternalPixelType = float;
    using OutputPixelType = float;
    using InputImageType = itk::Image<InputPixelType, 4>;
    using InternalImageType = itk::Image<InternalPixelType, 4>;
    using OutputImageType = itk::Image<OutputPixelType, 4>;
    InputImageType::Pointer inputImage;
    ReadImage<InputImageType>(inputImage, argv[2]);
    const InputImageType::SpacingType & inputSpacing = inputImage->GetSpacing();
    OutputImageType::SpacingType        spacing;
    for (int i = 0; i < 4; i++)
    {
      spacing[i] = inputSpacing[i];
    }
    std::cout << " spacing " << spacing << " dim " << 4 << std::endl;
    bool dosmooth = true;
    if (argc > 4)
    {
      spacing[0] = atof(argv[4]);
    }
    if (argc > 5)
    {
      spacing[1] = atof(argv[5]);
    }
    if (argc > 6)
    {
      spacing[2] = atof(argv[6]);
    }
    if (argc > 7)
    {
      spacing[3] = atof(argv[7]);
    }
    if (argc > 8)
    {
      dosmooth = std::stoi(argv[8]);
    }
    int addvox = 0;
    if (argc > 9)
    {
      addvox = std::stoi(argv[9]);
    }
    bool nn = false;
    if (argc > 9)
    {
      nn = std::stoi(argv[10]);
    }
    std::cout << " spacing2 " << spacing << std::endl;
    InternalImageType::Pointer smoothedImage = inputImage;
    if (dosmooth)
    {
      for (int sm = 0; sm < 4; sm++)
      {
        using GaussianFilterType = itk::RecursiveGaussianImageFilter<OutputImageType, OutputImageType>;
        GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
        smootherX->SetInput(smoothedImage);
        const float sig = atof(argv[4 + sm]) / inputSpacing[sm] - 1.0;
        std::cout << " smoothing by : " << sig << " dir " << sm << std::endl;
        smootherX->SetSigma(sig);
        smootherX->SetDirection(sm);
        smootherX->SetNormalizeAcrossScale(false);
        if (sig > 0 && dosmooth)
        {
          try
          {
            smootherX->Update();
          }
          catch (const itk::ExceptionObject & excep)
          {
            std::cout << "Exception catched !" << std::endl;
            std::cout << excep << std::endl;
          }
          smoothedImage = smootherX->GetOutput();
        }
      }
    }
    // InternalImageType::ConstPointer smoothedImage = smootherY->GetOutput();
    // InternalImageType::ConstPointer smoothedImage = reader->GetOutput();
    // smoothedImage =SmoothImage<ImageType>(reader->GetOutput() , );
    using ResampleFilterType = itk::ResampleImageFilter<InternalImageType, OutputImageType>;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    using TransformType = itk::IdentityTransform<double, 4>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<InternalImageType, double>;
    using InterpolatorType2 = itk::NearestNeighborInterpolateImageFunction<InternalImageType, double>;
    InterpolatorType::Pointer  interpolator = InterpolatorType::New();
    InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();
    resampler->SetInterpolator(interpolator);
    if (nn == 1)
    {
      resampler->SetInterpolator(interpolator2);
    }
    InternalImageType::IndexType ind;
    ind.Fill(1);
    resampler->SetDefaultPixelValue(inputImage->GetPixel(ind)); // zero regions without source
    // Use the inputImage as initial template
    resampler->SetOutputParametersFromImage(inputImage);
    // Reset spacing by explicit specification
    std::cout << " out space " << spacing << std::endl;
    resampler->SetOutputSpacing(spacing);
    InputImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
    using SizeValueType = InputImageType::SizeType::SizeValueType;
    InputImageType::SizeType size;
    for (int i = 0; i < 4; i++)
    {
      size[i] = static_cast<SizeValueType>(inputSize[i] * inputSpacing[i] / spacing[i] + addvox);
    }
    std::cout << " output size " << size << " spc " << spacing << std::endl;
    resampler->SetSize(size);
    resampler->SetInput(smoothedImage);
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampler->SetTransform(transform);
    resampler->Update();
    ANTs::WriteImage<OutputImageType>(resampler->GetOutput(), argv[3]);
  }

  return EXIT_SUCCESS;
}
} // namespace ants
