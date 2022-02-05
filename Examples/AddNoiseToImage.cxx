#include "antsAllocImage.h"
#include "antsCommandLineParser.h"
#include "antsUtilities.h"

#include "ReadWriteData.h"

#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkSaltAndPepperNoiseImageFilter.h"
#include "itkShotNoiseImageFilter.h"
#include "itkSpeckleNoiseImageFilter.h"

#include "ANTsVersion.h"

namespace ants
{


template <typename TFilter>
class CommandProgressUpdate : public itk::Command
{
public:
  using Self = CommandProgressUpdate<TFilter>;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<CommandProgressUpdate<TFilter>>;
  itkNewMacro(CommandProgressUpdate);

protected:
  CommandProgressUpdate() = default;
  ;

  using FilterType = TFilter;

  unsigned int m_CurrentProgress{ 0 };

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    auto * po = dynamic_cast<itk::ProcessObject *>(caller);
    if (!po)
      return;
    //    std::cout << po->GetProgress() << std::endl;
    if (typeid(event) == typeid(itk::ProgressEvent))
    {
      if (this->m_CurrentProgress < 99)
      {
        this->m_CurrentProgress++;
        if (this->m_CurrentProgress % 10 == 0)
        {
          std::cout << this->m_CurrentProgress << std::flush;
        }
        else
        {
          std::cout << "*" << std::flush;
        }
      }
    }
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    auto * po = dynamic_cast<itk::ProcessObject *>(const_cast<itk::Object *>(object));
    if (!po)
      return;

    if (typeid(event) == typeid(itk::ProgressEvent))
    {
      if (this->m_CurrentProgress < 99)
      {
        this->m_CurrentProgress++;
        if (this->m_CurrentProgress % 10 == 0)
        {
          std::cout << this->m_CurrentProgress << std::flush;
        }
        else
        {
          std::cout << "*" << std::flush;
        }
      }
    }
  }
};

template <unsigned int ImageDimension>
int
AddNoise(itk::ants::CommandLineParser * parser)
{
  using RealType = float;

  using OptionType = typename itk::ants::CommandLineParser::OptionType;

  bool                                                       verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  if (verbose)
  {
    std::cout << std::endl << "Running for " << ImageDimension << "-dimensional images." << std::endl << std::endl;
  }

  using ImageType = itk::Image<RealType, ImageDimension>;
  typename ImageType::Pointer inputImage = nullptr;

  // typedef itk::Image<RealType, ImageDimension> MaskImageType;
  // typename MaskImageType::Pointer maskImage = nullptr;

  typename OptionType::Pointer inputImageOption = parser->GetOption("input-image");
  if (inputImageOption && inputImageOption->GetNumberOfFunctions())
  {
    std::string inputFile = inputImageOption->GetFunction(0)->GetName();
    ReadImage<ImageType>(inputImage, inputFile.c_str());
    inputImage->Update();
    inputImage->DisconnectPipeline();
  }
  else
  {
    if (verbose)
    {
      std::cerr << "Input image not specified." << std::endl;
    }
    return EXIT_FAILURE;
  }

  typename ImageType::Pointer outputImage = nullptr;

  typename OptionType::Pointer noiseModelOption = parser->GetOption("noise-model");

  if (!noiseModelOption)
  {
    if (verbose)
    {
      std::cerr << "Input noise model not specified." << std::endl;
    }
    return EXIT_FAILURE;
  }
  if (noiseModelOption->GetNumberOfFunctions())
  {
    std::string noiseModel = noiseModelOption->GetFunction(0)->GetName();
    ConvertToLowerCase(noiseModel);

    if (std::strcmp(noiseModel.c_str(), "additivegaussian") == 0)
    {
      using NoiseFilterType = itk::AdditiveGaussianNoiseImageFilter<ImageType, ImageType>;
      typename NoiseFilterType::Pointer noiser = NoiseFilterType::New();
      noiser->SetInput(inputImage);

      RealType mean = 0.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        std::string meanString = noiseModelOption->GetFunction(0)->GetParameter(0);
        mean = parser->Convert<RealType>(meanString);
      }
      RealType sd = 1.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 1)
      {
        std::string sdString = noiseModelOption->GetFunction(0)->GetParameter(1);
        sd = parser->Convert<RealType>(sdString);
      }

      noiser->SetMean(mean);
      noiser->SetStandardDeviation(sd);

      if (verbose)
      {
        std::cout << "Noise model: Additive Gaussian" << std::endl;
        std::cout << "  mean = " << noiser->GetMean() << std::endl;
        std::cout << "  standard deviation = " << noiser->GetStandardDeviation() << std::endl;
      }

      outputImage = noiser->GetOutput();
      outputImage->Update();
      outputImage->DisconnectPipeline();
    }
    else if (std::strcmp(noiseModel.c_str(), "saltandpepper") == 0)
    {
      using NoiseFilterType = itk::SaltAndPepperNoiseImageFilter<ImageType, ImageType>;
      typename NoiseFilterType::Pointer noiser = NoiseFilterType::New();
      noiser->SetInput(inputImage);

      RealType probability = 0.01;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        std::string probabilityString = noiseModelOption->GetFunction(0)->GetParameter(0);
        probability = parser->Convert<RealType>(probabilityString);
      }
      RealType salt = 0.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 1)
      {
        std::string saltString = noiseModelOption->GetFunction(0)->GetParameter(1);
        salt = parser->Convert<RealType>(saltString);
        noiser->SetSaltValue(salt);
      }
      RealType pepper = 1.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 2)
      {
        std::string pepperString = noiseModelOption->GetFunction(0)->GetParameter(2);
        pepper = parser->Convert<RealType>(pepperString);
        noiser->SetPepperValue(pepper);
      }
      noiser->SetProbability(probability);

      if (verbose)
      {
        std::cout << "Noise model: Salt and pepper" << std::endl;
        std::cout << "  probability = " << noiser->GetProbability() << std::endl;
        std::cout << "  salt = " << noiser->GetSaltValue() << std::endl;
        std::cout << "  pepper = " << noiser->GetPepperValue() << std::endl;
      }

      outputImage = noiser->GetOutput();
      outputImage->Update();
      outputImage->DisconnectPipeline();
    }
    else if (std::strcmp(noiseModel.c_str(), "shot") == 0)
    {
      using NoiseFilterType = itk::ShotNoiseImageFilter<ImageType, ImageType>;
      typename NoiseFilterType::Pointer noiser = NoiseFilterType::New();
      noiser->SetInput(inputImage);

      RealType scale = 1.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        std::string scaleString = noiseModelOption->GetFunction(0)->GetParameter(0);
        scale = parser->Convert<RealType>(scaleString);
      }

      noiser->SetScale(scale);

      if (verbose)
      {
        std::cout << "Noise model: Shot" << std::endl;
        std::cout << "  scale = " << noiser->GetScale() << std::endl;
      }

      outputImage = noiser->GetOutput();
      outputImage->Update();
      outputImage->DisconnectPipeline();
    }
    else if (std::strcmp(noiseModel.c_str(), "speckle") == 0)
    {
      using NoiseFilterType = itk::SpeckleNoiseImageFilter<ImageType, ImageType>;
      typename NoiseFilterType::Pointer noiser = NoiseFilterType::New();
      noiser->SetInput(inputImage);

      RealType sd = 1.0;
      if (noiseModelOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        std::string sdString = noiseModelOption->GetFunction(0)->GetParameter(0);
        sd = parser->Convert<RealType>(sdString);
      }
      noiser->SetStandardDeviation(sd);

      if (verbose)
      {
        std::cout << "Noise model: Speckle" << std::endl;
        std::cout << "  standard deviation = " << noiser->GetStandardDeviation() << std::endl;
      }

      outputImage = noiser->GetOutput();
      outputImage->Update();
      outputImage->DisconnectPipeline();
    }
    else
    {
      if (verbose)
      {
        std::cerr << "Unrecognized noise model:  " << noiseModel << ".  See help menu." << std::endl;
      }
      return EXIT_FAILURE;
    }
  }

  /**
   * output
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");
  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    if (outputOption->GetFunction(0)->GetNumberOfParameters() == 0)
    {
      ANTs::WriteImage<ImageType>(outputImage, (outputOption->GetFunction(0)->GetName()).c_str());
    }
  }

  return EXIT_SUCCESS;
}

void
InitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  using OptionType = itk::ants::CommandLineParser::OptionType;

  {
    std::string description = std::string("This option forces the image to be treated as a specified-") +
                              std::string("dimensional image.  If not specified, the program tries to ") +
                              std::string("infer the dimensionality from the input image.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("image-dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3/4");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("A scalar image is expected as input for noise correction.  ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input-image");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputImageFilename");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Use different noise models each with its own (default) parameters. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("noise-model");
    option->SetShortName('n');
    option->SetUsageOption(0, "AdditiveGaussian[<mean=0.0>,<standardDeviation=1.0>]");
    option->SetUsageOption(1, "SaltAndPepper[<probability=0.01>,<saltValue=minPixelType>,<pepperValue=maxPixelType>]");
    option->SetUsageOption(2, "Shot[<scale=1.0>]");
    option->SetUsageOption(3, "Speckle[<standardDeviation=1.0>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("The output consists of the noise corrected version of the ") + std::string("input image.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "noiseCorruptedImage");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Get version information.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("version");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Verbose output.");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('v');
    option->SetLongName("verbose");
    option->SetUsageOption(0, "(0)/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
AddNoiseToImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "AddNoiseToImage");

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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand(argv[0]);

  std::string commandDescription = std::string("Add various types of noise to an image.");

  parser->SetCommandDescription(commandDescription);
  InitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc == 1)
  {
    parser->PrintMenu(std::cerr, 5, false);
    return EXIT_FAILURE;
  }
  else if (parser->GetOption("help")->GetFunction() &&
           parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_SUCCESS;
  }
  else if (parser->GetOption('h')->GetFunction() &&
           parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }
  // Show automatic version
  itk::ants::CommandLineParser::OptionType::Pointer versionOption = parser->GetOption("version");
  if (versionOption && versionOption->GetNumberOfFunctions())
  {
    std::string versionFunction = versionOption->GetFunction(0)->GetName();
    ConvertToLowerCase(versionFunction);
    if (versionFunction.compare("1") == 0 || versionFunction.compare("true") == 0)
    {
      // Print Version Information
      std::cout << ANTs::Version::ExtendedVersionString() << std::endl;
      return EXIT_SUCCESS;
    }
  }
  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption("image-dimensionality");
  if (dimOption && dimOption->GetNumberOfFunctions())
  {
    dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
  }
  else
  {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption = parser->GetOption("input-image");
    if (imageOption && imageOption->GetNumberOfFunctions() > 0)
    {
      if (imageOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        filename = imageOption->GetFunction(0)->GetParameter(0);
      }
      else
      {
        filename = imageOption->GetFunction(0)->GetName();
      }
    }
    else
    {
      std::cerr << "No input images were specified.  Specify an input image"
                << " with the -i option" << std::endl;
      return EXIT_FAILURE;
    }
    itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::IOFileModeEnum::ReadMode);
    dimension = imageIO->GetNumberOfDimensions();
  }

  switch (dimension)
  {
    case 2:
    {
      return AddNoise<2>(parser);
    }
    break;
    case 3:
    {
      return AddNoise<3>(parser);
    }
    break;
    case 4:
    {
      return AddNoise<4>(parser);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
