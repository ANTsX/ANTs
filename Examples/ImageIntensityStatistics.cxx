#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"

#include <vector>
#include <algorithm>
#include <iomanip>

namespace ants
{
template <unsigned int ImageDimension>
int
ImageIntensityStatistics(int argc, char * argv[])
{
  using LabelType = int;
  using RealType = float;

  using LabelImageType = itk::Image<LabelType, ImageDimension>;
  using RealImageType = itk::Image<RealType, ImageDimension>;

  typename RealImageType::Pointer intensityImage = RealImageType::New();
  ReadImage<RealImageType>(intensityImage, argv[2]);

  typename LabelImageType::Pointer labelImage = LabelImageType::New();
  if (argc > 3)
  {
    ReadImage<LabelImageType>(labelImage, argv[3]);
  }
  else
  {
    labelImage->CopyInformation(intensityImage);
    labelImage->SetRegions(intensityImage->GetLargestPossibleRegion());
    labelImage->Allocate();
    labelImage->FillBuffer(itk::NumericTraits<LabelType>::One);
  }

  std::vector<LabelType>                   labels;
  itk::ImageRegionIterator<LabelImageType> It(labelImage, labelImage->GetLargestPossibleRegion());

  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (It.Get() != 0 && std::find(labels.begin(), labels.end(), It.Get()) == labels.end())
    {
      labels.push_back(It.Get());
    }
  }
  std::sort(labels.begin(), labels.end());

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();

  itk::ImageRegionConstIterator<RealImageType> ItI(intensityImage, intensityImage->GetLargestPossibleRegion());
  for (ItI.GoToBegin(), It.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++It)
  {
    if (It.Get() == 0)
    {
      continue;
    }
    RealType value = ItI.Get();

    if (value < minValue)
    {
      minValue = value;
    }
    else if (value > maxValue)
    {
      maxValue = value;
    }
  }

  using HistogramGeneratorType = itk::LabelStatisticsImageFilter<RealImageType, LabelImageType>;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput(intensityImage);
  stats->SetLabelInput(labelImage);
  stats->UseHistogramsOn();
  stats->SetHistogramParameters(200, minValue, maxValue);
  stats->Update();

  // Calculate moments for skewness and kurtosis calculations

  std::vector<RealType> m3(labels.size(), 0.0);
  std::vector<RealType> m4(labels.size(), 0.0);
  std::vector<RealType> N(labels.size(), 0.0);
  for (ItI.GoToBegin(), It.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++It)
  {
    LabelType label = It.Get();
    if (label == 0)
    {
      continue;
    }
    auto it = std::find(labels.begin(), labels.end(), label);
    if (it == labels.end())
    {
      std::cerr << "Label not found.  Shouldn't get here." << std::endl;
      return EXIT_FAILURE;
    }
    RealType difference = ItI.Get() - static_cast<RealType>(stats->GetMean(label));

    unsigned long index = it - labels.begin();

    m3[index] += (difference * difference * difference);
    m4[index] += (m3[index] * difference);
    N[index] += itk::NumericTraits<RealType>::OneValue();
  }
  for (unsigned int n = 0; n < N.size(); n++)
  {
    m3[n] /= N[n];
    m4[n] /= N[n];
  }

  //   std::cout << "                                       "
  //             << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw(8) << "Label" << std::setw(14) << "Mean" << std::setw(14) << "Sigma" << std::setw(14)
            << "Skewness" << std::setw(14) << "Kurtosis" << std::setw(14) << "Entropy" << std::setw(14) << "Sum"
            << std::setw(14) << "5th%" << std::setw(14) << "95th%" << std::setw(14) << "Min" << std::setw(14) << "Max"
            << std::endl;

  std::vector<LabelType>::iterator it;
  for (it = labels.begin(); it != labels.end(); ++it)
  {
    using HistogramType = typename HistogramGeneratorType::HistogramType;
    const HistogramType * histogram = stats->GetHistogram(*it);

    RealType fifthPercentileValue = histogram->Quantile(0, 0.05);
    RealType ninetyFifthPercentileValue = histogram->Quantile(0, 0.95);
    RealType entropy = 0.0;
    for (unsigned int i = 0; i < histogram->Size(); i++)
    {
      RealType p =
        static_cast<RealType>(histogram->GetFrequency(i, 0)) / static_cast<RealType>(histogram->GetTotalFrequency());
      if (p > 0)
      {
        entropy += static_cast<RealType>(-p * std::log(p) / std::log(2.0f));
      }
    }

    auto          it2 = std::find(labels.begin(), labels.end(), *it);
    unsigned long index = it2 - labels.begin();

    RealType m2 = itk::Math::sqr(stats->GetSigma(*it));
    RealType k2 = (N[index]) / (N[index] - 1.0f) * m2;

    RealType prefactor3 = itk::Math::sqr(N[index]) / ((N[index] - 1.0f) * (N[index] - 2.0f));
    RealType k3 = prefactor3 * m3[index];

    RealType prefactor4 = itk::Math::sqr(N[index]) / ((N[index] - 1.0f) * (N[index] - 2.0f) * (N[index] - 3.0f));
    RealType k4 = prefactor4 * ((N[index] + 1.0f) * m4[index] - 3.0f * (N[index] - 1.0f) * itk::Math::sqr(m2));

    RealType skewness = k3 / std::sqrt(k2 * k2 * k2);
    RealType kurtosis = k4 / itk::Math::sqr(k2);

    std::cout << std::setw(8) << *it;
    std::cout << std::setw(14) << stats->GetMean(*it);
    std::cout << std::setw(14) << stats->GetSigma(*it);
    std::cout << std::setw(14) << skewness;
    std::cout << std::setw(14) << kurtosis;
    std::cout << std::setw(14) << entropy;
    std::cout << std::setw(14) << stats->GetSum(*it);
    std::cout << std::setw(14) << fifthPercentileValue;
    std::cout << std::setw(14) << ninetyFifthPercentileValue;
    std::cout << std::setw(14) << stats->GetMinimum(*it);
    std::cout << std::setw(14) << stats->GetMaximum(*it);
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}

int
ImageIntensityStatistics(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ImageIntensityStatistics");

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

  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage <labelImage>" << std::endl;
    exit(1);
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      return ImageIntensityStatistics<2>(argc, argv);
      break;
    case 3:
      return ImageIntensityStatistics<3>(argc, argv);
      break;
    case 4:
      return ImageIntensityStatistics<4>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
} // namespace ants
