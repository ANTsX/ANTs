#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include "itkAffineTransform.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCSVArray2DDataObject.h"
#include "itkCSVArray2DFileReader.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkImage.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMap.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkShapeLabelObject.h"
#include "itkStatisticsImageFilter.h"
#include "itkStatisticsLabelMapFilter.h"
#include "itkTransformFileWriter.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

namespace ants
{

template <unsigned int Dimension, typename CentroidType>
std::string formatCentroid(const CentroidType& centroid)
{
  std::ostringstream oss;
  oss << "[";
  for (unsigned int i = 0; i < Dimension; ++i) {
    oss << std::fixed << std::setprecision(4) << centroid[i];
    if (i < Dimension - 1) oss << ", ";
  }
  oss << "]";
  return oss.str();
}

template <unsigned int Dimension>
std::string formatAxesLengths(const std::vector<double>& axesLengths)
{
  std::ostringstream oss;
  oss << "[";
  for (unsigned int i = 0; i < Dimension; ++i)
  {
    oss << std::fixed << std::setprecision(4) << axesLengths[i];
    if (i < Dimension - 1) oss << ", ";
  }
  oss << "]";
  return oss.str();
}

template <unsigned int Dimension>
std::string formatBoundingBox(const itk::ImageRegion<Dimension>& region)
{
  std::ostringstream oss;
  auto index = region.GetIndex();
  auto size = region.GetSize();

  // Append the starting index and size of each dimension to the stream
  oss << "[";
  for (unsigned int i = 0; i < Dimension; ++i)
  {
      oss << index[i];
      if (i < Dimension - 1) oss << ", ";
  }
  oss << ", ";
  for (unsigned int i = 0; i < Dimension; ++i)
  {
    oss << index[i] + size[i] - 1;
    if (i < Dimension - 1) oss << ", ";
  }
  oss << "]";
  return oss.str();
}

template <unsigned int ImageDimension>
int
LabelGeometryMeasures(int argc, char * argv[])
{
  using LabelType = unsigned int;
  using LabelImageType = itk::Image<LabelType, ImageDimension>;
  // mask image is used to binarize the labels, used to mask the intensity image
  // to better initialize the histogram for stats
  using MaskType = unsigned char;
  using MaskImageType = itk::Image<MaskType, ImageDimension>;
  using RealType = float;
  using RealImageType = itk::Image<RealType, ImageDimension>;

  typename LabelImageType::Pointer labelImage = LabelImageType::New();
  ReadImage<LabelImageType>(labelImage, argv[2]);

  typename RealImageType::Pointer intensityImage = RealImageType::New();
  bool intensityImageUsed = false;
  if (argc > 3 && std::string(argv[3]) != "none" && std::string(argv[3]) != "na")
  {
    ReadImage<RealImageType>(intensityImage, argv[3]);
    intensityImageUsed = true;
  }

  bool writeCSV = false;

  if (argc > 4 && std::string(argv[4]) != "none" && std::string(argv[4]) != "na")
  {
    writeCSV = true;
  }

  using FilterType = itk::LabelImageToShapeLabelMapFilter<LabelImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetComputeOrientedBoundingBox(false);
  filter->SetComputePerimeter(true);
  filter->SetComputeFeretDiameter(false); // slow for large labels eg brain mask
  filter->SetInput(labelImage);
  filter->Update();

  using LabelMapType = typename FilterType::OutputImageType;
  typename LabelMapType::Pointer labelMap = filter->GetOutput();

  using StatisticsFilterType = itk::LabelStatisticsImageFilter<RealImageType, LabelImageType>;
  typename StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  if (intensityImageUsed)
  {

    // Threshold label image to create a binary mask
    auto thresholdFilter = itk::BinaryThresholdImageFilter<LabelImageType, MaskImageType>::New();
    thresholdFilter->SetInput(labelImage);
    thresholdFilter->SetLowerThreshold(1);
    thresholdFilter->SetUpperThreshold(itk::NumericTraits<LabelType>::max());
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
    auto binaryMask = thresholdFilter->GetOutput();

    // mask the intensity image with the binary mask
    auto maskFilter = itk::MaskImageFilter<RealImageType, MaskImageType>::New();
    maskFilter->SetInput(intensityImage);
    maskFilter->SetMaskImage(binaryMask);
    maskFilter->SetOutsideValue(0);
    maskFilter->Update();
    auto maskedScalarImage = maskFilter->GetOutput();

    // Iterate through the scalar image and mask
    itk::ImageRegionConstIterator<RealImageType> scalarIt(maskedScalarImage, maskedScalarImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<MaskImageType> maskIt(binaryMask, binaryMask->GetLargestPossibleRegion());

    float minValue = itk::NumericTraits<float>::max();
    float maxValue = itk::NumericTraits<float>::NonpositiveMin();

    while (!scalarIt.IsAtEnd())
    {
      if (maskIt.Get() > 0)  // Only consider voxels inside the mask
      {
        float value = scalarIt.Get();
        if (value < minValue)
          minValue = value;
        if (value > maxValue)
          maxValue = value;
      }
      ++scalarIt;
      ++maskIt;
   }

    statisticsFilter->SetInput(intensityImage);
    statisticsFilter->SetLabelInput(labelImage);
    statisticsFilter->SetUseHistograms(true);
    statisticsFilter->SetHistogramParameters(4096, minValue, maxValue);
    statisticsFilter->Update();
  }

  std::vector<std::string> columnHeaders = {"Label", "VolumeInVoxels", "VolumeInMillimeters", "SurfaceAreaInMillimetersSquared",
                                            "Eccentricity", "Elongation", "Roundness", "Flatness"};

  if (writeCSV)
  {
    if (ImageDimension == 2)
    {
      columnHeaders.push_back("Centroid_x");
      columnHeaders.push_back("Centroid_y");
      columnHeaders.push_back("AxesLength_x");
      columnHeaders.push_back("AxesLength_y");
      columnHeaders.push_back("BoundingBoxLower_x");
      columnHeaders.push_back("BoundingBoxLower_y");
      columnHeaders.push_back("BoundingBoxUpper_x");
      columnHeaders.push_back("BoundingBoxUpper_y");
    }
    else if (ImageDimension == 3)
    {
      columnHeaders.push_back("Centroid_x");
      columnHeaders.push_back("Centroid_y");
      columnHeaders.push_back("Centroid_z");
      columnHeaders.push_back("AxesLength_x");
      columnHeaders.push_back("AxesLength_y");
      columnHeaders.push_back("AxesLength_z");
      columnHeaders.push_back("BoundingBoxLower_x");
      columnHeaders.push_back("BoundingBoxLower_y");
      columnHeaders.push_back("BoundingBoxLower_z");
      columnHeaders.push_back("BoundingBoxUpper_x");
      columnHeaders.push_back("BoundingBoxUpper_y");
      columnHeaders.push_back("BoundingBoxUpper_z");
    }
  }
  else
  {
    columnHeaders.push_back("Centroid");
    columnHeaders.push_back("AxesLengths");
    columnHeaders.push_back("BoundingBox");
  }

  if (intensityImageUsed)
  {
    columnHeaders.push_back("MeanIntensity");
    columnHeaders.push_back("SigmaIntensity");
    columnHeaders.push_back("MedianIntensity");
    columnHeaders.push_back("MinIntensity");
    columnHeaders.push_back("MaxIntensity");
    columnHeaders.push_back("IntegratedIntensity");
  }

  std::vector<std::vector<std::string>> data;
  for (unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); ++i)
  {
    auto labelObject = labelMap->GetNthLabelObject(i);
    if (labelObject->GetLabel() == 0)
    {
      continue;  // Skip background
    }

    // Get principal moments and use them to calculate eccentricity and axes lengths
    auto principalMoments = labelObject->GetPrincipalMoments();

    // define in 3D such that it describes the ellipse with axes propotional to the
    // two largest principal moments. This is useful eg for cortex, where the thickness
    // is much smaller than the other two dimensions.
    //
    // Roundness gives a more general measure of deviation from a sphere, including all three dimensions.
    //
    double lambda1 = principalMoments[ImageDimension - 2];
    double lambdaN = principalMoments[ImageDimension - 1];
    double eccentricity = 0.0;

    if (!itk::Math::FloatAlmostEqual(lambda1, 0.0))
    {
      eccentricity = std::sqrt(1.0 - (lambda1 * lambda1) / (lambdaN * lambdaN));
    }

    // calculate axes lengths
    std::vector<double> axesLengths(ImageDimension, 0.0);

    for (unsigned int idx = 0; idx < ImageDimension; ++idx)
    {
      if (principalMoments[idx] > 0)
        axesLengths[idx] = 4.0 * std::sqrt(principalMoments[idx]);
      else
        axesLengths[idx] = 0.0;
    }
    // row is a vector of str
    std::vector<std::string> row = {
      std::to_string(labelObject->GetLabel()),
      std::to_string(labelObject->GetNumberOfPixels()),
      std::to_string(labelObject->GetPhysicalSize()),
      std::to_string(labelObject->GetPerimeter()),
      std::to_string(eccentricity),
      std::to_string(labelObject->GetElongation()),
      std::to_string(labelObject->GetRoundness()),
      std::to_string(labelObject->GetFlatness())
    };

    if (writeCSV)
    {
      if (ImageDimension == 2)
      {
       auto bbLowerX = labelObject->GetBoundingBox().GetIndex()[0];
       auto bbLowerY = labelObject->GetBoundingBox().GetIndex()[1];

       auto bbUpperX = bbLowerX + labelObject->GetBoundingBox().GetSize()[0] - 1;
       auto bbUpperY = bbLowerY + labelObject->GetBoundingBox().GetSize()[1] - 1;

       row.push_back(std::to_string(labelObject->GetCentroid()[0]));
       row.push_back(std::to_string(labelObject->GetCentroid()[1]));
       row.push_back(std::to_string(axesLengths[0]));
       row.push_back(std::to_string(axesLengths[1]));
       row.push_back(std::to_string(bbLowerX));
       row.push_back(std::to_string(bbLowerY));
       row.push_back(std::to_string(bbUpperX));
       row.push_back(std::to_string(bbUpperY));
      }
      else if (ImageDimension == 3)
      {
       auto bbLowerX = labelObject->GetBoundingBox().GetIndex()[0];
       auto bbLowerY = labelObject->GetBoundingBox().GetIndex()[1];
       auto bbLowerZ = labelObject->GetBoundingBox().GetIndex()[2];

       auto bbUpperX = bbLowerX + labelObject->GetBoundingBox().GetSize()[0] - 1;
       auto bbUpperY = bbLowerY + labelObject->GetBoundingBox().GetSize()[1] - 1;
       auto bbUpperZ = bbLowerZ + labelObject->GetBoundingBox().GetSize()[2] - 1;

       row.push_back(std::to_string(labelObject->GetCentroid()[0]));
       row.push_back(std::to_string(labelObject->GetCentroid()[1]));
       row.push_back(std::to_string(labelObject->GetCentroid()[2]));
       row.push_back(std::to_string(axesLengths[0]));
       row.push_back(std::to_string(axesLengths[1]));
       row.push_back(std::to_string(axesLengths[2]));
       row.push_back(std::to_string(bbLowerX));
       row.push_back(std::to_string(bbLowerY));
       row.push_back(std::to_string(bbLowerZ));
       row.push_back(std::to_string(bbUpperX));
       row.push_back(std::to_string(bbUpperY));
       row.push_back(std::to_string(bbUpperZ));
      }
    }
    else
    {
      row.push_back(formatCentroid<ImageDimension>(labelObject->GetCentroid()));
      row.push_back(formatAxesLengths<ImageDimension>(axesLengths));
      row.push_back(formatBoundingBox<ImageDimension>(labelObject->GetBoundingBox()));
    }

    if (intensityImageUsed)
    {
      row.push_back(std::to_string(statisticsFilter->GetMean(labelObject->GetLabel())));
      row.push_back(std::to_string(statisticsFilter->GetSigma(labelObject->GetLabel())));
      row.push_back(std::to_string(statisticsFilter->GetMedian(labelObject->GetLabel())));
      row.push_back(std::to_string(statisticsFilter->GetMinimum(labelObject->GetLabel())));
      row.push_back(std::to_string(statisticsFilter->GetMaximum(labelObject->GetLabel())));
      row.push_back(std::to_string(statisticsFilter->GetSum(labelObject->GetLabel())));
    }

    data.push_back(row);
  }

  std::ostream* outStream = &std::cout; // Default to std::cout
  std::string delimiter = "\t"; // Default to tab delimiter

  std::ofstream outFile;

  if (writeCSV)
  {
    outFile.open(argv[4]);
    if (!outFile.is_open())
    {
      std::cerr << "Error opening output file" << std::endl;
      return EXIT_FAILURE;
    }
    outStream = &outFile; // Point to the file stream
    delimiter = ","; // Use comma for files
  }

  // Write column headers
  for (size_t i = 0; i < columnHeaders.size(); ++i)
  {
    *outStream << columnHeaders[i];
    if (i < columnHeaders.size() - 1)
      *outStream << delimiter;
  }
  *outStream << "\n";

  // Write data
  for (const auto& row : data)
  {
    for (size_t j = 0; j < row.size(); ++j)
    {
      *outStream << row[j];
      if (j < row.size() - 1)
        *outStream << delimiter;
    }
    *outStream << "\n";
  }

  if (outFile.is_open())
  {
    outFile.close();
  }
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
LabelGeometryMeasures(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "LabelGeometryMeasures");

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
    std::cout << "Usage 1: " << argv[0] << " imageDimension labelImage [intensityImage] [csvFile]" << std::endl;
    std::cout << "To output to csvFile without an intensity image, use \"none\" or \"na\" as the third arg" << std::endl;
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
      LabelGeometryMeasures<2>(argc, argv);
    }
    break;
    case 3:
    {
      LabelGeometryMeasures<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


} // namespace ants
