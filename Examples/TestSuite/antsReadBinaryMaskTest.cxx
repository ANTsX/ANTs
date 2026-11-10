#include "ReadWriteData.h"
#include "itkArray.h"
#include "itkExtractImageFilter.h"
#include "itkPointSet.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

void
Require(bool condition, const char * message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

template <unsigned int Dimension>
void
CheckMask(const std::string & directory)
{
  using FloatImage = itk::Image<float, Dimension>;
  using MaskImage = itk::Image<unsigned char, Dimension>;
  const std::array<float, 9> values = { 0.0f, -0.0f, 0.25f, -0.25f, 1.0f, -1.0f,
                                      256.0f, -256.0f, 512.0f };
  auto input = FloatImage::New();
  typename FloatImage::SizeType size;
  size.Fill(8);
  input->SetRegions(size);
  typename FloatImage::PointType origin;
  typename FloatImage::SpacingType spacing;
  typename FloatImage::DirectionType direction;
  direction.SetIdentity();
  direction[0][0] = -1;
  for (unsigned int d = 0; d < Dimension; ++d)
  {
    origin[d] = 10.0 + 3.0 * d;
    spacing[d] = 1.0 + 0.5 * d;
  }
  input->SetOrigin(origin);
  input->SetSpacing(spacing);
  input->SetDirection(direction);
  input->AllocateInitialized();
  std::copy(values.begin(), values.end(), input->GetBufferPointer());
  const auto filename = directory + "/mask" + std::to_string(Dimension) + ".mha";
  auto writer = itk::ImageFileWriter<FloatImage>::New();
  writer->SetFileName(filename);
  writer->SetInput(input);
  writer->Update();

  typename MaskImage::Pointer mask;
  Require(ReadBinaryMask<MaskImage>(mask, filename.c_str()), "Mask read failed");
  Require(mask->GetOrigin() == origin && mask->GetSpacing() == spacing && mask->GetDirection() == direction &&
            mask->GetLargestPossibleRegion() == input->GetLargestPossibleRegion(), "Mask geometry changed");
  for (size_t i = 0; i < mask->GetLargestPossibleRegion().GetNumberOfPixels(); ++i)
  {
    const unsigned char expected = i < values.size() && values[i] != 0.0f ? 1 : 0;
    Require(mask->GetBufferPointer()[i] == expected, "Incorrect nonzero mask conversion");
  }

  // Exercise the intensity-point-set reader with the same file, not just the
  // shared conversion helper. Values that used to truncate/wrap to zero must
  // contribute points at their original physical coordinates.
  if constexpr (Dimension == 2)
  {
    using PointSet = itk::PointSet<itk::Array<float>, Dimension>;
    input->FillBuffer(1.0f);
    const auto intensityFilename = directory + "/intensity.mha";
    writer->SetFileName(intensityFilename);
    writer->Update();
    typename PointSet::Pointer points;
    Require((ReadImageIntensityPointSet<FloatImage, MaskImage, PointSet>(
              points, intensityFilename.c_str(), filename.c_str(), std::vector<unsigned int>(Dimension, 0), 1.0)),
            "Intensity point set read failed");
    Require(points->GetNumberOfPoints() == values.size() - 2, "Incorrect masked point count");
    for (size_t i = 2; i < values.size(); ++i)
    {
      typename FloatImage::PointType expected;
      input->TransformIndexToPhysicalPoint(input->ComputeIndex(i), expected);
      typename PointSet::PointType actual;
      Require(points->GetPoint(i - 2, &actual), "Missing masked point");
      for (unsigned int d = 0; d < Dimension; ++d)
      {
        Require(std::abs(actual[d] - expected[d]) < 1e-5, "Masked point geometry changed");
      }
    }
  }
  const std::array<float, 3> invalidValues = { std::numeric_limits<float>::quiet_NaN(),
                                             std::numeric_limits<float>::infinity(),
                                             -std::numeric_limits<float>::infinity() };
  for (size_t i = 0; i < invalidValues.size(); ++i)
  {
    input->FillBuffer(1.0f);
    const auto lastIndex = input->ComputeIndex(input->GetLargestPossibleRegion().GetNumberOfPixels() - 1);
    input->SetPixel(lastIndex, invalidValues[i]);
    input->Modified();
    const auto invalidFilename = directory + "/invalid" + std::to_string(Dimension) + "-" + std::to_string(i) + ".mha";
    writer->SetFileName(invalidFilename);
    writer->Update();
    mask = MaskImage::New();
    bool threw = false;
    try
    {
      ReadBinaryMask<MaskImage>(mask, invalidFilename.c_str());
    }
    catch (const itk::ExceptionObject & error)
    {
      threw = true;
      const std::string description = error.GetDescription();
      Require(description.find(invalidFilename) != std::string::npos &&
                description.find("NaN or Inf") != std::string::npos,
              "Nonfinite mask exception lacks filename or reason");
    }
    Require(threw, "Nonfinite mask did not raise an exception");
    Require(mask.IsNull(), "Nonfinite mask left an output image");
  }
  Require(!ReadBinaryMask<MaskImage>(mask, (directory + "/missing.mha").c_str()) && mask.IsNull(),
          "Failed mask read did not clear output");
}

int
main(int argc, char ** argv)
{
  try
  {
    Require(argc == 2, "Expected a test output directory");
    std::filesystem::create_directories(argv[1]);
    CheckMask<2>(argv[1]);
    CheckMask<3>(argv[1]);
    CheckMask<4>(argv[1]);
  }
  catch (const std::exception & error)
  {
    std::cerr << error.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
