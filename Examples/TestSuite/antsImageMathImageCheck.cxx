#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkVectorImage.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

// Write the script's textual fixture values directly, avoiding ImageIO's
// interpretation of ASCII multicomponent MetaImages.
template <unsigned int Dimension>
int
WriteFixture(const char * file)
{
  using ImageType = itk::VectorImage<float, Dimension>;
  std::ifstream input(file);
  typename ImageType::SizeType size;
  size.Fill(0);
  unsigned int components = 1;
  std::string line;
  while (std::getline(input, line))
  {
    if (line.find("DimSize = ") == 0)
    {
      std::istringstream values(line.substr(10));
      for (unsigned int d = 0; d < Dimension; ++d)
      {
        values >> size[d];
      }
    }
    else if (line.find("ElementNumberOfChannels = ") == 0)
    {
      components = std::stoul(line.substr(26));
    }
    else if (line == "ElementDataFile = LOCAL")
    {
      break;
    }
  }
  auto image = ImageType::New();
  image->SetRegions(size);
  image->SetNumberOfComponentsPerPixel(components);
  image->Allocate();
  const auto count = image->GetLargestPossibleRegion().GetNumberOfPixels() * components;
  for (size_t i = 0; i < count; ++i)
  {
    if (!(input >> image->GetBufferPointer()[i]))
    {
      throw std::runtime_error(std::string("Missing fixture values: ") + file);
    }
  }
  input.close();
  auto writer = itk::ImageFileWriter<ImageType>::New();
  writer->SetFileName(file);
  writer->SetInput(image);
  writer->Update();
  return EXIT_SUCCESS;
}

// Compare ImageMath outputs to independently constructed expected images,
// including their geometry and every scalar, RGB, or tensor component.
template <unsigned int Dimension>
int
Check(const char * actualFile, const char * expectedFile, double tolerance)
{
  using ImageType = itk::VectorImage<double, Dimension>;
  auto actualReader = itk::ImageFileReader<ImageType>::New();
  auto expectedReader = itk::ImageFileReader<ImageType>::New();
  actualReader->SetFileName(actualFile);
  expectedReader->SetFileName(expectedFile);
  actualReader->Update();
  expectedReader->Update();
  const auto * actual = actualReader->GetOutput();
  const auto * expected = expectedReader->GetOutput();
  if (!actual->IsSameImageGeometryAs(expected) ||
      actual->GetNumberOfComponentsPerPixel() != expected->GetNumberOfComponentsPerPixel())
  {
    std::cerr << "Image geometry or component count differs: " << actualFile << std::endl;
    return EXIT_FAILURE;
  }
  itk::ImageRegionConstIterator<ImageType> a(actual, actual->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<ImageType> e(expected, expected->GetLargestPossibleRegion());
  for (a.GoToBegin(), e.GoToBegin(); !a.IsAtEnd(); ++a, ++e)
  {
    for (unsigned int component = 0; component < actual->GetNumberOfComponentsPerPixel(); ++component)
    {
      const double value = a.Get()[component];
      const double target = e.Get()[component];
      if (!std::isfinite(value) || !std::isfinite(target) || std::abs(value - target) > tolerance)
      {
        std::cerr << actualFile << " at " << a.GetIndex() << ", component " << component
                  << ": expected " << target << ", got " << value << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}

int
main(int argc, char ** argv)
{
  try
  {
    if (argc == 4 && std::string(argv[1]) == "--write-fixture")
    {
      switch (std::stoi(argv[2]))
      {
        case 2: return WriteFixture<2>(argv[3]);
        case 3: return WriteFixture<3>(argv[3]);
        case 4: return WriteFixture<4>(argv[3]);
        default: return EXIT_FAILURE;
      }
    }
    if (argc != 5)
    {
      std::cerr << "Usage: dimension actualImage expectedImage tolerance" << std::endl;
      return EXIT_FAILURE;
    }
    switch (std::stoi(argv[1]))
    {
      case 2:
        return Check<2>(argv[2], argv[3], std::stod(argv[4]));
      case 3:
        return Check<3>(argv[2], argv[3], std::stod(argv[4]));
      case 4:
        return Check<4>(argv[2], argv[3], std::stod(argv[4]));
      default:
        return EXIT_FAILURE;
    }
  }
  catch (const std::exception & error)
  {
    std::cerr << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
