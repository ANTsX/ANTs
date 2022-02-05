/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkLabeledPointSetFileReader_hxx
#define itkLabeledPointSetFileReader_hxx


#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkByteSwapper.h"

#include <fstream>
#include <cstdio>
#include <string>

namespace itk
{
//
// Constructor
//
template <typename TOutputMesh>
LabeledPointSetFileReader<TOutputMesh>::LabeledPointSetFileReader()
  : m_ExtractBoundaryPoints(false)
  , m_FileName()
  , m_RandomPercentage(1.0)
  , m_LabelSet()
  , m_MultiComponentScalars(nullptr)
  , m_Lines(nullptr)
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::GenerateData()
{
  if (this->m_FileName == "")
  {
    itkExceptionMacro("No input FileName");
    return;
  }

  //
  // Read input file
  //
  std::ifstream inputFile(m_FileName.c_str());

  if (!inputFile.is_open())
  {
    itkExceptionMacro("Unable to open file\n"
                      "inputFilename= "
                      << m_FileName);
    return;
  }
  else
  {
    inputFile.close();
  }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind(".");
  std::string            extension(this->m_FileName, pos + 1, this->m_FileName.length() - 1);

  if (extension == "txt")
  {
    this->ReadPointsFromAvantsFile();
  }
  else if (extension == "vtk")
  {
    this->ReadVTKFile();
  }
  else // try reading the file as an image
  {
    this->ReadPointsFromImageFile();
  }

  if (this->m_RandomPercentage < 1.0)
  {
    typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    typename GeneratorType::Pointer                           generator = GeneratorType::New();
    generator->SetSeed();

    typename OutputMeshType::Pointer output = OutputMeshType::New();
    output->Initialize();

    if (this->GetOutput()->GetNumberOfPoints() > 0)
    {
      typename OutputMeshType::PointsContainerIterator It = this->GetOutput()->GetPoints()->Begin();

      unsigned long count = 0;
      while (It != this->GetOutput()->GetPoints()->End())
      {
        if (generator->GetVariateWithClosedRange() <= this->m_RandomPercentage)
        {
          output->SetPoint(count, It.Value());
          PixelType label;
          bool      elementExists = this->GetOutput()->GetPointData()->GetElementIfIndexExists(It.Index(), &label);
          if (elementExists)
          {
            output->SetPointData(count, label);
          }
          count++;
        }
        ++It;
      }
    }
    this->GraftOutput(output);
  }
  this->m_LabelSet.clear();

  /**
   * If the number of points does not match the number of
   * point data, fill the point data with zeros
   */
  if (!this->GetOutput()->GetPointData() ||
      this->GetOutput()->GetPointData()->Size() != this->GetOutput()->GetPoints()->Size())
  {
    itkWarningMacro("Number of points does not match number of labels. "
                    << "Filling point data with label zero.");
    typename OutputMeshType::PointsContainerIterator It = this->GetOutput()->GetPoints()->Begin();

    while (It != this->GetOutput()->GetPoints()->End())
    {
      this->GetOutput()->SetPointData(It.Index(), NumericTraits<PixelType>::ZeroValue());
      ++It;
    }
  }

  if (this->GetOutput()->GetNumberOfPoints() > 0)
  {
    typename OutputMeshType::PointDataContainerIterator ItD = this->GetOutput()->GetPointData()->Begin();
    while (ItD != this->GetOutput()->GetPointData()->End())
    {
      if (find(this->m_LabelSet.begin(), this->m_LabelSet.end(), ItD.Value()) == this->m_LabelSet.end())
      {
        this->m_LabelSet.push_back(ItD.Value());
      }
      ++ItD;
    }
  }
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadPointsFromAvantsFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile(m_FileName.c_str());

  unsigned long count = 0;
  while (!inputFile.eof())
  {
    PointType point;
    PixelType label;

    if (Dimension == 2)
    {
      float trash;
      inputFile >> point >> trash >> label;
    }
    else // Dimension == 3
    {
      inputFile >> point >> label;
    }
    if ((point.GetVectorFromOrigin()).GetSquaredNorm() > 0.0 || label != 0)
    {
      outputMesh->SetPointData(count, label);
      outputMesh->SetPoint(count, point);
      count++;
    }
  }

  inputFile.close();
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadVTKFile()
{
  this->ReadPointsFromVTKFile();
  this->ReadScalarsFromVTKFile();
  this->ReadLinesFromVTKFile();
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadPointsFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile(this->m_FileName.c_str());

  std::string line;

  bool isBinary = false;

  while (!inputFile.eof())
  {
    std::getline(inputFile, line);

    if (line.find("BINARY") != std::string::npos)
    {
      isBinary = true;
    }

    if (line.find("POINTS") != std::string::npos)
    {
      break;
    }
  }

  itkDebugMacro("POINTS line" << line);

  std::string pointLine(line, strlen("POINTS "), line.length());
  itkDebugMacro("pointLine " << pointLine);

  int numberOfPoints = -1;

  if (sscanf(pointLine.c_str(), "%d", &numberOfPoints) != 1)
  {
    itkExceptionMacro("ERROR: Failed to read numberOfPoints\n"
                      "       pointLine = "
                      << pointLine);
    return;
  }

  itkDebugMacro("numberOfPoints = " << numberOfPoints);

  if (numberOfPoints < 1)
  {
    itkExceptionMacro("numberOfPoints < 1"
                      << "       numberOfPoints = " << numberOfPoints);
    return;
  }

  outputMesh->GetPoints()->Reserve(numberOfPoints);

  //
  // Load the point coordinates into the itk::Mesh
  //
  PointType point;

  if (isBinary)
  {
    itkDebugMacro("Data is binary");

    float * ptData = new float[numberOfPoints * 3];
    inputFile.read(reinterpret_cast<char *>(ptData), static_cast<std::size_t>(3) * numberOfPoints * sizeof(float));
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData, numberOfPoints * static_cast<std::size_t>(3));
    for (long i = 0; i < numberOfPoints; i++)
    {
      for (long j = 0; j < Dimension; j++)
      {
        point[j] = ptData[i * 3 + j];
      }
      outputMesh->SetPoint(i, point);
    }

    delete[] ptData;
  }
  else
  {
    for (long i = 0; i < numberOfPoints; i++)
    {
      if (Dimension == 2)
      {
        float trash;
        inputFile >> point >> trash;
      }
      else // Dimension = 3
      {
        inputFile >> point;
      }
      outputMesh->SetPoint(i, point);
    }
  }

  inputFile.close();
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadScalarsFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile(this->m_FileName.c_str());

  std::string line;

  bool isBinary = false;

  //
  // Find the labels associated with each pixel
  //
  while (!inputFile.eof())
  {
    std::getline(inputFile, line);

    if (line.find("BINARY") != std::string::npos)
    {
      isBinary = true;
    }

    if (line.find("SCALARS") != std::string::npos)
    {
      break;
    }
  }

  if (inputFile.eof())
  {
    inputFile.close();
    return;
  }

  std::string::size_type pos = line.rfind(" ");

  std::string temp = std::string(line, pos + 1, line.length() - 1);

  unsigned int numberOfComponents = static_cast<unsigned int>(std::atoi(temp.c_str()));

  std::getline(inputFile, line);

  if (isBinary)
  {
    int   numberOfValues = outputMesh->GetNumberOfPoints() * numberOfComponents;
    float p;
    int * scalarData = new int[numberOfValues];
    inputFile.read(reinterpret_cast<char *>(scalarData), numberOfComponents * sizeof(p));
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(scalarData, numberOfValues);

    if (numberOfComponents == 1)
    {
      // PixelType label;
      for (unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++)
      {
        outputMesh->SetPointData(i, scalarData[i]);
      }
      //    itkExceptionMacro( "Only single label components are readable" );
    }
    else
    {
      this->m_MultiComponentScalars = MultiComponentScalarSetType::New();
      this->m_MultiComponentScalars->Initialize();
      for (unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++)
      {
        MultiComponentScalarType scalar;
        scalar.SetSize(numberOfComponents);
        for (unsigned int d = 0; d < numberOfComponents; d++)
        {
          scalar[d] = scalarData[i * numberOfComponents + d];
        }
        this->m_MultiComponentScalars->InsertElement(i, scalar);
      }
    }

    delete[] scalarData;
  }
  else
  {
    if (numberOfComponents == 1)
    {
      PixelType label;
      for (unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++)
      {
        inputFile >> label;
        outputMesh->SetPointData(i, label);
      }
      //    itkExceptionMacro( "Only single label components are readable" );
    }
    else
    {
      this->m_MultiComponentScalars = MultiComponentScalarSetType::New();
      this->m_MultiComponentScalars->Initialize();
      for (unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++)
      {
        MultiComponentScalarType scalar;
        scalar.SetSize(numberOfComponents);
        for (unsigned int d = 0; d < numberOfComponents; d++)
        {
          inputFile >> scalar[d];
        }
        this->m_MultiComponentScalars->InsertElement(i, scalar);
      }
    }
  }

  inputFile.close();
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadLinesFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile(this->m_FileName.c_str());

  std::string line;

  bool isBinary = false;

  //
  // Find the labels associated with each pixel
  //
  while (!inputFile.eof())
  {
    std::getline(inputFile, line);

    if (line.find("BINARY") != std::string::npos)
    {
      isBinary = true;
    }

    if (line.find("LINES") != std::string::npos)
    {
      break;
    }
  }

  if (inputFile.eof())
  {
    inputFile.close();
    return;
  }

  std::string::size_type pos = line.rfind(" ");

  std::string  temp = std::string(line, 6, pos - 1);
  unsigned int numberOfLines = static_cast<unsigned int>(std::atoi(temp.c_str()));

  temp = std::string(line, pos, line.length() - 1);
  unsigned int numberOfValues = static_cast<unsigned int>(std::atoi(temp.c_str()));

  this->m_Lines = LineSetType::New();
  this->m_Lines->Initialize();

  if (isBinary)
  {
    int   p;
    int * lineData = new int[numberOfValues];
    inputFile.read(reinterpret_cast<char *>(lineData), numberOfValues * sizeof(p));
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(lineData, numberOfValues);

    unsigned long valueId = 0;
    unsigned long lineId = 0;
    while (valueId < numberOfValues)
    {
      itk::SizeValueType lineLength = lineData[valueId];
      ++valueId;

      LineType polyLine;
      polyLine.SetSize(static_cast<itk::SizeValueType>(lineLength));
      for (itk::SizeValueType i = 0; i < lineLength; i++)
      {
        polyLine[i] = lineData[valueId];
        ++valueId;
      }
      this->m_Lines->InsertElement(lineId, polyLine);
      ++lineId;
    }

    delete[] lineData;
  }
  else
  {
    for (unsigned int i = 0; i < numberOfLines; i++)
    {
      LineType     line_loc;
      unsigned int numberOfPoints;
      inputFile >> numberOfPoints;
      line_loc.SetSize(numberOfPoints);
      for (unsigned int d = 0; d < numberOfPoints; d++)
      {
        inputFile >> line_loc[d];
      }
      this->m_Lines->InsertElement(i, line_loc);
    }
  }

  inputFile.close();
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::ReadPointsFromImageFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  typedef ImageFileReader<LabeledPointSetImageType> ImageReaderType;
  typename ImageReaderType::Pointer                 imageReader = ImageReaderType::New();
  imageReader->SetFileName(this->m_FileName.c_str());
  imageReader->Update();

  if (!this->m_ExtractBoundaryPoints)
  {
    ImageRegionIteratorWithIndex<LabeledPointSetImageType> It(imageReader->GetOutput(),
                                                              imageReader->GetOutput()->GetLargestPossibleRegion());

    unsigned long count = 0;
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      PixelType label = It.Get();
      if (label != NumericTraits<PixelType>::ZeroValue())
      {
        typename LabeledPointSetImageType::PointType imagePoint;
        imageReader->GetOutput()->TransformIndexToPhysicalPoint(It.GetIndex(), imagePoint);

        PointType point;
        point.CastFrom(imagePoint);
        outputMesh->SetPoint(count, point);
        outputMesh->SetPointData(count, label);
        count++;
      }
    }
  }
  else
  {
    typedef LabelContourImageFilter<LabeledPointSetImageType, LabeledPointSetImageType> ContourFilterType;
    typename ContourFilterType::Pointer contourFilter = ContourFilterType::New();
    contourFilter->SetInput(imageReader->GetOutput());
    contourFilter->SetFullyConnected(true);
    contourFilter->SetBackgroundValue(0);
    contourFilter->Update();

    this->m_LabelSet.clear();

    ImageRegionIteratorWithIndex<LabeledPointSetImageType> It(contourFilter->GetOutput(),
                                                              contourFilter->GetOutput()->GetLargestPossibleRegion());
    unsigned long                                          count = 0;
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      PixelType label = It.Get();
      if (It.Get() > 0)
      {
        typename LabeledPointSetImageType::PointType imagePoint;
        contourFilter->GetOutput()->TransformIndexToPhysicalPoint(It.GetIndex(), imagePoint);

        PointType point;
        point.CastFrom(imagePoint);
        outputMesh->SetPoint(count, point);
        outputMesh->SetPointData(count, It.Get());
        count++;

        if (find(this->m_LabelSet.begin(), this->m_LabelSet.end(), label) == this->m_LabelSet.end())
        {
          this->m_LabelSet.push_back(label);
        }
      }
    }
  }
}

template <typename TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
  os << indent << "RandomPercentage: " << this->m_RandomPercentage << std::endl;
}
} // end of namespace itk

#endif
