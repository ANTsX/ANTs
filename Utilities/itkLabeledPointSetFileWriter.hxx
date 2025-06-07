/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkLabeledPointSetFileWriter_hxx
#define itkLabeledPointSetFileWriter_hxx


#include "itkBoundingBox.h"
#include "itkImageFileWriter.h"

#include <fstream>

namespace itk
{
//
// Constructor
//
template <typename TInputMesh>
LabeledPointSetFileWriter<TInputMesh>::LabeledPointSetFileWriter()
{
  this->m_Input = nullptr;
  this->m_FileName = "";
  this->m_MultiComponentScalars = nullptr;
  this->m_Lines = nullptr;

  this->m_ImageSize.Fill(0);
}

//
// Destructor
//
template <typename TInputMesh>
LabeledPointSetFileWriter<TInputMesh>::~LabeledPointSetFileWriter() = default;

//
// Set the input mesh
//
template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::Write()
{
  this->GenerateData();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::GenerateData()
{
  if (this->m_FileName == "")
  {
    itkExceptionMacro("No FileName");
    return;
  }

  if (this->m_ImageSize[0] == 0)
  {
    typedef BoundingBox<unsigned long,
                        Dimension,
                        typename TInputMesh::CoordRepType,
                        typename TInputMesh::PointsContainer>
                                      BoundingBoxType;
    typename BoundingBoxType::Pointer bbox = BoundingBoxType::New();
    bbox->SetPoints(this->m_Input->GetPoints());
    bbox->ComputeBoundingBox();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      this->m_ImageSpacing[d] = static_cast<double>(bbox->GetMaximum()[d] - bbox->GetMinimum()[d]) /
                                static_cast<double>(this->m_ImageSize[d] + 1);
    }
    this->m_ImageDirection.SetIdentity();
  }

  //
  // Read output file
  //
  std::ofstream outputFile(m_FileName.c_str());

  if (!outputFile.is_open())
  {
    itkExceptionMacro("Unable to open file\n"
                      "outputFilename= "
                      << m_FileName);
    return;
  }
  else
  {
    outputFile.close();
  }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind(".");
  std::string            extension(this->m_FileName, pos + 1, this->m_FileName.length() - 1);

  if (extension == "txt")
  {
    this->WritePointsToAvantsFile();
  }
  else if (extension == "vtk")
  {
    this->WriteVTKFile();
  }
  else
  {
    try
    {
      this->WritePointsToImageFile();
    }
    catch (...)
    {
      itkExceptionMacro("Unknown extension: " << extension);
    }
  }
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WriteVTKFile()
{
  this->WritePointsToVTKFile();
  this->WriteLinesToVTKFile();
  this->WriteScalarsToVTKFile();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WritePointsToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile(this->m_FileName.c_str());

  outputFile << "#vtk DataFile Version 2.0" << std::endl;
  outputFile << "File written by itkLabeledPointSetFileWriter" << std::endl;
  outputFile << "ASCII" << std::endl;
  outputFile << "DATASET POLYDATA" << std::endl;

  // POINTS go first

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();
  outputFile << "POINTS " << numberOfPoints << " float" << std::endl;

  typename InputMeshType::PointsContainerIterator pointIterator = this->m_Input->GetPoints()->Begin();
  typename InputMeshType::PointsContainerIterator pointEnd = this->m_Input->GetPoints()->End();
  while (pointIterator != pointEnd)
  {
    PointType point = pointIterator.Value();
    outputFile << point[0] << " " << point[1];
    if (Dimension == 2)
    {
      outputFile << " 0 " << std::endl;
    }
    else if (Dimension == 3)
    {
      outputFile << " " << point[2] << " " << std::endl;
    }
    pointIterator++;
  }

  outputFile.close();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WriteScalarsToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile(this->m_FileName.c_str(), std::ios::app);

  // No point data conditions
  if (!this->m_Input->GetPointData())
  {
    return;
  }
  if (this->m_Input->GetPointData()->Size() == 0)
  {
    return;
  }

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();

  outputFile << std::endl;
  outputFile << "POINT_DATA " << numberOfPoints << std::endl;

  std::string type = std::string("float");

  if (!this->m_MultiComponentScalars)
  {
    outputFile << "SCALARS pointLabels " << type << " 1" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;

    typename InputMeshType::PointDataContainerIterator pointDataIterator = this->m_Input->GetPointData()->Begin();
    typename InputMeshType::PointDataContainerIterator pointDataEnd = this->m_Input->GetPointData()->End();

    while (pointDataIterator != pointDataEnd)
    {
      outputFile << pointDataIterator.Value() << " ";
      pointDataIterator++;
    }

    outputFile << std::endl;
  }
  else
  {
    MultiComponentScalarType scalar = this->m_MultiComponentScalars->GetElement(0);
    unsigned int             numberOfComponents = scalar.GetSize();

    outputFile << "SCALARS scalars " << type << numberOfComponents << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;

    typename MultiComponentScalarSetType::Iterator It = this->m_MultiComponentScalars->Begin();
    typename MultiComponentScalarSetType::Iterator ItEnd = this->m_MultiComponentScalars->End();

    while (It != ItEnd)
    {
      outputFile << It.Value() << " ";
      ++It;
    }

    outputFile << std::endl;
  }
  outputFile.close();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WriteLinesToVTKFile()
{
  if (this->m_Lines)
  {
    std::ofstream outputFile(this->m_FileName.c_str(), std::ios::app);

    unsigned int numberOfLines = this->m_Lines->Size();
    unsigned int totalSize = 0;

    typename LineSetType::Iterator It = this->m_Lines->Begin();
    typename LineSetType::Iterator ItEnd = this->m_Lines->End();

    while (It != ItEnd)
    {
      totalSize += (It.Value()).Size();
      totalSize++;
      ++It;
    }

    outputFile << "LINES " << numberOfLines << " " << totalSize << std::endl;

    It = this->m_Lines->Begin();
    while (It != ItEnd)
    {
      unsigned int numberOfPoints = (It.Value()).Size();
      outputFile << numberOfPoints << " ";
      for (unsigned int d = 0; d < numberOfPoints; d++)
      {
        outputFile << (It.Value())[d] << " ";
      }
      outputFile << std::endl;
      ++It;
    }

    outputFile << std::endl;
    outputFile.close();
  }
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WritePointsToAvantsFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile(this->m_FileName.c_str());

  outputFile << "0 0 0 0" << std::endl;

  if (this->m_Input->GetNumberOfPoints() > 0)
  {
    typename InputMeshType::PointsContainerIterator pointIterator = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator = this->m_Input->GetPointData()->Begin();

    while (pointIterator != pointEnd)
    {
      PointType point = pointIterator.Value();
      outputFile << point[0] << " " << point[1];
      if (Dimension == 2)
      {
        outputFile << " 0 ";
      }
      else if (Dimension == 3)
      {
        outputFile << " " << point[2] << " ";
      }
      outputFile << pointDataIterator.Value() << std::endl;
      pointIterator++;
      pointDataIterator++;
    }
  }

  outputFile << "0 0 0 0" << std::endl;

  outputFile.close();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::WritePointsToImageFile()
{
  typename LabeledPointSetImageType::Pointer outputImage = LabeledPointSetImageType::New();
  outputImage->SetDirection(this->m_ImageDirection);
  outputImage->SetRegions(this->m_ImageSize);
  outputImage->SetOrigin(this->m_ImageOrigin);
  outputImage->SetSpacing(this->m_ImageSpacing);
  outputImage->AllocateInitialized();

  if (this->m_Input->GetNumberOfPoints() > 0)
  {
    typename InputMeshType::PointsContainerIterator pointIterator = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator = this->m_Input->GetPointData()->Begin();

    while (pointIterator != pointEnd)
    {
      PointType point = pointIterator.Value();
      PixelType label = pointDataIterator.Value();

      typename LabeledPointSetImageType::IndexType index;
      typename LabeledPointSetImageType::PointType ipoint;
      ipoint.CastFrom(point);
      if (outputImage->TransformPhysicalPointToIndex(ipoint, index))
      {
        outputImage->SetPixel(index, label);
      }
      pointIterator++;
      pointDataIterator++;
    }
  }

  typedef ImageFileWriter<LabeledPointSetImageType> WriterType;
  typename WriterType::Pointer                      writer = WriterType::New();
  writer->SetFileName(this->m_FileName.c_str());
  writer->SetInput(outputImage);
  writer->Update();
}

template <typename TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}
} // end of namespace itk

#endif
