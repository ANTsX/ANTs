/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkLabeledPointSetFileReader_h
#define itkLabeledPointSetFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"

#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

#include <vector>

namespace itk
{
/** \class LabeledPointSetFileReader
 * \brief
 * Reads a file and creates an itkMesh.
 *
 */
template <typename TOutputMesh>
class LabeledPointSetFileReader final : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef LabeledPointSetFileReader Self;
  typedef MeshSource<TOutputMesh>   Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from the output mesh. */
  static constexpr unsigned int Dimension = TOutputMesh::PointType::Dimension;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(LabeledPointSetFileReader);

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                                     OutputMeshType;
  typedef typename OutputMeshType::MeshTraits             MeshTraits;
  typedef typename OutputMeshType::Superclass             PointSetType;
  typedef typename OutputMeshType::PointType              PointType;
  typedef typename MeshTraits::PixelType                  PixelType;
  typedef Array<PixelType>                                MultiComponentScalarType;
  typedef Array<unsigned long>                            LineType;
  typedef VectorContainer<long, MultiComponentScalarType> MultiComponentScalarSetType;
  typedef VectorContainer<long, LineType>                 LineSetType;

  typedef Image<PixelType, Self::Dimension> LabeledPointSetImageType;

  typedef std::vector<PixelType> LabelSetType;

  /** Set/Get the name of the file to be read. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  itkSetMacro(ExtractBoundaryPoints, bool);
  itkGetMacro(ExtractBoundaryPoints, bool);
  itkBooleanMacro(ExtractBoundaryPoints);

  /**
   * Percentage of points selected randomnly
   */
  itkSetClampMacro(RandomPercentage, double, 0.0, 1.0);
  itkGetConstMacro(RandomPercentage, double);

  LabelSetType *
  GetLabelSet()
  {
    return &this->m_LabelSet;
  }

  unsigned int
  GetNumberOfLabels() const
  {
    return this->m_LabelSet.size();
  }

  MultiComponentScalarSetType *
  GetMultiComponentScalars()
  {
    return this->m_MultiComponentScalars.GetPointer();
  }

  LineSetType *
  GetLines()
  {
    return this->m_Lines.GetPointer();
  }

protected:
  LabeledPointSetFileReader();
  ~LabeledPointSetFileReader() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Reads the file */
  void
  GenerateData() override;

  bool m_ExtractBoundaryPoints;

  std::string                                   m_FileName;
  double                                        m_RandomPercentage;
  LabelSetType                                  m_LabelSet;
  typename MultiComponentScalarSetType::Pointer m_MultiComponentScalars;
  typename LineSetType::Pointer                 m_Lines;

private:
  LabeledPointSetFileReader(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  void
  ReadPointsFromImageFile();

  void
  ReadPointsFromAvantsFile();

  void
  ReadVTKFile();

  void
  ReadPointsFromVTKFile();

  void
  ReadScalarsFromVTKFile();

  void
  ReadLinesFromVTKFile();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLabeledPointSetFileReader.hxx"
#endif

#endif // _itkLabeledPointSetFileReader_h
