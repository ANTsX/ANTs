/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkLabeledPointSetFileWriter_h
#define itkLabeledPointSetFileWriter_h

#include "itkMesh.h"

#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

namespace itk
{
/** \class LabeledPointSetFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <typename TInputMesh>
class LabeledPointSetFileWriter final : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef LabeledPointSetFileWriter Self;
  typedef Object                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void
  Update();

  void
  Write();

  /** Extract dimension from the output mesh. */
  static constexpr unsigned int Dimension = TInputMesh::PointType::Dimension;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(LabeledPointSetFileWriter);

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                                          InputMeshType;
  typedef typename TInputMesh::Pointer                        InputMeshPointer;
  typedef typename InputMeshType::MeshTraits                  MeshTraits;
  typedef typename InputMeshType::Superclass                  PointSetType;
  typedef typename InputMeshType::PointType                   PointType;
  typedef typename MeshTraits::PixelType                      PixelType;
  typedef Array<PixelType>                                    MultiComponentScalarType;
  typedef Array<unsigned long>                                LineType;
  typedef VectorContainer<long, MultiComponentScalarType>     MultiComponentScalarSetType;
  typedef VectorContainer<long, LineType>                     LineSetType;
  typedef Image<PixelType, Self::Dimension> LabeledPointSetImageType;
  typedef typename LabeledPointSetImageType::SizeType         ImageSizeType;
  typedef typename LabeledPointSetImageType::PointType        ImageOriginType;
  typedef typename LabeledPointSetImageType::SpacingType      ImageSpacingType;
  typedef typename LabeledPointSetImageType::DirectionType    ImageDirectionType;

  /** Set the Input */
  void
  SetInput(InputMeshType * input);

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  /** Specify other attributes */
  itkSetMacro(Lines, typename LineSetType::Pointer);

  itkSetMacro(MultiComponentScalars, typename MultiComponentScalarSetType::Pointer);

  /** Specify image attributes if output is an image. */
  itkSetMacro(ImageSize, ImageSizeType);
  itkGetConstMacro(ImageSize, ImageSizeType);

  itkSetMacro(ImageOrigin, ImageOriginType);
  itkGetConstMacro(ImageOrigin, ImageOriginType);

  itkSetMacro(ImageSpacing, ImageSpacingType);
  itkGetConstMacro(ImageSpacing, ImageSpacingType);

  itkSetMacro(ImageDirection, ImageDirectionType);
  itkGetConstMacro(ImageDirection, ImageDirectionType);

protected:
  LabeledPointSetFileWriter();
  ~LabeledPointSetFileWriter() override;

  virtual void
  GenerateData();

  std::string      m_FileName;
  InputMeshPointer m_Input;

  typename MultiComponentScalarSetType::Pointer m_MultiComponentScalars;
  typename LineSetType::Pointer                 m_Lines;

  /**
   * If output is an image type, the attributes must be specified.
   */
  ImageSizeType      m_ImageSize;
  ImageSpacingType   m_ImageSpacing;
  ImageOriginType    m_ImageOrigin;
  ImageDirectionType m_ImageDirection;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  LabeledPointSetFileWriter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  void
  WritePointsToAvantsFile();

  void
  WritePointsToImageFile();

  void
  WriteVTKFile();

  void
  WritePointsToVTKFile();

  void
  WriteScalarsToVTKFile();

  void
  WriteLinesToVTKFile();
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLabeledPointSetFileWriter.hxx"
#endif

#endif
