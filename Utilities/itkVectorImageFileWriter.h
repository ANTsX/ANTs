/*=========================================================================

  Program:   Advanced Normalization Tools
  Date:      $$
  Version:   $ $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorImageFileWriter_h
#define __itkVectorImageFileWriter_h

#include "itkProcessObject.h"
#include "itkImageIOBase.h"
#include "itkExceptionObject.h"
#include "itkSize.h"
#include "itkImageIORegion.h"

namespace itk
{
/** \brief Base exception class for IO problems during writing. */
class VectorImageFileWriterException : public ExceptionObject
{
public:
  /** Run-time information. */
  itkOverrideGetNameOfClassMacro(VectorImageFileWriterException);

  /** Constructor. */
  VectorImageFileWriterException(const char * file,
                                 unsigned int line,
                                 const char * message = "Error in IO",
                                 const char * loc = "Unknown")
    : ExceptionObject(file, line, message, loc)
  {}

  /** Constructor. */
  VectorImageFileWriterException(const std::string & file,
                                 unsigned int        line,
                                 const char *        message = "Error in IO",
                                 const char *        loc = "Unknown")
    : ExceptionObject(file, line, message, loc)
  {}
};

/** \class VectorImageFileWriter
 * \brief Writes the deformation field as component images files.
 *
 * \sa VectorImageFileWriter
 * \sa ImageSeriesReader
 * \sa ImageIOBase
 *
 * \ingroup IOFilters
 */
template <typename TVectorImage, typename TImage>
class VectorImageFileWriter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef VectorImageFileWriter    Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(VectorImageFileWriter);

  /** Some convenient typedefs. */
  typedef TVectorImage                         VectorImageType;
  typedef typename VectorImageType::Pointer    VectorImagePointer;
  typedef typename VectorImageType::RegionType VectorImageRegionType;
  typedef typename VectorImageType::PixelType  VectorImagePixelType;
  typedef TImage                               ImageType;
  typedef typename ImageType::Pointer          ImagePointer;
  typedef typename ImageType::RegionType       ImageRegionType;
  typedef typename ImageType::PixelType        ImagePixelType;

  /** Set/Get the image input of this writer.  */
  void
  SetInput(const VectorImageType * input);

  const VectorImageType *
  GetInput();

  const VectorImageType *
  GetInput(unsigned int idx);

  /** Specify the name of the output file to write. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  /** Set/Get the ImageIO helper class. Usually this is created via the object
   * factory mechanism that determines whether a particular ImageIO can
   * write a certain file. This method provides a way to get the ImageIO
   * instance that is created, or one can be manually set where the
   * IO factory mechanism may not work (for example, raw image files or
   * image files with non-standard filename suffix's.
   * If the user specifies the ImageIO, we assume she makes the
   * correct choice and will allow a file to be created regardless of
   * the file extension. If the factory has set the ImageIO, the
   * extension must be supported by the specified ImageIO. */
  void
  SetImageIO(ImageIOBase * io)
  {
    if (this->m_ImageIO != io)
    {
      this->Modified();
      this->m_ImageIO = io;
    }
    m_FactorySpecifiedImageIO = false;
  }

  itkGetModifiableObjectMacro(ImageIO, ImageIOBase);

  /** A special version of the Update() method for writers.  It
   * invokes start and end events and handles releasing data. It
   * eventually calls GenerateData() which does the actual writing.
   * Note: the write method will write data specified by the
   * IORegion. If not set, then then the whole image is written.  Note
   * that the region will be cropped to fit the input image's
   * LargestPossibleRegion. */
  virtual void
  Write();

  /** Specify the region to write. If left NULL, then the whole image
   * is written. */
  void
  SetIORegion(const ImageIORegion & region);

  itkGetConstReferenceMacro(IORegion, ImageIORegion);

  /** Aliased to the Write() method to be consistent with the rest of the
   * pipeline. */
  virtual void
  Update()
  {
    this->Write();
  }

  /** Set the compression On or Off */
  itkSetMacro(UseCompression, bool);
  itkGetConstReferenceMacro(UseCompression, bool);
  itkBooleanMacro(UseCompression);

  /** Set the Avants' naming convention On or Off */
  itkSetMacro(UseAvantsNamingConvention, bool);
  itkGetConstReferenceMacro(UseAvantsNamingConvention, bool);
  itkBooleanMacro(UseAvantsNamingConvention);

  /** Set the Avants' naming convention On or Off */
  itkSetMacro(UseZhangNamingConvention, bool);
  itkGetConstReferenceMacro(UseZhangNamingConvention, bool);
  itkBooleanMacro(UseZhangNamingConvention);

  /** By default the MetaDataDictionary is taken from the input image and
   *  passed to the ImageIO. In some cases, however, a user may prefer to
   *  introduce her/his own MetaDataDictionary. This is often the case of
   *  the ImageSeriesWriter. This flag defined whether the MetaDataDictionary
   *  to use will be the one from the input image or the one already set in
   *  the ImageIO object. */
  itkSetMacro(UseInputMetaDataDictionary, bool);
  itkGetConstReferenceMacro(UseInputMetaDataDictionary, bool);
  itkBooleanMacro(UseInputMetaDataDictionary);

protected:
  VectorImageFileWriter();
  ~VectorImageFileWriter();
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** Does the real work. */
  void
  GenerateData();

private:
  VectorImageFileWriter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  std::string m_FileName;
  std::string m_ComponentImageFileName;
  bool        m_UseAvantsNamingConvention;
  bool        m_UseZhangNamingConvention;

  ImagePointer m_Image;

  ImageIOBase::Pointer m_ImageIO;
  bool                 m_UserSpecifiedImageIO; // track whether the ImageIO is user specified

  ImageIORegion m_IORegion;
  bool          m_UserSpecifiedIORegion; //
                                         // track whether the region is user specified
  bool m_FactorySpecifiedImageIO;        // track whether the factory mechanism set the ImageIO
  bool m_UseCompression;
  bool m_UseInputMetaDataDictionary; // whether to use the MetaDataDictionary from the input or not.
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkVectorImageFileWriter.hxx"
#endif

#endif // __itkVectorImageFileWriter_h
