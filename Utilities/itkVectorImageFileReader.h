/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorImageFileReader_h
#define __itkVectorImageFileReader_h

#include "itkImageIOBase.h"
#include "itkImageSource.h"
#include "itkExceptionObject.h"
#include "itkSize.h"
#include "itkImageRegion.h"
#include "itkDefaultConvertPixelTraits.h"

namespace itk
{
/** \brief Base exception class for IO conflicts. */
class VectorImageFileReaderException : public ExceptionObject
{
public:
  /** Run-time information. */
  itkOverrideGetNameOfClassMacro(VectorImageFileReaderException);

  /** Constructor. */
  VectorImageFileReaderException(const char * file,
                                 unsigned int line,
                                 const char * message = "Error in IO",
                                 const char * loc = "Unknown")
    : ExceptionObject(file, line, message, loc)
  {}

  /** Constructor. */
  VectorImageFileReaderException(const std::string & file,
                                 unsigned int        line,
                                 const char *        message = "Error in IO",
                                 const char *        loc = "Unknown")
    : ExceptionObject(file, line, message, loc)
  {}
};

/** \brief Data source that reads image data from a single file.
 *
 * This source object is a general filter to read data from
 * a variety of file formats. It works with a ImageIOBase subclass
 * to actually do the reading of the data. Object factory machinery
 * can be used to automatically create the ImageIOBase, or the
 * ImageIOBase can be manually created and set. Note that this
 * class reads data from a single file; if you wish to read data
 * from a series of files use ImageSeriesReader.
 *
 * TImage is the type expected by the external users of the
 * filter. If data stored in the file is stored in a different format
 * then specified by TImage, than this filter converts data
 * between the file type and the external expected type.  The
 * ConvertTraits template argument is used to do the conversion.
 *
 * A Pluggable factory pattern is used this allows different kinds of readers
 * to be registered (even at run time) without having to modify the
 * code in this class. Normally just setting the FileName with the
 * appropriate suffix is enough to get the reader to instantiate the
 * correct ImageIO and read the file properly. However, some files (like
 * raw binary format) have no accepted suffix, so you will have to
 * manually create the ImageIO instance of the write type.
 *
 * \sa ImageSeriesReader
 * \sa ImageIOBase
 *
 * \ingroup IOFilters
 *
 */
template <typename TImage,
          typename TVectorImage,
          typename ConvertPixelTraits = DefaultConvertPixelTraits<ITK_TYPENAME TImage::IOPixelType>>
class VectorImageFileReader : public ImageSource<TVectorImage>
{
public:
  /** Standard class typedefs. */
  typedef VectorImageFileReader     Self;
  typedef ImageSource<TVectorImage> Superclass;
  typedef SmartPointer<Self>        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(VectorImageFileReader);

  /** Image types */
  typedef typename TImage::RegionType        ImageRegionType;
  typedef typename TImage::InternalPixelType ImagePixelType;

  /** Deformation field types */
  typedef typename TVectorImage::RegionType        VectorImageRegionType;
  typedef typename TVectorImage::InternalPixelType VectorImagePixelType;

  /** Specify the file to read. This is forwarded to the IO instance. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  /** Set the Avants' naming convention On or Off */
  itkSetMacro(UseAvantsNamingConvention, bool);
  itkGetConstReferenceMacro(UseAvantsNamingConvention, bool);
  itkBooleanMacro(UseAvantsNamingConvention);

  /** Set/Get the ImageIO helper class. Often this is created via the object
   * factory mechanism that determines whether a particular ImageIO can
   * read a certain file. This method provides a way to get the ImageIO
   * instance that is created. Or you can directly specify the ImageIO
   * to use to read a particular file in case the factory mechanism will
   * not work properly (e.g., unknown or unusual extension). */
  void
  SetImageIO(ImageIOBase * imageIO);

  itkGetModifiableObjectMacro(ImageIO, ImageIOBase);

  /** Prepare the allocation of the output image during the first back
   * propagation of the pipeline. */
  virtual void
  GenerateOutputInformation();

  /** Give the reader a chance to indicate that it will produce more
   * output than it was requested to produce. VectorImageFileReader cannot
   * currently read a portion of an image (since the ImageIO objects
   * cannot read a portion of an image), so the VectorImageFileReader must
   * enlarge the RequestedRegion to the size of the image on disk. */
  virtual void
  EnlargeOutputRequestedRegion(DataObject * output);

protected:
  VectorImageFileReader();
  ~VectorImageFileReader();
  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** Convert a block of pixels from one type to another. */
  void
  DoConvertBuffer(void * buffer, unsigned long numberOfPixels);

  /** Test whether the given filename exist and it is readable,
      this is intended to be called before attempting to use
      ImageIO classes for actually reading the file. If the file
      doesn't exist or it is not readable, and exception with an
      approriate message will be thrown. */
  void
  TestFileExistanceAndReadability();

  /** Does the real work. */
  virtual void
  GenerateData();

  ImageIOBase::Pointer m_ImageIO;
  bool                 m_UserSpecifiedImageIO; // keep track whether the ImageIO is user specified

  /** The file to be read. */
  std::string m_FileName;

private:
  VectorImageFileReader(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  std::string m_ExceptionMessage;

  typename TImage::Pointer m_Image;
  bool                     m_UseAvantsNamingConvention;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkVectorImageFileReader.hxx"
#endif

#endif // __itkVectorImageFileReader_h
