/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkVectorImageFileWriter_hxx
#define _itkVectorImageFileWriter_hxx

#include "itkImageFileWriter.h"
#include "itkDataObject.h"
#include "itkObjectFactoryBase.h"
#include "itkImageIOFactory.h"
#include "itkCommand.h"
#include "vnl/vnl_vector.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk
{
// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
VectorImageFileWriter<TVectorImage, TImage>::VectorImageFileWriter()
  : m_FileName("")
  , m_ImageIO(0)
  , m_UserSpecifiedImageIO(false)
  , m_UserSpecifiedIORegion(false)
{
  m_UseCompression = false;
  m_UseInputMetaDataDictionary = true;
  m_FactorySpecifiedImageIO = false;
  m_UseAvantsNamingConvention = true;
  m_UseZhangNamingConvention = false;
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
VectorImageFileWriter<TVectorImage, TImage>::~VectorImageFileWriter()
{}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
void
VectorImageFileWriter<TVectorImage, TImage>::SetInput(const VectorImageType * input)
{
  this->ProcessObject::SetNthInput(0, const_cast<TVectorImage *>(input));
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
const typename VectorImageFileWriter<TVectorImage, TImage>::VectorImageType *
VectorImageFileWriter<TVectorImage, TImage>::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
  {
    return 0;
  }

  return static_cast<TVectorImage *>(this->ProcessObject::GetInput(0));
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
const typename VectorImageFileWriter<TVectorImage, TImage>::VectorImageType *
VectorImageFileWriter<TVectorImage, TImage>::GetInput(unsigned int idx)
{
  return static_cast<TVectorImage *>(this->ProcessObject::GetInput(idx));
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
void
VectorImageFileWriter<TVectorImage, TImage>::SetIORegion(const ImageIORegion & region)
{
  itkDebugMacro("setting IORegion to " << region);
  if (m_IORegion != region)
  {
    m_IORegion = region;
    this->Modified();
    m_UserSpecifiedIORegion = true;
  }
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
void
VectorImageFileWriter<TVectorImage, TImage>::GenerateData(void)
{
  itkDebugMacro(<< "Writing file: " << m_ComponentImageFileName);

  // Make sure that the image is the right type and no more than
  // four components.
  typedef typename ImageType::PixelType ScalarType;

  if (strcmp(m_Image->GetNameOfClass(), "VectorImage") == 0)
  {
    typedef typename ImageType::InternalPixelType VectorImageScalarType;
    // m_ImageIO->SetPixelTypeInfo( typeid(VectorImageScalarType) );

    typedef typename ImageType::AccessorFunctorType AccessorFunctorType;
    m_ImageIO->SetNumberOfComponents(AccessorFunctorType::GetVectorLength(m_Image));
  }
  else
  {
    // Set the pixel and component type; the number of components.
    // m_ImageIO->SetPixelTypeInfo(typeid(ScalarType));
  }

  // Setup the image IO for writing.
  //
  m_ImageIO->SetFileName(m_ComponentImageFileName.c_str());

  // okay, now extract the data as a raw buffer pointer
  const void * dataPtr = (const void *)this->m_Image->GetBufferPointer();

  m_ImageIO->Write(dataPtr);
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
void
VectorImageFileWriter<TVectorImage, TImage>::Write()
{
  typedef VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> SelectorType;
  typename SelectorType::Pointer                                          selector = SelectorType::New();
  selector->SetInput(this->GetInput());

  unsigned int dimension = VectorImageType::PixelType::VectorDimension;

  std::string            filename = this->m_FileName;
  std::string::size_type pos = this->m_FileName.rfind(".");
  std::string            extension(this->m_FileName, pos, this->m_FileName.length() - 1);

  std::string gzExtension("");
  if (extension == std::string(".gz"))
  {
    gzExtension = extension;
    filename = std::string(filename, 0, pos);
    pos = filename.rfind(".");
    extension = std::string(filename, pos, this->m_FileName.length() - 1);
  }
  for (unsigned int i = 0; i < dimension; i++)
  {
    selector->SetIndex(i);
    selector->Update();

    std::string filename(this->m_FileName, 0, pos);

    if (this->m_UseAvantsNamingConvention && dimension <= 3)
    {
      switch (i)
      {
        case 0:
        {
          filename += std::string("xvec");
        }
        break;
        case 1:
        {
          filename += std::string("yvec");
        }
        break;
        case 2:
        {
          filename += std::string("zvec");
        }
        break;
        default:
        {
          filename += std::string("you_are_screwed_vec");
        }
        break;
      }
    }
    else if (this->m_UseZhangNamingConvention && dimension == 6)
    {
      switch (i)
      {
        case 0:
        {
          filename += std::string("xx");
        }
        break;
        case 1:
        {
          filename += std::string("yx");
        }
        break;
        case 2:
        {
          filename += std::string("yy");
        }
        break;
        case 3:
        {
          filename += std::string("zx");
        }
        break;
        case 4:
        {
          filename += std::string("zy");
        }
        break;
        case 5:
        {
          filename += std::string("zz");
        }
        break;
        default:
        {
          filename += std::string("you_are_screwed");
        }
        break;
      }
    }
    else
    {
      std::ostringstream buf;
      buf << i;
      filename += (std::string(".") + std::string(buf.str().c_str()));
    }
    filename += extension;
    if (!gzExtension.empty())
    {
      filename += std::string(".gz");
    }

    m_ComponentImageFileName = filename;

    ImageType * input = selector->GetOutput();
    this->m_Image = selector->GetOutput();

    itkDebugMacro(<< "Writing an image file");

    // Make sure input is available
    if (input == 0)
    {
      itkExceptionMacro(<< "No input to writer!");
    }

    // Make sure that we can write the file given the name
    //
    if (filename == "")
    {
      itkExceptionMacro(<< "No filename was specified");
    }

    if (m_ImageIO.IsNull()) // try creating via factory
    {
      itkDebugMacro(<< "Attempting factory creation of ImageIO for file: " << filename);
      m_ImageIO = ImageIOFactory::CreateImageIO(filename.c_str(), ImageIOFactory::WriteMode);
      m_FactorySpecifiedImageIO = true;
    }
    else
    {
      if (m_FactorySpecifiedImageIO && !m_ImageIO->CanWriteFile(filename.c_str()))
      {
        itkDebugMacro(<< "ImageIO exists but doesn't know how to write file:" << m_FileName);
        itkDebugMacro(<< "Attempting creation of ImageIO with a factory for file:" << m_FileName);
        m_ImageIO = ImageIOFactory::CreateImageIO(filename.c_str(), ImageIOFactory::WriteMode);
        m_FactorySpecifiedImageIO = true;
      }
    }

    if (m_ImageIO.IsNull())
    {
      ImageFileWriterException e(__FILE__, __LINE__);
      std::ostringstream       msg;
      msg << " Could not create IO object for file " << filename.c_str() << std::endl;
      msg << "  Tried to create one of the following:" << std::endl;
      std::list<LightObject::Pointer> allobjects = ObjectFactoryBase::CreateAllInstance("itkImageIOBase");
      for (std::list<LightObject::Pointer>::iterator i = allobjects.begin(); i != allobjects.end(); ++i)
      {
        ImageIOBase * io = dynamic_cast<ImageIOBase *>(i->GetPointer());
        msg << "    " << io->GetNameOfClass() << std::endl;
      }
      msg << "  You probably failed to set a file suffix, or" << std::endl;
      msg << "    set the suffix to an unsupported type." << std::endl;
      e.SetDescription(msg.str().c_str());
      e.SetLocation(ITK_LOCATION);
      throw e;
    }

    // NOTE: this const_cast<> is due to the lack of const-correctness
    // of the ProcessObject.
    ImageType * nonConstImage = input;

    typedef typename TImage::RegionType RegionType;

    if (!m_UserSpecifiedIORegion)
    {
      // Make sure the data is up-to-date.
      if (nonConstImage->GetSource())
      {
        nonConstImage->GetSource()->UpdateLargestPossibleRegion();
      }
      // Write the whole image
      ImageIORegion ioRegion(TImage::ImageDimension);
      RegionType    region = this->m_Image->GetLargestPossibleRegion();
      for (unsigned int i = 0; i < TVectorImage::ImageDimension; i++)
      {
        ioRegion.SetSize(i, region.GetSize(i));
        ioRegion.SetIndex(i, region.GetIndex(i));
      }
      m_IORegion = ioRegion; // used by GenerateData
    }
    else
    {
      nonConstImage->Update();
    }

    // Setup the ImageIO
    //
    m_ImageIO->SetNumberOfDimensions(TImage::ImageDimension);
    RegionType                             region = this->m_Image->GetLargestPossibleRegion();
    const typename TImage::SpacingType &   spacing = this->m_Image->GetSpacing();
    const typename TImage::PointType &     origin = this->m_Image->GetOrigin();
    const typename TImage::DirectionType & direction = this->m_Image->GetDirection();
    for (unsigned int k = 0; k < TVectorImage::ImageDimension; k++)
    {
      m_ImageIO->SetDimensions(k, region.GetSize(k));
      m_ImageIO->SetSpacing(k, spacing[k]);
      m_ImageIO->SetOrigin(k, origin[k]);
      vnl_vector<double> axisDirection(TVectorImage::ImageDimension);
      // Please note: direction cosines are stored as columns of the
      // direction matrix
      for (unsigned int j = 0; j < TImage::ImageDimension; j++)
      {
        axisDirection[j] = direction[j][k];
      }
      m_ImageIO->SetDirection(k, axisDirection);
    }

    m_ImageIO->SetUseCompression(m_UseCompression);
    m_ImageIO->SetIORegion(m_IORegion);
    if (m_UseInputMetaDataDictionary)
    {
      m_ImageIO->SetMetaDataDictionary(input->GetMetaDataDictionary());
    }
    // Notify start event observers
    this->InvokeEvent(StartEvent());

    // Actually do something
    this->GenerateData();

    // Notify end event observers
    this->InvokeEvent(EndEvent());

    // Release upstream data if requested
    if (input->ShouldIReleaseData())
    {
      nonConstImage->ReleaseData();
    }
  }
}

// ---------------------------------------------------------
template <typename TVectorImage, typename TImage>
void
VectorImageFileWriter<TVectorImage, TImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (m_FileName.data() ? m_FileName.data() : "(none)") << std::endl;

  os << "Number of vector components (i.e. number of images created): " << VectorImageType::PixelType::VectorDimension
     << std::endl;

  os << indent << "Image IO: ";
  if (m_ImageIO.IsNull())
  {
    os << "(none)\n";
  }
  else
  {
    os << m_ImageIO << "\n";
  }

  os << indent << "IO Region: " << m_IORegion << "\n";

  if (m_UseCompression)
  {
    os << indent << "Compression: On\n";
  }
  else
  {
    os << indent << "Compression: Off\n";
  }

  if (m_UseInputMetaDataDictionary)
  {
    os << indent << "UseInputMetaDataDictionary: On\n";
  }
  else
  {
    os << indent << "UseInputMetaDataDictionary: Off\n";
  }

  if (m_FactorySpecifiedImageIO)
  {
    os << indent << "FactorySpecifiedmageIO: On\n";
  }
  else
  {
    os << indent << "FactorySpecifiedmageIO: Off\n";
  }
}
} // end namespace itk

#endif
