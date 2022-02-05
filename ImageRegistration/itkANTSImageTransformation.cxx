/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkANTSImageTransformation_hxx_
#define _itkANTSImageTransformation_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#  pragma warning(disable : 4786)
#endif

#include "ANTS_affine_registration2.h"
#include "itkANTSImageTransformation.h"

#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkDisplacementFieldTransform.h"

#include "itkMath.h"

namespace itk
{
template <unsigned int TDimension, typename TReal>
ANTSImageTransformation<TDimension, TReal>::ANTSImageTransformation()
{
  this->m_DisplacementField = nullptr;
  this->m_InverseDisplacementField = nullptr;
  this->m_AffineTransform = nullptr;
  this->m_WriteComponentImages = false;
  m_DeformationRegionOfInterestSize.Fill(0);
  m_DeformationRegionSpacing.Fill(1);
  m_DeformationRegionOfInterestCenter.Fill(0);
}

template <unsigned int TDimension, typename TReal>
void
ANTSImageTransformation<TDimension, TReal>::Compose()
{}

template <unsigned int TDimension, typename TReal>
void
ANTSImageTransformation<TDimension, TReal>::Write()
{
  std::cout << " begin writing " << m_NamingConvention << std::endl;

  std::string            filePrefix = this->m_NamingConvention;
  std::string::size_type pos = filePrefix.rfind(".");
  std::string            extension;
  std::string            gzExtension("");

  if (pos != std::string::npos && std::string(filePrefix, pos, pos + 2) != std::string("./"))
  {
    bool is_supported = false;
    filePrefix = std::string(filePrefix, 0, pos);
    extension = std::string(this->m_NamingConvention, pos, this->m_NamingConvention.length() - 1);
    if (extension != std::string(".xfm"))
    {

      if (extension == std::string(".gz"))
      {
        gzExtension = std::string(".gz");
        pos = filePrefix.rfind(".");
        extension = std::string(filePrefix, pos, filePrefix.length() - 1);
        filePrefix = std::string(filePrefix, 0, pos);
      }
      // GetSupportedWriteExtensions
      typedef itk::ImageIOBase                           IOBaseType;
      typedef std::list<itk::LightObject::Pointer>       ArrayOfImageIOType;
      typedef typename IOBaseType::ArrayOfExtensionsType ArrayOfExtensionsType;
      ArrayOfImageIOType           allobjects = itk::ObjectFactoryBase::CreateAllInstance("itkImageIOBase");
      ArrayOfImageIOType::iterator itr = allobjects.begin();
      while (itr != allobjects.end())
      {
        IOBaseType * io = dynamic_cast<IOBaseType *>(itr->GetPointer());
        if (!io)
        {
          std::cout << "Got a null pointer in the array" << std::endl;
        }
        else
        {
          const ArrayOfExtensionsType &         writeExtensions = io->GetSupportedWriteExtensions();
          ArrayOfExtensionsType::const_iterator writeItr = writeExtensions.begin();
          while (writeItr != writeExtensions.end())
          {
            std::string test_ext = *writeItr;
            //  std::cout <<" compare " << extension << " to " << test_ext << std::endl;
            if (extension == test_ext)
            {
              is_supported = true;
            }
            //        else std::cout <<" not the same " << std::endl;
            ++writeItr;
          }
        }
        ++itr;
      }

      if (!is_supported)
      {
        std::cout << " WARNING! we are guessing you want .nii.gz for your image output instead of: " << extension
                  << std::endl;
        filePrefix = this->m_NamingConvention;
        extension = std::string(".nii.gz");
      }
    }
  }
  else
  {
    extension = std::string(".nii.gz");
  }

  if (extension == std::string(".xfm"))
  {
    typedef itk::DisplacementFieldTransform<TReal, TDimension> DisplacementFieldTransform;

    itk::TransformFactory<itk::DisplacementFieldTransform<TReal, TDimension>>::RegisterTransform();
    itk::TransformFactory<itk::MatrixOffsetTransformBase<TReal, TDimension, TDimension>>::RegisterTransform();

    std::string fw_filename = filePrefix + extension;
    std::string bw_filename = filePrefix + "_inverse" + extension;

    typedef itk::TransformFileWriterTemplate<TReal> TransformFileWriterFloat;

    typename TransformFileWriterFloat::Pointer fw_writer = TransformFileWriterFloat::New();

    if (this->m_AffineTransform)
      fw_writer->AddTransform(this->m_AffineTransform);
    if (this->m_DisplacementField)
    {
      typename DisplacementFieldTransform::Pointer disp = DisplacementFieldTransform::New();
      disp->SetDisplacementField(this->m_DisplacementField);
      fw_writer->AddTransform(disp);
    }

    fw_writer->SetFileName(fw_filename);
    fw_writer->Update();

    typename TransformFileWriterFloat::Pointer bw_writer = TransformFileWriterFloat::New();

    if (this->m_InverseDisplacementField)
    {
      typename DisplacementFieldTransform::Pointer disp = DisplacementFieldTransform::New();
      disp->SetDisplacementField(this->m_InverseDisplacementField);

      bw_writer->AddTransform(disp);
    }

    typename AffineTransformType::Pointer tmp = AffineTransformType::New();
    if (this->m_AffineTransform)
    {
      this->m_AffineTransform->GetInverse(tmp);
      bw_writer->AddTransform(tmp);
    }

    bw_writer->SetFileName(bw_filename);
    bw_writer->Update();
  }
  else
  {
    // Added by songgang
    if (this->m_AffineTransform)
    {
      std::cout << " writing " << filePrefix << " affine " << std::endl;
      const std::string filename = filePrefix + std::string("Affine.txt");
      WriteAffineTransformFile<AffineTransformType>(this->m_AffineTransform, filename);
    }

    if (this->m_DisplacementField)
    {
      std::cout << " writing " << filePrefix << " def " << std::endl;
      if (extension != std::string(".mha"))
      {
        std::string filename = filePrefix + std::string("Warp") + extension + gzExtension;

        typedef ImageFileWriter<DisplacementFieldType> WriterType;
        typename WriterType::Pointer                   writer = WriterType::New();
        std::cout << "filename " << filename << std::endl;
        writer->SetFileName(filename);
        //            writer->SetUseAvantsNamingConvention( true );
        writer->SetInput(this->m_DisplacementField);
        writer->Update();
      }
      else
      {
        std::string                                    filename = filePrefix + std::string("Warp.nii.gz");
        typedef ImageFileWriter<DisplacementFieldType> WriterType;
        typename WriterType::Pointer                   writer = WriterType::New();
        writer->SetFileName(filename);
        writer->SetInput(this->m_DisplacementField);
        writer->Update();
      }
    }

    if (this->m_InverseDisplacementField)
    {
      if (extension != std::string(".mha"))
      {
        std::string filename = filePrefix + std::string("InverseWarp") + extension + gzExtension;
        typedef ImageFileWriter<DisplacementFieldType> WriterType;
        typename WriterType::Pointer                   writer = WriterType::New();
        writer->SetFileName(filename);
        //            writer->SetUseAvantsNamingConvention( true );
        writer->SetInput(this->m_InverseDisplacementField);
        writer->Update();
      }
      else
      {
        std::string                                    filename = filePrefix + std::string("InverseWarp.mha");
        typedef ImageFileWriter<DisplacementFieldType> WriterType;
        typename WriterType::Pointer                   writer = WriterType::New();
        writer->SetFileName(filename);
        writer->SetInput(this->m_InverseDisplacementField);
        writer->Update();
      }
    }
  }
}

/**
 * Standard "PrintSelf" method
 */
template <unsigned int TDimension, typename TReal>
void
ANTSImageTransformation<TDimension, TReal>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk
#endif
