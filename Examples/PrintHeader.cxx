/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: PrintHeader.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/05 20:09:47 $
  Version:   $Revision: 1.4 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <sys/stat.h>

#include <fstream>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkSpatialOrientation.h"

using namespace std;
/** below code from Paul Yushkevich's c3d */
template <class AnyType>
bool
try_print_metadata(itk::MetaDataDictionary & mdd, std::string key)
{
  AnyType value = 0;

  if( itk::ExposeMetaData<AnyType>(mdd, key, value) )
    {
    cout << "    " << key << " = " << value << endl;
    return true;
    }
  else
    {
    return false;
    }
}

string
get_rai_code(itk::SpatialOrientation::ValidCoordinateOrientationFlags code)
{
  std::map<itk::SpatialOrientation::ValidCoordinateOrientationFlags, string> m_CodeToString;

  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP] = "RIP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP] = "LIP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP] = "RSP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP] = "LSP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA] = "RIA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA] = "LIA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA] = "RSA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA] = "LSA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP] = "IRP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP] = "ILP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP] = "SRP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP] = "SLP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA] = "IRA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA] = "ILA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA] = "SRA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA] = "SLA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI] = "RPI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI] = "LPI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI] = "RAI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI] = "LAI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS] = "RPS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS] = "LPS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS] = "RAS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS] = "LAS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI] = "PRI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI] = "PLI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI] = "ARI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI] = "ALI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS] = "PRS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS] = "PLS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS] = "ARS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS] = "ALS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR] = "IPR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR] = "SPR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR] = "IAR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR] = "SAR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL] = "IPL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL] = "SPL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL] = "IAL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL] = "SAL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR] = "PIR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR] = "PSR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR] = "AIR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR] = "ASR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL] = "PIL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL] = "PSL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL] = "AIL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL] = "ASL";
  return m_CodeToString[code];
}

template <unsigned int ImageDimension>
int PrintHeader(int argc, char *argv[])
{
  typedef  float                                     outPixelType;
  typedef  float                                     floatPixelType;
  typedef  float                                     inPixelType;
  typedef itk::Image<inPixelType, ImageDimension>    ImageType;
  typedef itk::Image<floatPixelType, ImageDimension> IntermediateType;
  typedef itk::Image<outPixelType, ImageDimension>   OutImageType;
  typedef itk::ImageFileReader<ImageType>            readertype;
  typedef itk::ImageFileWriter<OutImageType>         writertype;

  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  std::cout << " Spacing " << reader->GetOutput()->GetSpacing() << std::endl;
  std::cout << " Origin " << reader->GetOutput()->GetOrigin() << std::endl;
  std::cout << " Direction " << std::endl << reader->GetOutput()->GetDirection() << std::endl;
  if( ImageDimension == 1 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " " <<   std::endl;
    }
  else if( ImageDimension == 2 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << " " << std::endl;
    }
  else if( ImageDimension == 3 )
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << " " <<  " "
              << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    }
  else
    {
    std::cout << " Size : " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
    }
//  std::cout << " Orientation " << reader->GetOutput()->GetOrientation() << std::endl;

  unsigned int VDim = ImageDimension;
  // Get the input image
  typename ImageType::Pointer image = reader->GetOutput();

  // Compute the bounding box
  vnl_vector<double> bb0, bb1, ospm;
  bb0.set_size(VDim);
  bb1.set_size(VDim);
  ospm.set_size(VDim);
  for( size_t i = 0; i < VDim; i++ )
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    ospm[i] = -image->GetOrigin()[i] / image->GetSpacing()[i];
    }

  // Compute the intensity range of the image
  size_t n = image->GetBufferedRegion().GetNumberOfPixels();
  float *vox = image->GetBufferPointer();
  double iMax = vox[0], iMin = vox[0], iMean = vox[0];
  for( size_t i = 1; i < n; i++ )
    {
    iMax = (iMax > vox[i]) ? iMax : vox[i];
    iMin = (iMin < vox[i]) ? iMin : vox[i];
    iMean += vox[i];
    }
  iMean /= n;

  // Short or long?
  bool full = true;
  if( !full )
    {
    cout << " dim = " << image->GetBufferedRegion().GetSize() << "; ";
    cout << " bb = {[" << bb0 << "], [" << bb1 << "]}; ";
    cout << " vox = " << image->GetSpacing() << "; ";
    cout << " range = [" << iMin << ", " << iMax << "]; ";
    cout << endl;
    }
  else
    {
    cout << endl;
    cout << "  Image Dimensions   : " << image->GetBufferedRegion().GetSize() << endl;
    cout << "  Bounding Box       : " << "{[" << bb0 << "], [" << bb1 << "]}" << endl;
    cout << "  Voxel Spacing      : " << image->GetSpacing() << endl;
    cout << "  Intensity Range    : [" << iMin << ", " << iMax << "]" << endl;
    cout << "  Mean Intensity     : " << iMean << endl;
    cout << "  Direction Cos Mtx. : " << endl;
    std::cout << image->GetDirection().GetVnlMatrix() << std::endl;
    // Print NIFTI s-form matrix (check against freesurfer's MRIinfo)
    cout << "  Voxel->RAS x-form  : " << endl;
    //    image->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
    //    std::cout << image->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix() << std::endl;

    //
    // Print metadata
    cout << "  Image Metadata: " << endl;
    itk::MetaDataDictionary &              mdd = image->GetMetaDataDictionary();
    itk::MetaDataDictionary::ConstIterator itMeta;
    for( itMeta = mdd.Begin(); itMeta != mdd.End(); ++itMeta )
      {
      // Get the metadata as a generic object
      string                                                   key = itMeta->first, v_string;
      itk::SpatialOrientation::ValidCoordinateOrientationFlags v_oflags;

      if( itk::ExposeMetaData<string>(mdd, key, v_string) )
        {
        // For some weird reason, some of the strings returned by this method
        // contain '\0' characters. We will replace them by spaces
        std::ostringstream sout("");
        for( unsigned int i = 0; i < v_string.length(); i++ )
          {
          if( v_string[i] >= ' ' )
            {
            sout << v_string[i];
            }
          }
        v_string = sout.str();

        // Make sure the value has more than blanks
        if( v_string.find_first_not_of(" ") != v_string.npos )
          {
          cout << "    " << key << " = " << v_string << endl;
          }
        }
      else if( itk::ExposeMetaData(mdd, key, v_oflags) )
        {
        cout << "    " << key << " = " << get_rai_code(v_oflags) << endl;
        }
      else
        {
        bool rc = false;
        if( !rc )
          {
          rc |= try_print_metadata<double>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<float>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<int>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<unsigned int>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<long>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<unsigned long>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<short>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<unsigned short>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<char>(mdd, key);
          }
        if( !rc )
          {
          rc |= try_print_metadata<unsigned char>(mdd, key);
          }

        if( !rc )
          {
          cout << "    " << key << " of unsupported type "
               << itMeta->second->GetMetaDataObjectTypeName() << endl;
          }
        }
      }
    }

  return 1;
}

bool FileExists(string strFilename)
{
  struct stat stFileInfo;
  bool        blnReturn;
  int         intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(), &stFileInfo);
  if( intStat == 0 )
    {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
    }
  else
    {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
    }

  return blnReturn;
}

int main(int argc, char *argv[])
{
  if( argc < 2  || ( (argc == 2) && strcmp(argv[1], "--help") == 0) )
    {
    std::cout << "Usage:  " << argv[0] << " image.ext " << std::endl;
    return 1;
    }
  // Get the image dimension
  std::string fn = std::string(argv[1]);
  if( !FileExists(fn) )
    {
    std::cout << " file " << fn << " does not exist . " << std::endl;  return 1;
    }
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(
      fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  try
    {
    imageIO->ReadImageInformation();
    }
  catch( ... )
    {
    std::cout << " cant read " << fn << std::endl;
    return 1;
    }

  switch( imageIO->GetNumberOfDimensions() )
    {
    case 2:
      {
      PrintHeader<2>(argc, argv);
      }
      break;
    case 3:
      {
      PrintHeader<3>(argc, argv);
      }
      break;
    case 4:
      {
      PrintHeader<4>(argc, argv);
      }
      break;
    default:
      std::cerr << "Unsupported dimension " <<  imageIO->GetNumberOfDimensions() << std::endl;
      exit( EXIT_FAILURE );
    }

  return 0;
}
