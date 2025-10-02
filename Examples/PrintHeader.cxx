/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include <algorithm>

#include <iomanip>
#include <iostream>
#include <sys/stat.h>

#include <fstream>
#include <cstdio>
#include <cstring>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkSpatialOrientation.h"
#include "itkImageIOFactory.h"
#include "itkImageIOBase.h"

#include <map>
#include <sstream>
#include <type_traits>
#include <cmath>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_diag_matrix.h>

namespace ants
{
using namespace std;

// Aliases for metadata composite types commonly encountered
// Matrices
using Mat22d = itk::Matrix<double, 2, 2>;
using Mat33d = itk::Matrix<double, 3, 3>;
using Mat44d = itk::Matrix<double, 4, 4>;
using Mat55d = itk::Matrix<double, 5, 5>;
using Mat22f = itk::Matrix<float,  2, 2>;
using Mat33f = itk::Matrix<float,  3, 3>;
using Mat44f = itk::Matrix<float,  4, 4>;  // qto_xyz, sto_xyz, qto_ijk, sto_ijk

// Fixed arrays (spacing/origin sometimes appear this way)
using Fix1d  = itk::FixedArray<double, 1>;
using Fix2d  = itk::FixedArray<double, 2>;
using Fix3d  = itk::FixedArray<double, 3>;
using Fix4d  = itk::FixedArray<double, 4>;
using Fix5d  = itk::FixedArray<double, 5>;

/** below code (metadata printing) inspired by Paul Yushkevich's c3d */

// ---------- Helpers: metadata printers ----------

// Generic helper to pretty-print a vnl_matrix
template <typename TVnlMatrix>
void
print_vnl_matrix(const TVnlMatrix & mat,
                 int precision = 5,
                 int width = 11,
                 int indent = 8)
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(precision);

  std::string spaces(indent, ' ');

  for (unsigned int i = 0; i < mat.rows(); ++i)
  {
    std::cout << spaces;
    for (unsigned int j = 0; j < mat.cols(); ++j)
    {
      std::cout << std::setw(width) << mat(i, j);
    }
    std::cout << std::endl;
  }
}

template <typename AnyType>
bool
try_print_metadata(const itk::MetaDataDictionary & mdd, const std::string & key)
{
  AnyType value{};
  if (itk::ExposeMetaData<AnyType>(mdd, key, value))
  {
    std::cout << "    " << key << " = " << value << std::endl;
    return true;
  }
  return false;
}

template <typename T>
bool
try_print_metadata_std_vector(const itk::MetaDataDictionary & mdd, const std::string & key)
{
  std::vector<T> v;
  if (itk::ExposeMetaData<std::vector<T>>(mdd, key, v))
  {
    std::cout << "    " << key << " = [";
    for (size_t i = 0; i < v.size(); ++i)
    {
      std::cout << v[i] << (i + 1 < v.size() ? ", " : "");
    }
    std::cout << "]" << std::endl;
    return true;
  }
  return false;
}

template <typename TFixedArray>
bool try_print_metadata_fixed_array(const itk::MetaDataDictionary& mdd, const std::string& key)
{
  TFixedArray arr;
  if (itk::ExposeMetaData<TFixedArray>(mdd, key, arr))
  {
    std::cout << "    " << key << " = [";
    for (unsigned int i = 0; i < arr.Size(); ++i)
    {
      std::cout << arr[i];
      if (i + 1 < arr.Size()) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    return true;
  }
  return false;
}


template <typename TMatrix>
bool
try_print_metadata_matrix(const itk::MetaDataDictionary & mdd,
                          const std::string & key,
                          int precision = 5,
                          int width = 11,
                          int indent = 8)
{
  TMatrix mat;
  if (itk::ExposeMetaData<TMatrix>(mdd, key, mat))
  {
    std::cout << "    " << key << " =" << std::endl;
    auto vnl_mat = mat.GetVnlMatrix();
    print_vnl_matrix(vnl_mat, precision, width, indent);
    return true;
  }
  return false;
}

// ---------- Math helpers ----------

// Copy top-left MxN block from a dynamic vnl_matrix into a fixed one
template <unsigned int M, unsigned int N>
vnl_matrix_fixed<double, M, N>
TopLeftBlock(const vnl_matrix<double> & A)
{
  vnl_matrix_fixed<double, M, N> B;
  B.fill(0.0);
  const unsigned int mr = std::min<unsigned int>(M, A.rows());
  const unsigned int nc = std::min<unsigned int>(N, A.cols());
  for (unsigned int i = 0; i < mr; ++i)
  {
    for (unsigned int j = 0; j < nc; ++j)
    {
      B(i, j) = A(i, j);
    }
  }
  return B;
}

// Dimension-agnostic Voxel(index)->RAS(world) 4×4 affine from ITK LPS geometry
template <typename TImage>
vnl_matrix_fixed<double, 4, 4>
GetVoxelSpaceToRASPhysicalSpaceMatrix(const TImage * image)
{
  static_assert(TImage::ImageDimension >= 1, "ImageDimension must be >= 1");
  constexpr unsigned int Dim = TImage::ImageDimension;
  const unsigned int S = Dim >= 3 ? 3u : Dim; // spatial axes considered

  // Direction embedded to 3×3
  const vnl_matrix<double> Dfull = image->GetDirection().GetVnlMatrix().as_matrix();
  vnl_matrix_fixed<double, 3, 3> D3;
  D3.set_identity();
  if (S == 3)
  {
    D3 = TopLeftBlock<3, 3>(Dfull);
  }
  else if (S == 2)
  {
    D3 = TopLeftBlock<3, 3>(Dfull);
    D3(2, 2) = 1.0;
  }
  else // S == 1
  {
    D3.set_identity();
  }

  // Spacing padded to 3
  vnl_vector_fixed<double, 3> s3(1.0, 1.0, 1.0);
  for (unsigned int i = 0; i < S; ++i)
  {
    s3[i] = image->GetSpacing()[i];
  }

  // Origin padded to 3
  vnl_vector_fixed<double, 3> o3(0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < S; ++i)
  {
    o3[i] = image->GetOrigin()[i];
  }

  // LPS->RAS = diag(-1, -1, +1)
  vnl_diag_matrix<double> LPS2RAS(3);
  LPS2RAS.fill(1.0);
  LPS2RAS[0] = -1.0;
  LPS2RAS[1] = -1.0;

  // M = LPS->RAS * (D3 * diag(s))
  vnl_diag_matrix<double> Sdiag(s3.as_ref());
  vnl_matrix<double> M = LPS2RAS * (D3.as_ref() * Sdiag.as_ref());

  // t = LPS->RAS * origin
  vnl_vector<double> t = LPS2RAS * o3.as_ref();

  // Assemble homogeneous 4×4
  vnl_matrix_fixed<double, 4, 4> H;
  H.set_identity();
  H.update(M);
  H(0, 3) = t[0];
  H(1, 3) = t[1];
  H(2, 3) = t[2];
  return H;
}

// RAI code from direction matrix (supports 1D/2D/3D)
template <unsigned int VDim>
std::string
GetRAICodeFromDirectionMatrix(const vnl_matrix_fixed<double, VDim, VDim> & dir)
{
  static_assert(VDim >= 1 && VDim <= 3, "VDim must be 1..3");
  const char codes[3][2] = { { 'R', 'L' }, { 'A', 'P' }, { 'I', 'S' } };

  char out[VDim + 1];
  out[VDim] = '\0';
  bool oblique = false;

  for (unsigned int i = 0; i < VDim; ++i)
  {
    double       amax = 0.0;
    unsigned int argmax = 0;
    unsigned int sgn = 0;
    for (unsigned int j = 0; j < VDim; ++j)
    {
      const double a = std::abs(dir(j, i));
      if (a > amax)
      {
        amax = a;
        argmax = j;
        sgn = (dir(j, i) >= 0.0) ? 0 : 1;
      }
    }
    if (std::abs(amax - 1.0) > 1e-6)
    {
      oblique = true;
    }
    out[i] = codes[argmax][sgn];
  }

  if (oblique)
  {
    std::ostringstream sout;
    sout << "Oblique, closest to " << out;
    return sout.str();
  }
  return std::string(out);
}



// ---------- Core templated PrintHeader ----------

template <unsigned int ImageDimension>
int
PrintHeader(int argc, char * argv[])
{
  using inPixelType = float;
  using ImageType = itk::Image<inPixelType, ImageDimension>;
  using readertype = itk::ImageFileReader<ImageType>;

  typename readertype::Pointer reader = readertype::New();
  if (argc < 2)
  {
    std::cout << "missing input image name" << std::endl;
    throw;
  }
  reader->SetFileName(argv[1]);
  reader->Update();

  // Print only specific header information
  if (argc > 2)
  {
    switch (std::stoi(argv[2]))
    {
      case 0: // origin
      {
        for (int d = 0; d < static_cast<int>(ImageDimension) - 1; d++)
        {
          std::cout << reader->GetOutput()->GetOrigin()[d] << 'x';
        }
        std::cout << reader->GetOutput()->GetOrigin()[static_cast<int>(ImageDimension) - 1] << std::endl;
        break;
      }
      case 1: // spacing
      {
        for (int d = 0; d < static_cast<int>(ImageDimension) - 1; d++)
        {
          std::cout << reader->GetOutput()->GetSpacing()[d] << 'x';
        }
        std::cout << reader->GetOutput()->GetSpacing()[static_cast<int>(ImageDimension) - 1] << std::endl;
        break;
      }
      case 2: // size
      {
        for (int d = 0; d < static_cast<int>(ImageDimension) - 1; d++)
        {
          std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] << 'x';
        }
        std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[static_cast<int>(ImageDimension) - 1]
                  << std::endl;
        break;
      }
      case 3: // index
      {
        for (int d = 0; d < static_cast<int>(ImageDimension) - 1; d++)
        {
          std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[d] << 'x';
        }
        std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[static_cast<int>(ImageDimension) - 1]
                  << std::endl;
        break;
      }
      case 4: // direction
      {
        for (int di = 0; di < static_cast<int>(ImageDimension); di++)
        {
          for (int dj = 0; dj < static_cast<int>(ImageDimension); dj++)
          {
            std::cout << reader->GetOutput()->GetDirection()[di][dj];
            if (di == dj && di == static_cast<int>(ImageDimension) - 1)
            {
              std::cout << std::endl;
            }
            else
            {
              std::cout << 'x';
            }
          }
        }
        break;
      }
      default:
        break;
    }
    return EXIT_SUCCESS;
  }

  // else print out entire header information

  typename ImageType::Pointer image = reader->GetOutput();
  unsigned int VDim = ImageDimension;

  // Compute the bounding box and origin/spacing per mm
  vnl_vector<double> bb0, bb1, ospm;
  bb0.set_size(VDim);
  bb1.set_size(VDim);
  ospm.set_size(VDim);
  for (size_t i = 0; i < VDim; i++)
  {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    ospm[i] = -image->GetOrigin()[i] / image->GetSpacing()[i];
  }

  // Compute the intensity range of the image (assumes float pixel type here)
  size_t  n = image->GetBufferedRegion().GetNumberOfPixels();
  float * vox = image->GetBufferPointer();
  double  iMax = vox[0], iMin = vox[0], iMean = vox[0];
  for (size_t i = 1; i < n; i++)
  {
    iMax = (iMax > static_cast<double>(vox[i])) ? iMax : static_cast<double>(vox[i]);
    iMin = (iMin < static_cast<double>(vox[i])) ? iMin : static_cast<double>(vox[i]);
    iMean += static_cast<double>(vox[i]);
  }
  iMean /= static_cast<double>(n);

  // Full header dump
  bool full = true;
  if (!full)
  {
    cout << " dims = " << image->GetBufferedRegion().GetSize() << "; ";
    cout << " bb = {[" << bb0 << "], [" << bb1 << "]}; ";
    cout << " spc = " << image->GetSpacing() << "; ";
    cout << " orig = " << image->GetOrigin() << "; ";
    cout << " range = [" << iMin << ", " << iMax << "]; ";
    cout << endl;
  }
  else
  {
    cout << endl;
    cout << "  Image Dimensions   : " << image->GetBufferedRegion().GetSize() << endl;
    cout << "  Bounding Box       : "
         << "{[" << bb0 << "], [" << bb1 << "]}" << endl;
    cout << "  Origin             : " << image->GetOrigin() << endl;
    cout << "  Voxel Spacing      : " << image->GetSpacing() << endl;
    cout << "  Intensity Range    : [" << iMin << ", " << iMax << "]" << endl;
    cout << "  Mean Intensity     : " << iMean << endl;

    // Canonical orientation (RAI-like string) from direction
    {
      constexpr unsigned int Dim = ImageDimension;
      const unsigned int S = Dim >= 3 ? 3u : Dim;
      if (S == 3)
      {
        auto D3 = TopLeftBlock<3, 3>(image->GetDirection().GetVnlMatrix().as_matrix());
        cout << "  Canon. FROM Orient : " << GetRAICodeFromDirectionMatrix<3>(D3) << endl;
      }
      else if (S == 2)
      {
        auto D2 = TopLeftBlock<2, 2>(image->GetDirection().GetVnlMatrix().as_matrix());
        cout << "  Canon. FROM Orient : " << GetRAICodeFromDirectionMatrix<2>(D2) << endl;
      }
      else /* S == 1 */
      {
        vnl_matrix_fixed<double, 1, 1> D1;
        D1(0, 0) = 1.0;
        cout << "  Canon. FROM Orient : " << GetRAICodeFromDirectionMatrix<1>(D1) << endl;
      }
    }

    cout << "  Direction Cos Mtx. : " << endl;
    print_vnl_matrix(image->GetDirection().GetVnlMatrix());

    // Print Voxel->RAS 4x4 x-form (derived from ITK LPS geometry)
    cout << "  Voxel->RAS+ x-form  : " << endl;
    print_vnl_matrix(GetVoxelSpaceToRASPhysicalSpaceMatrix(image.GetPointer()));

    //
    // Print metadata
    cout << "  Image Metadata: " << endl;
    const itk::MetaDataDictionary &              mdd = image->GetMetaDataDictionary();
    itk::MetaDataDictionary::ConstIterator itMeta;
    for (itMeta = mdd.Begin(); itMeta != mdd.End(); ++itMeta)
    {
      // Get the metadata as a generic object
      std::string                                               key = itMeta->first, v_string;
      if (itk::ExposeMetaData<std::string>(mdd, key, v_string))
      {
        // For some weird reason, some of the strings returned by this method
        // contain '\0' characters. We will replace them by spaces
        std::ostringstream sout("");
        for (char c : v_string)
        {
          if (c >= ' ')
          {
            sout << c;
          }
        }
        v_string = sout.str();

        // Make sure the value has more than blanks
        if (v_string.find_first_not_of(" ") != v_string.npos)
        {
          cout << "    " << key << " = " << v_string << endl;
        }
      }
      else
      {
        bool rc = false;
        // order here depends on actual type in the ITK metadata dictionary, rather than
        // a successful attempt to parse the data to a specific type
        if (!rc) { rc |= try_print_metadata<double>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<float>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<int>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<unsigned int>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<long>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<unsigned long>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<short>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<unsigned short>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<char>(mdd, key); }
        if (!rc) { rc |= try_print_metadata<unsigned char>(mdd, key); }
        // matrix metadata (handle spatial-only and full-dim cases)
        if (!rc) { rc |= try_print_metadata_matrix<Mat22d>(mdd, key); }  // 2D original direction
        if (!rc) { rc |= try_print_metadata_matrix<Mat33d>(mdd, key); }  // 3D original direction
        if (!rc) { rc |= try_print_metadata_matrix<Mat44d>(mdd, key); }  // 4D original direction
        if (!rc) { rc |= try_print_metadata_matrix<Mat22f>(mdd, key); }  // float variants
        if (!rc) { rc |= try_print_metadata_matrix<Mat33f>(mdd, key); }
        if (!rc) { rc |= try_print_metadata_matrix<Mat44f>(mdd, key); }  // qto_xyz / sto_xyz (NIfTI)
        // fixed-array variants of spacing/origin/etc. (varies by writer/reader)
        if (!rc) { rc |= try_print_metadata_fixed_array<Fix1d>(mdd, key); }
        if (!rc) { rc |= try_print_metadata_fixed_array<Fix2d>(mdd, key); }
        if (!rc) { rc |= try_print_metadata_fixed_array<Fix3d>(mdd, key); }
        if (!rc) { rc |= try_print_metadata_fixed_array<Fix4d>(mdd, key); }
        if (!rc) { rc |= try_print_metadata_fixed_array<Fix5d>(mdd, key); }

        // dynamic vector fallback (e.g., ITK_original_spacing)
        if (!rc) { rc |= try_print_metadata_std_vector<double>(mdd, key); }


        if (!rc)
        {
          cout << "    " << key << " of unsupported type " << itMeta->second->GetMetaDataObjectTypeName() << endl;
        }
      }
    }
  }
  return EXIT_SUCCESS;
}

// ---------- Utility: FileExists ----------

bool
FileExists(std::string strFilename)
{
  struct stat stFileInfo;
  int         intStat = stat(strFilename.c_str(), &stFileInfo);
  return (intStat == 0);
}

// ---------- Entry point wrapper used by ANTs tools ----------

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
PrintHeader(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "PrintHeader");

  int     argc = static_cast<int>(args.size());
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;

  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}

    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  if (argc < 2 || ((argc == 2) && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)))
  {
    std::cout << "Usage:  " << argv[0] << " image.ext [whatInformation]" << std::endl;
    std::cout << "  whatInformation:  " << std::endl;
    std::cout << "    0 = origin" << std::endl;
    std::cout << "    1 = spacing" << std::endl;
    std::cout << "    2 = size" << std::endl;
    std::cout << "    3 = index" << std::endl;
    std::cout << "    4 = direction" << std::endl;
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  // Get the image dimension
  std::string fn = std::string(argv[1]);
  if (!FileExists(fn))
  {
    std::cout << " file " << fn << " does not exist . " << std::endl;
    return EXIT_FAILURE;
  }
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::IOFileModeEnum::ReadMode);
  imageIO->SetFileName(fn.c_str());
  try
  {
    imageIO->ReadImageInformation();
  }
  catch (...)
  {
    std::cout << " cant read " << fn << std::endl;
    return EXIT_FAILURE;
  }

  switch (imageIO->GetNumberOfDimensions())
  {
    case 1:
    {
      PrintHeader<1>(argc, argv);
    }
    break;
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
    case 5:
    {
      PrintHeader<5>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension " << imageIO->GetNumberOfDimensions() << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

} // namespace ants
