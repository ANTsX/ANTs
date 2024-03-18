/*=========================================================================1

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

#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkMacro.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkCommonEnums.h>

#include <itkMetaDataDictionary.h>

#include <iostream>
#include <sys/stat.h>

#include <fstream>
#include <cstdio>

#include <cstring>
#include <sstream>


namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
DeNrrd(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "DeNrrd");

  const int argc = args.size();
  char **   argv = new char *[args.size() + 1];
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

  // antscout->set_stream( out_stream );

  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " inImage.nrrd outImage.nii gradients.txt" << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  const char * const input_image_filename = argv[1];
  const char * const output_image_filename = argv[2];
  const char * const output_gradients_filename = argv[3];

  using PixelType = float;
  using DiffusionImageType = itk::VectorImage<PixelType, 3>;

  using FileReaderType = itk::ImageFileReader<DiffusionImageType, itk::DefaultConvertPixelTraits<PixelType>>;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(input_image_filename);
  reader->Update();

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  // io->SetNrrdVectorType( nrrdKindList );
  io->SetFileType(itk::CommonEnums::IOFile::ASCII);

  // std::cout  << "Number of dwi values " << reader->GetOutput()->GetNumberOfComponentsPerPixel() << std::endl;

  std::string v_string;
  std::string b_string;

  itk::MetaDataDictionary & mdd = reader->GetOutput()->GetMetaDataDictionary();

  if (mdd.HasKey("modality"))
  {
    std::string modality;
    itk::ExposeMetaData<std::string>(mdd, "modality", modality);
    if (modality.compare("DWMRI") != 0)
    {
      std::cerr << "Data in not DWMRI" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << "Error: No b-value found" << std::endl;
    return EXIT_FAILURE;
  }


  if (mdd.HasKey("DWMRI_b-value"))
  {
    itk::ExposeMetaData<std::string>(mdd, "DWMRI_b-value", b_string);
    // std::cout << "BValue = " << b_string << std::endl;
  }
  else
  {
    std::cerr << "Error: No b-value found" << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream gradientfile;
  gradientfile.open(output_gradients_filename);
  gradientfile << "VERSION: 2" << std::endl;

  for (unsigned short i = 0; i < static_cast<unsigned short>(reader->GetOutput()->GetNumberOfComponentsPerPixel()); i++)
  {
    char gradKey[40];
    sprintf(gradKey, "DWMRI_gradient_%04hu", i);
    itk::ExposeMetaData<std::string>(mdd, gradKey, v_string);
    // std::cout << "Gradient = " << v_string << std::endl;

    std::istringstream iss(v_string);
    double             x, y, z;
    iss >> x >> y >> z;

    if (itk::Math::FloatAlmostEqual(x * x + y * y + z * z, itk::NumericTraits<double>::ZeroValue()))
    {
      gradientfile << v_string << " 0" << std::endl;
    }
    else
    {
      gradientfile << v_string << " " << b_string << std::endl;
    }
  }

  gradientfile.close();

  /*
  OutputImageType::Pointer outImage = OutputImageType::New();
  OutputImageType::SizeType outSize;
  outSize[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  outSize[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  outSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
  outSize[3] = reader->GetOutput()->GetNumberOfComponentsPerPixel();
  OutputImageType::RegionType outRegion;
  outRegion.SetSize( outSize );
  outImage->SetRegions( outRegion );
  outImage->AllocateInitialized();

  OutputImageType::SpacingType spacing;
  spacing[0] = reader->GetOutput()->GetSpacing()[0];
  spacing[1] = reader->GetOutput()->GetSpacing()[1];
  spacing[2] = reader->GetOutput()->GetSpacing()[2];
  spacing[3] = 1.0;
  outImage->SetSpacing(spacing);

  OutputImageType::PointType origin;
  origin[0] = reader->GetOutput()->GetOrigin()[0];
  origin[1] = reader->GetOutput()->GetOrigin()[1];
  origin[2] = reader->GetOutput()->GetOrigin()[2];
  outImage->SetOrigin(origin);

  OutputImageType::DirectionType dirMat;
  for ( unsigned int i=0; i<3; i++ )
    for ( unsigned int j=0; j<3; j++)
      dirMat(i,j) = reader->GetOutput()->GetDirection()(i,j);
  dirMat(3,3) = 1.0;
  outImage->SetDirection(dirMat);

  for ( unsigned int x=0; x<outSize[0]; x++ )
    for ( unsigned int y=0; y<outSize[1]; y++ )
      for ( unsigned int z=0; z<outSize[2]; z++ )
        for ( unsigned int t=0; t<outSize[3]; t++ )
          {
          OutputImageType::IndexType oIdx;
          oIdx[0] = x;
          oIdx[1] = y;
          oIdx[2] = z;
          oIdx[3] = t;

          DiffusionImageType::IndexType dIdx;
          dIdx[0] = x;
          dIdx[1] = y;
          dIdx[2] = z;

          outImage->SetPixel(oIdx,reader->GetOutput()->GetPixel(dIdx)[t]);

          }
*/


  using FileWriterType = itk::ImageFileWriter<DiffusionImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName(output_image_filename);
  writer->SetInput(reader->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}
} // namespace ants
