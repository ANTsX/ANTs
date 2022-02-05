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
#include "antsUtilities.h"

#include "ReadWriteData.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

namespace ants
{
static bool
ReorientTensorImage_ParseInput(int              argc,
                               char **          argv,
                               char *&          moving_image_filename,
                               char *&          output_image_filename,
                               TRAN_OPT_QUEUE & opt_queue)
{
  opt_queue.clear();
  opt_queue.reserve(argc - 2);

  moving_image_filename = argv[0];
  output_image_filename = argv[1];

  int ind = 2;

  while (ind < argc)
  {
    TRAN_OPT opt;
    opt.filename = argv[ind];
    opt.file_type = CheckFileType(opt.filename.c_str());
    opt.do_affine_inv = false;

    if (strcmp(argv[ind], "-i") == 0)
    {
      std::cout << "ERROR - inverse transforms not yet supported\n" << std::endl;
      return false;
    }
    else
    {
      bool set_current_affine_inv = false;
      if (opt.file_type == AFFINE_FILE)
      {
        SetAffineInvFlag(opt, set_current_affine_inv);
      }
      else
      {
        if (opt.file_type == DEFORMATION_FILE && set_current_affine_inv)
        {
          std::cout << "Ignore inversion of non-affine file type! " << std::endl;
          std::cout << "opt.do_affine_inv:" << opt.do_affine_inv << std::endl;
        }
      }

      opt_queue.push_back(opt);
      DisplayOpt(opt);
    }
    ++ind;
  }

  return true;
}

template <int ImageDimension>
void
ReorientTensorImage(char * moving_image_filename, char * output_image_filename, TRAN_OPT_QUEUE & opt_queue)
{
  using PixelType = itk::DiffusionTensor3D<double>;
  using TensorImageType = itk::Image<PixelType, ImageDimension>;
  using ImageType = itk::Image<float, ImageDimension>;
  using VectorType = itk::Vector<double, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  using ImageFileReaderType = itk::ImageFileReader<ImageType>;
  typename TensorImageType::Pointer img_mov;

  // No reason to use log-euclidean space
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  typename ImageType::Pointer img_ref = nullptr;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  using TranReaderType = itk::TransformFileReader;
  using FieldReaderType = itk::ImageFileReader<DisplacementFieldType>;
  typename DisplacementFieldType::Pointer field = nullptr;
  typename AffineTransformType::Pointer   aff = nullptr;

  const int kOptQueueSize = opt_queue.size();

  if (kOptQueueSize > 1)
  {
    std::cout << "ERROR: Only 1 input transform is permitted" << std::endl;
    return;
  }

  using PPDReorientType =
    itk::PreservationOfPrincipalDirectionTensorReorientationImageFilter<TensorImageType, DisplacementFieldType>;
  typename PPDReorientType::Pointer reo = PPDReorientType::New();
  reo->SetInput(img_mov);

  const TRAN_OPT & opt = opt_queue[0];

  switch (opt.file_type)
  {
    case AFFINE_FILE:
    {
      typename TranReaderType::Pointer tran_reader = TranReaderType::New();
      tran_reader->SetFileName(opt.filename);
      tran_reader->Update();
      aff = dynamic_cast<AffineTransformType *>((tran_reader->GetTransformList())->front().GetPointer());
      reo->SetAffineTransform(aff);
      std::cout << "Affine transform" << std::endl;
    }
    break;
    case DEFORMATION_FILE:
    {
      typename FieldReaderType::Pointer field_reader = FieldReaderType::New();
      field_reader->SetFileName(opt.filename);
      field_reader->Update();
      // field = field_reader->GetOutput();
      reo->SetDisplacementField(field_reader->GetOutput());
      std::cout << "Warp transform" << std::endl;
    }
    break;
    default:
    {
      std::cout << "Unknown file type!" << std::endl;
    }
  }
  reo->Update();

  typename TensorImageType::Pointer img_output = reo->GetOutput();
  // No reason to use log-euclidean space here
  WriteTensorImage<TensorImageType>(img_output, output_image_filename, false);
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ReorientTensorImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ReorientTensorImage");

  int     argc = args.size();
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
  // antscout->set_stream( out_stream );

  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <warp.nii/affine.txt> " << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  TRAN_OPT_QUEUE opt_queue;
  char *         moving_image_filename = nullptr;
  char *         output_image_filename = nullptr;

  bool is_parsing_ok = false;
  int  dim = std::stoi(argv[1]);

  if (dim != 3)
  {
    std::cout << "ReorientTensorImage only supports 3D image volumes" << std::endl;
    return EXIT_FAILURE;
  }

  is_parsing_ok =
    ReorientTensorImage_ParseInput(argc - 2, argv + 2, moving_image_filename, output_image_filename, opt_queue);

  if (is_parsing_ok)
  {
    std::cout << "moving_image_filename: " << moving_image_filename << std::endl;
    std::cout << "output_image_filename: " << output_image_filename << std::endl;
    DisplayOptQueue(opt_queue);

    ReorientTensorImage<3>(moving_image_filename, output_image_filename, opt_queue);
  }
  else
  {
    std::cout << "Input error!" << std::endl;
    return EXIT_FAILURE;
  }
  // ReorientTensorImage<3>(argc,argv);
  // WarpImageForward(argc,argv);
  return EXIT_SUCCESS;
}
} // namespace ants
