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
#include "itkantsReadWriteTransform.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

namespace ants
{

template <int ImageDimension>
void
ReorientTensorImage(char * moving_image_filename, char * output_image_filename, char * transform_filename)
{
  if (ImageDimension != 3)
  {
    std::cout << "ReorientTensorImage only supports 3D image volumes" << std::endl;
    return;
  }

  using RealType = double;
  using TensorPixelType = itk::DiffusionTensor3D<RealType>;
  using TensorImageType = itk::Image<TensorPixelType, ImageDimension>;
  using ImageType = itk::Image<float, ImageDimension>;

  using AffineTransformType = itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  using ImageFileReaderType = itk::ImageFileReader<ImageType>;
  typename TensorImageType::Pointer img_mov;

  // No reason to use log-euclidean space
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  typename ImageType::Pointer img_ref = nullptr;

  typename ImageFileReaderType::Pointer reader_img_ref = ImageFileReaderType::New();

  using PPDReorientType =
    itk::PreservationOfPrincipalDirectionTensorReorientationImageFilter<TensorImageType>;
  typename PPDReorientType::Pointer reo = PPDReorientType::New();
  reo->SetInput(img_mov);

  using TransformType = itk::Transform<RealType, ImageDimension, ImageDimension>;
  using CompositeTransformType = itk::CompositeTransform<RealType, ImageDimension>;

  // Read the transform from the file
  typename TransformType::Pointer transform = itk::ants::ReadTransform<RealType, ImageDimension>(transform_filename, false);

  // Check if the transform is a composite transform
  if (transform->GetNameOfClass() == "CompositeTransform")
  {
    typename CompositeTransformType::Pointer composite_transform =
        dynamic_cast<CompositeTransformType *>(transform.GetPointer());
    if (composite_transform)
    {
        reo->SetCompositeTransform(composite_transform);
    }
    else
    {
        std::cerr << "Error: Failed to cast to CompositeTransformType." << std::endl;
        return;
    }
  }
  else
  {
    // Create a new composite transform and add the individual transform we just read
    typename CompositeTransformType::Pointer composite_transform = CompositeTransformType::New();
    composite_transform->AddTransform(transform);
    reo->SetCompositeTransform(composite_transform);
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

  if (argc != 5)
  {
    std::cout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <composite.h5/warp.nii.gz/affine.mat/affine.txt> "
        << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  char *         moving_image_filename = nullptr;
  char *         output_image_filename = nullptr;
  char *         transform_filename = nullptr;

  int dim = std::stoi(argv[1]);

  moving_image_filename = argv[2];
  output_image_filename = argv[3];
  transform_filename = argv[4];

  if (dim != 3)
  {
    std::cout << "ReorientTensorImage only supports 3D image volumes" << std::endl;
    return EXIT_FAILURE;
  }

  ReorientTensorImage<3>(moving_image_filename, output_image_filename, transform_filename);

  return EXIT_SUCCESS;
}
} // namespace ants
