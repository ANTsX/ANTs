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

#include "ReadWriteData.h"
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpTensorImageMultiTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

namespace ants
{
// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
RebaseTensorImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "RebaseTensorImage");

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

  if (argc < 5)
  {
    std::cout << "Usage: " << argv[0] << " Dimension infile.nii outfile.nii <PHYSICAL/LOCAL/reference.nii.gz> "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  const int          dim = std::stoi(argv[1]);
  const char * const moving_image_filename = argv[2];
  const char * const output_image_filename = argv[3];

  if (dim != 3)
  {
    std::cout << "RebaseTensorImage only supports 3D image volumes" << std::endl;
    return EXIT_FAILURE;
  }

  using PixelType = itk::DiffusionTensor3D<double>;
  using TensorImageType = itk::Image<PixelType, 3>;
  using ImageType = itk::Image<float, 3>;


  // No reason to use log-euclidean space
  TensorImageType::Pointer img_mov;
  ReadTensorImage<TensorImageType>(img_mov, moving_image_filename, false);

  TensorImageType::DirectionType::InternalMatrixType direction = img_mov->GetDirection().GetVnlMatrix();
  direction.set_identity();

  std::cout << "Transforming space of " << moving_image_filename;

  // std::cout << i << " = " << argv[i-1] << std::endl;
  char * convert = argv[4];

  if (strcmp(convert, "PHYSICAL") == 0)
  {
    std::cout << " -> physical space";
    direction = img_mov->GetDirection().GetVnlMatrix();
  }
  else if (strcmp(convert, "LOCAL") == 0)
  {
    std::cout << " -> local space";
    direction = img_mov->GetDirection().GetTranspose();
  }
  else
  {
    std::cout << " -> " << convert << " space";
    ImageType::Pointer target;
    ReadImage<ImageType>(target, convert);
    // converting from LOCAL to a reference LOCAL space
    direction = target->GetDirection().GetTranspose() * img_mov->GetDirection().GetVnlMatrix();
  }

  // direction = direction.transpose(); // to accomodate for how
  // eigenvectors are stored

  std::cout << std::endl;
  std::cout << "Final rebasing matrix: " << std::endl << direction << std::endl;

  if (!direction.is_identity(0.00001))
  {
    itk::ImageRegionIteratorWithIndex<TensorImageType> it(img_mov, img_mov->GetLargestPossibleRegion());
    while (!it.IsAtEnd())
    {
      /*
      PixelType dt = it.Value();
      PixelType::EigenValuesArrayType evalues;
      PixelType::EigenVectorsMatrixType evectors;
      dt.ComputeEigenAnalysis( evalues, evectors );

      evectors = evectors * direction;

      PixelType::EigenVectorsMatrixType emat;
      emat.Fill( 0.0 );
      for (unsigned int i=0; i<3; i++)
        {
        emat(i,i) = evalues[i];
        }

      PixelType::EigenVectorsMatrixType::InternalMatrixType matrixDT = evectors.GetTranspose() * emat.GetVnlMatrix() *
      evectors.GetVnlMatrix();
      */

      PixelType::EigenVectorsMatrixType::InternalMatrixType dt;
      dt(0, 0) = it.Value()[0];
      dt(0, 1) = dt(1, 0) = it.Value()[1];
      dt(0, 2) = dt(2, 0) = it.Value()[2];
      dt(1, 1) = it.Value()[3];
      dt(1, 2) = dt(2, 1) = it.Value()[4];
      dt(2, 2) = it.Value()[5];

      if ((it.Value()[0] + it.Value()[3] + it.Value()[5]) > 0.00001)
      {
        PixelType::EigenVectorsMatrixType::InternalMatrixType matrixDT = direction * dt * direction.transpose();

        PixelType outDT;
        outDT[0] = matrixDT(0, 0);
        outDT[1] = matrixDT(1, 0);
        outDT[2] = matrixDT(2, 0);
        outDT[3] = matrixDT(1, 1);
        outDT[4] = matrixDT(2, 1);
        outDT[5] = matrixDT(2, 2);
        it.Set(outDT);
      }

      ++it;
    }
  }
  else
  {
    std::cout << "Identity transform detected.. image unmodified" << std::endl;
  }

  // No reason to use log-euclidean space here
  WriteTensorImage<TensorImageType>(img_mov, output_image_filename, false);

  return EXIT_SUCCESS;
}
} // namespace ants
