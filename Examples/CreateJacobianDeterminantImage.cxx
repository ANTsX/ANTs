#include "antsUtilities.h"
#include <algorithm>
#include "ReadWriteData.h"
#include "itkDeformationFieldGradientTensorImageFilter.h"
#include "itkDeterminantTensorImageFilter.h"
#include "itkGeometricJacobianDeterminantImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkMaximumImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension>
int
CreateJacobianDeterminantImage(int argc, char * argv[])
{
  using RealType = double;
  using ImageType = itk::Image<RealType, ImageDimension>;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using VectorImageType = itk::Image<VectorType, ImageDimension>;

  /**
   * Read in vector field
   */
  using ReaderType = itk::ImageFileReader<VectorImageType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  reader->Update();

  typename ImageType::Pointer jacobian = nullptr;

  typename ImageType::Pointer minimumConstantImage = ImageType::New();
  minimumConstantImage->CopyInformation(reader->GetOutput());
  minimumConstantImage->SetRegions(reader->GetOutput()->GetRequestedRegion());
  minimumConstantImage->Allocate();
  minimumConstantImage->FillBuffer(0.001);

  bool calculateLogJacobian = false;
  if (argc > 4)
  {
    calculateLogJacobian = static_cast<bool>(std::stoi(argv[4]));
  }

  bool calculateGeometricJacobian = false;
  if (argc > 5)
  {
    calculateGeometricJacobian = static_cast<bool>(std::stoi(argv[5]));
  }

  bool returnDeformationGradient = false;
  if (argc > 6)
  {
    returnDeformationGradient = static_cast<bool>(std::stoi(argv[6]));
  }

  if (calculateGeometricJacobian)
  {
    using JacobianFilterType = itk::GeometricJacobianDeterminantImageFilter<VectorImageType, RealType, ImageType>;
    typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
    jacobianFilter->SetInput(reader->GetOutput());

    jacobian = jacobianFilter->GetOutput();
    jacobian->Update();
    jacobian->DisconnectPipeline();
  }
  else
  {
    using JacobianFilterType = itk::DeformationFieldGradientTensorImageFilter<VectorImageType, RealType>;
    typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
    jacobianFilter->SetInput(reader->GetOutput());
    jacobianFilter->SetCalculateJacobian(true);
    jacobianFilter->SetUseImageSpacing(true);
    jacobianFilter->SetOrder(2);
    jacobianFilter->SetUseCenteredDifference(true);

    if ( returnDeformationGradient ) {
      jacobianFilter->Update();
      ANTs::WriteImage<typename JacobianFilterType::OutputImageType>(jacobianFilter->GetOutput(), argv[3]);
      return EXIT_SUCCESS;
    }

    using DeterminantFilterType =
      itk::DeterminantTensorImageFilter<typename JacobianFilterType::OutputImageType, RealType>;
    typename DeterminantFilterType::Pointer determinantFilter = DeterminantFilterType::New();
    determinantFilter->SetInput(jacobianFilter->GetOutput());
    determinantFilter->Update();

    minimumConstantImage->FillBuffer(0.0);

    using MaxFilterType = itk::MaximumImageFilter<ImageType, ImageType, ImageType>;
    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1(determinantFilter->GetOutput());
    maxFilter->SetInput2(minimumConstantImage);

    jacobian = maxFilter->GetOutput();
    jacobian->Update();
    jacobian->DisconnectPipeline();
  }

  if (calculateLogJacobian)
  {
    minimumConstantImage->FillBuffer(0.001);

    using MaxFilterType = itk::MaximumImageFilter<ImageType, ImageType, ImageType>;
    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1(jacobian);
    maxFilter->SetInput2(minimumConstantImage);

    using LogFilterType = itk::LogImageFilter<ImageType, ImageType>;
    typename LogFilterType::Pointer logFilter = LogFilterType::New();
    logFilter->SetInput(maxFilter->GetOutput());
    logFilter->Update();
    ANTs::WriteImage<ImageType>(logFilter->GetOutput(), argv[3]);
  }
  else
  {
    ANTs::WriteImage<ImageType>(jacobian, argv[3]);
  }
  return EXIT_SUCCESS;
}


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
CreateJacobianDeterminantImage(std::vector<std::string> args, std::ostream * itkNotUsed(out_stream))
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "CreateJacobianDeterminantImage");

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

  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0]
              << " imageDimension deformationField outputImage [doLogJacobian=0] [useGeometric=0] [deformationGradient=0]" << std::endl;
    std::cout << "deformationGradient cannot be written to nifti: try mhd or nrrd." << std::endl;
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return CreateJacobianDeterminantImage<2>(argc, argv);
    }
    break;
    case 3:
    {
      return CreateJacobianDeterminantImage<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
