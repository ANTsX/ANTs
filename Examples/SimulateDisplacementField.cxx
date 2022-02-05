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
#include "ReadWriteData.h"

#include "itkCastImageFilter.h"
#include "itkSimulatedBSplineDisplacementFieldSource.h"
#include "itkSimulatedExponentialDisplacementFieldSource.h"

#include "string.h"

namespace ants
{

template <unsigned int ImageDimension>
int
SimulateDisplacementField(int argc, char * argv[])
{
  using RealType = float;
  using ImageType = itk::Image<RealType, ImageDimension>;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;

  typename ImageType::Pointer domainImage;
  ReadImage<ImageType>(domainImage, argv[3]);

  if (argc < 5)
  {
    std::cout << argv[0] << " imageDimension displacementFieldType domainImage outputField ";
    std::cout
      << "<numberOfRandomPoints=1000> <standardDeviationDisplacementField=10> <enforceStationaryBoundary=1> [options]"
      << std::endl;
    std::cout
      << "Usage 1 [options] (displacementFieldType = \"BSpline\"): <numberOfFittingLevels=4> <numberOfControlPoints=4>"
      << std::endl;
    std::cout << "Usage 2 [options] (displacementFieldType = \"Exponential\"): <smoothingStandardDeviation=4>"
              << std::endl;
  }

  itk::SizeValueType numberOfRandomPoints = 1000;
  if (argc > 5)
  {
    numberOfRandomPoints = static_cast<itk::SizeValueType>(std::stoi(argv[5]));
  }

  RealType standardDeviationDisplacementField = 10.0;
  if (argc > 6)
  {
    standardDeviationDisplacementField = static_cast<RealType>(atof(argv[6]));
  }

  bool enforceStationaryBoundary = true;
  if (argc > 7)
  {
    enforceStationaryBoundary = static_cast<bool>(std::stoi(argv[7]));
  }

  if (strcmp(argv[2], "BSpline") == 0 || strcmp(argv[2], "bspline") == 0)
  {
    using BSplineSimulatorType = itk::SimulatedBSplineDisplacementFieldSource<DisplacementFieldType>;

    itk::SizeValueType numberOfFittingLevels = 4;
    if (argc > 8)
    {
      numberOfFittingLevels = static_cast<itk::SizeValueType>(std::stoi(argv[8]));
    }

    typename BSplineSimulatorType::ArrayType numberOfControlPoints;

    numberOfControlPoints.Fill(4);
    if (argc > 9)
    {
      std::vector<RealType> ncps = ConvertVector<RealType>(std::string(argv[9]));
      if (ncps.size() == 1)
      {
        numberOfControlPoints.Fill(ncps[0]);
      }
      else if (ncps.size() == ImageDimension)
      {
        for (itk::SizeValueType d = 0; d < ImageDimension; d++)
        {
          numberOfControlPoints[d] = ncps[d];
        }
      }
      else
      {
        std::cerr << "Incorrect specification of the number of control points." << std::endl;
        return EXIT_FAILURE;
      }
    }

    using RealImageType = typename BSplineSimulatorType::RealImageType;
    using CastImageFilterType = itk::CastImageFilter<ImageType, RealImageType>;
    typename CastImageFilterType::Pointer caster = CastImageFilterType::New();
    caster->SetInput(domainImage);
    caster->Update();

    typename BSplineSimulatorType::Pointer bsplineSimulator = BSplineSimulatorType::New();
    bsplineSimulator->SetDisplacementFieldDomainFromImage(caster->GetOutput());
    bsplineSimulator->SetNumberOfRandomPoints(numberOfRandomPoints);
    bsplineSimulator->SetEnforceStationaryBoundary(enforceStationaryBoundary);
    bsplineSimulator->SetDisplacementNoiseStandardDeviation(standardDeviationDisplacementField);
    bsplineSimulator->SetNumberOfFittingLevels(numberOfFittingLevels);
    bsplineSimulator->SetNumberOfControlPoints(numberOfControlPoints);
    bsplineSimulator->Update();

    bsplineSimulator->Print(std::cout, 5);

    ANTs::WriteImage<DisplacementFieldType>(bsplineSimulator->GetOutput(), argv[4]);
  }
  else if (strcmp(argv[2], "Exponential") == 0 || strcmp(argv[2], "exponential") == 0)
  {
    RealType standardDeviationSmoothing = 10.0;
    if (argc > 8)
    {
      standardDeviationSmoothing = static_cast<RealType>(std::stoi(argv[8]));
    }

    using ExponentialSimulatorType = itk::SimulatedExponentialDisplacementFieldSource<DisplacementFieldType>;

    using RealImageType = typename ExponentialSimulatorType::RealImageType;
    using CastImageFilterType = itk::CastImageFilter<ImageType, RealImageType>;
    typename CastImageFilterType::Pointer caster = CastImageFilterType::New();
    caster->SetInput(domainImage);
    caster->Update();

    typename ExponentialSimulatorType::Pointer exponentialSimulator = ExponentialSimulatorType::New();
    exponentialSimulator->SetDisplacementFieldDomainFromImage(caster->GetOutput());
    exponentialSimulator->SetNumberOfRandomPoints(numberOfRandomPoints);
    exponentialSimulator->SetEnforceStationaryBoundary(enforceStationaryBoundary);
    exponentialSimulator->SetDisplacementNoiseStandardDeviation(standardDeviationDisplacementField);
    exponentialSimulator->SetSmoothingStandardDeviation(standardDeviationSmoothing);
    exponentialSimulator->Update();
    ANTs::WriteImage<DisplacementFieldType>(exponentialSimulator->GetOutput(), argv[4]);
  }
  else
  {
    std::cerr << "Unknown displacementFieldType." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
SimulateDisplacementField(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "SimulateDisplacementField");

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

  if (argc < 5)
  {
    std::cout << argv[0] << " imageDimension displacementFieldType domainImage outputField ";
    std::cout
      << "<numberOfRandomPoints=1000> <standardDeviationDisplacementField=10> <enforceStationaryBoundary=1> [options]"
      << std::endl;
    std::cout
      << "Usage 1 [options] (displacementFieldType = \"BSpline\"): <numberOfFittingLevels=4> <numberOfControlPoints=4>"
      << std::endl;
    std::cout << "Usage 2 [options] (displacementFieldType = \"Exponential\"): <smoothingStandardDeviation=4>"
              << std::endl;

    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      return SimulateDisplacementField<2>(argc, argv);
      break;
    case 3:
      return SimulateDisplacementField<3>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

} // namespace ants
