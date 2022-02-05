#include "antsUtilities.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkDisplacementFieldToBSplineImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"

#include "ReadWriteData.h"

#include <cstdio>
#include <vector>
#include <fstream>
#include <string>

template <typename TValue>
std::vector<TValue>
ConvertDelimitedArray(std::string optionString)
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find(',', 0);

  if (crosspos == std::string::npos)
  {
    values.push_back(ants::Convert<TValue>(optionString));
  }
  else
  {
    std::string element = optionString.substr(0, crosspos);

    TValue             value;
    std::istringstream iss(element);
    iss >> value;
    values.push_back(value);
    while (crosspos != std::string::npos)
    {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find(',', crossposfrom + 1);
      if (crosspos == std::string::npos)
      {
        element = optionString.substr(crossposfrom + 1, optionString.length());
      }
      else
      {
        element = optionString.substr(crossposfrom + 1, crosspos - (crossposfrom + 1));
      }

      std::istringstream iss2(element);
      iss2 >> value;
      values.push_back(value);
    }
  }
  return values;
}

namespace ants
{

template <unsigned int PointDimension>
int
FitBSplineWarpFieldToPoints(unsigned int argc, char * argv[])
{
  using RealType = float;

  using VectorType = itk::Vector<RealType, PointDimension>;
  using DisplacementFieldType = itk::Image<VectorType, PointDimension>;

  using BSplineFilterType = itk::DisplacementFieldToBSplineImageFilter<DisplacementFieldType>;
  using BSplinePointSetType = typename BSplineFilterType::InputPointSetType;

  typename BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();

  typename BSplinePointSetType::Pointer pointSet = BSplinePointSetType::New();
  pointSet->Initialize();

  // Read in points

  std::ifstream inputPointFile(argv[2]);
  std::string   line;

  unsigned int count = 0;
  if (inputPointFile.is_open())
  {
    while (std::getline(inputPointFile, line))
    {
      typename BSplinePointSetType::PointType point;
      float                                   weight = 1.0;

      std::vector<RealType> pointAndWeight = ConvertDelimitedArray<RealType>(line);

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        point[d] = pointAndWeight[d];
      }
      weight = pointAndWeight[pointAndWeight.size() - 1];

      pointSet->SetPoint(count, point);
      weights->InsertElement(count, weight);
      count++;
    }
  }

  unsigned int numberOfPoints = count;

  // Read in displacements
  std::ifstream inputDisplacementsFile(argv[3]);

  count = 0;
  if (inputDisplacementsFile.is_open())
  {
    while (std::getline(inputDisplacementsFile, line))
    {
      VectorType vector;

      std::vector<RealType> displacement = ConvertDelimitedArray<RealType>(line);

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        vector[d] = displacement[d];
      }

      pointSet->SetPointData(count, vector);
      count++;
    }
  }

  unsigned int numberOfPixels = count;

  if (numberOfPixels != numberOfPoints)
  {
    std::cerr << "The number of data does not equal the number of points." << std::endl;
    return EXIT_FAILURE;
  }

  using ImageType = itk::Image<RealType, PointDimension>;

  typename ImageType::Pointer domainImage = nullptr;

  ReadImage<ImageType>(domainImage, argv[4]);
  if (!domainImage)
  {
    std::cerr << "Cannot read image file." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>(std::string(argv[6]));

  constexpr unsigned int splineOrder = 3;

  typename BSplineFilterType::ArrayType ncps;
  if (meshSize.size() == 1)
  {
    ncps.Fill(meshSize[0] + splineOrder);
  }
  else if (meshSize.size() == PointDimension)
  {
    for (unsigned int d = 0; d < PointDimension; d++)
    {
      ncps[d] = meshSize[d] + splineOrder;
    }
  }
  else
  {
    std::cerr << "Invalid ncps format." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int numberOfLevels = 1;
  if (argc > 7)
  {
    numberOfLevels = Convert<unsigned int>(std::string(argv[7]));
  }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetPointSet(pointSet);
  bspliner->SetBSplineDomain(domainImage->GetOrigin(),
                             domainImage->GetSpacing(),
                             domainImage->GetLargestPossibleRegion().GetSize(),
                             domainImage->GetDirection());
  bspliner->SetNumberOfControlPoints(ncps);
  bspliner->SetSplineOrder(splineOrder);
  bspliner->SetNumberOfFittingLevels(numberOfLevels);
  bspliner->SetEnforceStationaryBoundary(true);
  bspliner->SetEstimateInverse(false);
  bspliner->Update();

  ANTs::WriteImage<DisplacementFieldType>(bspliner->GetOutput(), argv[5]);

  return EXIT_SUCCESS;
}

template <unsigned int PointDimension>
int
FitBSplineCurveToPoints(unsigned int argc, char * argv[])
{
  using RealType = float;

  using VectorType = itk::Vector<RealType, PointDimension>;
  using CurveImageType = itk::Image<VectorType, 1>;

  using PointSetType = itk::PointSet<VectorType, 1>;
  typename PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  using FilterType = itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, CurveImageType>;
  typename FilterType::Pointer filter = FilterType::New();

  typename FilterType::WeightsContainerType::Pointer weights = FilterType::WeightsContainerType::New();

  RealType totalDistance = 0.0;

  std::ifstream file(argv[2]);
  std::string   line;

  unsigned int count = 0;
  if (file.is_open())
  {
    while (std::getline(file, line))
    {
      VectorType vector(0.0);
      float      weight = 1.0;

      std::vector<RealType> vectorAndWeight = ConvertDelimitedArray<RealType>(line);

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        vector[d] = vectorAndWeight[d];
      }
      weight = vectorAndWeight[vectorAndWeight.size() - 1];

      pointSet->SetPointData(count, vector);

      if (count > 0)
      {
        VectorType previous(0.0);
        pointSet->GetPointData(count - 1, &previous);
        totalDistance += static_cast<RealType>((previous - vector).GetNorm());
      }

      typename PointSetType::PointType point;
      point[0] = 0.0;
      pointSet->SetPoint(count, point);

      weights->InsertElement(count, weight);
      count++;
    }
  }

  RealType cumSum = 0.0;
  for (unsigned int i = 1; i < pointSet->GetNumberOfPoints(); i++)
  {
    VectorType vector(0.0), previous(0.0);
    pointSet->GetPointData(i, &vector);
    pointSet->GetPointData(i - 1, &previous);

    cumSum += static_cast<RealType>((vector - previous).GetNorm());
    typename PointSetType::PointType point;
    point[0] = cumSum / totalDistance;

    pointSet->SetPoint(i, point);
  }

  filter->SetInput(pointSet);
  filter->SetGenerateOutputImage(true);

  typename CurveImageType::PointType origin;
  origin.Fill(0.0);
  filter->SetOrigin(origin);
  typename CurveImageType::SpacingType spacing;
  spacing[0] = 0.001;
  if (argc > 6)
  {
    spacing[0] = atof(argv[6]);
  }

  filter->SetSpacing(spacing);
  typename CurveImageType::SizeType size;
  size[0] = static_cast<unsigned int>(1.0 / spacing[0] + 1);
  filter->SetSize(size);
  typename FilterType::ArrayType order;
  order[0] = 3;
  if (argc > 3)
  {
    order[0] = std::stoi(argv[3]);
  }
  filter->SetSplineOrder(order);
  typename FilterType::ArrayType ncps;
  ncps[0] = order[0] + 1;
  if (argc > 5)
  {
    ncps[0] = std::stoi(argv[5]);
  }
  filter->SetNumberOfControlPoints(ncps);
  typename FilterType::ArrayType nlevels;
  nlevels[0] = 5;
  if (argc > 4)
  {
    nlevels[0] = std::stoi(argv[4]);
  }
  filter->SetNumberOfLevels(nlevels);
  typename FilterType::ArrayType close;
  close[0] = false;
  if (argc > 7)
  {
    close[0] = std::stoi(argv[7]);
  }
  filter->SetCloseDimension(close);

  filter->Update();

  itk::ImageRegionIterator<CurveImageType> It(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    VectorType vector = It.Get();
    for (unsigned int d = 0; d < PointDimension - 1; d++)
    {
      std::cout << vector[d] << ",";
    }
    std::cout << vector[PointDimension - 1] << std::endl;
  }

  //   {
  //   std::string filename = std::string( argv[2] ) + std::string( "_cps.txt" );
  //   std::ofstream ostr( filename.c_str() );
  //   ostr << "0 0 0 0" << std::endl;
  //
  //   itk::ImageRegionIterator<CurveImageType> It(
  //     filter->GetPhiLattice(), filter->GetPhiLattice()->GetLargestPossibleRegion() );
  //   for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
  //     {
  //     ostr << It.Get()[0] << " " << It.Get()[1] << " " << It.Get()[2] << " 1" << std::endl;
  //     }
  //   ostr << "0 0 0 0" << std::endl;
  //   ostr.close();
  //   }

  return EXIT_SUCCESS;
}


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
FitBSplineToPoints(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "FitBSplineToPoints");

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
    std::cout << "Usage: " << argv[0] << " pointDimension inputLandmarksFile "
              << " [order=3] [nlevels=10] "
              << " [numberOfControlPoints=4] [sampleSpacing=0.001] [closed?=0]" << std::endl;
    std::cout << "Usage 2: " << argv[0] << " pointDimension inputPointFile inputDisplacementFile "
              << " domainImage outputDisplacementField controlPointMeshSize [nlevels=1]" << std::endl;
    std::cout << "  Note:  1. Points are assumed to be parametrically ordered for fitting to a curve. " << std::endl
              << "         2. The last column (pointDimension+1) is used for weights." << std::endl
              << "         3. To specify a warp field, add a 'w' to the dimension argument, e.g. 2w." << std::endl;
    return EXIT_FAILURE;
  }

  std::string imageDimensionString(argv[1]);
  if (imageDimensionString.length() == 2)
  {
    if (imageDimensionString[1] == 'w')
    {
      switch (std::stoi(&imageDimensionString[0]))
      {
        case 2:
          return FitBSplineWarpFieldToPoints<2>(argc, argv);
          break;
        case 3:
          return FitBSplineWarpFieldToPoints<3>(argc, argv);
          break;
        default:
          std::cerr << "Unsupported dimension" << std::endl;
          return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Incorrect dimension specification.  See help." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    switch (std::stoi(argv[1]))
    {
      case 1:
        return FitBSplineCurveToPoints<2>(argc, argv);
        break;
      case 2:
        return FitBSplineCurveToPoints<2>(argc, argv);
        break;
      case 3:
        return FitBSplineCurveToPoints<3>(argc, argv);
        break;
      case 4:
        return FitBSplineCurveToPoints<4>(argc, argv);
        break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
} // namespace ants
