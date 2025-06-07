/** ANTS Landmarks used to initialize an b-spline displacement field ... */

#include "antsUtilities.h"

#include "itkAffineTransform.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkRigid2DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkTransformFileWriter.h"
#include "itkFloatTypes.h"

#include <string>
#include <vector>
#include <algorithm>

namespace ants
{

template <typename ImageType, typename PointSetType>
void
ReadLabeledPointSetFromImage(typename ImageType::Pointer                  image,
                             typename PointSetType::Pointer               pointSet,
                             std::vector<typename ImageType::PixelType> & labels)
{
  labels.clear();

  typename PointSetType::Pointer allPoints = PointSetType::New();
  allPoints->Initialize();

  itk::ImageRegionIteratorWithIndex<ImageType> It(image, image->GetLargestPossibleRegion());

  unsigned int count = 0;
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    typename ImageType::PixelType value = It.Get();
    if (value != 0)
    {
      if (std::find(labels.begin(), labels.end(), value) == labels.end())
      {
        labels.push_back(value);
      }
      typename PointSetType::PointType point;
      image->TransformIndexToPhysicalPoint(It.GetIndex(), point);
      allPoints->SetPointData(count, value);
      allPoints->SetPoint(count++, point);
    }
  }
  std::sort(labels.begin(), labels.end());

  pointSet->Initialize();

  for (unsigned int n = 0; n < labels.size(); n++)
  {
    typename ImageType::PixelType    currentLabel = labels[n];
    typename PointSetType::PointType center;
    center.Fill(0);
    float                                               N = 0;
    typename PointSetType::PointsContainerConstIterator ItP = allPoints->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator   ItD = allPoints->GetPointData()->Begin();
    while (ItP != allPoints->GetPoints()->End())
    {
      if (ItD.Value() == currentLabel)
      {
        typename PointSetType::PointType point = ItP.Value();
        for (unsigned int d = 0; d < ImageType::ImageDimension; d++)
        {
          center[d] += point[d];
        }
        N += 1.0f;
      }
      ++ItP;
      ++ItD;
    }
    for (unsigned int d = 0; d < ImageType::ImageDimension; d++)
    {
      center[d] /= N;
    }
    pointSet->SetPoint(n, center);
    pointSet->SetPointData(n, currentLabel);
  }
}

template <unsigned int ImageDimension, typename TransformType>
int
InitializeLinearTransform(int itkNotUsed(argc), char * argv[])
{
  using LabelType = unsigned int;
  using LabelImageType = itk::Image<LabelType, ImageDimension>;

  using TransformInitializerType =
    itk::LandmarkBasedTransformInitializer<TransformType, LabelImageType, LabelImageType>;
  using LandmarkContainerType = typename TransformInitializerType::LandmarkPointContainer;

  using PointSetType = itk::PointSet<LabelType, ImageDimension>;

  //
  // Read in the fixed and moving images and convert to point sets
  //

  using ImageReaderType = itk::ImageFileReader<LabelImageType>;
  typename ImageReaderType::Pointer fixedReader = ImageReaderType::New();
  fixedReader->SetFileName(argv[2]);
  fixedReader->Update();
  typename LabelImageType::Pointer fixedImage = fixedReader->GetOutput();

  typename PointSetType::Pointer fixedPoints = PointSetType::New();
  std::vector<LabelType>         fixedLabels;
  ReadLabeledPointSetFromImage<LabelImageType, PointSetType>(fixedImage, fixedPoints, fixedLabels);

  typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
  movingReader->SetFileName(argv[3]);
  movingReader->Update();
  typename LabelImageType::Pointer movingImage = movingReader->GetOutput();

  typename PointSetType::Pointer movingPoints = PointSetType::New();
  std::vector<LabelType>         movingLabels;
  ReadLabeledPointSetFromImage<LabelImageType, PointSetType>(movingImage, movingPoints, movingLabels);


  LandmarkContainerType                               fixedLandmarks;
  typename PointSetType::PointsContainerConstIterator ItF = fixedPoints->GetPoints()->Begin();
  while (ItF != fixedPoints->GetPoints()->End())
  {
    fixedLandmarks.push_back(ItF.Value());
    ++ItF;
  }


  LandmarkContainerType                               movingLandmarks;
  typename PointSetType::PointsContainerConstIterator ItM = movingPoints->GetPoints()->Begin();
  while (ItM != movingPoints->GetPoints()->End())
  {
    movingLandmarks.push_back(ItM.Value());
    ++ItM;
  }

  if (fixedLandmarks.size() != movingLandmarks.size())
  {
    std::cerr << "The number of fixed points and moving points must be the same." << std::endl;
    return EXIT_FAILURE;
  }

  typename std::vector<LabelType>::const_iterator itf;
  for (itf = fixedLabels.begin(); itf != fixedLabels.end(); ++itf)
  {
    if (std::find(movingLabels.begin(), movingLabels.end(), *itf) == movingLabels.end())
    {
      std::cerr << "Labels do not match." << std::endl;
      return EXIT_FAILURE;
    }
  }

  //
  // Calculate the transform
  //

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetFixedLandmarks(fixedLandmarks);
  initializer->SetMovingLandmarks(movingLandmarks);
  initializer->SetTransform(transform);
  initializer->InitializeTransform();

  //
  // Write the transform
  //

  typename itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
  transformWriter->SetFileName(argv[5]);
  transformWriter->SetInput(transform);
#if ITK_VERSION_MAJOR >= 5
  transformWriter->SetUseCompression(true);
#endif

  try
  {
    transformWriter->Update();
  }
  catch (const itk::ExceptionObject & itkNotUsed(err))
  {
    std::cerr << "Exception in writing transform file: " << argv[5] << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int
InitializeBSplineTransform(int argc, char * argv[])
{
  using RealType = float;
  using LabelType = unsigned int;
  using LabelImageType = itk::Image<LabelType, ImageDimension>;

  using PointSetType = itk::PointSet<LabelType, ImageDimension>;

  //
  // Read in the fixed and moving images and convert to point sets
  //

  using ImageReaderType = itk::ImageFileReader<LabelImageType>;
  typename ImageReaderType::Pointer fixedReader = ImageReaderType::New();
  fixedReader->SetFileName(argv[2]);
  fixedReader->Update();
  typename LabelImageType::Pointer fixedImage = fixedReader->GetOutput();

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItF(fixedImage, fixedImage->GetLargestPossibleRegion());

  typename PointSetType::Pointer  fixedLandmarks = PointSetType::New();
  typename std::vector<LabelType> fixedLabels;

  ReadLabeledPointSetFromImage<LabelImageType, PointSetType>(fixedImage, fixedLandmarks, fixedLabels);

  typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
  movingReader->SetFileName(argv[3]);
  movingReader->Update();
  typename LabelImageType::Pointer movingImage = movingReader->GetOutput();

  typename PointSetType::Pointer  movingLandmarks = PointSetType::New();
  typename std::vector<LabelType> movingLabels;
  ReadLabeledPointSetFromImage<LabelImageType, PointSetType>(movingImage, movingLandmarks, movingLabels);

  if (fixedLandmarks->GetNumberOfPoints() != movingLandmarks->GetNumberOfPoints())
  {
    std::cerr << "The number of fixed points and moving points must be the same." << std::endl;
    return EXIT_FAILURE;
  }

  typename std::vector<LabelType>::const_iterator itf;
  for (itf = fixedLabels.begin(); itf != fixedLabels.end(); ++itf)
  {
    if (std::find(movingLabels.begin(), movingLabels.end(), *itf) == movingLabels.end())
    {
      std::cerr << "Labels do not match." << std::endl;
      return EXIT_FAILURE;
    }
  }

  //
  // Calculate the transform
  //

  typename LabelImageType::DirectionType fixedDirection = fixedImage->GetDirection();
  typename LabelImageType::DirectionType fixedDirectionInverse(fixedDirection.GetInverse());

  typename LabelImageType::DirectionType identityDirection;
  identityDirection.SetIdentity();

  const typename LabelImageType::RegionType & bufferedRegion = fixedImage->GetBufferedRegion();
  const itk::SizeValueType                    numberOfPixels = bufferedRegion.GetNumberOfPixels();

  const bool filterHandlesMemory = false;

  using ImporterType = itk::ImportImageFilter<LabelType, ImageDimension>;
  typename ImporterType::Pointer importer = ImporterType::New();
  importer->SetImportPointer(
    const_cast<LabelType *>(fixedImage->GetBufferPointer()), numberOfPixels, filterHandlesMemory);
  importer->SetRegion(fixedImage->GetBufferedRegion());
  importer->SetOrigin(fixedImage->GetOrigin());
  importer->SetSpacing(fixedImage->GetSpacing());
  importer->SetDirection(identityDirection);
  importer->Update();

  const typename ImporterType::OutputImageType * parametricInputImage = importer->GetOutput();

  using VectorType = itk::Vector<RealType, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;

  // Read in the optional label weights

  std::vector<float>     labelWeights;
  std::vector<LabelType> userLabels;

  bool useWeights = false;

  unsigned int labelCount = 0;
  if (argc > 10)
  {
    useWeights = true;

    std::fstream labelStr(argv[10]);

    if (labelStr.is_open())
    {
      while (!labelStr.eof())
      {
        char line[256];
        labelStr.getline(line, 256);

        std::string lineString = std::string(line);
        std::size_t pos = lineString.find(',');

        RealType value;
        if (pos == std::string::npos)
        {
          std::istringstream iss(lineString);
          iss >> value;
          labelWeights.push_back(value);
          userLabels.push_back(movingLabels[labelCount++]);
        }
        else
        {
          unsigned int localLabel;

          std::string        element = lineString.substr(0, pos);
          std::istringstream iss(element);
          iss >> localLabel;
          userLabels.push_back(localLabel);

          element = lineString.substr(pos + 1, lineString.length());
          std::istringstream iss2(element);
          iss2 >> value;
          labelWeights.push_back(value);
        }
      }

      labelStr.close();
    }
    else
    {
      std::cerr << "File " << argv[10] << " cannot be opened." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Now match up the center points

  using DisplacementFieldPointSetType = itk::PointSet<VectorType, ImageDimension>;
  using BSplineFilterType =
    itk::BSplineScatteredDataPointSetToImageFilter<DisplacementFieldPointSetType, DisplacementFieldType>;
  using WeightsContainerType = typename BSplineFilterType::WeightsContainerType;

  typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
  weights->Initialize();
  const typename WeightsContainerType::Element boundaryWeight = 1.0e10;
  typename WeightsContainerType::Element       weight = 1.0;

  typename DisplacementFieldPointSetType::Pointer fieldPoints = DisplacementFieldPointSetType::New();
  fieldPoints->Initialize();
  unsigned long count = 0;

  typename PointSetType::PointsContainerConstIterator mIt = movingLandmarks->GetPoints()->Begin();
  typename PointSetType::PointDataContainerIterator   mItD = movingLandmarks->GetPointData()->Begin();

  while (mItD != movingLandmarks->GetPointData()->End())
  {
    typename PointSetType::PointsContainerConstIterator fIt = fixedLandmarks->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator   fItD = fixedLandmarks->GetPointData()->Begin();

    while (fItD != fixedLandmarks->GetPointData()->End())
    {
      if (fItD.Value() == mItD.Value())
      {
        typename PointSetType::PointType fpoint = fIt.Value();
        typename PointSetType::PointType mpoint = mIt.Value();

        VectorType vector;

        typename LabelImageType::PointType fixedPhysicalPoint;
        for (unsigned int i = 0; i < ImageDimension; i++)
        {
          fixedPhysicalPoint[i] = fpoint[i];
          vector[i] = mpoint[i] - fpoint[i];
        }

        itk::ContinuousIndex<double, ImageDimension> fixedCidx;
        fixedCidx = fixedImage->
                template TransformPhysicalPointToContinuousIndex<double, itk::SpacePrecisionType>(fixedPhysicalPoint);

        typename DisplacementFieldType::PointType fieldPoint;
        fieldPoint = parametricInputImage->
                template TransformContinuousIndexToPhysicalPoint<double, itk::SpacePrecisionType>(fixedCidx);

        fieldPoints->SetPoint(count, fieldPoint);
        fieldPoints->SetPointData(count, vector);

        if (useWeights)
        {
          auto it = std::find(userLabels.begin(), userLabels.end(), mItD.Value());
          if (it != userLabels.end())
          {
            weights->InsertElement(count, labelWeights[it - userLabels.begin()]);
          }
          else
          {
            std::cerr << "Unspecified label " << mItD.Value() << " in specified user label weights." << std::endl;
            return EXIT_FAILURE;
          }
        }
        else
        {
          weights->InsertElement(count, weight);
        }

        count++;

        break;
      }
      ++fItD;
      ++fIt;
    }

    ++mItD;
    ++mIt;
  }

  bool enforceStationaryBoundary = true;
  if (argc > 9)
  {
    enforceStationaryBoundary = static_cast<bool>(std::stoi(argv[9]));
  }
  if (enforceStationaryBoundary)
  {
    typename LabelImageType::IndexType startIndex2 = fixedImage->GetLargestPossibleRegion().GetIndex();

    typename LabelImageType::SizeType inputSize2 = fixedImage->GetLargestPossibleRegion().GetSize();
    for (ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF)
    {
      typename LabelImageType::IndexType index = ItF.GetIndex();

      bool isOnStationaryBoundary = false;
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        if (index[d] == startIndex2[d] || index[d] == startIndex2[d] + static_cast<int>(inputSize2[d]) - 1)
        {
          isOnStationaryBoundary = true;
          break;
        }
      }

      if (isOnStationaryBoundary)
      {
        VectorType vector;

        vector.Fill(0.0);

        typename PointSetType::PointType fixedPoint;
        parametricInputImage->TransformIndexToPhysicalPoint(index, fixedPoint);

        fieldPoints->SetPoint(count, fixedPoint);
        fieldPoints->SetPointData(count, vector);
        weights->InsertElement(count, boundaryWeight);
        count++;
      }
    }
  }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int numberOfLevels = 4;
  if (argc > 7)
  {
    numberOfLevels = std::stoi(argv[7]);
  }

  unsigned int splineOrder = 3;
  if (argc > 8)
  {
    splineOrder = std::stoi(argv[8]);
  }


  typename BSplineFilterType::ArrayType ncps;
  ncps.Fill(1 + splineOrder);

  if (argc > 6)
  {
    std::vector<unsigned int> meshSize = ConvertVector<unsigned int>(std::string(argv[6]));

    if (meshSize.size() == 1)
    {
      ncps.Fill(meshSize[0] + splineOrder);
    }
    else if (meshSize.size() == ImageDimension)
    {
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        ncps[d] = meshSize[d] + splineOrder;
      }
    }
    else
    {
      std::cerr << "Invalid meshSize format." << std::endl;
    }
  }

  bspliner->SetOrigin(fixedImage->GetOrigin());
  bspliner->SetSpacing(fixedImage->GetSpacing());
  bspliner->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  bspliner->SetDirection(fixedImage->GetDirection());
  bspliner->SetGenerateOutputImage(true);
  bspliner->SetNumberOfLevels(numberOfLevels);
  bspliner->SetSplineOrder(splineOrder);
  bspliner->SetNumberOfControlPoints(ncps);
  bspliner->SetInput(fieldPoints);
  bspliner->SetPointWeights(weights);
  bspliner->Update();

  using WriterType = itk::ImageFileWriter<DisplacementFieldType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[5]);
  writer->SetInput(bspliner->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsLandmarkBasedTransformInitializer(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsLandmarkBasedTransformInitializer");
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

  if (argc < 6)
  {
    std::cerr << "Usage: " << argv[0] << " dimension fixedImage movingImage transformType "
              << "outputTransform [meshSize[0]xmeshSize[1]x...=1x1x1] [numberOfLevels=4] [order=3] "
              << "[enforceStationaryBoundaries=1] [landmarkWeights] " << std::endl;
    std::cerr << std::endl << "Notes:" << std::endl;
    std::cerr << R"( 1) transformType = "rigid", "affine", or "bspline".)" << std::endl;
    std::cerr << " 2) The optional arguments only apply to the bspline transform." << std::endl;
    std::cerr << " 3) The landmark weights are read from a text file where each row is either:" << std::endl;
    std::cerr << "     \"label,labelWeight\" or " << std::endl;
    std::cerr << "     \"labelWeight\" " << std::endl;
    std::cerr
      << "    If the latter format is used, the label weights are assumed to be arranged in ascending order by label."
      << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  const unsigned int dimension = static_cast<unsigned int>(std::stoi(argv[1]));

  using Rigid2DTransformType = itk::Rigid2DTransform<double>;
  using Rigid3DTransformType = itk::VersorRigid3DTransform<double>;
  using AffineTransform2DType = itk::AffineTransform<double, 2>;
  using AffineTransform3DType = itk::AffineTransform<double, 3>;

  std::string transformType(argv[4]);
  ConvertToLowerCase(transformType);

  switch (dimension)
  {
    case 2:
    {
      if (std::strcmp(transformType.c_str(), "affine") == 0)
      {
        return InitializeLinearTransform<2, AffineTransform2DType>(argc, argv);
      }
      else if (std::strcmp(transformType.c_str(), "rigid") == 0)
      {
        return InitializeLinearTransform<2, Rigid2DTransformType>(argc, argv);
      }
      else if (std::strcmp(transformType.c_str(), "bspline") == 0)
      {
        return InitializeBSplineTransform<2>(argc, argv);
      }
      else
      {
        std::cerr << "Unrecognized transform:  " << transformType.c_str() << std::endl;
      }
    }
    break;
    case 3:
    {
      if (std::strcmp(transformType.c_str(), "affine") == 0)
      {
        return InitializeLinearTransform<3, AffineTransform3DType>(argc, argv);
      }
      else if (std::strcmp(transformType.c_str(), "rigid") == 0)
      {
        return InitializeLinearTransform<3, Rigid3DTransformType>(argc, argv);
      }
      else if (std::strcmp(transformType.c_str(), "bspline") == 0)
      {
        return InitializeBSplineTransform<3>(argc, argv);
      }
      else
      {
        std::cerr << "Unrecognized transform:  " << transformType.c_str() << std::endl;
      }
    }
    break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
