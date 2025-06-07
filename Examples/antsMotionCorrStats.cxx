/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "ReadWriteData.h"
#include "antsCommandLineParser.h"

#include "itkCSVNumericObjectFileWriter.h"
#include "itkCSVArray2DFileReader.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"
#include <sstream>

namespace ants
{


int
ants_motion_stats(itk::ants::CommandLineParser * parser)
{

  constexpr unsigned int ImageDimension = 3;

  using RealType = double;

  using ImageType = itk::Image<RealType, ImageDimension>;
  using vMatrix = vnl_matrix<RealType>;
  vMatrix param_values;
  using TimeSeriesImageType = itk::Image<RealType, 4>;
  using TimeSeriesRegionType = TimeSeriesImageType::RegionType;
  using TimeSeriesIndexType = TimeSeriesImageType::IndexType;
  using TimeSeriesSizeType = TimeSeriesImageType::SizeType;
  using TimeSeriesIteratorType = itk::ImageRegionIterator<TimeSeriesImageType>;

  using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
  using AffineTransformType = itk::AffineTransform<RealType, ImageDimension>;
  using RigidTransformType = itk::Euler3DTransform<RealType>;
  using TransformWriterType = itk::TransformFileWriterTemplate<RealType>;

  using ParserType = itk::ants::CommandLineParser;
  using OptionType = ParserType::OptionType;

  using ParameterValueType = double;
  using MocoReaderType = itk::CSVArray2DFileReader<ParameterValueType>;
  using MocoDataArrayType = MocoReaderType::Array2DDataObjectType;

  using WriterType = itk::CSVNumericObjectFileWriter<ParameterValueType>;
  using WriterMatrixType = WriterType::vnlMatrixType;


  std::string   outputName = "";
  std::string   mocoName = "";
  std::string   spatialName = "";
  std::string   timeseriesDisplacementName = "";
  unsigned long transformIndex = 0;
  bool          writeMap = false;
  bool          writeTransform = false;
  bool          timeseriesDisplacement = false;

  OptionType::Pointer outputOption = parser->GetOption("output");
  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    outputName = outputOption->GetFunction(0)->GetName();
    std::cout << "Output: " << outputName << std::endl;
  }
  else
  {
    std::cerr << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
  }

  OptionType::Pointer spatialOption = parser->GetOption("spatial-map");
  if (spatialOption && spatialOption->GetNumberOfFunctions())
  {
    spatialName = spatialOption->GetFunction(0)->GetName();
    std::cout << "Spatial map output: " << spatialName << std::endl;
    writeMap = true;
  }

  OptionType::Pointer timeseriesDisplacementOption = parser->GetOption("timeseries-displacement");
  if (timeseriesDisplacementOption && timeseriesDisplacementOption->GetNumberOfFunctions())
  {
    timeseriesDisplacementName = timeseriesDisplacementOption->GetFunction(0)->GetName();
    std::cout << "Time-series displacement map input 4d image: " << timeseriesDisplacementName << std::endl;
    timeseriesDisplacement = true;
  }
  OptionType::Pointer transformOption = parser->GetOption("transform");
  if (transformOption && transformOption->GetNumberOfFunctions())
  {
    transformIndex = std::stoi(transformOption->GetFunction(0)->GetName().c_str());
    std::cout << "Index of transform to output: " << transformIndex << std::endl;
    writeTransform = true;
  }

  OptionType::Pointer mocoOption = parser->GetOption("moco");
  if (mocoOption && mocoOption->GetNumberOfFunctions())
  {
    mocoName = mocoOption->GetFunction(0)->GetName();
    std::cout << "Moco file: " << mocoName << std::endl;
  }
  else
  {
    std::cerr << "Motion parameter file not specified" << std::endl;
    return EXIT_FAILURE;
  }

  OptionType::Pointer maskOption = parser->GetOption("mask");

  ImageType::Pointer mask = ImageType::New();
  if (maskOption && maskOption->GetNumberOfFunctions())
  {
    ReadImage<ImageType>(mask, maskOption->GetFunction(0)->GetName().c_str());
  }
  else
  {
    if (!writeTransform)
    {
      std::cerr << "Must use mask image" << std::endl;
      return EXIT_FAILURE;
    }
  }

  ImageType::Pointer map = ImageType::New();
  if (writeMap)
  {
    map->SetRegions(mask->GetLargestPossibleRegion());
    map->AllocateInitialized();

    map->SetOrigin(mask->GetOrigin());
    map->SetSpacing(mask->GetSpacing());
    map->SetDirection(mask->GetDirection());
  }


  TimeSeriesImageType::Pointer timeseriesImage = TimeSeriesImageType::New();
  TimeSeriesImageType::Pointer timeseriesDisplacementImage = TimeSeriesImageType::New();
  if (timeseriesDisplacement)
  {
    ReadImage<TimeSeriesImageType>(timeseriesImage, timeseriesDisplacementName.c_str());
    timeseriesDisplacementImage->SetRegions(timeseriesImage->GetLargestPossibleRegion());
    timeseriesDisplacementImage->SetOrigin(timeseriesImage->GetOrigin());
    timeseriesDisplacementImage->SetSpacing(timeseriesImage->GetSpacing());
    timeseriesDisplacementImage->SetDirection(timeseriesImage->GetDirection());
    timeseriesDisplacementImage->AllocateInitialized();
  }


  bool doFramewise = false;
  doFramewise = parser->Convert<bool>(parser->GetOption("framewise")->GetFunction()->GetName());
  std::cout << "Framewise = " << doFramewise << std::endl;

  MocoReaderType::Pointer mocoReader = MocoReaderType::New();
  mocoReader->SetFileName(mocoName.c_str());
  mocoReader->SetFieldDelimiterCharacter(',');
  mocoReader->HasColumnHeadersOn();
  mocoReader->HasRowHeadersOff();
  mocoReader->Update();

  MocoDataArrayType::Pointer mocoDataArray = mocoReader->GetOutput();
  std::cout << "Read motion correction data of size: " << mocoDataArray->GetMatrix().rows() << " x "
            << mocoDataArray->GetMatrix().cols() << std::endl;

  unsigned int nTransformParams = mocoDataArray->GetMatrix().cols() - 2;
  // std::cout << "# Transform parameters = " << nTransformParams << std::endl;

  WriterMatrixType dataMatrix(mocoDataArray->GetMatrix().rows(), 2);

  // Extract a single 3D transform and write to file
  if (writeTransform)
  {

    AffineTransformType::Pointer        affineTransform1 = AffineTransformType::New();
    AffineTransformType::ParametersType params1;
    params1.SetSize(nTransformParams);

    for (unsigned int t = 0; t < nTransformParams; t++)
    {
      params1[t] = mocoDataArray->GetMatrix()(transformIndex, t + 2);
    }

    // If rigid motion corretion
    if (nTransformParams == 6)
    {
      RigidTransformType::Pointer rigid1 = RigidTransformType::New();
      rigid1->SetParameters(params1);
      affineTransform1->SetMatrix(rigid1->GetMatrix());
      affineTransform1->SetTranslation(rigid1->GetTranslation());
    }
    else if (nTransformParams == 12)
    {
      affineTransform1->SetParameters(params1);
    }
    else
    {
      std::cout << "Unknown transform type! - Exiting" << std::endl;
      return EXIT_FAILURE;
    }

    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetInput(affineTransform1);
    transformWriter->SetFileName(outputName.c_str());
#if ITK_VERSION_MAJOR >= 5
    transformWriter->SetUseCompression(true);
#endif
    transformWriter->Update();

    return EXIT_SUCCESS;
  }

  for (unsigned int i = 0; i < mocoDataArray->GetMatrix().rows(); i++)
  {

    AffineTransformType::Pointer        affineTransform1 = AffineTransformType::New();
    AffineTransformType::Pointer        affineTransform2 = AffineTransformType::New();
    AffineTransformType::ParametersType params1;
    AffineTransformType::ParametersType params2;
    params1.SetSize(nTransformParams);
    params2.SetSize(nTransformParams);

    for (unsigned int t = 0; t < nTransformParams; t++)
    {
      params1[t] = mocoDataArray->GetMatrix()(i, t + 2);
      if (i < (mocoDataArray->GetMatrix().rows() - 1))
      {
        params2[t] = mocoDataArray->GetMatrix()(i + 1, t + 2);
      }
    }

    // If rigid motion correction
    if (nTransformParams == 6)
    {
      RigidTransformType::Pointer rigid1 = RigidTransformType::New();
      RigidTransformType::Pointer rigid2 = RigidTransformType::New();
      rigid1->SetParameters(params1);
      rigid2->SetParameters(params2);

      affineTransform1->SetMatrix(rigid1->GetMatrix());
      affineTransform1->SetTranslation(rigid1->GetTranslation());

      affineTransform2->SetMatrix(rigid2->GetMatrix());
      affineTransform2->SetTranslation(rigid2->GetTranslation());
    }
    else if (nTransformParams == 12)
    {
      affineTransform1->SetParameters(params1);
      affineTransform2->SetParameters(params2);
    }
    else
    {
      std::cout << "Unknown transform type! - Exiting" << std::endl;
      return EXIT_FAILURE;
    }


    double meanDisplacement = 0.0;
    double maxDisplacement = 0.0;
    double count = 0;

    // we iterate over the time-series one time volume at a time
    TimeSeriesRegionType timeSeriesRegion;
    TimeSeriesIndexType  timeSeriesIndex;
    TimeSeriesSizeType   timeSeriesSize;
    if (timeseriesDisplacement)
    {
      for (int jj = 0; jj < 3; jj++)
      {
        timeSeriesIndex[jj] = mask->GetLargestPossibleRegion().GetIndex()[jj];
        timeSeriesSize[jj] = mask->GetLargestPossibleRegion().GetSize()[jj];
      }
      timeSeriesIndex[3] = i;
      timeSeriesSize[3] = 1;
      timeSeriesRegion.SetIndex(timeSeriesIndex);
      timeSeriesRegion.SetSize(timeSeriesSize);
    }
    TimeSeriesIteratorType timeseriesIterator(timeseriesDisplacementImage, timeSeriesRegion);
    IteratorType           it(mask, mask->GetLargestPossibleRegion());
    while (!it.IsAtEnd())
    {

      if (it.Value() > 0)
      {

        ImageType::IndexType idx = it.GetIndex();
        ImageType::PointType pt;
        mask->TransformIndexToPhysicalPoint(idx, pt);

        ImageType::PointType pt1 = affineTransform1->TransformPoint(pt);

        double dist = 0;
        if (doFramewise && (i < (mocoDataArray->GetMatrix().rows() - 1)))
        {
          ImageType::PointType pt2 = affineTransform2->TransformPoint(pt);
          dist = pt1.EuclideanDistanceTo(pt2);
          if (writeMap)
          {
            map->SetPixel(idx, map->GetPixel(idx) + dist / (mocoDataArray->GetMatrix().rows() - 1));
          }
        }
        else
        {
          dist = pt.EuclideanDistanceTo(pt1);
          if (writeMap)
          {
            map->SetPixel(idx, map->GetPixel(idx) + dist / mocoDataArray->GetMatrix().rows());
          }
        }

        if (doFramewise && (i == mocoDataArray->GetMatrix().rows() - 1))
        {
          dist = 0.0;
        }


        if (dist > maxDisplacement)
        {
          maxDisplacement = dist;
        }
        meanDisplacement += dist;
        ++count;
        if (timeseriesDisplacement)
        {
          timeseriesIterator.Set(dist);
        }
      }
      ++it;
      if (timeseriesDisplacement)
      {
        ++timeseriesIterator;
      }
    }

    meanDisplacement /= count;
    dataMatrix(i, 0) = meanDisplacement;
    dataMatrix(i, 1) = maxDisplacement;
    // std::cout << i << "," << maxDisplacement << "," << meanDisplacement << std::endl;
  }

  // Write summary stats to output file
  WriterType::Pointer writer = WriterType::New();
  writer->ColumnHeadersPushBack("Mean");
  writer->ColumnHeadersPushBack("Max");
  writer->SetInput(&dataMatrix);
  writer->SetFileName(outputOption->GetFunction(0)->GetName().c_str());
  writer->Write();

  if (writeMap)
  {
    ANTs::WriteImage<ImageType>(map, spatialOption->GetFunction(0)->GetName().c_str());
  }

  if (timeseriesDisplacement)
  {
    std::string tsOutImageName = outputName;
    std::size_t lastdot = outputName.find_last_of(".");
    if (lastdot != std::string::npos)
      tsOutImageName = outputName.substr(0, lastdot);
    tsOutImageName += ".nii.gz";
    ANTs::WriteImage<TimeSeriesImageType>(timeseriesDisplacementImage, tsOutImageName.c_str());
  }

  return EXIT_SUCCESS;
}

void
antsMotionCorrStatsInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  using OptionType = itk::ants::CommandLineParser::OptionType;

  {
    std::string description = std::string("Mask image - compute displacements within mask.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask");
    option->SetShortName('x');
    option->SetDescription(description);
    option->SetUsageOption(0, "mask.nii.gz");
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("motion correction parameters from antsMotionCorr.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("moco");
    option->SetShortName('m');
    option->SetDescription(description);
    option->SetUsageOption(0, "MOCOparams.csv");
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Specify the output file");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetDescription(description);
    option->SetUsageOption(0, "corrected.csv");
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Specify the index for a 3D transform to output");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transform");
    option->SetShortName('t');
    option->SetDescription(description);
    // option->SetUsageOption( 0, "0" );
    parser->AddOption(option);
  }

  {
    std::string description = std::string("do framewise summarywise stats");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("framewise");
    option->SetShortName('f');
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("output image of displacement magnitude");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("spatial-map");
    option->SetShortName('s');
    option->SetDescription(description);
    // option->AddFunction( std::string( "0" ) );
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("output 4d time-series image of displacement magnitude");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("timeseries-displacement");
    option->SetShortName('d');
    option->SetDescription(description);
    parser->AddOption(option);
  }
  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsMotionCorrStats(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsMotionCorrStats");

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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand(argv[0]);

  std::string commandDescription = std::string(
    "antsMotionCorrStats - create summary measures of the parameters that are output by antsMotionCorr. Currently only "
    "works for linear transforms. Outputs the mean and max displacements for the voxels within a provided mask, at "
    "each time point. By default the displacements are relative to the reference space, but the framewise option may "
    "be used to provide displacements between consecutive time points");
  parser->SetCommandDescription(commandDescription);
  antsMotionCorrStatsInitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc < 2 || parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, false);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  std::cout << std::endl << "Running " << argv[0] << std::endl << std::endl;

  return ants_motion_stats(parser);
}

} // namespace ants
