/*=========================================================================

  Program:   ITK General

  Author: Jeffrey T. Duda (jtduda@seas.upenn.edu)
  Institution: PICSL

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include <cstdio>

#include "itkBinaryThresholdImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "ReadWriteData.h"

namespace ants
{
/* FlipScalarVolume
 * This program takes a volume and flips it along the
 * indicated axes
 */

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
StackSlices(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "StackSlices");

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

  // Pixel and Image typedefs
  using PixelType = float;

  using ImageSeriesType = itk::Image<PixelType, 4>;
  using ImageType = itk::Image<PixelType, 3>;
  using SliceType = itk::Image<PixelType, 2>;
  using LabelSliceType = itk::Image<unsigned int, 2>;

  using ReaderType = itk::ImageFileReader<ImageType>;
  using Reader4DType = itk::ImageFileReader<ImageSeriesType>;

  using ExtractFilterType = itk::ExtractImageFilter<ImageType, SliceType>;
  using ExtractFilterType2 = itk::ExtractImageFilter<ImageSeriesType, SliceType>;

  using SliceIt = itk::ImageRegionIteratorWithIndex<SliceType>;

  // Check for valid input parameters
  if (argc < 5)
  {
    std::cout << "Usage: " << argv[0] << " outputvolume x y z inputvolume(s)" << std::endl;
    std::cout << "  The specific slice is chosen by specifying the index for x, y, xor z." << std::endl;
    std::cout << R"(  For example, an "x y z" selection of "30 -1 -1" will stack slice 30 )" << std::endl;
    std::cout << "  along the first dimension.  Also note that input 4-D volumes are treated " << std::endl;
    std::cout << "  as a series of 3-D volumes." << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  char * stackName = argv[1];
  int    dimVars[3];
  dimVars[0] = std::stoi(argv[2]);
  dimVars[1] = std::stoi(argv[3]);
  dimVars[2] = std::stoi(argv[4]);

  int dim = -1;
  int slice = -1;
  // which dim to extract slice from
  for (unsigned int i = 0; i < 3; i++)
  {
    if (dimVars[i] > -1)
    {
      if ((dim > -1) || (slice > -1))
      {
        std::cout << "Can only choose slice from 1 dimension" << std::endl;
        return EXIT_FAILURE;
      }
      dim = i;
      slice = dimVars[i];
    }
  }

  SliceType::SizeType size;
  size.Fill(0);

  ImageType::RegionType region3D;

  ImageSeriesType::Pointer imageSeries = nullptr;

  unsigned long nSlices = 0;

  bool inputIsA4DImage = true;
  if (argc > 6)
  {
    inputIsA4DImage = false;
    std::cout << "  Input is a set of 3-D volumes." << std::endl;
  }
  else
  {
    std::cout << "  Input is a 4-D image." << std::endl;
  }

  if (!inputIsA4DImage) // input is a set of 3-D volumes
  {
    nSlices = argc - 5;

    // std::cout << nSlices << std::endl;

    ReaderType::Pointer firstReader = ReaderType::New();
    firstReader->SetFileName(argv[5]);
    firstReader->Update();
    if (dim == 0)
    {
      size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
      size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
    }
    if (dim == 1)
    {
      size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
      size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
    }
    if (dim == 2)
    {
      size[0] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
      size[1] = firstReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
    }

    region3D = firstReader->GetOutput()->GetLargestPossibleRegion();
    region3D.SetSize(dim, nSlices);
    region3D.SetIndex(dim, 0);
  }
  else
  {
    Reader4DType::Pointer reader4D = Reader4DType::New();
    reader4D->SetFileName(argv[5]);
    imageSeries = reader4D->GetOutput();
    imageSeries->Update();
    imageSeries->DisconnectPipeline();

    nSlices = imageSeries->GetLargestPossibleRegion().GetSize()[3];
    ImageSeriesType::IndexType index4D = imageSeries->GetLargestPossibleRegion().GetIndex();

    if (dim == 0)
    {
      size[0] = imageSeries->GetLargestPossibleRegion().GetSize()[1];
      size[1] = imageSeries->GetLargestPossibleRegion().GetSize()[2];
      region3D.SetSize(0, nSlices);
      region3D.SetSize(1, size[0]);
      region3D.SetSize(2, size[1]);
      region3D.SetIndex(0, 0);
      region3D.SetIndex(1, index4D[1]);
      region3D.SetIndex(2, index4D[2]);
    }
    if (dim == 1)
    {
      size[0] = imageSeries->GetLargestPossibleRegion().GetSize()[0];
      size[1] = imageSeries->GetLargestPossibleRegion().GetSize()[2];
      region3D.SetSize(0, size[0]);
      region3D.SetSize(1, nSlices);
      region3D.SetSize(2, size[1]);
      region3D.SetIndex(0, index4D[0]);
      region3D.SetIndex(1, 0);
      region3D.SetIndex(2, index4D[2]);
    }
    if (dim == 2)
    {
      size[0] = imageSeries->GetLargestPossibleRegion().GetSize()[0];
      size[1] = imageSeries->GetLargestPossibleRegion().GetSize()[1];
      region3D.SetSize(0, size[0]);
      region3D.SetSize(1, size[1]);
      region3D.SetSize(2, nSlices);
      region3D.SetIndex(0, index4D[0]);
      region3D.SetIndex(1, index4D[1]);
      region3D.SetIndex(2, 0);
    }
  }

  std::cout << "  Output region size = " << region3D.GetSize() << std::endl;

  ImageType::Pointer stack = AllocImage<ImageType>(region3D);

  // Start stacking the slices while normalizing by the mean at each slice.

  for (unsigned int i = 0; i < nSlices; i++)
  {
    SliceType::Pointer stackSlice = nullptr;

    if (!inputIsA4DImage)
    {
      std::cout << " Slice " << i << " :: " << std::string(argv[5 + i]) << std::endl;

      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(argv[5 + i]);
      reader->Update();

      ImageType::RegionType extractRegion = reader->GetOutput()->GetLargestPossibleRegion();
      extractRegion.SetSize(dim, 0);
      extractRegion.SetIndex(dim, slice);

      ExtractFilterType::Pointer extracter = ExtractFilterType::New();
      extracter->SetInput(reader->GetOutput());
      extracter->SetDirectionCollapseToIdentity();
      extracter->SetExtractionRegion(extractRegion);

      stackSlice = extracter->GetOutput();
      stackSlice->Update();
      stackSlice->DisconnectPipeline();
    }
    else
    {
      std::cout << " Slice " << i << " :: " << std::endl;

      ImageSeriesType::RegionType extractRegion = imageSeries->GetLargestPossibleRegion();
      extractRegion.SetSize(dim, 0);
      extractRegion.SetIndex(dim, slice);
      extractRegion.SetSize(3, 0);
      extractRegion.SetIndex(3, i);

      ExtractFilterType2::Pointer extracter2 = ExtractFilterType2::New();
      extracter2->SetInput(imageSeries);
      extracter2->SetDirectionCollapseToIdentity();
      extracter2->SetExtractionRegion(extractRegion);

      stackSlice = extracter2->GetOutput();
      stackSlice->Update();
      stackSlice->DisconnectPipeline();
    }

    using ThresholderType = itk::BinaryThresholdImageFilter<SliceType, LabelSliceType>;
    ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput(stackSlice);
    thresholder->SetInsideValue(0);
    thresholder->SetOutsideValue(1);
    thresholder->SetLowerThreshold(itk::NumericTraits<PixelType>::NonpositiveMin());
    thresholder->SetUpperThreshold(0);

    using StatsFilterType = itk::LabelStatisticsImageFilter<SliceType, LabelSliceType>;
    StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetInput(stackSlice);
    stats->SetLabelInput(thresholder->GetOutput());
    stats->Update();

    PixelType sliceMean = stats->GetMean(1);

    SliceIt It(stackSlice, stackSlice->GetLargestPossibleRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      PixelType value = It.Get();

      ImageType::IndexType index;
      index.Fill(0);

      index[dim] = i;

      if (dim == 0)
      {
        index[1] = It.GetIndex()[0];
        index[2] = It.GetIndex()[1];
      }
      if (dim == 1)
      {
        index[0] = It.GetIndex()[0];
        index[2] = It.GetIndex()[1];
      }
      if (dim == 2)
      {
        index[0] = It.GetIndex()[0];
        index[1] = It.GetIndex()[1];
      }
      stack->SetPixel(index, value / sliceMean);
    }
  }

  ANTs::WriteImage<ImageType>(stack, stackName);

  // Input parameters
  //   char * inputName  = argv[1];
  //   unsigned int flip_x = std::stoi( argv[2] );
  //   unsigned int flip_y = std::stoi( argv[3] );
  //   unsigned int flip_z = std::stoi( argv[4] );
  //   char * outputName = argv[5];

  //   ReaderType::Pointer reader = ReaderType::New();
  //   reader->SetFileName( inputName );
  //   reader->Update();

  //   // flip in desired directions for correct display
  //   FlipFilterType::Pointer flip = FlipFilterType::New();
  //   FlipFilterType::FlipAxesArrayType flipOver;
  //   flipOver[0] = flip_x;
  //   flipOver[1] = flip_y;
  //   flipOver[2] = flip_z;
  //   flip->SetInput( reader->GetOutput() );
  //   flip->SetFlipAxes( flipOver );
  //   flip->Update();

  //   // write output
  //   WriterType::Pointer writer = WriterType::New();
  //   writer->SetInput( flip->GetOutput() );
  //   writer->SetFileName( outputName );
  //   writer->Update();

  return EXIT_SUCCESS;
}
} // namespace ants
