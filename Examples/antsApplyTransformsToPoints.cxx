#include "itkCSVNumericObjectFileWriter.h"
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkantsRegistrationHelper.h"
#include "itkCSVArray2DFileReader.h"
#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
template <unsigned int Dimension, typename RealType>
int
antsApplyTransformsToPoints(itk::ants::CommandLineParser::Pointer & parser)
{
  using MatrixType = vnl_matrix<RealType>;
  MatrixType points_out;
  MatrixType points_in;
  using ImageType = itk::Image<RealType, 2>;
  using ReaderType = itk::CSVArray2DFileReader<RealType>;
  using DataFrameObjectType = itk::CSVArray2DDataObject<RealType>;
  using StringVectorType = typename DataFrameObjectType::StringVectorType;
  StringVectorType            colheadernames;
  typename ImageType::Pointer pointimage = nullptr;

  itk::ants::CommandLineParser::OptionType::Pointer antsrOption = parser->GetOption("forantsr");
  unsigned int                                      forANTsR = 0;
  if (antsrOption && antsrOption->GetNumberOfFunctions() > 0)
  {
    forANTsR = parser->Convert<unsigned int>(antsrOption->GetFunction(0)->GetName());
  }

  /**
   * Input object option
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption("input");
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");
  if (inputOption && inputOption->GetNumberOfFunctions() > 0)
  {
    std::size_t lengthInputFileName = std::strlen(inputOption->GetFunction(0)->GetName().c_str());
    std::string ext = (inputOption->GetFunction(0)->GetName()).substr(lengthInputFileName - 4);

    if (strcmp(ext.c_str(), ".csv") == 0)
    {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName((inputOption->GetFunction(0)->GetName()).c_str());
      reader->SetFieldDelimiterCharacter(',');
      reader->SetStringDelimiterCharacter('"');
      reader->HasColumnHeadersOn();
      reader->HasRowHeadersOff();
      //    reader->UseStringDelimiterCharacterOff();
      try
      {
        reader->Update();
      }
      catch (const itk::ExceptionObject & exp)
      {
        std::cerr << "Exception caught!" << std::endl;
        std::cerr << exp << std::endl;
      }
      typename DataFrameObjectType::Pointer dfo = reader->GetOutput();
      colheadernames = dfo->GetColumnHeaders();
      if (colheadernames.size() < Dimension)
      {
        std::cerr
          << "Input csv file must have column names such as x,y,z,t,label - where there are a minimum of "
             "N-Spatial-Dimensions names e.g. x,y in 2D.  ***Or pass in a 2D mha (meta format) binary image file."
          << std::endl;
        return EXIT_FAILURE;
      }
      points_in = dfo->GetMatrix();
      points_out.set_size(points_in.rows(), points_in.cols());
    }
    else if (strcmp(ext.c_str(), ".mha") == 0 || forANTsR)
    {
      std::string fn1 = inputOption->GetFunction(0)->GetName();
      ReadImage<ImageType>(pointimage, fn1.c_str());
      typename ImageType::IndexType ind;
      ind.Fill(0);
      typename ImageType::SizeType sz;
      sz.Fill(0);
      sz = pointimage->GetLargestPossibleRegion().GetSize();
      points_in.set_size(sz[0], sz[1]);
      points_out.set_size(points_in.rows(), points_in.cols());
      for (unsigned int d = 0; d < sz[0]; d++)
      {
        for (unsigned int dd = 0; dd < sz[1]; dd++)
        {
          ind[0] = d;
          ind[1] = dd;
          points_in(d, dd) = pointimage->GetPixel(ind);
        }
      }
    }
    else
    {
      std::cerr << "An input csv or mha file is required." << std::endl;
      return EXIT_FAILURE;
    }

    if (points_in.cols() < Dimension)
    {
      std::cerr << "The number of columns in the input point set is fewer than " << Dimension << " Exiting."
                << std::endl;
      return EXIT_FAILURE;
    }

    if (outputOption && outputOption->GetNumberOfFunctions() > 0)
    {
      if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) == 0)
      {
        std::cerr << "An input csv file is required." << std::endl;
        return EXIT_FAILURE;
      }
    }

    /**
     * Transform option
     */
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    using MatrixOffsetTransformType = itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension>;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
    using MatrixOffsetTransformType = itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension>;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

    /**
     * Load an identity transform in case no transforms are loaded.
     */
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    using AffineTransformType = itk::AffineTransform<RealType, Dimension>;
    typename AffineTransformType::Pointer aff = AffineTransformType::New();
    aff->SetIdentity();

    using CompositeTransformType = itk::CompositeTransform<RealType, Dimension>;
    typename CompositeTransformType::InputPointType            point_in;
    typename CompositeTransformType::OutputPointType           point_out;
    typename itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption("transform");

    std::vector<bool>                        isDerivedTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<RealType, Dimension>(parser, transformOption, isDerivedTransform, forANTsR);

    if (compositeTransform->GetNumberOfTransforms() == 0)
      compositeTransform->AddTransform(aff);

    if (compositeTransform.IsNull())
    {
      return EXIT_FAILURE;
    }
    for (unsigned int pointct = 0; pointct < points_in.rows(); pointct++)
    {
      point_in.Fill(0);
      point_out.Fill(0);
      for (unsigned int p = 0; p < Dimension; p++)
      {
        point_in[p] = points_in(pointct, p);
      }
      point_out = compositeTransform->TransformPoint(point_in);
      for (unsigned int p = 0; p < Dimension; p++)
      {
        points_out(pointct, p) = point_out[p];
      }
      for (unsigned int p = Dimension; p < points_in.cols(); p++)
      {
        points_out(pointct, p) = points_in(pointct, p);
      }
    }
    /**
     * output
     */
    if (outputOption && outputOption->GetNumberOfFunctions() > 0)
    {
      std::string outputFileName = "";
      if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>(outputOption->GetFunction(0)->GetParameter(1)) == 0)
      {
        outputFileName = outputOption->GetFunction(0)->GetParameter(0);
      }
      else
      {
        outputFileName = outputOption->GetFunction(0)->GetName();
      }
      std::size_t lengthOutputFileName = std::strlen(outputFileName.c_str());
      std::string exto = outputFileName.substr(lengthOutputFileName - 4);

      if (strcmp(exto.c_str(), ".csv") == 0)
      {
        StringVectorType ColumnHeaders = colheadernames;
        using WriterType = itk::CSVNumericObjectFileWriter<RealType, 1, 1>;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(outputFileName);
        writer->SetInput(&points_out);
        writer->SetColumnHeaders(ColumnHeaders);
        try
        {
          writer->Write();
        }
        catch (const itk::ExceptionObject & exp)
        {
          std::cerr << "Exception caught!" << std::endl;
          std::cerr << exp << std::endl;
          return EXIT_FAILURE;
        }
      }
      if ((strcmp(exto.c_str(), ".mha") == 0 || forANTsR) && (!pointimage.IsNull()))
      {
        typename ImageType::IndexType ind;
        ind.Fill(0);
        typename ImageType::SizeType sz;
        sz.Fill(0);
        sz = pointimage->GetLargestPossibleRegion().GetSize();
        if (sz[0] != points_out.rows() || sz[1] != points_out.cols())
        {
          std::cout << " the size of points_out must match the input pointimage" << std::endl;
          return EXIT_FAILURE;
        }
        for (unsigned int d = 0; d < sz[0]; d++)
          for (unsigned int dd = 0; dd < sz[1]; dd++)
          {
            ind[0] = d;
            ind[1] = dd;
            pointimage->SetPixel(ind, points_out(d, dd));
          }
        ANTs::WriteImage<ImageType>(pointimage, outputFileName.c_str());
      }
    }
  }

  return EXIT_SUCCESS;
}

static void
antsApplyTransformsToPointsInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{

  {
    std::string description =
      std::string("This option forces the points to be treated as a specified-") + std::string("dimensionality.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("use-double-precision");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("precision");
    option->SetShortName('p');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("set true for ANTsR IO");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("forantsr");
    option->SetShortName('f');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Currently, the only input supported is a csv file with ") +
      std::string("columns including x,y,z,t (all 4) column headers. ") +
      std::string("if you dont have 4D data, still supply 4D filling in extra places with zero. ") +
      std::string("The points should be defined in physical space. ") +
      std::string("Points are transformed in the OPPOSITE direction of images, therefore ") +
      std::string("you should pass the inverse of what is needed to warp the images. ") +
      std::string("Eg if the image is warped by  Affine.mat, you should pass the inverse of Affine.mat ") +
      std::string("to transform points defined in the same space as the image. ") +
      std::string("If in doubt how to convert coordinates from your files to the space ") +
      std::string("required by antsApplyTransformsToPoints try creating/drawing a simple ") +
      std::string("label volume with only one voxel set to 1 and all others set to 0. ") +
      std::string("Write down the voxel coordinates. Then use ImageMaths LabelStats to find ") +
      std::string("out what coordinates for this voxel antsApplyTransformsToPoints is ") +
      std::string("expecting.  ITK uses a LPS coordinate system.  See "
                  "http://sourceforge.net/p/advants/discussion/840261/thread/2a1e9307/") +
      std::string(" ***Or pass in a 2D mha (meta format) binary image file.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("One can output the warped points to a csv file.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "warpedOutputFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Several transform options are supported including all ") +
                              std::string("those defined in the ITK library in addition to ") +
                              std::string("a deformation field transform.  The ordering of ") +
                              std::string("the transformations follows the ordering specified ") +
                              std::string("on the command line.  An identity transform is pushed ") +
                              std::string("onto the transformation stack. Each new transform ") +
                              std::string("encountered on the command line is also pushed onto ") +
                              std::string("the transformation stack. Then, to warp the input object, ") +
                              std::string("each point comprising the input object is warped first ") +
                              std::string("according to the last transform pushed onto the stack ") +
                              std::string("followed by the second to last transform, etc. until ") +
                              std::string("the last transform encountered which is the identity ") +
                              std::string("transform. ") +
                              std::string("Also, it should be noted that the inverse transform can ") +
                              std::string("be accommodated with the usual caveat that such an inverse ") +
                              std::string("must be defined by the specified transform class ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transform");
    option->SetShortName('t');
    option->SetUsageOption(0, "transformFileName");
    option->SetUsageOption(1, "[transformFileName,useInverse]");
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
antsApplyTransformsToPoints(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsApplyTransformsToPoints");
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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand(argv[0]);

  std::string examplestring =
    std::string("reads in a csv file with the first D columns defining the spatial location where the spatial location "
                "is defined in physical coordinates.    the csv file should have a header row.   here is an example") +
    std::string("\n") + std::string("cat chicken-3.csv ") + std::string("x,y,z,t,label,comment") + std::string("\n") +
    std::string("82.5,116.5,0,0,1,this is the breast") + std::string("\n") +
    std::string("137.5,35.5,0,0,2,this is the beak") + std::string("\n") +
    std::string("antsApplyTransformsToPoints -d 2 -i chicken-3.csv -o test.csv -t [chicken3to4.mat ,1 ]") +
    std::string("\n") + std::string("cat test.csv ") + std::string("\n") + std::string("x,y,z,t,label,comment") +
    std::string("\n") + std::string("10.8945447481644,162.082675013049,0,0,1,nan") + std::string("\n") +
    std::string("7.5367085472988,52.099713111629,0,0,2,nan") + std::string("\n") +
    std::string("the nan appears in the last column until the ITK CSV I/O can handle mixed numeric / string types.  if "
                "your input is fully numeric, all is well.");

  std::string mhastring =
    std::string("\n\n**** We now can also read / write .mha files.") + std::string("\n") +
    std::string("This is a simple binary format (Meta format - look it up!) that is much faster to read/write than csv "
                "format.\n Note: To write a mha file, you must also pass an mha file as input.\n");

  std::string commandDescription =
    std::string("antsApplyTransformsToPoints, applied to an input image, transforms it ") +
    std::string("according to a reference image and a transform ") + std::string("(or a set of transforms).  ") +
    examplestring + mhastring;

  parser->SetCommandDescription(commandDescription);
  antsApplyTransformsToPointsInitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc < 2 ||
      (parser->GetOption("help") && (parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))))
  {
    parser->PrintMenu(std::cout, 5, false);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (parser->GetOption('h') && (parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName())))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  unsigned int                                      dimension = 3;
  itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption("dimensionality");
  if (dimOption && dimOption->GetNumberOfFunctions() > 0)
  {
    dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
  }
  else
  {
    std::cerr << "No -d ( dimensionality ) option is specified.  Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  itk::ants::CommandLineParser::OptionType::Pointer precOption = parser->GetOption("precision");
  unsigned int                                      myprecision = 0;
  if (precOption && precOption->GetNumberOfFunctions() > 0)
  {
    myprecision = parser->Convert<unsigned int>(precOption->GetFunction(0)->GetName());
  }

  if (myprecision == 1)
  {
    switch (dimension)
    {
      case 2:
      {
        return antsApplyTransformsToPoints<2, double>(parser);
      }
      break;
      case 3:
      {
        return antsApplyTransformsToPoints<3, double>(parser);
      }
      break;
      case 4:
      {
        return antsApplyTransformsToPoints<4, double>(parser);
      }
      break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
  }
  else
  {
    switch (dimension)
    {
      case 2:
      {
        return antsApplyTransformsToPoints<2, float>(parser);
      }
      break;
      case 3:
      {
        return antsApplyTransformsToPoints<3, float>(parser);
      }
      break;
      case 4:
      {
        return antsApplyTransformsToPoints<4, float>(parser);
      }
      break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

} // namespace ants
