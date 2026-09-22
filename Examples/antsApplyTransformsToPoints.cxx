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
    std::string         description = std::string("Use double-precision. Default = 0 (single precision).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("precision");
    option->SetShortName('p');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Set true for ANTsR IO. Default = 0.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("forantsr");
    option->SetShortName('f');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
        "Input can either be a CSV file or a 2D binary meta image (.mha). CSV input should have "
        "at least D columns where D is the spatial dimensionality of the transform. The first D "
        "columns should have column headers x,y,z,t. Additional numerical columns are passed "
        "through to the output file, but not transformed. "
        "MHA input should be 2D with the first dimension being the point index and the second "
        "dimension being the point coordinates. "
        " "
        "The input points should be defined in LPS+ physical space as defined by ITK. ";

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Output file name. Output format is the same as the input format.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "warpedOutputFileName");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "An ANTs transforms to apply. Use multiple times to chain transforms. Use [transformFile,1] "
      "to apply the inverse, for transforms that define an explicit inverse (eg affine transforms)."
      " "
      "Note on transform direction: The required 'forward' or 'inverse' warps for points "
      "are the OPPOSITE of those used to resample images. For warps from antsRegistration with a given "
      "'fixed' and 'moving' image: to warp a surface defined in the moving-image space into the "
      "fixed-image space, use the same transforms you would use with antsApplyTransforms to warp "
      "the fixed image into moving space. See "
      "https://github.com/ANTsX/ANTs/wiki/Applying-transforms-to-point-data";

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

  std::string commandDescription =
    "antsApplyTransformsToPoints transforms points in ITK LPS+ physical space by applying a set of "
    "transforms supplied on the command line. "

    "Input and output can be text or binary. Output format is the same as the input."
    "Text format is a csv file with at least D columns where D is the spatial dimensionality of the transform. "
    "Additional numeric columns (eg, labels) are passed through to the output file, but not transformed. "
    "Non-numeric data are not supported. The first D columns should be physical coordinate in ITK LPS+ space. "
    "Binary format is a 2D meta image (.mha) with the first dimension being the point index and the second "
    "dimension being the point coordinates, again in ITK LPS+ space, in units of mm. "

    "Note on transforms: The required 'forward' or 'inverse' warps for points are the OPPOSITE of those used "
    "to resample images in the same direction. ";

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
