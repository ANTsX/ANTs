/*
 * antsApplyTransformsToTRX
 *
 * Applies a series of ANTs transforms (warps, affines, composite, etc.) to the
 * streamline coordinates of a TRX tractography file (.trx), producing a warped
 * output in TRX format.
 *
 * This tool is the TRX equivalent of antsApplyTransformsToPoints
 *
 * Coordinate system handling
 * --------------------------
 * TRX files store streamline coordinates in RAS+ physical space, whereas ITK
 * (and all ANTs transforms) operate in LPS+ space.  This tool automatically
 * converts between the two conventions by negating the X and Y coordinates
 * immediately before applying the composite transform and again immediately
 * after, equivalent to the itk_lps=True flag used when converting to/from CSV.
 *
 *
 * Example
 * -------
 * Apply a SyN warp and affine to a tractogram (transforms applied last-first,
 * same ordering as antsApplyTransforms for images):
 *
 *   antsApplyTransformsToTRX \
 *     -i tractogram.trx \
 *     -o tractogram.warped.trx \
 *     -t Warp.nii.gz \
 *     -t Affine.mat
 *
 * Note on transform direction: streamline vertices are transformed in the
 * OPPOSITE direction to images.  To warp a tractogram defined in moving-image
 * space into fixed-image space, pass the same transforms you would use for
 * antsApplyTransforms on the moving image (not the inverses).
 *
 * See also: antsApplyTransforms, antsApplyTransformsToPoints,
 *           antsApplyTransformsToGifti
 */

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"

#include "itkTrxFileReader.h"
#include "itkTrxStreamWriter.h"
#include "itkTrxStreamlineData.h"

#include "vnl/vnl_matrix.h"

#include <string>
#include <vector>

namespace ants
{

template <typename RealType>
int
antsApplyTransformsToTRX(itk::ants::CommandLineParser::Pointer & parser)
{
  constexpr unsigned int Dimension = 3;

  using AffineTransformType = itk::AffineTransform<RealType, Dimension>;
  using CompositeTransformType = itk::CompositeTransform<RealType, Dimension>;

  // Register the matrix offset transform base class so that ANTs affine files
  // (.mat) are recognized by the transform factory.
  using MatrixOffsetTransformType = itk::MatrixOffsetTransformBase<RealType, Dimension, Dimension>;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

  // Identity fallback if no transforms are specified.
  typename AffineTransformType::Pointer identityAff = AffineTransformType::New();
  identityAff->SetIdentity();

  // -----------------------------------------------------------------------
  // Parse required options
  // -----------------------------------------------------------------------
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption("input");
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");

  if (!inputOption || inputOption->GetNumberOfFunctions() == 0)
  {
    std::cerr << "No input TRX file specified (use -i / --input)." << std::endl;
    return EXIT_FAILURE;
  }
  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    std::cerr << "No output TRX file specified (use -o / --output)." << std::endl;
    return EXIT_FAILURE;
  }

  std::string inputFile = inputOption->GetFunction(0)->GetName();
  std::string outputFile = outputOption->GetFunction(0)->GetName();

  // -----------------------------------------------------------------------
  // Build composite ITK transform from --transform options
  // -----------------------------------------------------------------------
  typename itk::ants::CommandLineParser::OptionType::Pointer transformOption =
    parser->GetOption("transform");

  std::vector<bool> isDerivedTransform;
  typename CompositeTransformType::Pointer compositeTransform =
    GetCompositeTransformFromParserOption<RealType, Dimension>(parser, transformOption, isDerivedTransform);

  if (compositeTransform.IsNull())
  {
    return EXIT_FAILURE;
  }
  if (compositeTransform->GetNumberOfTransforms() == 0)
  {
    compositeTransform->AddTransform(identityAff);
  }

  // -----------------------------------------------------------------------
  // Read TRX file
  // -----------------------------------------------------------------------
  itk::TrxFileReader::Pointer reader = itk::TrxFileReader::New();
  reader->SetFileName(inputFile);
  try
  {
    reader->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Failed to read TRX file: " << inputFile << "\n" << e << std::endl;
    return EXIT_FAILURE;
  }

  itk::TrxStreamlineData::Pointer trxData = reader->GetOutput();

  std::cout << "Read " << trxData->GetNumberOfStreamlines() << " streamlines ("
            << trxData->GetNumberOfVertices() << " vertices) from: " << inputFile << std::endl;

  // -----------------------------------------------------------------------
  // Stream-transform to output.
  //
  // TransformToWriterChunkedReuseVnlBuffer iterates the backing TRX handle
  // in streamline-sized chunks (lazy, no full load into RAM), converts
  // RAS->LPS internally per chunk, applies the composite transform, then
  // pushes to TrxStreamWriter which converts LPS->RAS for storage.
  // -----------------------------------------------------------------------
  itk::TrxStreamWriter::Pointer writer = itk::TrxStreamWriter::New();
  writer->SetFileName(outputFile);

  if (trxData->HasVoxelToRasMatrix())
  {
    writer->SetVoxelToRasMatrix(trxData->GetVoxelToRasMatrix());
  }
  if (trxData->HasVoxelToLpsMatrix())
  {
    writer->SetVoxelToLpsMatrix(trxData->GetVoxelToLpsMatrix());
  }
  if (trxData->HasDimensions())
  {
    writer->SetDimensions(trxData->GetDimensions());
  }

  vnl_matrix<double> buffer;
  try
  {
    trxData->TransformToWriterChunkedReuseVnlBuffer(compositeTransform.GetPointer(), writer.GetPointer(), buffer);
    writer->Finalize();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Failed to write TRX file: " << outputFile << "\n" << e << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


static void
antsApplyTransformsToTRXInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  {
    std::string description =
      "Use double-precision floating point for transform computation (0 = float,"
      " 1 = double).  Float is faster; double may improve accuracy for very large"
      " deformation fields.  Streamline coordinates are stored at their original"
      " precision in the output TRX regardless of this setting.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("precision");
    option->SetShortName('p');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "Input TRX tractography file (.trx).";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "tractogram.trx");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Output TRX tractography file (.trx).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "tractogram.warped.trx");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "One or more ANTs transforms to apply, specified in the same order as"
      " antsApplyTransformsToPoints.  Transforms are"
      " applied last-specified first.  Use [transformFile,1] to apply the inverse"
      " of a transform."
      "\n\n"
      "Note on transform direction: streamline vertices move in the OPPOSITE direction "
      "to images.  Given warps from antsRegistration with a given 'fixed' and 'moving' image: "
      "to warp a tractogram defined in the moving-image space into the fixed-image space "
      "use the same transforms you would use with antsApplyTransforms to warp the fixed image "
      "into moving space. See https://github.com/ANTsX/ANTs/wiki/Applying-transforms-to-point-data";
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


// entry point for the library; parameter 'args' is equivalent to 'argv' in
// (argc,argv) of commandline parameters to 'main()'
int
antsApplyTransformsToTRX(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // Put the arguments into standard (argc, argv) format expected by the parser.
  args.insert(args.begin(), "antsApplyTransformsToTRX");
  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;

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
        delete[] argv[i];
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
    "Apply ANTs transforms to the streamline coordinates of a TRX tractography file.\n\n"

    "TRX files store coordinates in RAS+ space.  The RAS<->LPS conversion is\n"
    "performed automatically before and after the ANTs transform is applied.\n"
    "\n"
    "Example:\n"
    "  antsApplyTransformsToTRX -i tractogram.trx -o tractogram.warped.trx \\\n"
    "    -t Warp.nii.gz -t Affine.mat";

  parser->SetCommandDescription(commandDescription);
  antsApplyTransformsToTRXInitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc < 2 ||
      (parser->GetOption("help") &&
       (parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))))
  {
    parser->PrintMenu(std::cout, 5, false);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (parser->GetOption('h') &&
           (parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName())))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  // Select floating-point precision for transform computation.
  itk::ants::CommandLineParser::OptionType::Pointer precOption = parser->GetOption("precision");
  unsigned int myprecision = 0;
  if (precOption && precOption->GetNumberOfFunctions() > 0)
  {
    myprecision = parser->Convert<unsigned int>(precOption->GetFunction(0)->GetName());
  }

  if (myprecision == 1)
  {
    return antsApplyTransformsToTRX<double>(parser);
  }
  else
  {
    return antsApplyTransformsToTRX<float>(parser);
  }
}

} // namespace ants
