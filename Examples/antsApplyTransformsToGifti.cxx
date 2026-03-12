/*
 * antsApplyTransformsToGifti
 *
 * Applies a series of ANTs transforms (warps, affines, composite, etc.) to the
 * vertex coordinates of a GIFTI surface file (.gii), producing a warped output
 * surface in GIFTI format.
 *
 * This tool is the GIFTI-native equivalent of antsApplyTransformsToPoints and
 * replaces the multi-step workflow previously used in niworkflows
 * (https://github.com/nipreps/niworkflows/blob/master/niworkflows/interfaces/surf.py):
 *
 *   NormalizeSurf -> GiftiToCSV(itk_lps=True)
 *     -> antsApplyTransformsToPoints -> CSVToGifti(itk_lps=True)
 *
 *
 * Coordinate system handling
 * --------------------------
 * GIFTI files store vertex coordinates in RAS+ physical space, whereas ITK (and
 * all ANTs transforms) operate in LPS+ space.  This tool automatically converts
 * between the two conventions by negating the X and Y coordinates immediately
 * before applying the composite transform and again immediately after,
 * equivalent to the itk_lps=True flag in the niworkflows GiftiToCSV and
 * CSVToGifti interfaces.
 *
 *
 * FreeSurfer C_RAS correction
 * ---------------------------
 * When FreeSurfer converts a surface to GIFTI using mris_convert, vertex
 * coordinates are stored in "native surface RAS" space, which is offset from
 * world RAS by a translation vector called C_RAS.  This offset is embedded in
 * each NIFTI_INTENT_POINTSET data array's metadata under three keys:
 *
 *   VolGeomC_R   (C_RAS offset in the R direction)
 *   VolGeomC_A   (C_RAS offset in the A direction)
 *   VolGeomC_S   (C_RAS offset in the S direction)
 *
 * Before applying any ITK transform, this tool detects and adds the C_RAS
 * offset to bring coordinates into world RAS space, matching the behavior of
 * niworkflows' NormalizeSurf interface (see normalize_surfs() in surf.py)
 * After the transform is applied, the VolGeomC_R/A/S metadata
 * fields are zeroed out in the output file so that downstream tools do not
 * apply the offset a second time.
 *
 * Header fields modified in the output
 * ------------------------------------
 *   - NIFTI_INTENT_POINTSET vertex coordinates (the float32 data buffer) are
 *     updated in place; face topology and all other data arrays are preserved.
 *   - VolGeomC_R, VolGeomC_A, VolGeomC_S are set to "0.000000" if they were
 *     non-zero (C_RAS offset has been baked into the coordinates).
 *   - All other metadata, including the coordsys/xform matrix, label tables,
 *     and per-file metadata, are written through unchanged.
 *
 *
 * Example
 * -------
 * Apply a SyN warp and affine to a cortical surface (transforms applied last
 * first, same ordering as antsApplyTransforms for images):
 *
 *   antsApplyTransformsToGifti \
 *     -i lh.pial.surf.gii \
 *     -o lh.pial.warped.surf.gii \
 *     -t [Warp.nii.gz,0] \
 *     -t [Affine.mat,0]
 *
 * Note on transform direction: surface vertices are transformed in the OPPOSITE
 * direction to images.  To warp a surface defined in moving-image space into
 * fixed-image space, pass the same transforms you would use for
 * antsApplyTransforms on the moving image (not the inverses).
 *
 * See also: antsApplyTransforms, antsApplyTransformsToPoints
 */

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"

// gifticlib is a private dependency of ITK's MeshGifti module; its headers are
// installed alongside the other ITK headers during the ANTs superbuild.
#include "gifti_io.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace ants
{

template <typename RealType>
int
antsApplyTransformsToGifti(itk::ants::CommandLineParser::Pointer & parser)
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
    std::cerr << "No input GIFTI file specified (use -i / --input)." << std::endl;
    return EXIT_FAILURE;
  }
  if (!outputOption || outputOption->GetNumberOfFunctions() == 0)
  {
    std::cerr << "No output GIFTI file specified (use -o / --output)." << std::endl;
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
  // Read GIFTI file with all data.
  // Silence gifticlib's validator: HCP- and fMRIPrep-produced files commonly
  // carry a coordsys on the TRIANGLE array, which is valid in practice even
  // though the GIFTI spec reserves coordsys for POINTSET arrays only.
  // -----------------------------------------------------------------------
  const int savedVerb = gifti_get_verb();
  gifti_set_verb(0);
  gifti_image * gim = gifti_read_image(inputFile.c_str(), /*read_data=*/1);
  gifti_set_verb(savedVerb);
  if (!gim)
  {
    std::cerr << "Failed to read GIFTI file: " << inputFile << std::endl;
    return EXIT_FAILURE;
  }

  // -----------------------------------------------------------------------
  // Locate the NIFTI_INTENT_POINTSET data array (vertex coordinates)
  // -----------------------------------------------------------------------
  int pointsetIdx = -1;
  for (int i = 0; i < gim->numDA; ++i)
  {
    if (gim->darray[i]->intent == NIFTI_INTENT_POINTSET)
    {
      pointsetIdx = i;
      break;
    }
  }
  if (pointsetIdx < 0)
  {
    std::cerr << "No NIFTI_INTENT_POINTSET data array found in: " << inputFile << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }

  giiDataArray * da = gim->darray[pointsetIdx];

  if (da->datatype != NIFTI_TYPE_FLOAT32)
  {
    std::cerr << "POINTSET data array must use NIFTI_TYPE_FLOAT32. "
              << "Got datatype " << da->datatype << " in: " << inputFile << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }
  if (da->num_dim < 2 || da->dims[1] != 3)
  {
    std::cerr << "POINTSET data array must be Nx3. "
              << "Got dimensions " << da->dims[0] << "x" << da->dims[1]
              << " in: " << inputFile << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }

  const int    nVerts = da->dims[0];
  float *      coords = static_cast<float *>(da->data);

  // -----------------------------------------------------------------------
  // Detect and read the FreeSurfer C_RAS offset from data array metadata.
  //
  // When mris_convert produces a GIFTI, vertex coordinates are in FreeSurfer's
  // "native surface RAS" space.  World RAS = native surface RAS + C_RAS.
  // The offset is stored in the POINTSET data array's metadata under the keys
  // VolGeomC_R, VolGeomC_A, VolGeomC_S.
  //
  // Reference: niworkflows normalize_surfs() in interfaces/surf.py
  // (https://github.com/nipreps/niworkflows)
  // -----------------------------------------------------------------------
  const char * cRasKeys[3] = { "VolGeomC_R", "VolGeomC_A", "VolGeomC_S" };
  double       cRas[3] = { 0.0, 0.0, 0.0 };
  bool         hasCRas = false;

  for (int d = 0; d < 3; ++d)
  {
    char * val = gifti_get_meta_value(&da->meta, cRasKeys[d]);
    if (val)
    {
      cRas[d] = std::atof(val);
      if (cRas[d] != 0.0)
      {
        hasCRas = true;
      }
    }
  }

  if (hasCRas)
  {
    std::cout << "Detected FreeSurfer C_RAS offset: "
              << "R=" << cRas[0] << "  A=" << cRas[1] << "  S=" << cRas[2]
              << "  -- applying before transform." << std::endl;
  }

  // -----------------------------------------------------------------------
  // Apply transforms to each vertex
  // -----------------------------------------------------------------------
  for (int v = 0; v < nVerts; ++v)
  {
    float * pt = coords + v * 3;

    // 1. Apply the FreeSurfer C_RAS offset: native surface RAS -> world RAS.
    //    This matches the addition performed in niworkflows normalize_surfs()
    //    before applying any further transforms.
    double x = static_cast<double>(pt[0]) + cRas[0];
    double y = static_cast<double>(pt[1]) + cRas[1];
    double z = static_cast<double>(pt[2]) + cRas[2];

    // 2. RAS+ -> LPS+: negate X and Y to convert from the GIFTI/neuroimaging
    //    convention to ITK's internal coordinate convention.
    typename CompositeTransformType::InputPointType itkPt;
    itkPt[0] = static_cast<RealType>(-x);
    itkPt[1] = static_cast<RealType>(-y);
    itkPt[2] = static_cast<RealType>(z);

    // 3. Apply the ANTs composite transform in LPS+ space.
    typename CompositeTransformType::OutputPointType itkPtOut =
      compositeTransform->TransformPoint(itkPt);

    // 4. LPS+ -> RAS+: negate X and Y back before storing.
    pt[0] = static_cast<float>(-itkPtOut[0]);
    pt[1] = static_cast<float>(-itkPtOut[1]);
    pt[2] = static_cast<float>(itkPtOut[2]);
  }

  // -----------------------------------------------------------------------
  // Zero out C_RAS metadata in the output so downstream tools do not apply
  // the offset a second time.
  // -----------------------------------------------------------------------
  if (hasCRas)
  {
    int dalist[1] = { pointsetIdx };
    for (int d = 0; d < 3; ++d)
    {
      gifti_set_DA_meta(gim, cRasKeys[d], "0.000000", dalist, 1, /*replace=*/1);
    }
  }

  // -----------------------------------------------------------------------
  // Write output GIFTI.  gifti_write_image preserves all data arrays (face
  // topology, shape data, etc.), the coordsys/xform matrix, label tables, and
  // all file- and array-level metadata.  Only the POINTSET coordinate buffer
  // and the zeroed-out VolGeomC_* fields differ from the input.
  //
  // -----------------------------------------------------------------------
  if (gifti_write_image(gim, outputFile.c_str(), /*write_data=*/1) != 0)
  {
    std::cerr << "Failed to write output GIFTI: " << outputFile << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }

  gifti_free_image(gim);
  return EXIT_SUCCESS;
}


static void
antsApplyTransformsToGiftiInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  {
    std::string description =
      "Use double-precision floating point for transform computation (0 = float,"
      " 1 = double).  Float is faster; double may improve accuracy for very large"
      " deformation fields.  Vertex coordinates are always stored as float32 in the"
      " output GIFTI regardless of this setting.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("precision");
    option->SetShortName('p');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "Input GIFTI surface file (.gii).  Both standard GIFTI files (e.g., produced"
      " by HCP pipelines or Connectome Workbench) and FreeSurfer-generated GIFTI"
      " files (produced by mris_convert) are supported.  Native FreeSurfer binary"
      " surfaces (.pial, .white, .inflated, etc.) must first be converted to GIFTI"
      " using mris_convert -- see the tool description for the exact commands.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputSurface.gii");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Output GIFTI surface file (.gii).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "outputSurface.gii");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "One or more ANTs transforms to apply, specified in the same order as"
      " antsApplyTransforms and antsApplyTransformsToPoints.  Transforms are"
      " applied last-specified first.  Use [transformFile,1] to apply the inverse"
      " of a transform."
      "\n\n"
      " Note on transform direction: surface vertices move in the OPPOSITE direction"
      " to images.  To warp a surface defined in moving-image space into fixed-image"
      " space, pass the same transforms you would supply to antsApplyTransforms for"
      " the moving image (not their inverses).";
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
antsApplyTransformsToGifti(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // Put the arguments into standard (argc, argv) format expected by the parser.
  args.insert(args.begin(), "antsApplyTransformsToGifti");
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
    "Apply ANTs transforms to the vertex coordinates of a GIFTI surface file.\n\n"

    "Preparing FreeSurfer surfaces:\n"
    "  mris_convert $SUBJECTS_DIR/$SUBJECT/surf/lh.pial  lh.pial.surf.gii\n"
    "  mris_convert $SUBJECTS_DIR/$SUBJECT/surf/rh.pial  rh.pial.surf.gii\n"
    "  mris_convert $SUBJECTS_DIR/$SUBJECT/surf/lh.white lh.white.surf.gii\n"
    "  mris_convert $SUBJECTS_DIR/$SUBJECT/surf/rh.white rh.white.surf.gii\n"
    "\n"
    "The FreeSurfer C_RAS offset (VolGeomC_R/A/S) is applied automatically\n"
    "and zeroed out in the output.  No separate normalization step is needed.\n"
    "\n"
    "Example:\n"
    "  antsApplyTransformsToGifti -i lh.pial.surf.gii -o lh.pial.warped.surf.gii \\\n"
    "    -t Warp.nii.gz -t Affine.mat";

  parser->SetCommandDescription(commandDescription);
  antsApplyTransformsToGiftiInitializeCommandLineOptions(parser);

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
    return antsApplyTransformsToGifti<double>(parser);
  }
  else
  {
    return antsApplyTransformsToGifti<float>(parser);
  }
}

} // namespace ants
