/*
 * antsApplyTransformsToGifti
 *
 * Applies ANTs transforms to the vertices of a GIFTI surface (.gii).
 *
 * Input coordinate contract
 * -------------------------
 * Stored vertices must be physical RAS+ coordinates in millimeters in the
 * scanner/world frame expected by the transforms. Neither FreeSurfer VolGeom*
 * metadata (including C_RAS) nor GIFTI coordinate-system matrices are applied.
 * The POINTSET must have coordinate-system metadata, with every DataSpace set
 * to NIFTI_XFORM_SCANNER_ANAT. Otherwise an exception is thrown. Users certain
 * that the stored coordinates are already scanner RAS can bypass this check
 * with --input-is-scanner-ras, for example for older, incorrectly marked files.
 * This override does not convert coordinates or apply any header transform.
 * The tool converts RAS+ to ITK LPS+ by negating X and Y, applies the composite
 * transform, then converts back to RAS+ in the target image's physical frame.
 *
 * Preparing FreeSurfer surfaces
 * ----------------------------
 * Convert native surfaces explicitly to scanner RAS before using this tool:
 *
 *   mris_convert --to-scanner \
 *     $SUBJECTS_DIR/$SUBJECT/surf/lh.pial lh.pial.scanner.surf.gii
 *
 * Plain mris_convert output is normally in tkregister RAS and is not suitable.
 * Unlike earlier versions of this tool, no automatic C_RAS correction is made.
 * --to-scanner handles the full surface-to-scanner mapping, including volume
 * orientations for which adding C_RAS alone would be insufficient. FreeSurfer
 * 7.4.0 or newer is recommended for conversion and viewing: older converters
 * can write scanner coordinates without the correct GIFTI DataSpace marking.
 *
 * Output metadata
 * ---------------
 * The transformed POINTSET has one CoordinateSystemTransformMatrix with
 * DataSpace = NIFTI_XFORM_SCANNER_ANAT, TransformedSpace = NIFTI_XFORM_UNKNOWN,
 * and an identity matrix. Here scanner RAS denotes the target physical frame;
 * the transform does not identify a particular anatomical template.
 * Original VolGeom* entries are removed from file and data-array metadata.
 * With -r/--reference-image, the output POINTSET receives volume geometry from
 * that 3D target image for FreeSurfer/Freeview compatibility. Without a reference,
 * no VolGeom* entries are written. The reference changes metadata only, never
 * vertex coordinates. ITK LPS geometry is converted to RAS; C_RAS is evaluated
 * at continuous voxel index dimensions/2, following FreeSurfer's convention.
 * Coordinate-system blocks are removed from TRIANGLE arrays, where the GIFTI
 * specification does not allow them. Face topology, other data buffers, label
 * tables, and unrelated metadata are preserved.
 *
 * Transform direction and example
 * -------------------------------
 * Point mapping runs in the opposite direction to image resampling. To move
 * surface vertices from the moving image to the fixed image, use the transforms
 * for resampling the fixed image into moving space. For standard antsRegistration
 * outputs with prefix subjectToTemplate_ (fixed = template, moving = subject):
 *
 *   antsApplyTransformsToGifti \
 *     -i lh.pial.scanner.surf.gii \
 *     -o lh.pial.template.surf.gii \
 *     -r template.nii.gz \
 *     -t '[subjectToTemplate_0GenericAffine.mat,1]' \
 *     -t subjectToTemplate_1InverseWarp.nii.gz
 *
 * Transforms are applied last-specified first, as in antsApplyTransformsToPoints.
 * With no transforms, coordinates are unchanged and output metadata is updated.
 *
 * See also: antsApplyTransforms, antsApplyTransformsToPoints
 */

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkImageIOFactory.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"

// gifticlib is a private dependency of ITK's MeshGifti module; its headers are
// installed alongside the other ITK headers during the ANTs superbuild.
#include "gifti_io.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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

  // Read only the reference header. Its geometry describes the output frame;
  // it does not participate in transforming the vertices.
  std::vector<std::pair<std::string, std::string>> referenceGeometry;
  auto referenceOption = parser->GetOption("reference-image");
  if (referenceOption && referenceOption->GetNumberOfFunctions() > 0)
  {
    const std::string referenceFile = referenceOption->GetFunction(0)->GetName();
    auto imageIO = itk::ImageIOFactory::CreateImageIO(referenceFile.c_str(), itk::IOFileModeEnum::ReadMode);
    if (!imageIO)
    {
      throw std::runtime_error("Cannot read reference image: " + referenceFile);
    }
    imageIO->SetFileName(referenceFile);
    imageIO->ReadImageInformation();
    if (imageIO->GetNumberOfDimensions() != Dimension)
    {
      throw std::runtime_error("Reference image must be 3D: " + referenceFile);
    }
    const auto number = [](double value) {
      if (!std::isfinite(value))
      {
        throw std::runtime_error("Reference image contains non-finite geometry.");
      }
      std::ostringstream stream;
      stream << std::setprecision(17) << value;
      return stream.str();
    };
    referenceGeometry.emplace_back("VolGeomFname", referenceFile);
    const char * dimensions[3] = { "VolGeomWidth", "VolGeomHeight", "VolGeomDepth" };
    const char * axes[3] = { "X", "Y", "Z" };
    const char * ras[3] = { "R", "A", "S" };
    double center[3] = { imageIO->GetOrigin(0), imageIO->GetOrigin(1), imageIO->GetOrigin(2) };
    for (unsigned int axis = 0; axis < Dimension; ++axis)
    {
      const auto size = imageIO->GetDimensions(axis);
      const double spacing = imageIO->GetSpacing(axis);
      if (size == 0 || !std::isfinite(spacing) || spacing <= 0.0)
      {
        throw std::runtime_error("Reference image must have positive dimensions and spacing: " + referenceFile);
      }
      referenceGeometry.emplace_back(dimensions[axis], std::to_string(size));
      referenceGeometry.emplace_back(std::string("VolGeom") + axes[axis] + "size", number(spacing));
      const auto direction = imageIO->GetDirection(axis);
      for (unsigned int component = 0; component < Dimension; ++component)
      {
        const double sign = component < 2 ? -1.0 : 1.0;
        referenceGeometry.emplace_back(std::string("VolGeom") + axes[axis] + "_" + ras[component],
                                       number(sign * direction[component]));
        center[component] += direction[component] * spacing * static_cast<double>(size) / 2.0;
      }
    }
    for (unsigned int component = 0; component < Dimension; ++component)
    {
      const double sign = component < 2 ? -1.0 : 1.0;
      referenceGeometry.emplace_back(std::string("VolGeomC_") + ras[component], number(sign * center[component]));
    }
  }

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
  // Some input files carry a coordsys on the TRIANGLE array. Suppress the
  // validator warning on read; these nonstandard blocks are removed on output.
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

  typename itk::ants::CommandLineParser::OptionType::Pointer scannerRasOption =
    parser->GetOption("input-is-scanner-ras");
  const bool inputIsScannerRas = scannerRasOption && scannerRasOption->GetNumberOfFunctions() > 0 &&
    parser->Convert<bool>(scannerRasOption->GetFunction(0)->GetName());
  if (!inputIsScannerRas)
  {
    // Multiple coordinate-system blocks must agree on the stored data's frame.
    // TransformedSpace describes the result of a header matrix, not the vertices.
    bool scannerDataSpace = da->numCS > 0 && da->coordsys;
    for (int cs = 0; scannerDataSpace && cs < da->numCS; ++cs)
    {
      const giiCoordSystem * csys = da->coordsys[cs];
      scannerDataSpace = csys && csys->dataspace &&
        std::strcmp(csys->dataspace, "NIFTI_XFORM_SCANNER_ANAT") == 0;
    }
    if (!scannerDataSpace)
    {
      gifti_free_image(gim);
      throw std::runtime_error(
        "Input GIFTI POINTSET must declare DataSpace = NIFTI_XFORM_SCANNER_ANAT in every "
        "coordinate-system block: " + inputFile +
        ". Convert native FreeSurfer surfaces with mris_convert --to-scanner. "
        "Use --input-is-scanner-ras only if you are certain the stored vertices are already "
        "in scanner RAS+ physical coordinates; this override does not convert coordinates.");
    }
  }

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
  // Apply transforms to each vertex
  // -----------------------------------------------------------------------
  for (int v = 0; v < nVerts; ++v)
  {
    float * pt = coords + v * 3;

    // RAS+ -> LPS+: the stored vertices are already in scanner/world space.
    typename CompositeTransformType::InputPointType itkPt;
    itkPt[0] = static_cast<RealType>(-pt[0]);
    itkPt[1] = static_cast<RealType>(-pt[1]);
    itkPt[2] = static_cast<RealType>(pt[2]);

    // Apply the ANTs composite transform in LPS+ space.
    typename CompositeTransformType::OutputPointType itkPtOut =
      compositeTransform->TransformPoint(itkPt);

    // LPS+ -> RAS+: negate X and Y back before storing.
    pt[0] = static_cast<float>(-itkPtOut[0]);
    pt[1] = static_cast<float>(-itkPtOut[1]);
    pt[2] = static_cast<float>(itkPtOut[2]);
  }

  // Drop stale FreeSurfer volume geometry at both file and array level.
  // Compact the owned name/value arrays without disturbing unrelated metadata.
  const auto removeVolGeom = [](giiMetaData & meta) {
    int kept = 0;
    for (int i = 0; i < meta.length; ++i)
    {
      if (meta.name[i] && std::strncmp(meta.name[i], "VolGeom", 7) == 0)
      {
        free(meta.name[i]);
        free(meta.value[i]);
      }
      else
      {
        meta.name[kept] = meta.name[i];
        meta.value[kept] = meta.value[i];
        ++kept;
      }
    }
    for (int i = kept; i < meta.length; ++i)
    {
      meta.name[i] = nullptr;
      meta.value[i] = nullptr;
    }
    meta.length = kept;
  };
  removeVolGeom(gim->meta);
  for (int i = 0; i < gim->numDA; ++i)
  {
    removeVolGeom(gim->darray[i]->meta);
    if (gim->darray[i]->intent == NIFTI_INTENT_TRIANGLE)
    {
      gifti_free_CS_list(gim->darray[i]);
    }
  }

  for (const auto & entry : referenceGeometry)
  {
    if (gifti_add_to_meta(&da->meta, entry.first.c_str(), entry.second.c_str(), 1) != 0)
    {
      gifti_free_image(gim);
      throw std::runtime_error("Failed to write reference-image volume geometry to GIFTI metadata.");
    }
  }

  // Replace all old mappings with a single identity mapping describing the
  // stored target-space RAS coordinates. This also handles missing coordsys.
  gifti_free_CS_list(da);
  if (gifti_add_empty_CS(da) != 0)
  {
    std::cerr << "Failed to allocate output GIFTI coordinate system." << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }
  giiCoordSystem * csys = da->coordsys[0];
  csys->dataspace = strdup("NIFTI_XFORM_SCANNER_ANAT");
  csys->xformspace = strdup("NIFTI_XFORM_UNKNOWN");
  if (!csys->dataspace || !csys->xformspace)
  {
    std::cerr << "Failed to allocate output GIFTI coordinate-space names." << std::endl;
    gifti_free_image(gim);
    return EXIT_FAILURE;
  }
  memset(csys->xform, 0, sizeof(csys->xform));
  for (int i = 0; i < 4; ++i)
  {
    csys->xform[i][i] = 1.0;
  }

  // Write the transformed coordinates and updated metadata, preserving all
  // other data buffers and the label table.
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
      "Use double-precision floating point for transform computation (0 = float, "
      "1 = double).  Float is faster; double may improve accuracy for very large "
      "deformation fields.  Vertex coordinates are always stored as float32 in the "
      "output GIFTI regardless of this setting.  Default = 0 (float). ";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("precision");
    option->SetShortName('p');
    option->SetUsageOption(0, "0/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      "Input GIFTI surface file (.gii) with stored vertices in scanner/world RAS+ "
      "physical coordinates (millimeters). Convert native FreeSurfer surfaces "
      "to this format by using mris_convert --to-scanner. "
      "POINTSET DataSpace must be NIFTI_XFORM_SCANNER_ANAT; missing or conflicting "
      "declarations cause an exception unless --input-is-scanner-ras is used. "
      "No C_RAS offset or GIFTI coordinate-system matrix is applied.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputSurface.gii");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input-is-scanner-ras");
    option->SetUsageOption(0, "[0/1]");
    option->SetDescription(
      "Bypasses the POINTSET DataSpace check. Use only when you are certain the stored "
      "vertices are already scanner RAS+ physical coordinates, despite missing or "
      "incorrect metadata (for example, older FreeSurfer GIFTIs). This flag does not "
      "convert coordinates or apply header transforms. Default = 0; using the flag by "
      "itself or with 1 bypasses the check.");
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description =
      "Output GIFTI surface file (.gii), with vertices in target physical RAS+. "
      "POINTSET DataSpace is SCANNER_ANAT, TransformedSpace is UNKNOWN, and the "
      "coordinate-system matrix is identity. Original VolGeom* metadata and TRIANGLE "
      "coordinate-system blocks are removed. Use --reference-image to populate "
      "output VolGeom* from the target image.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "outputSurface.gii");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("reference-image");
    option->SetShortName('r');
    option->SetUsageOption(0, "targetImage.nii.gz");
    option->SetDescription(
      "Optional 3D reference image in the output surface's target physical space. "
      "Populate POINTSET VolGeom* metadata from its dimensions, spacing, RAS direction "
      "cosines, center RAS, and filename for FreeSurfer/Freeview compatibility. "
      "Only the header is read; this option does not change vertex coordinates or "
      "the applied transforms. Without a reference, no VolGeom* metadata is written. ");
    parser->AddOption(option);
  }

  {
    std::string description =
      "An ANTs transforms to apply. Use multiple times to chain transforms, specified in "
      "the same order as antsApplyTransformsToPoints. Use [transformFile,1] to apply the "
      "inverse, for transforms that define an explicit inverse (eg affine transforms)."
      " "
      "Note on transform direction: The required 'forward' or 'inverse' warps for surface vertices "
      "are the OPPOSITE of those used to resample images. For warps from antsRegistration with a given "
      "'fixed' and 'moving' image: to warp a surface defined in the moving-image space into the "
      "fixed-image space, use the same transforms you would use with antsApplyTransforms to warp "
      "the fixed image into moving space. See https://github.com/ANTsX/ANTs/wiki/Applying-transforms-to-point-data";
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
    "Apply ANTs transforms to scanner/world RAS+ vertices of a GIFTI surface file. "
    "Input vertices must already be in scanner/world RAS+ physical coordinates. "
    "The POINTSET DataSpace must be NIFTI_XFORM_SCANNER_ANAT. If this metadata is missing "
    "or incorrect, but you are certain that the mesh points are in RAS+ physical coordinates, "
    "you can force the program to proceed with the --input-is-scanner-ras option. "
    "FreeSurfer 7.4.0 or newer is recommended for conversion and viewing; older "
    "converters may not mark scanner coordinates correctly in DataSpace. "
    "In the output, original VolGeom* metadata is removed. Supply -r/--reference-image "
    "to write geometry from the target image for freeview compatibility. "
    "Without a reference, no VolGeom* metadata is written.";

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
    std::cout << R"(
PREPARING FREESURFER SURFACES:
     mris_convert --to-scanner \
       "$SUBJECTS_DIR/$SUBJECT/surf/lh.pial" lh.pial.surf.gii

EXAMPLE:
     Move a subject surface to the fixed image 'template.nii.gz' using registration outputs:

     antsApplyTransformsToGifti \
       -i lh.pial.surf.gii \
       -o lh.pial.warped.surf.gii \
       -r template.nii.gz \
       -t [ subjectToTemplate_0GenericAffine.mat, 1 ] \
       -t subjectToTemplate_1InverseWarp.nii.gz

)";

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
