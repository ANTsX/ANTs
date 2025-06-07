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

#include <iostream>
#include <vector>
#include <string>
#include "antsRegistrationTemplateHeader.h"

#include "ANTsVersion.h"

namespace ants
{

static void
antsRegistrationInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{

  // short names in use-  a:b:c:d:f:g:h:i:j:k:l:m:n:o:q:r:s:t:u:v:w:x:z
  {
    const std::string description = std::string("Get Version Information.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("version");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("This option forces the image to be treated as a specified-") +
                              std::string("dimensional image.  If not specified, we try to ") +
                              std::string("infer the dimensionality from the input image.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3/4");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the output transform prefix (output format is .nii.gz ). ") +
      std::string("Optionally, one can choose to warp the moving image to the fixed space and, if the ") +
      std::string("inverse transform exists, one can also output the warped fixed image.  Note that ") +
      std::string("only the images specified in the first metric call are warped.  Use antsApplyTransforms ") +
      std::string(
        "to warp other images using the resultant transform(s). When a composite transform is not specified, ") +
      std::string("linear transforms are specified with a \'.mat\' suffix and displacement fields with a ") +
      std::string("\'Warp.nii.gz\' suffix (and \'InverseWarp.nii.gz\', when applicable.  In addition, for ") +
      std::string("velocity-based transforms, the full velocity field is written to file (\'VelocityField.nii.gz\') as "
                  "long as the ") +
      std::string("collapse transforms flag is turned off (\'-z 0\').");
    ;

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "outputTransformPrefix");
    option->SetUsageOption(1, "[outputTransformPrefix,<outputWarpedImage>,<outputInverseWarpedImage>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the output file for the current state of the registration. ") +
      std::string("The state file is written to an hdf5 composite file. It is specially usefull if ") +
      std::string("we want to save the current state of a SyN registration to the disk, so ") +
      std::string("we can load and restore that later to continue the next registration process ") +
      std::string("directly started from the last saved state. ") +
      std::string("The output file of this flag is the same as the write-composite-transform, ") +
      std::string("unless the last transform is a SyN transform. In that case, the inverse ") +
      std::string("displacement field of the SyN transform is also added to the output composite transform. ") +
      std::string("Again notice that this file cannot be treated as a transform, and restore-state option ") +
      std::string("must be used to load the written file by this flag.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("save-state");
    option->SetShortName('j');
    option->SetUsageOption(0, "saveStateAsTransform");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the initial state of the registration which get immediately ") +
      std::string("used to directly initialize the registration process. ") +
      std::string("The flag is mutually exclusive with other intialization flags.") +
      std::string("If this flag is used, none of the initial-moving-transform and initial-fixed-transform ") +
      std::string("cannot be used.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("restore-state");
    option->SetShortName('k');
    option->SetUsageOption(0, "restoreStateAsATransform");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Boolean specifying whether or not the ") +
                              std::string("composite transform (and its inverse, if it exists) should ") +
                              std::string("be written to an hdf5 composite file.  This is false by default ") +
                              std::string("so that only the transform for each stage is written to file.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("write-composite-transform");
    option->SetShortName('a');
    option->SetUsageOption(0, "1/(0)");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  // This is currently not functioning properly for all linear transforms.  If I
  // restrict the linear transforms to rigid transforms, then it seems to work.
  // I think there's something in working with images that don't work properly
  // with a generic affine transform in the header.  You can certainly store it
  // and read it from the header but perhaps this interferes with something fundamental
  // like transforming indices to physical coordinates.  I'll have to investigate
  // in the future.

  //     {
  //     std::string description = std::string( "Collapse initial linear transforms " )
  //       + std::string( "to the fixed image header.  This should speed up subsequent " )
  //       + std::string( "nonlinear transform optimizations." );
  //     OptionType::Pointer option = OptionType::New();
  //     option->SetLongName( "collapse-linear-transforms-to-fixed-image-header" );
  //     option->SetShortName( 'b' );
  //     option->SetUsageOption( 0, "1/(0)" );
  //     option->SetDescription( description );
  //     option->AddFunction( std::string( "0" ) );
  //     parser->AddOption( option );
  //     }

  {
    std::string description =
      std::string("Prints out the CC similarity metric measure ") +
      std::string("between the full-size input fixed and the transformed moving images at each iteration ") +
      std::string("a value of 0 (the default) indicates that the full scale computation should not take place") +
      std::string("any value greater than 0 represents the interval of full scale metric computation.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("print-similarity-measure-interval");
    option->SetShortName('p');
    option->SetUsageOption(0, "<unsignedIntegerValue>");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Writes out the output volume at each iteration. It helps to present the registration process as a "
                  "short movie ") +
      std::string("a value of 0 (the default) indicates that this option should not take place") +
      std::string("any value greater than 0 represents the interval between the iterations which outputs are written "
                  "to the disk.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("write-interval-volumes");
    //  option->SetShortName( 'v' ); // BUG! dont set as v because v is verbose
    option->SetUsageOption(0, "<unsignedIntegerValue>");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Collapse output transforms. ") +
      std::string("Specifically, enabling this option combines all adjacent transforms where") +
      std::string("possible.  All adjacent linear transforms are written to disk in the form") +
      std::string("an itk affine transform (called xxxGenericAffine.mat).  Similarly, all ") +
      std::string("adjacent displacement field transforms are combined when written to disk ") +
      std::string("(e.g. xxxWarp.nii.gz and xxxInverseWarp.nii.gz (if available)).") +
      std::string("Also, an output composite transform including the collapsed transforms is ") +
      std::string("written to the disk (called outputCollapsed(Inverse)Composite).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("collapse-output-transforms");
    option->SetShortName('z');
    option->SetUsageOption(0, "(1)/0");
    option->SetDescription(description);
    option->AddFunction(std::string("1"));
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Initialize linear transforms from the previous stage. ") +
      std::string("By enabling this option, the current linear stage transform is directly intialized ") +
      std::string("from the previous stage's linear transform; this allows multiple linear stages to be run ") +
      std::string("where each stage directly updates the estimated linear transform from the previous stage. ") +
      std::string("(e.g. Translation -> Rigid -> Affine).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initialize-transforms-per-stage");
    option->SetShortName('i');
    option->SetUsageOption(0, "(1)/0");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Several interpolation options are available in ITK. ") +
                              std::string("These have all been made available.  Currently the interpolator ") +
                              std::string("choice is only used to warp (and possibly inverse warp) the final ") +
                              std::string("output image(s).");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("interpolation");
    option->SetShortName('n');
    option->SetUsageOption(0, "Linear");
    option->SetUsageOption(1, "NearestNeighbor");
    option->SetUsageOption(2, "MultiLabel[<sigma=imageSpacing>,<alpha=4.0>]");
    option->SetUsageOption(3, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]");
    option->SetUsageOption(4, "BSpline[<order=3>]");
    option->SetUsageOption(5, "CosineWindowedSinc");
    option->SetUsageOption(6, "WelchWindowedSinc");
    option->SetUsageOption(7, "HammingWindowedSinc");
    option->SetUsageOption(8, "LanczosWindowedSinc");
    option->SetUsageOption(9, "GenericLabel[<interpolator=Linear>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("This option allows the user to restrict the ") +
      std::string("optimization of the displacement field, translation, rigid or affine ") +
      std::string("transform on a per-component basis.  For example, if one wants to limit ") +
      std::string("the deformation or rotation of 3-D volume to the first two dimensions, ") +
      std::string("this is possible by specifying a weight vector of \'1x1x0\' for a ") +
      std::string("deformation field or  \'1x1x0x1x1x0\' for a rigid transformation. ") +
      std::string("Low-dimensional restriction only works if there are no preceding transformations.") +
      std::string("All stages up to and including the desired stage must have this option specified,") +
      std::string("even if they should not be restricted (in which case specify 1x1x1...)");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("restrict-deformation");
    option->SetShortName('g');
    option->SetUsageOption(0, "PxQxR");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the initial fixed transform(s) which get immediately ") +
      std::string("incorporated into the composite transform.  The order of the ") +
      std::string("transforms is stack-esque in that the last transform specified on ") +
      std::string("the command line is the first to be applied.  In addition to initialization ") +
      std::string("with ITK transforms, the user can perform an initial translation alignment ") +
      std::string("by specifying the fixed and moving images and selecting an initialization ") +
      std::string("feature.  These features include using the geometric center of the images (=0), ") +
      std::string("the image intensities (=1), or the origin of the images (=2).");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initial-fixed-transform");
    option->SetShortName('q');
    option->SetUsageOption(0, "initialTransform");
    option->SetUsageOption(1, "[initialTransform,<useInverse>]");
    option->SetUsageOption(2, "[fixedImage,movingImage,initializationFeature]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the initial moving transform(s) which get immediately ") +
      std::string("incorporated into the composite transform.  The order of the ") +
      std::string("transforms is stack-esque in that the last transform specified on ") +
      std::string("the command line is the first to be applied.  In addition to initialization ") +
      std::string("with ITK transforms, the user can perform an initial translation alignment ") +
      std::string("by specifying the fixed and moving images and selecting an initialization ") +
      std::string("feature.  These features include using the geometric center of the images (=0), ") +
      std::string("the image intensities (=1), or the origin of the images (=2).");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("initial-moving-transform");
    option->SetShortName('r');
    option->SetUsageOption(0, "initialTransform");
    option->SetUsageOption(1, "[initialTransform,<useInverse>]");
    option->SetUsageOption(2, "[fixedImage,movingImage,initializationFeature]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("These image metrics are available--- ") +
      std::string("CC:  ANTS neighborhood cross correlation, MI:  Mutual information, ") +
      std::string("Demons: (Thirion), MeanSquares, and GC: Global Correlation. ") +
      std::string("The \"metricWeight\" variable is used to modulate the per stage weighting of the metrics.  ") +
      std::string("The metrics can also employ a sampling strategy defined by a ") +
      std::string("sampling percentage. The sampling strategy defaults to \'None\' (aka a dense sampling of ") +
      std::string("one sample per voxel), otherwise it defines a point set over which to optimize the metric. ") +
      std::string("The point set can be on a regular lattice or a random lattice of points slightly ") +
      std::string("perturbed to minimize aliasing artifacts. samplingPercentage defines the ") +
      std::string("fraction of points to select from the domain. useGradientFilter specifies whether a smoothing") +
      std::string("filter is applied when estimating the metric gradient.") +
      std::string("In addition, three point set metrics are available:  Euclidean ") +
      std::string("(ICP), Point-set expectation (PSE), and Jensen-Havrda-Charvet-Tsallis (JHCT).");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("metric");
    option->SetShortName('m');
    option->SetUsageOption(0,
                           "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={None,Regular,Random}>,<"
                           "samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(1,
                           "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,Random}"
                           ">,<samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(2,
                           "Mattes[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(3,
                           "MeanSquares[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(4,
                           "Demons[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(5,
                           "GC[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,Random}>,<"
                           "samplingPercentage=[0,1]>,<useGradientFilter=false>]");
    option->SetUsageOption(
      6, "ICP[fixedPointSet,movingPointSet,metricWeight,<samplingPercentage=[0,1]>,<boundaryPointsOnly=0>]");
    option->SetUsageOption(7,
                           "PSE[fixedPointSet,movingPointSet,metricWeight,<samplingPercentage=[0,1]>,<"
                           "boundaryPointsOnly=0>,<pointSetSigma=1>,<kNeighborhood=50>]");
    option->SetUsageOption(
      8,
      "JHCT[fixedPointSet,movingPointSet,metricWeight,<samplingPercentage=[0,1]>,<boundaryPointsOnly=0>,<pointSetSigma="
      "1>,<kNeighborhood=50>,<alpha=1.1>,<useAnisotropicCovariances=1>]");
    option->SetUsageOption(9,
                           "IGDM[fixedImage,movingImage,metricWeight,fixedMask,movingMask,<neighborhoodRadius=0x0>,<"
                           "intensitySigma=0>,<distanceSigma=0>,<kNeighborhood=1>,<gradientSigma=1>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Several transform options are available.  The gradientStep or ") +
      std::string("learningRate characterizes the gradient descent optimization and is scaled appropriately ") +
      std::string("for each transform using the shift scales estimator.  Subsequent parameters are ") +
      std::string("transform-specific and can be determined from the usage. For the B-spline transforms ") +
      std::string("one can also specify the smoothing in terms of spline distance (i.e. knot spacing). ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transform");
    option->SetShortName('t');
    option->SetUsageOption(0, "Rigid[gradientStep]");
    option->SetUsageOption(1, "Affine[gradientStep]");
    option->SetUsageOption(2, "CompositeAffine[gradientStep]");
    option->SetUsageOption(3, "Similarity[gradientStep]");
    option->SetUsageOption(4, "Translation[gradientStep]");
    option->SetUsageOption(5, "BSpline[gradientStep,meshSizeAtBaseLevel]");
    option->SetUsageOption(
      6, "GaussianDisplacementField[gradientStep,updateFieldVarianceInVoxelSpace,totalFieldVarianceInVoxelSpace]");
    option->SetUsageOption(7,
                           "BSplineDisplacementField[gradientStep,updateFieldMeshSizeAtBaseLevel,<"
                           "totalFieldMeshSizeAtBaseLevel=0>,<splineOrder=3>]");
    option->SetUsageOption(8,
                           "TimeVaryingVelocityField[gradientStep,numberOfTimeIndices,updateFieldVarianceInVoxelSpace,"
                           "updateFieldTimeVariance,totalFieldVarianceInVoxelSpace,totalFieldTimeVariance]");
    option->SetUsageOption(9,
                           "TimeVaryingBSplineVelocityField[gradientStep,velocityFieldMeshSize,<"
                           "numberOfTimePointSamples=4>,<splineOrder=3>]");
    option->SetUsageOption(10,
                           "SyN[gradientStep,<updateFieldVarianceInVoxelSpace=3>,<totalFieldVarianceInVoxelSpace=0>]");
    option->SetUsageOption(
      11, "BSplineSyN[gradientStep,updateFieldMeshSizeAtBaseLevel,<totalFieldMeshSizeAtBaseLevel=0>,<splineOrder=3>]");
    option->SetUsageOption(12,
                           "Exponential[gradientStep,updateFieldVarianceInVoxelSpace,velocityFieldVarianceInVoxelSpace,"
                           "<numberOfIntegrationSteps>]");
    option->SetUsageOption(13,
                           "BSplineExponential[gradientStep,updateFieldMeshSizeAtBaseLevel,<"
                           "velocityFieldMeshSizeAtBaseLevel=0>,<numberOfIntegrationSteps>,<splineOrder=3>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Convergence is determined from the number of iterations per level ") +
                              std::string("and is determined by fitting a line to the normalized energy ") +
                              std::string("profile of the last N iterations (where N is specified by ") +
                              std::string("the window size) and determining the slope which is then ") +
                              std::string("compared with the convergence threshold.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("convergence");
    option->SetShortName('c');
    option->SetUsageOption(0, "MxNxO");
    option->SetUsageOption(1, "[MxNxO,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the sigma of gaussian smoothing at each level.  ") +
      std::string(R"(Units are given in terms of voxels ('vox') or physical spacing ('mm'). )") +
      std::string(R"(Example usage is '4x2x1mm' and '4x2x1vox' where no units implies voxel spacing.)");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("smoothing-sigmas");
    option->SetShortName('s');
    option->SetUsageOption(0, "MxNxO...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the shrink factor for the virtual domain (typically the fixed image) at each level.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("shrink-factors");
    option->SetShortName('f');
    option->SetUsageOption(0, "MxNxO...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Histogram match the images before registration.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("use-histogram-matching");
    option->SetShortName('u');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  // {
  //   std::string description =
  //     std::string("Turn on the option that lets you estimate the learning rate step size only at the beginning of each "
  //                 "level.  This is useful as a second stage of fine-scale registration.");
  //   OptionType::Pointer option = OptionType::New();
  //   option->SetLongName("use-estimate-learning-rate-once");
  //   option->SetShortName('l');
  //   option->SetDescription(description);
  //   parser->AddOption(option);
  // }

  {
    std::string description = std::string("Winsorize data based on specified quantiles.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("winsorize-image-intensities");
    option->SetShortName('w');
    option->SetUsageOption(0, "[lowerQuantile,upperQuantile]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Image masks to limit voxels considered by the metric. ") +
                              std::string("Two options are allowed for mask specification:  1) Either ") +
                              std::string("the user specifies a single mask to be used for all stages or ") +
                              std::string("2) the user specifies a mask for each stage.  With the latter ") +
                              std::string("one can select to which stages masks are applied by supplying ") +
                              std::string("valid file names.  If the file does not exist, a mask will not ") +
                              std::string("be used for that stage.  Note that we handle the fixed and moving ") +
                              std::string("masks separately to enforce this constraint.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("masks");
    option->SetShortName('x');
    option->SetUsageOption(0, "[fixedImageMask,movingImageMask]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Use 'float' instead of 'double' for computations.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("float");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Use MINC file formats for transformations.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("minc");
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Use a fixed seed for random number generation. ") +
                              std::string("By default, the system clock is used to initialize the seeding. ") +
                              std::string("The fixed seed can be any nonzero int value.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("random-seed");
    option->SetUsageOption(0, "seedValue");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Verbose output.");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('v');
    option->SetLongName("verbose");
    option->SetUsageOption(0, "(0)/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu.  Will also print values ") +
                              std::string("used on the current command line call.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'

int
antsRegistration(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  try
  {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert(args.begin(), "antsRegistration");
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

    //    // antscout->set_stream( out_stream );

    ParserType::Pointer parser = ParserType::New();

    parser->SetCommand(argv[0]);

    std::string commandDescription =
      std::string("This program is a user-level ") +
      std::string("registration application meant to utilize classes in ITK v4.0 and later. The user can specify ") +
      std::string("any number of \"stages\" where a stage consists of a transform; an image metric; ") +
      std::string("and iterations, shrink factors, and smoothing sigmas for each level.  ") +
      std::string("Note that explicitly setting the dimensionality, metric, transform, output, ") +
      std::string("convergence, shrink-factors, and smoothing-sigmas parameters is mandatory.");

    parser->SetCommandDescription(commandDescription);
    antsRegistrationInitializeCommandLineOptions(parser);

    if (parser->Parse(argc, argv) == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
    }

    bool                                              verbose = false;
    itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
    if (verboseOption && verboseOption->GetNumberOfFunctions())
    {
      verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
    }

    OptionType::Pointer collapseOutputTransformsOption = parser->GetOption("collapse-output-transforms");
    OptionType::Pointer compositeOutputOption = parser->GetOption("write-composite-transform");
    OptionType::Pointer initializePerStageOption = parser->GetOption("initialize-transforms-per-stage");
    OptionType::Pointer saveStateOption = parser->GetOption("save-state");

    const bool writeCompositeTransform = parser->Convert<bool>(compositeOutputOption->GetFunction(0)->GetName());
    const bool shouldInitializePerStage = parser->Convert<bool>(initializePerStageOption->GetFunction(0)->GetName());
    if (shouldInitializePerStage && (!writeCompositeTransform))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  --initialize-transforms-per-stage requires --write-composite-transform" << std::endl;
        std::cerr << "        because the initializizing transform is collapsed into each stage for optimization"
                  << std::endl;
      }
      return EXIT_FAILURE;
    }
    if ((saveStateOption && saveStateOption->GetNumberOfFunctions()) && (!writeCompositeTransform))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  --save-state requires --write-composite-transform" << std::endl;
        std::cerr << "        because the the output transform will contain the this processes initializer"
                  << std::endl;
      }
      return EXIT_FAILURE;
    }

    if (verbose)
    {
      std::cout << "All_Command_lines_OK" << std::endl;
    }

    if (argc == 1)
    {
      parser->PrintMenu(std::cout, 5, false);
      return EXIT_FAILURE;
    }
    else if (parser->GetOption("help")->GetFunction() &&
             parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
    {
      parser->PrintMenu(std::cout, 5, false);
      return EXIT_SUCCESS;
    }
    else if (parser->GetOption('h')->GetFunction() &&
             parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
    {
      parser->PrintMenu(std::cout, 5, true);
      return EXIT_SUCCESS;
    }
    ParserType::OptionType::Pointer versionOption = parser->GetOption("version");
    if (versionOption && versionOption->GetNumberOfFunctions())
    {
      std::string versionFunction = versionOption->GetFunction(0)->GetName();
      ConvertToLowerCase(versionFunction);
      if (versionFunction.compare("1") == 0 || versionFunction.compare("true") == 0)
      {
        // Print Version Information
        std::cout << ANTs::Version::ExtendedVersionString() << std::endl;
        return EXIT_SUCCESS;
      }
    }

    unsigned int dimension = 3;

    ParserType::OptionType::Pointer dimOption = parser->GetOption("dimensionality");
    if (dimOption && dimOption->GetNumberOfFunctions())
    {
      dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
    }
    else
    {
      if (verbose)
      {
        std::cerr << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
      }
      return EXIT_FAILURE;
    }

    std::string precisionType;
    ParserType::OptionType::Pointer typeOption = parser->GetOption("float");
    if (typeOption && typeOption->GetNumberOfFunctions() && parser->Convert<bool>(typeOption->GetFunction(0)->GetName()))
    {
      if (verbose)
      {
        std::cout << "Using single precision for computations." << std::endl;
      }
      precisionType = std::string("float");
    }
    else
    {
      if (verbose)
      {
        std::cout << "Using double precision for computations." << std::endl;
      }
      precisionType = std::string("double");
    }

    switch (dimension)
    {
      case 2:
      {
        if (strcmp(precisionType.c_str(), "float") == 0)
        {
          return antsRegistration2DFloat(parser);
        }
        else
        {
          return antsRegistration2DDouble(parser);
        }
      }
      case 3:
      {
        if (strcmp(precisionType.c_str(), "float") == 0)
        {
          return antsRegistration3DFloat(parser);
        }
        else
        {
          return antsRegistration3DDouble(parser);
        }
      }
      case 4:
      {
        if (strcmp(precisionType.c_str(), "float") == 0)
        {
          return antsRegistration4DFloat(parser);
        }
        else
        {
          return antsRegistration4DDouble(parser);
        }
      }
      default:
      {
        if (verbose)
        {
          std::cerr << "bad image dimension " << dimension << std::endl;
        }
        return EXIT_FAILURE;
      }
    }
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "Exception Object caught: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
