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
#include "itkantsRegistrationHelper.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  // short names in use-  a:b:c:d:f:h:l:m:n:o:q:r:s:t:u::w:x:z

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, we try to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the output transform prefix (output format is .nii.gz ). " )
      + std::string( "Optionally, one can choose to warp the moving image to the fixed space and, if the " )
      + std::string( "inverse transform exists, one can also output the warped fixed image.  Note that " )
      + std::string( "only the images specified in the first metric call are warped.  Use antsApplyTransforms " )
      + std::string( "to warp other images using the resultant transform(s)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "outputTransformPrefix" );
    option->SetUsageOption( 1, "[outputTransformPrefix,<outputWarpedImage>,<outputInverseWarpedImage>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Boolean specifying whether or not the " )
      + std::string( "composite transform (and its inverse, if it exists) should " )
      + std::string( "be written to an hdf5 composite file.  This is false by default " )
      + std::string( "so that only the transform for each stage is written to file." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "write-composite-transform" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "1/(0)" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Several interpolation options are available in ITK. " )
      + std::string( "These have all been made available.  Currently the interpolator " )
      + std::string( "choice is only used to warp (and possibly inverse warp) the final " )
      + std::string( "output image(s)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "interpolation" );
    option->SetShortName( 'n' );
    option->SetUsageOption( 0, "Linear" );
    option->SetUsageOption( 1, "NearestNeighbor" );
    option->SetUsageOption( 2, "MultiLabel[<sigma=imageSpacing>,<alpha=4.0>]" );
    option->SetUsageOption( 3, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]" );
    option->SetUsageOption( 4, "BSpline[<order=3>]" );
    option->SetUsageOption( 5, "CosineWindowedSinc" );
    option->SetUsageOption( 6, "WelchWindowedSinc" );
    option->SetUsageOption( 7, "HammingWindowedSinc" );
    option->SetUsageOption( 8, "LanczosWindowedSinc" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

  // Although this option will eventually be used, it is not needed now.

    {
    std::string description = std::string( "Specify the initial fixed transform(s) which get immediately " )
      + std::string( "incorporated into the composite transform.  The order of the " )
      + std::string( "transforms is stack-esque in that the last transform specified on " )
      + std::string( "the command line is the first to be applied.  In addition to specification " )
      + std::string( "of ITK transforms, the user can perform an initial translation alignment " )
      + std::string( "by specifying the fixed and moving images.  These images can then be aligned " )
      + std::string( "using the image intensities, i.e. \"center of mass\", or the geometric center " )
      + std::string( "of the images." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initial-fixed-transform" );
    option->SetShortName( 'q' );
    option->SetUsageOption( 0, "initialTransform" );
    option->SetUsageOption( 1, "[initialTransform,<useInverse>]" );
    option->SetUsageOption( 2, "[fixedImage,movingImage,useCenterOfMass]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

// This is currently not functioning properly for all linear transforms.  If I
// restrict the linear transforms to rigid transforms, then it seems to work.
// I think there's something in working with images that doesn't work properly
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
    std::string description = std::string( "Prints out the CC similarity metric measure " )
      + std::string( "between the full-size input fixed and the transformed moving images at each iteration " )
      + std::string( "a value of 0 (the default) indicates that the full scale computation should not take place")
      + std::string( "any value greater than 0 represents the interval of full scale metric computation.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "print-similarity-measure-interval" );
    option->SetShortName( 'p' );
    option->SetUsageOption( 0, "<unsigned integer value>" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string(
        "Writes out the output volume at each iteration. It helps to present the registration process as a short movie " )
      + std::string( "a value of 0 (the default) indicates that this option should not take place")
      + std::string(
        "any value greater than 0 represents the interval between the iterations which outputs are written to the disk.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "write-interval-volumes" );
    option->SetShortName( 'v' );
    option->SetUsageOption( 0, "<unsigned integer value>" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Collapse output transforms. " )
      + std::string( "Specifically, enabling this option combines all adjacent linear transforms " )
      + std::string( "and composes all adjacent displacement field transforms before writing the " )
      + std::string( "results to disk in the form of an itk affine transform (called xxxGenericAffine.mat). " );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "collapse-output-transforms" );
    option->SetShortName( 'z' );
    option->SetUsageOption( 0, "(1)/0" );
    option->SetDescription( description );
    option->AddFunction( std::string( "1" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the initial moving transform(s) which get immediately " )
      + std::string( "incorporated into the composite transform.  The order of the " )
      + std::string( "transforms is stack-esque in that the last transform specified on " )
      + std::string( "the command line is the first to be applied.  In addition to specification " )
      + std::string( "of ITK transforms, the user can perform an initial translation alignment " )
      + std::string( "by specifying the fixed and moving images.  These images can then be aligned " )
      + std::string( "using the image intensities, i.e. \"center of mass\", or the geometric center " )
      + std::string( "of the images." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initial-moving-transform" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "initialTransform" );
    option->SetUsageOption( 1, "[initialTransform,<useInverse>]" );
    option->SetUsageOption( 2, "[fixedImage,movingImage,useCenterOfMass]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "These image metrics are available--- " )
      + std::string( "CC:  ANTS neighborhood cross correlation, MI:  Mutual information, " )
      + std::string( "Demons: (Thirion), MeanSquares, and " )
      + std::string( "GC: Global Correlation. " )
      + std::string( "Note that the metricWeight is currently not used.  " )
      + std::string( "Rather, it is a temporary place holder until multivariate metrics " )
      + std::string( "are available for a single stage. " )
      + std::string( "The metrics can also employ a sampling strategy defined by a " )
      + std::string( "sampling percentage. The sampling strategy defaults to dense, otherwise " )
      + std::string( "it defines a point set over which to optimize the metric. " )
      + std::string( "The point set can be on a regular lattice or a random lattice of points slightly " )
      + std::string( "perturbed to minimize aliasing artifacts. samplingPercentage defines the " )
      + std::string( "fraction of points to select from the domain. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "metric" );
    option->SetShortName( 'm' );
    option->SetUsageOption(
      0,
      "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      1,
      "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      2,
      "Mattes[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      3,
      "MeanSquares[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      4,
      "Demons[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      5,
      "GC[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Several transform options are available.  The gradientStep or " )
      + std::string( "learningRate characterizes the gradient descent optimization and is scaled appropriately " )
      + std::string( "for each transform using the shift scales estimator.  Subsequent parameters are " )
      + std::string( "transform-specific and can be determined from the usage. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "transform" );
    option->SetShortName( 't' );
    option->SetUsageOption(  0, "Rigid[gradientStep]" );
    option->SetUsageOption(  1, "Affine[gradientStep]" );
    option->SetUsageOption(  2, "CompositeAffine[gradientStep]" );
    option->SetUsageOption(  3, "Similarity[gradientStep]" );
    option->SetUsageOption(  4, "Translation[gradientStep]" );
    option->SetUsageOption(  5, "BSpline[gradientStep,meshSizeAtBaseLevel]" );
    option->SetUsageOption(
      6, "GaussianDisplacementField[gradientStep,updateFieldVarianceInVoxelSpace,totalFieldVarianceInVoxelSpace]" );
    option->SetUsageOption(
      7,
      "BSplineDisplacementField[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
    option->SetUsageOption(
      8,
      "TimeVaryingVelocityField[gradientStep,numberOfTimeIndices,updateFieldVarianceInVoxelSpace,updateFieldTimeVariance,totalFieldVarianceInVoxelSpace,totalFieldTimeVariance]" );
    option->SetUsageOption(
      9,
      "TimeVaryingBSplineVelocityField[gradientStep,velocityFieldMeshSize,<numberOfTimePointSamples=4>,<splineOrder=3>]" );
    option->SetUsageOption( 10, "SyN[gradientStep,updateFieldVarianceInVoxelSpace,totalFieldVarianceInVoxelSpace]" );
    option->SetUsageOption(
      11, "BSplineSyN[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
    option->SetUsageOption(
      12,
      "Exponential[gradientStep,updateFieldVarianceInVoxelSpace,velocityFieldVarianceInVoxelSpace,<numberOfIntegrationSteps>]" );
    option->SetUsageOption(
      13,
      "BSplineExponential[gradientStep,updateFieldMeshSizeAtBaseLevel,velocityFieldMeshSizeAtBaseLevel,<numberOfIntegrationSteps>,<splineOrder=3>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Convergence is determined from the number of iterations per level " )
      + std::string( "and is determined by fitting a line to the normalized energy " )
      + std::string( "profile of the last N iterations (where N is specified by " )
      + std::string( "the window size) and determining the slope which is then " )
      + std::string( "compared with the convergence threshold." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "convergence" );
    option->SetShortName( 'c' );
    option->SetUsageOption( 0, "MxNxO" );
    option->SetUsageOption( 1, "[MxNxO,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the sigma of gaussian smoothing at each level.  " )
      + std::string( "Units are given in terms of voxels (\'vox\') or physical spacing (\'mm\'). " )
      + std::string( "Example usage is \'4x2x1mm\' and \'4x2x1vox\' where no units implies voxel spacing." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "smoothing-sigmas" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "MxNxO..." );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string(
        "Specify the shrink factor for the virtual domain (typically the fixed image) at each level." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "shrink-factors" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "MxNxO..." );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Histogram match the images before registration." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "use-histogram-matching" );
    option->SetShortName( 'u' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string(
        "turn on the option that lets you estimate the learning rate step size only at the beginning of each level.  * useful as a second stage of fine-scale registration." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "use-estimate-learning-rate-once" );
    option->SetShortName( 'l' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Winsorize data based on specified quantiles." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "winsorize-image-intensities" );
    option->SetShortName( 'w' );
    option->SetUsageOption( 0, "[lowerQuantile,upperQuantile]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string         description = "Image masks to limit voxels considered by the metric.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "masks" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "[fixedImageMask,movingImageMask]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }
}

static
const char *
RegTypeToFileName(const std::string & type, bool & writeInverse, bool & writeVelocityField)
{
  std::string str(type);

  ConvertToLowerCase(str);
  if( str == "syn" ||
      str == "symmetricnormalization" ||
      str == "bsplinesyn" ||
      str == "timevaryingbsplinevelocityfield" ||
      str == "tvdmffd" ||
      str == "timevaryingvelocityfield" ||
      str == "tvf"  ||
      str == "exponential" ||
      str == "bsplineexponential" )
    {
    writeInverse = true;
    }
  else
    {
    writeInverse = false;
    }

  if( str == "timevaryingvelocityfield" ||
      str == "tvf" )
    {
    writeVelocityField = true;
    }
  else
    {
    writeVelocityField = false;
    }

  if( str == "rigid" )
    {
    return "Rigid.mat";
    }
  else if( str == "affine" ||
           str == "compositeaffine" || str == "compaff" )
    {
    return "Affine.mat";
    }
  else if( str == "similarity" )
    {
    return "Similarity.mat";
    }
  else if( str == "translation" )
    {
    return "Translation.mat";
    }
  else if( str == "bspline" ||
           str == "ffd" )
    {
    return "BSpline.txt";
    }
  else if( str == "genericaffine" )
    {
    return "GenericAffine.mat";
    }
  else if( str == "gaussiandisplacementfield" ||
           str == "gdf" ||
           str == "bsplinedisplacementfield" ||
           str == "dmffd" ||
           str == "syn" ||
           str == "symmetricnormalization" ||
           str == "bsplinesyn" ||
           str == "timevaryingvelocityfield" ||
           str == "tvf" ||
           str == "timevaryingbsplinevelocityfield" ||
           str == "tvdmffd" ||
           str == "exp" ||
           str == "exponential" ||
           str == "bsplineexponential" )
    {
    return "Warp.nii.gz";
    }
  return "BOGUS.XXXX";
}

template <unsigned VImageDimension>
int
DoRegistration(typename ParserType::Pointer & parser)
{
  typedef typename ants::RegistrationHelper<VImageDimension>      RegistrationHelperType;
  typedef typename RegistrationHelperType::ImageType              ImageType;
  typedef typename RegistrationHelperType::CompositeTransformType CompositeTransformType;

  typename RegistrationHelperType::Pointer regHelper =
    RegistrationHelperType::New();

  OptionType::Pointer transformOption = parser->GetOption( "transform" );

  OptionType::Pointer metricOption = parser->GetOption( "metric" );

  OptionType::Pointer convergenceOption = parser->GetOption( "convergence" );

  OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );

  OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );

  OptionType::Pointer outputOption = parser->GetOption( "output" );

  OptionType::Pointer maskOption = parser->GetOption( "masks" );

  OptionType::Pointer compositeOutputOption = parser->GetOption( "write-composite-transform" );

  OptionType::Pointer collapseOutputTransformsOption = parser->GetOption( "collapse-output-transforms" );

  if( !outputOption || outputOption->GetNumberOfFunctions() == 0 )
    {
    antscout << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
    }

  OptionType::Pointer collapseLinearTransforms =
    parser->GetOption( "collapse-linear-transforms-to-fixed-image-header" );
  if( collapseLinearTransforms && parser->Convert<bool>( collapseLinearTransforms->GetFunction( 0 )->GetName() ) )
    {
    regHelper->SetApplyLinearTransformsToFixedImageHeader( true );
    }
  else
    {
    regHelper->SetApplyLinearTransformsToFixedImageHeader( false );
    }

  OptionType::Pointer printSimilarityMeasureInterval = parser->GetOption( "print-similarity-measure-interval" );
  if( printSimilarityMeasureInterval && printSimilarityMeasureInterval->GetNumberOfFunctions() )
    {
    unsigned int intervalLength = parser->Convert<unsigned int>( printSimilarityMeasureInterval->GetFunction(
                                                                   0 )->GetName() );
    regHelper->SetPrintSimilarityMeasureInterval(intervalLength);
    }
  else
    {
    regHelper->SetPrintSimilarityMeasureInterval( 0 );
    }

  OptionType::Pointer writeIntervalVolumes = parser->GetOption( "write-interval-volumes" );
  if( writeIntervalVolumes && writeIntervalVolumes->GetNumberOfFunctions() )
    {
    unsigned int LengthOfIntervals = parser->Convert<unsigned int>( writeIntervalVolumes->GetFunction( 0 )->GetName() );
    regHelper->SetWriteIntervalVolumes(LengthOfIntervals);
    }
  else
    {
    regHelper->SetPrintSimilarityMeasureInterval( 0 );
    }

  std::string outputPrefix = outputOption->GetFunction( 0 )->GetName();
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
    {
    outputPrefix = outputOption->GetFunction( 0 )->GetParameter( 0 );
    }
  std::string outputWarpedImageName;
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
    {
    outputWarpedImageName = outputOption->GetFunction( 0 )->GetParameter( 1 );
    }

  std::string outputInverseWarpedImageName;
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
    {
    outputInverseWarpedImageName = outputOption->GetFunction( 0 )->GetParameter( 2 );
    }

  ParserType::OptionType::Pointer initialMovingTransformOption = parser->GetOption( "initial-moving-transform" );

  if( initialMovingTransformOption && initialMovingTransformOption->GetNumberOfFunctions() )
    {
    std::vector<bool> isDerivedInitialMovingTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<VImageDimension>( parser, initialMovingTransformOption,
                                                              isDerivedInitialMovingTransform );

    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    regHelper->SetMovingInitialTransform( compositeTransform );

    // Write out initial derived transforms only if we're not collapsing them in the output
    if( !parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      for( unsigned int n = 0; n < isDerivedInitialMovingTransform.size(); n++ )
        {
        std::stringstream curFileName;
        curFileName << outputPrefix << n << "DerivedInitialMovingTranslation.mat";

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer curTransform =
          compositeTransform->GetNthTransform( n );
        itk::ants::WriteTransform<VImageDimension>( curTransform, curFileName.str() );
        }
      }
    }

  ParserType::OptionType::Pointer initialFixedTransformOption = parser->GetOption( "initial-fixed-transform" );

  if( initialFixedTransformOption && initialFixedTransformOption->GetNumberOfFunctions() )
    {
    std::vector<bool> isDerivedInitialFixedTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<VImageDimension>( parser, initialFixedTransformOption,
                                                              isDerivedInitialFixedTransform );
    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    regHelper->SetFixedInitialTransform( compositeTransform );

    // Write out initial derived transforms only if we're not collapsing them in the output
    if( !parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      for( unsigned int n = 0; n < isDerivedInitialFixedTransform.size(); n++ )
        {
        std::stringstream curFileName;
        curFileName << outputPrefix << n << "DerivedInitialFixedTranslation.mat";

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer curTransform =
          compositeTransform->GetNthTransform( n );
        itk::ants::WriteTransform<VImageDimension>( curTransform, curFileName.str() );
        }
      }
    }

  if( maskOption.IsNotNull() && maskOption->GetNumberOfFunctions() )
    {
    typedef typename RegistrationHelperType::MaskImageType MaskImageType;
    typedef itk::ImageFileReader<MaskImageType>            ImageReaderType;
    for( unsigned m = 0; m < maskOption->GetFunction( 0 )->GetNumberOfParameters(); m++ )
      {
      std::string fname = maskOption->GetFunction( 0 )->GetParameter( m );

      typename MaskImageType::Pointer maskImage;
      typename ImageReaderType::Pointer reader = ImageReaderType::New();

      reader->SetFileName( fname.c_str() );
      try
        {
        reader->Update();
        maskImage = reader->GetOutput();
        }
      catch( itk::ExceptionObject & err )
        {
        antscout << "Can't read specified mask image " << fname.c_str() << std::endl;
        antscout << "Exception Object caught: " << std::endl;
        antscout << err << std::endl;
        return EXIT_FAILURE;
        }
      if( m == 0 )
        {
        regHelper->SetFixedImageMask( maskImage );
        }
      else if( m == 1 )
        {
        regHelper->SetMovingImageMask( maskImage );
        }
      }
    }

  // The misc. options
  //  * winsorize image intensities
  //  * use histogram matching
  //  * estimate learning rate
  // are currently specified once on the command line and then apply to all
  // stages.  Advanced parameter specification might require us to rewrite
  // this in the future.

  float lowerQuantile = 0.0;
  float upperQuantile = 1.0;

  bool doWinsorize = false;

  OptionType::Pointer winsorizeOption = parser->GetOption( "winsorize-image-intensities" );
  if( winsorizeOption && winsorizeOption->GetNumberOfFunctions() )
    {
    doWinsorize = true;
    if( winsorizeOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      lowerQuantile = parser->Convert<float>( winsorizeOption->GetFunction( 0 )->GetParameter( 0 ) );
      }
    if( winsorizeOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      upperQuantile = parser->Convert<float>( winsorizeOption->GetFunction( 0 )->GetParameter( 1 ) );
      }
    }
  regHelper->SetWinsorizeImageIntensities( doWinsorize, lowerQuantile, upperQuantile );

  bool doHistogramMatch = false;

  OptionType::Pointer histOption = parser->GetOption( "use-histogram-matching" );
  if( histOption && histOption->GetNumberOfFunctions() )
    {
    std::string histFunction = histOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( histFunction );
    if( histFunction.compare( "1" ) == 0 || histFunction.compare( "true" ) == 0 )
      {
      doHistogramMatch = true;
      }
    }
  regHelper->SetUseHistogramMatching( doHistogramMatch);

  bool doEstimateLearningRateAtEachIteration = true;

  OptionType::Pointer rateOption = parser->GetOption( "use-estimate-learning-rate-once" );
  if( rateOption && rateOption->GetNumberOfFunctions() )
    {
    std::string rateFunction = rateOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( rateFunction );
    if( rateFunction.compare( "1" ) == 0 || rateFunction.compare( "true" ) == 0 )
      {
      doEstimateLearningRateAtEachIteration = false;
      }
    }
  regHelper->SetDoEstimateLearningRateAtEachIteration( doEstimateLearningRateAtEachIteration );

  // We find both the number of transforms and the number of metrics

  unsigned int numberOfTransforms = transformOption->GetNumberOfFunctions();
  if( transformOption.IsNull() || numberOfTransforms == 0 )
    {
    antscout << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<std::vector<unsigned int> > iterationList;
  std::vector<double>                     convergenceThresholdList;
  std::vector<unsigned int>               convergenceWindowSizeList;
  std::vector<std::vector<unsigned int> > shrinkFactorsList;
  std::vector<std::vector<float> >        smoothingSigmasList;
  std::vector<bool>                       smoothingSigmasAreInPhysicalUnitsList;
  std::deque<std::string>                 TransformTypeNames;
  // Each image registration "stage" is characterized by
  //   * a transform
  //   * a set of convergence criteria (number of iterations, convergence threshold,
  //     and/or convergence window)
  //   * a set of shrink factors
  //   * a set of smoothing factors (specified in physical space or voxel space)
  // Note that the set of the number of iterations, the set of shrink factors, and the
  // set of smoothing factors imply the number of levels that will be used for that stage
  // and thus they all need to be the same vector length, e.g. "-c 100x50x10 -f 4x2x1 -s 2x1x0".
  // We use the number of transforms to implicitly guess how many stages are being used in
  // current registration call.  We don't add the metrics in this loop as there could be
  // more than one metric per stage.
  //
  // Also, we iterate backwards because the command line options are stored as a stack (first
  // in last out).
  for( int currentStage = numberOfTransforms - 1; currentStage >= 0; currentStage-- )
    {
    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations;
    double                    convergenceThreshold = 1e-6;
    unsigned int              convergenceWindowSize = 10;
    if( convergenceOption.IsNotNull() && convergenceOption->GetNumberOfFunctions() )
      {
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() == 0 )
        {
        iterations = parser->ConvertVector<unsigned int>( convergenceOption->GetFunction( currentStage )->GetName() );
        }
      else if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 0 )
        {
        iterations =
          parser->ConvertVector<unsigned int>( convergenceOption->GetFunction( currentStage )->GetParameter( 0 ) );
        }
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 1 )
        {
        convergenceThreshold = parser->Convert<double>( convergenceOption->GetFunction( currentStage )->GetParameter(
                                                          1 ) );
        }
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
        {
        convergenceWindowSize = parser->Convert<unsigned int>( convergenceOption->GetFunction(
                                                                 currentStage )->GetParameter( 2 ) );
        const unsigned int minAllowedconvergenceWindowSize = 2; // The BSplineScatteredDataPoints requires at least 2
                                                                // points for interpolation.
        if( convergenceWindowSize < minAllowedconvergenceWindowSize )
          {
          antscout << "Convergence Window Size must be greater than or equal to " << minAllowedconvergenceWindowSize
                   << std::endl;
          }
        }
      }
    else
      {
      antscout << "No convergence criteria are specified." << std::endl;
      return EXIT_FAILURE;
      }

    iterationList.push_back( iterations );
    convergenceThresholdList.push_back( convergenceThreshold );
    convergenceWindowSizeList.push_back( convergenceWindowSize );

    unsigned int numberOfLevels = iterations.size();
    antscout << "  number of levels = " << numberOfLevels << std::endl;

    // Get shrink factors

    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetFunction( currentStage )->GetName() );
    shrinkFactorsList.push_back( factors );

    // Get smoothing sigmas
    std::string  smoothingSigmasString = smoothingSigmasOption->GetFunction( currentStage )->GetName();
    const size_t mmPosition = smoothingSigmasString.find( "mm" );
    const size_t voxPosition = smoothingSigmasString.find( "vox" );
    if( mmPosition != std::string::npos )
      {
      smoothingSigmasString.replace( mmPosition, 2, "" );
      smoothingSigmasAreInPhysicalUnitsList.push_back( true );
      }
    else if( voxPosition != std::string::npos )
      {
      smoothingSigmasString.replace( voxPosition, 3, "" );
      smoothingSigmasAreInPhysicalUnitsList.push_back( false );
      }
    else
      {
      smoothingSigmasAreInPhysicalUnitsList.push_back( false );
      }

    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasString );
    if( sigmas.size() == 1 )
      {
      sigmas.resize( numberOfLevels, sigmas[0] );
      }
    smoothingSigmasList.push_back( sigmas );

    // Set up the optimizer.  To change the iteration number for each level we rely
    // on the command observer.

    float learningRate = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter( 0 ) );

    std::string whichTransform = transformOption->GetFunction( currentStage )->GetName();
    ConvertToLowerCase( whichTransform );

    TransformTypeNames.push_back( whichTransform );

    typename RegistrationHelperType::XfrmMethod xfrmMethod = regHelper->StringToXfrmMethod( whichTransform );

    switch( xfrmMethod )
      {
      case RegistrationHelperType::Affine:
        {
        regHelper->AddAffineTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Rigid:
        {
        regHelper->AddRigidTransform( learningRate );
        }
        break;
      case RegistrationHelperType::CompositeAffine:
        {
        regHelper->AddCompositeAffineTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Similarity:
        {
        regHelper->AddSimilarityTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Translation:
        {
        regHelper->AddTranslationTransform( learningRate );
        }
        break;
      case RegistrationHelperType::GaussianDisplacementField:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 2 ) );
        regHelper->AddGaussianDisplacementFieldTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
        break;
      case RegistrationHelperType::BSplineDisplacementField:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );

        std::vector<unsigned int> meshSizeForTheTotalField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          meshSizeForTheTotalField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
          }
        else
          {
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            meshSizeForTheTotalField.push_back( 0 );
            }
          }

        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }

        regHelper->AddBSplineDisplacementFieldTransform( learningRate, meshSizeForTheUpdateField,
                                                         meshSizeForTheTotalField,
                                                         splineOrder );
        }
        break;
      case RegistrationHelperType::BSpline:
        {
        std::vector<unsigned int> meshSizeAtBaseLevel =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
        regHelper->AddBSplineTransform( learningRate, meshSizeAtBaseLevel );
        }
        break;
      case RegistrationHelperType::TimeVaryingVelocityField:
        {
        unsigned int numberOfTimeIndices = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                            0 )->GetParameter( 1 ) );

        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 2 ) );
        const float varianceForUpdateFieldTime = parser->Convert<float>( transformOption->GetFunction(
                                                                           currentStage )->GetParameter( 3 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 4 ) );
        const float varianceForTotalFieldTime = parser->Convert<float>( transformOption->GetFunction(
                                                                          currentStage )->GetParameter( 5 ) );
        regHelper->AddTimeVaryingVelocityFieldTransform( learningRate,
                                                         numberOfTimeIndices,
                                                         varianceForUpdateField,
                                                         varianceForUpdateFieldTime,
                                                         varianceForTotalField,
                                                         varianceForTotalFieldTime );
        }
        break;
      case RegistrationHelperType::TimeVaryingBSplineVelocityField:
        {
        std::vector<unsigned int> meshSize = parser->ConvertVector<unsigned int>( transformOption->GetFunction(
                                                                                    0 )->GetParameter( 1 ) );

        unsigned int numberOfTimePointSamples = 4;

        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          numberOfTimePointSamples = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 2 ) );
          }
        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }
        regHelper->AddTimeVaryingBSplineVelocityFieldTransform( learningRate,
                                                                meshSize,
                                                                numberOfTimePointSamples,
                                                                splineOrder );
        }
        break;
      case RegistrationHelperType::SyN:
        {
        float varianceForUpdateField = 3.0;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 1 )
          {
          varianceForUpdateField =
            parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
          }
        float varianceForTotalField = 0.0;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          varianceForTotalField = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter(
                                                            2 ) );
          }
        regHelper->AddSyNTransform( learningRate, varianceForUpdateField, varianceForTotalField );
        }
        break;
      case RegistrationHelperType::BSplineSyN:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
        std::vector<unsigned int> meshSizeForTheTotalField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );

        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }

        regHelper->AddBSplineSyNTransform( learningRate, meshSizeForTheUpdateField,
                                           meshSizeForTheTotalField,
                                           splineOrder );
        }
        break;
      case RegistrationHelperType::Exponential:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 1 ) );
        const float varianceForVelocityField = parser->Convert<float>( transformOption->GetFunction(
                                                                         currentStage )->GetParameter( 2 ) );
        unsigned int numberOfIntegrationSteps = 0;  // If the number of integration steps = 0, compute steps
                                                    // automatically
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          numberOfIntegrationSteps = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 3 ) );
          }
        regHelper->AddExponentialTransform( learningRate, varianceForUpdateField, varianceForVelocityField,
                                            numberOfIntegrationSteps );
        }
        break;
      case RegistrationHelperType::BSplineExponential:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );

        std::vector<unsigned int> meshSizeForTheVelocityField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          meshSizeForTheVelocityField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
          }
        else
          {
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            meshSizeForTheVelocityField.push_back( 0 );
            }
          }

        unsigned int numberOfIntegrationSteps = 0;  // If the number of integration steps = 0, compute steps
                                                    // automatically
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          numberOfIntegrationSteps = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 3 ) );
          }

        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 4 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 4 ) );
          }

        regHelper->AddBSplineExponentialTransform( learningRate, meshSizeForTheUpdateField,
                                                   meshSizeForTheVelocityField,
                                                   numberOfIntegrationSteps, splineOrder );
        }
        break;
      default:
        {
        antscout << "Unknown registration method " << "\"" << whichTransform << "\"" << std::endl;
        }
        break;
      }
    }

  // set the vector-vector parameters accumulated
  regHelper->SetIterations( iterationList );
  regHelper->SetConvergenceWindowSizes( convergenceWindowSizeList );
  regHelper->SetConvergenceThresholds( convergenceThresholdList );
  regHelper->SetSmoothingSigmas( smoothingSigmasList );
  regHelper->SetSmoothingSigmasAreInPhysicalUnits( smoothingSigmasAreInPhysicalUnitsList );
  regHelper->SetShrinkFactors( shrinkFactorsList );

  // We iterate through each of the metric "functions" specified on the command
  // line and add it the registration helper.  We also need to assign the stage
  // ID to the added metric.  Multiple metrics for a single stage are specified
  // on the command line by being specified adjacently.

  unsigned int numberOfMetrics = metricOption->GetNumberOfFunctions();
  for( int currentMetricNumber = numberOfMetrics - 1; currentMetricNumber >= 0; currentMetricNumber-- )
    {
    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
    std::string movingImageFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 1 );
    antscout << "  fixed image: " << fixedImageFileName << std::endl;
    antscout << "  moving image: " << movingImageFileName << std::endl;

    typename ImageType::Pointer fixedImage;
    typename ImageType::Pointer movingImage;
    ReadImage<ImageType>( fixedImage,  fixedImageFileName.c_str() );
    ReadImage<ImageType>( movingImage, movingImageFileName.c_str() );
    fixedImage->DisconnectPipeline();
    movingImage->DisconnectPipeline();

    // Get the stage ID
    unsigned int stageID = metricOption->GetFunction( currentMetricNumber )->GetStageID();

    // We check the last stage ID (first iteration) to ensure that the number of stages
    // (as determined by the number of transforms) is equal to the number of stages (as
    // determined by the metrics command line specification).
    if( currentMetricNumber == static_cast<int>( numberOfMetrics - 1 ) )
      {
      if( stageID != numberOfTransforms - 1 )
        {
        ::ants::antscout << "\n\n\n"
                         << "Error:  The number of stages does not match up with the metrics." << std::endl
                         << "The number of transforms is " << numberOfTransforms << " and the last stage ID "
                         << " as determined by the metrics is " << stageID << "." << std::endl;
        return EXIT_FAILURE;
        }
      }

    std::string whichMetric = metricOption->GetFunction( currentMetricNumber )->GetName();
    ConvertToLowerCase( whichMetric );

    float metricWeighting = 1.0;
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 2 )
      {
      metricWeighting = parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 2 ) );
      }

    float samplingPercentage = 1.0;
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 5 )
      {
      samplingPercentage =
        parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 5 ) );
      }

    std::string strategy = "none";
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 4 )
      {
      strategy = metricOption->GetFunction( currentMetricNumber )->GetParameter( 4 );
      }
    ConvertToLowerCase( strategy );

    typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::invalid;
    if( strategy == "random" )
      {
      samplingStrategy = RegistrationHelperType::random;
      }
    else if( strategy == "regular" )
      {
      samplingStrategy = RegistrationHelperType::regular;
      }
    else if( ( strategy == "none" ) || ( strategy == "" ) )
      {
      samplingStrategy = RegistrationHelperType::none;
      }
    else
      {
      samplingStrategy = RegistrationHelperType::invalid;
      std::cout << "ERROR: invalid sampling strategy specified: " << strategy << std::endl;
      return EXIT_FAILURE;
      }

    typename RegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType( whichMetric );

    switch( curMetric )
      {
      case RegistrationHelperType::CC:
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                     currentMetricNumber )->GetParameter( 3 ) );
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              1,
                              radiusOption,
                              samplingPercentage );
        }
        break;
      case RegistrationHelperType::GC:
      case RegistrationHelperType::MeanSquares:
      case RegistrationHelperType::Demons:
        {
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              1,
                              1,
                              samplingPercentage );
        }
        break;
      case RegistrationHelperType::Mattes:
      case RegistrationHelperType::MI:
        {
        unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                  currentMetricNumber )->GetParameter( 3 ) );
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              binOption,
                              1,
                              samplingPercentage );
        }
        break;
      default:
        antscout << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
        return EXIT_FAILURE;
      }
    }

  // Perform the registration

  if( regHelper->DoRegistration() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }
  //
  // write out transforms stored in the composite
  typename CompositeTransformType::Pointer resultTransform = regHelper->GetModifiableCompositeTransform();

  if( parser->Convert<bool>( compositeOutputOption->GetFunction( 0 )->GetName() ) )
    {
    std::string compositeTransformFileName = outputPrefix + std::string( "Composite.h5" );
    std::string inverseCompositeTransformFileName = outputPrefix + std::string( "InverseComposite.h5" );

    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer compositeTransform =
      resultTransform.GetPointer();
    itk::ants::WriteTransform<VImageDimension>( compositeTransform, compositeTransformFileName.c_str() );

    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer inverseCompositeTransform =
      compositeTransform->GetInverseTransform();
    if( inverseCompositeTransform.IsNotNull() )
      {
      itk::ants::WriteTransform<VImageDimension>( inverseCompositeTransform,
                                                  inverseCompositeTransformFileName.c_str() );
      }
    }
  unsigned int numTransforms = resultTransform->GetNumberOfTransforms();

  // write out transforms actually computed, so skip any initial transforms unless
  // we're collapsing the output transforms.

  typedef typename RegistrationHelperType::CompositeTransformType         CompositeTransformType;
  typedef typename CompositeTransformType::Pointer                        CompositeTransformPointer;
  typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;
  typedef typename RegistrationHelperType::TransformType                  TransformType;

  unsigned int              startIndex = initialMovingTransformOption->GetNumberOfFunctions();
  CompositeTransformPointer collapsedResultTransform;
  if( parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
    {
    collapsedResultTransform = regHelper->CollapseCompositeTransform( resultTransform );
    numTransforms = collapsedResultTransform->GetNumberOfTransforms();
    startIndex = 0;
    TransformTypeNames.clear();
    for( unsigned int i = 0; i < numTransforms; i++ )
      {
      if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() == TransformType::Linear )
        {
        // The output type must be Affine, not matrixoffset!  TransformTypeNames.push_back( "matrixoffset" );
        TransformTypeNames.push_back( "genericaffine" );
        }
      else if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() ==
               TransformType::DisplacementField )
        {
        typename DisplacementFieldTransformType::Pointer nthTransform =
          dynamic_cast<DisplacementFieldTransformType *>( collapsedResultTransform->GetNthTransform( i ).GetPointer() );

        // We don't know what set of displacement field transforms were optimized.
        // All we know is whether or not an inverse displacement field exists.  If so,
        // we simply pass a transform name which either does have an inverse or does
        // not.
        if( nthTransform && nthTransform->GetInverseDisplacementField() )
          {
          TransformTypeNames.push_back( "syn" );
          }
        else
          {
          TransformTypeNames.push_back( "gdf" );
          }
        }
      else if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() == TransformType::BSpline )
        {
        TransformTypeNames.push_back( "bspline" );
        }
      }
    }
  for( unsigned int i = startIndex; i < numTransforms; ++i )
    {
    typename CompositeTransformType::TransformTypePointer curTransform;
    if( parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      curTransform = collapsedResultTransform->GetNthTransform( i );
      }
    else
      {
      curTransform = resultTransform->GetNthTransform( i );
      }

    //
    // only registrations not part of the initial transforms in the
    // TransformTypeNames list.
    const std::string curTransformType = TransformTypeNames.front();
    TransformTypeNames.pop_front();

    bool writeInverse;
    bool writeVelocityField;

    std::string transformTemplateName = RegTypeToFileName( curTransformType, writeInverse, writeVelocityField );

    std::stringstream curFileName;
    curFileName << outputPrefix << i << transformTemplateName;

    // WriteTransform will spit all sorts of error messages if it
    // fails, and we want to keep going even if it does so ignore its
    // return value.
    itk::ants::WriteTransform<VImageDimension>( curTransform, curFileName.str() );

    typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::DisplacementFieldType  DisplacementFieldType;
    typename DisplacementFieldTransformType::Pointer dispTransform =
      dynamic_cast<DisplacementFieldTransformType *>(curTransform.GetPointer() );
    // write inverse transform file
    if( writeInverse && dispTransform.IsNotNull() )
      {
      typename DisplacementFieldType::ConstPointer inverseDispField = dispTransform->GetInverseDisplacementField();
      if( inverseDispField.IsNotNull() )
        {
        std::stringstream curInverseFileName;
        curInverseFileName << outputPrefix << i << "InverseWarp.nii.gz";
        typedef itk::ImageFileWriter<DisplacementFieldType> InverseWriterType;
        typename InverseWriterType::Pointer inverseWriter = InverseWriterType::New();
        inverseWriter->SetInput( dispTransform->GetInverseDisplacementField() );
        inverseWriter->SetFileName( curInverseFileName.str().c_str() );
        try
          {
          inverseWriter->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          antscout << "Can't write transform file " << curInverseFileName.str().c_str() << std::endl;
          antscout << "Exception Object caught: " << std::endl;
          antscout << err << std::endl;
          }
        }
      }
    if( writeVelocityField )
      {
      // write velocity field (if applicable)
      typedef typename RegistrationHelperType::TimeVaryingVelocityFieldTransformType
        VelocityFieldTransformType;

      typedef itk::Image<itk::Vector<double, VImageDimension>, VImageDimension + 1> VelocityFieldType;
      typename VelocityFieldTransformType::Pointer velocityFieldTransform =
        dynamic_cast<VelocityFieldTransformType *>(curTransform.GetPointer() );
      if( !velocityFieldTransform.IsNull() )
        {
        std::stringstream curVelocityFieldFileName;
        curVelocityFieldFileName << outputPrefix << i << "VelocityField.nii.gz";

        typedef itk::ImageFileWriter<VelocityFieldType> VelocityFieldWriterType;
        typename VelocityFieldWriterType::Pointer velocityFieldWriter = VelocityFieldWriterType::New();
        velocityFieldWriter->SetInput( velocityFieldTransform->GetTimeVaryingVelocityField() );
        velocityFieldWriter->SetFileName( curVelocityFieldFileName.str().c_str() );
        try
          {
          velocityFieldWriter->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          antscout << "Can't write velocity field transform file " << curVelocityFieldFileName.str().c_str()
                   << std::endl;
          antscout << "Exception Object caught: " << std::endl;
          antscout << err << std::endl;
          }
        }
      }
    }

  typedef double RealType;
  std::string whichInterpolator( "linear" );
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfFunctions() )
    {
    whichInterpolator = interpolationOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( whichInterpolator );
    }

  typename ImageType::SpacingType cache_spacing_for_smoothing_sigmas;
  if( !std::strcmp( whichInterpolator.c_str(), "gaussian" )
      ||   !std::strcmp( whichInterpolator.c_str(), "multilabel" )
      )
    {
#if 1
    // HACK:: This can just be cached when reading the fixedImage from above!!
    //
    std::string fixedImageFileName = metricOption->GetFunction( numberOfTransforms - 1 )->GetParameter( 0 );

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();

    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
#endif
    cache_spacing_for_smoothing_sigmas = fixedImage->GetSpacing();
    }

#include "make_interpolator_snip.tmpl"
  regHelper->SetInterpolator( interpolator );

  typename ImageType::Pointer warpedImage = regHelper->GetWarpedImage();
  typedef itk::ImageFileWriter<ImageType> WarpedImageWriterType;
  if( !outputWarpedImageName.empty() )
    {
    WriteImage<ImageType>( warpedImage, outputWarpedImageName.c_str()  );
    }

  if( !outputInverseWarpedImageName.empty() )
    {
    typename ImageType::Pointer inverseWarpedImage = regHelper->GetInverseWarpedImage();
    if( inverseWarpedImage.IsNotNull() )
      {
      WriteImage<ImageType>( inverseWarpedImage, outputInverseWarpedImageName.c_str()  );
      }
    }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsRegistration( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  try
    {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert( args.begin(), "antsRegistration" );
    int     argc = args.size();
    char* * argv = new char *[args.size() + 1];
    for( unsigned int i = 0; i < args.size(); ++i )
      {
      // allocate space for the string plus a null character
      argv[i] = new char[args[i].length() + 1];
      std::strncpy( argv[i], args[i].c_str(), args[i].length() );
      // place the null character in the end
      argv[i][args[i].length()] = '\0';
      }
    argv[argc] = 0;
    // class to automatically cleanup argv upon destruction
    class Cleanup_argv
    {
public:
      Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
      {
      }

      ~Cleanup_argv()
      {
        for( unsigned int i = 0; i < argc_plus_one; ++i )
          {
          delete[] argv[i];
          }
        delete[] argv;
      }

private:
      char* *      argv;
      unsigned int argc_plus_one;
    };
    Cleanup_argv cleanup_argv( argv, argc + 1 );

    antscout->set_stream( out_stream );

    ParserType::Pointer parser = ParserType::New();

    parser->SetCommand( argv[0] );

    std::string commandDescription = std::string( "This program is a user-level " )
      + std::string( "registration application meant to utilize ITKv4-only classes. The user can specify " )
      + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
      + std::string( "and iterations, shrink factors, and smoothing sigmas for each level." );

    parser->SetCommandDescription( commandDescription );
    InitializeCommandLineOptions( parser );

    parser->Parse( argc, argv );

    if( argc < 2 || parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
      {
      parser->PrintMenu( antscout, 5, false );
      if( argc < 2 )
        {
        return EXIT_FAILURE;
        }
      return EXIT_SUCCESS;
      }
    else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
      {
      parser->PrintMenu( antscout, 5, true );
      return EXIT_SUCCESS;
      }
    unsigned int dimension = 3;

    ParserType::OptionType::Pointer dimOption = parser->GetOption( "dimensionality" );
    if( dimOption && dimOption->GetNumberOfFunctions() )
      {
      dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
      }
    else
      {
      antscout << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
      return EXIT_FAILURE;
      }

    switch( dimension )
      {
      case 2:
        {
        return DoRegistration<2>(parser);
        }
      case 3:
        {
        return DoRegistration<3>(parser);
        }
      default:
        antscout << "bad image dimension " << dimension << std::endl;
        return EXIT_FAILURE;
      }
    }
  catch( itk::ExceptionObject & err )
    {
    antscout << "Exception Object caught: " << std::endl;
    antscout << err << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
