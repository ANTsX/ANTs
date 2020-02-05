
#include "antsUtilities.h"
#include <algorithm>

#include "antsCommandLineParser.h"
#include "antsAllocImage.h"
#include "itkArray2D.h"
#include "itkDecomposeTensorFunction.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTimeProbe.h"
#include "itkVariableSizeMatrix.h"

#include <itksys/SystemTools.hxx>

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include <iomanip>

#include <string>
#include <algorithm>
#include <vector>
#include <fstream>

namespace ants
{
template <typename TensorType>
double CalculateFractionalAnisotropy( TensorType tensor )
{
  typename TensorType::EigenValuesArrayType eigenvalues;
  typename TensorType::EigenVectorsMatrixType eigenvectors;

  tensor.ComputeEigenAnalysis( eigenvalues, eigenvectors );

  if( eigenvalues[0] < 0 )
    {
    eigenvalues[0] = eigenvalues[1];
    }
  if( TensorType::Dimension == 3 && eigenvalues[2] < 0 )
    {
    eigenvalues[2] = eigenvalues[1];
    }

  double fa = 0.0;
  double mean = tensor.GetTrace() / static_cast<double>(
      TensorType::Dimension );

  double numerator = 0.0;
  double denominator = 0.0;
  for( unsigned int d = 0; d < TensorType::Dimension; d++ )
    {
    numerator += static_cast<double>( itk::Math::sqr( static_cast<double>( eigenvalues[d] ) - mean ) );
    denominator += static_cast<double>( itk::Math::sqr( eigenvalues[d] ) );
    }
  fa = std::sqrt( ( 3.0 * numerator ) / ( 2.0 * denominator ) );

  return fa;
}

template <typename TensorType>
double CalculateMeanDiffusivity( TensorType tensor )
{
  typename TensorType::EigenValuesArrayType eigenvalues;
  typename TensorType::EigenVectorsMatrixType eigenvectors;

  tensor.ComputeEigenAnalysis( eigenvalues, eigenvectors );

  if( eigenvalues[0] < 0 )
    {
    eigenvalues[0] = eigenvalues[1];
    }
  if( TensorType::Dimension == 3 && eigenvalues[2] < 0 )
    {
    eigenvalues[2] = eigenvalues[1];
    }

  double mean = tensor.GetTrace() / static_cast<double>(
      TensorType::Dimension );

  return mean;
}

template <unsigned int ImageDimension>
int CreateDTICohort( itk::ants::CommandLineParser *parser )
{
  typedef float                                                    RealType;
  typedef itk::SymmetricSecondRankTensor<RealType, ImageDimension> TensorType;
  typedef itk::VariableSizeMatrix<typename TensorType::ValueType>  MatrixType;

  typedef itk::Image<RealType, ImageDimension> ImageType;

  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage = nullptr;

  typedef itk::Image<TensorType, ImageDimension> TensorImageType;
  typename TensorImageType::Pointer inputAtlas = nullptr;

  typedef itk::ImageFileReader<TensorImageType> TensorReaderType;
  typename TensorReaderType::Pointer reader = TensorReaderType::New();

  typedef itk::DecomposeTensorFunction<MatrixType,
                                       typename MatrixType::ValueType, MatrixType> DecomposerType;
  typename DecomposerType::Pointer decomposer = DecomposerType::New();

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
  typename RandomizerType::Pointer randomizer = RandomizerType::New();
  randomizer->Initialize();

  //
  // Get the input DTI atlas
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer inputAtlasOption =
    parser->GetOption( "dti-atlas" );
  if( inputAtlasOption && inputAtlasOption->GetNumberOfFunctions() )
    {
    std::string inputFile = inputAtlasOption->GetFunction( 0 )->GetName();
    reader->SetFileName( inputFile.c_str() );

    inputAtlas = reader->GetOutput();
    inputAtlas->Update();
    inputAtlas->DisconnectPipeline();
    }
  else
    {
    std::cout << "ERROR:  Input DTI atlas not specified." << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Get the number of output images and duplicate the atlas for each cohort.
  //
  std::string  outputDirectory( "./" );
  std::string  rootOutputFileName( "outputDWI" );
  unsigned int numberOfControls = 10;
  unsigned int numberOfExperimentals = 10;

  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      outputDirectory = outputOption->GetFunction( 0 )->GetParameter( 0 );
      outputDirectory += std::string( "/" );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      rootOutputFileName = outputOption->GetFunction( 0 )->GetParameter( 1 );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      numberOfControls = parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 2 ) );
      }
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
      {
      numberOfExperimentals = parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 3 ) );
      }
    }
  else
    {
    std::cout << "ERROR:  No output specified." << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Get the label mask.  If not specified, create one from the DTI atlas.
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer maskImageOption =
    parser->GetOption( "label-mask-image" );
  RealType lowerThresholdFunction = 0.2;
  if( maskImageOption && maskImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = maskImageOption->GetFunction( 0 )->GetName();
    typedef itk::ImageFileReader<MaskImageType> ReaderType;
    typename ReaderType::Pointer maskreader = ReaderType::New();
    maskreader->SetFileName( inputFile.c_str() );
    try
      {
      maskImage = maskreader->GetOutput();
      maskImage->Update();
      maskImage->DisconnectPipeline();
      }
    catch( ... )
      {
      lowerThresholdFunction = parser->Convert<RealType>( maskImageOption->GetFunction( 0 )->GetName() );
      }
    }
  if( maskImage.IsNull() )
    {
    std::cout << "Mask not read.  Creating mask by thresholding "
             << "the FA of the DTI atlas at >= " << lowerThresholdFunction
             << "." << std::endl << std::endl;

    typename ImageType::Pointer faImage =
      AllocImage<ImageType>(inputAtlas, 0.0);

    itk::ImageRegionIterator<TensorImageType> ItA( inputAtlas,
                                                   inputAtlas->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> ItF( faImage,
                                             faImage->GetLargestPossibleRegion() );
    for( ItA.GoToBegin(), ItF.GoToBegin(); !ItA.IsAtEnd(); ++ItA, ++ItF )
      {
      ItF.Set( CalculateFractionalAnisotropy<TensorType>( ItA.Get() ) );
      }

    typedef itk::BinaryThresholdImageFilter<ImageType, MaskImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( faImage );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->SetLowerThreshold( lowerThresholdFunction );
    thresholder->SetUpperThreshold( 2.0 );

    maskImage = thresholder->GetOutput();
    maskImage->Update();
    maskImage->DisconnectPipeline();
    }

  //
  // Get label information for pathology option
  //
  typedef itk::LabelGeometryImageFilter<MaskImageType, ImageType>
    LabelGeometryFilterType;
  typename LabelGeometryFilterType::Pointer labelGeometry =
    LabelGeometryFilterType::New();
  labelGeometry->SetInput( maskImage );
  labelGeometry->CalculatePixelIndicesOff();
  labelGeometry->CalculateOrientedBoundingBoxOff();
  labelGeometry->CalculateOrientedLabelRegionsOff();
  labelGeometry->Update();

  typename LabelGeometryFilterType::LabelsType labels =
    labelGeometry->GetLabels();
  std::sort( labels.begin(), labels.end() );

  unsigned int totalMaskVolume = 0;
  for( unsigned int n = 0; n < labels.size(); n++ )
    {
    totalMaskVolume += static_cast<unsigned int>( labelGeometry->GetVolume(
                                                    labels[n] ) );
    }

  // Fill in default values per label:
  //    column 1:  percentage change in longitudinal eigenvalue
  //    column 2:  percentage change in average of transverse eigenvalue(s)
  //    column 3:  percentage of affected voxels
  itk::Array2D<RealType> pathologyParameters( labels.size(), 3 );
  for( unsigned int i = 0; i < labels.size(); i++ )
    {
    pathologyParameters(i, 0) = 0.0;
    pathologyParameters(i, 1) = 0.0;
    pathologyParameters(i, 2) = 0.0;
    }

  /**
   * labels
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer pathologyOption =
    parser->GetOption( "pathology" );
  if( pathologyOption && pathologyOption->GetNumberOfFunctions() )
    {
    if( pathologyOption->GetNumberOfFunctions() == 1 &&
        ( pathologyOption->GetFunction( 0 )->GetName() ).empty() )
      {
      float pathologyDeltaEig1 = 0.0;
      float pathologyDeltaEig2_Eig3 = 0.0;
      float percentageVoxels = 0.0;

      if( pathologyOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        pathologyDeltaEig1 = parser->Convert<float>(
            pathologyOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      if( pathologyOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        pathologyDeltaEig2_Eig3 = parser->Convert<float>(
            pathologyOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      if( pathologyOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        percentageVoxels = parser->Convert<float>(
            pathologyOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      for( unsigned int n = 0; n < labels.size(); n++ )
        {
        RealType percentage = percentageVoxels;
        if( percentage > itk::NumericTraits<RealType>::OneValue() )
          {
          percentage /= static_cast<RealType>( labelGeometry->GetVolume( labels[n] ) );
          }

        pathologyParameters(n, 0) = pathologyDeltaEig1;
        pathologyParameters(n, 1) = pathologyDeltaEig2_Eig3;
        pathologyParameters(n, 2) = percentage;
        }
      }
    else
      {
      for( unsigned int n = 0; n < pathologyOption->GetNumberOfFunctions(); n++ )
        {
        auto whichClass = parser->Convert<LabelType>( pathologyOption->GetFunction( n )->GetName() );

        auto it = std::find( labels.begin(), labels.end(), whichClass );

        if( it == labels.end() )
          {
          continue;
          }

        float pathologyDeltaEig1 = 0.1;
        float pathologyDeltaEig2_Eig3 = 0.1;
        float percentageVoxels = 1.0;

        if( pathologyOption->GetFunction( n )->GetNumberOfParameters() > 0 )
          {
          pathologyDeltaEig1 = parser->Convert<float>( pathologyOption->GetFunction( n )->GetParameter( 0 ) );
          }
        if( pathologyOption->GetFunction( n )->GetNumberOfParameters() > 1 )
          {
          pathologyDeltaEig2_Eig3 = parser->Convert<float>( pathologyOption->GetFunction( n )->GetParameter( 1 ) );
          }
        if( pathologyOption->GetFunction( n )->GetNumberOfParameters() > 2 )
          {
          percentageVoxels = parser->Convert<float>( pathologyOption->GetFunction( n )->GetParameter( 2 ) );
          }

        RealType percentage = percentageVoxels;
        if( percentage > itk::NumericTraits<RealType>::OneValue() )
          {
          percentage /= static_cast<RealType>( labelGeometry->GetVolume( *it ) );
          }
        pathologyParameters(it - labels.begin(), 0) = pathologyDeltaEig1;
        pathologyParameters(it - labels.begin(), 1) = pathologyDeltaEig2_Eig3;
        pathologyParameters(it - labels.begin(), 2) = percentage;
        }
      }
    }

  /**
   * Perform PCA decomposition on input registered population
   */
  bool applyISV = false;

  typename MatrixType::InternalMatrixType ISV(1, 1);
  typename itk::ants::CommandLineParser::OptionType::Pointer populationOption = parser->GetOption(
      "registered-population" );
  if( populationOption && populationOption->GetNumberOfFunctions() )
    {
    std::cout << "--- Modeling intersubject variability ---" << std::endl
             << std::endl;

    std::vector<std::string> imageNames;

    std::string filename = populationOption->GetFunction( 0 )->GetName();

    std::string imageFile;

    std::fstream str( filename.c_str() );
    while( str >> imageFile )
      {
      imageNames.push_back( imageFile );
      }

    str.close();

    MatrixType M;
    MatrixType Mt;
    MatrixType E;
    MatrixType Lambda;
    M.SetSize( imageNames.size(), 2 * totalMaskVolume );
    M.Fill( 0 );
    for( unsigned int k = 0; k < imageNames.size(); k++ )
      {
      std::cout << "Processing " << imageNames[k] << " (" << k + 1 << " of "
               << imageNames.size() << ")." << std::endl;
      typename TensorReaderType::Pointer tensorReader = TensorReaderType::New();
      tensorReader->SetFileName( imageNames[k].c_str() );
      tensorReader->Update();

      unsigned int count = 0;

      itk::ImageRegionIterator<TensorImageType> It( tensorReader->GetOutput(),
                                                    tensorReader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
                                                   maskImage->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItM )
        {
        if( ItM.Get() != 0 )
          {
          TensorType tensor = It.Get();

          typename TensorType::EigenValuesArrayType eigenvalues;
          typename TensorType::EigenVectorsMatrixType eigenvectors;
          tensor.ComputeEigenAnalysis( eigenvalues, eigenvectors );

          if( eigenvalues[0] < 0 )
            {
            eigenvalues[0] = eigenvalues[1];
            }
          if( ImageDimension == 3 && eigenvalues[2] < 0 )
            {
            eigenvalues[2] = eigenvalues[1];
            }

          if( ImageDimension == 2 )
            {
            M(k, count) = eigenvalues[1];
            M(k, totalMaskVolume + count) = eigenvalues[0];
            }
          else
            {
            M(k, count) = eigenvalues[2];
            M(k, totalMaskVolume + count) = static_cast<RealType>( 0.5 )
              * ( eigenvalues[0] + eigenvalues[1] );
            }
          ++count;
          }
        }
      }

    std::cout << std::endl;
    // Now that the matrix M has been calculated, we need to subtract out
    // the longitudinal mean before performing PCA
    for( unsigned int i = 0; i < M.Cols(); i++ )
      {
      RealType columnAverage = 0.0;
      for( unsigned int j = 0; j < M.Rows(); j++ )
        {
        columnAverage += M(j, i);
        }
      columnAverage /= static_cast<RealType>( M.Rows() );
      for( unsigned int j = 0; j < M.Rows(); j++ )
        {
        M(j, i) -= columnAverage;
        }
      }
    // Perform PCA decomposition

    MatrixType MMt = M;
    MMt *= M.GetTranspose();
    decomposer->EvaluateSymmetricEigenDecomposition( MMt, Lambda, E );

    ISV = ( M.GetTranspose() * E.GetVnlMatrix() )
      / std::sqrt( static_cast<float>( imageNames.size() ) );

    applyISV = true;
    }

  //
  // Get DWI parameters
  //
  typename ImageType::Pointer b0Image = nullptr;
  unsigned int                       numberOfDirections = 0;
  std::vector<vnl_vector<RealType> > directions;
  std::vector<RealType>              bvalues;

  vnl_vector<RealType> direction( ImageDimension );

  // Add a B0 value at direction 0
  direction.fill( 0.0 );
  directions.push_back( direction );
  bvalues.push_back( 0 );

  typename itk::ants::CommandLineParser::OptionType::Pointer dwiOption =
    parser->GetOption( "dwi-parameters" );
  if( dwiOption && dwiOption->GetNumberOfFunctions() &&
      dwiOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName( dwiOption->GetFunction( 0 )->GetParameter( 0 ) );
    reader2->Update();
    b0Image = reader2->GetOutput();
    b0Image->DisconnectPipeline();

    std::string  directionsFileName = dwiOption->GetFunction( 0 )->GetParameter( 1 );
    std::fstream str( directionsFileName.c_str() );

    if( dwiOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      bvalues.push_back(
        parser->Convert<RealType>( dwiOption->GetFunction( 0 )->GetParameter( 2 ) ) );

      str >> numberOfDirections;
      }

    RealType x = 0.0;

    unsigned int count = 0;
    while( str >> x )
      {
      direction[count % ImageDimension] = x;
      ++count;
      if( count % ImageDimension == 0 )
        {
        directions.push_back( direction );
        if( dwiOption->GetFunction( 0 )->GetNumberOfParameters() < 3 )
          {
          str >> x;
          bvalues.push_back( x );
          }
        else
          {
          bvalues.push_back( bvalues[1] );
          }
        }
      }

    if( dwiOption->GetFunction( 0 )->GetNumberOfParameters() < 3 )
      {
      if( bvalues.size() != directions.size() )
        {
        std::cout << "ERROR:  Number of bvalues does not match the number of directions."
                 << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      if( numberOfDirections != directions.size() - 1 )
        {
        std::cout << "ERROR:  Number of directions does not match the data file."
                 << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  else
    {
    std::cout << "ERROR:  No DWI parameters specified." << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Get Rician noise parameter
  //
  RealType noiseSigma = 0;

  typename itk::ants::CommandLineParser::OptionType::Pointer noiseOption =
    parser->GetOption( "noise-sigma" );
  if( noiseOption && noiseOption->GetNumberOfFunctions() )
    {
    noiseSigma = parser->Convert<RealType>( noiseOption->GetFunction()->GetName() );
    }

  //
  // Create the simulated diffusion-weighted images.  For each image, we
  // perform the following steps:
  //   1. Copy the atlas
  //   2. Construct new DTI
  //     2a. Apply pathology (only for the experimentals).
  //     2b. Introduce subject intervariability
  //   3. For each direction, write new DWI
  //     3a. Use DTI from 2 to reconstruct DWI in current direction
  //     3b. Add Rician noise
  //
  itksys::SystemTools::MakeDirectory( outputDirectory.c_str() );

  itk::Array2D<RealType> meanFAandMD( labels.size(), 5 );
  meanFAandMD.Fill( 0.0 );
  for( unsigned n = 0; n <= numberOfControls + numberOfExperimentals; n++ )
    {
    if( n == 0 )
      {
      std::cout << "--- Calculating regional average FA and MD values (original and "
               << "pathology + intersubject variability) ---" << std::endl << std::endl;
      }
    else if( n <= numberOfControls )
      {
      if( n == 1 )
        {
        std::cout << std::endl << "--- Writing images ---" << std::endl << std::endl;
        }
      std::cout << "Writing control " << n
               << " (of " << numberOfControls << ") DWI images." << std::endl;
      }
    else
      {
      std::cout << "Writing experimental " << n - numberOfControls
               << " (of " << numberOfExperimentals << ") DWI images." << std::endl;
      }

    // copy atlas
    typedef itk::ImageDuplicator<TensorImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( inputAtlas );
    duplicator->Update();

    typename TensorImageType::Pointer dti = duplicator->GetOutput();

    // If we are to apply intersubject variability, we calculate random
    // projection.
    vnl_vector<RealType> eigenISVProjection( 1 );
    if( applyISV )
      {
      vnl_vector<RealType> R( ISV.cols() );
      for(float & d : R)
        {
        d = randomizer->GetNormalVariate( 0.0, 1.0 );
        }
      eigenISVProjection = ISV * R;
      }

    //
    // Iterate through the image to apply pathology and inter-subject variability
    //
    unsigned long count = 0;

    itk::ImageRegionIterator<TensorImageType> It( dti,
                                                  dti->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
                                                 maskImage->GetLargestPossibleRegion() );
    for( It.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItM )
      {
      LabelType  label = ItM.Get();
      TensorType tensor = It.Get();

      typename TensorType::EigenValuesArrayType eigenvalues;
      typename TensorType::EigenVectorsMatrixType eigenvectors;
      tensor.ComputeEigenAnalysis( eigenvalues, eigenvectors );

      if( eigenvalues[0] < 0 )
        {
        eigenvalues[0] = eigenvalues[1];
        }
      if( ImageDimension == 3 && eigenvalues[2] < 0 )
        {
        eigenvalues[2] = eigenvalues[1];
        }

      auto it = std::find( labels.begin(),
                                                             labels.end(), label );
      if( it == labels.end() )
        {
        std::cout << "ERROR:  unknown label." << std::endl;
        }
      unsigned int labelIndex = it - labels.begin();

      typename TensorType::EigenValuesArrayType newEigenvalues;

      //
      // Only apply pathology to a certain fraction of the voxels for a
      // particular label.  We "throw the dice" to determine whether or not
      // to apply to the current voxel.
      //
      RealType pathologyLongitudinalChange = 0.0;
      RealType pathologyTransverseChange = 0.0;
      if( ( n == 0 || n > numberOfControls ) && static_cast<RealType>(
            randomizer->GetUniformVariate( 0.0, 1.0 ) ) <= pathologyParameters(labelIndex, 2) )
        {
        pathologyLongitudinalChange = pathologyParameters(labelIndex, 0);
        pathologyTransverseChange = pathologyParameters(labelIndex, 1);
        }

      //
      // Apply intersubject variability
      //
      RealType isvLongitudinalProjection = 0.0;
      RealType isvTransverseProjection = 0.0;
      if( label != 0 && applyISV )
        {
        isvLongitudinalProjection = eigenISVProjection(count);
        isvTransverseProjection = eigenISVProjection(totalMaskVolume + count);
        count++;
        }

      //
      // Reconstruct the tensor
      //
      if( ImageDimension == 2 )
        {
        newEigenvalues[1] = eigenvalues[1]
          + eigenvalues[1] * pathologyLongitudinalChange
          + isvLongitudinalProjection;
        newEigenvalues[0] = eigenvalues[0] * ( itk::NumericTraits<RealType>::OneValue() + eigenvalues[0] )
          * pathologyTransverseChange + isvTransverseProjection;
        if( newEigenvalues[0] >= newEigenvalues[1] )
          {
          newEigenvalues[0] = newEigenvalues[1] - static_cast<RealType>( 1.0e-6 );
          }
        }
      else
        {
        newEigenvalues[2] = eigenvalues[2]
          + eigenvalues[2] * pathologyLongitudinalChange
          + isvLongitudinalProjection;
        RealType eigenAverage = static_cast<RealType>( 0.5 ) * ( eigenvalues[1] + eigenvalues[0] );
        newEigenvalues[1] = ( static_cast<RealType>( 2.0 ) * eigenAverage
                              * ( itk::NumericTraits<RealType>::OneValue() + pathologyTransverseChange )
                              + isvTransverseProjection ) / ( eigenvalues[0] / eigenvalues[1] + itk::NumericTraits<RealType>::OneValue() );
        if( newEigenvalues[1] >= newEigenvalues[2] )
          {
          newEigenvalues[1] = newEigenvalues[2] - static_cast<RealType>( 1.0e-6 );
          }
        newEigenvalues[0] = ( eigenvalues[0] / eigenvalues[1] )
          * newEigenvalues[1];
        }
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        if( std::isnan( newEigenvalues[d] ) )
          {
          newEigenvalues[d] = 0.0;
          }
        }

      if( newEigenvalues[0] < 0 )
        {
        newEigenvalues[0] = newEigenvalues[1];
        }
      if( ImageDimension == 3 && newEigenvalues[2] < itk::NumericTraits<RealType>::ZeroValue() )
        {
        newEigenvalues[2] = newEigenvalues[1];
        }

      typename TensorType::MatrixType eigenvalueMatrix;
      eigenvalueMatrix.Fill( 0.0 );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        eigenvalueMatrix(d, d) = newEigenvalues[d];
        }

      typename TensorType::MatrixType D( eigenvectors.GetTranspose() );
      D *= eigenvalueMatrix;
      D *= eigenvectors;

      TensorType newTensor;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for( unsigned int j = i; j < ImageDimension; j++ )
          {
          newTensor(i, j) = D(i, j);
          }
        }

      if( label != 0 && n == 0 )
        {
        meanFAandMD(labelIndex, 0) += static_cast<RealType>(
          CalculateFractionalAnisotropy<TensorType>( tensor ) );
        meanFAandMD(labelIndex, 1) += static_cast<RealType>(
          CalculateMeanDiffusivity<TensorType>( tensor ) );
        meanFAandMD(labelIndex, 2) += static_cast<RealType>(
          CalculateFractionalAnisotropy<TensorType>( newTensor ) );
        meanFAandMD(labelIndex, 3) += static_cast<RealType>(
          CalculateMeanDiffusivity<TensorType>( newTensor ) );
        meanFAandMD(labelIndex, 4)++;
        }
      else if( n != 0 )
        {
        It.Set( newTensor );
        }
      }

    if( n == 0 )
      {
      std::cout << "   " << std::left << std::setw( 7 ) << "Region"
               << std::left << std::setw( 15 ) << "FA (original)"
               << std::left << std::setw( 15 ) << "FA (path+isv)"
               << std::left << std::setw( 15 ) << "FA (% change)"
               << std::left << std::setw( 15 ) << "MD (original)"
               << std::left << std::setw( 15 ) << "MD (path+isv)"
               << std::left << std::setw( 15 ) << "MD (% change)"
               << std::endl;
      for( unsigned int l = 1; l < labels.size(); l++ )
        {
        std::cout << "   " << std::left << std::setw( 7 ) << labels[l]
                 << std::left << std::setw( 15 ) << meanFAandMD(l, 0) / meanFAandMD(l, 4)
                 << std::left << std::setw( 15 ) << meanFAandMD(l, 2) / meanFAandMD(l, 4)
                 << std::left << std::setw( 15 )
                 << ( meanFAandMD(l, 2) - meanFAandMD(l, 0) ) / meanFAandMD(l, 0)
                 << std::left << std::setw( 15 ) << meanFAandMD(l, 1) / meanFAandMD(l, 4)
                 << std::left << std::setw( 15 ) << meanFAandMD(l, 3) / meanFAandMD(l, 4)
                 << std::left << std::setw( 15 )
                 << ( meanFAandMD(l, 3) - meanFAandMD(l, 1) ) / meanFAandMD(l, 1)
                 << std::endl;
        }
      }
    else
      {
      std::string which;
      if( n <= numberOfControls )
        {
        which = std::string( "Control" );
        }
      else
        {
        which = std::string( "Experimental" );
        }

      std::stringstream istream;
      if( n <= numberOfControls )
        {
        istream << n;
        }
      else
        {
        istream << ( n - numberOfControls );
        }
      std::string dwiSeriesFileNames = outputDirectory + which + istream.str()
        + rootOutputFileName + std::string( "Direction%03d.nii.gz" );

      itk::NumericSeriesFileNames::Pointer dwiFileNamesCreator =
        itk::NumericSeriesFileNames::New();
      dwiFileNamesCreator->SetStartIndex( 0 );
      dwiFileNamesCreator->SetEndIndex( directions.size() - 1 );
      dwiFileNamesCreator->SetSeriesFormat( dwiSeriesFileNames.c_str() );
      std::vector<std::string> dwiImageNames = dwiFileNamesCreator->GetFileNames();
      for( unsigned int d = 0; d < directions.size(); d++ )
        {
        vnl_vector<RealType> bk = directions[d];
        RealType             bvalue = bvalues[d];

        std::cout << "  Applying direction " << d << " (of "
                 << directions.size() - 1 << "): [" << bk << "]"
                 << ", bvalue = " << bvalue << std::endl;

        typename ImageType::Pointer dwi =
          AllocImage<ImageType>(dti, 0);

        itk::ImageRegionConstIterator<ImageType> ItB( b0Image,
                                                      b0Image->GetLargestPossibleRegion() );
        itk::ImageRegionIterator<ImageType> ItD( dwi,
                                                 dwi->GetLargestPossibleRegion() );
        for( It.GoToBegin(), ItB.GoToBegin(), ItD.GoToBegin(); !It.IsAtEnd();
             ++It, ++ItB, ++ItD )
          {
          TensorType tensor = It.Get();
          for( unsigned int i = 0; i < tensor.GetNumberOfComponents(); i++ )
            {
            if( std::isnan( tensor[i] ) )
              {
              tensor[i] = 0.0;
              }
            }

          vnl_matrix<RealType> D(ImageDimension, ImageDimension);
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            for( unsigned int j = 0; j < ImageDimension; j++ )
              {
              D(i, j) = tensor(i, j);
              }
            }

          vnl_vector<RealType> bkD = bk * D;

          RealType signal = ItB.Get() * std::exp( -bvalue * inner_product( bkD, bk ) );

          // Add Rician noise
          RealType realNoise = 0.0;
          RealType imagNoise = 0.0;
          if( noiseSigma > itk::NumericTraits<RealType>::ZeroValue() )
            {
            realNoise = randomizer->GetNormalVariate( 0.0,
                                                      itk::Math::sqr( noiseSigma ) );
            imagNoise = randomizer->GetNormalVariate( 0.0,
                                                      itk::Math::sqr( noiseSigma ) );
            }
          RealType realSignal = signal + realNoise;
          RealType imagSignal = imagNoise;

          std::complex<RealType> noisySignal( realSignal, imagSignal );

          RealType finalSignal = std::sqrt( std::norm( noisySignal ) );

          if( signal <= ItB.Get() )
            {
            ItD.Set( finalSignal );
            }
          }
        typedef itk::ImageFileWriter<ImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( dwiImageNames[d].c_str() );
        writer->SetInput( dwi );
        writer->Update();
        }
      }
    }
  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, the program tries to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "image-dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "A diffusion tensor atlas image is required input for " )
      + std::string( "creating the cohort. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dti-atlas" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "inputDTIAtlasFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "A mask image can be specified which determines the region(s). " )
      + std::string( "to which the simulated pathology operations are applied. " )
      + std::string( "See also the option '--pathology'.  If no mask is specified " )
      + std::string( "one is created by thresholding the atlas FA map at 0.2 unless  " )
      + std::string( "a lower threshold is specified." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "label-mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskImageFileName" );
    option->SetUsageOption( 1, "lowerThresholdFunction" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "This parameter characterizes the Rician noise in the original DWI" )
      + std::string( "images.  Van Hecke uses the noise-estimation method of Sijbers et " )
      + std::string( "al. \"Automatic estimation of the noise variance from the " )
      + std::string( "histogram of a magnetic resonance image\", Phys. Med. Biol. " )
      + std::string( "52:1335-1348, 2007." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "noise-sigma" );
    option->SetShortName( 'n' );
    option->SetUsageOption( 0, "<noiseSigma=18>" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The user can specify the simulated pathology in a given " )
      + std::string( "area using a label mask. If no label is prepended to " )
      + std::string( "parameters, the specified parameters are applied to all labels." )
      + std::string( "Pathology is simulated by changing the eigenvalues. Typically " )
      + std::string( "this involves a decrease in the largest eigenvalue and an " )
      + std::string( "increase in the average of the remaining eigenvalues. " )
      + std::string( "Change is specified as a percentage of the current eigenvalues. " )
      + std::string( "However, care is taken " )
      + std::string( "to ensure that diffusion direction does not change. " )
      + std::string( "Additionally, one can specify the number of voxels affected " )
      + std::string( "in each region or one can specify the percentage of voxels " )
      + std::string( "affected.  Default is to change all voxels.  Note that the " )
      + std::string( "percentages must be specified in the range [0,1]. For " )
      + std::string( "dimension=3 where the average transverse diffusion eigenvalues " )
      + std::string( "are altered, this change is propagated to the distinct eigenvalues " )
      + std::string( "by forcing the ratio to be the same before the change. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "pathology" );
    option->SetUsageOption(
      0,
      "label[<percentageChangeEig1=-0.05>,<percentageChangeAvgEig2andEig3=0.05>,<numberOfVoxels=all or percentageOfvoxels>]" );
    option->SetShortName( 'p' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "This option specifies the parameters of the output " )
      + std::string( "diffusion-weighted images including the directions and " )
      + std::string( "b-values.  The directions are specified using a direction " )
      + std::string( "file which has as its first line the number of directions." )
      + std::string( "Each successive three lines contains the x, y, and z " )
      + std::string( "directions, respectively, and a single b-value. " )
      + std::string( "Note that several direction files of this format are " )
      + std::string( "distributed with the Camino DTI toolkit " )
      + std::string( "(http://web4.cs.ucl.ac.uk/research/medic/camino/pmwiki/pmwiki.php).  " )
      + std::string( "Alternatively, one can specify a scheme file where each direction " )
      + std::string( "is specified followed by a b-value for that direction, i.e. " )
      + std::string( "<x1> <y1> <z1> <bvalue1> ... <xN><yN><zN><bvalueN>." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dwi-parameters" );
    option->SetShortName( 'w' );
    option->SetUsageOption( 0, "[B0Image,directionFile,bvalue]" );
    option->SetUsageOption( 1, "[B0Image,schemeFile]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "If one wants to introduce inter-subject variability" )
      + std::string( "a registered DTI population to the DTI atlas is " )
      + std::string( "required.  This variability is modeled by a PCA " )
      + std::string( "decomposition on a combination of the first eigenvalue " )
      + std::string( "image and the average of the second and third eigenvalues." )
      + std::string( "The registered image file names are specified using " )
      + std::string( "a text file " )
      + std::string( "where each line is the name of an individual DTI." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "registered-population" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "textFileWithFileNames.txt" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of a set of diffusion-weighted images " )
      + std::string( "for each subject.  Each file name is prepended with the " )
      + std::string( "word 'Control' or 'Experimental'.  The number of control " )
      + std::string( "and experimental subjects can be also be specified on the " )
      + std::string( "command line.  Default is 10 for each group." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0,
                            "[outputDirectory,fileNameSeriesRootName,<numberOfControls=10>,<numberOfExperimentals=10>]" );
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

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int CreateDTICohort( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CreateDTICohort" );

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
  argv[argc] = nullptr;
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

  // antscout->set_stream( out_stream );

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "CreateDTICohort implements the work of Van Hecke et al. (" )
    + std::string( "On the construction of a ground truth framework for " )
    + std::string( "evaluating voxl-based diffusion tensor MRI analysis " )
    + std::string( "methods, Neuroimage 46:692-707, 2009) to create " )
    + std::string( "simulated DTI data sets.  The only " )
    + std::string( "difference is that all registrations (both for the input " )
    + std::string( "population and for the output population) are assumed to " )
    + std::string( "take place outside of this program." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc < 2 || parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction()->GetName() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "dti-atlas" );
    if( imageOption && imageOption->GetNumberOfFunctions() )
      {
      if( imageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        filename = imageOption->GetFunction( 0 )->GetParameter( 0 );
        }
      else
        {
        filename = imageOption->GetFunction( 0 )->GetName();
        }
      }
    else
      {
      std::cout << "No input atlas was specified.  Specify a dti atlas"
               << " with the -a option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  std::cout << std::endl << "Creating DTI cohort for "
           << dimension << "-dimensional images." << std::endl << std::endl;

  switch( dimension )
    {
    case 2:
      {
      return CreateDTICohort<2>( parser );
      }
      break;
    case 3:
      {
      return CreateDTICohort<3>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
