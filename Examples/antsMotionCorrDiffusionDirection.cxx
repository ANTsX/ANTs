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
#include "itkImageRegistrationMethodv4.h"
#include "itkSyNImageRegistrationMethod.h"
#include "itkDisplacementFieldTransform.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkImageToHistogramFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"

#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransform.h"
#include "itkExtractImageFilter.h"

#include "itkBSplineTransformParametersAdaptor.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingVelocityFieldTransformParametersAdaptor.h"

#include "itkGradientDescentOptimizerv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkQuasiNewtonOptimizerv4.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMacro.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkResampleImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"

#include <sstream>

namespace ants
{

template <typename T>
inline std::string ants_moco_to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template <typename T>
struct ants_moco_index_cmp
  {
  ants_moco_index_cmp(const T _arr) : arr(_arr)
  {
  }

  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }

  const T arr;
  };



// Transform traits to generalize the rigid transform
//
template <unsigned int ImageDimension>
class RigidTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<double, ImageDimension> TransformType;
};

template <>
class RigidTransformTraits<2>
{
public:
  typedef itk::Euler2DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<3>
{
public:
  // typedef itk::VersorRigid3DTransform<double> TransformType;
  // typedef itk::QuaternionRigidTransform<double>  TransformType;
  typedef itk::Euler3DTransform<double> TransformType;
};

template <unsigned int ImageDimension>
class SimilarityTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<double, ImageDimension> TransformType;
};

template <>
class SimilarityTransformTraits<2>
{
public:
  typedef itk::Similarity2DTransform<double> TransformType;
};

template <>
class SimilarityTransformTraits<3>
{
public:
  typedef itk::Similarity3DTransform<double> TransformType;
};


int ants_motion_directions( itk::ants::CommandLineParser *parser )
{

  constexpr unsigned int ImageDimension = 3;

  typedef double                                    RealType;
  typedef itk::Image<RealType, ImageDimension>      FixedImageType;
  typedef itk::ImageFileReader<FixedImageType>      ImageReaderType;
  typedef vnl_matrix<RealType>                      vMatrix;
  vMatrix param_values;
  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;

  std::vector<CompositeTransformType::Pointer> CompositeTransformVector;

  typedef itk::Euler3DTransform<RealType>                   RigidTransformType;
  typedef itk::AffineTransform<RealType, ImageDimension>    AffineTransformType;

  typedef itk::ants::CommandLineParser ParserType;
  typedef ParserType::OptionType       OptionType;

  typedef double                                        ParameterValueType;
  typedef itk::CSVArray2DFileReader<ParameterValueType> MocoReaderType;
  typedef MocoReaderType::Array2DDataObjectType         MocoDataArrayType;
  typedef itk::Array2D<ParameterValueType>              DirectionArrayType;

  std::string outputName = "";
  std::string mocoName = "";
  std::string physicalName = "";


  OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    outputName = outputOption->GetFunction(0)->GetName();
    std::cout << "Output: " << outputName << std::endl;
    }
  else
    {
    std::cerr << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
    }

  OptionType::Pointer mocoOption = parser->GetOption( "moco" );
  if( mocoOption && mocoOption->GetNumberOfFunctions() )
    {
    mocoName = mocoOption->GetFunction(0)->GetName();
    std::cout << "Moco file: " << mocoName << std::endl;
    }
  else
    {
    std::cerr << "Motion parameter file not specified" << std::endl;
    return EXIT_FAILURE;
    }

  OptionType::Pointer physicalOption = parser->GetOption( "physical" );
  if( physicalOption && physicalOption->GetNumberOfFunctions() )
    {
    physicalName = physicalOption->GetFunction(0)->GetName();
    std::cout << "Physical space image: " << physicalName << std::endl;
    }
  else
    {
    std::cerr << "Physical space image not specified" << std::endl;
    return EXIT_FAILURE;
    }


  OptionType::Pointer schemeOption = parser->GetOption( "scheme" );
  OptionType::Pointer bvecOption = parser->GetOption( "bvec" );

  if ( bvecOption->GetNumberOfFunctions() && schemeOption->GetNumberOfFunctions() )
    {
    std::cerr << "Must use either scheme or bvec input, not both" << std::endl;
    return EXIT_FAILURE;
    }

  DirectionArrayType directionArray;

  if( schemeOption && schemeOption->GetNumberOfFunctions() )
    {
    std::string schemeName = schemeOption->GetFunction(0)->GetName();
    //std::cout << "Scheme file: " << schemeName << std::endl;

    std::cout << "Scheme input format still a work in progress" << std::endl;
    return EXIT_FAILURE;

    MocoReaderType::Pointer schemeReader = MocoReaderType::New();
    schemeReader->SetFileName( schemeName.c_str() ) ;
    schemeReader->SetFieldDelimiterCharacter( ' ' );
    schemeReader->HasColumnHeadersOff();
    schemeReader->HasRowHeadersOff();
    schemeReader->Parse();

    // Manually skip first 2 lines of data array
    MocoDataArrayType::MatrixType schemeMatrix = schemeReader->GetOutput()->GetMatrix();
    //std::cout << "scheme data array size = " << schemeMatrix.rows() << " x " << schemeMatrix.cols() << std::endl;

    directionArray.SetSize( schemeMatrix.rows()-2, schemeMatrix.cols() );

    for ( unsigned int i=2; i < schemeMatrix.rows(); i++ )
      {
        for ( unsigned int j=0; j < 3; j++ )
          {
            directionArray(i-2,j) = schemeMatrix(i-2,j);
          }
      }

    }

  bool transposeArray = true;
  if( bvecOption && bvecOption->GetNumberOfFunctions() )
    {
    std::string bvecName = bvecOption->GetFunction(0)->GetName();
    //std::cout << "bvec file: " << bvecName << std::endl;

    MocoReaderType::Pointer bvecReader = MocoReaderType::New();
    bvecReader->SetFileName( bvecName.c_str() ) ;
    bvecReader->SetFieldDelimiterCharacter( ' ' );
    bvecReader->HasColumnHeadersOff();
    bvecReader->HasRowHeadersOff();
    bvecReader->Update();
    MocoDataArrayType::MatrixType bvecMatrix = bvecReader->GetOutput()->GetMatrix();

    //std::cout << "BVEC data array size = " << bvecMatrix.rows() << " x " << bvecMatrix.cols() << std::endl;
    if ( bvecMatrix.cols() == 3 )
      {
      transposeArray = false;
      std::cout << "Column based format" << std::endl;
      }

    // Transpose the array
    if ( transposeArray )
      {
      directionArray.SetSize( bvecMatrix.cols(), bvecMatrix.rows() );
      }
    else
      {
      directionArray.SetSize( bvecMatrix.rows(), bvecMatrix.cols() );
      }

    if ( transposeArray )
      {
      for ( unsigned int i=0; i < bvecMatrix.cols(); i++ )
        {
        for ( unsigned int j=0; j < bvecMatrix.rows(); j++ )
          {
          directionArray(i,j) = bvecMatrix(j,i);
          }
        }
      }
    else
      {
      for ( unsigned int i=0; i < bvecMatrix.cols(); i++ )
        {
        for ( unsigned int j=0; j < bvecMatrix.rows(); j++ )
          {
          directionArray(j,i) = bvecMatrix(j,i);
          }
        }
      }
    }

  std::cout << "Read direction data of size: " << directionArray.rows() << " x "
              << directionArray.cols() << std::endl;

  // itkImageFileReader will set direction to identity if the image being read has more dimensions than the class template
  // eg if you pass a 4D image file name to a ReaderType whose dimension is 3
  //
  // Therefore check reference image is 3D, and fail if not
  //
 itk::ImageIOBase::Pointer imageIO =
   itk::ImageIOFactory::CreateImageIO(physicalName.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
 imageIO->SetFileName(physicalName.c_str() );
 try
   {
     imageIO->ReadImageInformation();
   }
 catch( ... )
   {
     std::cout << "Can't read reference image " << physicalName << std::endl;
     return EXIT_FAILURE;
   }
 if (imageIO->GetNumberOfDimensions() != ImageDimension)
   {
     std::cout << "Reference image must be 3D " << std::endl;
     return EXIT_FAILURE;
   }

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(physicalName.c_str());
  imageReader->Update();

  DirectionArrayType outputDirectionArray;
  outputDirectionArray.SetSize( directionArray.rows(), directionArray.cols() );

  MocoReaderType::Pointer mocoReader = MocoReaderType::New();
  mocoReader->SetFileName( mocoName.c_str() );
  mocoReader->SetFieldDelimiterCharacter( ',' );
  mocoReader->HasColumnHeadersOn();
  mocoReader->HasRowHeadersOff();
  mocoReader->Update();

  MocoDataArrayType::Pointer mocoDataArray = mocoReader->GetOutput();
  std::cout << "Read motion correction data of size: " << mocoDataArray->GetMatrix().rows() << " x "
            << mocoDataArray->GetMatrix().cols() << std::endl;

 // [-0.231374, -0.0543334, 0.677384]

  unsigned int nTransformParams = mocoDataArray->GetMatrix().cols() - 2;
  //std::cout << "# Transform parameters = " << nTransformParams << std::endl;

  AffineTransformType::Pointer toPhysical = AffineTransformType::New();
  toPhysical->SetIdentity();
  AffineTransformType::Pointer toIndex = AffineTransformType::New();
  toIndex->SetIdentity();

  AffineTransformType::MatrixType toPhysicalMatrix = toPhysical->GetMatrix();
  for ( unsigned int i=0; i < ImageDimension; i++)
    {
    for ( unsigned int j=0; j < ImageDimension; j++)
      {
      toPhysicalMatrix(i,j) = imageReader->GetOutput()->GetDirection()(i,j);
      }
    }
  toPhysical->SetMatrix(toPhysicalMatrix);
  toPhysical->GetInverse(toIndex);

  for ( unsigned int i=0; i<directionArray.rows(); i++ )
    {

    AffineTransformType::Pointer affineTransform = AffineTransformType::New();
    AffineTransformType::Pointer directionTransform = AffineTransformType::New();
    AffineTransformType::ParametersType params;

    if (nTransformParams == 6) {
      RigidTransformType::Pointer rigid = RigidTransformType::New();
      RigidTransformType::ParametersType rParams;

      rParams.SetSize( nTransformParams );
      for ( unsigned int t=0; t<nTransformParams; t++ )
        {
          rParams[t] = mocoDataArray->GetMatrix()(i,t+2);
        }
      rigid->SetParameters( rParams );
      affineTransform->SetMatrix( rigid->GetMatrix() );
      affineTransform->SetTranslation( rigid->GetTranslation() );
    }
    else if (nTransformParams == 12) {
      params.SetSize( nTransformParams );
      for ( unsigned int t=0; t<nTransformParams; t++ )
        {
          params[t] = mocoDataArray->GetMatrix()(i,t+2);
        }
      affineTransform->SetParameters( params );
      affineTransform->GetInverse( directionTransform );
    }
    else {
      // Not rigid (6 params) or affine (12), something is wrong
      return EXIT_FAILURE;
    }

    //std::cout << affineTransform->GetTranslation() << std::endl;

    AffineTransformType::InputVectorType dir;
    AffineTransformType::OutputVectorType rotatedDir;



    for ( unsigned int j=0; j<directionArray.cols(); j++ )
      {
      dir[j] = directionArray(i,j);
      }
    if ( dir.GetNorm() > 0 )
      {
      dir.Normalize();
      rotatedDir = toPhysical->TransformVector( dir );
      rotatedDir = affineTransform->TransformVector( rotatedDir );
      rotatedDir.Normalize();
      rotatedDir = toIndex->TransformVector( rotatedDir );

      //rotatedDir = affineTransform->TransformVector( dir );
      //rotatedDir.Normalize();

      }
    else
      {
      for ( unsigned int j=0; j<directionArray.cols(); j++)
        {
        rotatedDir[j] = dir[j];
        }
      }

    for ( unsigned int j=0; j<directionArray.cols(); j++ )
      {
      outputDirectionArray(i,j) = rotatedDir[j];
      }

    std::cout << dir << " -> " << rotatedDir << std::endl;
    }


  // Write new directions to output file
  if ( outputName.find( ".bvec" ) != std::string::npos )
    {
    //std::cout << "Writing bvec file " << outputName << std::endl;
    std::ofstream outfile( outputName.c_str() );


    if ( transposeArray )
      {
      for ( unsigned int i=0; i<outputDirectionArray.cols(); i++ )
        {
        for ( unsigned int j=0; j<outputDirectionArray.rows(); j++)
          {
          outfile << outputDirectionArray(j,i) << " ";
          }
        outfile << std::endl;
        }
      }
    else
      {
      for ( unsigned int i=0; i<outputDirectionArray.rows(); i++ )
        {
        for ( unsigned int j=0; j<outputDirectionArray.cols(); j++)
          {
          outfile << outputDirectionArray(i,j) << " ";
          }
        outfile << std::endl;
        }
      }


    outfile.close();

    }
  else
    {
    std::cout << "WARNING - only bvec output currently supported" << std::endl;
    }

  return EXIT_SUCCESS;
}

void antsMotionCorrDiffusionDirectionInitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description = std::string( "Camino scheme file specificy acquisition parameters." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "scheme" );
    option->SetShortName( 's' );
    option->SetDescription( description );
    option->SetUsageOption( 0, "dwi.scheme" );
    parser->AddOption( option );
    }

    {
    std::string description =       std::string( "bvec image specifying diffusion directions.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "bvec" );
    option->SetShortName( 'b' );
    option->SetDescription( description );
    option->SetUsageOption( 0, "dwi.bvec" );
    parser->AddOption( option );
    }

    {
    std::string description =       std::string( "3D image in dwi space");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "physical" );
    option->SetShortName( 'p' );
    option->SetDescription( description );
    option->SetUsageOption( 0, "ref.nii.gz" );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "motion correction parameters from antsMotionCorr." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "moco" );
    option->SetShortName( 'm' );
    option->SetDescription( description );
    option->SetUsageOption( 0, "MOCOparams.csv" );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the output file for corrected directions." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetDescription( description );
    option->SetUsageOption( 0, "corrected.scheme" );
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
int antsMotionCorrDiffusionDirection( std::vector<std::string> args, std::ostream * /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsMotionCorrDiffusionDirection" );

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

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription = std::string( "antsMotionCorrDiffusionDirection" );
  parser->SetCommandDescription( commandDescription );
  antsMotionCorrDiffusionDirectionInitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

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

  std::cout << std::endl << "Running " << argv[0] << "  for 3-dimensional images." << std::endl
            << std::endl;

  return ants_motion_directions( parser );

}

} // namespace ants
