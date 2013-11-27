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
#include "itkTransformToDisplacementFieldSource.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
template <unsigned int Dimension>
int antsApplyTransformsToPoints( itk::ants::CommandLineParser::Pointer & parser )
{
  typedef double                RealType;
  typedef double                PixelType;
  typedef vnl_matrix<PixelType> MatrixType;
  MatrixType points_out;
  MatrixType points_in;
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  typedef itk::CSVArray2DDataObject<double> DataFrameObjectType;
  typedef typename DataFrameObjectType::StringVectorType StringVectorType;
  StringVectorType colheadernames;
  //  bool input_points_are_indices = false;

  /**
   * Input object option - for now, we're limiting this to images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption( "input" );
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( inputOption && inputOption->GetNumberOfFunctions() > 0 )
    {
    std::cout << "Input csv file: " << inputOption->GetFunction( 0 )->GetName() << std::endl;
    std::string ext =
      itksys::SystemTools::GetFilenameExtension(  ( inputOption->GetFunction( 0 )->GetName() ).c_str()  );
    if( strcmp(ext.c_str(), ".csv") == 0 )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(  ( inputOption->GetFunction( 0 )->GetName() ).c_str()  );
      reader->SetFieldDelimiterCharacter( ',' );
      reader->SetStringDelimiterCharacter( '"' );
      reader->HasColumnHeadersOn();
      reader->HasRowHeadersOff();
      //    reader->UseStringDelimiterCharacterOff();
      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject& exp )
        {
        std::cout << "Exception caught!" << std::endl;
        std::cout << exp << std::endl;
        }
      DataFrameObjectType::Pointer dfo = reader->GetOutput();
      colheadernames = dfo->GetColumnHeaders();
      if ( colheadernames.size() < Dimension ) 
	{
	std::cout << "Input csv file must have column names such as x,y,z,t,label - where there are a minimum of N-Spatial-Dimensions names e.g. x,y in 2D." << std::endl;
	return EXIT_FAILURE;
	}
      points_in = dfo->GetMatrix();
      points_out.set_size( points_in.rows(),  points_in.cols() );
      }
    else
      {
      std::cout << "An input csv file is required." << std::endl;
      return EXIT_FAILURE;
      }

    if(  points_in.cols() < Dimension )
      {
      std::cout << "The number of columns in the input point set is fewer than " << Dimension << " Exiting."
               << std::endl;
      return EXIT_FAILURE;
      }

    if( outputOption && outputOption->GetNumberOfFunctions() > 0 )
      {
      if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) == 0 )
        {
        std::cout << "An input csv file is required." << std::endl;
        return EXIT_FAILURE;
        }
      }

    /**
     * Transform option
     */
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> MatrixOffsetTransformType;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
    typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> MatrixOffsetTransformType;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

    /**
     * Load an identity transform in case no transforms are loaded.
     */
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    typedef itk::MatrixOffsetTransformBase<double, 2, 2> MatrixOffsetTransformType2D;
    itk::TransformFactory<MatrixOffsetTransformType2D>::RegisterTransform();
    typedef itk::MatrixOffsetTransformBase<double, 3, 3> MatrixOffsetTransformType3D;
    itk::TransformFactory<MatrixOffsetTransformType3D>::RegisterTransform();

    typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;
    typename CompositeTransformType::InputPointType point_in;
    typename CompositeTransformType::OutputPointType point_out;
    typename itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption( "transform" );

    std::vector<bool> isDerivedTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<RealType, Dimension>( parser, transformOption, isDerivedTransform );
    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    for( unsigned int pointct = 0; pointct < points_in.rows(); pointct++ )
      {
      point_in.Fill( 0 );
      point_out.Fill( 0 );
      for( unsigned int p = 0; p < Dimension; p++ )
        {
        point_in[p] = points_in( pointct, p );
        }
      point_out = compositeTransform->TransformPoint( point_in );
      for( unsigned int p = 0; p < Dimension; p++ )
        {
        points_out( pointct, p ) = point_out[p];
        }
      for( unsigned int p = Dimension; p < points_in.cols(); p++ )
        {
	points_out( pointct, p ) = points_in( pointct, p );
        }
      }

    /**
     * output
     */
    if( outputOption && outputOption->GetNumberOfFunctions() > 0 )
      {
      std::string outputFileName = "";
      if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 &&
          parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) == 0 )
        {
        outputFileName = outputOption->GetFunction( 0 )->GetParameter( 0 );
        }
      else
        {
        outputFileName = outputOption->GetFunction( 0 )->GetName();
        }
      std::cout << "Output warped points to csv file: " << outputFileName << std::endl;
      StringVectorType ColumnHeaders = colheadernames;
      typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( outputFileName );
      writer->SetInput( &points_out );
      writer->SetColumnHeaders( ColumnHeaders );
      try
        {
        writer->Write();
        }
      catch( itk::ExceptionObject& exp )
        {
        std::cout << "Exception caught!" << std::endl;
        std::cout << exp << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  return EXIT_SUCCESS;
}

static void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "This option forces the points to be treated as a specified-" )
      + std::string( "dimensionality." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Currently, the only input supported is a csv file with " )
      + std::string( "columns including x,y (2D), x,y,z (3D) or x,y,z,t (4D) column headers." )
      + std::string( "The points should be defined in physical space." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "inputFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "One can output the warped points to a csv file.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "warpedOutputFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Several transform options are supported including all " )
      + std::string( "those defined in the ITK library in addition to " )
      + std::string( "a deformation field transform.  The ordering of " )
      + std::string( "the transformations follows the ordering specified " )
      + std::string( "on the command line.  An identity transform is pushed " )
      + std::string( "onto the transformation stack. Each new transform " )
      + std::string( "encountered on the command line is also pushed onto " )
      + std::string( "the transformation stack. Then, to warp the input object, " )
      + std::string( "each point comprising the input object is warped first " )
      + std::string( "according to the last transform pushed onto the stack " )
      + std::string( "followed by the second to last transform, etc. until " )
      + std::string( "the last transform encountered which is the identity " )
      + std::string( "transform. " )
      + std::string( "Also, it should be noted that the inverse transform can " )
      + std::string( "be accommodated with the usual caveat that such an inverse " )
      + std::string( "must be defined by the specified transform class " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "transform" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "transformFileName" );
    option->SetUsageOption( 1, "[transformFileName,useInverse]" );
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
int antsApplyTransformsToPoints( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsApplyTransformsToPoints" );
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

  // antscout->set_stream( out_stream );

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string examplestring = std::string( "reads in a csv file with the first D columns defining the spatial location where the spatial location is defined in physical coordinates.    the csv file should have a header row.   here is an example") + std::string("\n") + std::string("cat chicken-3.csv ") + std::string("x,y,z,t,label,comment")+std::string("\n")+std::string("82.5,116.5,0,0,1,this is the breast")+std::string("\n")+std::string("137.5,35.5,0,0,2,this is the beak")+std::string("\n")+std::string("antsApplyTransformsToPoints -d 2 -i chicken-3.csv -o test.csv -t [chicken3to4.mat ,1 ]")+std::string("\n")+std::string("cat test.csv ")+std::string("\n")+std::string("x,y,z,t,label,comment")+std::string("\n")+std::string("10.8945447481644,162.082675013049,0,0,1,nan")+std::string("\n")+std::string("7.5367085472988,52.099713111629,0,0,2,nan")+std::string("\n")+std::string("the nan appears in the last column until the ITK CSV I/O can handle mixed numeric / string types.  if your input is fully numeric, all is well.");
  std::string commandDescription =
    std::string( "antsApplyTransformsToPoints, applied to an input image, transforms it " )
    + std::string( "according to a reference image and a transform " )
    + std::string( "(or a set of transforms).  " ) + examplestring;

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || ( parser->GetOption( "help" ) &&
                    ( parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) ) ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' ) &&
           ( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

#if 0 // HACK This makes no sense here, filename is never used.
  //Perhaps the "input" option is not needed in this program
  //but is a copy/paste error from another program.
  // Read in the first intensity image to get the image dimension.
  std::string filename;
  itk::ants::CommandLineParser::OptionType::Pointer inputOption =
    parser->GetOption( "input" );
  if( inputOption && inputOption->GetNumberOfFunctions() > 0 )
    {
    if( inputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      filename = inputOption->GetFunction( 0 )->GetParameter( 0 );
      }
    else
      {
      filename = inputOption->GetFunction( 0 )->GetName();
      }
    }
  else
    {
    std::cout << "No csv file point set was specified." << std::endl;
    return EXIT_FAILURE;
    }
#endif

  unsigned int dimension = 3;
  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    std::cout << "No -d ( dimensionality ) option is specified.  Exiting." << std::endl;
    return EXIT_FAILURE;
    }

  switch( dimension )
    {
    case 2:
      {
      antsApplyTransformsToPoints<2>( parser );
      }
      break;
    case 3:
      {
      antsApplyTransformsToPoints<3>( parser );
      }
      break;
    case 4:
      {
      antsApplyTransformsToPoints<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
