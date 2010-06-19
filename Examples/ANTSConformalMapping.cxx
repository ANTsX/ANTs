#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "itkFEM.h"
#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMElement3DC0LinearTriangularLaplaceBeltrami.h"
#include "itkFEMElement3DC0LinearTriangularMembrane.h"
#include "itkFEMDiscConformalMap.h"

#include <string>
#include <algorithm>
#include <vector>

template <class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate()
  {
  };
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter =
      dynamic_cast<const TFilter *>( object );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    std::cout << "Iteration " << filter->GetElapsedIterations()
              << " (of " << filter->GetMaximumNumberOfIterations() << "): ";
    std::cout << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
  }
};

/** inflation
  std::cout << " read " << std::string(filename) << " write " << outfn << std::endl;
  vtkPolyDataReader *fltReader = vtkPolyDataReader::New();
  fltReader->SetFileName(filename);
  fltReader->Update();
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetInput(fltReader->GetOutput());
  smoother->SetNumberOfIterations( (int) param );
  smoother->BoundarySmoothingOn();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(180.0);
  smoother->SetEdgeAngle(180.0);
  smoother->SetPassBand(1.e-3); // smaller values increase smoothing
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOff();
  smoother->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(	smoother->GetOutput() );
  writer->SetFileName(outfn.c_str());
  writer->SetFileTypeToBinary();
  writer->Update();
  std::cout << " done writing ";
*/

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), (int (*)(int) )tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
}

template <unsigned int ImageDimension>
int ANTSConformalMapping( itk::CommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  // we define the options in the InitializeCommandLineOptions function
  // and then use them here ...
  typedef vtkPolyData                                   MeshType;
  typedef itk::FEMDiscConformalMap<MeshType, ImageType> ParamType;
  typename ParamType::Pointer flattener = ParamType::New();
  flattener->SetDebug(false);
  flattener->SetSigma(1);
  //  flattener->SetSurfaceMesh(vtkmesh);

  /**
   * Initialization
   */
  typename itk::CommandLineParser::OptionType::Pointer initializationOption =
    parser->GetOption( "initialization" );
  if( initializationOption
      && initializationOption->GetNumberOfParameters() < 1 )
    {
    std::cerr << "Incorrect initialization option specification." << std::endl;
    std::cerr << "   " << initializationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    /**
     * intensity images
     */
    typename itk::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "intensity-image" );
    if( imageOption && imageOption->GetNumberOfValues() > 0 )
      {
      for( unsigned int n = 0; n < imageOption->GetNumberOfValues(); n++ )
        {
        typedef itk::ImageFileReader<ImageType> ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();
        if( imageOption->GetNumberOfParameters( n ) > 0 )
          {
          reader->SetFileName( imageOption->GetParameter( n, 0 ) );
          }
        else
          {
          reader->SetFileName( imageOption->GetValue( n ) );
          }
        reader->Update();
        }
      }
    else
      {
      std::cerr << "No input images were specified.  Specify an input image"
                << " with the -a option" << std::endl;
      }

    /**
     * output
     */
    typename itk::CommandLineParser::OptionType::Pointer outputOption =
      parser->GetOption( "output" );
    if( outputOption && outputOption->GetNumberOfValues() > 0 )
      {
      if( outputOption->GetNumberOfParameters() == 0 )
        {
        // WriteImage<ImageType>(img,outputOption->GetValue() ).c_str() );
        }
      if( outputOption->GetNumberOfParameters() > 0 )
        {
        //        WriteImage<ImageType>(img,outputOption->GetParameter( 0 ) ).c_str()  );
        }
      }
    }

  bool paramws = parser->template Convert<bool>( parser->GetOption( "param-while-searching" )->GetValue() );
  flattener->SetParamWhileSearching(paramws);

  std::string canonicaldomain = parser->GetOption( "canonical-domain" )->GetValue();
  std::cout << " you will map to a " << canonicaldomain << std::endl;
  if( canonicaldomain == std::string("circle") )
    {
    flattener->SetMapToCircle();
    }
  else if( canonicaldomain == std::string("square")  )
    {
    flattener->SetMapToSquare();
    }
  else
    {
    std::cout << " that domain is not an option " << std::endl;  return 1;
    }

  std::string boundaryparam = parser->GetOption( "boundary-param" )->GetValue();
  // do stuff -- but not implemented yet

  unsigned int labeltoflatten = parser->template Convert<unsigned int>(
      parser->GetOption( "label-to-flatten" )->GetValue() );
  std::cout << " you will flatten " << labeltoflatten << std::endl;
  flattener->SetLabelToFlatten(labeltoflatten);

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::CommandLineParser *parser )
{
  typedef itk::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "One or more scalar images is specified as input " );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "intensity-image" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "[intensityImage,<ScalarWeight>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of ...");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "[someImage,<some_Other_Optional_Images>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Map to a canonical domain : pass square or circle ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "canonical-domain" );
    option->SetShortName( 'c' );
    option->SetUsageOption( 0, "[domain]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string(
        "Parameterize the boundary while searching for the boundary (a bit slow and not guaranteed to be doable). " )
      + std::string( "If false, we try to parameterize after the searching is done. " )
      + std::string( "This option is meaningless if you pass the boundary in as an option. " );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "param-while-searching" );
    option->SetShortName( 'p' );
    option->SetUsageOption( 0, "[true/false]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The name of a boundary parameterization file (not implemented)." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "boundary-param" );
    option->SetShortName( 'b' );
    option->SetUsageOption( 0, "[filename.vtk]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cout << "Usage: " << argv[0]
              << " imageDimension args" << std::endl;
    exit( 1 );
    }

  itk::CommandLineParser::Pointer parser = itk::CommandLineParser::New();
  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "A tool for conformal mapping to various canonical coordinate systems: disc, square " )
    + std::string( " operates on 3D vtk triangulated meshes.")
    + std::string(
      " Open problems include computation of, consistent orientation of and parameterization of the boundary-condition defining loop.   Should we use curve matching ?  Knot points?  Min distortion?  " );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>(
        parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    exit( EXIT_FAILURE );
    }
  else if( parser->Convert<bool>(
             parser->GetOption( 'h' )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    exit( EXIT_FAILURE );
    }

  ANTSConformalMapping<3>( parser );
}
