#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkFEM.h"
#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMElement3DC0LinearTriangularLaplaceBeltrami.h"
#include "itkFEMElement3DC0LinearTriangularMembrane.h"
#include "itkFEMDiscConformalMap.h"

#include "vtkCallbackCommand.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include <vtkSmartPointer.h>
#include <vtkWindowedSincPolyDataFilter.h>

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

  // first find out if the user wants to inflate the mesh ...
  unsigned int inflate_iterations = 0;
  float        inflate_param = 0;
  typename itk::CommandLineParser::OptionType::Pointer infOption =
    parser->GetOption( "inflate" );
  if( infOption && infOption->GetNumberOfValues() > 0 )
    {
    if( infOption->GetNumberOfParameters() == 2 )
      {
      inflate_param = parser->Convert<float>(infOption->GetParameter( 0 ) );
      inflate_iterations = parser->Convert<unsigned int>(infOption->GetParameter( 1 ) );
      std::cout << " you will inflate before flattening with params " << inflate_param << " applied over  "
                << inflate_iterations << " iterations. " <<  std::endl;
      }
    else
      {
      std::cerr << " wrong params for inflation. ignoring. " << std::endl;
      std::cerr << "   " << infOption->GetDescription() << std::endl;
      return EXIT_FAILURE;
      }
    }

  typename itk::CommandLineParser::OptionType::Pointer displayOption = parser->GetOption( "display-mesh" );
  if( displayOption && displayOption->GetNumberOfValues() > 0 )
    {
    if( displayOption->GetNumberOfParameters() > 0 )
      {
      std::string dispm = displayOption->GetParameter( 0 );
      std::cout << " render " << dispm << std::endl;
      // read the vtk file ...
      vtkPolyDataReader *fltReader = vtkPolyDataReader::New();
      fltReader->SetFileName(dispm.c_str() );
      fltReader->Update();

      vtkRenderer*     ren1 = vtkRenderer::New();
      vtkRenderWindow* renWin = vtkRenderWindow::New();
      renWin->AddRenderer(ren1);
      vtkRenderWindowInteractor* inter = vtkRenderWindowInteractor::New();
      inter->SetRenderWindow(renWin);
      vtkCallbackCommand *cbc = vtkCallbackCommand::New();
      ren1->AddObserver(vtkCommand::KeyPressEvent, cbc);

      vtkDataSetMapper* mapper = vtkDataSetMapper::New();
      mapper->SetInput(fltReader->GetOutput() );
      mapper->SetScalarRange(0, 255);
      vtkActor* actor = vtkActor::New();
      actor->SetMapper(mapper);
      ren1->SetViewport(0.0, 0.0, 1.0, 1.0);
      ren1->AddActor(actor);

      renWin->Render();
      inter->Start();

      mapper->Delete();
      actor->Delete();
      ren1->Delete();
      renWin->Delete();
      inter->Delete();
      return 0;
      }
    }

  exit(1);
  /**
   * Initialization
   */
  typename itk::CommandLineParser::OptionType::Pointer initializationOption =
    parser->GetOption( "input-mesh" );
  if( initializationOption
      && initializationOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << "Incorrect input mesh specification." << std::endl;
    std::cerr << "   " << initializationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    /**
     * intensity images
     */
    typename itk::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "input-mesh" );
    if(  imageOption->GetNumberOfValues() > 0 )
      {
/** inflation
  std::cout << " read " << std::string(filename) << " write " << outfn << std::endl;
  vtkPolyDataReader *fltReader = vtkPolyDataReader::New();
  fltReader->SetFileName(filename);
  fltReader->Update();
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetInput(fltReader->GetOutput());
  smoother->SetNumberOfIterations( (int) inflate_iterations );
  smoother->BoundarySmoothingOn();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(180.0);
  smoother->SetEdgeAngle(180.0);
  smoother->SetPassBand( inflate_param ); // smaller values increase smoothing
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
      if( outputOption->GetNumberOfParameters() > 0 )
        {
        //        WriteImage<ImageType>(img,outputOption->GetParameter( 0 ) ).c_str()  );
        for( unsigned int p = 0; p < outputOption->GetNumberOfParameters(); p++ )
          {
          if( p == 0 )
            {
            vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
            writer->SetInput(flattener->m_DiskSurfaceMesh);
            std::string outnm = outputOption->GetParameter( p );
            std::cout << " writing " << outnm << std::endl;
            writer->SetFileName(outnm.c_str() );
            writer->SetFileTypeToBinary();
            if( flattener->m_DiskSurfaceMesh )
              {
              writer->Update();
              }
            }
          if( p == 1 )
            {
            vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
            writer->SetInput(flattener->m_ExtractedSurfaceMesh);
            std::string outnm = outputOption->GetParameter( 1 );
            std::cout << " writing " << outnm << std::endl;
            writer->SetFileName(outnm.c_str() );
            writer->SetFileTypeToBinary();
            if( flattener->m_DiskSurfaceMesh )
              {
              writer->Update();
              }
            }
          }
        }
      }
    }

  bool paramws = parser->template Convert<bool>( parser->GetOption( "param-while-searching" )->GetValue() );
  flattener->SetParamWhileSearching(paramws);

  std::string canonicaldomain = parser->GetOption( "canonical-domain" )->GetValue();
  ConvertToLowerCase( canonicaldomain );
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
      std::string(
        "Two mesh images are specified as input - 1. defines the label mesh.  2. defines the feature mesh.  we put the 2nd mesh's values into the flat space." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input-mesh" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "[InputMesh1.vtk,<InputMesh2.vtk>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Display the mesh." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "display-mesh" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "[InputMesh1.vtk]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Inflation --- two params : \n 1. BandPass (smaller increases smoothing) -- e.g. 0.001. \n " )
      + std::string( "2. number of iterations --- higher increases smoothing. " );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "inflate" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "[<InverseSmoothingFactor=0.001>,<iterations=150>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string(
        "The output consists of one (or two) meshes ... 1. the flattened mesh with features mapped.  2. the extracted extrinisic mesh. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "[MyFlatVTKImage.vtk,<MyOptionalExtrinsicVTKImage.vtk>]" );
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
