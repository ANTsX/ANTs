
#include "antsUtilities.h"
#include <algorithm>

#include "antsCommandLineParser.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

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

namespace ants
{
template <typename TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate(){};

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event)
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) != typeid(itk::IterationEvent))
    {
      return;
    }

    std::cout << "Iteration " << filter->GetElapsedIterations() << " (of " << filter->GetMaximumNumberOfIterations()
              << "): ";
    std::cout << filter->GetCurrentConvergenceMeasurement() << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
  }
};

void
InitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
    std::string description =
      std::string("Two mesh images are specified as input - 1. defines the label mesh. must have a scalar attached to "
                  "vertices named 'Label'. 2. defines the feature mesh.  scalar name is 'Feature'.  we put the 2nd "
                  "mesh's values into the flat space.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input-mesh");
    option->SetShortName('i');
    option->SetUsageOption(0, "[InputMesh1.vtk,<InputMesh2.vtk>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Display the mesh.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("display-mesh");
    option->SetShortName('d');
    option->SetUsageOption(0, "[InputMesh1.vtk]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Inflation --- two params : \n 1. BandPass (smaller increases smoothing) -- e.g. 0.001. \n ") +
      std::string("2. number of iterations --- higher increases smoothing. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("inflate");
    option->SetShortName('f');
    option->SetUsageOption(0, "[<InverseSmoothingFactor=0.001>,<iterations=150>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string(
        "SegmentationCost --- 3 params : \n 1. Float-Max-Cost : controls the size of the output region.  \n ") +
      std::string("2. Float-Weight for distance cost =  edge_length*W_d \n ") +
      std::string("3. Float-Weight for label cost = H(fabs( desiredLabel - localLabel ))*W_l*MaxCost \n where H is the "
                  "heaviside function.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("segmentation-cost");
    option->SetShortName('s');
    option->SetUsageOption(0, "e.g. [40,1,0]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("The output consists of one (or more) meshes ... 1. the flattened mesh with "
                                          "features mapped.  2. the extracted extrinisic mesh.  3. inflated mesh. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "[MyFlatMesh.vtk,<MyOptionalExtrinsicMesh.vtk>,<MyOptionalInflatedMesh.vtk>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Map to a canonical domain : pass square or circle ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("canonical-domain");
    option->SetShortName('c');
    option->SetUsageOption(0, "[domain]");
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

  {
    std::string description =
      std::string(
        "Parameterize the boundary while searching for the boundary (a bit slow and not guaranteed to be doable). ") +
      std::string("If false, we try to parameterize after the searching is done. ") +
      std::string("This option is meaningless if you pass the boundary in as an option. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("param-while-searching");
    option->SetShortName('p');
    option->SetUsageOption(0, "0 / 1 ");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Which label (unsigned integer) to flatten. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("label-to-flatten");
    option->SetShortName('l');
    option->SetUsageOption(0, "1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("The name of a boundary parameterization file (not implemented).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("boundary-param");
    option->SetShortName('b');
    option->SetUsageOption(0, "[filename.vtk]");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

template <unsigned int ImageDimension>
int
ANTSConformalMapping(itk::ants::CommandLineParser * parser)
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  // we define the options in the InitializeCommandLineOptions function
  // and then use them here ...
  typedef vtkPolyData                                   MeshType;
  vtkSmartPointer<MeshType>                             labelmesh = nullptr;
  vtkSmartPointer<MeshType>                             featuremesh = nullptr;
  vtkSmartPointer<MeshType>                             inflatedmesh = nullptr;
  typedef itk::FEMDiscConformalMap<MeshType, ImageType> ParamType;
  typename ParamType::Pointer                           flattener = ParamType::New();
  flattener->SetDebug(false);
  flattener->SetSigma(1);
  //  flattener->SetSurfaceMesh(vtkmesh);

  // first find out if the user wants to inflate the mesh ...
  unsigned int                                               inflate_iterations = 0;
  float                                                      inflate_param = 0;
  typename itk::ants::CommandLineParser::OptionType::Pointer infOption = parser->GetOption("inflate");
  if (infOption && infOption->GetNumberOfFunctions())
  {
    if (infOption->GetFunction(0)->GetNumberOfParameters() == 2)
    {
      inflate_param = parser->Convert<float>(infOption->GetParameter(0));
      inflate_iterations = parser->Convert<unsigned int>(infOption->GetParameter(1));
      std::cout << " you will inflate before flattening with params " << inflate_param << " applied over  "
                << inflate_iterations << " iterations. " << std::endl;
    }
    else
    {
      std::cerr << " wrong params for inflation. ignoring. " << std::endl;
      std::cerr << "   " << infOption->GetDescription() << std::endl;
      return EXIT_FAILURE;
    }
  }

  float                                                      maxCost = 40, distCostW = 1, labelCostW = 0;
  typename itk::ants::CommandLineParser::OptionType::Pointer costOption = parser->GetOption("segmentation-cost");
  if (costOption && costOption->GetNumberOfFunctions())
  {
    if (costOption->GetFunction(0)->GetNumberOfParameters() == 3)
    {
      maxCost = parser->Convert<float>(costOption->GetParameter(0));
      distCostW = parser->Convert<float>(costOption->GetParameter(1));
      labelCostW = parser->Convert<float>(costOption->GetParameter(2));
    }
    else
    {
      std::cerr << " wrong params for cost weights. " << std::endl;
      std::cerr << "   " << costOption->GetDescription() << std::endl;
      return EXIT_FAILURE;
    }
  }

  typename itk::ants::CommandLineParser::OptionType::Pointer displayOption = parser->GetOption("display-mesh");
  if (displayOption && displayOption->GetNumberOfFunctions())
  {
    if (displayOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      std::string dispm = displayOption->GetParameter(0);
      std::cout << " render " << dispm << std::endl;
      // read the vtk file ...
      vtkPolyDataReader * fltReader = vtkPolyDataReader::New();
      fltReader->SetFileName(dispm.c_str());
      fltReader->Update();

      vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
      normalGenerator->SetInput(fltReader->GetOutput());
      normalGenerator->Update();

      vtkRenderer *     ren1 = vtkRenderer::New();
      vtkRenderWindow * renWin = vtkRenderWindow::New();
      renWin->AddRenderer(ren1);
      vtkRenderWindowInteractor * inter = vtkRenderWindowInteractor::New();
      inter->SetRenderWindow(renWin);
      vtkCallbackCommand * cbc = vtkCallbackCommand::New();
      ren1->AddObserver(vtkCommand::KeyPressEvent, cbc);

      vtkDataSetMapper * mapper = vtkDataSetMapper::New();
      mapper->SetInput(normalGenerator->GetOutput());
      mapper->SetScalarRange(0, 255);
      vtkActor * actor = vtkActor::New();
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
      fltReader->Delete();
      return 0;
    }
  }

  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");
  /**
   * Initialization
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inOption = parser->GetOption("input-mesh");
  if (inOption && inOption->GetFunction(0)->GetNumberOfParameters() == 2)
  {
    std::string                        innm = inOption->GetParameter(0);
    vtkSmartPointer<vtkPolyDataReader> labReader = vtkSmartPointer<vtkPolyDataReader>::New();
    labReader->SetFileName(innm.c_str());
    labReader->Update();
    labelmesh = vtkSmartPointer<vtkPolyData>(labReader->GetOutput());

    vtkDataArray * labels = labelmesh->GetPointData()->GetArray("Label");
    if (!labels)
    {
      std::cerr << "  Cannot find vtk Array named 'Label' in " << innm << std::endl;
      std::cerr << "  This could cause problems " << std::endl;
      //      std::cout <<" exiting " << std::endl;
      // throw std::exception();
    }
    innm = inOption->GetParameter(1);
    vtkSmartPointer<vtkPolyDataReader> fltReader = vtkSmartPointer<vtkPolyDataReader>::New();
    fltReader->SetFileName(innm.c_str());
    fltReader->Update();
    featuremesh = vtkSmartPointer<vtkPolyData>(fltReader->GetOutput());
    vtkDataArray * feats = featuremesh->GetPointData()->GetArray("Feature");
    if (!feats)
    {
      std::cerr << "  Cannot find vtk Array named 'Feature' in " << innm << std::endl;
      std::cerr << " continuing " << std::endl;
    }

    /** inflation */
    if (inflate_iterations > 0)
    {
      vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
      smoother->SetInput(labelmesh);
      smoother->SetNumberOfIterations((int)inflate_iterations);
      smoother->BoundarySmoothingOn();
      smoother->FeatureEdgeSmoothingOff();
      smoother->SetFeatureAngle(180.0);
      smoother->SetEdgeAngle(180.0);
      smoother->SetPassBand(inflate_param); // smaller values increase smoothing
      smoother->NonManifoldSmoothingOn();
      smoother->NormalizeCoordinatesOff();
      smoother->Update();
      inflatedmesh = vtkSmartPointer<vtkPolyData>(smoother->GetOutput());
      std::cout << " done smoothing " << std::endl;
      flattener->SetSurfaceMesh(inflatedmesh);
      if (outputOption->GetFunction(0)->GetNumberOfParameters() > 0)
      {
        for (unsigned int p = 0; p < outputOption->GetFunction(0)->GetNumberOfParameters(); p++)
        {
          if (p == 2 && inflatedmesh)
          {
            vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
            writer->SetInput(inflatedmesh);
            std::string outnm = outputOption->GetParameter(2);
            std::cout << " writing " << outnm << std::endl;
            writer->SetFileName(outnm.c_str());
            writer->SetFileTypeToBinary();
            writer->Update();
          }
        }
      }
    }
    else
    {
      flattener->SetSurfaceMesh(labelmesh);
    }
    flattener->SetSurfaceFeatureMesh(featuremesh);
  }

  bool paramws = parser->template Convert<bool>(parser->GetOption("param-while-searching")->GetFunction());
  flattener->SetParamWhileSearching(paramws);

  unsigned int labeltoflatten =
    parser->template Convert<unsigned int>(parser->GetOption("label-to-flatten")->GetFunction());
  flattener->SetLabelToFlatten(labeltoflatten);
  std::string canonicaldomain = parser->GetOption("canonical-domain")->GetFunction();
  //    canonicaldomain=ConvertToLowerCase( canonicaldomain );
  std::cout << " you will map label " << labeltoflatten << " to a " << canonicaldomain << std::endl;
  if (canonicaldomain == std::string("circle"))
  {
    flattener->SetMapToCircle();
  }
  else if (canonicaldomain == std::string("square"))
  {
    flattener->SetMapToSquare();
  }
  else
  {
    std::cerr << " that domain is not an option -- exiting. " << std::endl;
    return EXIT_FAILURE;
  }

  // do stuff -- but not implemented yet
  //   flattener->SetDiscBoundaryList(nullptr);

  std::cout << " you will flatten " << labeltoflatten << ".  param while searching? " << paramws << std::endl;
  flattener->SetSigma(1);

  flattener->SetMaxCost(maxCost);
  flattener->SetDistanceCostWeight(distCostW);
  flattener->SetLabelCostWeight(labelCostW);
  std::cout << " MC " << maxCost << " DW " << distCostW << " LW " << labelCostW << std::endl;
  flattener->ExtractSurfaceDisc();
  std::cout << " begin conformal mapping ";
  flattener->ConformalMap();

  /**
   * output
   */
  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    if (outputOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      for (unsigned int p = 0; p < outputOption->GetFunction(0)->GetNumberOfParameters(); p++)
      {
        if (p == 0)
        {
          vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
          writer->SetInput(flattener->m_DiskSurfaceMesh);
          std::string outnm = outputOption->GetParameter(p);
          std::cout << " writing " << outnm << std::endl;
          writer->SetFileName(outnm.c_str());
          writer->SetFileTypeToBinary();
          if (flattener->m_DiskSurfaceMesh)
          {
            writer->Update();
          }
        }
        if (p == 1)
        {
          vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
          writer->SetInput(flattener->m_ExtractedSurfaceMesh);
          std::string outnm = outputOption->GetParameter(1);
          std::cout << " writing " << outnm << std::endl;
          writer->SetFileName(outnm.c_str());
          writer->SetFileTypeToBinary();
          if (flattener->m_ExtractedSurfaceMesh)
          {
            writer->Update();
          }
        }
        if (p == 2 && inflatedmesh)
        {
          vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
          writer->SetInput(inflatedmesh);
          std::string outnm = outputOption->GetParameter(2);
          std::cout << " writing " << outnm << std::endl;
          writer->SetFileName(outnm.c_str());
          writer->SetFileTypeToBinary();
          writer->Update();
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ANTSConformalMapping(std::vector<std::string> args, std::ostream * out_stream = nullptr)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ANTSConformalMapping");
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
  argv[argc] = 0;
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

  // antscout->set_stream( out_stream );

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();
  parser->SetCommand(argv[0]);

  std::string commandDescription =
    std::string("A tool for conformal mapping to various canonical coordinate systems: disc, square ") +
    std::string(" operates on 3D vtk triangulated meshes.") +
    std::string(" Open problems include computation of, consistent orientation of and parameterization of the "
                "boundary-condition defining loop.   Should we use curve matching ?  Knot points?  Min distortion?  ");

  parser->SetCommandDescription(commandDescription);
  InitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc < 2 || parser->Convert<bool>(parser->GetOption("help")->GetFunction()))
  {
    parser->PrintMenu(std::cout, 5, false);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (parser->Convert<bool>(parser->GetOption('h')->GetFunction()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  ANTSConformalMapping<3>(parser);
}
} // namespace ants
