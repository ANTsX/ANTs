


#include "antsUtilities.h"
#include <algorithm>

#include <string>

#include <cmath>
#include <ctime>
#include "itkVTKPolyDataWriter.h"
#include "vtkSTLWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMesh.h"
#include "itkSphereMeshSource.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkExtractEdges.h"
#include "vtkPolyDataReader.h"

#include "BinaryImageToMeshFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkDataSetMapper.h"
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkGraphicsFactory.h>
#include "ReadWriteData.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vtkDelaunay2D.h"
#include "vtkFloatArray.h"
#include <vtkSmartPointer.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include "vtkVolume16Reader.h"
#include "vtkImageReader2.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkOutlineFilter.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkPolyData.h"
#include "vtkPolyVertex.h"
#include "vtkPointData.h"
#include "vtkExtractEdges.h"
#include "vtkPolyDataNormals.h"
#include "vtkMarchingCubes.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkDecimatePro.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataConnectivityFilter.h"

#include "vtkSmoothPolyDataFilter.h"
#include "vtkSTLWriter.h"
#include "itkSurfaceImageCurvature.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>


namespace ants
{

void
Display(vtkUnstructuredGrid * vtkgrid, std::string offscreen, bool secondwin = false, bool delinter = true)
{

  vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
  normalGenerator->SetInputData(vtkgrid);
  normalGenerator->ComputePointNormalsOn();
  normalGenerator->ComputeCellNormalsOff();
  normalGenerator->Update();

  vtkSmartPointer<vtkGraphicsFactory> graphics_factory = vtkSmartPointer<vtkGraphicsFactory>::New();
  if (offscreen.length() > 4)
    graphics_factory->SetOffScreenOnlyMode(1);
  graphics_factory->SetUseMesaClasses(1);

  vtkRenderer *     ren1 = vtkRenderer::New();
  vtkRenderer *     ren2 = vtkRenderer::New();
  vtkRenderWindow * renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);
  if (secondwin)
  {
    renWin->AddRenderer(ren2);
  }
  vtkRenderWindowInteractor * inter = vtkRenderWindowInteractor::New();
  inter->SetRenderWindow(renWin);
  vtkCallbackCommand * cbc = vtkCallbackCommand::New();
  ren1->AddObserver(vtkCommand::KeyPressEvent, cbc);

  vtkDataSetMapper * mapper = vtkDataSetMapper::New();
  mapper->SetInputData(normalGenerator->GetOutput());

  // Create a lookup table to map cell data to colors
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  int                             tableSize = std::max(10, 10);
  lut->SetNumberOfTableValues(tableSize);
  lut->Build();
  lut->SetTableRange(0, 1);
  lut->SetTableValue(0, 1, 1, 1, 1);                // Black
  lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
  lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
  lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
  lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
  lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
  lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
  lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
  lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
  lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
  mapper->SetLookupTable(lut);
  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(0, 1);
  //  mapper->SetScalarModeToUsePointData();
  //  mapper->SetColorModeToMapScalars();

  vtkActor * actor = vtkActor::New();
  actor->SetMapper(mapper);
  vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(mapper->GetLookupTable());
  scalarBar->SetTitle("F(x)");
  scalarBar->SetNumberOfLabels(4);
  // Create a lookup table to share between the mapper and the scalarbar
  vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
  hueLut->SetTableRange(0, 1);
  hueLut->SetHueRange(0, 1);
  hueLut->SetSaturationRange(1, 1);
  hueLut->SetValueRange(1, 1);
  hueLut->Build();

  //  mapper->SetLookupTable( hueLut );
  scalarBar->SetLookupTable(lut);
  vtkDataSetMapper * mapper2 = vtkDataSetMapper::New();
  if (secondwin)
  {
    mapper2->SetInputData(vtkgrid);
  }
  if (secondwin)
  {
    mapper2->SetScalarRange(0, 255);
  }
  vtkActor * actor2 = vtkActor::New();
  if (secondwin)
  {
    actor2->SetMapper(mapper2);
  }
  if (secondwin)
  {
    ren1->SetViewport(0.0, 0.0, 0.5, 1.0);
    ren2->SetViewport(0.5, 0.0, 1.0, 1.0);
  }
  else
  {
    ren1->SetViewport(0.0, 0.0, 1.0, 1.0);
  }
  if (secondwin)
  {
    ren2->AddActor(actor2);
  }
  // add the actor and start the render loop
  // see http://www.vtk.org/doc/nightly/html/classvtkProperty.html
  actor->GetProperty()->SetInterpolationToFlat();
  actor->GetProperty()->SetInterpolationToGouraud();
  actor->GetProperty()->ShadingOff();
  ren1->AddActor(actor);
  ren1->SetBackground(1, 1, 1); // Background color
  ren1->AddActor2D(scalarBar);
  renWin->Render();

  if (offscreen.length() > 4)
  {
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renWin);
    windowToImageFilter->SetScale(4);
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(offscreen.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
  }
  if (offscreen.length() < 5)
    inter->Start();
  mapper->Delete();
  actor->Delete();
  ren1->Delete();
  mapper2->Delete();
  actor2->Delete();
  ren2->Delete();
  renWin->Delete();
  if (delinter)
  {
    inter->Delete();
  }
}

float
ComputeGenus(vtkPolyData * pd1)
{
  vtkExtractEdges * edgeex = vtkExtractEdges::New();

  edgeex->SetInputData(pd1);
  edgeex->Update();
  vtkPolyData * edg1 = edgeex->GetOutput();

  vtkIdType nedg = edg1->GetNumberOfCells();
  vtkIdType vers = pd1->GetNumberOfPoints();
  int       nfac = pd1->GetNumberOfPolys();

  float g = 0.5 * (2.0 - vers + nedg - nfac);
  std::cout << " Genus " << g << std::endl;

  std::cout << " face " << nfac << " edg " << nedg << " vert " << vers << std::endl;

  return g;
}


float
vtkComputeTopology(vtkPolyData * pd)
{
  // Marching cubes
  //    std::cout << " Marching Cubes ";
  //    vtkMarchingCubes *marchingCubes = vtkMarchingCubes::New();
  //    vtkContourFilter *marchingCubes = vtkContourFilter::New();
  //    vtkKitwareContourFilter *marchingCubes = vtkKitwareContourFilter::New();
  //    marchingCubes->SetInput((vtkDataSet*) vds);
  //    marchingCubes->SetValue(0, hithresh);
  //    int nc;
  //    std::cout << " Input #conts "; std::cin >> nc;
  //    marchingCubes->SetNumberOfContours(2);
  //    marchingCubes->SetComputeScalars(false);
  //    marchingCubes->SetComputeGradients(false);
  //    marchingCubes->SetComputeNormals(false);
  vtkPolyDataConnectivityFilter * con = vtkPolyDataConnectivityFilter::New();

  con->SetExtractionModeToLargestRegion();
  //    con->SetInput(marchingCubes->GetOutput());
  con->SetInputData(pd);
  con->Update();
  float g = ComputeGenus(con->GetOutput());
  return g;
#if 0
  int inputNumberOfPoints = con->GetOutput()->GetNumberOfPoints();
  int inputNumberOfPolys = con->GetOutput()->GetNumberOfPolys();

  vtkPolyDataConnectivityFilter *polyDataConnectivityFilter =
    vtkPolyDataConnectivityFilter::New();
  polyDataConnectivityFilter->SetInputData( con->GetOutput() );
  polyDataConnectivityFilter->SetExtractionModeToAllRegions();
  polyDataConnectivityFilter->SetExtractionModeToLargestRegion();
  polyDataConnectivityFilter->Update();

  int connectivityNumberOfExtractedRegions = polyDataConnectivityFilter->
    GetNumberOfExtractedRegions();

  polyDataConnectivityFilter->Delete();

  vtkExtractEdges *extractEdges = vtkExtractEdges::New();
  extractEdges->SetInputData( con->GetOutput() );
  extractEdges->Update();

  int extractNumberOfLines = extractEdges->GetOutput()->GetNumberOfLines();

  extractEdges->Delete();

  int EulerCharacteristic = inputNumberOfPoints - extractNumberOfLines
    + inputNumberOfPolys;

  double genus = 0.5 * ( 2 * connectivityNumberOfExtractedRegions
                         - EulerCharacteristic );

  std::cout << "EulerCharacteristic " << EulerCharacteristic << std::endl;
  std::cout << "genus " << genus << std::endl;

  return genus;
#endif
}

template <typename TImage>
void
GetValueMesh(typename TImage::Pointer image,
             typename TImage::Pointer image2,
             std::string              outfn,
             const char *             paramname,
             float                    scaledata,
             float                    aaParm,
             std::string              offscreen,
             unsigned int             inflate)
{
  //  std::cout << " parname " << std::string(paramname) << std::endl;
  using ImageType = TImage;

  using FilterType = BinaryImageToMeshFilter<ImageType>;
  typename FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput(image);
  fltMesh->SetAntiAliasMaxRMSError(aaParm); // to do nothing, set negative
  fltMesh->SetSmoothingIterations(0);
  fltMesh->SetDecimateFactor(0.0);
  fltMesh->Update();
  vtkPolyData * vtkmesh = fltMesh->GetMesh();
  // assign scalars to the original surface mesh
  //  std::string offsc=std::string("");
  //  Display((vtkUnstructuredGrid*)vtkmesh, offsc );
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetInputData(vtkmesh);
  smoother->SetNumberOfIterations(25);
  smoother->BoundarySmoothingOff();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(120.0);
  smoother->SetPassBand(.01);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->Update();
  vtkmesh = smoother->GetOutput();

  std::cout << " Genus " << vtkComputeTopology(vtkmesh) << std::endl;

  vtkPoints * vtkpoints = vtkmesh->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  float       mx = 0, mn = 9.e9, meank = 0;
  for (int i = 0; i < numPoints; i++)
  {
    typename ImageType::IndexType index;
    typename ImageType::PointType point;
    for (int j = 0; j < 3; j++)
    {
      point[j] = (vtkpoints->GetPoint(i)[j]);
    }
    image2->TransformPhysicalPointToIndex(point, index);
    float temp = image2->GetPixel(index);
    if (fabs(temp) > mx)
    {
      mx = fabs(temp);
    }
    if (fabs(temp) < mn && mn > 0)
    {
      mn = fabs(temp);
    }
    meank += fabs(temp);
  }
  std::cout << " max kap " << mx << " mn k " << mn << std::endl;
  meank /= numPoints;
  //  mx=1.3;
  //  mx=2.0;
  float localscaledata = scaledata;
  localscaledata = 1;
  // while (!done)
  {
    vtkFloatArray * param = vtkFloatArray::New();
    param->SetName(paramname);
    float       dif = (mx - mn) * localscaledata;
    const float mx2 = meank + dif;
    const float mn2 = meank - dif;
    dif = mx2 - mn2;
    for (int i = 0; i < (numPoints); i++)
    {
      typename ImageType::IndexType index;
      typename ImageType::PointType point;
      for (int j = 0; j < 3; j++)
      {
        point[j] = (vtkpoints->GetPoint(i)[j]);
      }
      image2->TransformPhysicalPointToIndex(point, index);
      float temp = image2->GetPixel(index);
      //    param->InsertNextValue(temp);
      //    float temp=surfk->CurvatureAtIndex(index);
      if (i % 1000 == 0)
      {
        std::cout << " kappa " << temp << std::endl;
      }
      // =fabs(manifoldIntegrator->GetGraphNode(i)->GetTotalCost());

      temp = std::fabs(temp);
      float vvv = (temp - mn2) * 255.0f / dif;
      vvv = (temp - mn) / dif;
      /*
      if (vvv > 128)
        {
      float dif=255-vvv;
      vvv = 128 - dif;
        }
      else
        {
      float dif=128-vvv;
      vvv = 128 + dif;
      }*/
      param->InsertNextValue(vvv);
    }
    vtkmesh->GetPointData()->SetScalars(param);
  }
  std::cout << " Now display to " << offscreen << std::endl;
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> inflater = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  inflater->SetInputData(vtkmesh);
  inflater->SetNumberOfIterations(inflate);
  inflater->BoundarySmoothingOn();
  inflater->FeatureEdgeSmoothingOff();
  inflater->SetFeatureAngle(180.0);
  inflater->SetEdgeAngle(180.0);
  inflater->SetPassBand(0.001); // smaller values increase smoothing
  inflater->NonManifoldSmoothingOn();
  inflater->NormalizeCoordinatesOff();
  if (inflate > 0)
  {
    inflater->Update();
    vtkmesh = inflater->GetOutput();
  }
  if (offscreen.length() > 2)
    Display((vtkUnstructuredGrid *)vtkmesh, offscreen);
  std::cout << " done with mesh map ";
  /*
  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInputData( vtkmesh );
  writer->SetFileName( outfn.c_str() );
  writer->Update();
  */
  vtkSTLWriter * writer = vtkSTLWriter::New();
  writer->SetInputData(vtkmesh);
  std::cout << " writing " << outfn << std::endl;
  writer->SetFileName(outfn.c_str());
  writer->Write();
  std::cout << " done writing ";
  inflater->Delete();
  smoother->Delete();
  std::cout << " done writing2 ";
}

template <typename TImage>
float
GetImageTopology(typename TImage::Pointer image)
{
  using ImageType = TImage;

  double aaParm = 0.024;
  using FilterType = BinaryImageToMeshFilter<ImageType>;
  typename FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput(image);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->SetAntiAliasMaxRMSError(-1000.0); // to do nothing
  fltMesh->Update();
  vtkPolyData * vtkmesh = fltMesh->GetMesh();
  // assign scalars to the original surface mesh
  //  Display((vtkUnstructuredGrid*)vtkmesh);

  float genus = vtkComputeTopology(vtkmesh);
  std::cout << " Genus " << genus << std::endl;

  return genus;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
GetMeshAndTopology(std::vector<std::string> args, std::ostream *)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "GetMeshAndTopology");

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

  // antscout->set_stream( out_stream );

  if (argc < 2)
  {
    std::cout
      << argv[0]
      << " binaryimage valueimage  out paramname ValueScale AntiaAliasParm=0.001 offscreen.png  inflation-interations "
      << std::endl;
    std::cout << " outputs vtk version of input image -- assumes object is defined by non-zero values " << std::endl;
    std::cout << " mesh is colored by the value of the image voxel " << std::endl;
    std::cout << " the AntiaAliasParm could cause topo problems but makes nicer meshes  " << std::endl;
    std::cout << " the offscreen param will render to screen if set to win, 0 means no rendering " << std::endl;
    std::cout << " ValueScale controls contrast in image appearance - lower increaseses contrast -- should be <= 1 "
              << std::endl;
    return EXIT_FAILURE;
  }

  // Define the dimension of the images
  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  // Declare the types of the output images
  using ImageType = itk::Image<PixelType, Dimension>;

  // Declare the type of the Mesh

  std::string outfn = std::string(argv[3]);

  ImageType::Pointer image;
  ReadImage<ImageType>(image, argv[1]);
  ImageType::Pointer image2;
  ReadImage<ImageType>(image2, argv[2]);

  ImageType::DirectionType fmat = image->GetDirection();
  fmat.SetIdentity();
  image->SetDirection(fmat);
  image2->SetDirection(fmat);

  // Save the mesh
  float        aaParm = 0.03;
  const char * paramname = std::string("f(x)").c_str();
  if (argc > 4)
  {
    paramname = (argv[4]);
  }
  float scaledata = 0.5;
  if (argc > 5)
  {
    scaledata = atof(argv[5]);
  }
  if (argc > 6)
  {
    aaParm = atof(argv[6]);
  }
  std::cout << "aaParm " << aaParm << std::endl;
  std::string offscreen = "win";
  if (argc > 7)
  {
    offscreen = std::string(argv[7]);
  }
  unsigned int inflate = 0;
  if (argc > 8)
  {
    inflate = std::stoi(argv[8]);
  }
  GetValueMesh<ImageType>(image, image2, outfn, paramname, scaledata, aaParm, offscreen, inflate);
  //  GetImageTopology<ImageType>(image);
  return EXIT_SUCCESS;
}
} // namespace ants
