



#include "antsUtilities.h"
#include <algorithm>

#include <string>

#include <math.h>
#include <time.h>

#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMesh.h"
#include "itkSphereMeshSource.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "vtkPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkExtractEdges.h"
#include "vtkPolyDataReader.h"

// #include "itkFEMConformalMap.h"
// #include "itkFEMDiscConformalMap.h"

#include "BinaryImageToMeshFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "ReadWriteImage.h"
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
// #include "vtkKitwareContourFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSTLWriter.h"
#include "vtkUnstructuredGridToPolyDataFilter.h"
#include "itkSurfaceImageCurvature.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"

namespace ants
{
template <class TImage>
typename TImage::Pointer BinaryThreshold(
  typename TImage::PixelType bkg,
  typename TImage::PixelType foreground,
  typename TImage::PixelType replaceval, typename TImage::Pointer input)
{
  typedef typename TImage::PixelType PixelType;
  // Begin Threshold Image
  typedef itk::BinaryThresholdImageFilter<TImage, TImage>
    InputThresholderType;
  typename InputThresholderType::Pointer inputThresholder =
    InputThresholderType::New();

  inputThresholder->SetInput( input );
  inputThresholder->SetInsideValue(  replaceval );
  int outval = 0;
  if( (float) replaceval == (float) -1 )
    {
    outval = 1;
    }
  inputThresholder->SetOutsideValue( outval );

  float low = bkg;
  float high = foreground;
  if( high < low )
    {
    high = 255;
    }
  inputThresholder->SetLowerThreshold( (PixelType) low );
  inputThresholder->SetUpperThreshold( (PixelType) high);
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

float ComputeGenus(vtkPolyData* pd1)
{
  vtkExtractEdges* edgeex = vtkExtractEdges::New();

  edgeex->SetInput(pd1);
  edgeex->Update();
  vtkPolyData* edg1 = edgeex->GetOutput();

  vtkIdType nedg = edg1->GetNumberOfCells();
  vtkIdType vers = pd1->GetNumberOfPoints();
  int       nfac = pd1->GetNumberOfPolys();

  float g = 0.5 * (2.0 - vers + nedg - nfac);
  antscout << " Genus " << g << std::endl;

  antscout << " face " << nfac << " edg " << nedg <<  " vert " << vers << std::endl;

  return g;
}

float vtkComputeTopology(vtkPolyData* pd)
{
  // Marching cubes
//    antscout << " Marching Cubes ";
//    vtkMarchingCubes *marchingCubes = vtkMarchingCubes::New();
//    vtkContourFilter *marchingCubes = vtkContourFilter::New();
//    vtkKitwareContourFilter *marchingCubes = vtkKitwareContourFilter::New();
//    marchingCubes->SetInput((vtkDataSet*) vds);
//    marchingCubes->SetValue(0, hithresh);
//    int nc;
//    antscout << " Input #conts "; std::cin >> nc;
//    marchingCubes->SetNumberOfContours(2);
//    marchingCubes->SetComputeScalars(false);
//    marchingCubes->SetComputeGradients(false);
//    marchingCubes->SetComputeNormals(false);
  vtkPolyDataConnectivityFilter* con = vtkPolyDataConnectivityFilter::New();

  con->SetExtractionModeToLargestRegion();
//    con->SetInput(marchingCubes->GetOutput());
  con->SetInput(pd);
  con->Update();
  float g = ComputeGenus(con->GetOutput() );
  return g;
//    vtkUnstructuredGridToPolyDataFilter* gp = vtkUnstructuredGridToPolyDataFilter::New();
//    gp->SetInput(con->GetOutput());

//    marchingCubes->Delete();
#if 0
  int inputNumberOfPoints = con->GetOutput()->GetNumberOfPoints();
  int inputNumberOfPolys = con->GetOutput()->GetNumberOfPolys();

  vtkPolyDataConnectivityFilter *polyDataConnectivityFilter =
    vtkPolyDataConnectivityFilter::New();
  polyDataConnectivityFilter->SetInput( con->GetOutput() );
  polyDataConnectivityFilter->SetExtractionModeToAllRegions();
  polyDataConnectivityFilter->SetExtractionModeToLargestRegion();
  polyDataConnectivityFilter->Update();

  int connectivityNumberOfExtractedRegions = polyDataConnectivityFilter->
    GetNumberOfExtractedRegions();

  polyDataConnectivityFilter->Delete();

  vtkExtractEdges *extractEdges = vtkExtractEdges::New();
  extractEdges->SetInput( con->GetOutput() );
  extractEdges->Update();

  int extractNumberOfLines = extractEdges->GetOutput()->GetNumberOfLines();

  extractEdges->Delete();

  int EulerCharacteristic = inputNumberOfPoints - extractNumberOfLines
    + inputNumberOfPolys;

  double genus = 0.5 * ( 2 * connectivityNumberOfExtractedRegions
                         - EulerCharacteristic );

  antscout << "EulerCharacteristic " << EulerCharacteristic << std::endl;
  antscout << "genus " << genus << std::endl;

  return genus;
#endif
}

template <class TImage>
void GetValueMesh(typename TImage::Pointer image, typename TImage::Pointer image2,  std::string outfn,
                  const char* paramname, float scaledata,
                  float aaParm )
{
  //  antscout << " parname " << std::string(paramname) << std::endl;
  typedef TImage      ImageType;
  typedef ImageType   itype;
  typedef vtkPolyData MeshType;

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(0.8);
  filter->SetMaximumError(.01f);
  filter->SetUseImageSpacingOn();
  filter->SetInput(image);
  filter->Update();

  typedef BinaryImageToMeshFilter<ImageType> FilterType;
  typename  FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput( image );
  fltMesh->SetAntiAliasMaxRMSError( aaParm ); // to do nothing, set negative
  fltMesh->SetSmoothingIterations( 0 );
  fltMesh->Update();
  vtkPolyData* vtkmesh = fltMesh->GetMesh();
// assign scalars to the original surface mesh
//  Display((vtkUnstructuredGrid*)vtkmesh);

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetInput(vtkmesh);
  smoother->SetNumberOfIterations(25);
  smoother->BoundarySmoothingOff();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(120.0);
  smoother->SetPassBand(.01);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->Update();
  vtkmesh = smoother->GetOutput();

  antscout << " Genus " << vtkComputeTopology(vtkmesh) << std::endl;

  typename itype::SpacingType spacing = image->GetSpacing();

  vtkPoints* vtkpoints = vtkmesh->GetPoints();
  int        numPoints = vtkpoints->GetNumberOfPoints();
  float      mx = 0, mn = 9.e9, meank = 0;
  for( int i = 0; i < numPoints; i++ )
    {
    typename ImageType::IndexType index;
    typename ImageType::PointType point;
    for( int j = 0; j < 3; j++ )
      {
      point[j] = (vtkpoints->GetPoint(i)[j]);
      }
    image2->TransformPhysicalPointToIndex(point, index);
    float temp = image2->GetPixel(index);
    if( fabs(temp) > mx )
      {
      mx = fabs(temp);
      }
    if( fabs(temp) < mn && mn > 0 )
      {
      mn = fabs(temp);
      }
    meank += fabs(temp);
    }
  antscout << " max kap " << mx << " mn k " << mn <<  std::endl;
  meank /= numPoints;
//  mx=1.3;
//  mx=2.0;

  vtkFloatArray* param;

  // while (!done)
    {
    param = vtkFloatArray::New();
    param->SetName(paramname);
    float dif = (mx - mn) * scaledata;
    float mx2 = meank + dif;
    float mn2 = meank - dif;
    dif = mx2 - mn2;
    for( int i = 0; i < numPoints; i++ )
      {
      typename ImageType::IndexType index;
      typename ImageType::PointType point;
      for( int j = 0; j < 3; j++ )
        {
        point[j] = (vtkpoints->GetPoint(i)[j]);
        }
      image2->TransformPhysicalPointToIndex(point, index);
      float temp = image2->GetPixel(index);
      //    param->InsertNextValue(temp);
      //    float temp=surfk->CurvatureAtIndex(index);
      if( i % 1000 == 0 )
        {
        antscout << " kappa " << temp << std::endl;
        }
      // =fabs(manifoldIntegrator->GetGraphNode(i)->GetTotalCost());

      temp = fabs(temp);
      float vvv = (temp - mn2) * 255. / dif;
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
    //  Display((vtkUnstructuredGrid*)vtkmesh);
//  antscout<<"DOne? "; std::cin >> done;
    }
  antscout << " done with mesh map ";
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(vtkmesh);
  antscout << " writing " << outfn << std::endl;
  // outnm="C:\\temp\\mesh.vtk";
  writer->SetFileName(outfn.c_str() );
  writer->SetFileTypeToBinary();
  writer->Update();
  antscout << " done writing ";
  return;
}

template <class TImage>
float GetImageTopology(typename TImage::Pointer image)
{
  typedef TImage      ImageType;
  typedef vtkPolyData MeshType;

  double aaParm = 0.024;
  typedef BinaryImageToMeshFilter<ImageType> FilterType;
  typename  FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput(image);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->SetAntiAliasMaxRMSError( -1000.0 );   // to do nothing
  fltMesh->Update();
  vtkPolyData* vtkmesh = fltMesh->GetMesh();
// assign scalars to the original surface mesh
//  Display((vtkUnstructuredGrid*)vtkmesh);

  float genus =  vtkComputeTopology(vtkmesh);
  antscout << " Genus " << genus << std::endl;

  return genus;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int GetMeshAndTopology( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "GetMeshAndTopology" );

  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  if( argc < 2 )
    {
    antscout << argv[0] << " binaryimage valueimage  out paramname ValueScale AntiaAliasParm=0.001" << std::endl;
    antscout << " outputs vtk version of input image -- assumes object is defined by non-zero values " << std::endl;
    antscout << " mesh is colored by the value of the image voxel " << std::endl;
    antscout <<  " the AA-Param could cause topo problems but makes nicer meshes  " << std::endl;
    antscout << " ValueScale controls contrast in image appearance - lower increaseses -- should be <= 1 "
             << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int Dimension = 3;
  typedef float PixelType;
  // Declare the types of the output images
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::Image<PixelType, 2>         Image2DType;

  // Declare the type of the index,size and region to initialize images
  typedef itk::Index<Dimension>               IndexType;
  typedef itk::Size<Dimension>                SizeType;
  typedef itk::ImageRegion<Dimension>         RegionType;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  // Declare the type of the Mesh

  std::string outfn = std::string(argv[3]);

  ImageType::Pointer image2;
  ImageType::Pointer image;
  ReadImage<ImageType>(image, argv[1]);
  ReadImage<ImageType>(image2, argv[2]);

  ImageType::DirectionType fmat = image->GetDirection();
  fmat.SetIdentity();
  image->SetDirection(fmat);
  image2->SetDirection(fmat);

  // Save the mesh
  float       aaParm = 0.03;
  const char* paramname = std::string("f(x)").c_str();
  if( argc > 4 )
    {
    paramname = (argv[4]);
    }
  float scaledata = 0.5;
  if( argc > 5 )
    {
    scaledata = atof(argv[5]);
    }
  if( argc > 6 )
    {
    aaParm = atof(argv[6]);
    }
  antscout << "aaParm " << aaParm << std::endl;
  GetValueMesh<ImageType>(image, image2, outfn, paramname, scaledata, aaParm);
  //  GetImageTopology<ImageType>(image);

  return 0;
}
} // namespace ants
