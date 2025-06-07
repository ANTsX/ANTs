/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include <algorithm>
#include <string>

#include <cmath>
#include <ctime>

#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMesh.h"
#include "itkSphereMeshSource.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkFEMConformalMap.h"
#include "itkFEMDiscConformalMap.h"

#include "BinaryImageToMeshFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "ReadWriteData.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkSurfaceImageCurvature.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"
#include <vtkSmartPointer.h>
#include <vtkWindowedSincPolyDataFilter.h>

namespace ants
{
/*

OPEN QUESTIONS:

  LOCAL FROM GLOBAL FOR TRIANGLE IN 3D (MARCELO)
    USEFUL FOR INTERPOLATING STEREOGRAPHIC COORDS INTO IMAGE.
    CURRENTLY JUST DENSELY INTERPOLATE ACROSS EACH TRIANGLE. -- Not necessary for triangle

  TRIANGLE AREA IN 3D (D0NE : USING CROSS PRODUCT )

  - need to insure we use the hypotneuse as the AB edge : so theta is <= 1.0;

  - try to apply bcs to even out the map

  - what is the effect of material param E?

  */

void
Display(vtkUnstructuredGrid * vtkgrid, bool secondwin = false, bool delinter = true)
{
  // Create the renderer and window stuff
  std::cout << " second win " << secondwin << std::endl;
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

  //   vtkCellPicker* cpicker = vtkCellPicker::New();

  vtkCallbackCommand * cbc = vtkCallbackCommand::New();
  /*-------------------------------------------------------------------------
  void vtkCompositeManagerExitInteractor(vtkObject *vtkNotUsed(o),
                                         unsigned long vtkNotUsed(event),
                                         void *clientData, void *)
  {
    vtkCompositeManager *self = (vtkCompositeManager *)clientData;

    self->ExitInteractor();
  }*/

  //  cbc=vectraEventHandler; //-> SetCallback(vectraEventHandler);
  ren1->AddObserver(vtkCommand::KeyPressEvent, cbc);

  vtkDataSetMapper * mapper = vtkDataSetMapper::New();
  mapper->SetInput(vtkgrid);
  mapper->SetScalarRange(0, 255);
  vtkActor * actor = vtkActor::New();
  actor->SetMapper(mapper);
  vtkDataSetMapper * mapper2 = vtkDataSetMapper::New();
  if (secondwin)
  {
    mapper2->SetInput(vtkgrid);
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
  ren1->AddActor(actor);

  //  this->InteractionPicker->Pick(x, y, 0.0, this->CurrentRenderer);
  //  this->InteractionPicker->GetPickPosition(this->DownPt);

  renWin->Render();
  inter->Start();
  /*
        int X, Y;
        char keypressed = *(inter -> GetKeySym());
        inter -> GetMousePosition(&X, &Y);
        std::cout <<" X Y " << X << " " << Y << std::endl;
    renWin->Render();
    inter->Start();*/

  /*
    vtkPointPicker* picker=vtkPointPicker::New();
  //  vtkPoints *GetPickedPositions() {return this->PickedPositions;};
    float x = inter->GetEventPosition()[0];
    float y = inter->GetEventPosition()[1];
    float z = 0.0;
    picker->Pick(x,y,z,ren1);
    std::cout <<" picked " << (*picker->GetPickedPositions()) << std::endl;*/

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

template <typename TImage>
void
MapToSphere(typename TImage::Pointer image, int fixdir, float e)
{
  typedef TImage                                           ImageType;
  typedef unsigned char                                    PixelType;
  typedef itk::Mesh<float>                                 MeshType;
  typedef itk::BinaryMask3DMeshSource<ImageType, MeshType> MeshSourceType;

  PixelType                        internalValue = 1;
  typename MeshSourceType::Pointer meshSource = MeshSourceType::New();
  meshSource->SetBinaryImage(image);
  meshSource->SetObjectValue(internalValue);
  try
  {
    meshSource->Update();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cout << "Exception thrown during Update() " << std::endl;
    std::cout << exp << std::endl;
    return;
  }
  meshSource->GetOutput();
  std::cout << meshSource->GetNumberOfNodes() << std::endl;
  std::cout << meshSource->GetNumberOfCells() << std::endl;

  typedef itk::FEMConformalMap<MeshType, ImageType> ParamType;
  typename ParamType::Pointer                       Parameterizer = ParamType::New();

  Parameterizer->SetDebug(false);
  //  Parameterizer->SetDebug(true);
  Parameterizer->SetReadFromFile(false);
  Parameterizer->SetParameterFileName("");
  Parameterizer->SetImage(image);
  std::cout << " fixdir " << fixdir << " e " << e << " best e ~ 3.e-3 " << std::endl;

  Parameterizer->SetNorthPole(fixdir);
  Parameterizer->SetSigma(e);
  Parameterizer->SetSurfaceMesh(meshSource->GetOutput());
  //  Parameterizer->GenerateSystemFromSurfaceMesh();
  std::cout << std::endl;

  Parameterizer->ConformalMap();
  Parameterizer->ComputeStereographicCoordinates();

  if (Parameterizer->GetImage())
  {
    std::cout << " writing param images " << std::endl;
    {
      Parameterizer->MapCheckerboardToImage(0.05);
      typename itk::ImageFileWriter<ImageType>::Pointer writer;
      writer = itk::ImageFileWriter<ImageType>::New();
      std::string fn = "checker.nii";
      writer->SetFileName(fn.c_str());
      writer->SetInput(Parameterizer->GetImage());
      writer->Write();
    }
  }
}

template <typename TImage>
typename TImage::Pointer
SmoothImage(typename TImage::Pointer input, double var)
{
  typedef TImage                        ImageType;
  typedef typename ImageType::PixelType PixelType;
  enum
  {
    ImageDimension = ImageType::ImageDimension
  };
  typedef itk::Image<float, ImageDimension> realImageType;

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;

  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance(var);
  filter->SetMaximumError(.01f);
  filter->SetUseImageSpacing(true);
  filter->SetInput(input);
  filter->Update();

  typedef itk::CastImageFilter<ImageType, ImageType> CasterType1;
  typename CasterType1::Pointer                      caster1 = CasterType1::New();
  caster1->SetInput(filter->GetOutput());
  caster1->Update();

  return caster1->GetOutput();
}

float
ComputeGenus(vtkPolyData * pd1)
{
  vtkExtractEdges * edgeex = vtkExtractEdges::New();

  edgeex->SetInput(pd1);
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
  con->SetInput(pd);

  //    vtkUnstructuredGridToPolyDataFilter* gp = vtkUnstructuredGridToPolyDataFilter::New();
  //    gp->SetInput(con->GetOutput());

  float g = ComputeGenus(con->GetOutput());
  //    marchingCubes->Delete();

  return g;
}

template <typename TImage>
void
GetMeshAndCurvature(typename TImage::Pointer image, float e, const char * filename)
{
  typedef TImage      ImageType;
  typedef vtkPolyData MeshType;

  double                                     aaParm = 0.024;
  typedef BinaryImageToMeshFilter<ImageType> FilterType;
  typename FilterType::Pointer               fltMesh = FilterType::New();
  fltMesh->SetInput(image);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->SetAntiAliasMaxRMSError(-1000.0); // to do nothing
  fltMesh->Update();
  vtkPolyData * vtkmesh = fltMesh->GetMesh();
  // assign scalars to the original surface mesh
  //  Display((vtkUnstructuredGrid*)vtkmesh);

  std::cout << " Genus " << vtkComputeTopology(vtkmesh) << std::endl;

  typedef itk::SurfaceImageCurvature<ImageType> surfktype;
  typename surfktype::Pointer                   surfk = surfktype::New();
  //  kappa=;

  typedef itk::Image<float, 3> FloatImageType;
  surfk->SetInput(image); // SmoothImage<ImageType>(image,1.0));
  surfk->SetNeighborhoodRadius(1.5);
  surfk->SetSigma(1.);
  surfk->SetUseLabel(false);
  surfk->SetUseGeodesicNeighborhood(false);
  surfk->SetkSign(1.0);
  surfk->ComputeFrameOverDomain(3);
  typename FloatImageType::Pointer kappa = surfk->GetFunctionImage();

  typedef itk::Image<unsigned char, 3>                            itype;
  typedef itk::RescaleIntensityImageFilter<FloatImageType, itype> CastFilterType;
  typename CastFilterType::Pointer                                caster = CastFilterType::New();
  caster->SetInput(image);
  caster->SetOutputMinimum(0);
  caster->SetOutputMaximum(255);
  std::string fn = std::string(filename);
  fn = fn.substr(0, fn.length() - 4) + "kappa.nii";
  typedef itk::ImageFileWriter<itype> writertype;
  typename writertype::Pointer        w = writertype::New();
  typename itype::Pointer             kimg = caster->GetOutput();
  w->SetInput(kimg);
  w->SetFileName(fn.c_str());
  w->Write();

  typename itype::SpacingType spacing = image->GetSpacing();

  vtkPoints * vtkpoints = vtkmesh->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  float       mx = 0, mn = 9.e9, meank = 0;
  for (int i = 0; i < numPoints; i++)
  {
    typename ImageType::IndexType index;
    for (int j = 0; j < 3; j++)
    {
      index[j] = (int)(vtkpoints->GetPoint(i)[j] / spacing[j] + 0.5);
    }
    float temp = kimg->GetPixel(index);
    //    float temp=image->GetPixel(index);
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

  //   bool done=false;
  // while (!done)
  {
    vtkFloatArray * param = vtkFloatArray::New();
    param->SetName("angle");
    float       dif = (mx - mn) * 0.25;
    const float mx2 = meank + dif;
    const float mn2 = meank - dif;
    dif = mx2 - mn2;
    for (int i = 0; i < numPoints; i++)
    {
      typename ImageType::IndexType index;
      for (int j = 0; j < 3; j++)
      {
        index[j] = (int)(vtkpoints->GetPoint(i)[j] / spacing[j] + 0.5);
      }
      float temp = kimg->GetPixel(index);
      //    float temp=surfk->CurvatureAtIndex(index);
      if (i % 1000 == 0)
      {
        std::cout << " kappa " << temp << std::endl;
      }
      // =fabs(manifoldIntegrator->GetGraphNode(i)->GetTotalCost());
      temp = fabs(temp);
      param->InsertNextValue((temp - mn2) * 255. / dif);
    }
    vtkmesh->GetPointData()->SetScalars(param);
    //  Display((vtkUnstructuredGrid*)vtkmesh);
    //  std::cout<<"DOne? "; std::cin >> done;
  }
  std::cout << " done with curvature map ";
  vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
  writer->SetInput(vtkmesh);
  std::string outnm = std::string(filename);
  outnm = outnm.substr(0, outnm.length() - 4) + ".vtk";
  std::cout << " writing " << outnm << std::endl;
  // outnm="C:\\temp\\mesh.vtk";
  writer->SetFileName(outnm.c_str());
  writer->SetFileTypeToBinary();
  writer->Update();
  std::cout << " done writing ";
  return;
}

template <typename TImage>
float
GetImageTopology(typename TImage::Pointer image, float e, const char * filename)
{
  typedef TImage      ImageType;
  typedef vtkPolyData MeshType;

  double                                     aaParm = 0.024;
  typedef BinaryImageToMeshFilter<ImageType> FilterType;
  typename FilterType::Pointer               fltMesh = FilterType::New();
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

template <typename TImage>
void
MapToDisc(vtkPolyData * vtkmesh, float e, std::string outfn)
{
  typedef TImage                                        ImageType;
  typedef vtkPolyData                                   MeshType;
  typedef itk::FEMDiscConformalMap<MeshType, ImageType> ParamType;
  typename ParamType::Pointer                           Parameterizer = ParamType::New();
  Parameterizer->SetDebug(false);
  //  Parameterizer->SetDebug(true);
  Parameterizer->SetReadFromFile(false);
  Parameterizer->SetParameterFileName("");
  Parameterizer->SetSigma(e);
  Parameterizer->SetSurfaceMesh(vtkmesh);

  ComputeGenus((vtkPolyData *)Parameterizer->GetSurfaceMesh());
  Parameterizer->ExtractSurfaceDisc();
  ComputeGenus((vtkPolyData *)Parameterizer->GetSurfaceMesh());
  //  Display((vtkUnstructuredGrid*)Parameterizer->GetSurfaceMesh());

  std::cout << " begin conformal mapping ";
  Parameterizer->ConformalMap();

  std::cout << " display patch ";
  //  ComputeGenus((vtkPolyData*)Parameterizer->m_ExtractedSurfaceMesh);
  // Display((vtkUnstructuredGrid*)Parameterizer->m_ExtractedSurfaceMesh);

  //  std::cout << " display flattened patch ";

  //   float step  = 0.1;
  //   float maxt=0.0;
  //  for (float tt = 0.0; tt<=maxt; tt=tt+step)
  {
    //      std::cout <<" Building at : " << tt << std::endl;
    // /  Parameterizer->BuildOutputMeshes(tt);

    //      if (tt == 0.)
    {
      vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
      writer->SetInput(Parameterizer->m_ExtractedSurfaceMesh);
      std::string outnm; // =std::string(filename);
      // outnm=outnm.substr(0,outnm.length()-4)+".vtk";
      outnm = outfn + "mapspace.vtk";
      std::cout << " writing " << outnm << std::endl;
      writer->SetFileName(outnm.c_str());
      writer->SetFileTypeToBinary();
      writer->Update();
    }

    {
      vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
      writer->SetInput(Parameterizer->m_DiskSurfaceMesh);
      std::string outnm;
      outnm = outfn + "mapflat.vtk";
      std::cout << " writing " << outnm << std::endl;
      writer->SetFileName(outnm.c_str());
      writer->SetFileTypeToBinary();
      writer->Update();
    }
  }
  /*
    for (float sig=0.4; sig<=1.0; sig=sig+0.1)
    {
    Parameterizer->SetSmooth(sig);
    Parameterizer->ConjugateHarmonic();
  //  ComputeGenus((vtkPolyData*)Parameterizer->m_DiskSurfaceMesh);
    Display((vtkUnstructuredGrid*)Parameterizer->m_DiskSurfaceMesh);
    }
  */

  typedef typename ParamType::FlatImageType      imtype;
  typename itk::ImageFileWriter<imtype>::Pointer writer;
  writer = itk::ImageFileWriter<imtype>::New();
  std::string fn = outfn + "flat.nii";
  writer->SetFileName(fn.c_str());
  writer->SetInput(Parameterizer->m_FlatImage);
  if (Parameterizer->m_FlatImage)
  {
    writer->Write();
  }

  std::cout << " done writing ";
}

/*
template <typename TImage>
void MeshToImage(vtkPolyData* vtkmesh, int imagesize, std::string outfn)
{
    constexpr unsigned int Dimension = 2;
    constexpr int splineOrder = 3;
    constexpr int component = 0;
    const char* dataname = "scalars";
    typedef TImage ImageType;
    typedef vtkPolyData MeshType;
    typedef typename ImageType::PixelType PixelType;
    typedef typename itk::Vector<PixelType, 1> FunctionalDataType;
    typedef typename itk::Image<FunctionalDataType, Dimension> FunctionalImageType;

    int numVertices = vtkmesh->GetNumberOfPoints();
    vtkmesh->GetPointData()->SetActiveScalars( dataname );
    std::cout << "#of vertices " << numVertices << " Fitting variable " <<
vtkmesh->GetPointData()->GetScalars()->GetName( ) << std::endl;

    typedef itk::PointSet<FunctionalDataType, Dimension> FunctionalMapType;
    typedef itk::BSplineScatteredDataPointSetToImageFilter
        <FunctionalMapType, FunctionalImageType> BSplineFilterType;

    typename FunctionalMapType::Pointer funcdataPoints =
        FunctionalMapType::New();
    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

    FunctionalDataType data, data1;
    typename FunctionalMapType::PointType point;
    double bounds[6];
    vtkPoints *meshpoints = vtkmesh->GetPoints();
    meshpoints->ComputeBounds();
    meshpoints->GetBounds( bounds );

    for (int ID=0; ID < numVertices; ++ID)
    {
        data[0] = vtkmesh->GetPointData()->GetScalars()->GetComponent(ID, component);
        for (int dir=0; dir < Dimension; ++dir)
            point[dir] = meshpoints->GetPoint(ID)[dir];
        //Debug
        //data[0] = 1*point[0] + 0*point[1];


        funcdataPoints->SetPointData( ID, data);
        funcdataPoints->SetPoint( ID, point);
        //std::cout << ID << " " << point[0] << " " << point[1] << " " << data[0] << std::endl;
    }

    std::cout << "Read all mesh data" << std::endl;

    typename FunctionalImageType::PointType origin;
    typename FunctionalImageType::SpacingType spacing;
    typename FunctionalImageType::SizeType size;
    float maxX, minX, maxY, minY;
       maxX = static_cast<float>(bounds[1]);
       minX = static_cast<float>(bounds[0]);
       maxY = static_cast<float>(bounds[3]);
       minY = static_cast<float>(bounds[2]);
    std::cout << "minX " << minX << " maxX " << maxX << " minY " << minY << " maxY " << maxY << std::endl;

    size[0] = imagesize;
    size[1] = imagesize;

    origin[0] = (minX + maxX)/2;
    origin[1] = (minY + maxY)/2;
    origin[0] = minX - (maxX-minX)*0.1;
    origin[1] = minY - (maxY-minY)*0.1;
    origin[0] = minX;
    origin[1] = minY;

        spacing[0] = 1.1*((((maxX-minX)>(maxY-minY))?(maxX-minX):(maxY-minY)))/imagesize;
    spacing[1] = spacing[0];

    std::cout << "size " << size << " origin " << origin << " spacing " << spacing << std::endl;
        typename BSplineFilterType::ArrayType ncps;
    ncps.Fill(  splineOrder + 1 );

    //const int numLevels = round((1/static_cast<double>(Dimension))*log(static_cast<double>(numVertices))/log(2.0));
    constexpr int numLevels = 13;

      typename BSplineFilterType::ArrayType close;
     close.Fill( false );

    std::cout << "No. of levels " << numLevels << " splineOrder " << splineOrder << " #control points " << ncps <<
std::endl;

        bspliner->SetOrigin( origin );
        bspliner->SetSpacing( spacing );
        bspliner->SetSize( size );
        bspliner->SetGenerateOutputImage( true );
        bspliner->SetNumberOfLevels( numLevels );
        bspliner->SetSplineOrder( splineOrder );
        bspliner->SetNumberOfControlPoints( ncps );
        bspliner->SetCloseDimension( close );
    bspliner->SetInput( funcdataPoints );
    bspliner->DebugOn();


    std::cout << "Entering BSpline" << std::endl;
    bspliner->Update();
    std::cout << "BSpline fitting done" << std::endl;

    // ORIENTATION ALERT  -- the original code here
    // set the region, spacing, and origin without setting directions.
    typename ImageType::Pointer outimage =
      AllocImage<ImageType>(bspliner->GetOutput());

    typename itk::ImageRegionIterator<ImageType>
        ItO( outimage, outimage->GetRequestedRegion() );
    typename itk::ImageRegionConstIterator<FunctionalImageType>
        ItB( bspliner->GetOutput(), bspliner->GetOutput()->GetRequestedRegion() );

    for ( ItB.GoToBegin(), ItO.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItO )
        {
            ItO.Set( static_cast<PixelType>( ItB.Get()[0] ) );
        }

      typename itk::ImageFileWriter<ImageType>::Pointer writer;
    writer = itk::ImageFileWriter<ImageType>::New();
      std::string fn=outfn+"flat.nii";
      writer->SetFileName(fn.c_str());
      writer->SetInput( outimage );
      writer->Write();
*/
/*** Evaluate quality of fit ***/
/*
typedef typename itk::LinearInterpolateImageFunction< ImageType, float > InterpolatorType;
typename InterpolatorType::PointType testpoint;
typename InterpolatorType::Pointer interp = InterpolatorType::New();
interp->SetInputImage( outimage );
for (int ID=0; ID < numVertices; ++ID)
{
data[0] = vtkmesh->GetPointData()->GetScalars()->GetComponent(ID, component);
for (int dir=0; dir < Dimension; ++dir)
{
    point[dir] = meshpoints->GetPoint(ID)[dir];
    testpoint[dir] = point[dir] - origin[dir];
}
vtkmesh->GetPointData()->GetScalars()->SetComponent(ID, component, static_cast<int>(data[0]) -
static_cast<int>(interp->Evaluate( testpoint )));

//vtkmesh->GetPointData()->GetScalars()->SetComponent(ID, component,  interp->Evaluate( testpoint ));
data1[0] = vtkmesh->GetPointData()->GetScalars()->GetComponent(ID, component);
cout << "error  " << data1[0] << " original " << data[0] << " at " << point << " interpolated " << interp->Evaluate(
testpoint ) << " at " << testpoint << endl;
}
vtkPolyDataWriter *vtkwriter = vtkPolyDataWriter::New();
vtkwriter->SetInput( vtkmesh );
vtkwriter->SetFileName( "newmesh3.vtk" );
vtkwriter->Write();
*/

// }

template <typename TImage>
typename TImage::Pointer
RemoveNaNs(typename TImage::Pointer image, float replaceval)
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator                                          vfIter(image, image->GetLargestPossibleRegion());
  for (vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter)
  {
    typename TImage::PixelType v1 = vfIter.Get();
    if (std::isnan(v1))
    {
      vfIter.Set(replaceval);
    }
  }
  return image;
}

/*
template <typename TImage>
void ImageToMesh(vtkPolyData* vtkmesh, typename TImage::Pointer image, std::string outfn)
{
    constexpr unsigned int Dimension = 2;
    constexpr int splineOrder = 3;
    constexpr int imagesize = 512;
    constexpr int component = 0;
    const char* dataname = "points";
    typedef TImage ImageType;
    typedef vtkPolyData MeshType;
    typedef typename ImageType::PixelType PixelType;
    typedef typename itk::Vector<PixelType, 1> FunctionalDataType;
    typedef typename itk::Image<FunctionalDataType, Dimension> FunctionalImageType;

    float replaceval = 5;
    image = RemoveNaNs<TImage>(image, replaceval);

    int numVertices = vtkmesh->GetNumberOfPoints();
    vtkmesh->GetPointData()->SetActiveScalars( dataname );
    std::cout << "#of vertices " << numVertices << " Fitting variable " <<
vtkmesh->GetPointData()->GetScalars()->GetName( ) << std::endl;

    typedef itk::PointSet<FunctionalDataType, Dimension> FunctionalMapType;
    typedef itk::BSplineScatteredDataPointSetToImageFilter
        <FunctionalMapType, FunctionalImageType> BSplineFilterType;

    typename FunctionalMapType::Pointer funcdataPoints =
        FunctionalMapType::New();
    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

    FunctionalDataType data;
    typename FunctionalMapType::PointType point;
    double bounds[6];
    vtkPoints *meshpoints = vtkmesh->GetPoints();
    meshpoints->ComputeBounds();
    meshpoints->GetBounds( bounds );

    for (int ID=0; ID < numVertices; ++ID)
    {
        data[0] = vtkmesh->GetPointData()->GetScalars()->GetComponent(ID, component);
        for (int dir=0; dir < Dimension; ++dir)
            point[dir] = meshpoints->GetPoint(ID)[dir];
        //Debug
        //data[0] = 1*point[0] + 0*point[1];


        funcdataPoints->SetPointData( ID, data);
        funcdataPoints->SetPoint( ID, point);
        //std::cout << ID << " " << point[0] << " " << point[1] << " " << data[0] << std::endl;
    }

    std::cout << "Read all mesh data" << std::endl;

    // Create new data array
    vtkFloatArray *newdata;
    newdata = vtkFloatArray::New();
    newdata->SetName("atlas");

    typename FunctionalImageType::PointType origin;
    typename FunctionalImageType::SpacingType spacing;
    typename FunctionalImageType::SizeType size;
    float maxX, minX, maxY, minY;
       maxX = static_cast<float>(bounds[1]);
       minX = static_cast<float>(bounds[0]);
       maxY = static_cast<float>(bounds[3]);
       minY = static_cast<float>(bounds[2]);
    std::cout << "minX " << minX << " maxX " << maxX << " minY " << minY << " maxY " << maxY << std::endl;

    size[0] = imagesize;
    size[1] = imagesize;

    origin[0] = (minX + maxX)/2;
    origin[1] = (minY + maxY)/2;
    origin[0] = minX - (maxX-minX)*0.1;
    origin[1] = minY - (maxY-minY)*0.1;
    origin[0] = minX;
    origin[1] = minY;

        spacing[0] = 1.1*((((maxX-minX)>(maxY-minY))?(maxX-minX):(maxY-minY)))/imagesize;
    spacing[1] = spacing[0];

    std::cout << "size " << size << " origin " << origin << " spacing " << spacing << std::endl;
*/
/* bspline code -- not doing b-spline now

    typename BSplineFilterType::ArrayType ncps;
ncps.Fill(  splineOrder + 1 );

//const int numLevels = round((1/static_cast<double>(Dimension))*log(static_cast<double>(numVertices))/log(2.0));
constexpr int numLevels = 13;

  typename BSplineFilterType::ArrayType close;
 close.Fill( false );

std::cout << "No. of levels " << numLevels << " splineOrder " << splineOrder << " #control points " << ncps <<
std::endl;

    bspliner->SetOrigin( origin );
    bspliner->SetSpacing( spacing );
    bspliner->SetSize( size );
    bspliner->SetGenerateOutputImage( true );
    bspliner->SetNumberOfLevels( numLevels );
    bspliner->SetSplineOrder( splineOrder );
    bspliner->SetNumberOfControlPoints( ncps );
    bspliner->SetCloseDimension( close );
bspliner->SetInput( funcdataPoints );
bspliner->DebugOn();


std::cout << "Entering BSpline" << std::endl;
bspliner->Update();
std::cout << "BSpline fitting done" << std::endl;

typename ImageType::Pointer outimage = ImageType::New();
outimage->SetSpacing( bspliner->GetOutput()->GetSpacing() );
outimage->SetOrigin( bspliner->GetOutput()->GetOrigin() );
outimage->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion() );
outimage->AllocateInitialized();

*/
/*
typename itk::ImageRegionIterator<ImageType>
    ItO( image, image->GetRequestedRegion() );



// Evaluate quality of fit
//typedef typename itk::LinearInterpolateImageFunction< ImageType, float > InterpolatorType;
typedef typename itk::NearestNeighborInterpolateImageFunction< ImageType, float > InterpolatorType;
typename InterpolatorType::PointType testpoint;
typedef typename itk::ContinuousIndex< float, Dimension > ContinuousIndexType;

typename InterpolatorType::Pointer interp = InterpolatorType::New();
interp->SetInputImage( image );
    ContinuousIndexType contind;


for (int ID=0; ID < numVertices; ++ID)
{
    for (int dir=0; dir < Dimension; ++dir)
        testpoint[dir] = meshpoints->GetPoint(ID)[dir] - origin[dir];
    newdata->InsertNextValue( interp->Evaluate( testpoint ) );
    //vtkmesh->GetPointData()->GetScalars()->SetComponent(ID, component, (interp->Evaluate( testpoint )));

    //vtkmesh->GetPointData()->GetScalars()->SetComponent(ID, component,  interp->Evaluate( testpoint ));
//    cout << vtkmesh->GetPointData()->GetScalars()->GetComponent(ID, component) << " " << data[0] << " " <<
interp->Evaluate( testpoint ) << endl;
}
vtkmesh->GetPointData()->AddArray( newdata );
vtkPolyDataWriter *vtkwriter = vtkPolyDataWriter::New();
vtkwriter->SetInput( vtkmesh );
  std::string fn=outfn+".vtk";
  vtkwriter->SetFileName(fn.c_str());
vtkwriter->Write();


}
*/

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ConformalMapping(std::vector<std::string> args, std::ostream * out_stream = nullptr)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ConformalMapping");
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

  // Define the dimension of the images
  constexpr unsigned int Dimension = 3;

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

  char *      filename;
  PixelType   loth = 1;
  PixelType   hith = 256;
  int         fixdir = 1;
  float       param = 2.e-4;
  std::string outfn;

  std::cout << "to get mesh: ConformalMapping   image.nii 1 255 1.e-3 0 outname " << std::endl;
  std::cout << "to flatten mesh: ConformalMapping   mesh.vtk 1 2 3 4 outname " << std::endl;
  std::cout << "to view mesh: ConformalMapping   mesh.vtk 1 2 3 1 outname " << std::endl;
  std::cout << "to get image topology: ConformalMapping   image.nii 1 255 1.e-3 5 outname " << std::endl;
  std::cout << "to convert flattened mesh to image : ConformalMapping   mesh.vtk 1 2 3 6 outname " << std::endl;
  std::cout << "to interpolate data in flattened image domain to original mesh: ConformalMapping   image.nii 1 2 3 7 "
               "outname originalflatmesh.vtk"
            << std::endl;
  std::cout << " to smooth a mesh --- ConformalMapping mesh.vtk 1 2 NumSmoothIts 8 outname " << std::endl;
  if (argc >= 2)
  {
    filename = argv[1];
    loth = atof(argv[2]);
    hith = atof(argv[3]);
    param = atof(argv[4]);
    fixdir = std::stoi(argv[5]);
    outfn = std::string(argv[6]);
  }
  else
  {
    filename = (char *)"image.nii";
    return 0;
  }

  if (fixdir == 0)
  {
    typedef itk::ImageFileReader<ImageType> FileSourceType;
    typedef ImageType::PixelType            PixType;
    FileSourceType::Pointer                 readfilter = FileSourceType::New();
    readfilter->SetFileName(filename);
    try
    {
      readfilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cout << "Exception caught during reference file reading " << std::endl;
      std::cout << e << std::endl;
      return -1;
    }

    ImageType::Pointer image = readfilter->GetOutput();

    ImageType::Pointer thresh = BinaryThreshold<ImageType>(loth, hith, hith, image);

    // Save the mesh
    GetMeshAndCurvature<ImageType>(thresh, param, filename);
    //  MapToSphere<ImageType>(thresh,fixdir,e);
  }
  else if (fixdir == 5)
  {
    typedef itk::ImageFileReader<ImageType> FileSourceType;
    typedef ImageType::PixelType            PixType;
    FileSourceType::Pointer                 readfilter = FileSourceType::New();
    readfilter->SetFileName(filename);
    try
    {
      readfilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cout << "Exception caught during reference file reading " << std::endl;
      std::cout << e << std::endl;
      return -1;
    }

    ImageType::Pointer image = readfilter->GetOutput();

    ImageType::Pointer thresh = BinaryThreshold<ImageType>(loth, hith, hith, image);

    // Save the mesh
    GetImageTopology<ImageType>(thresh, param, filename);
    //  MapToSphere<ImageType>(thresh,fixdir,e);
  }
  /*************** flat mesh to image ****************/
  else if (fixdir == 6)
  {
    vtkPolyDataReader * fltReader = vtkPolyDataReader::New();
    fltReader->SetFileName(filename);
    fltReader->Update();
    // MeshToImage<Image2DType>(fltReader->GetOutput(), 512, outfn);
  }
  /*************** flat mesh to image ****************/
  else if (fixdir == 7)
  /*************** flat map to mesh *****************/
  {
    char * const        refmeshname = argv[7];
    vtkPolyDataReader * fltReader = vtkPolyDataReader::New();
    fltReader->SetFileName(refmeshname);
    fltReader->Update();
    typedef itk::ImageFileReader<Image2DType> FileSourceType;
    FileSourceType::Pointer                   readfilter = FileSourceType::New();
    readfilter->SetFileName(filename);
    readfilter->Update();
    Image2DType::Pointer image = readfilter->GetOutput();
    // ImageToMesh<Image2DType>(fltReader->GetOutput(), image, outfn);
  }
  else if (fixdir == 8)
  /*************** flat map to mesh *****************/
  {
    std::cout << " read " << std::string(filename) << " write " << outfn << std::endl;
    vtkPolyDataReader * fltReader = vtkPolyDataReader::New();
    fltReader->SetFileName(filename);
    fltReader->Update();
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smoother->SetInput(fltReader->GetOutput());
    smoother->SetNumberOfIterations((int)param);
    smoother->BoundarySmoothingOn();
    smoother->FeatureEdgeSmoothingOff();
    smoother->SetFeatureAngle(180.0);
    smoother->SetEdgeAngle(180.0);
    smoother->SetPassBand(1.e-3); // smaller values increase smoothing
    smoother->NonManifoldSmoothingOn();
    smoother->NormalizeCoordinatesOff();
    smoother->Update();

    vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
    writer->SetInput(smoother->GetOutput());
    writer->SetFileName(outfn.c_str());
    writer->SetFileTypeToBinary();
    writer->Update();
    std::cout << " done writing ";
  }
  /*************** flat map to mesh ****************/
  else
  {
    vtkPolyDataReader * fltReader = vtkPolyDataReader::New();
    fltReader->SetFileName(filename);
    fltReader->Update();
    vtkPolyData * polydata = fltReader->GetOutput();

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInput(polydata);
    normalGenerator->Update();
    polydata = normalGenerator->GetOutput();
    if (fixdir == 1)
    {
      Display(static_cast<vtkUnstructuredGrid *>(polydata));
    }
    //    if ( fixdir == 1 )  Display((vtkUnstructuredGrid*)fltReader->GetOutput());
    std::cout << " m_Smooth " << param << std::endl;
    MapToDisc<ImageType>(fltReader->GetOutput(), param, outfn);
  }
  return 0;
}
} // namespace ants
