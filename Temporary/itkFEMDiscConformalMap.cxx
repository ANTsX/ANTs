/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMDiscConformalMap.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/07 15:30:39 $
  Version:   $Revision: 1.3 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detailm_Solver.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _FEMDiscConformalMap_txx
#define _FEMDiscConformalMap_txx

#include <vcl_cmath.h>
#include <vcl_iostream.h>
#include <vnl/vnl_real_polynomial.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_math.h>
#include "vtkDelaunay2D.h"
#include "vtkFloatArray.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkFEMDiscConformalMap.h"
#include "itkFEMLoadNode.h"
#include "itkSurfaceMeshCurvature.h"

namespace itk
{
template <typename TSurface, typename TImage, unsigned int TDimension>
FEMDiscConformalMap<TSurface, TImage, TDimension>
::FEMDiscConformalMap()
{
  m_Sigma = 2.0e-4;
  m_ParameterFileName = "";
  m_NorthPole = 0;
  m_SourceNodeNumber = 1;

  m_Pi = 3.14159265358979323846;
  m_ReadFromFile = true;
  m_Debug = false;
  m_FindingRealSolution = true;
  m_SurfaceMesh = NULL;
  for( int i = 0; i < 7; i++ )
    {
    m_PoleElementsGN[i] = 0;
    }

  m_Smooth = 2.0;

  this->m_FlatImage  = NULL;

  manifoldIntegrator = ManifoldIntegratorType::New();
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool  FEMDiscConformalMap<TSurface, TImage, TDimension>
::InBorder(typename FEMDiscConformalMap<TSurface, TImage, TDimension>::GraphSearchNodePointer g)
{
  float dist = g->GetTotalCost();

  if( dist > m_MaxCost && dist < m_MaxCost + 1.0 )
    {
    return true;
    }
  return false;

//	NodeLocationType loc;
//	loc=g->GetLocation()-manifoldIntegrator->GetGraphNode(0)->GetLocation();
//	float dist=loc.magnitude();
//	std::cout << " dist " << dist << std::endl;
//	if (dist > m_MaxCost && dist < m_MaxCost + 3) return true;

  //	if (g->GetTotalCost() > m_MaxCost) return true;
  return false;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool  FEMDiscConformalMap<TSurface, TImage, TDimension>
::InDisc(typename FEMDiscConformalMap<TSurface, TImage, TDimension>::GraphSearchNodePointer g)
{
  float d = g->GetTotalCost();

  if( d < m_MaxCost + 1. )
    {
    return true;
    }
  return false;

  if( g->GetPredecessor() )
    {
    return true;
    }
  else
    {
    return false;
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::FindSource(IndexType index)
{
  vtkPoints* vtkpoints = m_SurfaceMesh->GetPoints();
  int        numPoints = vtkpoints->GetNumberOfPoints();
  float      mindist = 1.e9;

  for( int i = 0; i < numPoints; i++ )
    {
    double* pt = vtkpoints->GetPoint(i);
    float   dist = 0.0;
    for( int j = 0; j < ImageDimension; j++ )
      {
      dist += (pt[j] - (float)index[j]) * (pt[j] - (float)index[j]);
      }
    dist = sqrt(dist);
    if( dist < mindist )
      {
      mindist = dist;
      m_SourceNodeNumber = i;
      }
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ExtractSurfaceDisc()
{
  std::cout << " set surface mesh " << std::endl;

  manifoldIntegrator->SetSurfaceMesh(m_SurfaceMesh);
  std::cout << " begin initializing graph " << std::endl;
  manifoldIntegrator->InitializeGraph3();
  m_SurfaceMesh = manifoldIntegrator->GetSurfaceMesh();
//   float frac=0;
  IndexType index;
  std::cout << " enter src node coords ";  std::cin >> index[0] >> index[1] >> index[2];
  std::cout << " entered index " << index;
  this->FindSource(index);
// 181,158,113  if (frac > 0.99) frac=0.99; else if (frac < 0) frac=0;
//  this->m_SourceNodeNumber=(int)(float)(manifoldIntegrator->GetGraphSize())*frac;
  manifoldIntegrator->SetSource(manifoldIntegrator->GetGraphNode(this->m_SourceNodeNumber) );
  manifoldIntegrator->InitializeQueue();
  std::cout << " done initializing graph " << std::endl;
  float mc = 0;
  std::cout << " Enter max cost ";  std::cin >> mc;
  m_MaxCost = mc;
  manifoldIntegrator->SetMaxCost(mc + 1.0);
  std::cout << " findpath in extractsurfacedisk ";
  manifoldIntegrator->FindPath();
  std::cout << " findpath in extractsurfacedisk done ";

// assign scalars to the original surface mesh

  typedef itk::SurfaceMeshCurvature<GraphSearchNodeType, GraphSearchNodeType> surfktype;
  typename surfktype::Pointer surfk = surfktype::New();

  vtkPoints*     vtkpoints = m_SurfaceMesh->GetPoints();
  int            numPoints = vtkpoints->GetNumberOfPoints();
  vtkFloatArray* param = vtkFloatArray::New();
  param->SetName("angle");
  for( int i = 0; i < numPoints; i++ )
    {
    float temp = fabs(manifoldIntegrator->GetGraphNode(i)->GetTotalCost() );
    if( temp > m_MaxCost )
      {
      temp = m_MaxCost;
      }
    param->InsertNextValue(temp * 255. / m_MaxCost);
    }

  std::cout << " extractsurfacedisk done ";
  //  m_SurfaceMesh->GetPointData()->SetScalars(param);
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool  FEMDiscConformalMap<TSurface, TImage, TDimension>
::GenerateSystemFromSurfaceMesh()
{
  if( !m_SurfaceMesh )
    {
    std::cout << " NO MESH"; return false;
    }
  std::cout << " Generate system from surface mesh " << std::endl;
  m_Smooth = m_Sigma;
  // create a material
  // Choose the material properties
  typename MaterialType::Pointer m;
  m = MaterialType::New();
  m->GN = 0;  // Global number of the material ///
  m->E = 4.0; // Young modulus -- used in the membrane ///
  m->A = 0.0; // Crossection area ///
  m->h = 0.0; // Crossection area ///
  m->I = 0.0; // Moment of inertia ///
  m->nu = 0.; // .0;    // poissons -- DONT CHOOSE 1.0!!///
  m->RhoC = 0.0;

  // Create the element type
  ElementType::Pointer e1 = ElementType::New();
  e1->m_mat = dynamic_cast<MaterialType *>( m );

  vtkPoints* vtkpoints = m_SurfaceMesh->GetPoints();
  int        numPoints = vtkpoints->GetNumberOfPoints();
  int        foundnum = 0;
  int        boundsz = 0;
  for( int i = 0; i < numPoints; i++ )
    {
    double* pt = vtkpoints->GetPoint(i);
    typename NodeType::Pointer n;
    n = new NodeType(pt[0], pt[1], pt[2]);
    if( this->InDisc(manifoldIntegrator->GetGraphNode(i) ) )
      {
      n->GN = i;
      m_Solver.node.push_back(itk::fem::FEMP<NodeType>(n) );
      foundnum++;
      }
    if( this->InBorder(manifoldIntegrator->GetGraphNode(i) ) )
      {
      boundsz++;
      }
    }

  typename NodeType::Pointer * narr = new NodeType::Pointer[numPoints];
  for( int i = 0; i < numPoints; i++ )
    {
    if( this->InDisc(manifoldIntegrator->GetGraphNode(i) ) || this->InBorder(manifoldIntegrator->GetGraphNode(i) ) )
      {
      narr[i] = m_Solver.node.Find( i );
      }
    }

  std::cout << " Found " << foundnum << " nodes " << std::endl;
  std::cout <<  " bound " << boundsz << std::endl;
  vtkCellArray* vtkcells = m_SurfaceMesh->GetPolys();

  vtkIdType     npts;
  vtkIdType*    pts;
  unsigned long i = 0;
//   unsigned long toti = vtkcells->GetNumberOfCells();
//   unsigned long rate = toti/50;
//  std::cout << " progress ";
  for( vtkcells->InitTraversal(); vtkcells->GetNextCell(npts, pts); )
    {
    //  if ( i % rate == 0 && i > rate ) std::cout << "  " <<  (float) i / (float) toti << " ";
    // turn the cell into an element
    //	  std::cout << " points ids a " << pts[0] << " b " << pts[1] << " c " << pts[2] << std::endl;
    bool                 eltok = true;
    ElementType::Pointer e;
    e = dynamic_cast<ElementType *>(e1->Clone() );

    if( this->InDisc(manifoldIntegrator->GetGraphNode(pts[0]) ) )
      {
      e->SetNode(2, narr[pts[0]]);
      }
    //		  e->SetNode(2,m_Solver.node.Find( pts[0] ));
    else
      {
      eltok = false;
      }

    if( this->InDisc(manifoldIntegrator->GetGraphNode(pts[1]) ) )
      {
      e->SetNode(1, narr[pts[1]]);
      }
    //		  e->SetNode(1,m_Solver.node.Find( pts[1] ));
    else
      {
      eltok = false;
      }

    if( this->InDisc(manifoldIntegrator->GetGraphNode(pts[2]) ) )
      {
      e->SetNode(0, narr[pts[2]]);
      }
    //		  e->SetNode(0,m_Solver.node.Find( pts[2] ));
    else
      {
      eltok = false;
      }

    if( eltok )
      {
      e->GN = i;
      if( pts[0] == this->m_SourceNodeNumber )
        {
        std::cout << "choosing pole " << i << std::endl;
        m_PoleElementsGN[0] = i;
        }
      m_Solver.el.push_back(itk::fem::FEMP<itk::fem::Element>(e) );
      i++;
      } // else std::cout <<" cannot find elt " << std::endl;
    }

  std::cout << " DONE: NUMBER OF CELLS " << i << std::endl;

  return true;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::FixPointsBeyondDisc()
{
//  if (dim > 3 || dim < -3 || dim == 0 ) return;

  unsigned int dofsperelt;

  itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin();
  unsigned int                           Nnodes = (*elt)->GetNumberOfNodes();
  unsigned int                           dofs = (*elt)->GetNumberOfDegreesOfFreedomPerNode();

  dofsperelt = dofs * Nnodes;

  int fixct = 0;
  int eltct = 0;
  while( elt != m_Solver.el.end() )
    {
    for( unsigned int i = 0; i < dofs; i++ )
      {
      int nodeid = (*elt)->GetNode(i)->GN;
      //      if( manifoldIntegrator->GetGraphNode(nodeid)->GetTotalCost() > 222 )
      if( this->InBorder(manifoldIntegrator->GetGraphNode(nodeid) ) )
        {
        itk::fem::LoadBC::Pointer l1;
        l1 = itk::fem::LoadBC::New();
        l1->m_element = (*elt);
        l1->m_dof = i;
        //	      l1->m_value=vnl_vector<double>(1,0.0); // for exp
        l1->m_value = vnl_vector<double>(1, 0.0);   // for direct rad
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );
        if( i == 0 )
          {
          fixct++;
          }
        //	      std::cout << " bord coord "<< manifoldIntegrator->GetGraphNode(nodeid)->GetLocation() << std::endl;
        }
      }
    elt++;
    eltct++;
    }

  std::cout << " Fixed elt number " << fixct << " of " << eltct << std::endl;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ApplyRealForces()
{
  itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[0]];
  // load node 0
    {
    itk::fem::LoadNode::Pointer ln1 = itk::fem::LoadNode::New();
    ln1->m_pt = 0;
    ln1->F.set_size(1);
    ln1->F.fill(-4.0 * m_Pi);
    ln1->m_element = (e);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln1) );

    itk::fem::LoadNode::Pointer ln2 = itk::fem::LoadNode::New();
    ln2->m_pt = 1;
    ln2->F.set_size(1);
    ln2->F.fill(-4.0 * m_Pi);
    ln2->m_element = (e);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln2) );

    itk::fem::LoadNode::Pointer ln3 = itk::fem::LoadNode::New();
    ln3->m_pt = 2;
    ln3->F.set_size(1);
    ln3->F.fill(-4.0 * m_Pi);
    ln3->m_element = (e);
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*ln3) );

    /*
    itk::fem::LoadBC::Pointer l1;
    l1=itk::fem::LoadBC::New();
    l1->m_element=(e);
    l1->m_dof=0;
    l1->m_value=vnl_vector<double>(1,-2./3.);// for exp
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );

    itk::fem::LoadBC::Pointer l2;
    l2=itk::fem::LoadBC::New();
    l2->m_element=(e);
    l2->m_dof=1;
    l2->m_value=vnl_vector<double>(1,1./3.);// for exp
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l2) );

    itk::fem::LoadBC::Pointer l3;
    l3=itk::fem::LoadBC::New();
    l3->m_element=(e);
    l3->m_dof=2;
    l3->m_value=vnl_vector<double>(1,1./3.);// for exp
    m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l3) );
    */
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ApplyImaginaryForces()
{
  itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[0]];

  itk::fem::LoadBC::Pointer l1;

  l1 = itk::fem::LoadBC::New();
  l1->m_element = (e);
  l1->m_dof = 0;
  l1->m_value = vnl_vector<double>(1, -1. / 3.); // for exp
  m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );

  itk::fem::LoadBC::Pointer l2;
  l2 = itk::fem::LoadBC::New();
  l2->m_element = (e);
  l2->m_dof = 1;
  l2->m_value = vnl_vector<double>(1, -1. / 3.); // for exp
  m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l2) );

  itk::fem::LoadBC::Pointer l3;
  l3 = itk::fem::LoadBC::New();
  l3->m_element = (e);
  l3->m_dof = 2;
  l3->m_value = vnl_vector<double>(1, 2. / 3.); // for exp
  m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l3) );

  /*
  unsigned int dofsperelt;

  itk::fem::Element::ArrayType::iterator elt=m_Solver.el.begin();
  unsigned int Nnodes= (*elt)->GetNumberOfNodes();
  unsigned int dofs=(*elt)->GetNumberOfDegreesOfFreedomPerNode();
  dofsperelt=dofs*Nnodes;

  while(elt!=m_Solver.el.end())
    {
      for (int i=0; i<dofs; i++)
  {
    int nodeid=(*elt)->GetNode(i)->GN;
    if( this->InBorder(manifoldIntegrator->GetGraphNode(nodeid)) )
      {

        itk::fem::LoadBC::Pointer l1;
        l1=itk::fem::LoadBC::New();
        l1->m_element=(*elt);
        l1->m_dof=0;
        l1->m_value=vnl_vector<double>(1,0.0);// for exp
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );

        itk::fem::LoadBC::Pointer l2;
        l2=itk::fem::LoadBC::New();
        l2->m_element=(*elt);
        l2->m_dof=1;
        l2->m_value=vnl_vector<double>(1,2.0*m_Pi);// for exp
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l2) );

        return;
      }
  }
      elt++;
    }
  */
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>::MakeFlatImage()
{
// first declare the flat image
  typename FlatImageType::RegionType region;
  typename FlatImageType::SizeType size;
  int sz = 256;
  size[0] = sz;
  size[1] = sz;
  region.SetSize( size );
  m_FlatImage = FlatImageType::New();
  m_FlatImage->SetRegions( region );
  m_FlatImage->Allocate();
  typename FlatImageType::IndexType index;

  std::cout << " Making flat image " << std::endl;
  int maxits = 100;
  for( int its = 0; its <= maxits; its++ )
    {
    for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
      {
      float temp = 255.0 - manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(3); // curvature
      //	  float temp=255.0*manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(2); // extrinsic dist
      index[0] =
        (long int)(0.5 + ( 1.0 + manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(0) ) * (float)(sz - 1) / 2.);
      index[1] =
        (long int)(0.5 + ( 1.0 + manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(1) ) * (float)(sz - 1) / 2.);
      // std::cout << " ind " << index << std::endl;
      m_FlatImage->SetPixel(index, temp);
      }
    typedef itk::DiscreteGaussianImageFilter<FlatImageType, FlatImageType> dgf;
    typename dgf::Pointer filter = dgf::New();
    filter->SetVariance(1.0);
    filter->SetUseImageSpacingOff();
    filter->SetMaximumError(.01f);
    filter->SetInput(m_FlatImage);
    filter->Update();
    m_FlatImage = filter->GetOutput();

    if( its < maxits )
      {
      int                                              center = (int)sz / 2;
      itk::ImageRegionIteratorWithIndex<FlatImageType> it( m_FlatImage, region );
      it.GoToBegin();
      typename FlatImageType::IndexType index;
      while( !it.IsAtEnd() )
        {
        index = it.GetIndex();
        float x = (float)index[0] - (float)sz / 2.0;
        float y = (float)index[1] - (float)sz / 2.0;
        float dist = sqrt(x * x + y * y);
        //	    std::cout << "center " << center <<  " index " << index << "dist " << dist ;
        if( dist > center )
          {
          it.Set( 0.0 );
          }
        ++it;
        }
      }
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>::BuildOutputMeshes(float tval)
{
  typedef GraphSearchNodeType::NodeLocationType loctype;
  // Get the number of points in the mesh
  int numPoints = m_Solver.node.size();

  // Create vtk polydata
  vtkPolyData* polydata1 = vtkPolyData::New();
  vtkPolyData* polydata2 = vtkPolyData::New();

  // Create the vtkPoints object and set the number of points
  vtkPoints* vpoints1 = vtkPoints::New();
  vtkPoints* vpoints2 = vtkPoints::New();
  vpoints1->SetNumberOfPoints(numPoints);
  vpoints2->SetNumberOfPoints(numPoints);

  std::cout << " start pts ";
  int idx = 0;

  vtkFloatArray* param = vtkFloatArray::New();
  param->SetName("curvature");

  vtkFloatArray* paramAngle = vtkFloatArray::New();
  paramAngle->SetName("angle");

  vtkFloatArray* paramDistance = vtkFloatArray::New();
  paramDistance->SetName("distance");

  vtkIdTypeArray* paramPoints = vtkIdTypeArray::New();
  paramPoints->SetName("points");
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    loctype loc = manifoldIntegrator->GetGraphNode( (*n)->GN)->GetLocation();
    float   pt1[3];
    pt1[0] = loc[0];
    pt1[1] = loc[1];
    pt1[2] = loc[2];
    float pt2[3];
    pt2[0] = manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(0) * 100. * (1.0 - tval) + tval * loc[0];
    pt2[1] = manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(1) * 100. * (1.0 - tval) + tval * loc[1];
    pt2[2] = 1. * (1.0 - tval) + tval * loc[2];
    vpoints1->SetPoint(idx, pt1);
    vpoints2->SetPoint(idx, pt2);

    //	float temp=(float)manifoldIntegrator->GetGraphNode((*n)->GN)->GetTotalCost();
    //	float temp=(float)manifoldIntegrator->GetGraphNode((*n)->GN)->GetValue(2);
    //	if (temp > m_MaxCost) temp=m_MaxCost;
    float temp = manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(3);        // for curvature
    float temp2 = manifoldIntegrator->GetGraphNode( (*n)->GN)->GetValue(2) * 255; // for length
    //    temp=m_RealSolution[(*n)->GN]*255.0;
    param->InsertNextValue(temp);
    paramDistance->InsertNextValue(temp2);
    paramPoints->InsertNextValue( (*n)->GN);
    paramAngle->InsertNextValue(temp); // curvature
    //    param->InsertNextValue(temp*255./m_MaxCost);

    (*n)->GN = idx;
    idx++;
    }

  std::cout << " done with pts " << std::endl;
  vtkCellArray* tris1 = vtkCellArray::New();
  vtkCellArray* tris2 = vtkCellArray::New();

  std::cout << " start with tris " << std::endl;
  for( ::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); n++ )
    {
    tris1->InsertNextCell(3);
    tris2->InsertNextCell(3);
    for( unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++ )
      {
      tris1->InsertCellPoint( (*n)->GetNode(i)->GN);
      tris2->InsertCellPoint( (*n)->GetNode(i)->GN);
      }
    }
  std::cout << " done with tris " << std::endl;
  // Assign points and cells
  polydata1->SetPoints(vpoints1);
  polydata2->SetPoints(vpoints2);
  polydata1->SetPolys(tris1);
  polydata2->SetPolys(tris2);
  polydata1->GetPointData()->SetScalars(param);
  polydata2->GetPointData()->SetScalars(param);
  polydata1->GetPointData()->AddArray(paramAngle);
  polydata2->GetPointData()->AddArray(paramAngle);
  polydata1->GetPointData()->AddArray(paramDistance);
  polydata2->GetPointData()->AddArray(paramDistance);
  polydata1->GetPointData()->AddArray(paramPoints);
  polydata2->GetPointData()->AddArray(paramPoints);

  m_ExtractedSurfaceMesh = polydata1; //

  vtkDelaunay2D* delny2 = vtkDelaunay2D::New();
  delny2->SetInput(polydata2);
  m_DiskSurfaceMesh = delny2->GetOutput();
  // m_DiskSurfaceMesh=polydata2;

  return;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ConformalMap()
{
  m_Solver.load.clear();
  m_Solver.node.clear();
  m_Solver.el.clear();

  /**
   * Open the file and assign it to stream object f
   */
  if( m_ReadFromFile )
    {
    const char* filename = m_ParameterFileName.c_str();
    std::cout << "Reading FEM problem from file: " << std::string(filename) << "\n";
    std::ifstream f;
    f.open(filename);
    if( !f )
      {
      std::cout << "File " << filename << " not found!\n";
      return;
      }

    try
      {
      m_Solver.Read(f);
      }
    catch( ::itk::fem::FEMException e )
      {
      std::cout << "Error reading FEM problem: " << filename << "!\n";
      e.Print(std::cout);
      return;
      }

    f.close();
    }
  else if( this->GenerateSystemFromSurfaceMesh() )
    {
    ;
    }
  else
    {
    return;
    }

  /**
   * Assign a unique id (global freedom number - GFN)
   * to every degree of freedom (DOF) in a system.
   */
  m_Solver.GenerateGFN();
  m_ImagSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom() );
  m_RealSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom() );
  m_Radius.set_size(m_Solver.GetNumberOfDegreesOfFreedom() );
  m_RealSolution.fill(0);
  m_ImagSolution.fill(0);
  m_Radius.fill(0);

  m_Debug = false;
  if( m_Debug )
    {
    for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
      {
      std::cout << "Node#: " << (*n)->GN << ": ";
      std::cout << " coord " << (*n)->GetCoordinates()
                << " coord2 " << manifoldIntegrator->GetGraphNode( (*n)->GN)->GetLocation() << std::endl;
      }
    for( ::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); n++ )
      {
      std::cout << "Elt#: " << (*n)->GN << ": has " << (*n)->GetNumberOfNodes() << " nodes ";
      for( unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++ )
        {
        std::cout << " coord " << (*n)->GetNode(i)->GetCoordinates() << std::endl;
        }
      }
    }

  unsigned int maxits = m_Solver.GetNumberOfDegreesOfFreedom(); // should be > twice ndofs
  // if (m_Debug)
  std::cout << " ndof " << maxits << std::endl;
  itpackWrapper.SetMaximumNumberIterations(maxits * 3);
  itpackWrapper.SetTolerance(1.e-1);
  itpackWrapper.SuccessiveOverrelaxation();
  // itpackWrapper.JacobianConjugateGradient();
  itpackWrapper.SetMaximumNonZeroValuesInMatrix(maxits * 10);
  m_Solver.SetLinearSystemWrapper(&itpackWrapper);

  this->FixPointsBeyondDisc();
  m_Solver.AssembleK();
  m_Solver.DecomposeK();
  this->ApplyRealForces();
  std::cout << " appl force ";
  m_Solver.AssembleF();
  //  for (int i=0; i<maxits; i++) if (m_Solver.GetVectorValue(i) != 0) m_Solver.SetVectorValue(i,1.0);
  std::cout << " b solve ";
  m_Solver.Solve();
  std::cout << " e solve ";
  m_Solver.UpdateDisplacements(); // copies solution to nodes
  unsigned long ct  =  0;
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    for( unsigned int d = 0, dof; (dof = (*n)->GetDegreeOfFreedom(d) ) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
         d++ )
      {
      m_RealSolution[dof] = m_Solver.GetSolution(dof);

      // FIXME m_MaxCost = actual radius value - 1  .....

      m_Radius[dof] = manifoldIntegrator->GetGraphNode( (*n)->GN )->GetTotalCost() / m_MaxCost;
//    double diff = m_Radius[dof] - m_MaxCost;
      //	  std::cout << " diff " << diff << "  " << m_Smooth << std::endl;
      // if (diff > 0) diff = 0 ;
      if( ct % 100 == 0 )
        {
        std::cout << " mrdof " <<  m_RealSolution[dof]  <<  " rad " << m_Radius[dof] << std::endl;
        }
      // m_RealSolution[dof]=m_Radius[dof];//exp(diff/m_Smooth);
      }
    ct++;
    }

  this->ConformalMap2();
  // this->ConformalMap3();
  // this->MapToSquare();
  this->ConjugateHarmonic();
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::MeasureLengthDistortion()
{
  // now measure the length distortion of the given solution
  manifoldIntegrator->EmptyQ();
  manifoldIntegrator->SetSearchFinished( false );
  manifoldIntegrator->m_PureDist = false;
  for( int i = 0; i < manifoldIntegrator->GetGraphSize(); i++ )
    {
    manifoldIntegrator->GetGraphNode(i)->SetTotalCost(vnl_huge_val(manifoldIntegrator->GetMaxCost() ) );
    manifoldIntegrator->GetGraphNode(i)->SetUnVisited();
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 1);
    manifoldIntegrator->GetGraphNode(i)->SetValue(0.0, 2);
    manifoldIntegrator->GetGraphNode(i)->SetPredecessor(NULL);
    }
  manifoldIntegrator->SetSource(manifoldIntegrator->GetGraphNode(this->m_SourceNodeNumber) );
  manifoldIntegrator->InitializeQueue();
  float mchere = 1.2;
  manifoldIntegrator->SetMaxCost(mchere);
  // here we want to find another path between the source and sink, next to the path
  // found above
  manifoldIntegrator->FindPath();

  // backtrack everywhere to set up forweard tracking
  float maxmanifolddist = 0;
  float meanmanifolddist = 0;
  float distDistortion = 0;

  unsigned int ct = 0;
  for( int i = 0; i < manifoldIntegrator->GetGraphSize(); i++ )
    {
    if( manifoldIntegrator->GetGraphNode(i) )
      {
      if( manifoldIntegrator->GetGraphNode(i)->GetTotalCost() <= manifoldIntegrator->GetMaxCost() )
        {
        float ttt = manifoldIntegrator->GetGraphNode(i)->GetValue(2);
        meanmanifolddist += ttt;
        if( ttt > maxmanifolddist )
          {
          maxmanifolddist = ttt;
          }
        ct++;
        }
      }
    }
  meanmanifolddist /= (float)ct;
  ct = 0;
  for( int i = 0; i < manifoldIntegrator->GetGraphSize(); i++ )
    {
    float manifolddist;
    float rad;
    if( manifoldIntegrator->GetGraphNode(i) )
      {
      if( manifoldIntegrator->GetGraphNode(i)->GetTotalCost() < manifoldIntegrator->GetMaxCost() )
        {
        rad = manifoldIntegrator->GetGraphNode(i)->GetValue(0);
        manifolddist = manifoldIntegrator->GetGraphNode(i)->GetValue(2) / maxmanifolddist;
        manifoldIntegrator->GetGraphNode(i)->SetValue(manifolddist, 2);
        distDistortion += fabs(rad - manifolddist);
        ct++;
        }
      }
    }
  std::cout << "  distDistortion/ct " << distDistortion / (float)ct << " maxmfd " << maxmanifolddist << std::endl;
  return;

  /*    m_ExtractedSurfaceMesh=polydata1;//
  //   m_DiskSurfaceMesh=delny2->GetOutput();
  float total=0.0;
  for (int i=0; i<m_SurfaceMesh->GetPoints()->GetNumberOfPoints(); i++)
  {
  float length1=0;
  float length2=0;
  float pt1[3];

  for (int ne=0; ne<manifoldIntegrator->GetGraphNode(i)->GetNumberOfNeighbors(); ne++)
  {
  }

  }
*/
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::FixThetaAlongBorder()
{
  /*
  unsigned int ct = 0;
  for( ::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n!=m_Solver.el.end(); n++)
    {
      for (int i=0; i<(*n)->GetNumberOfNodes(); i++)
  {
    unsigned int dof=(*n)->GetNode(i)->GetDegreeOfFreedom(0);
    if ( m_Radius[(*n)->GetNode(i)->GN] > 0 )
      {
        itk::fem::LoadBC::Pointer l1;
        l1=itk::fem::LoadBC::New();
        l1->m_element=(*n);
        l1->m_dof=0;
        l1->m_value=vnl_vector<double>(1, m_ImagSolution[dof]);
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );
        ct++;
      }
  }
    }
    std::cout << "MAP3 " << ct << std::endl;
*/
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::FixPointsAlongRadialLine()
{
  int bordernodeindex = 0;
  // first find a point that is in the border and has max dist
  float maxdist = 0;

  for( int i = 0; i < manifoldIntegrator->GetGraphSize(); i++ )
    {
    float cost = manifoldIntegrator->GetGraphNode(i)->GetTotalCost();
    if( this->InBorder(manifoldIntegrator->GetGraphNode(i) )  &&
        cost >= maxdist && manifoldIntegrator->GetGraphNode(i)->GetPredecessor() )
      {
      bordernodeindex = i;
      maxdist = cost;
      }
    }

  // now backtrack from here
  manifoldIntegrator->BackTrack( manifoldIntegrator->GetGraphNode(bordernodeindex) );
  int pathsz = manifoldIntegrator->GetPathSize();
  for( int j = pathsz - 1; j >= 0; j-- )
    { // just use some specific number to label these nodes
    manifoldIntegrator->GetPathAtIndex(j)->SetValue(-99.0, 2);
    if( j != (pathsz - 1) && j != 0 )
      {
      manifoldIntegrator->GetPathAtIndex(j)->SetUnVisitable();
      }
    }
  float d1 = manifoldIntegrator->GetPathAtIndex(0)->GetTotalCost();
  for( ::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); n++ )
    {
    for( unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++ )
      {
      if( manifoldIntegrator->GetGraphNode( (*n)->GetNode(i)->GN)->GetValue(2) == -99.0 )
        {
        manifoldIntegrator->GetGraphNode( (*n)->GetNode(i)->GN)->SetUnVisitable();
        itk::fem::LoadBC::Pointer l1;
        l1 = itk::fem::LoadBC::New();
        l1->m_element = (*n);
        l1->m_dof = 0;
        l1->m_value = vnl_vector<double>(1, 2.0 * m_Pi);
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );
        }
      }
    }

  // now basically repeat the search and the above code as well
  manifoldIntegrator->EmptyQ();
  for( int i = 0; i < manifoldIntegrator->GetGraphSize(); i++ )
    {
    manifoldIntegrator->GetGraphNode(i)->SetTotalCost(vnl_huge_val(manifoldIntegrator->GetMaxCost() ) );
    if( !manifoldIntegrator->GetGraphNode(i)->GetUnVisitable() )
      {
      manifoldIntegrator->GetGraphNode(i)->SetUnVisited();
      manifoldIntegrator->GetGraphNode(i)->SetPredecessor(NULL);
      }
    }
  manifoldIntegrator->SetSearchFinished( false );
  manifoldIntegrator->SetSource(manifoldIntegrator->GetGraphNode(this->m_SourceNodeNumber) );
  manifoldIntegrator->SetSink(manifoldIntegrator->GetGraphNode(bordernodeindex) );
  manifoldIntegrator->InitializeQueue();
  manifoldIntegrator->SetMaxCost(m_MaxCost * 2.0);
  // here we want to find another path between the source and sink, next to the path
  // found above
  std::cout << " findpath in fix points along radial ";
  manifoldIntegrator->FindPath();
  std::cout << " findpath in fix points along radial done ";

  // now backtrack from here
  if( !manifoldIntegrator->GetGraphNode(bordernodeindex)->GetPredecessor() )
    {
    bordernodeindex = 0;
    for( int ne = 0; ne < manifoldIntegrator->GetGraphNode(bordernodeindex)->GetNumberOfNeighbors(); ne++ )
      {
      if( manifoldIntegrator->GetGraphNode(bordernodeindex)->GetNeighbor(ne)->GetPredecessor() )
        {
        std::cout << " Found new bordernodeindex ";
        bordernodeindex = manifoldIntegrator->GetGraphNode(bordernodeindex)->GetNeighbor(ne)->GetIdentity();
        ne = manifoldIntegrator->GetGraphNode(bordernodeindex)->GetNumberOfNeighbors();
        }
      }
    }

  manifoldIntegrator->BackTrack( manifoldIntegrator->GetGraphNode(bordernodeindex) );
  int pathsz2 = manifoldIntegrator->GetPathSize();
  for( int j = pathsz2 - 1; j >= 0; j-- )
    { // just use some specific number to label these nodes
    manifoldIntegrator->GetPathAtIndex(j)->SetValue(-999.0, 2);
    }
  float d2 = manifoldIntegrator->GetPathAtIndex(0)->GetTotalCost();
  std::cout << " pathsz1 " << pathsz  << " len  " << d1
            << " pathsz2 " << pathsz2 << " len2 " << d2 << std::endl;
  for( ::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); n++ )
    {
    for( unsigned int i = 0; i < (*n)->GetNumberOfNodes(); i++ )
      {
      if( manifoldIntegrator->GetGraphNode( (*n)->GetNode(i)->GN)->GetValue(2) == -999.0
          && !manifoldIntegrator->GetGraphNode( (*n)->GetNode(i)->GN)->GetUnVisitable() )
        {
        itk::fem::LoadBC::Pointer l1;
        l1 = itk::fem::LoadBC::New();
        l1->m_element = (*n);
        l1->m_dof = 0;
        double epsilon = 0.01;
        l1->m_value = vnl_vector<double>(1, epsilon);
        m_Solver.load.push_back( itk::fem::FEMP<itk::fem::Load>(&*l1) );
        }
      }
    }

  //   std::cout << " Pathsize  2 ";
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ConformalMap2()
{
//   unsigned int maxits=m_Solver.GetNumberOfDegreesOfFreedom(); // should be > twice ndofs
  m_Solver.load.clear();
  //    m_Solver.node.clear();
  //     m_Solver.el.clear();
  this->FixPointsAlongRadialLine();
  m_Solver.AssembleK(); // need to reassemble b/c LoadBC's affect K
  // this->ApplyImaginaryForces();
  m_Solver.AssembleF();
  m_Solver.Solve();
  m_Solver.UpdateDisplacements(); // copies solution to nodes
  unsigned long ct = 0;
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    for( unsigned int d = 0, dof; (dof = (*n)->GetDegreeOfFreedom(d) ) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
         d++ )
      {
      m_ImagSolution[dof] = m_Solver.GetSolution(dof);
      if( ct % 100 == 0 )
        {
        std::cout << " midof " <<  m_ImagSolution[dof]  <<  " dof " << dof << std::endl;
        }
      }
    ct++;
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ConformalMap3()
{
  unsigned int maxits = m_Solver.GetNumberOfDegreesOfFreedom(); // should be > twice ndofs

  m_Solver.load.clear();
  this->FixThetaAlongBorder();
  m_Solver.AssembleK(); // need to reassemble b/c LoadBC's affect K
  m_Solver.AssembleF();
  m_Solver.Solve();
  m_Solver.UpdateDisplacements(); // copies solution to nodes
  unsigned long ct = 0;
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    for( unsigned int d = 0, dof; (dof = (*n)->GetDegreeOfFreedom(d) ) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
         d++ )
      {
      m_ImagSolution[dof] = m_Solver.GetSolution(dof);
      if( ct % 100 == 0 )
        {
        std::cout << " midof " <<  m_ImagSolution[dof]  <<  " dof " << dof << std::endl;
        }
      }
    ct++;
    }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::MapToSquare()
{
  // now backtrack from all pts to the source - this sets up the ancestors
  float minintval = 9.e9;
  float maxintval = -9.e9;

  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    int   dof = (*n)->GetDegreeOfFreedom(0);
    float U = m_RealSolution[dof];
    manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(U, 0);
    }
  this->MeasureLengthDistortion();

  int   tct = 0;
  float minu = 9.e9;
  float minv = 9.e9;
  float maxu = 0;
  float maxv = 0;
  typedef itk::SurfaceMeshCurvature<GraphSearchNodeType, GraphSearchNodeType> surfktype;
  typename surfktype::Pointer surfk = surfktype::New();

  float         totallengthdistortion = 0.;
  unsigned long ct  = 0;
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    int   dof = (*n)->GetDegreeOfFreedom(0);
    float U = m_RealSolution[dof];
    float V = m_ImagSolution[dof];  // manifoldIntegrator->GetGraphNode( (*n)->GN )->GetValue(2);
    if( U < minu )
      {
      minu = U;
      }
    if( V < minv )
      {
      minv = V;
      }
    if( U > maxu )
      {
      maxu = U;
      }
    if( V > maxv )
      {
      maxv = V;
      }
    }
  std::cout << " MINU " << minu << " MINV " << minv << std::endl;
  std::cout << " MaxU " << maxu << " MaxV " << maxv << std::endl;
  float urange = maxu - minu;
  float vrange = maxv - minv;
  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
    {
    ct++;
    int   dof = (*n)->GetDegreeOfFreedom(0);
    float U = m_RealSolution[dof];
    float V = m_ImagSolution[dof];  // manifoldIntegrator->GetGraphNode( (*n)->GN )->GetValue(2);

    float x;
    if( x > 0 )
      {
      x = (maxu - U) / urange;
      }
    else
      {
      x = (minu - U) / urange;
      }
    float y;
    if( y > 0 )
      {
      y = (maxv - V) / vrange;
      }
    else
      {
      y = (minv - V) / vrange;
      }
    manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(x, 0);
    manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(y, 1);
    }
  //  this->MakeFlatImage();
  //  this->BuildOutputMeshes(0.5);

  return;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void  FEMDiscConformalMap<TSurface, TImage, TDimension>
::ConjugateHarmonic()
{
  // now backtrack from all pts to the source - this sets up the ancestors
//   float minintval=9.e9;
//   float maxintval=-9.e9;
/*
for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n!=m_Solver.node.end(); n++)
  {
    int dof=(*n)->GetDegreeOfFreedom(0);
    float U = m_RealSolution[dof];
    manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(U,0);
  }
this->MeasureLengthDistortion();
//
  for (int i=0; i<manifoldIntegrator->GetGraphSize(); i++)
  {
  if (this->InDisc(manifoldIntegrator->GetGraphNode(i) ))
  {
  // backtracking gives us a curve along which to integrate
    manifoldIntegrator->BackTrack( manifoldIntegrator->GetGraphNode(i) );
  // then perform integration to the given point.  this set its V value.
  int pathsz=manifoldIntegrator->GetPathSize();
  float intval=0.0;
  typedef GraphSearchNodeType::NodeLocationType loctype;
  for (int j=pathsz-2; j>=0; j--)
  {
    loctype delt=manifoldIntegrator->GetPathAtIndex(j+1)->GetLocation()-
      manifoldIntegrator->GetPathAtIndex(j)->GetLocation();
    float dstarU=manifoldIntegrator->dstarUestimate(
      manifoldIntegrator->GetPathAtIndex(j))*delt.magnitude();
    if (j==0) dstarU*=0.5;
//		  dstarU=fabs(dstarU)*(-1.0);
      intval+=dstarU;
    manifoldIntegrator->GetPathAtIndex(j)->SetValue(intval,2);
  }
  if (intval < minintval) minintval=intval;
  if (intval > maxintval) maxintval=intval;
  }
}
std::cout << " MAX int VAL " << maxintval << std::endl;
std::cout << " MIN int VAL " << minintval << std::endl;
*/

// now print out the conf coords U+iV
  std::cout << " Conformal coordinates " << std::endl;
//   int tct = 0;
  float minu = 9.e9;
  float minv = 9.e9;
  float maxu = 0;
  float maxv = 0;

  typedef itk::SurfaceMeshCurvature<GraphSearchNodeType, GraphSearchNodeType> surfktype;
  typename surfktype::Pointer surfk = surfktype::New();

  //    for (int DER=0; DER<4; DER++)
    {
    float dorad = 0.; // smaller more length pres
    std::cout << " input dorad "; std::cin >> dorad;

    float         totallengthdistortion = 0.;
    unsigned long ct  = 0;
    for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); n++ )
      {
      ct++;
      int   dof = (*n)->GetDegreeOfFreedom(0);
      float U = m_RealSolution[dof];
      float V = m_ImagSolution[dof];           // manifoldIntegrator->GetGraphNode( (*n)->GN )->GetValue(2);
      float expu = exp( (double)U / m_Smooth); // radius
      float dp = m_Radius[dof];

      float r = (1. - dorad) * dp + dorad * expu;

      //	  std::cout << " rad " << r << " ct  " << ct << std::endl;

      float rat1 = dp / r; if( rat1 > 10 )
        {
        rat1 = 10;
        }
      float rat2 = r / dp; if( rat2 > 10 )
        {
        rat2 = 10;
        }
      if( rat1 > 1 )
        {
        totallengthdistortion += rat1;
        }
      if( rat2 > 1 )
        {
        totallengthdistortion += rat2;
        }

      if( r > 1.0 )
        {
        r = 1.0;
        }
      if( r < 0.0 )
        {
        r = 0.0;
        }

      float         expv = V; // exp(V/sig)*(2.0*m_Pi);//+fabs(minintval);//exp(V);
      unsigned long mult = (unsigned long)(expv / (2.0 * m_Pi) );
      float         theta = expv - ( (float)(mult) ) * (2.0 * m_Pi); // theta

      //      if (ct % 100 == 0) std::cout << " R "  << r << " expu " <<  expu << " Theta " << theta << " tld " <<
      // totallengthdistortion/(float)ct << " V " << V << " mult " << mult <<  std::endl;

      float x = r * cos(theta); // expu/denom;
      float y = r * sin(theta); // expv/denom;
//       float denom=sqrt(x*x+y*y);

      bool makesquare = false;
      //      if (theta <= m_Pi) makesquare = true;
      if( makesquare )
        {
        std::cout << " old x " << x << " old y " << y << " th " << theta << std::endl;
        // use linear interp
        if( theta >= 0 && theta < m_Pi / 4.0 )
          {
          y = tan(theta);
          x = (sqrt(2.) - sqrt(y * y + 1.) ) * r;
          }
        else if( theta >= m_Pi / 4.0 && theta < m_Pi / 2.0 )
          {
          x = tan(theta - m_Pi / 4.);
          y = (sqrt(2.) - sqrt(x * x + 1.) ) * r;
          }
        else if( theta >= m_Pi / 2.0 && theta < m_Pi * 3. / 4.0 )
          {
          x = tan(theta - m_Pi / 2.) * (-1.);
          y = (sqrt(2.) - sqrt(x * x + 1.) ) * r;
          }
        else if( theta >= m_Pi * 3. / 4.0 && theta < m_Pi )
          {
          y = tan(theta - m_Pi / 3.);
          x = -1.0 * (sqrt(2.) - sqrt(y * y + 1.) ) * r;
          }
        std::cout << " NEW x " << x << " NEW y " << y << std::endl;
        }

      manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(x, 0);
      manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(y, 1);
      manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(r, 2);
      //      manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue(r,3);

      if( r < minu )
        {
        minu = r;
        }
      if( theta < minv )
        {
        minv = theta;
        }
      if( r > maxu )
        {
        maxu = r;
        }
      if( theta > maxv )
        {
        maxv = theta;
        }
      }
    m_MaxCost = maxu;
    std::cout << " MINU " << minu << " MINV " << minv << std::endl;
    std::cout << " MaxU " << maxu << " MaxV " << maxv <<  " m_Smooth " << m_Smooth << std::endl;

    std::cout <<  "  totallengthdistortion " << totallengthdistortion / (float)ct << std::endl;
    }

  //  this->MakeFlatImage();
  this->BuildOutputMeshes();

  return;
}

// now just loop through the graph and set the value correctly
/*
  manifoldIntegrator->ResetMaxCost();
  manifoldIntegrator->m_PureDist=false;
  manifoldIntegrator->EmptyQ();

  for (int i=0; i<manifoldIntegrator->GetGraphSize(); i++)
  {
    manifoldIntegrator->GetGraphNode(i)->SetValue(manifoldIntegrator->GetMaxCost());
    if (manifoldIntegrator->GetGraphNode(i)->GetPredecessor())
      manifoldIntegrator->GetGraphNode(i)->SetTotalCost(manifoldIntegrator->GetMaxCost());
    else
      manifoldIntegrator->GetGraphNode(i)->SetTotalCost(0);
    manifoldIntegrator->GetGraphNode(i)->SetUnVisited();
    manifoldIntegrator->GetGraphNode(i)->SetPredecessor(NULL);
    //      if (manifoldIntegrator->GetGraphNode(i)->GetPredecessor() )
      //	    manifoldIntegrator->GetGraphNode(i)->SetValue(   );
      //	  else manifoldIntegrator->GetGraphNode(i)->SetValue(manifoldIntegrator->GetMaxCost());
  }

  for( ::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n!=m_Solver.node.end(); n++)
    {
    int dof=(*n)->GetDegreeOfFreedom(0);
      manifoldIntegrator->GetGraphNode( (*n)->GN )->SetValue( m_RealSolution[dof] );
    }

    for (int i=0; i<m_DiscBoundaryList.size(); i++)
  {
    manifoldIntegrator->SetSource(manifoldIntegrator->GetGraphNode(m_DiscBoundaryList[i]));
  }
    manifoldIntegrator->EmptyQ();
    manifoldIntegrator->SetSearchFinished( false );
    manifoldIntegrator->InitializeQueue();
    manifoldIntegrator->FindPath();

  for (int i=0; i<manifoldIntegrator->GetGraphSize(); i++)
  {
    float val=manifoldIntegrator->GetGraphNode(i)->GetTotalCost();
      if  ( manifoldIntegrator->GetGraphNode(i)->GetPredecessor()
      && val < 444 && val > 0.2)
      std::cout << " intval " << val << std::endl;
  }
*/
} // namespace itk

#endif
