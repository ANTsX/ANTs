/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detailm_Solver.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _FEMConformalMap_hxx
#define _FEMConformalMap_hxx

#include "vtkFeatureEdges.h"
#include "vtkPointLocator.h"
#include "vtkCellLocator.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"

#include <vcl_compiler.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <vnl/vnl_real_polynomial.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkMath.h>

#include "itkFEMConformalMap.h"
#include "itkFEMLoadNode.h"

namespace itk
{
vtkPolyData *
vtkGenSmoothMesh(vtkPolyData * vds)
{
  float                     fai = 45;
  vtkSmoothPolyDataFilter * polySmoother = vtkSmoothPolyDataFilter::New();

  //    polySmoother->SetInput(decimator->GetOutput());
  polySmoother->SetInput(vds);
  polySmoother->SetNumberOfIterations(100);
  polySmoother->SetRelaxationFactor(0.01);
  polySmoother->SetFeatureAngle(fai);
  polySmoother->FeatureEdgeSmoothingOff();
  polySmoother->BoundarySmoothingOff();
  polySmoother->Update();

  return polySmoother->GetOutput();
}

template <typename TSurface, typename TImage, unsigned int TDimension>
FEMConformalMap<TSurface, TImage, TDimension>::FEMConformalMap()
{
  m_Sigma = 2.0e-4;
  m_ParameterFileName = "";
  m_NorthPole = 0;
  m_SouthPole = 1;

  m_VtkSurfaceMesh = nullptr;
  m_Pi = 3.14159265358979323846;
  m_ReadFromFile = true;
  m_Debug = false;
  m_FindingRealSolution = true;
  m_SurfaceMesh = nullptr;
  m_Image = nullptr;
  m_SphereImage = nullptr;
  for (int i = 0; i < 7; i++)
  {
    m_PoleElementsGN[i] = 0;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool
FEMConformalMap<TSurface, TImage, TDimension>::GenerateSystemFromVtkSurfaceMesh()
{
  if (!m_VtkSurfaceMesh)
  {
    ::std::cout << " NO MESH";
    return false;
  }
  ::std::cout << " Generate system from surface mesh " << std::endl;

  vtkTriangleFilter * fltTriangle = vtkTriangleFilter::New();
  fltTriangle->SetInput(m_VtkSurfaceMesh);
  cout << "   converting mesh to triangles " << endl;
  fltTriangle->Update();
  // Clean the data
  vtkCleanPolyData * fltCleaner = vtkCleanPolyData::New();
  fltCleaner->SetInput(fltTriangle->GetOutput());
  fltCleaner->SetTolerance(0);
  fltCleaner->ConvertPolysToLinesOn();
  cout << "   cleaning up triangle mesh " << endl;
  fltCleaner->Update();
  // Go through and delete the cells that are of the wrong type
  vtkPolyData * clean = fltCleaner->GetOutput();
  for (vtkIdType i = clean->GetNumberOfCells(); i > 0; i--)
  {
    if (clean->GetCellType(i - 1) != VTK_TRIANGLE)
    {
      clean->DeleteCell(i - 1);
    }
  }
  clean->BuildCells();
  m_VtkSurfaceMesh = clean;

  // create a material
  // Choose the material properties
  typename MaterialType::Pointer m;
  m = MaterialType::New();
  m->GN = 0;      // Global number of the material ///
  m->E = m_Sigma; // Young modulus -- used in the membrane ///
  m->A = 0.0;     // Crossection area ///
  m->h = 0.0;     // Crossection area ///
  m->I = 0.0;     // Moment of inertia ///
  m->nu = 0.;     // .0;    // poissons -- DONT CHOOSE 1.0!!///
  m->RhoC = 0.0;

  // Create the element type
  ElementType::Pointer e1 = ElementType::New();
  e1->m_mat = dynamic_cast<MaterialType *>(m);

  vtkPoints * vtkpoints = m_VtkSurfaceMesh->GetPoints();
  int         numPoints = vtkpoints->GetNumberOfPoints();
  int         foundnum = 0;
  for (int i = 0; i < numPoints; i++)
  {
    double *                   pt = vtkpoints->GetPoint(i);
    typename NodeType::Pointer n;
    n = new NodeType(pt[0], pt[1], pt[2]);
    n->GN = i;
    m_Solver.node.push_back(itk::fem::FEMP<NodeType>(n));
    foundnum++;
  }

  ::std::cout << " foundnum " << foundnum << std::endl;

  typename NodeType::Pointer * narr = new NodeType::Pointer[numPoints];
  for (int i = 0; i < numPoints; i++)
  {
    narr[i] = m_Solver.node.Find(i);
  }

  vtkCellArray * vtkcells = m_VtkSurfaceMesh->GetPolys();

  vtkIdType     npts;
  vtkIdType *   pts;
  unsigned long i = 0;
  unsigned long toti = vtkcells->GetNumberOfCells();
  //  ::std::cout << " progress ";
  ::std::cout << " DONE: NUMBER OF CELLS start " << toti << std::endl;
  for (vtkcells->InitTraversal(); vtkcells->GetNextCell(npts, pts);)
  {
    ElementType::Pointer e;
    e = dynamic_cast<ElementType *>(e1->Clone());

    e->SetNode(2, narr[pts[0]]);
    e->SetNode(1, narr[pts[1]]);
    e->SetNode(0, narr[pts[2]]);

    e->GN = i;
    m_Solver.el.push_back(itk::fem::FEMP<itk::fem::Element>(e));
    i++;
  }

  ::std::cout << " DONE: NUMBER OF CELLS " << i << std::endl;

  return true;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
bool
FEMConformalMap<TSurface, TImage, TDimension>::GenerateSystemFromSurfaceMesh()
{
  if (!m_SurfaceMesh)
  {
    return false;
  }

  // create a material
  // Choose the material properties
  typename MaterialType::Pointer m;
  m = MaterialType::New();
  m->GN = 0;      // Global number of the material ///
  m->E = m_Sigma; // Young modulus -- used in the membrane ///
  m->A = 0.0;     // Crossection area ///
  m->h = 0.0;     // Crossection area ///
  m->I = 0.0;     // Moment of inertia ///
  m->nu = 0.;     // .0;    // poissons -- DONT CHOOSE 1.0!!///
  m->RhoC = 0.0;
  // 3.e-4 prolly good
  //  ::std::cout << " input E ";  std::cin >> m->E;

  // Create the element type
  ElementType::Pointer e1 = ElementType::New();
  e1->m_mat = dynamic_cast<MaterialType *>(m);

  PointType * pt_ptr;
  PointType   pt;
  pt_ptr = &pt;
  for (int i = 0; i < m_SurfaceMesh->GetNumberOfPoints(); i++)
  {
    m_SurfaceMesh->GetPoint(i, pt_ptr);
    if (m_Debug)
    {
      ::std::cout << "Point: " << i << " coords " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
    }

    // turn the point into a node
    {
      // Create nodes
      typename NodeType::Pointer n;
      n = new NodeType(pt[0], pt[1], pt[2]);
      n->GN = i;
      m_Solver.node.push_back(itk::fem::FEMP<NodeType>(n));
    }
  }

  InputCellsContainerPointer  cellList = m_SurfaceMesh->GetCells();
  InputCellsContainerIterator cells = cellList->Begin();
  for (int i = 0; i < m_SurfaceMesh->GetNumberOfCells(); i++)
  {
    const unsigned long *            tp;
    typename SurfaceType::CellType * cellPtr = cells.Value();
    tp = cellPtr->GetPointIds();
    if (m_Debug)
    {
      ::std::cout << " pt 1 id " << tp[0] << " pt 2 id " << tp[1] << " pt 3 id " << tp[2] << std::endl;
    }
    cells++;

    // turn the cell into an element
    {
      // Create elements
      ElementType::Pointer e;
      e = dynamic_cast<ElementType *>(e1->Clone());
      e->SetNode(0, m_Solver.node.Find(tp[0]));
      e->SetNode(1, m_Solver.node.Find(tp[1]));
      e->SetNode(2, m_Solver.node.Find(tp[2]));
      e->GN = i;
      m_Solver.el.push_back(itk::fem::FEMP<itk::fem::Element>(e));
    }
  }

  return true;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::MapImageToSphere(TImage * img, float rad)
{
  typename ImageType::SizeType imageSize = img->GetLargestPossibleRegion().GetSize();

  vnl_vector<float> Pos; // solution at the point
  vnl_vector<float> Sol; // solution at the local point
  vnl_vector<float> Gpt; // global position given by local point

  Sol.set_size(ImageDimension);
  Gpt.set_size(ImageDimension);
  Pos.set_size(ImageDimension);

  VectorType disp;
  disp.set_size(3);
  disp.fill(0.0);

  float WIDTH = rad * 2.1;
  float HEIGHT = rad * 2.1;
  float DEPTH = rad * 2.1;

  typename ImageType::SizeType size;
  size[0] = WIDTH;
  size[1] = HEIGHT;
  size[2] = DEPTH;

  typename ImageType::IndexType start;
  start.Fill(0);

  typename ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  m_SphereImage = ImageType::New();

  m_SphereImage->SetRegions(region);
  m_SphereImage->AllocateInitialized();

  int center = (int)WIDTH / 2;

  typename ImageType::PixelType backgroundValue = 0;

  itk::ImageRegionIteratorWithIndex<ImageType> it(m_SphereImage, region);

  it.GoToBegin();
  typename ImageType::IndexType rindex = it.GetIndex();
  typename ImageType::IndexType sindex = it.GetIndex();

  while (!it.IsAtEnd())
  {
    rindex = it.GetIndex();
    float dist = 0;
    for (int ii = 0; ii < ImageDimension; ii++)
    {
      dist += (float)(rindex[ii] - center) * (float)(rindex[ii] - center);
    }
    dist = sqrt(dist);
    if (dist > rad)
    {
      it.Set(backgroundValue);
    }
    else
    {
      it.Set(1);
    }
    ++it;
  }

  PointType * pt_ptr;
  PointType   pt;
  pt_ptr = &pt;

  bool inimage;
  for (itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin(); elt != m_Solver.el.end(); ++elt)
  {
    unsigned int Nnodes = (*elt)->GetNumberOfNodes();
    inimage = true;
    for (unsigned int nd = 0; nd < Nnodes; nd++)
    {
      for (unsigned int f = 0; f < ImageDimension; f++)
      {
        float x;
        if ((*elt)->GetNumberOfDegreesOfFreedomPerNode() == 3)
        {
          Sol[f] = m_Solver.GetLinearSystemWrapper()->GetSolutionValue((*elt)->GetNode(nd)->GetDegreeOfFreedom(f), 0);
          x = (*elt)->GetNode(nd)->GetCoordinates()[f];
        }
        else
        {
          Sol[f] = (*elt)->GetNode(nd)->GetCoordinates()(f);
          m_SurfaceMesh->GetPoint((*elt)->GetNode(nd)->GN, pt_ptr);
          x = pt[f];
        }

        //  get img index
        long int temp;
        if (x != 0)
        {
          temp = (long int)((x) + 0.5);
        }
        else
        {
          temp = 0; // round
        }
        rindex[f] = temp;

        // get sphr index
        disp[f] = (float)((float)center) + rad * Sol[f];
        x = disp[f];
        if (x != 0)
        {
          temp = (long int)((x) + 0.5);
        }
        else
        {
          temp = 0; // round
        }
        sindex[f] = temp;
      }
      if (inimage)
      {
        m_SphereImage->SetPixel(sindex, img->GetPixel(rindex) + 1.0); // *100. +110.0);
      }
    }
  } // end of elt array loop
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::MapCheckerboardToImage(float increment)
{
  typedef float Float;

  if (increment < 0)
  {
    increment = 0;
  }
  if (increment > 1)
  {
    increment = 1;
  }

  float phiinc = m_Pi * increment;
  float thetainc = m_Pi * increment;

  vnl_vector<Float> Sol; // solution at the local point

  typename ImageType::SpacingType spacing = m_Image->GetSpacing();
  typename ImageType::SizeType    imageSize = m_Image->GetLargestPossibleRegion().GetSize();
  ::std::cout << " size " << imageSize << std::endl;
  VectorType disp;
  disp.set_size(3);
  disp.fill(0.0);

  //  Float solval,posval;
  bool                                    inimage;
  float                                   phi;
  float                                   theta;
  float                                   s;
  float                                   t;
  float                                   gridval;
  ImageRegionIteratorWithIndex<ImageType> wimIter(m_Image, m_Image->GetLargestPossibleRegion());
  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    wimIter.Set(0.0);
  }

  typename ImageType::IndexType rindex = wimIter.GetIndex();

  Sol.set_size(ImageDimension);

  PointType * pt_ptr;
  PointType   pt;
  pt_ptr = &pt;
  for (itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin(); elt != m_Solver.el.end(); ++elt)
  {
    unsigned int Nnodes = (*elt)->GetNumberOfNodes();

    inimage = true;
    for (unsigned int nd = 0; nd < Nnodes; nd++)
    {
      for (unsigned int f = 0; f < ImageDimension; f++)
      {
        Float x;

        Sol[f] = (*elt)->GetNode(nd)->GetCoordinates()(f);
        m_SurfaceMesh->GetPoint((*elt)->GetNode(nd)->GN, pt_ptr);
        x = pt[f];

        long int temp;
        if (x != 0)
        {
          temp = (long int)((x) + 0.5);
        }
        else
        {
          temp = 0; // round
        }
        rindex[f] = temp / spacing[f];
        disp[f] = (Float)1.0 * Sol[f];
        if (temp < 0 || temp > (long int)imageSize[f] - 1)
        {
          inimage = false;
        }
      }
      // convert the cartesian coordinates to theta and phi spherical coords
      phi = acos(Sol[2]);
      theta = acos(Sol[0] / sin(phi));
      s = phi / phiinc;
      t = theta / thetainc;
      gridval;
      if ((((int)t) + ((int)s)) % 2 == 0)
      {
        gridval = 100;
      }
      else
      {
        gridval = 200;
      }

      if (inimage)
      {
        m_Image->SetPixel(rindex, gridval); // *100. +110.0);
      }
    }
  } // end of elt array loop
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::MapStereographicCoordinatesToImage(int cdim)
{
  typedef float Float;

  vnl_vector<Float> Pos; // solution at the point
  vnl_vector<Float> Sol; // solution at the local point
  vnl_vector<Float> Gpt; // global position given by local point

  typename ImageType::SpacingType spacing = m_Image->GetSpacing();
  typename ImageType::SizeType    imageSize = m_Image->GetLargestPossibleRegion().GetSize();
  VectorType                      disp;
  disp.set_size(3);
  disp.fill(0.0);

  bool inimage;

  ImageRegionIteratorWithIndex<ImageType> wimIter(m_Image, m_Image->GetLargestPossibleRegion());
  wimIter.GoToBegin();
  for (; !wimIter.IsAtEnd(); ++wimIter)
  {
    wimIter.Set(0.0);
  }

  typename ImageType::IndexType rindex = wimIter.GetIndex();

  Sol.set_size(ImageDimension);
  Gpt.set_size(ImageDimension);
  Pos.set_size(ImageDimension);

  PointType * pt_ptr;
  PointType   pt;
  pt_ptr = &pt;
  for (itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin(); elt != m_Solver.el.end(); ++elt)
  {
    unsigned int Nnodes = (*elt)->GetNumberOfNodes();

    inimage = true;
    for (unsigned int nd = 0; nd < Nnodes; nd++)
    {
      for (unsigned int f = 0; f < ImageDimension; f++)
      {
        Float x;
        if ((*elt)->GetNumberOfDegreesOfFreedomPerNode() == 3)
        {
          Sol[f] = m_Solver.GetLinearSystemWrapper()->GetSolutionValue((*elt)->GetNode(nd)->GetDegreeOfFreedom(f), 0);
          x = (*elt)->GetNode(nd)->GetCoordinates()[f];
        }
        else
        {
          Sol[f] = (*elt)->GetNode(nd)->GetCoordinates()(f);
          m_SurfaceMesh->GetPoint((*elt)->GetNode(nd)->GN, pt_ptr);
          x = pt[f];
        }

        long int temp;
        if (x != 0)
        {
          temp = (long int)((x) + 0.5);
        }
        else
        {
          temp = 0; // round
        }
        rindex[f] = temp / spacing[f];
        disp[f] = (Float)1.0 * Sol[f];
        if (temp < 0 || temp > (long int)imageSize[f] - 1)
        {
          inimage = false;
        }
      }
      if (inimage)
      {
        m_Image->SetPixel(rindex, disp[cdim] + 2.0); // *100. +110.0);
      }
    }
  } // end of elt array loop
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::FindPoles(int dim)
{
  // find node with min x coord and its elt -> use as north pole (apply forces there)
  // find node with max x coord and its elt -> use as south pole (fix solution to zero)

  if (dim < 0)
  {
    return;
  }

  vnl_vector<float> minx(3, 9.e9);
  vnl_vector<float> maxx(3, 0.0);

  itk::fem::Element::VectorType minv;
  itk::fem::Element::VectorType maxv;

  itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin();
  unsigned int                           Nnodes = (*elt)->GetNumberOfNodes();

  for (; elt != m_Solver.el.end(); ++elt)
  {
    for (unsigned int nd = 0; nd < Nnodes; nd++)
    {
      for (unsigned int dd = 0; dd < ImageDimension; dd++)
      {
        float x = ((*elt)->GetNode(nd)->GetCoordinates()[dd]);
        if (x <= minx[dd])
        {
          //      ::std::cout << " x " << x << " minx " << minx << std::endl;
          m_PoleElementsGN[dd * 2] = (*elt)->GN;
          minv = (*elt)->GetNode(nd)->GetCoordinates();
          minx[dd] = x;
        }
        else if (x >= maxx[dd])
        {
          //      ::std::cout << " x " << x << " maxx " << maxx << std::endl;
          m_PoleElementsGN[dd * 2 + 1] = (*elt)->GN;
          maxv = (*elt)->GetNode(nd)->GetCoordinates();
          maxx[dd] = x;
        }
      }
    }
  }
  for (int jj = 0; jj < 3; jj++)
  {
    ::std::cout << " Element " << m_PoleElementsGN[jj * 2] << " is north " << jj << std::endl;
    ::std::cout << " Element " << m_PoleElementsGN[jj * 2 + 1] << " is south " << jj << std::endl;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::FixPoles(int dim)
{
  //  if (dim > 3 || dim < -3 || dim == 0 ) return;

  if (dim == m_NorthPole)
  {
    return;
  }

  itk::fem::Element::ArrayType::iterator elt = m_Solver.el.begin();
  unsigned int                           dofs = (*elt)->GetNumberOfDegreesOfFreedomPerNode();

  {
    float                     fixval;
    itk::fem::LoadBC::Pointer l1;
    itk::fem::LoadBC::Pointer l2;
    itk::fem::LoadBC::Pointer l3;
    itk::fem::LoadBC::Pointer l4;
    itk::fem::LoadBC::Pointer l5;
    //  itk::fem::LoadBC::Pointer l6;

    int   absdim = abs(dim);
    float oc = 2.0;
    ::std::cout << " Fixing BC  " << absdim << std::endl;

    switch (absdim)
    {
      case 0:
      {
        for (int i = 0; i < dofs; i++)
        {
          fixval = 0;
          l1 = itk::fem::LoadBC::New();
          l1->m_element = m_Solver.el[m_PoleElementsGN[0]];
          l1->m_dof = i;
          l1->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l1));
        }
      }
      break;
      case 1:
      {
        for (int i = 0; i < dofs; i++)
        {
          fixval = 0;
          l1 = itk::fem::LoadBC::New();
          l1->m_element = m_Solver.el[m_PoleElementsGN[1]];
          l1->m_dof = i;
          l1->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l1));
        }
      }
      break;
      case 2:
      {
        for (int i = 0; i < dofs; i++)
        {
          if (m_FindingRealSolution)
          {
            fixval = 0;
          }
          else
          {
            fixval = oc;
          }
          l2 = itk::fem::LoadBC::New();
          l2->m_element = m_Solver.el[m_PoleElementsGN[2]];
          l2->m_dof = i;
          l2->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l2));
        }
      }
      break;
      case 3:
      {
        for (int i = 0; i < dofs; i++)
        {
          if (m_FindingRealSolution)
          {
            fixval = 0;
          }
          else
          {
            fixval = -1.0 * oc; // was -2.0
          }
          l3 = itk::fem::LoadBC::New();
          l3->m_element = m_Solver.el[m_PoleElementsGN[3]];
          l3->m_dof = i;
          l3->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l3));
        }
      }
      break;
      case 4:
      {
        for (int i = 0; i < dofs; i++)
        {
          if (m_FindingRealSolution)
          {
            fixval = oc;
          }
          else
          {
            fixval = 0.0;
          }
          l4 = itk::fem::LoadBC::New();
          l4->m_element = m_Solver.el[m_PoleElementsGN[4]];
          l4->m_dof = i;
          l4->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l4));
        }
      }
      break;
      case 5:
      {
        for (int i = 0; i < dofs; i++)
        {
          if (m_FindingRealSolution)
          {
            fixval = -1.0 * oc;
          }
          else
          {
            fixval = 0.0;
          }
          l5 = itk::fem::LoadBC::New();
          l5->m_element = m_Solver.el[m_PoleElementsGN[5]];
          l5->m_dof = i;
          l5->m_value = vnl_vector<double>(1, fixval);
          m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*l5));
        }
      }
      break;
    }
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::ApplyRealForces(int dim)
{


  if (dim > 6)
  {
    dim = 6;
  }
  // loop over all elements, applying forces to the nodes 0 and 1 (real forces are zero on the 3rd node)
  //  for( ::itk::fem::Solver::ElementArray::iterator e = m_Solver.el.begin(); e!=m_Solver.el.end(); e++)
  {
    itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[dim]];

    int na = 0;
    int nb = 1;
    int nc = 2;
    {
      VectorType A = (e)->GetNode(na)->GetCoordinates();
      VectorType B = (e)->GetNode(nb)->GetCoordinates();
      VectorType C = (e)->GetNode(nc)->GetCoordinates();
      VectorType BA = B - A;
      VectorType AC = A - C;
      VectorType CB = C - B;
      float      bamag = BA.magnitude();
      float      cbmag = CB.magnitude();
      float      acmag = AC.magnitude();

      if (bamag > cbmag && bamag > acmag)
      {
        na = 0;
        nb = 1;
        nc = 2;
      }
      if (cbmag > bamag && cbmag > acmag)
      {
        na = 1;
        nb = 2;
        nc = 0;
      }
      if (acmag > bamag && acmag > cbmag)
      {
        na = 2;
        nb = 0;
        nc = 1;
      }
    }

    VectorType A = (e)->GetNode(na)->GetCoordinates();
    VectorType B = (e)->GetNode(nb)->GetCoordinates();
    VectorType C = (e)->GetNode(nc)->GetCoordinates();
    VectorType BA = B - A;
    VectorType CA = C - A;
    VectorType CB = C - B;
    float      bamag = BA.magnitude();

    float theta = (CA[0] * BA[0] + CA[1] * BA[1] + CA[2] * BA[2]) / bamag * bamag;
    if (theta > 1)
    {
      ::std::cout << " theta " << theta << std::endl;
    }

    VectorType E = A + BA * theta;
    VectorType CE = C - E;

    // load node 0
    {
      itk::fem::LoadNode::Pointer load = itk::fem::LoadNode::New();
      load->m_pt = na;
      load->F.set_size(3);
      load->F.fill(-1.0 / bamag);
      load->m_element = (e);
      m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*load));
    }
    // load node 1
    {
      itk::fem::LoadNode::Pointer load = itk::fem::LoadNode::New();
      load->m_pt = nb;
      load->F.set_size(3);
      load->F.fill(1.0 / bamag);
      load->m_element = (e);
      m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*load));
    }
    return;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::ApplyImaginaryForces(int dim)
{
  if (dim > 6)
  {
    dim = 6;
  }
  // loop over all elements, applying forces to the nodes 0 and 1 (real forces are zero on the 3rd node)
  //  for( ::itk::fem::Solver::ElementArray::iterator e = m_Solver.el.begin(); e!=m_Solver.el.end(); e++)
  {
    itk::fem::Element::Pointer e = m_Solver.el[m_PoleElementsGN[dim]];

    /*
        VectorType A=(e)->GetNode(0)->GetCoordinates();
        VectorType B=(e)->GetNode(1)->GetCoordinates();
        VectorType C=(e)->GetNode(2)->GetCoordinates();
        VectorType BA =B-A;
        VectorType CA =C-A;
        float bamag=BA.magnitude();

        float theta=(CA[0]*BA[0]+CA[1]*BA[1]+CA[2]*BA[2])/bamag*bamag;
        if (theta > 1) ::std::cout << " theta " << theta << std::endl;

        VectorType E = A+BA*theta;
        VectorType CE =C-E;

        float cemag=CE.magnitude();
    */
    int na = 0;
    int nb = 1;
    int nc = 2;

    {
      VectorType A = (e)->GetNode(na)->GetCoordinates();
      VectorType B = (e)->GetNode(nb)->GetCoordinates();
      VectorType C = (e)->GetNode(nc)->GetCoordinates();
      VectorType BA = B - A;
      VectorType AC = A - C;
      VectorType CB = C - B;
      float      bamag = BA.magnitude();
      float      cbmag = CB.magnitude();
      float      acmag = AC.magnitude();

      if (bamag > cbmag && bamag > acmag)
      {
        na = 0;
        nb = 1;
        nc = 2;
      }
      if (cbmag > bamag && cbmag > acmag)
      {
        na = 1;
        nb = 2;
        nc = 0;
      }
      if (acmag > bamag && acmag > cbmag)
      {
        na = 2;
        nb = 0;
        nc = 1;
      }
    }

    VectorType A = (e)->GetNode(na)->GetCoordinates();
    VectorType B = (e)->GetNode(nb)->GetCoordinates();
    VectorType C = (e)->GetNode(nc)->GetCoordinates();
    VectorType BA = B - A;
    VectorType CA = C - A;
    VectorType CB = C - B;
    float      bamag = BA.magnitude();

    float theta = (CA[0] * BA[0] + CA[1] * BA[1] + CA[2] * BA[2]) / bamag * bamag;
    if (theta > 1 || theta < 0)
    {
      //      ::std::cout << " ERROR : theta from non-acute angle  " << theta << std::endl;
      //      throw;
      //      return;
    }
    VectorType E = A + BA * theta;
    VectorType CE = C - E;

    float cemag = CE.magnitude();

    // load node 0
    {
      itk::fem::LoadNode::Pointer load = itk::fem::LoadNode::New();
      load->m_pt = na;
      load->F.set_size(3);
      load->F.fill((1.0 - theta) / cemag);
      load->m_element = (e);
      m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*load));
    }
    // load node 1
    {
      itk::fem::LoadNode::Pointer load = itk::fem::LoadNode::New();
      load->m_pt = nb;
      load->F.set_size(3);
      load->F.fill((theta) / cemag);
      load->m_element = (e);
      m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*load));
    }
    // load node 2
    {
      itk::fem::LoadNode::Pointer load = itk::fem::LoadNode::New();
      load->m_pt = nc;
      load->F.set_size(3);
      load->F.fill(-1.0 / cemag);
      load->m_element = (e);
      m_Solver.load.push_back(itk::fem::FEMP<itk::fem::Load>(&*load));
    }

    return;
  }
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::ComputeStereographicCoordinates()
{
  ::std::cout << " coords on unit sphere : " << std::endl;
  float radsq;      // radius squared
  float c1, c2, c3; // stereographic coordinates

  unsigned int dof;

  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    dof = (*n)->GetDegreeOfFreedom(0);
    float x = m_RealSolution[dof];
    float y = m_ImagSolution[dof];
    radsq = x * x + y * y;
    c1 = 2.0 * x / (1.0 + radsq);
    c2 = 2.0 * y / (1.0 + radsq);
    c3 = 2.0 * radsq / (1.0 + radsq) - 1.0;
    //      ::std::cout << c1 << "   " << c2 << "   " << c3 << " mag " << sqrt(c1*c1+c2*c2+c3*c3) << std::endl;

    VectorType coord = (*n)->GetCoordinates();
    coord[0] = c1;
    coord[1] = c2;
    coord[2] = c3;
    // /      (*n)->SetCoordinates(coord);
  }
  ::std::cout << " coords on unit sphere done " << std::endl;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::ConformalMap()
{
  m_Solver.load.clear();
  m_Solver.node.clear();
  m_Solver.el.clear();

  /**
   * Open the file and assign it to stream object f
   */
  if (m_ReadFromFile)
  {
    const char * filename = m_ParameterFileName.c_str();
    ::std::cout << "Reading FEM problem from file: " << std::string(filename) << "\n";
    std::ifstream f;
    f.open(filename);
    if (!f)
    {
      ::std::cout << "File " << filename << " not found!\n";
      return;
    }

    try
    {
      m_Solver.Read(f);
    }
    catch (const ::itk::fem::FEMException & e)
    {
      ::std::cout << "Error reading FEM problem: " << filename << "!\n";
      e.Print(::std::cout);
      return;
    }

    f.close();
  }
  else // if (m_VtkSurfaceMesh)
  {
    this->GenerateSystemFromVtkSurfaceMesh();
  }
  //   else this->GenerateSystemFromSurfaceMesh();

  int dir = m_NorthPole;
  //  ::std::cout << " input dir to fix "; std::cin >> dir;
  this->FindPoles(dir);
  if (m_NorthPole != m_SouthPole)
  {
    this->FixPoles(m_SouthPole);
  }
  //  for (int oo=0; oo<6; oo++)
  //   this->FixPoles(oo);

  /**
   * Assign a unique id (global freedom number - GFN)
   * to every degree of freedom (DOF) in a system.
   */
  m_Solver.GenerateGFN();
  m_ImagSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom());
  m_RealSolution.set_size(m_Solver.GetNumberOfDegreesOfFreedom());
  m_RealSolution.fill(0);
  m_ImagSolution.fill(0);

  if (m_Debug)
  {
    for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
    {
      ::std::cout << "Node#: " << (*n)->GN << ": ";
      ::std::cout << " coord " << (*n)->GetCoordinates() << std::endl;
    }
    for (::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); ++n)
    {
      ::std::cout << "Elt#: " << (*n)->GN << ": has " << (*n)->GetNumberOfNodes() << " nodes ";
      for (int i = 0; i < (*n)->GetNumberOfNodes(); i++)
      {
        ::std::cout << " coord " << (*n)->GetNode(i)->GetCoordinates() << std::endl;
      }
    }
  }

  unsigned int maxits = m_Solver.GetNumberOfDegreesOfFreedom(); // should be > twice ndofs
  // if (m_Debug)
  ::std::cout << " ndof " << maxits << std::endl;
  itpackWrapper.SetMaximumNumberIterations(maxits * 4);
  itpackWrapper.SetTolerance(1.e-5);
  itpackWrapper.JacobianConjugateGradient();
  itpackWrapper.SetMaximumNonZeroValuesInMatrix(maxits * 100);
  m_Solver.SetLinearSystemWrapper(&itpackWrapper);

  /**
   * Assemble the master stiffness matrix. In order to do this
   * the GFN's should already be assigned to every DOF.
   */
  ::std::cout << " assemble k " << std::endl;
  m_Solver.AssembleK();

  ::std::cout << " assemble k done " << std::endl;
  if (m_Debug)
  {
    for (int i = 0; i < m_Solver.GetNumberOfDegreesOfFreedom(); i++)
    {
      float sum = 0;
      for (int j = 0; j < m_Solver.GetNumberOfDegreesOfFreedom(); j++)
      {
        float val = m_Solver.GetLinearSystemWrapper()->GetMatrixValue(i, j);
        if (i != j)
        {
          sum += val;
        }
        //    if (val != 0) ::std::cout << " entry i,j " << i << " , " << j << "  =  " <<
        // m_Solver.GetLinearSystemWrapper()->GetMatrixValue(i,j)<<std::endl;
      }
      ::std::cout << " sum #" << i << " = " << sum << std::endl;
    }
  }
  /**
   * Invert the master stiffness matrix
   */
  m_Solver.DecomposeK();
  /**
   * Assemble the master force vector (from the applied loads)
   */
  // run this once for imaginary and once for real
  for (int im = 0; im < 2; im++)
  {
    if (im == 1)
    {
      m_FindingRealSolution = false;
      m_Solver.load.clear();
      if (m_NorthPole != m_SouthPole)
      {
        this->FixPoles(m_SouthPole);
      }
      //       for (int oo=0; oo<6; oo++)
      //         this->FixPoles(oo);
    }

    if (im == 0)
    {
      this->ApplyRealForces(m_NorthPole);
    }
    else
    {
      this->ApplyImaginaryForces(m_NorthPole);
    }
    m_Solver.AssembleF();
    ::std::cout << " force done " << std::endl;
    m_Solver.Solve();
    ::std::cout << " solve done " << std::endl;
    m_Solver.UpdateDisplacements(); // copies solution to nodes
    m_Debug = true;
    if (m_Debug)
    {
      ::std::cout << "\nNodal displacements:\n";
    }
    unsigned int ct = 0;
    for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
    {
      if (m_Debug)
      {
        ::std::cout << "Node#: " << (*n)->GN << ": ";
      }
      if (m_Debug)
      {
        ::std::cout << " coord " << (*n)->GetCoordinates() << std::endl;
      }
      // For each DOF in the node... */
      for (unsigned int d = 0, dof;
           (dof = (*n)->GetDegreeOfFreedom(d)) != ::itk::fem::Element::InvalidDegreeOfFreedomID;
           d++)
      {
        if (m_Debug)
        {
          ::std::cout << m_Solver.GetSolution(dof);
        }
        if (m_Debug)
        {
          ::std::cout << ",  ";
        }
        //         if (d==0 && im == 0) {m_RealSolution[ct]=m_Solver.GetSolution(dof); ct++;}
        //         if (d==0 && im == 1) {m_ImagSolution[ct]=m_Solver.GetSolution(dof); ct++;}
        if (im == 0)
        {
          m_RealSolution[dof] = m_Solver.GetSolution(dof);
          ct++;
        }
        if (im == 1)
        {
          m_ImagSolution[dof] = m_Solver.GetSolution(dof);
          ct++;
        }
        //         ::std::cout << " dof " << dof << " ct " << ct << std::endl;
      }
      if (m_Debug)
      {
        ::std::cout << "\b\b\b \b\n";
      }
      m_Debug = false;
    }
  }

  return;
}

template <typename TSurface, typename TImage, unsigned int TDimension>
void
FEMConformalMap<TSurface, TImage, TDimension>::BuildOutputMeshes(typename TImage::Pointer image)
{
  // Get the number of points in the mesh
  int numPoints = m_Solver.node.size();

  // Create vtk polydata
  vtkPolyData * polydata1 = vtkPolyData::New();
  vtkPolyData * polydata2 = vtkPolyData::New();

  // Create the vtkPoints object and set the number of points
  vtkPoints * vpoints1 = vtkPoints::New();
  vtkPoints * vpoints2 = vtkPoints::New();

  vpoints1->SetNumberOfPoints(numPoints);
  vpoints2->SetNumberOfPoints(numPoints);

  vtkFloatArray * param = vtkFloatArray::New();
  param->SetName("angle");
  ::std::cout << " start pts ";
  int idx = 0;

  vtkDataArray * scs = nullptr;
  if (m_VtkSurfaceMesh)
  {
    vtkPointData * pd = m_VtkSurfaceMesh->GetPointData();
    scs = pd->GetScalars();
  }
  typename TImage::SpacingType spacing = image->GetSpacing();
  for (::itk::fem::Solver::NodeArray::iterator n = m_Solver.node.begin(); n != m_Solver.node.end(); ++n)
  {
    // extrinisic coords
    VectorType loc = (*n)->GetCoordinates();

    // spherical coords
    unsigned int dof = (*n)->GetDegreeOfFreedom(0);
    float        x = m_RealSolution[dof];
    float        y = m_ImagSolution[dof];
    float        radsq = x * x + y * y;
    float        c1 = 2.0 * x / (1.0 + radsq);
    float        c2 = 2.0 * y / (1.0 + radsq);
    float        c3 = 2.0 * radsq / (1.0 + radsq) - 1.0;
    if (idx % 1000 == 0)
    {
      ::std::cout << c1 << "   " << c2 << "   " << c3 << " mag " << sqrt(c1 * c1 + c2 * c2 + c3 * c3) << std::endl;
    }

    //
    float pt1[3];
    pt1[0] = c1;
    pt1[1] = c2;
    pt1[2] = c3;
    float pt2[3];
    pt2[0] = loc[0];
    pt2[1] = loc[1];
    pt2[2] = loc[2];
    vpoints1->SetPoint(idx, pt1);
    vpoints2->SetPoint(idx, pt2);

    typename TImage::IndexType index;
    index[0] = (long int)loc[0] / spacing[0];
    index[1] = (long int)loc[1] / spacing[1];
    index[2] = (long int)loc[2] / spacing[2];
    float temp;
    if (m_VtkSurfaceMesh)
    {
      temp = scs->GetTuple1(idx);
    }
    else
    {
      temp = (float)image->GetPixel(index);
    }

    param->InsertNextValue(temp);

    (*n)->GN = idx;
    idx++;
  }

  ::std::cout << " done with pts " << std::endl;
  vtkCellArray * tris1 = vtkCellArray::New();
  vtkCellArray * tris2 = vtkCellArray::New();

  ::std::cout << " start with tris " << std::endl;
  for (::itk::fem::Solver::ElementArray::iterator n = m_Solver.el.begin(); n != m_Solver.el.end(); ++n)
  {
    tris1->InsertNextCell(3);
    tris2->InsertNextCell(3);
    for (int i = 0; i < (*n)->GetNumberOfNodes(); i++)
    {
      tris1->InsertCellPoint((*n)->GetNode(i)->GN);
      tris2->InsertCellPoint((*n)->GetNode(i)->GN);
    }
  }
  ::std::cout << " done with tris " << std::endl;
  // Assign points and cells
  polydata1->SetPoints(vpoints1);
  polydata2->SetPoints(vpoints2);
  polydata1->SetPolys(tris1);
  polydata2->SetPolys(tris2);
  polydata1->GetPointData()->SetScalars(param);
  polydata2->GetPointData()->SetScalars(param);

  //  vtkDelaunay2D* delny1=vtkDelaunay2D::New();
  // delny1->SetInput(polydata1);
  // m_ExtractedSurfaceMesh=delny1->GetOutput();//polydata1;//
  m_ExtractedSurfaceMesh = vtkGenSmoothMesh(polydata1);

  //  vtkDelaunay2D* delny2=vtkDelaunay2D::New();
  // delny2->SetInput(polydata2);
  // m_DiskSurfaceMesh=delny2->GetOutput();
  if (!m_VtkSurfaceMesh)
  {
    m_VtkSurfaceMesh = vtkGenSmoothMesh(polydata2); // polydata2;
  }
  return;
}
} // namespace itk

#endif
